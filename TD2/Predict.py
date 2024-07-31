import os
import sys
import gzip
import copy
import time
import pandas
import pickle
import argparse
import warnings
import subprocess
import pandas
import numpy as np
from tqdm.auto import tqdm
from io import StringIO
from collections import defaultdict

from TD2.translator import Translator
from TD2.LongOrfs import load_fasta


def get_args():
    parser = argparse.ArgumentParser()
    
    # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("-t", dest="transcripts",  type=str, required=True, help="REQUIRED path to transcripts.fasta")
    
    # optional
    parser.add_argument("-P", dest="psauron_cutoff", type=float, required=False, help="minimum in-frame PSAURON score required to report ORF, assuming no homology hits (range: [0,1]; default: 0.25)", default=0.25)
    parser.add_argument("--retain_mmseqs_hits",  type=str, required=False, help="mmseqs output in '.m8' format. Any ORF with a MMseqs2 match will be retained in the final output.")
    parser.add_argument("--retain_blastp_hits",  type=str, required=False, help="blastp output in '-outfmt 6' format. Any ORF with a blastp match will be retained in the final output.")
    parser.add_argument("--retain_hmmer_hits",  type=str, required=False, help="domain table output file from running hmmer to search Pfam. Any ORF with a Pfam domain hit will be retained in the final output.")
    parser.add_argument("--retain_long_orfs_length",  type=int, required=False, help="retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence marks it as coding (default: 1000000, so essentially turned off by default.)", default=1000000)
    parser.add_argument("--single_best_only",  type=str, required=False, help="Retain only the single best orf per transcript (prioritized by homology then orf length)")
    parser.add_argument("-O", dest="output_dir", type=str, required=False, help="same output directory from LongOrfs", default="./transcripts.TD2_dir")
    parser.add_argument("-G", dest="genetic_code", type=int, required=False, help="genetic code a.k.a. translation table, NCBI integer codes, default=1", default=1)
    
    # TODO verbosity
    parser.add_argument("-v", "--verbose", action='store_true', help="set -v for verbose output with progress bars, default=False", default=False)


    # TODO start codon refinement with PWM? Is there anything better?

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help']) # prints help message if no args are provided by user
    return args

def find_encapsulated_intervals(intervals):
    # O(n log(n))
    # Assumes intervals is a list of tuples with names: e.g. [("foo", (1,100)), ("bar", (0,100))],
    # Sort intervals: first by low coord (ascending), then by high coord (descending)
    intervals.sort(key=lambda x: (x[1][0], -x[1][1]))
    
    # Set of encapsulated intervals
    encapsulated_intervals = set()
    
    # Traverse the sorted list
    bigboi_low = float('inf')
    bigboi_high = -float('inf')
    for x in intervals:
        name = x[0]
        smolboi_low, smolboi_high = x[1]
        if smolboi_high <= bigboi_high:
            # encapsulated
            encapsulated_intervals.add(name)
        else:
            # not encapsulated, new bigboi
            bigboi_low = smolboi_low
            bigboi_high = smolboi_high
            
    return encapsulated_intervals
    
def main():
    # supress annoying warnings
    warnings.filterwarnings('ignore')
    
    # parse command line arguments
    args = get_args()
    psauron_cutoff = args.psauron_cutoff
    
    # use absolute path of output
    output_dir = os.path.abspath(args.output_dir)    
    
    # run psauron to score ORFs
    print(f"Step 1: Running PSAURON", flush=True)
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    p_score = os.path.join(output_dir, "psauron_score.csv")
    if args.verbose:
        command_psauron = ["psauron", "-i", str(p_cds), "-o", str(p_score), "-m", "0", "--inframe", str(psauron_cutoff), "-v"]
        #result_psauron = subprocess.run(command_psauron, capture_output=False, text=True)
    else:
        command_psauron = ["psauron", "-i", str(p_cds), "-o", str(p_score), "-m", "0", "--inframe", str(psauron_cutoff)]
        #result_psauron = subprocess.run(command_psauron, capture_output=True, text=True)
    
    # load psauron results
    df_psauron = pandas.read_csv(p_score, skiprows=3)
    name_to_score = dict(zip([str(x.split(" ")[0]) for x in df_psauron["description"].tolist()], 
                              df_psauron["in_frame_score"]))
    
    # select transcripts based on psauron score
    df_psauron_selected = df_psauron[df_psauron.apply(lambda row: row['in_frame_score'] > psauron_cutoff and all(row[3:] < row['in_frame_score']), axis=1)]
    name_psauron_selected = set([str(x.split(" ")[0]) for x in df_psauron_selected["description"]])
    
    print(f"Done.")
    
    
    # integrate homology search results
    if any([args.retain_mmseqs_hits, args.retain_blastp_hits, args.retain_hmmer_hits]):
        print(f"Step 2: Integrating homology results", flush=True)
    else:
        print(f"Step 2: No homology search results provided, skipping", flush=True)
    
    # parse mmseqs
    hits_mmeseqs = set()
    if args.retain_mmseqs_hits:
        print(f"Loading MMseqs2 output", flush=True)
        p_mmseqs = args.retain_mmseqs_hits
        df_mmseqs = pandas.read_table(p_mmseqs, header=None)
        hits_mmeseqs = set(df_mmseqs[0])
        print(f"Found {int(len(hits_mmeseqs)):d} ORFs with MMseqs2 hits", flush=True)
    
    # parse blast
    hits_blastp = set()
    if args.retain_blastp_hits:
        print(f"Loading blastp output", flush=True)
        p_blastp = args.retain_blastp_hits
        df_blastp = pandas.read_table(p_blastp, header=None)
        hits_blastp = set(df_blastp[0])
        print(f"Found {int(len(hits_blastp)):d} ORFs with blastp hits", flush=True)   
        
    # parse hmmer
    hits_hmmer = set()
    if args.retain_hmmer_hits:
        print(f"Loading hmmer output", flush=True)
        p_hmmer = args.retain_hmmer_hits
        # remove lines that start with "#"
        with open(p_hmmer, 'r') as f:
            filtered_lines = [line for line in f if not line.startswith("#")]
        f_filtered = StringIO(''.join(filtered_lines))
        df_hmmer = pandas.read_table(f_filtered, header=None)
        hits_hmmer = set(df_hmmer[0])
        print(f"Found {int(len(hits_hmmer)):d} ORFs with hmmer hits", flush=True)  
    print(f"Done.")
    
    
    # generate final ORfs
    print(f"Step 3: Generating final ORFs", flush=True)
    
    # load LongOrfs output
    print(f"Loading LongOrfs output from {output_dir}", flush=True)
    p_pep = os.path.join(output_dir, "longest_orfs.pep")
    p_gff3 = os.path.join(output_dir, "longest_orfs.gff3")
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    
    # get full description lines to retrieve intervals
    with open(p_pep, "rt") as f:
        pep_description_list = [x for x in f.readlines() if x.startswith(">")]
    pep_name_list = []
    pep_transcript_list = []
    pep_strand_list = []
    pep_lowcoord_list = []
    pep_highcoord_list = []
    for d in pep_description_list:
        l = d.split(" ")
        name = str(l[0][1:])
        transcript = str(l[-1].split(":")[0])
        coords = str(l[-1].split(":")[1][:-4])
        coords_int = [int(x) for x in coords.split("-")]
        lowcoord = min(coords_int)
        highcoord = max(coords_int)
        
        pep_name_list.append(name)
        pep_transcript_list.append(transcript)
        pep_lowcoord_list.append(lowcoord)
        pep_highcoord_list.append(highcoord)
    
    # group intervals by transcript
    transcript_intervals = defaultdict(list)
    for i, transcript in enumerate(pep_transcript_list):
        transcript_intervals[transcript].append((pep_name_list[i],
                                                 (pep_lowcoord_list[i],
                                                  pep_highcoord_list[i])))
    
    # find fully encapsulated ORFs
    name_encapsulated = set()
    for transcript, intervals in transcript_intervals.items():
        encapsulated = find_encapsulated_intervals(intervals)
        name_encapsulated.update(encapsulated)
    print(f"Removing {len(name_encapsulated):d} encapsulated ORFs", flush=True)
    
    # write final outputs to current working directory
    basename = os.path.basename(args.transcripts)
    p_pep_final = basename + ".TD2.pep"
    p_gff3_final = basename + ".TD2.gff3"
    p_cds_final = basename + ".TD2.cds"
    p_bed_final = basename + ".TD2.bed" # TODO bed file, probably better just to write gff then convert
    
    # pep fasta
    with open(p_pep, "rt") as f_pep, open(p_pep_final, "wt") as f_pep_final:
        longorfs_pep = f_pep.readlines()
        description_list = []
        seq_list = []
        for line in longorfs_pep:
            if line.startswith(">"):
                description_list.append(line[1:])
            else:
                seq_list.append(line)
        # get ORF information
        for i, ORF in enumerate(description_list):
            s = ORF.split(" ")
            name = str(s[0])
            # discard encapsulated
            if name in name_encapsulated:
                continue
            # filter with psauron and homology
            if not ((name in name_psauron_selected) or (name in hits_mmeseqs) or (name in hits_blastp) or (name in hits_hmmer)):
                continue
            gene = str(".".join(s[0].split(".")[:-1])) # TODO do we need to get gene name from the tab delimited file?
            ORF_type = str(s[1].split(":")[1])
            psauron_score = '{:.3f}'.format(round(float(name_to_score[name]), 3))
            length = str(s[2].split(":")[1])
            location = s[-1].rstrip()
            strand = location[-3:]
            # match TransDecoder output format
            description_line_final = ">" + name + " GENE." + gene + "~~" + name + "  ORF type:" + ORF_type + " " + strand + ",psauron_score=" + psauron_score + " len:" + length + " " + location + "\n"
            f_pep_final.write(description_line_final)
            seq = seq_list[i]
            while len(seq) > 60:
                f_pep_final.write(str(seq[:60]) + "\n")
                seq = seq[60:]
            f_pep_final.write(str(seq))
            
        # gff3
            
            
    
    
    
    
    
    #with open(p_pep, "rt") as f_pep, open(p_gff3, "rt") as f_gff3, open(p_cds, "rt") as f_cds:
        # open final output files
     #   with open(p_pep_final, "wt") as f_pep_final, open(p_gff3_final, "wt") as f_gff3_final, open(p_cds_final, "wt") as f_cds_final:
     #       for 
     #       pass
    

    
if __name__ == "__main__":
    main()