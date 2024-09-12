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
from collections import OrderedDict

from TD2.translator import Translator
from TD2.LongOrfs import load_fasta


def get_args():
    parser = argparse.ArgumentParser()
    
    # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("-t", dest="transcripts",  type=str, required=True, help="REQUIRED path to transcripts.fasta")
    
    # optional
    parser.add_argument("-P", dest="psauron_cutoff", type=float, required=False, help="minimum in-frame PSAURON score required to report ORF assuming no homology hits, higher is less sensitive and more precise (range: [0,1]; default: 0.25)", default=0.25)
    parser.add_argument("--single-best-only", action='store_true', help="retain only the single best ORF per transcript (prioritized by homology then ORF length), default=False")
    parser.add_argument("--retain-mmseqs-hits", type=str, required=False, help="mmseqs output in '.m8' format. Complete ORFs with a MMseqs2 match will be retained in the final output.")
    parser.add_argument("--retain-blastp_hits", type=str, required=False, help="blastp output in '-outfmt 6' format. Complete ORFs with a blastp match will be retained in the final output.")
    parser.add_argument("--retain-hmmer_hits", type=str, required=False, help="domain table output file from running hmmer to search Pfam. Complete ORFs with a Pfam domain hit will be retained in the final output.")
    parser.add_argument("--retain-long-orfs-length", type=int, required=False, help="retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence marks it as coding (default: 1000000, so essentially turned off by default.)", default=1000000)
    parser.add_argument("--retain-encapsulated", action='store_true', help="retain ORFs that are fully contained within larger ORFs, default=False")
    parser.add_argument("--retain-partial", action='store_true', help="retain 5' and 3' partial ORFs (may cause correct complete ORFs to be missed), default=False")
    parser.add_argument("--psauron-all-frame", action='store_true', help="require ORF to have highest PSAURON score compared to all other reading frames, set this argument for less sensitive and more precise ORFs, default=False")

    parser.add_argument("-G", dest="genetic_code", type=int, required=False, help="genetic code a.k.a. translation table, NCBI integer codes, default=1", default=1)
    parser.add_argument("-O", dest="output_dir", type=str, required=False, help="same output directory from LongOrfs", default="./transcripts.TD2_dir")
    
    # TODO verbosity
    parser.add_argument("-v", "--verbose", action='store_true', help="verbose output with progress bars, default=False", default=False)

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
    if args.output_dir == "./transcripts.TD2_dir":
        p_transcripts = os.path.abspath(args.transcripts)
        output_dir = os.path.splitext(os.path.basename(p_transcripts))[0]
    else:
        output_dir = os.path.abspath(args.output_dir) 
    
    # run psauron to score ORFs
    print(f"Step 1: Running PSAURON", flush=True)
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    p_score = os.path.join(output_dir, "psauron_score.csv")
    if args.verbose:
        if args.psauron_all_frame:
            command_psauron = ["psauron", "-i", str(p_cds), "-o", str(p_score), "-m", "0", "--inframe", str(psauron_cutoff), "-v"]
        else:
            command_psauron = ["psauron", "-i", str(p_cds), "-o", str(p_score), "-m", "0", "--inframe", str(psauron_cutoff), "-v", "-s"]
        result_psauron = subprocess.run(command_psauron, capture_output=False, text=True)
    else:
        if args.psauron_all_frame:
            command_psauron = ["psauron", "-i", str(p_cds), "-o", str(p_score), "-m", "0", "--inframe", str(psauron_cutoff)]
        else:
            command_psauron = ["psauron", "-i", str(p_cds), "-o", str(p_score), "-m", "0", "--inframe", str(psauron_cutoff), "-s"]
        result_psauron = subprocess.run(command_psauron, capture_output=True, text=True)
    
    # load psauron results
    if args.psauron_all_frame:
        df_psauron = pandas.read_csv(p_score, skiprows=3)
        ID_to_score = dict(zip([str(x.split(" ")[0]) for x in df_psauron["description"].tolist()], 
                                df_psauron["in_frame_score"]))
        df_psauron_selected = df_psauron[df_psauron.apply(lambda row: row['in_frame_score'] > psauron_cutoff and all(row[3:] < row['in_frame_score']), axis=1)]
    else:
        df_psauron = pandas.read_csv(p_score, skiprows=2)
        ID_to_score = dict(zip([str(x.split(" ")[0]) for x in df_psauron["description"].tolist()], 
                                df_psauron["in-frame_score"])) # this "-" is a bug in the -s mode for psauron, TODO fix this in psauron
        df_psauron_selected = df_psauron[df_psauron.apply(lambda row: row['in-frame_score'] > psauron_cutoff, axis=1)]
    
    ID_psauron_selected = set([str(x.split(" ")[0]) for x in df_psauron_selected["description"]])
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
        if int(len(hits_mmeseqs)) == 1:
            print(f"Found {int(len(hits_mmeseqs)):d} ORF with MMseqs2 hit", flush=True)
        else:
            print(f"Found {int(len(hits_mmeseqs)):d} ORFs with MMseqs2 hits", flush=True)
    
    # parse blast
    hits_blastp = set()
    if args.retain_blastp_hits:
        print(f"Loading blastp output", flush=True)
        p_blastp = args.retain_blastp_hits
        df_blastp = pandas.read_table(p_blastp, header=None)
        hits_blastp = set(df_blastp[0])
        if int(len(hits_blastp)) == 1:
            print(f"Found {int(len(hits_blastp)):d} ORF with blastp hit", flush=True)
        else:
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
        if int(len(hits_hmmer)) == 1:
            print(f"Found {int(len(hits_hmmer)):d} ORF with hmmer hit", flush=True)  
        else:
            print(f"Found {int(len(hits_hmmer)):d} ORFs with hmmer hits", flush=True)  
    print(f"Done.")
    
    
    # generate final ORfs
    print(f"Step 3: Deciding upon final ORFs", flush=True)
    
    # load LongOrfs output
    print(f"Loading LongOrfs output from {output_dir}", flush=True)
    p_pep = os.path.join(output_dir, "longest_orfs.pep")
    p_gff3 = os.path.join(output_dir, "longest_orfs.gff3")
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    
    # get full description lines to retrieve intervals
    with open(p_pep, "rt") as f:
        pep_description_list = [x for x in f.readlines() if x.startswith(">")]
    pep_ID_list = []
    pep_transcript_list = []
    pep_strand_list = []
    pep_lowcoord_list = []
    pep_highcoord_list = []
    pep_ID_to_info = dict()
    ID_to_partial = dict()
    for d in pep_description_list:
        l = d.split(" ")
        ID = str(l[0][1:])
        transcript = str(":".join(l[-1].split(":")[:-1]))
        coords = str(l[-1].split(":")[-1][:-4])
        coords_int = [int(x) for x in coords.split("-")]
        lowcoord = min(coords_int)
        highcoord = max(coords_int)
        if ("5prime_partial" in d) or ("3prime_partial" in d):
            partial = True
        else:
            partial = False
        
        pep_ID_list.append(ID)
        pep_transcript_list.append(transcript)
        pep_lowcoord_list.append(lowcoord)
        pep_highcoord_list.append(highcoord)
        pep_ID_to_info[ID] = (transcript, lowcoord, highcoord)
        ID_to_partial[ID] = partial
    
    # write final outputs to current working directory
    basename = os.path.basename(args.transcripts)
    p_pep_final = basename + ".TD2.pep"
    p_gff3_final = basename + ".TD2.gff3"
    p_cds_final = basename + ".TD2.cds"
    p_bed_final = basename + ".TD2.bed" # TODO bed file, probably better just to write gff then convert
    
    # ORF selection
    ID_selected = set() # IDs of all transcripts that pass filters
    ID_to_info = dict()
    ID_to_description = dict()
    with open(p_pep, "rt") as f_pep:
        longorfs_pep = f_pep.readlines()
        description_list = []
        seq_list_pep = []
        for line in longorfs_pep:
            if line.startswith(">"):
                description_list.append(line[1:])
            else:
                seq_list_pep.append(line)
        # get ORF information
        for i, ORF in enumerate(description_list):
            s = ORF.split(" ")
            ID = str(s[0])
            ID_to_description[ID] = ORF
            # filter with psauron and homology
            if not ((ID in ID_psauron_selected) or (ID in hits_mmeseqs) or (ID in hits_blastp) or (ID in hits_hmmer)):
                continue
                
            # remove partial ORFs by default
            if ID_to_partial[ID] and not args.retain_partial:
                continue
            
            # get transcript info
            gene = str(".".join(s[0].split(".")[:-1])) # TODO do we need to get gene ID from the tab delimited file?
            ORF_type = str(s[1].split(":")[1])
            psauron_score = '{:.3f}'.format(round(float(ID_to_score[ID]), 3))
            length = str(s[2].split(":")[1])
            location = s[-1].rstrip()
            strand = location[-3:]
            
            # match TransDecoder output format
            description_line_final = ">" + ID + " GENE." + gene + "~~" + ID + "  ORF type:" + ORF_type + " " + strand + ",psauron_score=" + psauron_score + " len:" + length + " " + location + "\n"
            seq_pep = seq_list_pep[i].rstrip()
            
            # keep info of IDs that pass filters
            info = (description_line_final, seq_pep)
            ID_to_info[ID] = info
            ID_selected.add(ID)
    
    # read cds
    ID_to_cds = dict()
    with open(p_cds, "rt") as f_cds:
        longorfs_cds = f_cds.readlines()  
        for line in longorfs_cds:
            if line.startswith(">"):
                ID = line.split(" ")[0][1:]
            else:
                seq_cds = line.rstrip()
                ID_to_cds[ID] = seq_cds
    
    # remove encapsulated ORFs
    if not args.retain_encapsulated:
        # group intervals by transcript
        transcript_intervals = defaultdict(list)
        for ID in ID_selected:
            transcript, lowcoord, highcoord = pep_ID_to_info[ID]
            transcript_intervals[transcript].append((ID, (lowcoord, highcoord)))            
        
        # find fully encapsulated ORFs
        ID_encapsulated = set()
        for transcript, intervals in transcript_intervals.items():
            encapsulated = find_encapsulated_intervals(intervals)
            ID_encapsulated.update(encapsulated)
        if len(ID_encapsulated) == 1:
            print(f"Removing {len(ID_encapsulated):d} encapsulated ORF", flush=True)
        else:
            print(f"Removing {len(ID_encapsulated):d} encapsulated ORFs", flush=True)
        
        # remove from final set
        ID_selected -= ID_encapsulated    
    
    # keep only one ORF per transcript
    # prioritizes homology, then ORF length
    if args.single_best_only:
        transcript_to_ID_info = dict()
        ID_not_single_best = set()
        for ID in ID_selected:
            transcript, lowcoord, highcoord = pep_ID_to_info[ID]
            length = highcoord - lowcoord + 1
            homology = any([ID in hits_blastp,
                            ID in hits_hmmer,
                            ID in hits_mmeseqs])
            info = (ID, homology, length)
            # will keep longest ORF with any homology hit
            if transcript in transcript_to_ID_info:
                info_stored = transcript_to_ID_info[transcript]
                if info[1] and info_stored[1]: # both have homology
                    if info[2] > info_stored[2]: # keep longer ORF
                        transcript_to_ID_info[transcript] = info
                        ID_not_single_best.add(info_stored[0]) # remove shorter ORF
                    else:
                        ID_not_single_best.add(info[0]) # remove shorter ORF
                elif info[1] and not info_stored[1]: 
                    transcript_to_ID_info[transcript] = info # keep ORF with homology
                    ID_not_single_best.add(info_stored[0])# remove ORF with no homology
                elif not info[1] and  info_stored[1]:
                    ID_not_single_best.add(info[0])# remove ORF with no homology
                elif not info[1] and not info_stored[1]:
                    if info[2] > info_stored[2]: # keep longer ORF
                        transcript_to_ID_info[transcript] = info
                        ID_not_single_best.add(info_stored[0]) # remove shorter ORF
                    else:
                        ID_not_single_best.add(info[0]) # remove shorter ORF
            else:
                transcript_to_ID_info[transcript] = info
                
        # remove from final set
        ID_selected -= ID_not_single_best   
    
    print(f"Writing final output to current working directory", flush=True)           
    # pep and cds fasta
    with open(p_pep_final, "wt") as f_pep_final, open(p_cds_final, "wt") as f_cds_final:
        #for ID in ID_selected:
        for ID in pep_ID_list:
            if not ID in ID_selected:
                continue
                
            # pep 
            description, seq_pep = ID_to_info[ID]
            f_pep_final.write(str(description))
            while len(seq_pep) > 60:
                f_pep_final.write(str(seq_pep[:60]) + "\n")
                seq_pep = seq_pep[60:]
            f_pep_final.write(str(seq_pep) + "\n")
            
            # cds
            seq_cds = ID_to_cds[ID]
            f_cds_final.write(str(description))
            while len(seq_cds) > 60:
                f_cds_final.write(str(seq_cds[:60]) + "\n")
                seq_cds = seq_cds[60:]
            f_cds_final.write(str(seq_cds) + "\n")
     
    # gff3
    ID_to_transcript_length = dict()
    with open(p_gff3, "rt") as f_gff3, open(p_gff3_final, "wt") as f_gff3_final:
        longorfs_gff3 = f_gff3.readlines()
        ID_to_block = OrderedDict()
        # associate gff3 blocks with transcript ids
        ID = ""
        for line in longorfs_gff3:
            if ID == "":
                ID = str(line.split("~")[-1].split(";")[0])
                transcript_length = int(line.split("\t")[4])
                ID_to_transcript_length[ID] = transcript_length
                block_list = [line]
            else:
                block_list.append(line)
                if line == "\n":
                    block = "".join(block_list)
                    ID_to_block[ID] = block
                    ID = ""
        # write final gff3
        for ID, block in ID_to_block.items():
            if ID in ID_selected:
                f_gff3_final.write(block)
    
    # bed
    with open(p_bed_final, "wt") as f:
        line = "track name='" + basename + ".TD2.gff3'\n"
        f.write(line)
        for ID in ID_selected:
            description = ID_to_description[ID].rstrip()
            
            # standard bed file format
            chrom = ".".join(ID.split(".")[:-1])
            chromStart = "0"
            chromEnd = str(ID_to_transcript_length[ID])
            name = "ID=" + ID + ";" + ";".join(description.split(" ")[1:])
            score = str(int(ID_to_score[ID] * 1000)) # uses psauron score, range 0-1000 for bed file display
            strand = str(description[-2:-1])
            thickStart = str(description.split(":")[-1].split("-")[0])
            thickEnd = str(description.split(":")[-1].split("-")[1][:-3])
            itemRgb = "0"
            blockCount = "1"
            blockSizes = chromEnd
            blockStarts = "0"
            
            linelist = [chrom, chromStart, chromEnd,
                        name, score, strand,
                        thickStart, thickEnd, itemRgb,
                        blockCount, blockSizes, blockStarts]
            line = "\t".join(linelist)
            f.write(line + "\n")

    # TODO sort bed by default?
    # sort -k 1,1 -k2,2n a.bed
    

    
if __name__ == "__main__":
    main()