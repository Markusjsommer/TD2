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
    
    # TODO allow for hits from multiple mmseqs databases?
    
    parser.add_argument("--retain_long_orfs_length",  type=int, required=False, help="retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence marks it as coding (default: 1000000, so essentially turned off by default.)", default=1000000)
    parser.add_argument("--single_best_only",  type=str, required=False, help="Retain only the single best orf per transcript (prioritized by homology then orf length)")
    parser.add_argument("-O", dest="output_dir", type=str, required=False, help="same output directory from LongOrfs", default="./transcripts.TD2_dir")
    parser.add_argument("-G", dest="genetic_code", type=int, required=False, help="genetic code a.k.a. translation table, NCBI integer codes, default=1", default=1)
    
    # TODO verbosity
    parser.add_argument("-v", "--verbose", action='store_true', help="set -v for verbose output with progress bars, default=False", default=False)


    # TODO start codon refinement with PWM? Is there anything better?

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help']) # prints help message if no args are provided by user
    return args
    
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
    
    # select transcripts based on psauron score
    df_psauron_selected = df_psauron[df_psauron.apply(lambda row: row['in_frame_score'] > psauron_cutoff and all(row[3:] < row['in_frame_score']), axis=1)]
    print(f"Done.")
    
    # integrate homology search results
    if any([args.retain_mmseqs_hits, args.retain_blastp_hits, args.retain_hmmer_hits]):
        print(f"Step 2: Integrating homology results", flush=True)
    else:
        print(f"Step 2: No homology search results provided, skipping", flush=True)
    
    # parse mmseqs
    if args.retain_mmseqs_hits:
        print(f"Loading MMseqs2 output", flush=True)
        p_mmseqs = args.retain_mmseqs_hits
        df_mmseqs = pandas.read_table(p_mmseqs, header=None)
        hits_mmeseqs = set(df_mmseqs[0])
        print(f"Found {int(len(hits_mmeseqs)):d} ORFs with MMseqs2 hits", flush=True)
    
    # parse blast
    if args.retain_blastp_hits:
        print(f"Loading blastp output", flush=True)
        p_blastp = args.retain_blastp_hits
        df_blastp = pandas.read_table(p_blastp, header=None)
        hits_blastp = set(df_blastp[0])
        print(f"Found {int(len(hits_blastp)):d} ORFs with blastp hits", flush=True)   
        
    # parse hmmer
    if args.retain_hmmer_hits:
        print(f"Loading hmmer output", flush=True)
        p_hmmer = args.retain_hmmer_hits
    
        
    
    
    
    print(f"Done.")
    
    
    # remove internal ORFs
    
    
    # write final outputs to current working directory
    p_pep = os.path.join(output_dir, "longest_orfs.pep")
    p_gff3 = os.path.join(output_dir, "longest_orfs.gff3")
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    
    basename = os.path.basename(args.transcripts)
    p_pep_final = basename + ".TD2.pep"
    p_gff3_final = basename + ".TD2.gff3"
    p_cds_final = basename + ".TD2.cds"
    p_bed_final = basename + ".TD2.bed" # TODO bed file, probably better just to write gff then convert
    
    with open(p_pep_final, "wt") as f_pep_final, open(p_gff3_final, "wt") as f_gff3_final, open(p_cds_final, "wt") as f_cds_final:
        pass
        
    
if __name__ == "__main__":
    main()