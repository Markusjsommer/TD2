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

from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq

from transmark.translator import Translator
from transmark.LongOrfs import load_fasta


def get_args():
    parser = argparse.ArgumentParser()
    
    # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("-t", dest="transcripts",  type=str, required=True, help="REQUIRED path to transcripts.fasta")
    
    # optional
    parser.add_argument("--retain_mmseqs_hits",  type=str, required=False, help="mmseqs output in default '.m8' format. Any ORF with a MMseqs2 match will be retained in the final output.")
    parser.add_argument("--retain_blastp_hits",  type=str, required=False, help="blastp output in '-outfmt 6' format. Any ORF with a blast match will be retained in the final output.")
    parser.add_argument("--retain_pfam_hits",  type=str, required=False, help="domain table output file from running hmmscan to search Pfam. Any ORF with a pfam domain hit will be retained in the final output.")
    
    # TODO mmseqs for domain search. Hayden says she has something better?
    
    parser.add_argument("--retain_long_orfs_length",  type=int, required=False, help="retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence marks it as coding (default: 1000000) so essentially turned off by default.)", default=1000000)
    parser.add_argument("--single_best_only",  type=str, required=False, help="Retain only the single best orf per transcript (prioritized by homology then orf length)")
    parser.add_argument("-O", dest="output_dir", type=str, required=False, help="same output directory from LongOrfs", default="./transcripts.transmark_dir")
    parser.add_argument("-G", dest="genetic_code", type=int, required=False, help="genetic code a.k.a. translation table, NCBI integer codes, default=1", default=1)

    # TODO start codon refinement with PWM? Is there anything better?

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help']) # prints help message if no args are provided by user
    return args
    
def main():
    # supress annoying warnings
    warnings.filterwarnings('ignore')
    
    # parse command line arguments
    args = get_args()
    
    # use absolute path of output
    output_dir = os.path.abspath(args.output_dir)    
    
    # run psauron to score ORFs
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    p_score = os.path.join(output_dir, "psauron_score.csv")
    command_psauron = ["psauron", "-i"+str(p_cds), "-o" +str(p_score), "-m0"]
    result_psauron = subprocess.run(command_psauron, capture_output=False, text=True)
    
    # load psauron results
    df_psauron = pandas.read_csv(p_score, skiprows=3)
    
    # select transcripts based on psauron score
    # in-frame > 0.5
    # all out-of-frame < in-frame
    df_psauron_selected = df_psauron[df_psauron.apply(lambda row: row['in_frame_score'] > 0.5 and all(row[3:] < row['in_frame_score']), axis=1)]
    
    ################### TEMP
    df_psauron_selected.to_csv("~/desktop/foo.csv")
    

    # integrate homology search results
    # parse mmseqs 
    
    # parse blast
    
    # parse pfam
    
    
    
    # write intermediate output files to working directory
    p_scores = os.path.join(output_dir, "longest_orfs.cds.scores")
    with open(p_scores, "wt") as f:
    # TODO
        pass
    p_scores_selected = os.path.join(output_dir, "longest_orfs.cds.scores.selected")
    with open(p_scores_selected, "wt") as f:
    # TODO
        pass
        
    p_best_candidates = os.path.join(output_dir, "longest_orfs.cds.best_candidates.gff3")
    with open(p_best_candidates, "wt") as f:
    # TODO
        pass
    
    
    # write final outputs to current working directory 
    # TODO does TransDecoder really just write to CWD regardless of -O? test this behavior
    with open("transcripts.fasta.transdecoder.pep", "wt") as f:
    # TODO
        pass
    with open("transcripts.fasta.transdecoder.cds", "wt") as f:
    # TODO
        pass
    with open("transcripts.fasta.transdecoder.gff3", "wt") as f:
    # TODO
        pass
    with open("transcripts.fasta.transdecoder.bed", "wt") as f:
    # TODO
        pass
    
    
if __name__ == "__main__":
    main()