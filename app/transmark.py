import os
import sys
import gzip
import copy
import time
import pandas
import pickle
import argparse
import warnings
import numpy as np
from tqdm.auto import tqdm

import pkg_resources

from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq
from app.translator import Translator


def get_args():
    parser = argparse.ArgumentParser()
    
    # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("-t", dest="transcripts",  type=str, required=True, help="REQUIRED path to transcripts.fasta")
    
    # optional
    parser.add_argument("-O", dest="output_dir", type=str, required=False, help="path to output results, default=./", default="./")
    parser.add_argument("-m", dest="minimum_length", type=int, required=False, help="minimum protein length, default=100", default=100)
    parser.add_argument("-S", dest="strand_specific", action='store_true', required=False, help="set -S for strand-specific ORFs (only analyzes top strand), default=False", default=False)
    parser.add_argument("-G", dest="genetic_code", type=int, required=False, help="genetic code a.k.a. translation table, NCBI integer codes, default=1", default=1)
    parser.add_argument("-c", dest="complete_orfs_only", action='store_true', required=False, help="set -c to yield only complete ORFs (peps start with Met (M), end with stop (*)), default=False", default=False)
    parser.add_argument("-R", dest="is_rna", action='store_true', required=False, help="set -R if the transcripts are RNA, default=False", default=False)
    
    parser.add_argument("-@", "--threads", dest="threads", type=int, help="number of threads to use, default=1", default=1)
    # if annotation file provided -> use ORFanage to find ORFs
    parser.add_argument("-A", dest="annotation_file", type=str, required=False, help="path to annotation file transcripts.gff, default=None", default=None)
    parser.add_argument("-M", dest="use_mmseqs2", action='store_true', required=False, help="use mmseqs2 tool, default=False", default=False) # TODO: find better way to describe utility of mmseqs2
    
    parser.add_argument("-v", "--verbose", action='store_true', help="set -v for verbose output with progress bars, default=False", default=False)

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help']) # prints help message if no args are provided by user
    return args
    
def load_fasta(filepath):
    '''Loads a FASTA file and returns a list of descriptions and a list of sequences'''
    print("Loading transcripts at", filepath)
    if filepath[-3:].lower() == ".gz":
        f = gzip.open(filepath, "rt")
    else:
        f = open(filepath, "rt")
    
    description_list = []
    seq_list = []
    for name, seq, qual in readfq(f):
        description_list.append(name)
        seq_list.append(seq.upper())
    
    # convert all str to biopython seqs
    # TODO measure performance hit of this
    # TODO consider scikit-bio for translation
    # TODO alternatively, just implement a dict translation on str (dealing with all possible translation tables will be annoying)
    seq_list = [Seq(s) for s in seq_list]
    
    f.close()
    print("Loaded...\n")
    return description_list, seq_list
    
def readfq(fp): 
    """this is a generator function copied from Heng Li's readfq project https://github.com/lh3/readfq/blob/master/readfq.py"""
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
    
def reverse_complement(seq):
    '''Reverse complements DNA seq, returns N for all non-ATGC chars'''
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq_reverse_complement = ''.join([complement.get(base, 'N') for base in seq[::-1]])
    return seq_reverse_complement

def complement(seq):
    '''Complements DNA seq, returns N for all non-ATGC chars'''
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq_complement = ''.join([complement.get(base, 'N') for base in seq])
    return seq_complement
    
def find_probable_ORFs(seq, min_len_aa, strand_specific, genetic_code):
    '''Finds all open reading frames above minimum length threshold'''
    
    pass
    
def find_complete_ORFs(seq, min_len_aa, strand_specific, genetic_code):
    '''Finds all open reading frames above minimum length threshold'''
    
    # initialize translator
    translator = Translator(table=genetic_code)
    all_orf_list = []
    
    # find orfs in forward frames
    for i in range(3):
        sequence, orfs = translator.find_orfs(seq[i:])
        if len(sequence) < min_len_aa:
            continue
        all_orf_list = (f'+{i+1}', sequence, orfs)
            
        
    # do reverse strand if not strand-specific
    if not strand_specific:
        for i in range(3):
            sequence, orfs = translator.find_orfs(reverse_complement(seq)[i:])
            if len(sequence) < min_len_aa:
                continue
            all_orf_list = (f'-{i+1}', sequence, orfs)
    
    return all_orf_list
    
def main():
    # supress annoying warnings
    warnings.filterwarnings('ignore')
    
    # parse command line arguments
    args = get_args()
    min_len_aa = args.minimum_length
    strand_specific = args.strand_specific
    complete_orfs_only = args.complete_orfs_only
    genetic_code = args.genetic_code
    is_rna = args.is_rna
    
    annotation_file = args.annotation_file
    if annotation_file:
        use_orfanage = True
    else:
        use_orfanage = False
    use_mmseqs2 = args.use_mmseqs2
    
    
    verbose = args.verbose # work on this at the end

    print("Python", sys.version, "\n")
    
    # check args
    # TODO
    
    # OUTLINE -> transcript -> orf -> protein_score -> high_score_orfs
    # orfanage for finding orfs
    # psauron for scoring orfs
    # mmseqs2 for clustering orfs
    
    # load FASTA
    description_list, seq_list = load_fasta(args.transcripts)
    
    # find all possible ORFs
    print(f"Step 1: Finding all ORFs with protein length >= {min_len_aa}", flush=True)
    start_time = time.time()
    all_ORF_list = []
    for seq in seq_list:
        seq_ORF_list = find_all_ORFs(seq, min_len_aa, strand_specific, complete_orfs_only, genetic_code)
    
    print(f"Done. {time.time() - start_time:.2f} seconds", flush=True)

if __name__ == "__main__":
    main()