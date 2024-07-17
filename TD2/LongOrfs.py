import os
import sys
import gzip
import time
import argparse
import warnings
from TD2.translator import Translator

### INIT GLOBALS ###
complement_map = ('ACTGNactgnYRWSKMDVHBXyrwskmdvhbx',
                  'TGACNtgacnRYWSMKHBDVXrywsmkhbdvx')
complement_table = str.maketrans(complement_map[0], complement_map[1])


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
    
    parser.add_argument("-@", "--threads", dest="threads", type=int, help="number of threads to use, default=1", default=1)
    
    # TODO gene to transcript mapping file
    
    # TODO if annotation file provided -> use ORFanage to find ORFs
    parser.add_argument("-A", dest="annotation_file", type=str, required=False, help="path to annotation file transcripts.gff, default=None", default=None)
    
    parser.add_argument("-v", "--verbose", action='store_true', help="set -v for verbose output with progress bars, default=False", default=False)
    # TODO refactor m-start to alt-start and reverse logic, also test/copy default behavior of transdecoder with regard to alternate starts (alternate starts technically result in m in the protein)
    parser.add_argument("--m-start", dest="m_start", action='store_true', required=False, help="use only ATG as initiator codon, default=False", default=False)

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
        seq_list.append(seq)
    
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
    seq_complement = seq.translate(complement_table)
    return seq_complement[::-1]

def complement(seq):
    '''Complements DNA seq, returns N for all non-ATGC chars'''
    seq_complement = seq.translate(complement_table)
    return seq_complement
    
def find_ORFs(seq, translator, min_len_aa, strand_specific, complete_orfs_only):
    '''Finds all open reading frames above minimum length threshold'''
    
    all_orf_list = []
    
    # determine whether to allow partial ORFs
    if complete_orfs_only:
        five_prime_partial = False
        three_prime_partial = False
    else:
        five_prime_partial = True
        three_prime_partial = True
    
    # find orfs in forward frames
    for i in range(3):
        # print(seq[i:])
        sequence, orfs = translator.find_orfs(seq[i:], five_prime_partial=five_prime_partial, three_prime_partial=three_prime_partial)
        # print(sequence, orfs)
        filtered_orfs = [orf for orf in orfs if orf[1] - orf[0] >= min_len_aa]
        all_orf_list.append((sequence, filtered_orfs, '+', i+1))
            
    # do reverse strand if not strand-specific
    if not strand_specific:
        for i in range(3):
            # print(reverse_complement(seq)[i:])
            sequence, orfs = translator.find_orfs(reverse_complement(seq)[i:], five_prime_partial=five_prime_partial, three_prime_partial=three_prime_partial)
            # print(sequence, orfs)
            filtered_orfs = [orf for orf in orfs if orf[1] - orf[0] >= min_len_aa]
            all_orf_list.append((sequence, filtered_orfs, '-', i+1))
    
    return all_orf_list

def calculate_start_end(orf, length, strand, frame):
    '''Calculates start and end positions of ORF in genomic coordinates'''
    start = orf[0] * 3 + frame
    end = orf[1] * 3 + frame - 1
    if strand == '-':
        start, end = length - start + 1, length - end + 1
    return start, end
    
def main():
    # suppress annoying warnings
    warnings.filterwarnings('ignore')
    print("Python", sys.version, "\n")
    
    print(f"Initializing args...", flush=True)
    start_time = time.time()
    
    # parse command line arguments
    args = get_args()
    min_len_aa = args.minimum_length
    strand_specific = args.strand_specific
    complete_orfs_only = args.complete_orfs_only
    genetic_code = args.genetic_code
    m_start = args.m_start
    
    # use absolute path of output
    output_dir = os.path.abspath(args.output_dir)
    
    annotation_file = args.annotation_file
    if annotation_file:
        use_orfanage = True
    else:
        use_orfanage = False
    
    verbose = args.verbose # TODO: work on this at the end -> tqdm stuff

    # TODO: check args
    
    # create working directory
    working_base = "transcripts.transmark_dir"
    working_dir = os.path.join(output_dir, working_base)
    print("Writing to", working_dir, flush=True)
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    # define output filepaths
    p_pep = os.path.join(working_dir, "longest_orfs.pep")
    p_gff3 = os.path.join(working_dir, "longest_orfs.gff3")
    p_cds = os.path.join(working_dir, "longest_orfs.cds")
    p_cds_top500 = os.path.join(working_dir, "longest_orfs.cds.top_500_longest")
    
    print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)
    
    
    print(f"Step 1: Finding all ORFs with protein length >= {min_len_aa}", flush=True)
    start_time = time.time()
    
    # load FASTA
    description_list, seq_list = load_fasta(args.transcripts)
    
    # create translator object
    translator = Translator(table=genetic_code, m_start=m_start)
        
    # find all ORFs
    seq_ORF_list = [find_ORFs(seq, translator, min_len_aa, strand_specific, complete_orfs_only) for seq in seq_list]
     
    
    print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)
    
    print(f"Step 2: Writing results to file", flush=True)
    start_time = time.time() 
    
    with open(p_pep, "wt") as f:
        
        if genetic_code == 1 and m_start:
            gc_name = 'universal'
        else:
            gc_name = f'ncbi_table_{genetic_code}'
            
        for entries, desc, gene_seq in zip(seq_ORF_list, description_list, seq_list):
            count = 1
            for entry in entries:
                prot_seq = entry[0]
                orfs = entry[1]
                strand = entry[2]
                frame = entry[3]
                name = desc
                for orf in orfs:
                    start, end = calculate_start_end(orf, len(gene_seq), strand, frame)
                    orf_seq = prot_seq[orf[0]:orf[1]]
                    header = f'>{name}.p{count} type:{orf[2]} gc:{gc_name} {name}:{start}-{end}({strand})'
                    f.write(f'{header}\n{orf_seq}\n')
                    count += 1
        
    with open(p_gff3, "wt") as f:
        # TODO
        pass
    with open(p_cds, "wt") as f:
        # TODO
        pass
    with open(p_cds_top500, "wt") as f:
        # TODO
        pass
    
    print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)

if __name__ == "__main__":
    main()