import os
import sys
import gzip
import time
import argparse
import warnings
from TD2.translator import Translator
from concurrent.futures import ProcessPoolExecutor, as_completed

####################
### INIT GLOBALS ###
####################

complement_map = ('ACTGNactgnYRWSKMDVHBXyrwskmdvhbx',
                  'TGACNtgacnRYWSMKHBDVXrywsmkhbdvx')
complement_table = str.maketrans(complement_map[0], complement_map[1])
ncbi_table_mapping = {
    1: "Universal",
    2: "Mitochondrial-Vertebrates",
    3: "Mitochondrial-Yeast",
    4: "Mitochondrial-Protozoan",
    5: "Mitochondrial-Invertebrates",
    6: "Ciliate, Dasycladacean, Hexamita",
    9: "Mitochondrial-Echinoderm and Mitochondrial-Flatworm",
    10: "Mitochondrial-Ascidian",
    11: "Mitochondrial-Scenedesmus_obliquus",
    12: "Pachysolen_tannophilus",
    13: "Mitochondrial-Chlorophycean",
    14: "SR1_Gracilibacteria",
    15: "Mitochondrial-Thraustochytrium",
    16: "Mitochondrial-Trematode",
    21: "Mitochondrial-Pterobranchia",
    22: "Mitochondrial-Ctenophore",
    23: "Mitochondrial-Mesodinium",
    24: "Mitochondrial-Euplotid",
    25: "Mitochondrial-Peritrich",
    26: "Tetrahymena",
    27: "Candida",
    28: "Acetabularia",
    29: "Mitochondrial-Spirochaete",
    30: "Mitochondrial-Apicomplexa",
    31: "Mitochondrial-Plasmodium",
    33: "SR1_Gracilibacteria"
}

#############
## HELPERS ##
#############

def get_args():
    parser = argparse.ArgumentParser()
    
   # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("-t", dest="transcripts",  type=str, required=True, help="REQUIRED path to transcripts.fasta")
    
    # optional
    parser.add_argument("-O", "--output-dir", dest="output_dir", type=str, required=False, help="path to output results, default=./{transcripts}.TD2_dir", default="./transcripts.TD2_dir")
    parser.add_argument("-m", "--min-length", dest="minimum_length", type=int, required=False, help="minimum protein length, default=100", default=100)
    parser.add_argument("-M", "--absolute-min-length", dest="absolute_min", type=float, required=False, help="absolute minimum protein length for small proteins, default=25", default=25)
    parser.add_argument("-L", "--length-scale", dest="len_scale", type=float, required=False, help="allow short ORFs if they are at least a fraction of the total transcript length, default=0.8 (set to any number greater than 1 to disable this)", default=0.8)
    parser.add_argument("-S", "--strand-specific", dest="strand_specific", action='store_true', required=False, help="set -S for strand-specific ORFs (only analyzes top strand), default=False", default=False)
    parser.add_argument("-G", "--genetic-code", dest="genetic_code", type=int, required=False, help="genetic code a.k.a. translation table, NCBI integer codes, default=1 (universal)", default=1)
    parser.add_argument("-i", "--incomplete-orfs", dest="incomplete_orfs", action='store_true', required=False, help="set -i to yield also consider incomplete ORFs (5' partial, 3' partial, internal), default=False", default=False)
    parser.add_argument("--alt-start", dest="alt_start", action='store_true', required=False, help="include alternative initiator codons from provided table, default=False", default=False)
    parser.add_argument("--top", dest='top', type=int, required=False, help='set -top to also record the top N CDS transcripts by length, default=0', default=0)
    
    # TODO validate gene to transcript mapping file
    parser.add_argument("--gene-trans-map", dest="gene_trans_map", type=str, required=False, help="gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<newline>)")
    
    # TODO make verbosity better
    parser.add_argument("-v", "--verbose", action='store_true', help="set -v for verbose output with progress bars, default=False", default=False)
    
    parser.add_argument("-%", "--memory-threshold", dest="memory_threshold", type=float, required=False, help="percent available virtual memory to set as maximum before flushing to file, default=None", default=None)
    parser.add_argument("-@", "--threads", dest="threads", type=int, required=False, help="number of threads to use, default=1", default=1)
    
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
    print(f"Loaded {len(seq_list)} sequences.")
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

def filter_len(seq_len, orf_len, high, low, scale):
    '''Filters out sequences based on dual length threshold, with a scaling factor between thresholds'''
    if orf_len < low:
        return False
    elif orf_len >= high:
        return True
    else:
        return orf_len >= scale * seq_len / 3 # seq is nucleotides so needs to be divided by 3

def find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only):
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
        filtered_orfs = [orf for orf in orfs if filter_len(len(seq), orf[1] - orf[0], min_len_aa, abs_min_len_aa, len_scale)]
        all_orf_list.append((sequence, filtered_orfs, '+', i+1))
            
    # do reverse strand if not strand-specific
    if not strand_specific:
        for i in range(3):
            # print(reverse_complement(seq)[i:])
            sequence, orfs = translator.find_orfs(reverse_complement(seq)[i:], five_prime_partial=five_prime_partial, three_prime_partial=three_prime_partial)
            # print(sequence, orfs)
            filtered_orfs = [orf for orf in orfs if filter_len(len(seq), orf[1] - orf[0], min_len_aa, abs_min_len_aa, len_scale)]
            all_orf_list.append((sequence, filtered_orfs, '-', i+1))
    
    return all_orf_list

def find_ORFs_with_index(index, seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only):
    '''Finds all open reading frames above minimum length threshold and returns index with result (for multithreading)'''
    orfs = find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only)
    return index, orfs

def calculate_start_end(orf, length, strand, frame):
    '''Calculates start and end positions of ORF in genomic coordinates'''
    start = orf[0] * 3 + frame
    end = orf[1] * 3 + frame - 1
    if strand == '-':
        start, end = length - start + 1, length - end + 1
    return start, end

def get_genetic_code(table_num, alt_start):
    '''Returns the name of the genetic code based on NCBI table number'''
    if alt_start:
        return f'{ncbi_table_mapping[table_num]}_alt'.lower()
    else:
        return ncbi_table_mapping[table_num].lower()
    
def create_gff_block(gene_id, gene_length, prot_length, start, end, strand, count, orf_type, transcript_gene_map=None):
    '''
    gene -> mRNA -> exon -> CDS (5'UTR, 3'UTR)
    gene_id\tTD2\tcomponent\tstart\tend\t.\tstrand\t.\tattributes
    attributes:
        - gene ID=GENE.gene_id~~ORF_ID;Name=ORF_NAME (type, len, description)
        - mRNA ID=ORF_ID;Parent=GENE.gene_id~~ORF_ID;Name=ORF_NAME
        - exon ID=ORF_ID.exon;Parent=ORF_ID
        - CDS ID=cds.ORF_ID;Parent=ORF_ID
        - 5'UTR ID=ORF_ID.utr5p1;Parent=ORF_ID
        - 3'UTR ID=ORF_ID.utr3p1;Parent=ORF_ID
    separate handling based on complete, 5prime_partial, 3prime_partial, and internal
    '''
    if transcript_gene_map:
        gene_acc = transcript_gene_map.get(gene_id, f'GENE.{gene_id}')
    else:
        gene_acc = f'GENE.{gene_id}'
        
    
    
    gene_line = f'{gene_id}\tTD2\tgene\t1\t{gene_length}\t.\t{strand}\t.\tID={gene_acc}~~{gene_id}.p{count};Name={gene_id} type:{orf_type} len:{prot_length} ({strand})'
    mrna_line = f'{gene_id}\tTD2\tmRNA\t1\t{gene_length}\t.\t{strand}\t.\tID={gene_id}.p{count};Parent={gene_acc}~~{gene_id}.p{count};Name={gene_id} type:{orf_type} len:{prot_length} ({strand})'
    exon_line = f'{gene_id}\tTD2\texon\t1\t{gene_length}\t.\t{strand}\t.\tID={gene_id}.p{count}.exon1;Parent={gene_id}.p{count}'

    if orf_type == 'complete':
        if strand == '+':
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count}'
            five_UTR_line = f'{gene_id}\tTD2\tfive_prime_UTR\t1\t{start-1}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr5p1;Parent={gene_id}.p{count}'
            three_UTR_line = f'{gene_id}\tTD2\tthree_prime_UTR\t{end+1}\t{gene_length}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr3p1;Parent={gene_id}.p{count}'
        else:
            start, end = end, start
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count}'
            three_UTR_line = f'{gene_id}\tTD2\tthree_prime_UTR\t1\t{start-1}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr3p1;Parent={gene_id}.p{count}'
            five_UTR_line = f'{gene_id}\tTD2\tfive_prime_UTR\t{end+1}\t{gene_length}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr5p1;Parent={gene_id}.p{count}'
        block = '\n'.join([gene_line, mrna_line, five_UTR_line, exon_line, cds_line, three_UTR_line]) + '\n\n'

    elif orf_type == '5prime_partial':
        if strand == '+':
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count};5_prime_partial=true'
            three_UTR_line = f'{gene_id}\tTD2\tthree_prime_UTR\t{end+1}\t{gene_length}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr3p1;Parent={gene_id}.p{count}'
        else:
            start, end = end, start
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count};5_prime_partial=true'
            three_UTR_line = f'{gene_id}\tTD2\tthree_prime_UTR\t1\t{start-1}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr3p1;Parent={gene_id}.p{count}'
        block = '\n'.join([gene_line, mrna_line, exon_line, cds_line, three_UTR_line]) + '\n\n'

    elif orf_type == '3prime_partial':
        if strand == '+':
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count};3_prime_partial=true'
            five_UTR_line = f'{gene_id}\tTD2\tfive_prime_UTR\t1\t{start-1}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr5p1;Parent={gene_id}.p{count}'
        else:
            start, end = end, start
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count};3_prime_partial=true'
            five_UTR_line = f'{gene_id}\tTD2\tfive_prime_UTR\t{end+1}\t{gene_length}\t.\t{strand}\t.\tID={gene_id}.p{count}.utr5p1;Parent={gene_id}.p{count}'
        block = '\n'.join([gene_line, mrna_line, five_UTR_line, exon_line, cds_line]) + '\n\n'
        
    elif orf_type == 'internal':
        if strand == '+':
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count};5_prime_partial=true;3_prime_partial=true'
        else:
            start, end = end, start
            cds_line = f'{gene_id}\tTD2\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{gene_id}.p{count};Parent={gene_id}.p{count};5_prime_partial=true;3_prime_partial=true'        
        block = '\n'.join([gene_line, mrna_line, exon_line, cds_line]) + '\n\n'

    else:
        raise ValueError(f"Invalid ORF type: {orf_type}")
    
    return block

############
## DRIVER ##
############
    
def main():
    # suppress annoying warnings
    warnings.filterwarnings('ignore')
    print("Python", sys.version, "\n")
    
    print(f"Step 1: Initializing args and loading inputs...", flush=True)
    start_time = time.time()
    
    # parse command line arguments
    args = get_args()
    min_len_aa = args.minimum_length
    abs_min_len_aa = args.absolute_min
    len_scale = args.len_scale
    strand_specific = args.strand_specific
    complete_orfs_only = not args.incomplete_orfs
    genetic_code = args.genetic_code
    alt_start = args.alt_start
    gene_trans_map = args.gene_trans_map
    verbose = args.verbose # TODO: work on this at the end -> tqdm stuff
    threads = args.threads
    memory_threshold = args.memory_threshold
    top = args.top
    
    # create working dir and define output filepaths
    if args.output_dir == "./transcripts.TD2_dir":
        p_transcripts = os.path.abspath(args.transcripts)
        output_dir = os.path.splitext(os.path.basename(p_transcripts))[0]
    else:
        output_dir = os.path.abspath(args.output_dir)
        
    p_pep = os.path.join(output_dir, "longest_orfs.pep")
    p_gff3 = os.path.join(output_dir, "longest_orfs.gff3")
    p_cds = os.path.join(output_dir, "longest_orfs.cds")
    path_list = [p_pep, p_gff3, p_cds]
    
    if top:
        import heapq
        p_cds_top = os.path.join(output_dir, f"longest_orfs.cds.top_{top}_longest")
        path_list.append(p_cds_top)
        longest_cds_heap = []

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Writing to", output_dir, flush=True)
    else:
        if all(os.path.exists(path) for path in path_list):
            print("Output directory already exists. Exiting...", flush=True)
            sys.exit(0)

    # load FASTA
    description_list, seq_list = load_fasta(args.transcripts)
    
    # create translator object
    translator = Translator(table=genetic_code, alt_start=alt_start)
    
    # get transcript to gene mapping
    if gene_trans_map:
        with open(gene_trans_map, 'r') as f:
            transcript_gene_map = {}
            for line in f:
                transcript_id, gene_id = line.strip().split('\t')
                transcript_gene_map[transcript_id] = gene_id
    else:
        transcript_gene_map = None
    
    print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)
    
    #### DYNAMIC MEMORY MODE ####
    if memory_threshold:
        import psutil

        print(f"Step 2+: Finding all ORFs with protein length >= {min_len_aa} in dynamic memory mode with {memory_threshold}% virtual memory allocation...", flush=True)
        start_time = time.time()

        def flush_results(f_pep, f_gff3, f_cds, results):
            for result in results:
                pep_header, prot_seq, cds_header, orf_gene_seq, gff_block = result
                f_pep.write(f'{pep_header}\n{prot_seq}\n')
                f_cds.write(f'{cds_header}\n{orf_gene_seq}\n')
                f_gff3.write(gff_block)
        
        # clear files if they exist
        for path in path_list:
            open(path, 'w').close()

        # write files in append mode in batches
        with open(p_pep, 'a') as f_pep, open(p_gff3, 'a') as f_gff3, open(p_cds, 'a') as f_cds:
            
            gc_name = get_genetic_code(genetic_code, alt_start)
            results = []
            max_batch_size = len(seq_list) // 2
            batch_size = max_batch_size // 5
            batch_start = 0

            while batch_start < len(seq_list):

                # adjust batch size based on memory usage
                cur_memory = psutil.virtual_memory().percent
                distance = cur_memory - memory_threshold
                scale = abs(distance) / memory_threshold
                if distance > 0 and batch_size > 2: # reduce batch size
                    batch_size -= max(1, int(batch_size * scale))
                elif distance < 0 and batch_size < max_batch_size: # increase batch size
                    batch_size += min(int(batch_size * scale), max_batch_size)
                    
                print(f"Percent virtual memory utilization: {cur_memory}%, Batch size: {batch_size}.", end='\t', flush=True)
                
                batch_end = min(batch_start + batch_size, len(seq_list))
                batch_seqs = seq_list[batch_start:batch_end]
                batch_descriptions = description_list[batch_start:batch_end]

                if threads == 1:
                    batch_results = [find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only) for seq in batch_seqs]
                else:
                    batch_results = [None] * len(batch_seqs)
                    with ProcessPoolExecutor(max_workers=threads) as executor:
                        future_to_index = {executor.submit(find_ORFs_with_index, i, seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only): i for i, seq in enumerate(batch_seqs)}
                        for future in as_completed(future_to_index):
                            index, orfs = future.result()
                            batch_results[index] = orfs

                for frames, name, gene_seq in zip(batch_results, batch_descriptions, batch_seqs):
                    count = 1
                    gene_len = len(gene_seq)
                    for frame_info in frames:
                        prot_seq = frame_info[0]
                        orfs = frame_info[1]
                        strand = frame_info[2]
                        frame = frame_info[3]
                        
                        for orf in orfs:
                            start, end = calculate_start_end(orf, gene_len, strand, frame)
                            orf_prot_seq = prot_seq[orf[0]:orf[1]]
                            orf_gene_seq = gene_seq[start-1:end] if strand == '+' else reverse_complement(gene_seq[end-1:start])
                            orf_prot_len = len(orf_prot_seq)
                            orf_type = orf[2]
                            if alt_start and orf_type in ['complete', '3prime_partial'] and orf_prot_seq[0] != 'M':
                                orf_prot_seq = 'M' + orf_prot_seq[1:]

                            pep_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} gc:{gc_name} {name}:{start}-{end}({strand})'
                            cds_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} {name}:{start}-{end}({strand})'
                            gff_block = create_gff_block(name, gene_len, orf_prot_len, start, end, strand, count, orf_type, transcript_gene_map)

                            results.append((pep_header, orf_prot_seq, cds_header, orf_gene_seq, gff_block))

                            if top:
                                cds_length = end - start + 1
                                if len(longest_cds_heap) < top:
                                    heapq.heappush(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
                                else:
                                    heapq.heappushpop(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
                            count += 1
                
                # flush results to final files if memory usage exceeds threshold
                flush_results(f_pep, f_gff3, f_cds, results)
                results.clear()
                print(f"Processed {batch_end} transcripts. {time.time() - start_time:.3f} seconds", flush=True)
                batch_start = batch_end

    #### STANDARD MODE ####
    else:
        print(f"Step 2: Finding ORFs...", flush=True)
        start_time = time.time()
        
        # find all ORFs using single thread
        if threads == 1:
            seq_ORF_list = [find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only) for seq in seq_list]
        # use multithreading/multiprocessing
        else:
            seq_ORF_list = [None] * len(seq_list)  # list to hold results in order
            with ProcessPoolExecutor(max_workers=threads) as executor:
                future_to_index = {executor.submit(find_ORFs_with_index, i, seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only): i for i, seq in enumerate(seq_list)}
                for future in as_completed(future_to_index):
                    index, orfs = future.result()
                    seq_ORF_list[index] = orfs  # place result at the correct index
                    
            # ensure all results are collected
            assert None not in seq_ORF_list
        
        print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)
        
        print(f"Step 3: Writing results to file...", flush=True)
        start_time = time.time() 
        
        # write all peps and cds to file
        with open(p_pep, "wt") as f_pep, open(p_gff3, "wt") as f_gff3, open(p_cds, "wt") as f_cds:
            
            gc_name = get_genetic_code(genetic_code, alt_start)
                
            for frames, name, gene_seq in zip(seq_ORF_list, description_list, seq_list):
                count = 1
                gene_len = len(gene_seq)
                for frame_info in frames:
                    prot_seq = frame_info[0]
                    orfs = frame_info[1]
                    strand = frame_info[2]
                    frame = frame_info[3]
                    
                    for orf in orfs:
                        
                        start, end = calculate_start_end(orf, gene_len, strand, frame)
                        orf_prot_seq = prot_seq[orf[0]:orf[1]]
                        orf_gene_seq = gene_seq[start-1:end] if strand == '+' else reverse_complement(gene_seq[end-1:start])
                        orf_prot_len = len(orf_prot_seq)
                        orf_type = orf[2]
                        
                        # fix protein sequence to start with M if alt_start and complete orf
                        if alt_start and orf_type in ['complete', '3prime_partial'] and orf_prot_seq[0] != 'M':
                            orf_prot_seq = 'M' + orf_prot_seq[1:]

                        # write pep file
                        pep_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} gc:{gc_name} {name}:{start}-{end}({strand})'
                        f_pep.write(f'{pep_header}\n{orf_prot_seq}\n')

                        # write cds file
                        cds_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} {name}:{start}-{end}({strand})'
                        f_cds.write(f'{cds_header}\n{orf_gene_seq}\n')

                        # write gff file
                        f_gff3.write(create_gff_block(name, gene_len, orf_prot_len, start, end, strand, count, orf_type, transcript_gene_map))
                        
                        if top:
                            # keep track of the N(top) longest cds
                            cds_length = end - start + 1
                            if len(longest_cds_heap) < top:
                                heapq.heappush(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
                            else:
                                heapq.heappushpop(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))

                        count += 1
        
    # write longest cds in descending order
    if top:
        with open(p_cds_top, "wt") as f_cds_top:
            for _, cds_header, orf_gene_seq in sorted(longest_cds_heap, reverse=True, key=lambda x: x[0]):
                f_cds_top.write(f'{cds_header}\n{orf_gene_seq}\n')
        
    print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)

if __name__ == "__main__":
    main()