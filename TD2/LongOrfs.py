import os
import sys
import gzip
import time
import argparse
import warnings
from tqdm.auto import tqdm
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
    required.add_argument("-t", dest="transcripts",  type=str, required=True,
                          help="REQUIRED path to transcripts.fasta")
    
    # optional
    parser.add_argument("-O", "--output-dir", dest="output_dir", type=str, required=False,
                        help="path to output results, default=./{transcripts}.TD2_dir",
                        default="./transcripts.TD2_dir")
    parser.add_argument("-m", "--min-length", dest="minimum_length", type=int, required=False,
                        help="minimum protein length for proteins in long transcripts, default=100", default=100)
    parser.add_argument("-M", "--absolute-min-length", dest="absolute_min", type=int, required=False,
                        help="minimum protein length for proteins in short transcripts, default=50",
                        default=50)
    parser.add_argument("-L", "--length-scale", dest="len_scale", type=float, required=False,
                        help="allow short ORFs in short transcripts if the ORF is at least a fraction of the total transcript length, default=0.8",
                        default=0.8)
    parser.add_argument("-S", "--strand-specific", dest="strand_specific", action='store_true', required=False,
                        help="set -S for strand-specific ORFs (only analyzes top strand), default=False", default=False)
    parser.add_argument("-G", "--genetic-code", dest="genetic_code", type=int, required=False,
                        help="genetic code (NCBI integer code), default=1 (universal)", default=1)
    parser.add_argument("-i", "--incomplete-orfs", dest="incomplete_orfs", action='store_true', required=False,
                        help="set -i to consider incomplete ORFs, default=False", default=False)
    parser.add_argument("--alt-start", dest="alt_start", action='store_true', required=False,
                        help="include alternative initiator codons, default=False", default=False)
    parser.add_argument("--top", dest='top', type=int, required=False,
                        help="record the top N CDS transcripts by length, default=0", default=0)
    
    parser.add_argument("--gene-trans-map", dest="gene_trans_map", type=str, required=False,
                        help="gene-to-transcript mapping file (tab-delimited)")
    
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="set -v for verbose output, default=False", default=False)
    
    # Use all available threads by default.
    parser.add_argument("-@", "--threads", dest="threads", type=int, required=False,
                        help=f"number of threads to use, default={os.cpu_count()}",
                        default=os.cpu_count())
    
    parser.add_argument("-%", "--memory-threshold", dest="memory_threshold", type=float, required=False,
                        help="percent of available memory to use per batch, default=None", default=None)
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args

def load_fasta(filepath):
    '''Loads a FASTA file and returns a list of descriptions and sequences'''
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
    """Generator function for reading FASTA/FASTQ files. From Heng Li."""
    last = None
    while True:
        if not last:
            for l in fp:
                if l[0] in '>@':
                    last = l[:-1]
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':
            yield name, ''.join(seqs), None
            if not last: break
        else:
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):
                    last = None
                    yield name, seq, ''.join(seqs)
                    break
            if last:
                yield name, seq, None
                break
    
def reverse_complement(seq):
    '''Reverse complements a DNA sequence.'''
    seq_complement = seq.translate(complement_table)
    return seq_complement[::-1]

def complement(seq):
    '''Complements a DNA sequence.'''
    seq_complement = seq.translate(complement_table)
    return seq_complement

def filter_len(seq_len, orf_len, high, low, scale):
    '''Filters sequences based on dual length thresholds with scaling.'''
    if orf_len < low:
        return False
    elif orf_len >= high:
        return True
    else:
        return orf_len >= scale * seq_len / 3

def find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only):
    '''Finds all ORFs above the minimum length threshold.'''
    all_orf_list = []
    
    if complete_orfs_only:
        five_prime_partial = False
        three_prime_partial = False
    else:
        five_prime_partial = True
        three_prime_partial = True
    
    for i in range(3):
        sequence, orfs = translator.find_orfs(seq[i:], five_prime_partial=five_prime_partial, three_prime_partial=three_prime_partial)
        filtered_orfs = [orf for orf in orfs if filter_len(len(seq), orf[1] - orf[0], min_len_aa, abs_min_len_aa, len_scale)]
        all_orf_list.append((sequence, filtered_orfs, '+', i+1))
            
    if not strand_specific:
        for i in range(3):
            sequence, orfs = translator.find_orfs(reverse_complement(seq)[i:], five_prime_partial=five_prime_partial, three_prime_partial=three_prime_partial)
            filtered_orfs = [orf for orf in orfs if filter_len(len(seq), orf[1] - orf[0], min_len_aa, abs_min_len_aa, len_scale)]
            all_orf_list.append((sequence, filtered_orfs, '-', i+1))
    
    return all_orf_list

def find_ORFs_with_index(index, seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only):
    '''Wrapper to allow multithreading with an index.'''
    orfs = find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only)
    return index, orfs

def calculate_start_end(orf, length, strand, frame):
    '''Calculates start and end positions of an ORF in genomic coordinates.'''
    start = orf[0] * 3 + frame
    end = orf[1] * 3 + frame - 1
    if strand == '-':
        start, end = length - start + 1, length - end + 1
    return start, end

def get_genetic_code(table_num, alt_start):
    '''Returns the genetic code name based on the NCBI table number.'''
    if alt_start:
        return f'{ncbi_table_mapping[table_num]}_alt'.lower()
    else:
        return ncbi_table_mapping[table_num].lower()
    
def create_gff_block(gene_id, gene_length, prot_length, start, end, strand, count, orf_type, transcript_gene_map=None):
    '''Creates a GFF3 block for the given ORF.'''
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
    warnings.filterwarnings('ignore')
    print("Python", sys.version, "\n")
    
    print("Step 1: Initializing args and loading inputs...", flush=True)
    start_time = time.time()
    
    args = get_args()
    min_len_aa = args.minimum_length
    abs_min_len_aa = min(args.absolute_min, min_len_aa)
    len_scale = args.len_scale
    strand_specific = args.strand_specific
    complete_orfs_only = not args.incomplete_orfs
    genetic_code = args.genetic_code
    alt_start = args.alt_start
    gene_trans_map = args.gene_trans_map
    verbose = args.verbose
    threads = args.threads
    memory_threshold = args.memory_threshold
    top = args.top
    
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
    
    description_list, seq_list = load_fasta(args.transcripts)
    
    translator = Translator(table=genetic_code, alt_start=alt_start)
    
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
        print(f"Step 2+: Finding all ORFs with protein length >= {min_len_aa} in dynamic memory mode using {memory_threshold}% of available memory...", flush=True)
        start_time = time.time()
    
        def flush_results(f_pep, f_gff3, f_cds, results):
            for result in results:
                pep_header, prot_seq, cds_header, orf_gene_seq, gff_block = result
                f_pep.write(f'{pep_header}\n{prot_seq}\n')
                f_cds.write(f'{cds_header}\n{orf_gene_seq}\n')
                f_gff3.write(gff_block)
    
        # Clear output files if they exist
        for path in path_list:
            open(path, 'w').close()
    
        with open(p_pep, 'a') as f_pep, open(p_gff3, 'a') as f_gff3, open(p_cds, 'a') as f_cds:
            vm = psutil.virtual_memory()
            allowed_memory = vm.available * (memory_threshold / 100.0)
            avg_result_size = 200  # bytes per transcript result
            batch_size = max(1, int(allowed_memory / avg_result_size))
            print(f"Available memory: {vm.available} bytes, target usage: {allowed_memory:.0f} bytes. Batch size set to: {batch_size}.", flush=True)
    
            batch_start = 0
            while batch_start < len(seq_list):
                batch_end = min(batch_start + batch_size, len(seq_list))
                batch_seqs = seq_list[batch_start:batch_end]
                batch_descriptions = description_list[batch_start:batch_end]
    
                if threads == 1:
                    if verbose:
                        batch_results = [find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only)
                                         for seq in tqdm(batch_seqs, desc="Finding ORFs", leave=False)]
                    else:
                        batch_results = [find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only)
                                         for seq in batch_seqs]
                else:
                    batch_results = [None] * len(batch_seqs)
                    with ProcessPoolExecutor(max_workers=threads) as executor:
                        futures = {executor.submit(find_ORFs_with_index, i, seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only): i
                                   for i, seq in enumerate(batch_seqs)}
                        if verbose:
                            for future in tqdm(as_completed(futures), total=len(futures), desc="Finding ORFs", leave=False):
                                index, orfs = future.result()
                                batch_results[index] = orfs
                        else:
                            for future in as_completed(futures):
                                index, orfs = future.result()
                                batch_results[index] = orfs
    
                results = []
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
    
                            pep_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} gc:{get_genetic_code(genetic_code, alt_start)} {name}:{start}-{end}({strand})'
                            cds_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} {name}:{start}-{end}({strand})'
                            gff_block = create_gff_block(name, gene_len, orf_prot_len, start, end, strand, count, orf_type, transcript_gene_map)
    
                            results.append((pep_header, orf_prot_seq, cds_header, orf_gene_seq, gff_block))
    
                            if top:
                                import heapq
                                cds_length = end - start + 1
                                if len(longest_cds_heap) < top:
                                    heapq.heappush(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
                                else:
                                    heapq.heappushpop(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
                            count += 1
    
                flush_results(f_pep, f_gff3, f_cds, results)
                print(f"Processed {batch_end} transcripts. {time.time() - start_time:.3f} seconds", flush=True)
                batch_start = batch_end
    
    #### STANDARD MODE ####
    else:
        print("Step 2: Finding ORFs...", flush=True)
        start_time = time.time()
        if threads == 1:
            if verbose:
                seq_ORF_list = [find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only)
                                for seq in tqdm(seq_list, desc="Finding ORFs", unit="transcript")]
            else:
                seq_ORF_list = [find_ORFs(seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only)
                                for seq in seq_list]
        else:
            seq_ORF_list = [None] * len(seq_list)
            with ProcessPoolExecutor(max_workers=threads) as executor:
                futures = {executor.submit(find_ORFs_with_index, i, seq, translator, min_len_aa, abs_min_len_aa, len_scale, strand_specific, complete_orfs_only): i 
                           for i, seq in enumerate(seq_list)}
                if verbose:
                    for future in tqdm(as_completed(futures), total=len(futures), desc="Finding ORFs", unit="transcript"):
                        index, orfs = future.result()
                        seq_ORF_list[index] = orfs
                else:
                    for future in as_completed(futures):
                        index, orfs = future.result()
                        seq_ORF_list[index] = orfs
            assert None not in seq_ORF_list
    
        print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)
    
        print("Step 3: Writing results to file...", flush=True)
        start_time = time.time()
        with open(p_pep, "wt") as f_pep, open(p_gff3, "wt") as f_gff3, open(p_cds, "wt") as f_cds:
            if verbose:
                iterator = tqdm(zip(seq_ORF_list, description_list, seq_list), total=len(seq_list),
                                desc="Writing results", unit="transcript")
            else:
                iterator = zip(seq_ORF_list, description_list, seq_list)
    
            for frames, name, gene_seq in iterator:
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
    
                        pep_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} gc:{get_genetic_code(genetic_code, alt_start)} {name}:{start}-{end}({strand})'
                        f_pep.write(f'{pep_header}\n{orf_prot_seq}\n')
    
                        cds_header = f'>{name}.p{count} type:{orf_type} len:{orf_prot_len} {name}:{start}-{end}({strand})'
                        f_cds.write(f'{cds_header}\n{orf_gene_seq}\n')
    
                        f_gff3.write(create_gff_block(name, gene_len, orf_prot_len, start, end, strand, count, orf_type, transcript_gene_map))
                        
                        if top:
                            import heapq
                            cds_length = end - start + 1
                            if len(longest_cds_heap) < top:
                                heapq.heappush(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
                            else:
                                heapq.heappushpop(longest_cds_heap, (cds_length, cds_header, orf_gene_seq))
    
                        count += 1
    
    if top:
        with open(p_cds_top, "wt") as f_cds_top:
            for _, cds_header, orf_gene_seq in sorted(longest_cds_heap, reverse=True, key=lambda x: x[0]):
                f_cds_top.write(f'{cds_header}\n{orf_gene_seq}\n')
    
    print(f"Done. {time.time() - start_time:.3f} seconds", flush=True)

if __name__ == "__main__":
    main()

