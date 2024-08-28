import os
import argparse
from TD2.translator import Translator
from pyfaidx import Fasta
import random

def trans_and_print(headers, seqs, translator):
    for header in headers:
        print('███', header)
        name = header[0].split(':')[0]
        start = int(header[0].split(':')[1].split('-')[0])
        end = int(header[0].split(':')[1].split('-')[1].split('(')[0])
        strand = header[0].split('(')[1].split(')')[0]
        if strand == '-':
            sequence = seqs[name][:][end-1:start]
            sequence = sequence.reverse.complement
            sequence = str(sequence)
        else:
            sequence = seqs[name][:][start-1:end]
            sequence = str(sequence)
            
        print(sequence, '█', len(sequence)/3)
        prot_seq, starts, ends = translator.translate(sequence)
        print(prot_seq, '█', len(prot_seq))
        print(len(starts), len(ends))
    

def decode(diff1, diff2, fasta):    
    
    seqs = Fasta(fasta)
    translator = Translator(table=1)
    
    print('decoding orfs in ref but not in query')
    trans_and_print(diff1, seqs, translator)
    print('-'*120)
        
    print('decoding orfs in query but not in ref')
    trans_and_print(diff2, seqs, translator)
    print('-'*120)
    
def investigate(headers1, headers2, limit, fasta):
    
    seqs = Fasta(fasta)
    translator = Translator(table=1)
    
    random.shuffle(headers1)
    random.shuffle(headers2)
    
    print('investigating matched ref orfs')
    trans_and_print(headers1[:limit], seqs, translator)
    print('-'*120)
        
    print('investigating matched query orfs')
    trans_and_print(headers2[:limit], seqs, translator) 
    print('-'*120)

def internal(internal1, internal2, fasta):
    
    seqs = Fasta(fasta)
    translator = Translator(table=1)
    
    print('internal orfs in ref')
    print(internal1, len(internal1))
    for header in internal1:
        print('~~', header)
        name = header[0].split(':')[0]
        start = int(header[0].split(':')[1].split('-')[0])
        end = int(header[0].split(':')[1].split('-')[1].split('(')[0])
        strand = header[0].split('(')[1].split(')')[0]
        if strand == '-':
            sequence = seqs[name][:][end-1:start]
            sequence = sequence.reverse.complement
            sequence = str(sequence)
        else:
            sequence = seqs[name][:][start-1:end]
            sequence = str(sequence)
        original_seq = seqs[name][:]
        if strand == '-':
            original_seq = original_seq.reverse.complement
        print(original_seq, len(original_seq))
        print(sequence, len(sequence))
        print(translator.translate(sequence))
    
    print('-'*120)

    print('internal orfs in query')
    print(internal2, len(internal2))
    for header in internal2:
        print('~~', header)
        name = header[0].split(':')[0]
        start = int(header[0].split(':')[1].split('-')[0])
        end = int(header[0].split(':')[1].split('-')[1].split('(')[0])
        strand = header[0].split('(')[1].split(')')[0]
        if strand == '-':
            sequence = seqs[name][:][end-1:start]
            sequence = sequence.reverse.complement
            sequence = str(sequence)
        else:
            sequence = seqs[name][:][start-1:end]
            sequence = str(sequence)
        original_seq = seqs[name][:]
        if strand == '-':
            original_seq = original_seq.reverse.complement
        print(original_seq, len(original_seq))
        print(sequence, len(sequence))
        print(translator.translate(sequence))

    print('-'*120)

def main(args):

    pep1 = args.pep1
    pep2 = args.pep2
    
    headers1 = []
    headers2 = []

    with open(pep1, 'r') as file:
        for line in file:
            if line.startswith('>'):
                headers1.append((line.strip().split(' ')[-1], line.strip().split(' ')[1]))

    with open(pep2, 'r') as file:
        for line in file:
            if line.startswith('>'):
                headers2.append((line.strip().split(' ')[-1], line.strip().split(' ')[1]))
    
    headers1.sort()
    headers2.sort()
    
    print('-'*120)
    
    print('total orfs in ref:', len(headers1))
    print('number of internal orfs:', len([header for header in headers1 if 'type:internal' in header[1]]))
    print('number of 3\' orfs:', len([header for header in headers1 if 'type:3prime_partial' in header[1]]))
    print('number of 5\' orfs:', len([header for header in headers1 if 'type:5prime_partial' in header[1]]))
    print('number of complete orfs:', len([header for header in headers1 if 'type:complete' in header[1]]))
    print('-'*120)
    
    print('total orfs in query:', len(headers2))
    print('number of internal orfs:', len([header for header in headers2 if 'type:internal' in header[1]]))
    print('number of 3\' orfs:', len([header for header in headers2 if 'type:3prime_partial' in header[1]]))
    print('number of 5\' orfs:', len([header for header in headers2 if 'type:5prime_partial' in header[1]]))
    print('number of complete orfs:', len([header for header in headers2 if 'type:complete' in header[1]]))
    print('-'*120)
        
    diff12 = []
    for header in headers1:
        if header not in headers2:
            diff12.append(header)
    print('orfs in ref but not in query:', len(diff12))
    print(diff12)
    print('-'*120)
            
    diff21 = []
    for header in headers2:
        if header not in headers1:
            diff21.append(header)
    print('orfs in query but not in ref:', len(diff21))
    print(diff21)
    print('-'*120)
    
    if args.diff:
        decode(diff12, diff21, args.f)
    if args.matched:
        investigate(headers1, headers2, 5, args.f)
    if args.internal:
        internal([header for header in headers1 if 'type:internal' in header[1]], [header for header in headers2 if 'type:internal' in header[1]], args.f)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract headers from a FASTA file')
    parser.add_argument('pep1', type=str, help='Input peptide file 1 (reference)')
    parser.add_argument('pep2', type=str, help='Input peptide file 2 (query)')
    parser.add_argument('-f', type=str, help='path to source fasta', default='./tests/data/MANE_test_rna.fna')
    parser.add_argument('--diff', action='store_true', help='investigates problematic orfs in one group but not the other', default=False)
    parser.add_argument('--matched', action='store_true', help='investigates a random sample of the correct matched orfs', default=False)
    parser.add_argument('--internal', action='store_true', help='investigates internal orfs', default=False)
    
    args = parser.parse_args()
    main(args)