import os
import argparse

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
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract headers from a FASTA file')
    parser.add_argument('pep1', type=str, help='Input peptide file 1 (reference)')
    parser.add_argument('pep2', type=str, help='Input peptide file 2 (query)')
    
    args = parser.parse_args()
    main(args)