"""
Parse the NCBI Table that is, by default in app/tables/ncbi_tables.txt and create a json object
Adapted from https://github.com/linsalrob/genetic_codes/blob/main/pygenetic_code/ncbi_table_to_json.py
"""
import json
import os, sys
import argparse
import re
import itertools


def genetic_codes(input_file, output_dir, verbose=False, three_letter=False):
    """Parses genetic codes from NCBI table and stores each into a separate JSON object"""
    # ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    head = re.compile(r'transl_table=(\d+)')
    table = re.compile(r'\w\w\w\s+')

    trans_table = None
    initiators = None
    codons = {}
    
    with open(input_file, 'r', encoding='utf-8') as f:
        for l in f: # skip comments
            if l.startswith('#'):
                continue
            l = l.strip()
            if not l: # skip empty lines
                continue
            m = head.search(l) # check if the line is a header
            if m:
                if codons:
                    # save the current table as a JSON file
                    save_table(output_dir, trans_table, {'initiators': sorted(list(initiators)), 'codons': codons})

                # start a new translation table
                trans_table = int(m.groups()[0])
                initiators = set()
                codons = {}
            elif table.match(l):
                p = l.strip().split()
                i = 0
                while i < len(p):
                    codon = p[i]
                    i += 2 if three_letter else 1
                    if codon in codons:
                        if verbose:
                            print(
                                f"For {codon} in trans table {trans_table} we already had {codons[codon]} but now "
                                f"we found {p[i]}", file=sys.stderr)
                    codons[codon] = p[i]
                    i += 1 if three_letter else 2
                    if i < len(p) and p[i] == 'i':
                        initiators.add(codon)
                        i += 1
            else:
                print(f"ERROR: We don't know what {l} means?", file=sys.stderr)

        # save the last table
        if codons:
            save_table(output_dir, trans_table, {'initiators': sorted(list(initiators)), 'codons': codons})

def generate_all_ambiguous_codons(codons):
    """Generate all possible ambiguous codons for a given set of codons."""
    # Transpose the list of codons to get bases at each position
    bases_at_positions = [set(bases) for bases in zip(*codons)]
    print(bases_at_positions)
    
    # Generate all possible ambiguous codons
    ambiguous_combinations = [[]]
    for bases in bases_at_positions:
        new_combinations = []
        for combination in ambiguous_combinations:
            # Find ambiguous nucleotides that can represent the current set of bases
            for ambiguous, possible_bases in ambiguous_nucleotides.items():
                if isinstance(possible_bases, list) and bases <= set(possible_bases):
                    new_combinations.append(combination + [ambiguous])
        ambiguous_combinations = new_combinations
    
    # Convert each combination back into a codon string
    ambiguous_codons = [''.join(codon) for codon in ambiguous_combinations]
    return ambiguous_codons

def add_ambiguous_codons(data, verbose=False):
    """Update the translation table by adding possible ambiguous codons."""

    ambiguous_nucleotides = ambiguous_to_standard()

    updated_codons = data['codons'].copy()
    
    # group codons by amino acid
    amino_acid_to_codons = {}
    for codon, amino_acid in data['codons'].items():
        if amino_acid not in amino_acid_to_codons:
            amino_acid_to_codons[amino_acid] = []
        amino_acid_to_codons[amino_acid].append(codon)
    
    # determine possible ambiguous codons
    for amino_acid, codons in amino_acid_to_codons.items():
        for i in range(3):  # each position in the codon
            # collect nucleotides at the same position in all codons
            position_nucleotides = [set(c[i] for c in codons)]
            possible_ambiguous = list(itertools.product(*position_nucleotides))
            
            for ambiguous_codon in possible_ambiguous:
                ambiguous_codon = ''.join(ambiguous_codon)
                possible_codon = True
                
                for pos, char in enumerate(ambiguous_codon):
                    if char not in ambiguous_nucleotides:
                        possible_codon = False
                        break
                    valid_bases = ambiguous_nucleotides[char] if isinstance(ambiguous_nucleotides[char], list) else [ambiguous_nucleotides[char]]
                    if any(c[pos] not in valid_bases for c in codons):
                        possible_codon = False
                        break
                
                if possible_codon:
                    updated_codons[ambiguous_codon] = amino_acid

    data['codons'] = updated_codons
    return data

def save_table(directory, table_number, data):
    """Save a SINGLE translation table to JSON file."""
    file_name = f"table_{table_number}.json"
    file_path = os.path.join(directory, file_name)
    
    with open(file_path, 'w', encoding='utf-8') as json_file:
        json.dump(data, json_file, indent=2)

def load_tables(directory, tables=[]):
    """Loads genetic codes from JSON files in the specified directory."""
    genetic_codes = {}
    
    for file_name in os.listdir(directory):
        if file_name.endswith('.json'):
            try:
                table_number = int(os.path.splitext(file_name)[0].split('_')[-1])
            except ValueError: # skip invalid names
                continue
            if table_number not in tables: # only get specified tables
                continue

            file_path = os.path.join(directory, file_name)
            
            # extract the table number from the file name
            table_number = os.path.splitext(file_name)[0].split('_')[-1]
            
            # load the JSON data
            with open(file_path, 'r', encoding='utf-8') as json_file:
                data = json.load(json_file)
                genetic_codes[table_number] = data
                
    # sort the dictionary by table numbers (converted to integers for proper sorting)
    return dict(sorted(genetic_codes.items(), key=lambda x: int(x[0])))

def print_table(table_number, data, columns=4):
    """Pretty prints the translation table in a grid format."""
    print(f"Translation Table {table_number}:")
    print("=" * (17 * columns - 1))
    
    # print initiator codons
    print("Initiator Codons:")
    print(", ".join(sorted(data['initiators'])))
    print()
    
    # print codon to amino acid mapping in grid format
    print("Codon to Amino Acid Mapping:")
    codons = sorted(data['codons'].items())
    
    # determine the number of rows
    num_rows = (len(codons) + columns - 1) // columns
    
    # print the header for the grid
    header = "{:<10} {:<5} " * columns
    print(header.format(*(f"Codon", f"AA") * columns))
    print("-" * (17 * columns - 1))
    
    # print each row of the grid
    for row in range(num_rows):
        row_data = []
        for col in range(columns):
            index = row + num_rows * col
            if index < len(codons):
                codon, amino_acid = codons[index]
                row_data.extend([codon, amino_acid])
            else:
                row_data.extend(["", ""])  # fill empty cells if not enough data
        print(header.format(*row_data))
    
    print("=" * (17 * columns - 1))
    print()

def ambiguous_to_standard():
    return {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
        'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
        'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse genetic codes from NCBI table and store into JSON objects.")
    parser.add_argument('-f', type=str, help="Path to the input file containing genetic codes.", default='./tables/ncbi_tables.txt')
    parser.add_argument('-o', type=str, help="Directory where the JSON files will be saved.", default='./tables')
    parser.add_argument('-v', action='store_true', help="Enable verbose output.")
    parser.add_argument('-t', action='store_true', help="Use three-letter amino acid codes.")
    args = parser.parse_args()

    # extract translation tables
    genetic_codes(args.f, args.o, args.v, args.t)

    # pretty print example tables
    example_tables = [1, 11, 21, 31]
    genetic_code_data = load_tables(args.o, example_tables)
    for table_number, data in genetic_code_data.items():
        print_table(table_number, data)

    # add ambiguous codons to all tables
    genetic_code_data = load_tables(args.o)
    for table_number, data in genetic_code_data.items():
        updated_data = add_ambiguous_codons(data) 
        save_table(args.o, table_number, updated_data)
    
    # pretty print modified example tables
    genetic_code_data = load_tables(args.o, example_tables)
    for table_number, data in genetic_code_data.items():
        print_table(table_number, data)
       
