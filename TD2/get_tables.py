import json
import os, sys
import argparse
import re
from itertools import product


def genetic_codes(input_file, output_dir, verbose=False, three_letter=False):
    """Parses genetic codes from NCBI table and stores each into a separate JSON object (Adapted from https://github.com/linsalrob/genetic_codes/blob/main/pygenetic_code/ncbi_table_to_json.py)"""
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

def add_ambiguous_codons(data, verbose=False):
    """Update the translation table by adding possible ambiguous codons."""
    
    amb2std = ambiguous_to_standard()

    # generate all possible 3-letter codons using IUPAC codes
    iupac_bases = list(amb2std.keys())
    all_combinations = [''.join(codon) for codon in product(iupac_bases, repeat=3)]

    def expand_ambiguous_codon(codon):
        """Expand an ambiguous codon into all possible standard codons it can represent."""
        return [''.join(bases) for bases in product(*[amb2std[nuc] for nuc in codon])]

    def generate_ambiguous_translation_dict(codon_to_amino_acid):
        """Generate a dictionary of ambiguous codons that map to a single amino acid."""
        ambiguous_translation_dict = {}

        for ambiguous_codon in all_combinations:
            expanded_codons = expand_ambiguous_codon(ambiguous_codon)
            
            # find the corresponding amino acid for each expanded codon
            amino_acids = {codon_to_amino_acid[codon] for codon in expanded_codons if codon in codon_to_amino_acid}
            
            # if all expanded codons map to the same amino acid, add to the dictionary
            if len(amino_acids) == 1: 
                ambiguous_translation_dict[ambiguous_codon] = amino_acids.pop()
        
        return ambiguous_translation_dict

    # generate the final translation dictionary
    updated_data = generate_ambiguous_translation_dict(data['codons'])
    data['codons'] = updated_data

    return data

def generate_alternative_initiators(genetic_code_data, output_dir):
    alternative_initiators = { # MANUALLY COMPILED
        1: ['TTG', 'CTG'],
        2: ['ATA', 'ATT', 'ATC', 'GTG'],
        3: ['GTG'],
        4: ['TTA', 'TTG', 'CTG', 'ATT', 'ATA', 'ATC', 'GTG'],
        5: ['TTG', 'ATT', 'ATC', 'GTG', 'ATA'],
        6: [],
        9: [],
        10: [],
        11: ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'GTG'],
        12: ['CTG'],
        13: ['ATA', 'GTG', 'TTG'],
        14: [],
        15: [],
        16: [],
        21: [],
        22: [],
        23: [],
        24: ['TTG', 'CTG'],
        25: [],
        26: ['GTG', 'TTG'],  # These are not part of NCBI table yet
        27: [],
        28: [],
        29: [],
        30: [],
        31: [],
        33: []
    }

    # Process each genetic code table
    for table_number, data in genetic_code_data.items():
        complete_initiators = data.get('initiators', [])
        alt_initiators = alternative_initiators.get(int(table_number), [])
        
        # Check if the table number exists in alternative_initiators
        if int(table_number) not in alternative_initiators:
            print(f"Warning: Table number {table_number} does not exist in alternative codons. Please update the NCBI table versions.")
            continue
        
        # Add alternative initiators to the list if they are not already present
        for alt in alt_initiators:
            if alt not in complete_initiators:
                print(f"Warning: Alternative initiator {alt} does not exist in complete initiators for table {table_number}. Adding...")
                complete_initiators.append(alt)

        # Remove alternative initiators from complete_initiators
        standard_initiators = [initiator for initiator in complete_initiators if initiator not in alt_initiators]

        # Save the updated table
        save_table(output_dir, table_number, {'initiators': sorted(standard_initiators), 'complete_initiators': sorted(complete_initiators), 'codons': data['codons']})
  

def save_table(directory, table_number, data):
    """Save a SINGLE translation table to JSON file."""
    file_name = f"table_{table_number}.json"
    file_path = os.path.join(directory, file_name)
    
    with open(file_path, 'w', encoding='utf-8') as json_file:
        json.dump(data, json_file, indent=2)

def load_tables(directory, tables=None):
    """Loads genetic codes from JSON files in the specified directory."""
    genetic_codes = {}
    
    for file_name in os.listdir(directory):
        if file_name.endswith('.json'):
            try:
                table_number = int(os.path.splitext(file_name)[0].split('_')[-1])
            except ValueError: # skip invalid names
                continue
            if tables and table_number not in tables: # only get specified tables
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

    # print alternative initiators (if they exist)
    if 'complete_initiators' in data:
        print("Alternative Initiator Codons:")
        print(", ".join(sorted(data['complete_initiators'])))
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

    # print codon for each amino acid
    amino_acid_to_codons = {}
    for codon, amino_acid in data['codons'].items():
        if amino_acid not in amino_acid_to_codons:
            amino_acid_to_codons[amino_acid] = []
        amino_acid_to_codons[amino_acid].append(codon)
    
    print("Codons for each Amino Acid:")
    for amino_acid, codons in sorted(amino_acid_to_codons.items()):
        print(f"{amino_acid}: {', '.join(sorted(codons))}")
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
    parser.add_argument('--alt', action='store_true', help="Generate alternative codons.")
    args = parser.parse_args()

    # extract translation tables
    genetic_codes(args.f, args.o, args.v, args.t)

    # generate alternative initiators
    if args.alt:
        generate_alternative_initiators(load_tables(args.o), args.o)

    # pretty print example tables
    example_tables = [1]
    genetic_code_data = load_tables(args.o, example_tables)
    for table_number, data in genetic_code_data.items():
        print_table(table_number, data)

    # add ambiguous codons to all tables
    genetic_code_data = load_tables(args.o)
    print(genetic_code_data)
    for table_number, data in genetic_code_data.items():
        updated_data = add_ambiguous_codons(data) 
        print_table(table_number, updated_data)
        save_table(args.o, table_number, updated_data)
    
    # pretty print modified example tables
    genetic_code_data = load_tables(args.o, example_tables)
    for table_number, data in genetic_code_data.items():
        print_table(table_number, data)