"""
Parse the NCBI Table that is, by default in codes/ncbi_tables.txt and create a json object
Adapted from https://github.com/linsalrob/genetic_codes/blob/main/pygenetic_code/ncbi_table_to_json.py
"""
import json
import os, sys
import argparse
import re
# from app.translator import standard_to_ambiguous, ambiguous_to_standard

def genetic_codes(args):
    """Parses genetic codes from NCBI table and stores each into a separate JSON object"""
    input_file = args.f
    verbose = args.v
    output_dir = args.o  

    # Ensure the output directory exists
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
                    file_name = f"table_{trans_table}.json"
                    with open(os.path.join(output_dir, file_name), 'w', encoding='utf-8') as json_file:
                        json.dump({'initiators': list(initiators), 'codons': codons}, json_file, indent=2)

                # start a new translation table
                trans_table = int(m.groups()[0])
                initiators = set()
                codons = {}
            elif table.match(l):
                p = l.strip().split()
                i = 0
                while i < len(p):
                    codon = p[i]
                    i += 2
                    if codon in codons:
                        if verbose:
                            print(
                                f"For {codon} in trans table {trans_table} we already had {codons[codon]} but now "
                                f"we found {p[i]}", file=sys.stderr)
                    codons[codon] = p[i]
                    i += 1
                    if i < len(p) and p[i] == 'i':
                        initiators.add(codon)
                        i += 1
            else:
                print(f"ERROR: We don't know what {l} means?", file=sys.stderr)

        # save the last table
        if codons:
            file_name = f"table_{trans_table}.json"
            with open(os.path.join(output_dir, file_name), 'w', encoding='utf-8') as json_file:
                json.dump({'initiators': list(initiators), 'codons': codons}, json_file, indent=2)

def add_ambiguous_codons():
    pass

def load_tables(directory, tables=None):
    """Loads genetic codes from JSON files in the specified directory."""
    genetic_codes = {}
    
    for file_name in os.listdir(directory):
        if file_name.endswith('.json'):
            if tables:
                try:
                    table_number = int(os.path.splitext(file_name)[0].split('_')[-1])
                except ValueError: # invalid name
                    continue
                if table_number not in tables: # only get specified tables
                    continue

            file_path = os.path.join(directory, file_name)
            
            # Extract the table number from the file name
            table_number = os.path.splitext(file_name)[0].split('_')[-1]
            
            # Load the JSON data
            with open(file_path, 'r', encoding='utf-8') as json_file:
                data = json.load(json_file)
                genetic_codes[table_number] = data
                
    return genetic_codes

def print_table(table_number, data, columns=4):
    """Pretty prints the translation table in a grid format."""
    print(f"Translation Table {table_number}:")
    print("=" * (17 * columns - 1))
    
    # Print initiator codons
    print("Initiator Codons:")
    print(", ".join(sorted(data['initiators'])))
    print()
    
    # Print codon to amino acid mapping in grid format
    print("Codon to Amino Acid Mapping:")
    codons = sorted(data['codons'].items())
    
    # Determine the number of rows
    num_rows = (len(codons) + columns - 1) // columns
    
    # Print the header for the grid
    header = "{:<10} {:<5} " * columns
    print(header.format(*(f"Codon", f"AA") * columns))
    print("-" * (17 * columns - 1))
    
    # Print each row of the grid
    for row in range(num_rows):
        row_data = []
        for col in range(columns):
            index = row + num_rows * col
            if index < len(codons):
                codon, amino_acid = codons[index]
                row_data.extend([codon, amino_acid])
            else:
                row_data.extend(["", ""])  # Fill empty cells if not enough data
        print(header.format(*row_data))
    
    print("=" * (17 * columns - 1))
    print()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse genetic codes from NCBI table and store into JSON objects.")
    parser.add_argument('-f', type=str, help="Path to the input file containing genetic codes.", default='./tables/ncbi_tables.txt')
    parser.add_argument('-o', type=str, help="Directory where the JSON files will be saved.", default='./tables')
    parser.add_argument('-v', action='store_true', help="Enable verbose output.")
    args = parser.parse_args()

    genetic_codes(args)

    example_tables = [1, 11, 21, 31]
    genetic_code_data = load_tables(args.o, example_tables)
    for table_number in sorted(genetic_code_data.keys(), key=int):
        data = genetic_code_data[table_number]
        print_table(table_number, data)

    # add_ambiguous_codons()

    genetic_code_data = load_tables(args.o, example_tables)
    for table_number in sorted(genetic_code_data.keys(), key=int):
        data = genetic_code_data[table_number]
        print_table(table_number, data)

