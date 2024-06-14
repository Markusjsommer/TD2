import os
import json
import itertools


def ambiguous_to_standard():
    return {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
        'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
        'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }

def generate_all_ambiguous_codons(codons):
    """Generate all possible ambiguous codons for a given set of codons."""
    # Transpose the list of codons to get bases at each position
    bases_at_positions = [set(bases) for bases in zip(*codons)]
    print(bases_at_positions)
    ambiguous_nucleotides = ambiguous_to_standard()

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

    updated_codons = data['codons'].copy()
    
    # Group codons by amino acid
    amino_acid_to_codons = {}
    for codon, amino_acid in data['codons'].items():
        if amino_acid not in amino_acid_to_codons:
            amino_acid_to_codons[amino_acid] = []
        amino_acid_to_codons[amino_acid].append(codon)
    
    # Determine possible ambiguous codons
    for amino_acid, codons in amino_acid_to_codons.items():
        if len(codons) > 1:  # Only consider groups with multiple codons
            print('For amino acid: ', amino_acid)
            ambiguous_codons = generate_all_ambiguous_codons(codons)
            for ambiguous_codon in ambiguous_codons:
                if ambiguous_codon not in updated_codons:
                    updated_codons[ambiguous_codon] = amino_acid

    data['codons'] = updated_codons
    return data


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

# pretty print example tables
example_table = [1]
genetic_code_data = load_tables('./tables', example_table)
for table_number, data in genetic_code_data.items():
    print_table(table_number, data, 8)

# add ambiguous codons to all tables
print('Getting ambiguous codons...')
for table_number, data in genetic_code_data.items():
    updated_data = add_ambiguous_codons(data) 
    print_table(table_number, updated_data, 8)