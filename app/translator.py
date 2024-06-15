
import json
import argparse

class Translator:
    def __init__(self, table=1, rna=False, three_letter=False):
        if table not in legal_tables():
            raise ValueError(f"Table {table} is not a legal table")
        if table in stopless_tables():
            print(f"[WARNING] Table {table} does not contain stop codons")
        self.table = table
        self.rna = rna
        self.three_letter = three_letter
    
    def translate(self, sequence):
        '''Translate the given DNA sequence to protein sequence'''
        translation_dict, initiators = load_translation_table(self.table)
        dna_sequence = standardize_sequence(sequence, self.rna)
        if set(dna_sequence) - legal_letters(): 
            raise ValueError(f"Sequence {dna_sequence} contains illegal letters")

        protein_sequence = []
        initiator_positions = []

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            amino_acid = translation_dict.get(codon, 'X')
            if self.three_letter:
                amino_acid = one_to_three_letter().get(amino_acid, 'Unk')
            protein_sequence.append(amino_acid)

            if codon in initiators:
                initiator_positions.append(i // 3)

        protein_string = ''.join(protein_sequence)

        return protein_string, initiator_positions

def standardize_sequence(sequence, rna=False):
    '''Ensures that given sequence is valid DNA, upper case, and multiple of 3'''
    dna = sequence.upper().strip()
    if rna:
        dna = dna.replace('U', 'T')
    if len(dna) % 3 != 0:
        dna += 'N' * (3 - len(dna) % 3)
    return dna

def load_translation_table(table_num=1):
    '''Get the corresponding translation table for the genetic code'''
    with open(f'./tables/table_{table_num}.json', 'r') as file:
        data = json.load(file)
    return data['codons'], data['initiators']

def legal_letters():
    return {
        "A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"
    }

def legal_tables():
    return {
        1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33
    }

def stopless_tables():
    return {
        27, 28, 31
    }

def one_to_three_letter():
    return {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp",
        "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
        "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
        "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
        "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
        "*": "Ter", "X": "Unk"
    }

def main():
    parser = argparse.ArgumentParser(description='DNA to Protein Translator')
    parser.add_argument('sequence', type=str, help='DNA sequence to translate')
    parser.add_argument('-t', '--table', type=int, default=1, help='Translation table to use')
    parser.add_argument('-r', '--rna', action='store_true', help='Use RNA instead of DNA')
    parser.add_argument('-3', '--three-letter', action='store_true', help='Use three-letter amino acid codes')
    args = parser.parse_args()

    # Create an instance of the Translator class
    dna2prot = Translator(args.table, args.rna, args.three_letter)

    # Call the translate method
    translation, starts = dna2prot.translate(args.sequence)

    # Print the translated text
    print(translation, starts)

if __name__ == '__main__':
    main()