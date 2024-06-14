
import json
import argparse

class Translator:
    def __init__(self, sequence, table=1, rna=False):
        if table not in legal_tables():
            raise ValueError(f"Table {table} is not a legal table")
        self.table = table
        self.rna = rna
        self.sequence = standardize_sequence(sequence)
        if set(self.sequence) - legal_letters(): 
            raise ValueError(f"Sequence {sequence} contains illegal letters")
    
    def translate(self, text):
        # Implement the translation logic here
        pass

    def standardize_sequence(sequence):
        '''Ensures that given sequence is valid DNA, upper case, and multiple of 3'''
        dna = sequence.upper()
        if self.rna:
            dna.replace('U', 'T')
        if len(dna) % 3 != 0:
            dna += 'N' * (3 - len(dna) % 3)
        return dna

    def legal_letters():
        return {
            "A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"
        }

    def legal_tables():
        return {
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33
        }

    def load_translation_table(table_num=1):
        '''Get the corresponding translation table for the genetic code'''
        return json.loads(open(f'{table_num}.json').read())

    def one_to_three_letter():
        return {
            "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp",
            "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
            "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
            "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
            "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
            "*": "Ter"
        }

# TODO: translate the sequence -> string, and then mark where start codons are (using initiators) -> list[int]

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='DNA to Protein Translator')
    parser.add_argument('sequence', type=str, help='DNA sequence to translate')
    parser.add_argument('-t', '--table', type=int, default=1, help='Translation table to use')
    parser.add_argument('-r', '--rna', action='store_true', help='Use RNA instead of DNA')
    args = parser.parse_args()

    # Create an instance of the Translator class
    translator = Translator()

    # Call the translate method
    translation = translator.translate(args.text)

    # Print the translated text
    print(translation)

if __name__ == '__main__':
    main()