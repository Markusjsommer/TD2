import json
import pkg_resources

class Translator:
    def __init__(self, table=1, rna=False, three_letter=False, m_start=False):
        if table not in legal_tables():
            raise ValueError(f"Table {table} is not a legal table")
        if table in stopless_tables():
            print(f"[WARNING] Table {table} does not contain stop codons")
        self.table_num = table
        self.table, self.initiators = load_translation_table(table)
        self.rna = rna
        self.three_letter = three_letter
        if m_start: # sets M as the sole initiator
            self.initiators = ['ATG']
    
    def translate(self, sequence):
        '''
        Translate the given sequence to protein sequence
        Parameters: sequence (str): DNA/RNA sequence to translate
        Returns: Tuple[str, List[int], List[int]]: Translated protein sequence, (sorted) initiator positions, and (sorted) stop positions
        '''
        dna_sequence = standardize_sequence(sequence, self.rna)
        if set(dna_sequence) - legal_letters(): 
            raise ValueError(f"Sequence {dna_sequence} contains illegal letters")

        protein_sequence = []
        start_positions = []
        end_positions = []

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            amino_acid = self.table.get(codon, 'X')
            if codon in self.initiators:
                start_positions.append(i // 3)
            elif amino_acid == '*':
                end_positions.append(i // 3)
            if self.three_letter:
                amino_acid = one_to_three_letter().get(amino_acid, 'Unk')
            protein_sequence.append(amino_acid)

        protein_string = ''.join(protein_sequence)

        return protein_string, start_positions, end_positions
    
    def translate_three_frames(self, sequence, strand='+'):
        '''
        Translate the given sequence in all three frames of given strand
        Parameters: sequence (str): DNA/RNA sequence to translate
        Returns: List[Tuple[str, List[int]]]: List of translated protein sequences and their initiator positions
        '''
        translations = []
        for frame in range(3):
            translated_sequence, initiator_positions, end_positions = self.translate(sequence[frame:])
            translations.append((f'{strand}{frame+1}', translated_sequence, initiator_positions, end_positions))
        return translations
    
    def find_orfs(self, sequence, five_prime_partials=False, three_prime_partials=False):
        '''
        Find all open reading frames in the given sequence
        Parameters: sequence (str): DNA/RNA sequence to analyze
        Returns: str, List[Tuple[int, int]]: Translated protein sequence and list of ORFs (start, end) positions
        '''
        protein_sequence, start_positions, end_positions = self.translate(sequence)
        orfs = []
        cur_pos = -1
        start_index = 0
        end_index = 0
        
        if five_prime_partials:
            end_positions = end_positions + [len(protein_sequence)]
        if three_prime_partials:
            start_positions = [0] + start_positions
               
        while start_index < len(start_positions) and end_index < len(end_positions):
            if start_positions[start_index] <= cur_pos:
                start_index += 1
            elif end_positions[end_index] <= start_positions[start_index]:
                end_index += 1
            else:
                orfs.append((start_positions[start_index], end_positions[end_index]))
                cur_pos = end_positions[end_index]
                start_index += 1
                end_index += 1
            
        return protein_sequence, orfs


def standardize_sequence(sequence, rna=False):
    '''Ensures that given sequence is valid DNA, upper case, and multiple of 3'''
    dna = sequence.upper().strip()
    dna = dna[:-(len(dna) % 3)] # truncates to last full codon
    if rna:
        dna = dna.replace('U', 'T')        
    return dna

def load_translation_table(table_num=1):
    '''Get the corresponding translation table for the genetic code'''
    path_table = pkg_resources.resource_filename(__name__, f'tables/table_{table_num}.json')
    with open(path_table, 'r') as file:
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
    import argparse

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