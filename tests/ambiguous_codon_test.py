from itertools import product
# import app.get_tables
# TODO: make this an actual test

# IUPAC nucleotide codes including ambiguous bases
iupac_codes = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

# Translation dictionary: codon to amino acid
codon_to_amino_acid = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
    'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
}

# print original codon dict
amino_acid_to_codons = {}
for codon, amino_acid in codon_to_amino_acid.items():
    if amino_acid not in amino_acid_to_codons:
        amino_acid_to_codons[amino_acid] = []
    amino_acid_to_codons[amino_acid].append(codon)

print("ORIGINAL Codons for each Amino Acid:")
for amino_acid, codons in sorted(amino_acid_to_codons.items()):
    print(f"{amino_acid}: {', '.join(sorted(codons))}")
print()

# Generate all possible 3-letter codons using IUPAC codes
iupac_bases = list(iupac_codes.keys())
all_combinations = [''.join(codon) for codon in product(iupac_bases, repeat=3)]
# print(len(all_combinations))

def expand_ambiguous_codon(codon):
    """Expand an ambiguous codon into all possible standard codons it can represent."""
    return [''.join(bases) for bases in product(*[iupac_codes[nuc] for nuc in codon])]

def generate_ambiguous_translation_dict(codon_to_amino_acid):
    """Generate a dictionary of ambiguous codons that map to a single amino acid."""
    ambiguous_translation_dict = {}

    for ambiguous_codon in all_combinations:
        expanded_codons = expand_ambiguous_codon(ambiguous_codon)
        # print(ambiguous_codon, len(expanded_codons))
        
        # Find the amino acid for each expanded codon
        amino_acids = {codon_to_amino_acid[codon] for codon in expanded_codons if codon in codon_to_amino_acid}
        
        # If all expanded codons map to the same amino acid, add to the dictionary
        if len(amino_acids) == 1:
            ambiguous_translation_dict[ambiguous_codon] = amino_acids.pop()
    
    return ambiguous_translation_dict

# Generate the final translation dictionary
ambiguous_translation_dict = generate_ambiguous_translation_dict(codon_to_amino_acid)

# Print the resulting ambiguous codon dictionary
amino_acid_to_codons = {}
for codon, amino_acid in ambiguous_translation_dict.items():
    if amino_acid not in amino_acid_to_codons:
        amino_acid_to_codons[amino_acid] = []
    amino_acid_to_codons[amino_acid].append(codon)

print("NEW Codons for each Amino Acid:")
for amino_acid, codons in sorted(amino_acid_to_codons.items()):
    print(f"{amino_acid}: {', '.join(sorted(codons))}")

'''
Codons for each Amino Acid:
*: TAA, TAG, TAR, TGA, TRA
A: GCA, GCB, GCC, GCD, GCG, GCH, GCK, GCM, GCN, GCR, GCS, GCT, GCV, GCW, GCY
C: TGC, TGT, TGY
D: GAC, GAT, GAY
E: GAA, GAG, GAR
F: TTC, TTT, TTY
G: GGA, GGB, GGC, GGD, GGG, GGH, GGK, GGM, GGN, GGR, GGS, GGT, GGV, GGW, GGY
H: CAC, CAT, CAY
I: ATA, ATC, ATH, ATM, ATT, ATW, ATY
K: AAA, AAG, AAR
L: CTA, CTB, CTC, CTD, CTG, CTH, CTK, CTM, CTN, CTR, CTS, CTT, CTV, CTW, CTY, TTA, TTG, TTR, YTA, YTG, YTR
M: ATG
N: AAC, AAT, AAY
P: CCA, CCB, CCC, CCD, CCG, CCH, CCK, CCM, CCN, CCR, CCS, CCT, CCV, CCW, CCY
Q: CAA, CAG, CAR
R: AGA, AGG, AGR, CGA, CGB, CGC, CGD, CGG, CGH, CGK, CGM, CGN, CGR, CGS, CGT, CGV, CGW, CGY, MGA, MGG, MGR
S: AGC, AGT, AGY, TCA, TCB, TCC, TCD, TCG, TCH, TCK, TCM, TCN, TCR, TCS, TCT, TCV, TCW, TCY
T: ACA, ACB, ACC, ACD, ACG, ACH, ACK, ACM, ACN, ACR, ACS, ACT, ACV, ACW, ACY
V: GTA, GTB, GTC, GTD, GTG, GTH, GTK, GTM, GTN, GTR, GTS, GTT, GTV, GTW, GTY
W: TGG
Y: TAC, TAT, TAY
'''