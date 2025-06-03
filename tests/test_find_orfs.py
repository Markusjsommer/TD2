from TD2.translator import Translator

def test_orf_short():
    sequence = "ATGGATTCATGATGTGATCGTATGCTAG"
    translator = Translator()
    protein_sequence, orfs = translator.find_orfs(sequence)

    print("Input DNA Sequence:", sequence)
    print("Translated Protein Sequence:", protein_sequence)
    print("ORFs (Start, End Positions):", orfs)

def test_orf_long():
    sequence = """actctgacctcaggtgatacacctgcctcggcctcccaaagtgctgggattacaggtgtgagccaccatgcctacctaGTTCT\
AGCTCTCTTAATtcccacaagagctggttgttaacaAGAGCCTGGCACAAACCCCTCTCTCTCGCCacgtgatctctgca\
catgccagcttcccttccccttctgccatgagtggaaacaGCCTAacgccctcaccagaagcaaatggtggcaccatgct\
tcttgcacaccttcagaactgtgagccaaataaacctctcttctttaaaattattcagcctctggtattcctttataaca\
acacacacacacacacacacacacatacacacacacgcaaaagCAGACTAAAACAGGAACTAATTAGAAATGGTGATGCA\
CCGAGGGATTGGCACCGAGGCTCCCCAACAGGAACTGAGGTCATGGATAGAAGGAcacattcatgttatttttttctaat\
ggttAAGTAATTATTTGCTCTTACTCTCAAAATTTCTGCCAAGGCCTCCCATGGACCAAACTCAACTAGAATCTAGGAAG\
CAGAGAACCTGAGTGTTGCATTCAGCAGAAGTCAGCTTCCTAGGGAATCTTGCAGGAAGGGTGAAGGTAGAGAATCTGGT\
GGGGAAGCAAGCAAATGCCCATCACATGCACTTTCCTCCAACAGAGCGACTCAGATGCTATAAAACTTGCTAACACAGTC\
TCAGGGTCTGATCACAGTAACATACAATCCAGGTTTTAATCATCAGAAATCACAGTCCTATTGTCTTCTGCACAGACCCA\
AACACACTTGGAGGTCATGTTCAATATGAATACCtcacagagaaggaaatttaCACGCGAGAAGTACATCTGCAGAAAGC\
"""
    translator = Translator()
    protein_sequence, orfs = translator.find_orfs(sequence)

    print("Input DNA Sequence:", sequence)
    print("Translated Protein Sequence:", protein_sequence)
    print("ORFs (Start, End Positions):", orfs)
    
def test_orf_long_partials():
    sequence = """actctgacctcaggtgatacacctgcctcggcctcccaaagtgctgggattacaggtgtgagccaccatgcctacctaGTTCT\
AGCTCTCTTAATtcccacaagagctggttgttaacaAGAGCCTGGCACAAACCCCTCTCTCTCGCCacgtgatctctgca\
catgccagcttcccttccccttctgccatgagtggaaacaGCCTAacgccctcaccagaagcaaatggtggcaccatgct\
tcttgcacaccttcagaactgtgagccaaataaacctctcttctttaaaattattcagcctctggtattcctttataaca\
acacacacacacacacacacacacatacacacacacgcaaaagCAGACTAAAACAGGAACTAATTAGAAATGGTGATGCA\
CCGAGGGATTGGCACCGAGGCTCCCCAACAGGAACTGAGGTCATGGATAGAAGGAcacattcatgttatttttttctaat\
ggttAAGTAATTATTTGCTCTTACTCTCAAAATTTCTGCCAAGGCCTCCCATGGACCAAACTCAACTAGAATCTAGGAAG\
CAGAGAACCTGAGTGTTGCATTCAGCAGAAGTCAGCTTCCTAGGGAATCTTGCAGGAAGGGTGAAGGTAGAGAATCTGGT\
GGGGAAGCAAGCAAATGCCCATCACATGCACTTTCCTCCAACAGAGCGACTCAGATGCTATAAAACTTGCTAACACAGTC\
TCAGGGTCTGATCACAGTAACATACAATCCAGGTTTTAATCATCAGAAATCACAGTCCTATTGTCTTCTGCACAGACCCA\
AACACACTTGGAGGTCATGTTCAATATGAATACCtcacagagaaggaaatttaCACGCGAGAAGTACATCTGCAGAAAGC\
"""
    translator = Translator()
    sequence, start_positions, end_positions = translator.translate(sequence)
    protein_sequence, orfs = translator.find_orfs(sequence, five_prime_partials=True, three_prime_partials=True)

    print("Input DNA Sequence:", sequence)
    print("Translated Protein Sequence:", protein_sequence)
    print("ORFs (Start, End Positions):", orfs)
    print("Start Positions:", start_positions)
    print("End Positions:", end_positions)
    
if __name__ == "__main__":
    test_orf_short()
    test_orf_long()
    test_orf_long_partials()
    # run with python -m test.test_find_orfs