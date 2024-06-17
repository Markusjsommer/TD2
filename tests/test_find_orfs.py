from app.translator import Translator

def test_orf_short():
    sequence = "ATGGATTCATGATGTGATCGTATGCTAG"
    translator = Translator()
    protein_sequence, orfs = translator.find_orfs(sequence)

    print("Input DNA Sequence:", sequence)
    print("Translated Protein Sequence:", protein_sequence)
    print("ORFs (Start, End Positions):", orfs)

def test_orf_long():
    sequence = """ctgacctcaggtgatacacctgcctcggcctcccaaagtgctgggattacaggtgtgagccaccatgcctacctaGTTCT\
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
    
if __name__ == "__main__":
    test_orf_short()
    test_orf_long()
    # run with python -m test.test_find_orfs