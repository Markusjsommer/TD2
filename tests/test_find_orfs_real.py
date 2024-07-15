from transmark.LongOrfs import find_ORFs, Translator, calculate_start_end

def test_orf_NM_000014_6():
    sequence = """GGGACCAGATGGATTGTAGGGAGTAGGGTACAATACAGTCTGTTCTCCTCCAGCTCCTTC\
TTTCTGCAACATGGGGAAGAACAAACTCCTTCATCCAAGTCTGGTTCTTCTCCTCTTGGT\
CCTCCTGCCCACAGACGCCTCAGTCTCTGGAAAACCGCAGTATATGGTTCTGGTCCCCTC\
CCTGCTCCACACTGAGACCACTGAGAAGGGCTGTGTCCTTCTGAGCTACCTGAATGAGAC\
AGTGACTGTAAGTGCTTCCTTGGAGTCTGTCAGGGGAAACAGGAGCCTCTTCACTGACCT\
GGAGGCGGAGAATGACGTACTCCACTGTGTCGCCTTCGCTGTCCCAAAGTCTTCATCCAA\
TGAGGAGGTAATGTTCCTCACTGTCCAAGTGAAAGGACCAACCCAAGAATTTAAGAAGCG\
GACCACAGTGATGGTTAAGAACGAGGACAGTCTGGTCTTTGTCCAGACAGACAAATCAAT\
CTACAAACCAGGGCAGACAGTGAAATTTCGTGTTGTCTCCATGGATGAAAACTTTCACCC\
CCTGAATGAGTTGATTCCACTAGTATACATTCAGGATCCCAAAGGAAATCGCATCGCACA\
ATGGCAGAGTTTCCAGTTAGAGGGTGGCCTCAAGCAATTTTCTTTTCCCCTCTCATCAGA\
GCCCTTCCAGGGCTCCTACAAGGTGGTGGTACAGAAGAAATCAGGTGGAAGGACAGAGCA\
CCCTTTCACCGTGGAGGAATTTGTTCTTCCCAAGTTTGAAGTACAAGTAACAGTGCCAAA\
GATAATCACCATCTTGGAAGAAGAGATGAATGTATCAGTGTGTGGCCTATACACATATGG\
GAAGCCTGTCCCTGGACATGTGACTGTGAGCATTTGCAGAAAGTATAGTGACGCTTCCGA\
CTGCCACGGTGAAGATTCACAGGCTTTCTGTGAGAAATTCAGTGGACAGCTAAACAGCCA\
TGGCTGCTTCTATCAGCAAGTAAAAACCAAGGTCTTCCAGCTGAAGAGGAAGGAGTATGA\
AATGAAACTTCACACTGAGGCCCAGATCCAAGAAGAAGGAACAGTGGTGGAATTGACTGG\
AAGGCAGTCCAGTGAAATCACAAGAACCATAACCAAACTCTCATTTGTGAAAGTGGACTC\
ACACTTTCGACAGGGAATTCCCTTCTTTGGGCAGGTGCGCCTAGTAGATGGGAAAGGCGT\
CCCTATACCAAATAAAGTCATATTCATCAGAGGAAATGAAGCAAACTATTACTCCAATGC\
TACCACGGATGAGCATGGCCTTGTACAGTTCTCTATCAACACCACCAATGTTATGGGTAC\
CTCTCTTACTGTTAGGGTCAATTACAAGGATCGTAGTCCCTGTTACGGCTACCAGTGGGT\
GTCAGAAGAACACGAAGAGGCACATCACACTGCTTATCTTGTGTTCTCCCCAAGCAAGAG\
CTTTGTCCACCTTGAGCCCATGTCTCATGAACTACCCTGTGGCCATACTCAGACAGTCCA\
GGCACATTATATTCTGAATGGAGGCACCCTGCTGGGGCTGAAGAAGCTCTCCTTCTATTA\
TCTGATAATGGCAAAGGGAGGCATTGTCCGAACTGGGACTCATGGACTGCTTGTGAAGCA\
GGAAGACATGAAGGGCCATTTTTCCATCTCAATCCCTGTGAAGTCAGACATTGCTCCTGT\
CGCTCGGTTGCTCATCTATGCTGTTTTACCTACCGGGGACGTGATTGGGGATTCTGCAAA\
ATATGATGTTGAAAATTGTCTGGCCAACAAGGTGGATTTGAGCTTCAGCCCATCACAAAG\
TCTCCCAGCCTCACACGCCCACCTGCGAGTCACAGCGGCTCCTCAGTCCGTCTGCGCCCT\
CCGTGCTGTGGACCAAAGCGTGCTGCTCATGAAGCCTGATGCTGAGCTCTCGGCGTCCTC\
GGTTTACAACCTGCTACCAGAAAAGGACCTCACTGGCTTCCCTGGGCCTTTGAATGACCA\
GGACAATGAAGACTGCATCAATCGTCATAATGTCTATATTAATGGAATCACATATACTCC\
AGTATCAAGTACAAATGAAAAGGATATGTACAGCTTCCTAGAGGACATGGGCTTAAAGGC\
ATTCACCAACTCAAAGATTCGTAAACCCAAAATGTGTCCACAGCTTCAACAGTATGAAAT\
GCATGGACCTGAAGGTCTACGTGTAGGTTTTTATGAGTCAGATGTAATGGGAAGAGGCCA\
TGCACGCCTGGTGCATGTTGAAGAGCCTCACACGGAGACCGTACGAAAGTACTTCCCTGA\
GACATGGATCTGGGATTTGGTGGTGGTAAACTCAGCAGGTGTGGCTGAGGTAGGAGTAAC\
AGTCCCTGACACCATCACCGAGTGGAAGGCAGGGGCCTTCTGCCTGTCTGAAGATGCTGG\
ACTTGGTATCTCTTCCACTGCCTCTCTCCGAGCCTTCCAGCCCTTCTTTGTGGAGCTCAC\
AATGCCTTACTCTGTGATTCGTGGAGAGGCCTTCACACTCAAGGCCACGGTCCTAAACTA\
CCTTCCCAAATGCATCCGGGTCAGTGTGCAGCTGGAAGCCTCTCCCGCCTTCCTAGCTGT\
CCCAGTGGAGAAGGAACAAGCGCCTCACTGCATCTGTGCAAACGGGCGGCAAACTGTGTC\
CTGGGCAGTAACCCCAAAGTCATTAGGAAATGTGAATTTCACTGTGAGCGCAGAGGCACT\
AGAGTCTCAAGAGCTGTGTGGGACTGAGGTGCCTTCAGTTCCTGAACACGGAAGGAAAGA\
CACAGTCATCAAGCCTCTGTTGGTTGAACCTGAAGGACTAGAGAAGGAAACAACATTCAA\
CTCCCTACTTTGTCCATCAGGTGGTGAGGTTTCTGAAGAATTATCCCTGAAACTGCCACC\
AAATGTGGTAGAAGAATCTGCCCGAGCTTCTGTCTCAGTTTTGGGAGACATATTAGGCTC\
TGCCATGCAAAACACACAAAATCTTCTCCAGATGCCCTATGGCTGTGGAGAGCAGAATAT\
GGTCCTCTTTGCTCCTAACATCTATGTACTGGATTATCTAAATGAAACACAGCAGCTTAC\
TCCAGAGATCAAGTCCAAGGCCATTGGCTATCTCAACACTGGTTACCAGAGACAGTTGAA\
CTACAAACACTATGATGGCTCCTACAGCACCTTTGGGGAGCGATATGGCAGGAACCAGGG\
CAACACCTGGCTCACAGCCTTTGTTCTGAAGACTTTTGCCCAAGCTCGAGCCTACATCTT\
CATCGATGAAGCACACATTACCCAAGCCCTCATATGGCTCTCCCAGAGGCAGAAGGACAA\
TGGCTGTTTCAGGAGCTCTGGGTCACTGCTCAACAATGCCATAAAGGGAGGAGTAGAAGA\
TGAAGTGACCCTCTCCGCCTATATCACCATCGCCCTTCTGGAGATTCCTCTCACAGTCAC\
TCACCCTGTTGTCCGCAATGCCCTGTTTTGCCTGGAGTCAGCCTGGAAGACAGCACAAGA\
AGGGGACCATGGCAGCCATGTATATACCAAAGCACTGCTGGCCTATGCTTTTGCCCTGGC\
AGGTAACCAGGACAAGAGGAAGGAAGTACTCAAGTCACTTAATGAGGAAGCTGTGAAGAA\
AGACAACTCTGTCCATTGGGAGCGCCCTCAGAAACCCAAGGCACCAGTGGGGCATTTTTA\
CGAACCCCAGGCTCCCTCTGCTGAGGTGGAGATGACATCCTATGTGCTCCTCGCTTATCT\
CACGGCCCAGCCAGCCCCAACCTCGGAGGACCTGACCTCTGCAACCAACATCGTGAAGTG\
GATCACGAAGCAGCAGAATGCCCAGGGCGGTTTCTCCTCCACCCAGGACACAGTGGTGGC\
TCTCCATGCTCTGTCCAAATATGGAGCAGCCACATTTACCAGGACTGGGAAGGCTGCACA\
GGTGACTATCCAGTCTTCAGGGACATTTTCCAGCAAATTCCAAGTGGACAACAACAACCG\
CCTGTTACTGCAGCAGGTCTCATTGCCAGAGCTGCCTGGGGAATACAGCATGAAAGTGAC\
AGGAGAAGGATGTGTCTACCTCCAGACATCCTTGAAATACAATATTCTCCCAGAAAAGGA\
AGAGTTCCCCTTTGCTTTAGGAGTGCAGACTCTGCCTCAAACTTGTGATGAACCCAAAGC\
CCACACCAGCTTCCAAATCTCCCTAAGTGTCAGTTACACAGGGAGCCGCTCTGCCTCCAA\
CATGGCGATCGTTGATGTGAAGATGGTCTCTGGCTTCATTCCCCTGAAGCCAACAGTGAA\
AATGCTTGAAAGATCTAACCATGTGAGCCGGACAGAAGTCAGCAGCAACCATGTCTTGAT\
TTACCTTGATAAGGTGTCAAATCAGACACTGAGCTTGTTCTTCACGGTTCTGCAAGATGT\
CCCAGTAAGAGATCTGAAACCAGCCATAGTGAAAGTCTATGATTACTACGAGACGGATGA\
GTTTGCAATTGCTGAGTACAATGCTCCTTGCAGCAAAGATCTTGGAAATGCTTGAAGACC\
ACAAGGCTGAAAAGTGCTTTGCTGGAGTCCTGTTCTCAGAGCTCCACAGAAGACACGTGT\
TTTTGTATCTTTAAAGACTTGATGAATAAACACTTTTTCTGGTCAATGTC\
"""

    complete_orfs = False
    translator = Translator(table=1, m_start=True)
    frame_orfs = find_ORFs(sequence, translator, min_len_aa=100, strand_specific=False, complete_orfs_only=complete_orfs)


    for entry in frame_orfs:
        prot_seq = entry[0]
        orfs = entry[1]
        strand = entry[2]
        frame = entry[3]
        name = 'NM_000014.6'
        print(name, orfs, len(sequence), strand, frame)
        count = 1
        for orf in orfs:
            start, end = calculate_start_end(orf, len(sequence), strand, frame)
            orf_seq = prot_seq[orf[0]:orf[1]]
            header = f'>{name}.p{count} type:{orf[2]} gc:universal {name}:{start}-{end}({strand})'
            # print(header)
            # print(orf_seq)
            count += 1

    
if __name__ == "__main__":
    test_orf_NM_000014_6()
    # run with python -m test.test_find_orfs_real