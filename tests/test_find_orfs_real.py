"""Exercises LongOrfs.find_ORFs on a real human mRNA transcript."""
from TD2.LongOrfs import find_ORFs, Translator, calculate_start_end, filter_len


NM_000014_6 = (
    # 2.5 kb tail of NM_000014.6 (A2M) -- enough to contain ORFs in multiple frames.
    "GGGACCAGATGGATTGTAGGGAGTAGGGTACAATACAGTCTGTTCTCCTCCAGCTCCTTC"
    "TTTCTGCAACATGGGGAAGAACAAACTCCTTCATCCAAGTCTGGTTCTTCTCCTCTTGGT"
    "CCTCCTGCCCACAGACGCCTCAGTCTCTGGAAAACCGCAGTATATGGTTCTGGTCCCCTC"
    "CCTGCTCCACACTGAGACCACTGAGAAGGGCTGTGTCCTTCTGAGCTACCTGAATGAGAC"
    "AGTGACTGTAAGTGCTTCCTTGGAGTCTGTCAGGGGAAACAGGAGCCTCTTCACTGACCT"
    "GGAGGCGGAGAATGACGTACTCCACTGTGTCGCCTTCGCTGTCCCAAAGTCTTCATCCAA"
    "TGAGGAGGTAATGTTCCTCACTGTCCAAGTGAAAGGACCAACCCAAGAATTTAAGAAGCG"
    "GACCACAGTGATGGTTAAGAACGAGGACAGTCTGGTCTTTGTCCAGACAGACAAATCAAT"
    "CTACAAACCAGGGCAGACAGTGAAATTTCGTGTTGTCTCCATGGATGAAAACTTTCACCC"
    "CCTGAATGAGTTGATTCCACTAGTATACATTCAGGATCCCAAAGGAAATCGCATCGCACA"
    "ATGGCAGAGTTTCCAGTTAGAGGGTGGCCTCAAGCAATTTTCTTTTCCCCTCTCATCAGA"
    "GCCCTTCCAGGGCTCCTACAAGGTGGTGGTACAGAAGAAATCAGGTGGAAGGACAGAGCA"
    "CCCTTTCACCGTGGAGGAATTTGTTCTTCCCAAGTTTGAAGTACAAGTAACAGTGCCAAA"
    "GATAATCACCATCTTGGAAGAAGAGATGAATGTATCAGTGTGTGGCCTATACACATATGG"
    "GAAGCCTGTCCCTGGACATGTGACTGTGAGCATTTGCAGAAAGTATAGTGACGCTTCCGA"
    "CTGCCACGGTGAAGATTCACAGGCTTTCTGTGAGAAATTCAGTGGACAGCTAAACAGCCA"
    "TGGCTGCTTCTATCAGCAAGTAAAAACCAAGGTCTTCCAGCTGAAGAGGAAGGAGTATGA"
)


def test_find_ORFs_returns_six_frames_with_complete_a2m_start():
    translator = Translator(table=1)
    frame_orfs = find_ORFs(
        NM_000014_6,
        translator,
        min_len_aa=100,
        abs_min_len_aa=100,
        len_scale=1.1,
        strand_specific=False,
        complete_orfs_only=False,
        all_stopless=False,
    )

    # Always six frames: 3 on the + strand, 3 on the - strand.
    assert len(frame_orfs) == 6
    strands = [f[2] for f in frame_orfs]
    assert strands.count("+") == 3 and strands.count("-") == 3

    # The canonical A2M ORF (starting with MGKNKLLHPSLVLLLLVLLPTDASVSGKPQYMVLVPSLLHTET...)
    # lives on the + strand in this slice; recover it from one of the frames.
    matches = []
    for prot_seq, orfs, strand, frame in frame_orfs:
        if strand != "+":
            continue
        for orf in orfs:
            sub = prot_seq[orf[0]:orf[1]].rstrip("*")
            if sub.startswith("MGKNKLLHPSLVLLLLVLLPTDA"):
                matches.append((sub, strand, frame, orf[2]))
    assert matches, "expected the canonical A2M signal-peptide ORF"


def test_strand_specific_skips_minus_strand():
    translator = Translator(table=1)
    frame_orfs = find_ORFs(
        NM_000014_6,
        translator,
        min_len_aa=100,
        abs_min_len_aa=100,
        len_scale=1.1,
        strand_specific=True,
        complete_orfs_only=False,
        all_stopless=False,
    )
    assert len(frame_orfs) == 3
    assert all(f[2] == "+" for f in frame_orfs)


def test_filter_len_uses_high_low_and_scale():
    # orf_len < low: rejected.
    assert filter_len(seq_len=1000, orf_len=50, high=100, low=80, scale=1.1) is False
    # orf_len >= high: accepted unconditionally.
    assert filter_len(seq_len=1000, orf_len=200, high=100, low=80, scale=1.1) is True
    # In between: accept iff orf_len >= scale * seq_len / 3.
    # 1000 * 1.1 / 3 ~= 366 codons, so a 90-codon ORF in a 1000 nt transcript fails the scale gate.
    assert filter_len(seq_len=1000, orf_len=90, high=100, low=80, scale=1.1) is False
    # ... but in a short transcript the same length passes the scale gate.
    assert filter_len(seq_len=200, orf_len=90, high=100, low=80, scale=1.1) is True


def test_calculate_start_end_plus_strand():
    # orf (0-based codons) [10, 20) in frame 1 (1-based), + strand, transcript length 1000.
    start, end = calculate_start_end((10, 20, "complete"), length=1000, strand="+", frame=1)
    assert start == 31  # 10 * 3 + 1
    assert end == 60    # 20 * 3 + 1 - 1


def test_calculate_start_end_minus_strand_flips_coordinates():
    start, end = calculate_start_end((10, 20, "complete"), length=1000, strand="-", frame=1)
    # On the minus strand, coords are reflected over the transcript length.
    assert start == 1000 - 31 + 1
    assert end == 1000 - 60 + 1
    assert start > end
