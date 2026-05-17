from TD2.translator import Translator


def test_complete_orf_returned():
    # AAA ATG AAA AAA TAA AAA -> K M K K * K
    sequence = "AAA" "ATG" "AAA" "AAA" "TAA" "AAA"
    protein, orfs = Translator().find_orfs(sequence)
    assert protein == "KMKK*K"
    assert orfs == [(1, 5, "complete")]


def test_no_orf_without_partials_when_no_start():
    sequence = "AAA" "AAA" "TAA" "AAA"  # no ATG
    _, orfs = Translator().find_orfs(sequence)
    assert orfs == []


def test_five_prime_partial_emitted():
    # First stop precedes any ATG -> 5' partial spanning to that first stop.
    sequence = "AAA" "AAA" "TAA" "AAA" "ATG" "AAA" "TAA"
    protein, orfs = Translator().find_orfs(
        sequence, five_prime_partial=True, three_prime_partial=False
    )
    assert protein == "KK*KMK*"
    assert (0, 3, "5prime_partial") in orfs
    assert (4, 7, "complete") in orfs


def test_three_prime_partial_emitted():
    # ATG without a downstream stop -> 3' partial to end of protein.
    sequence = "AAA" "ATG" "AAA" "AAA"  # protein KMKK, no stop
    protein, orfs = Translator().find_orfs(
        sequence, five_prime_partial=False, three_prime_partial=True
    )
    assert protein == "KMKK"
    assert orfs == [(1, 4, "3prime_partial")]


def test_internal_orf_when_no_starts_and_no_stops():
    # No start, no stop, both partials requested -> single internal ORF.
    sequence = "AAA" "AAA" "AAA"
    protein, orfs = Translator().find_orfs(
        sequence, five_prime_partial=True, three_prime_partial=True
    )
    assert protein == "KKK"
    assert orfs == [(0, 3, "internal")]


def test_internal_when_start_after_protein_start_no_stop():
    # Start exists but no stop, both partials requested -> internal spanning everything.
    sequence = "AAA" "ATG" "AAA"  # K M K, no stop
    _, orfs = Translator().find_orfs(
        sequence, five_prime_partial=True, three_prime_partial=True
    )
    assert (0, 3, "internal") in orfs


def test_all_stopless_treats_every_non_stop_codon_as_start():
    sequence = "AAA" "AAA" "TAA" "AAA"  # K K * K
    _, orfs = Translator().find_orfs(sequence, all_stopless=True)
    # Every codon before the stop is a potential start; the function should
    # emit a single ORF using the smallest "start" not after the stop.
    assert any(orf[2] == "complete" and orf[1] == 3 for orf in orfs)


def test_overlapping_complete_orfs_share_no_stop():
    # Two non-overlapping complete ORFs in one frame.
    sequence = ("ATG" "AAA" "TAA"   # ORF 1: M K *
                "AAA"                # spacer K
                "ATG" "GGG" "TAA")   # ORF 2: M G *
    _, orfs = Translator().find_orfs(sequence)
    completes = [o for o in orfs if o[2] == "complete"]
    assert len(completes) == 2
    # Both ORFs end at their own stop codon.
    ends = {o[1] for o in completes}
    assert ends == {3, 7}
