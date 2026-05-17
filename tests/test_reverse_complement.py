import pytest
from TD2.LongOrfs import reverse_complement, complement


def test_reverse_complement_default():
    assert reverse_complement("GATACA") == "TGTATC"


def test_reverse_complement_N_padding():
    assert reverse_complement("NNNNNNA") == "TNNNNNN"


def test_reverse_complement_preserves_case():
    # Case is preserved per the complement_table mapping.
    assert reverse_complement("aTGc") == "gCAt"


def test_reverse_complement_palindrome():
    # GAATTC (EcoRI site) is its own reverse complement.
    assert reverse_complement("GAATTC") == "GAATTC"


def test_reverse_complement_empty():
    assert reverse_complement("") == ""


@pytest.mark.parametrize("base,expected", [
    ("R", "Y"),  # A/G <-> T/C
    ("Y", "R"),
    ("S", "S"),  # G/C self-complementary
    ("W", "W"),  # A/T self-complementary
    ("K", "M"),  # G/T <-> A/C
    ("M", "K"),
    ("B", "V"),  # not-A <-> not-T
    ("V", "B"),
    ("D", "H"),  # not-C <-> not-G
    ("H", "D"),
])
def test_complement_handles_iupac_ambiguity(base, expected):
    assert complement(base) == expected


def test_double_reverse_complement_is_identity():
    seq = "ACGTRYSWKMBDHVN"
    assert reverse_complement(reverse_complement(seq)) == seq
