import pytest
from TD2.translator import (
    Translator,
    standardize_sequence,
    legal_letters,
    legal_tables,
    stopless_tables,
    one_to_three_letter,
    load_translation_table,
)


def test_translate_short():
    t = Translator()
    assert t.translate("ATG") == ("M", [0], [])


def test_translate_returns_three_lists():
    # Triple: protein string, start positions (codons), stop positions (codons).
    protein, starts, stops = Translator().translate("ATGAAATAA")
    assert protein == "MK*"
    assert starts == [0]
    assert stops == [2]


def test_translate_rna_input_treated_as_dna():
    # Translator should accept U and convert to T.
    assert Translator().translate("AUGAAA") == ("MK", [0], [])


def test_translate_lowercase_is_normalized():
    assert Translator().translate("atgaaataa") == ("MK*", [0], [2])


def test_translate_truncates_partial_codon():
    # Trailing 1-2 bases that don't form a full codon are dropped.
    protein, _, _ = Translator().translate("ATGAAAT")
    assert protein == "MK"


def test_translate_unknown_codon_is_X():
    # NNN doesn't disambiguate to a single amino acid, so the table maps it to X.
    protein, _, _ = Translator().translate("ATGNNNTAA")
    assert protein == "MX*"


def test_translate_three_letter_codes():
    t = Translator(three_letter=True)
    protein, _, _ = t.translate("ATGAAATAA")
    assert protein == "MetLysTer"


def test_translate_illegal_letter_raises():
    with pytest.raises(ValueError):
        Translator().translate("ATGZZZAAA")


def test_translate_alt_start_table1():
    # Table 1 alt initiators include CTG and TTG. Default initiators do not.
    default = Translator(table=1)
    alt = Translator(table=1, alt_start=True)
    # CTG is leucine in both translations, but only counts as a start under alt_start.
    _, default_starts, _ = default.translate("CTGAAATAA")
    _, alt_starts, _ = alt.translate("CTGAAATAA")
    assert default_starts == []
    assert alt_starts == [0]


def test_translate_three_frames_all_returned():
    seq = "ATGAAATAA" + "G"  # 10 nt, all three frames have content
    translations = Translator().translate_three_frames(seq)
    assert [x[0] for x in translations] == ["+1", "+2", "+3"]
    # Frame +1 is the in-frame ORF.
    assert translations[0][1] == "MK*"


def test_table_2_mitochondrial_AGA_is_stop():
    # Vertebrate mitochondrial code: AGA/AGG are stops, ATA initiates.
    t = Translator(table=2)
    protein, _, stops = t.translate("ATGAAAAGA")
    assert protein.endswith("*")
    assert stops == [2]


def test_table_2_ATA_is_alt_initiator():
    # Table 2 has ATG only by default; ATA/ATC/ATT/GTG come in under alt_start.
    default_starts = Translator(table=2).translate("ATAAAATAA")[1]
    alt_starts = Translator(table=2, alt_start=True).translate("ATAAAATAA")[1]
    assert default_starts == []
    assert 0 in alt_starts


def test_invalid_table_raises():
    with pytest.raises(ValueError):
        Translator(table=999)


def test_stopless_tables_warn(capsys):
    # Tables 27/28/31 are documented as stopless; constructing one should print a warning.
    Translator(table=27)
    out = capsys.readouterr().out
    assert "does not contain stop codons" in out


def test_legal_tables_matches_paper():
    # Paper text: "26 current NCBI translation tables".
    assert len(legal_tables()) == 26


def test_one_to_three_letter_covers_all_aas():
    table = one_to_three_letter()
    # 20 standard AAs + stop (*) + unknown (X).
    standard = set("ACDEFGHIKLMNPQRSTVWY")
    assert standard.issubset(table.keys())
    assert table["*"] == "Ter"
    assert table["X"] == "Unk"


def test_standardize_sequence_truncation_and_uppercase():
    assert standardize_sequence("atgaaat") == "ATGAAA"
    assert standardize_sequence("AU") == ""  # below codon length
    assert standardize_sequence("AUG") == "ATG"


def test_legal_letters_set_includes_iupac():
    letters = legal_letters()
    for n in "ACGTNRYSWKMBDHV":
        assert n in letters


@pytest.mark.parametrize("table_num", sorted(legal_tables()))
def test_every_listed_table_loads(table_num):
    # Exercises the importlib.resources.files() resource path that replaced pkg_resources.
    codons, initiators = load_translation_table(table_num)
    assert len(codons) >= 64  # ambiguous expansions add many more
    assert all(len(c) == 3 for c in codons)
    assert all(len(c) == 3 for c in initiators)
    assert codons["ATG"] in {"M", "I"}  # standard or special table 2 (ATG initiator semantics)
