"""Sanity-check the GFF3 emission paths in LongOrfs.create_gff_block."""
import pytest

from TD2.LongOrfs import create_gff_block, get_genetic_code


def _parse(block: str):
    return [line.split("\t") for line in block.strip().splitlines()]


def test_complete_orf_plus_strand_has_gene_mrna_exon_cds_and_utrs():
    block = create_gff_block(
        gene_id="tx1",
        gene_length=300,
        prot_length=80,
        start=10,
        end=249,
        strand="+",
        count=1,
        orf_type="complete",
    )
    rows = _parse(block)
    feature_types = [row[2] for row in rows]
    assert feature_types == ["gene", "mRNA", "five_prime_UTR", "exon", "CDS", "three_prime_UTR"]
    cds = next(row for row in rows if row[2] == "CDS")
    assert cds[3] == "10" and cds[4] == "249" and cds[6] == "+"


def test_complete_orf_minus_strand_swaps_start_and_end():
    block = create_gff_block(
        gene_id="tx1",
        gene_length=300,
        prot_length=80,
        start=249,
        end=10,
        strand="-",
        count=1,
        orf_type="complete",
    )
    rows = _parse(block)
    cds = next(row for row in rows if row[2] == "CDS")
    # On the minus strand the caller passes start>end; the GFF block must re-order
    # the coordinates so that CDS start < CDS end (GFF spec).
    assert int(cds[3]) < int(cds[4])
    assert cds[6] == "-"


def test_partial_orfs_emit_5_and_3_prime_partial_flags():
    block_5p = create_gff_block("tx", 300, 50, start=1, end=150, strand="+", count=1, orf_type="5prime_partial")
    block_3p = create_gff_block("tx", 300, 50, start=10, end=300, strand="+", count=1, orf_type="3prime_partial")
    assert "5_prime_partial=true" in block_5p
    assert "3_prime_partial=true" in block_3p


def test_internal_orf_marks_both_partial_flags():
    block = create_gff_block("tx", 300, 100, start=1, end=300, strand="+", count=1, orf_type="internal")
    assert "5_prime_partial=true" in block
    assert "3_prime_partial=true" in block
    # No UTRs for internal ORFs by spec.
    assert "five_prime_UTR" not in block
    assert "three_prime_UTR" not in block


def test_invalid_orf_type_raises():
    with pytest.raises(ValueError):
        create_gff_block("tx", 300, 100, 1, 100, "+", 1, "bogus")


def test_get_genetic_code_uses_alt_suffix():
    assert get_genetic_code(1, False) == "the standard code"
    assert get_genetic_code(1, True) == "the standard code_alt"
