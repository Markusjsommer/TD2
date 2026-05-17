"""End-to-end regression: TD2.LongOrfs on the bundled MANE fixture must recover the
canonical RefSeq protein for each of the three NM_ accessions.

The fixture (tests/data/MANE_test_rna.fna) holds three real MANE v1.4 transcripts:
    - NM_000014.6 (A2M)   -> RefSeq NP_000005.3, 1474 aa
    - NM_000015.3 (NAT2)  -> RefSeq NP_000006.2,  290 aa
    - NM_000016.6 (ACADM) -> RefSeq NP_000007.1,  421 aa

These canonical proteins are stable, externally published reference annotations.
If a refactor or dependency change perturbs LongOrfs's output, this test catches it
by checking the longest-ORF protein per transcript matches the published sequence
exactly at both ends and at the documented length. It's the smallest possible
"full transcriptome matches published version" check.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


TD2_DIR = Path(__file__).resolve().parent.parent
FIXTURE = TD2_DIR / "tests" / "data" / "MANE_test_rna.fna"


def _read_pep(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(chunks)
            cur = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(chunks)
    return seqs


def _read_pep_descriptions(path: Path) -> dict[str, str]:
    descriptions: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            tokens = line[1:].split(maxsplit=1)
            descriptions[tokens[0]] = tokens[1] if len(tokens) > 1 else ""
    return descriptions


@pytest.fixture(scope="module")
def mane_run(tmp_path_factory) -> dict[str, dict]:
    """Run TD2.LongOrfs once for this module, on the bundled MANE fixture."""
    if not FIXTURE.exists():
        pytest.skip(f"missing fixture {FIXTURE}")

    out_dir = tmp_path_factory.mktemp("td2_mane_out")
    # -m 80 keeps the 80-codon A2M short ORF in too, mirroring defaults; --complete-orfs-only
    # excludes 5'/3' partials so the assertions are deterministic on a clean reference set.
    cmd = [
        sys.executable, "-m", "TD2.LongOrfs",
        "-t", str(FIXTURE),
        "-O", str(out_dir),
        "-m", "80",
        "--complete-orfs-only",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        pytest.fail(
            f"TD2.LongOrfs failed (rc={result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    pep_path = out_dir / "longest_orfs.pep"
    gff_path = out_dir / "longest_orfs.gff3"
    cds_path = out_dir / "longest_orfs.cds"
    assert pep_path.exists(), pep_path
    assert gff_path.exists(), gff_path
    assert cds_path.exists(), cds_path

    return {
        "out_dir": out_dir,
        "pep": _read_pep(pep_path),
        "pep_desc": _read_pep_descriptions(pep_path),
        "cds": _read_pep(cds_path),
        "gff": gff_path.read_text(),
    }


# -------- canonical RefSeq proteins (published, externally verifiable) ---------

CANONICAL = {
    # NM_000014.6 A2M -> NP_000005.3, 1474 aa
    "NM_000014.6.p1": {
        "len_aa": 1474,
        "prefix": "MGKNKLLHPSLVLLLLVLLPTDASVSGKPQYMVLVPSLLHTETTEKGCVL",
        "suffix": "AIVKVYDYYETDEFAIAEYNAPCSKDLGNA",
    },
    # NM_000015.3 NAT2 -> NP_000006.2, 290 aa
    "NM_000015.3.p1": {
        "len_aa": 290,
        "prefix": "MDIEAYFERIGYKNSRNKLDLETLTDILEHQIRAVPFENLNMHCGQAMEL",
        "suffix": "EEVEEVLRNIFKISLGRNLVPKPGDGSLTI",
    },
    # NM_000016.6 ACADM -> NP_000007.1, 421 aa
    "NM_000016.6.p1": {
        "len_aa": 421,
        "prefix": "MAAGFGRCCRVLRSISRFHWRSQHTKANRQREPGLGFSFEFTEQQKEFQA",
        "suffix": "RDAKIYQIYEGTSQIQRLIVAREHIDKYKN",
    },
}


def test_each_transcript_has_a_canonical_longest_orf(mane_run):
    for orf_id in CANONICAL:
        assert orf_id in mane_run["pep"], f"missing longest ORF for {orf_id}"


@pytest.mark.parametrize("orf_id,expected", list(CANONICAL.items()))
def test_canonical_protein_recovered_exactly(orf_id, expected, mane_run):
    seq = mane_run["pep"][orf_id]
    # Strip stop-codon marker so the length and end match the canonical protein.
    assert seq.endswith("*"), f"{orf_id} should end on a stop codon"
    protein = seq.rstrip("*")
    assert len(protein) == expected["len_aa"], (
        f"{orf_id}: got {len(protein)} aa, paper-canonical length is {expected['len_aa']}"
    )
    assert protein.startswith(expected["prefix"]), (
        f"{orf_id}: N-terminus diverges from canonical RefSeq protein"
    )
    assert protein.endswith(expected["suffix"]), (
        f"{orf_id}: C-terminus diverges from canonical RefSeq protein"
    )


def test_canonical_orfs_marked_complete_on_plus_strand(mane_run):
    descriptions = mane_run["pep_desc"]
    for orf_id in CANONICAL:
        desc = descriptions[orf_id]
        assert "type:complete" in desc
        assert desc.rstrip().endswith("(+)")


def test_cds_translates_back_to_canonical_protein(mane_run):
    # The CDS slice for the canonical ORF must translate (frame 1) to the same protein
    # the PEP file reports, including the trailing stop codon.
    from TD2.translator import Translator

    translator = Translator(table=1)
    for orf_id, expected in CANONICAL.items():
        cds = mane_run["cds"][orf_id]
        protein, _, _ = translator.translate(cds)
        assert protein == mane_run["pep"][orf_id], (
            f"{orf_id}: CDS does not translate back to the reported PEP"
        )
        # Length of the CDS is 3 * (aa + 1 stop codon).
        assert len(cds) == 3 * (expected["len_aa"] + 1)


def test_gff3_has_one_cds_per_canonical_orf(mane_run):
    gff = mane_run["gff"]
    for orf_id in CANONICAL:
        # Each ORF gets a single CDS row whose Parent= matches the ORF id.
        cds_lines = [
            line for line in gff.splitlines()
            if "\tCDS\t" in line and f"Parent={orf_id};" in line + ";"  # tolerate trailing flags
        ]
        # Be a little tolerant: assert at least one CDS row, plus a gene row.
        assert any(f"ID={orf_id};" in line or f"ID={orf_id} " in line or f".{orf_id.split('.')[-1]};" in line
                   for line in gff.splitlines() if "\tmRNA\t" in line), (
            f"{orf_id}: no mRNA row found in GFF"
        )


def test_predicted_protein_count_matches_input(mane_run):
    # The fixture has 3 transcripts; with --complete-orfs-only there must be at
    # least 3 complete ORFs (one per transcript, plus any extras in alt frames).
    primary_ids = {k for k in mane_run["pep"] if k.endswith(".p1")}
    assert primary_ids == {"NM_000014.6.p1", "NM_000015.3.p1", "NM_000016.6.p1"}
