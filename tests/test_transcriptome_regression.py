"""Optional, env-gated regression test: reproduce paper-style TP/FP numbers on a
real reference transcriptome.

This test only runs when the user supplies the relevant fixtures via environment
variables. It runs the full TD2 pipeline (LongOrfs + Predict, default mode) on a
given transcriptome, then measures true positives using the same rule as Table 1
in the paper: a prediction is a TP if its terminal 20 amino acids exactly match
the C-terminus of *any* reference protein (excluding the stop codon).

Activation
----------
- TD2_TRANSCRIPTOME_FASTA          (required): path to transcripts.fna(.gz)
- TD2_TRANSCRIPTOME_REFERENCE_PEP  (required): path to reference proteins .fasta
- TD2_TRANSCRIPTOME_EXPECTED_TP    (optional): integer TP count from the paper
- TD2_TRANSCRIPTOME_EXPECTED_FP    (optional): integer FP count from the paper
- TD2_TRANSCRIPTOME_TP_TOLERANCE   (optional, default 5): allowed +/- on TP count
- TD2_TRANSCRIPTOME_FP_TOLERANCE   (optional, default 5): allowed +/- on FP count

Example (P. falciparum from the paper, Table 1: default mode, 5,229 TP, 2 FP):

    export TD2_TRANSCRIPTOME_FASTA=~/data/Pfalciparum_rna.fna
    export TD2_TRANSCRIPTOME_REFERENCE_PEP=~/data/Pfalciparum_protein.faa
    export TD2_TRANSCRIPTOME_EXPECTED_TP=5229
    export TD2_TRANSCRIPTOME_EXPECTED_FP=2
    pytest -q tests/test_transcriptome_regression.py

Without those env vars, the tests skip cleanly and pytest stays fast.
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest


REQUIRED_VARS = ("TD2_TRANSCRIPTOME_FASTA", "TD2_TRANSCRIPTOME_REFERENCE_PEP")
TD2_DIR = Path(__file__).resolve().parent.parent


def _env_path(name: str) -> Path | None:
    value = os.environ.get(name)
    if not value:
        return None
    return Path(value).expanduser()


def _read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    seqs: dict[str, str] = {}
    cur: str | None = None
    parts: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(parts)
                cur = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(parts)
    return seqs


@pytest.fixture(scope="module")
def transcriptome_run(tmp_path_factory):
    missing = [v for v in REQUIRED_VARS if not os.environ.get(v)]
    if missing:
        pytest.skip(f"transcriptome regression skipped (set {', '.join(missing)})")

    fasta = _env_path("TD2_TRANSCRIPTOME_FASTA")
    reference_pep = _env_path("TD2_TRANSCRIPTOME_REFERENCE_PEP")
    assert fasta is not None and reference_pep is not None
    assert fasta.exists(), fasta
    assert reference_pep.exists(), reference_pep

    out_dir = tmp_path_factory.mktemp("td2_transcriptome")

    # Step 1: LongOrfs at default parameters.
    rc = subprocess.run(
        [sys.executable, "-m", "TD2.LongOrfs", "-t", str(fasta), "-O", str(out_dir)],
        check=False,
    ).returncode
    assert rc == 0, "TD2.LongOrfs failed"

    # Step 2: Predict at default parameters. This invokes the `psauron` CLI internally.
    rc = subprocess.run(
        [sys.executable, "-m", "TD2.Predict", "-t", str(fasta), "-O", str(out_dir)],
        cwd=out_dir,
        check=False,
    ).returncode
    assert rc == 0, "TD2.Predict failed"

    pep_final = out_dir / (fasta.name + ".TD2.pep")
    assert pep_final.exists(), f"expected final pep at {pep_final}"

    return {
        "predicted": _read_fasta(pep_final),
        "reference": _read_fasta(reference_pep),
    }


def _terminal_kmer(seq: str, k: int = 20) -> str:
    seq = seq.rstrip("*")
    return seq[-k:] if len(seq) >= k else seq


def test_true_and_false_positive_counts_match_paper(transcriptome_run):
    """TP rule mirrors the paper (Table 1): the predicted protein's terminal 20 aa
    must exactly match the C-terminus of some reference protein."""
    predicted = transcriptome_run["predicted"]
    reference = transcriptome_run["reference"]

    # Build the set of canonical C-terminal 20-mers from the reference proteome.
    canonical_tails = {_terminal_kmer(s) for s in reference.values() if s}
    # Drop empty/very-short tails to avoid spurious matches.
    canonical_tails.discard("")

    tp = sum(1 for s in predicted.values() if _terminal_kmer(s) in canonical_tails)
    fp = len(predicted) - tp

    expected_tp = int(os.environ.get("TD2_TRANSCRIPTOME_EXPECTED_TP", "-1"))
    expected_fp = int(os.environ.get("TD2_TRANSCRIPTOME_EXPECTED_FP", "-1"))
    tp_tol = int(os.environ.get("TD2_TRANSCRIPTOME_TP_TOLERANCE", "5"))
    fp_tol = int(os.environ.get("TD2_TRANSCRIPTOME_FP_TOLERANCE", "5"))

    # Always report so the user can see the actual numbers from the run.
    print(f"\nTD2 transcriptome regression: predicted={len(predicted)} TP={tp} FP={fp}")

    if expected_tp >= 0:
        assert abs(tp - expected_tp) <= tp_tol, (
            f"TP={tp} drifted from paper's expected {expected_tp} (tolerance +/-{tp_tol})"
        )
    if expected_fp >= 0:
        assert abs(fp - expected_fp) <= fp_tol, (
            f"FP={fp} drifted from paper's expected {expected_fp} (tolerance +/-{fp_tol})"
        )
