"""Tests for the FDR math (Eqs. 3-8 in the paper) and other Predict helpers."""
import math

import numpy as np
import pytest
from scipy import stats

from TD2.Predict import (
    ER,
    VR,
    find_encapsulated_intervals,
    gc_content,
    length_at_fdr,
)


# --------------------------------------------------------------------- ER / VR


def test_ER_matches_paper_equation_3():
    # Paper Eq. 3: E(LSS) = log_{1/p}(n*q) + gamma / ln(1/p) - 1/2
    # (with r1, eps1 ~ 0).
    n, p = 1000.0, 0.61
    q = 1 - p
    expected = math.log(n * q, 1 / p) + np.euler_gamma / math.log(1 / p) - 0.5
    assert ER(n, p) == pytest.approx(expected)


def test_VR_matches_paper_equation_4():
    # Paper Eq. 4: Var(LSS) = pi^2 / (6 ln^2(1/p)) + 1/12
    n, p = 1000.0, 0.61
    expected = math.pi**2 / (6 * math.log(1 / p) ** 2) + 1 / 12
    assert VR(n, p) == pytest.approx(expected)


def test_ER_returns_zero_for_degenerate_inputs():
    assert ER(0, 0.5) == 0
    assert ER(100, 1.0) == 0  # p>=1 short-circuit
    assert ER(100, 1.1) == 0  # also p>=1


def test_VR_returns_zero_for_n_zero():
    assert VR(0, 0.5) == 0


# ------------------------------------------------------------- length_at_fdr


def test_length_at_fdr_zero_GC_uses_AT_only_stops():
    # p_stop = (0.5*(1-GC))^3 + 2*(0.5*(1-GC))^2 * (0.5*GC). At GC=0, p_stop = 1/8.
    # Translates to a finite cutoff much smaller than 2**31.
    L = length_at_fdr(transcript_length=1000, GC=0.0, fdr=0.01)
    assert 0 < L < 1000


def test_length_at_fdr_increases_with_GC_content():
    # Higher GC -> rarer stop codons -> longer expected LSS -> larger cutoff.
    low_gc = length_at_fdr(1000, GC=0.30, fdr=0.05)
    high_gc = length_at_fdr(1000, GC=0.70, fdr=0.05)
    assert high_gc > low_gc


def test_length_at_fdr_increases_with_transcript_length():
    # Longer transcript -> more "coin flips" -> longer expected longest stopless run.
    short = length_at_fdr(500, GC=0.5, fdr=0.05)
    long_ = length_at_fdr(50_000, GC=0.5, fdr=0.05)
    assert long_ > short


def test_length_at_fdr_decreases_with_stricter_fdr():
    # A stricter (smaller) FDR demands a higher cutoff to be rare-enough-by-chance.
    lenient = length_at_fdr(10_000, GC=0.5, fdr=0.10)
    strict = length_at_fdr(10_000, GC=0.5, fdr=0.001)
    assert strict > lenient


def test_length_at_fdr_returns_huge_for_GC_one():
    # GC = 1.0 -> no stop codons can ever appear -> cutoff saturates to the documented sentinel.
    assert length_at_fdr(10_000, GC=1.0, fdr=0.05) == 2**31


def test_length_at_fdr_matches_inverse_gumbel_construction():
    # Verifies the formula in Predict.length_at_fdr against the paper's explicit construction
    # (Eqs. 3-8): we should land on the same nt-length threshold.
    transcript_length = 5000
    GC = 0.42
    fdr = 0.05

    p_TAA = (0.5 * (1 - GC)) ** 3
    p_TAG = (0.5 * (1 - GC)) ** 2 * (0.5 * GC)
    p_TGA = (0.5 * (1 - GC)) ** 2 * (0.5 * GC)
    p_stop = p_TAA + p_TAG + p_TGA
    coin_flips = transcript_length / 3
    mean = ER(coin_flips, 1 - p_stop)
    var = VR(coin_flips, 1 - p_stop)
    std = math.sqrt(var)
    fdr_multiframe = 1 - (1 - fdr) ** (1 / 6)
    beta = std * math.sqrt(6) / math.pi
    mu = mean - np.euler_gamma * beta
    expected_nt = stats.gumbel_r.ppf(1 - fdr_multiframe, loc=mu, scale=beta) * 3

    assert length_at_fdr(transcript_length, GC, fdr) == pytest.approx(expected_nt)


# --------------------------------------------------------- gc_content


@pytest.mark.parametrize("seq,expected", [
    ("ATAT", 0.0),
    ("GCGC", 1.0),
    ("ACGT", 0.5),
    ("acgt", 0.5),
    ("GGCCAAA", 4 / 7),
])
def test_gc_content_known_values(seq, expected):
    assert gc_content(seq) == pytest.approx(expected)


# ------------------------------------------------- find_encapsulated_intervals


def test_no_encapsulation_when_intervals_are_disjoint():
    intervals = [("a", (10, 20)), ("b", (30, 40)), ("c", (50, 60))]
    assert find_encapsulated_intervals(intervals) == set()


def test_inner_interval_is_marked_encapsulated():
    intervals = [("outer", (0, 100)), ("inner", (10, 90))]
    assert find_encapsulated_intervals(intervals) == {"inner"}


def test_only_strict_containment_counts_in_general_order():
    # Two intervals sharing one endpoint: the second is still contained in the first.
    intervals = [("outer", (0, 100)), ("flush_right", (50, 100))]
    assert find_encapsulated_intervals(intervals) == {"flush_right"}


def test_intervals_with_identical_bounds_handled():
    intervals = [("a", (10, 20)), ("b", (10, 20))]
    # Sorting tie-break (low asc, high desc) keeps the first as bigboi;
    # the duplicate ends up encapsulated.
    encapsulated = find_encapsulated_intervals(intervals)
    assert len(encapsulated) == 1
