from pathlib import Path

import pytest

from sequence_alignment.sequences import load_sequences
from sequence_alignment.task1_distances import (
    compute_distances,
    format_task1_report,
    hamming_distance,
    pairwise_p_distance,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_hamming_equal_and_unequal_lengths():
    assert hamming_distance("ACGT", "AGGT") == 1
    with pytest.raises(ValueError, match="equal length"):
        hamming_distance("ACGT", "ACG")


def test_p_distance_excludes_gaps_and_handles_zero_compared_sites():
    res = pairwise_p_distance("s1", "ACGTACGT", "s2", "ACGTTCGT")
    assert (res.mismatches, res.compared_sites, res.p_distance) == (1, 8, 0.125)
    na_res = pairwise_p_distance("n1", "NNNN", "n2", "ACGT")
    assert na_res.compared_sites == 0
    assert na_res.p_distance is None


def test_p_distance_is_symmetric_on_tied_alignments():
    recs = {r.accession: r for r in load_sequences(DATA_DIR)}
    a, b = recs["NM_011640.3"], recs["NM_131327.2"]
    fwd = pairwise_p_distance(a.accession, a.seq, b.accession, b.seq)
    rev = pairwise_p_distance(b.accession, b.seq, a.accession, a.seq)
    assert (
        (fwd.mismatches, fwd.compared_sites)
        == (rev.mismatches, rev.compared_sites)
        == (648, 1755)
    )


def test_closest_pair_regions_expose_gap_deletion_bias():
    summary = compute_distances(load_sequences(DATA_DIR))
    assert summary.closest_regions["CDS"] == (1159, 198, 28)
    assert summary.closest_regions["3'UTR"] == (451, 35, 737)
    assert "3'UTR" in format_task1_report(summary)


def test_tp53_upper_triangular_matrix_and_closest_pair():
    recs = load_sequences(DATA_DIR)
    summary = compute_distances(recs)
    assert summary.closest_pair.id1 == "NM_000546.6"
    assert summary.closest_pair.id2 == "NM_011640.3"
    assert summary.closest_pair.p_distance == pytest.approx(287 / 1749, abs=1e-4)
    mat = summary.upper_matrix
    assert mat.columns == ["accession", "NM_000546.6", "NM_011640.3", "NM_131327.2"]
    assert mat.row(0) == ("NM_000546.6", "0.0000", "0.1641", "0.3809")
    assert mat.row(1)[1] == "" and mat.row(1)[2] == "0.0000"
    assert mat.row(2)[1:3] == ("", "")
    rendered = format_task1_report(summary)
    assert "0.1641" in rendered and "NM_000546.6" in rendered
