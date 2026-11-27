import csv

import pytest

from phylogenetics.dataset import PROJECT_ROOT, load_or_fetch_dataset
from phylogenetics.distances import (
    compute_distances,
    jukes_cantor,
    pair_distance,
    write_matrix_csv,
    write_pairs_csv,
)


def test_jukes_cantor_values():
    assert jukes_cantor(0.0) == 0.0
    assert jukes_cantor(0.1) == pytest.approx(0.107326, abs=1e-6)
    with pytest.raises(ValueError):
        jukes_cantor(0.75)


def test_pair_distance_toy_and_symmetric():
    d = pair_distance("a", "ACGTACGTAC", "b", "ACGTTCGTAC")
    assert (d.compared_sites, d.mismatches, d.p_distance) == (10, 1, 0.1)
    assert pair_distance("b", "ACGTTCGTAC", "a", "ACGTACGTAC") == d


def test_n_columns_not_compared():
    d = pair_distance("a", "ACGTACGTAC", "b", "ACGTNCGTAC")
    assert (d.compared_sites, d.mismatches) == (9, 0)


def test_matrix_from_committed_dataset(tmp_path):
    res = compute_distances(load_or_fetch_dataset(PROJECT_ROOT / "data/tp53_cds.fasta"))
    n = len(res.labels)
    assert n == 12 and len(res.pairs) == 66
    assert all(
        res.matrix[i][i] == 0 and res.matrix[i][j] == res.matrix[j][i]
        for i in range(n)
        for j in range(n)
    )
    write_matrix_csv(res, tmp_path / "m.csv")
    rows = list(csv.reader((tmp_path / "m.csv").open()))
    assert rows[0][1:] == res.labels and [r[0] for r in rows[1:]] == res.labels
    write_pairs_csv(res, tmp_path / "p.csv")
    header, *body = list(csv.reader((tmp_path / "p.csv").open()))
    assert header == [
        "id1",
        "id2",
        "compared_sites",
        "mismatches",
        "p_distance",
        "jc_distance",
    ]
    assert len(body) == 66
