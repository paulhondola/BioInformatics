from io import StringIO

import pytest
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylogenetics.msa import conserved_regions, robinson_foulds


def _tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def _aln(rows):
    return MultipleSeqAlignment([SeqRecord(Seq(s), id=k) for k, s in rows.items()])


CORE = "ATG" * 13 + "A"  # 40 nt, identical in every row


def _rows(b_core=CORE):
    # 60 columns; cols 10..49 identical; human row starts with 3 gaps.
    return {
        "Homo_sapiens": "---" + "C" * 7 + CORE + "G" * 10,
        "B": "T" * 10 + b_core + "C" * 10,
        "C": "G" * 10 + CORE + "T" * 10,
    }


def test_rf_identical_and_one_swap():
    t1 = _tree("((A,B),C,(D,E));")
    t2 = _tree("((A,C),B,(D,E));")
    assert robinson_foulds(t1, t1).rf == 0
    r = robinson_foulds(t1, t2)
    assert (r.rf, r.normalized) == (2, 0.5)
    assert r.shared == {frozenset("DE")}


def test_rf_rejects_different_leaves():
    with pytest.raises(ValueError):
        robinson_foulds(_tree("((A,B),C,(D,E));"), _tree("((A,B),C,(D,F));"))


def test_region_found_trimmed_and_mapped():
    s = conserved_regions(_aln(_rows()))
    [r] = s.regions
    assert (r.aln_start, r.aln_end, r.identity) == (10, 50, 1.0)
    assert r.ref_nt == (8, 47) and r.ref_codons == (3, 16)
    assert (s.columns, s.identical_columns, s.longest_run) == (60, 40, 40)


def test_gapped_column_not_identical():
    b_core = CORE[:20] + "-" + CORE[21:]  # column 30 gapped in B
    s = conserved_regions(_aln(_rows(b_core)))
    [r] = s.regions
    assert (r.aln_start, r.aln_end) == (10, 50)
    assert r.identity == pytest.approx(39 / 40)


def test_no_region_falls_back_to_best_window():
    rows = {"Homo_sapiens": "ACGT" * 15, "B": "CGTA" * 15, "C": "GTAC" * 15}
    for k in ("B", "C"):  # make column 5 the only identical one (Homo has C there)
        rows[k] = rows[k][:5] + "C" + rows[k][6:]
    s = conserved_regions(_aln(rows))
    assert s.regions == [] and s.identical_columns == 1
    assert s.best_window.aln_end - s.best_window.aln_start == 30
    assert s.best_window.aln_start == 0
