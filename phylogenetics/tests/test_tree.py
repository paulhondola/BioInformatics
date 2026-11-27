from Bio import Phylo

from phylogenetics.tree import build_nj_tree, monophyly, splits, write_newick

# Additive toy tree ((A:1,B:2):1,(C:1,D:2)) -> pairwise path lengths.
LABELS = ["A", "B", "C", "D"]
M = [[0, 3, 3, 4], [3, 0, 4, 5], [3, 4, 0, 3], [4, 5, 3, 0]]


def _root_sides(tree):
    return {frozenset(t.name for t in c.get_terminals()) for c in tree.root.clades}


def test_nj_recovers_toy_split():
    res = build_nj_tree(LABELS, M, outgroup=("C", "D"))
    assert splits(res.tree) == {frozenset({"C", "D"})}  # normalized away from "A"
    assert res.rooting == "outgroup"
    assert _root_sides(res.tree) == {frozenset("AB"), frozenset("CD")}


def test_outgroup_not_a_split_falls_back_to_midpoint(capsys):
    res = build_nj_tree(LABELS, M, outgroup=("A", "C"))
    assert res.rooting == "midpoint"
    assert "warning" in capsys.readouterr().out.lower()


def test_monophyly_uses_unrooted_splits():
    res = build_nj_tree(LABELS, M, outgroup=("C", "D"))
    groups = {"ab": frozenset("AB"), "ac": frozenset("AC"), "solo": frozenset("A")}
    assert monophyly(res.tree, groups) == {"ab": True, "ac": False}


def test_saved_newick_reloads_clean(tmp_path):
    res = build_nj_tree(LABELS, M, outgroup=("C", "D"))
    write_newick(res.tree, tmp_path / "t.nwk")
    assert "Inner" not in (tmp_path / "t.nwk").read_text()
    reloaded = Phylo.read(tmp_path / "t.nwk", "newick")
    assert sorted(t.name for t in reloaded.get_terminals()) == LABELS
