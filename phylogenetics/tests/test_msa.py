import pytest

from phylogenetics import msa
from phylogenetics.dataset import PROJECT_ROOT, CdsRecord, load_or_fetch_dataset
from phylogenetics.msa import (
    load_or_fetch_clustalo,
    parse_validated_msa,
    parse_validated_tree,
)

RECS = [
    CdsRecord("Homo_sapiens", "X1", "ATGAAATAA"),
    CdsRecord("B", "X2", "ATGAAGTAA"),
    CdsRecord("C", "X3", "ATGTAA"),
]


def _clustal(rows):
    body = "\n".join(f"{k:<14}{v}" for k, v in rows.items())
    return f"CLUSTAL O(1.2.4) multiple sequence alignment\n\n\n{body}\n"


GOOD_MSA = _clustal({"Homo_sapiens": "ATGAAATAA", "B": "ATGAAGTAA", "C": "ATG---TAA"})
GOOD_TREE = "((Homo_sapiens:0.1,B:0.2):0.05,C:0.3);"


def test_msa_parses_in_record_order():
    aln = parse_validated_msa(GOOD_MSA, list(reversed(RECS)))
    assert [r.id for r in aln] == ["C", "B", "Homo_sapiens"]


def test_msa_rejects_wrong_ids():
    bad = _clustal({"Homo_sapiens": "ATGAAATAA", "B": "ATGAAGTAA", "D": "ATG---TAA"})
    with pytest.raises(ValueError):
        parse_validated_msa(bad, RECS)


def test_msa_rejects_altered_sequence():
    bad = _clustal({"Homo_sapiens": "ATGAAATAA", "B": "ATGAAATAA", "C": "ATG---TAA"})
    with pytest.raises(ValueError):
        parse_validated_msa(bad, RECS)


def test_tree_rejects_mismatched_labels():
    with pytest.raises(ValueError):
        parse_validated_tree("((Homo_sapi:0.1,B:0.2):0.05,C:0.3);", RECS)


def test_tree_accepts_matching_labels():
    tree = parse_validated_tree(GOOD_TREE, RECS)
    assert sorted(t.name for t in tree.get_terminals()) == ["B", "C", "Homo_sapiens"]


def test_bad_refresh_keeps_both_caches(tmp_path, monkeypatch):
    msa_cache, tree_cache = tmp_path / "m.aln", tmp_path / "t.nwk"
    msa_cache.write_text(GOOD_MSA)
    tree_cache.write_text(GOOD_TREE)
    monkeypatch.setattr(
        msa, "fetch_clustalo", lambda records, email: (GOOD_MSA, "(Homo_sapiens,B);")
    )
    with pytest.raises(ValueError):
        load_or_fetch_clustalo(RECS, msa_cache, tree_cache, refresh=True)
    assert msa_cache.read_text() == GOOD_MSA
    assert tree_cache.read_text() == GOOD_TREE


def test_committed_caches_validate():
    msa_cache = PROJECT_ROOT / "data/clustalo_msa.aln"
    tree_cache = PROJECT_ROOT / "data/clustalo_phylotree.nwk"
    # Must never fall through to a live EBI fetch.
    assert msa_cache.is_file() and tree_cache.is_file(), "run --refresh-msa once"
    recs = load_or_fetch_dataset(PROJECT_ROOT / "data/tp53_cds.fasta")
    aln, tree = load_or_fetch_clustalo(recs, msa_cache, tree_cache)
    assert len(aln) == 12 and len(tree.get_terminals()) == 12
