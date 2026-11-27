import pytest

from phylogenetics import dataset
from phylogenetics.dataset import (
    GROUPS,
    PROJECT_ROOT,
    SPECIES,
    load_or_fetch_dataset,
    parse_cds_fasta,
    validate_cds,
)

GOOD = "ATG" + "GCA" * 10 + "TAA"


def test_validate_accepts_good_cds_with_n():
    validate_cds("X", "ATG" + "GNA" + "TGA")


@pytest.mark.parametrize(
    "seq",
    ["GCA" + GOOD[3:], GOOD + "A", "ATG" + "GCA" * 3, "ATG" + "GXA" + "TAA", ""],
    ids=["no-atg", "not-codons", "no-stop", "bad-char", "empty"],
)
def test_validate_rejects_bad_cds(seq):
    with pytest.raises(ValueError):
        validate_cds("X", seq)


def test_parse_requires_all_species():
    text = "".join(f">{s.label} {s.accession}\n{GOOD}\n" for s in SPECIES[:-1])
    with pytest.raises(ValueError):
        parse_cds_fasta(text)


def test_committed_dataset_parses_in_species_order():
    recs = parse_cds_fasta((PROJECT_ROOT / "data/tp53_cds.fasta").read_text())
    assert [r.label for r in recs] == [s.label for s in SPECIES]
    assert len(next(r for r in recs if r.label == "Homo_sapiens").seq) == 1182


def test_bad_refresh_keeps_cache(tmp_path, monkeypatch):
    cache = tmp_path / "tp53_cds.fasta"
    good = "".join(f">{s.label} {s.accession}\n{GOOD}\n" for s in SPECIES)
    cache.write_text(good)
    monkeypatch.setattr(
        dataset, "fetch_cds_fasta", lambda email: good.replace(GOOD, "GCA", 1)
    )
    with pytest.raises(ValueError):
        load_or_fetch_dataset(cache, refresh=True)
    assert cache.read_text() == good


def test_groups_cover_mammals_and_tetrapods():
    assert (
        GROUPS["Mammals"]
        == GROUPS["Primates"] | GROUPS["Rodents"] | GROUPS["Laurasiatheria"]
    )
    assert "Danio_rerio" not in GROUPS["Tetrapods"]
    assert "Gallus_gallus" in GROUPS["Tetrapods"]
