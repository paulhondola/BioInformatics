from pathlib import Path

import pytest

from fastq_and_ngs import vcf_pubmed
from fastq_and_ngs.pubmed_query import SearchResult
from fastq_and_ngs.vcf_pubmed import (
    Variant,
    VariantHit,
    build_query,
    format_variants_report,
    load_or_fetch_hits,
    parse_vcf,
    position_dropped,
)

VCF = Path(__file__).resolve().parents[1] / "data/tp53_variants.vcf"
EXPECTED = [
    ("17", 7676154, "rs1042522", "G", "C"),
    ("17", 7675088, "rs28934578", "C", "T"),
    ("17", 7674220, "rs11540652", "C", "T"),
    ("17", 7675091, None, "C", "G"),
    ("17", 7674224, None, "T", "G"),
]
RS = Variant("17", 1, "rs5", "A", "T")
NO_ID = Variant("17", 7675091, None, "C", "G")
DROPPED = '"chr17"[All Fields] AND "tp53"[All Fields]'
KEPT = '"chr17:7675091"[All Fields] AND "tp53"[All Fields]'
WARNING = "hits are about chr17 and TP53, not this variant"
EMAIL = "test@example.org"


def _articles_xml(pmids):
    body = "".join(
        f"<PubmedArticle><MedlineCitation><PMID>{p}</PMID><Article><ArticleTitle>Title {p}</ArticleTitle>"
        "</Article></MedlineCitation></PubmedArticle>"
        for p in pmids
    )
    return f"<PubmedArticleSet>{body}</PubmedArticleSet>"


def _stub_network(monkeypatch):
    searched = []

    def search(term, retmax, email):
        searched.append((term, retmax))
        return SearchResult(10, ["1", "2", "3"], f"translated {term}")

    monkeypatch.setattr(vcf_pubmed, "search_pubmed", search)
    monkeypatch.setattr(vcf_pubmed, "fetch_articles_xml", lambda pmids, email: _articles_xml(pmids))
    return searched


def _forbid_network(monkeypatch):
    def boom(*args):
        raise AssertionError("network used")

    monkeypatch.setattr(vcf_pubmed, "search_pubmed", boom)
    monkeypatch.setattr(vcf_pubmed, "fetch_articles_xml", boom)


def test_committed_vcf_matches_spec():
    assert [(v.chrom, v.pos, v.vid, v.ref, v.alt) for v in parse_vcf(VCF)] == EXPECTED
    text = VCF.read_text()
    assert "##reference=GRCh38" in text and "ClinVar" in text
    for line in (l for l in text.splitlines() if not l.startswith("##")):
        assert line.count("\t") >= 7


def test_parse_accepts_space_aligned_vcf(tmp_path):
    path = tmp_path / "v.vcf"
    path.write_text("##x\n#CHROM POS ID REF ALT\n17   10  rs1  G C . PASS X\n17   20  .    A T . PASS X\n")
    assert [v.vid for v in parse_vcf(path)] == ["rs1", None]


def test_parse_rejects_fewer_than_two_variants(tmp_path):
    path = tmp_path / "v.vcf"
    path.write_text("17\t10\trs1\tG\tC\n")
    with pytest.raises(ValueError, match="at least 2"):
        parse_vcf(path)


def test_parse_rejects_short_line(tmp_path):
    path = tmp_path / "v.vcf"
    path.write_text("##x\n17\t10\trs1\tG\tC\n17\t20\n")
    with pytest.raises(ValueError, match="line 3"):
        parse_vcf(path)


def test_build_query():
    assert build_query(RS) == "rs5"
    assert build_query(NO_ID) == "chr17:7675091 AND TP53"
    assert build_query(Variant("chr17", 7675091, None, "C", "G")) == "chr17:7675091 AND TP53"
    assert build_query(Variant("17", 1, "rs5;rs6", "A", "T")) == "rs5 OR rs6"


def test_position_dropped():
    assert position_dropped(NO_ID, DROPPED) is True
    assert position_dropped(NO_ID, KEPT) is False
    assert position_dropped(RS, '"rs5"[All Fields]') is False


def test_report_shows_translation_zero_hits_and_warning():
    dropped = VariantHit(NO_ID, "chr17:7675091 AND TP53", DROPPED, 13, (("7", "A paper"),))
    kept = VariantHit(NO_ID, "chr17:7675091 AND TP53", KEPT, 1, ())
    zero = VariantHit(RS, "rs5", '"rs5"[All Fields]', 0, ())
    text = format_variants_report([dropped, zero])
    for needed in ("17:7675091 C>G", "chr17:7675091 AND TP53", DROPPED, "Hits: 13", "PMID 7: A paper", "Hits: 0", WARNING):
        assert needed in text
    assert text.count("PMID ") == 1
    assert WARNING not in format_variants_report([kept, zero])


def test_cache_avoids_network_and_caps_titles(tmp_path, monkeypatch):
    searched = _stub_network(monkeypatch)
    cache = tmp_path / "sub" / "hits.json"
    first = load_or_fetch_hits([RS, NO_ID], cache, email=EMAIL)
    assert searched == [("rs5", 3), ("chr17:7675091 AND TP53", 3)]
    assert first[0].count == 10
    assert first[0].articles == (("1", "Title 1"), ("2", "Title 2"), ("3", "Title 3"))
    _forbid_network(monkeypatch)
    assert load_or_fetch_hits([RS, NO_ID], cache, email=EMAIL) == first


def test_zero_hits_skip_article_fetch(tmp_path, monkeypatch):
    def fetch(pmids, email):
        raise AssertionError("fetched with no PMIDs")

    monkeypatch.setattr(vcf_pubmed, "search_pubmed", lambda t, r, e: SearchResult(0, [], "x"))
    monkeypatch.setattr(vcf_pubmed, "fetch_articles_xml", fetch)
    hit = load_or_fetch_hits([RS], tmp_path / "h.json", email=EMAIL)[0]
    assert (hit.count, hit.articles) == (0, ())


def test_new_vcf_record_fetches_only_new_query(tmp_path, monkeypatch):
    cache = tmp_path / "hits.json"
    _stub_network(monkeypatch)
    load_or_fetch_hits([RS], cache, email=EMAIL)
    searched = _stub_network(monkeypatch)
    hits = load_or_fetch_hits([RS, NO_ID], cache, email=EMAIL)
    assert searched == [("chr17:7675091 AND TP53", 3)]
    assert len(hits) == 2


def test_refresh_failure_keeps_cache(tmp_path, monkeypatch):
    cache = tmp_path / "hits.json"
    _stub_network(monkeypatch)
    load_or_fetch_hits([RS, NO_ID], cache, email=EMAIL)
    before = cache.read_bytes()
    calls = []

    def flaky(term, retmax, email):
        calls.append(term)
        if len(calls) == 2:
            raise OSError("network down")
        return SearchResult(99, [], "changed")

    monkeypatch.setattr(vcf_pubmed, "search_pubmed", flaky)
    with pytest.raises(OSError):
        load_or_fetch_hits([RS, NO_ID], cache, refresh=True, email=EMAIL)
    assert cache.read_bytes() == before
