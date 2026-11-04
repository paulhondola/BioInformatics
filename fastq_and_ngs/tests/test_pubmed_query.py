import pytest

from fastq_and_ngs import pubmed_query
from fastq_and_ngs.pubmed_query import (
    SearchResult,
    format_pubmed_report,
    load_or_fetch_articles,
    parse_articles,
)

XML = """<?xml version="1.0"?>
<PubmedArticleSet>
<PubmedArticle><MedlineCitation><PMID>1</PMID><Article>
  <ArticleTitle><i>TP53</i> in cancer</ArticleTitle>
  <Abstract>
    <AbstractText Label="BACKGROUND">first.</AbstractText>
    <AbstractText Label="RESULTS">second.</AbstractText>
  </Abstract>
  <AuthorList>
    <Author><LastName>Smith</LastName><Initials>J</Initials></Author>
    <Author><CollectiveName>TP53 Consortium</CollectiveName></Author>
  </AuthorList>
</Article></MedlineCitation></PubmedArticle>
<PubmedArticle><MedlineCitation><PMID>2</PMID><Article>
  <ArticleTitle>Bare article</ArticleTitle>
</Article></MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""
EMAIL = "test@example.org"


def _stub_network(monkeypatch, xml=XML):
    calls = []

    def search(term, retmax, email):
        calls.append((term, retmax, email))
        return SearchResult(2, ["1", "2"], "translated")

    monkeypatch.setattr(pubmed_query, "search_pubmed", search)
    monkeypatch.setattr(pubmed_query, "fetch_articles_xml", lambda pmids, email: xml)
    return calls


def _forbid_network(monkeypatch):
    def boom(*args):
        raise AssertionError("network used")

    monkeypatch.setattr(pubmed_query, "search_pubmed", boom)
    monkeypatch.setattr(pubmed_query, "fetch_articles_xml", boom)


def test_parse_handles_real_world_shapes():
    a1, a2 = parse_articles(XML)
    assert a1.pmid == "1"
    assert a1.title == "TP53 in cancer"
    assert a1.authors == ("Smith J", "TP53 Consortium")
    assert a1.abstract == "BACKGROUND: first. RESULTS: second."
    assert (a2.authors, a2.abstract) == ((), None)


def test_report_lists_title_authors_abstract():
    text = format_pubmed_report(parse_articles(XML))
    for needed in (
        "PMID: 1",
        "Title: TP53 in cancer",
        "Authors: Smith J, TP53 Consortium",
        "Abstract: BACKGROUND",
        "Authors: (none listed)",
        "Abstract: (none available)",
    ):
        assert needed in text


def test_uses_exact_query_and_cap_then_cache(tmp_path, monkeypatch):
    calls = _stub_network(monkeypatch)
    cache = tmp_path / "sub" / "cache.xml"
    first = load_or_fetch_articles(cache, email=EMAIL)
    assert calls == [("TP53 AND cancer", 5, EMAIL)]
    _forbid_network(monkeypatch)
    assert load_or_fetch_articles(cache, email=EMAIL) == first


def test_refresh_refetches(tmp_path, monkeypatch):
    cache = tmp_path / "cache.xml"
    cache.write_text(
        "<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>9</PMID></MedlineCitation></PubmedArticle></PubmedArticleSet>"
    )
    calls = _stub_network(monkeypatch)
    articles = load_or_fetch_articles(cache, refresh=True, email=EMAIL)
    assert len(calls) == 1
    assert [a.pmid for a in articles] == ["1", "2"]
    assert cache.read_text() == XML


@pytest.mark.parametrize(
    "reply", ["", "<not xml", "<PubmedArticleSet></PubmedArticleSet>"]
)
def test_bad_reply_keeps_cache(tmp_path, monkeypatch, reply):
    cache = tmp_path / "cache.xml"
    cache.write_text(XML)
    _stub_network(monkeypatch, xml=reply)
    with pytest.raises(ValueError):
        load_or_fetch_articles(cache, refresh=True, email=EMAIL)
    assert cache.read_text() == XML
