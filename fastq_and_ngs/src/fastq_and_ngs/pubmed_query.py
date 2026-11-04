"""Task 1: Query PubMed for TP53 AND cancer and save title, authors and abstract."""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez

QUERY: str = "TP53 AND cancer"
MAX_ARTICLES: int = 5
DEFAULT_EMAIL: str = "paulhondola@gmail.com"
HANDLE: str = "paulhondola"
PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Article:
    """One PubMed record."""

    pmid: str
    title: str
    authors: tuple[str, ...]
    abstract: str | None


@dataclass(frozen=True)
class SearchResult:
    """One esearch reply, including how PubMed actually read the query."""

    count: int
    pmids: list[str]
    query_translation: str


def search_pubmed(term: str, retmax: int, email: str) -> SearchResult:
    """Run an esearch against PubMed."""
    Entrez.email = email
    with Entrez.esearch(db="pubmed", term=term, retmax=retmax) as handle:
        reply = Entrez.read(handle)
    return SearchResult(
        count=int(reply["Count"]),
        pmids=[str(pmid) for pmid in reply["IdList"]],
        query_translation=str(reply.get("QueryTranslation", "")),
    )


def fetch_articles_xml(pmids: list[str], email: str) -> str:
    """Fetch full PubMed records as XML text."""
    Entrez.email = email
    with Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml") as handle:
        raw = handle.read()
    return raw.decode("utf-8") if isinstance(raw, bytes) else raw


def _text(element: ET.Element | None) -> str:
    # itertext keeps inline markup such as <i>TP53</i> in titles.
    return "".join(element.itertext()).strip() if element is not None else ""


def _author_name(author: ET.Element) -> str:
    collective = _text(author.find("CollectiveName"))
    if collective:
        return collective
    last = _text(author.find("LastName"))
    initials = _text(author.find("Initials")) or _text(author.find("ForeName"))
    return f"{last} {initials}".strip()


def _abstract(article: ET.Element) -> str | None:
    parts = []
    for section in article.findall(".//Abstract/AbstractText"):
        label = section.get("Label")
        body = _text(section)
        parts.append(f"{label}: {body}" if label else body)
    return " ".join(part for part in parts if part) or None


def parse_articles(xml: str) -> list[Article]:
    """Parse efetch XML; raises ValueError when the text is not valid XML."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise ValueError(f"PubMed reply is not valid XML: {exc}") from exc
    articles = []
    for record in root.iter("PubmedArticle"):
        names = (_author_name(a) for a in record.findall(".//AuthorList/Author"))
        articles.append(
            Article(
                pmid=_text(record.find(".//MedlineCitation/PMID")),
                title=_text(record.find(".//ArticleTitle")),
                authors=tuple(name for name in names if name),
                abstract=_abstract(record),
            )
        )
    return articles


def load_or_fetch_articles(
    cache_path: Path, refresh: bool = False, email: str = DEFAULT_EMAIL
) -> list[Article]:
    """Load cached PubMed XML, or fetch it; a bad reply never replaces the cache."""
    cache = Path(cache_path)
    if cache.exists() and not refresh:
        return parse_articles(cache.read_text(encoding="utf-8"))
    result = search_pubmed(QUERY, MAX_ARTICLES, email)
    xml = fetch_articles_xml(result.pmids, email)
    articles = parse_articles(xml)
    if not articles:
        raise ValueError(f"PubMed returned no articles for {QUERY!r}")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(xml, encoding="utf-8")
    return articles


def format_pubmed_report(articles: list[Article]) -> str:
    """Render one block per article."""
    blocks = [f"PubMed query: {QUERY} | Articles: {len(articles)}"]
    for index, article in enumerate(articles, start=1):
        authors = ", ".join(article.authors) or "(none listed)"
        blocks.append(
            f"[{index}] PMID: {article.pmid}\n"
            f"Title: {article.title}\n"
            f"Authors: {authors}\n"
            f"Abstract: {article.abstract or '(none available)'}"
        )
    return "\n\n".join(blocks) + "\n"


def main(argv: list[str] | None = None) -> int:
    """Write the PubMed report."""
    parser = argparse.ArgumentParser(
        description="Save TP53 AND cancer articles from PubMed"
    )
    parser.add_argument(
        "--refresh", action="store_true", help="Query PubMed instead of using the cache"
    )
    parser.add_argument(
        "--email", default=DEFAULT_EMAIL, help="Email sent to NCBI Entrez"
    )
    args = parser.parse_args(argv)

    articles = load_or_fetch_articles(
        PROJECT_ROOT / "data/pubmed_tp53_cancer.xml",
        refresh=args.refresh,
        email=args.email,
    )
    out = PROJECT_ROOT / f"results/pubmed_{HANDLE}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(format_pubmed_report(articles), encoding="utf-8")
    print(f"[OK] PubMed report -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
