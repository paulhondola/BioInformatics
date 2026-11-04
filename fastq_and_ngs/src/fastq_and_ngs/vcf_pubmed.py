"""Task 4: Read a VCF and search PubMed for each variant."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path

from fastq_and_ngs.pubmed_query import (
    DEFAULT_EMAIL,
    HANDLE,
    PROJECT_ROOT,
    fetch_articles_xml,
    parse_articles,
    search_pubmed,
)

MAX_LISTED: int = 3
MIN_VARIANTS: int = 2


@dataclass(frozen=True)
class Variant:
    """One VCF record; ``vid`` is None when the ID column is '.'."""

    chrom: str
    pos: int
    vid: str | None
    ref: str
    alt: str


@dataclass(frozen=True)
class VariantHit:
    """PubMed result for one variant."""

    variant: Variant
    query: str
    query_translation: str
    count: int
    articles: tuple[tuple[str, str], ...]


def parse_vcf(path: Path) -> list[Variant]:
    """Read VCF records, splitting on any whitespace so space-aligned files also work."""
    variants: list[Variant] = []
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 5:
            raise ValueError(f"VCF line {number} has fewer than 5 columns: {line!r}")
        chrom, pos, vid, ref, alt = fields[:5]
        variants.append(Variant(chrom, int(pos), None if vid == "." else vid, ref, alt))
    if len(variants) < MIN_VARIANTS:
        raise ValueError(f"VCF must hold at least {MIN_VARIANTS} variants, found {len(variants)}")
    return variants


def _chrom_label(variant: Variant) -> str:
    return f"chr{variant.chrom.removeprefix('chr')}"


def build_query(variant: Variant) -> str:
    """Search the ID when there is one, otherwise chr<CHROM>:<POS> AND TP53."""
    if variant.vid:
        return " OR ".join(variant.vid.split(";"))
    return f"{_chrom_label(variant)}:{variant.pos} AND TP53"


def position_dropped(variant: Variant, translation: str) -> bool:
    """True when PubMed ignored the position of a no-ID variant's query."""
    return variant.vid is None and str(variant.pos) not in translation


def _fetch_hit(query: str, email: str) -> dict:
    result = search_pubmed(query, MAX_LISTED, email)
    articles = parse_articles(fetch_articles_xml(result.pmids, email)) if result.pmids else []
    if result.pmids and not articles:
        # ponytail: also rejects a result made only of PubMed book records; parse those if it ever happens.
        raise ValueError(f"PubMed returned no articles for PMIDs {result.pmids} (query {query!r})")
    return {
        "count": result.count,
        "query_translation": result.query_translation,
        "articles": [[article.pmid, article.title] for article in articles],
    }


def load_or_fetch_hits(
    variants: list[Variant],
    cache_path: Path,
    refresh: bool = False,
    email: str = DEFAULT_EMAIL,
) -> list[VariantHit]:
    """Search each variant, reusing cached queries; a failure part-way never changes the cache."""
    cache_file = Path(cache_path)
    cache: dict[str, dict] = (
        json.loads(cache_file.read_text(encoding="utf-8")) if cache_file.exists() else {}
    )
    queries = [build_query(variant) for variant in variants]
    missing = [query for query in dict.fromkeys(queries) if refresh or query not in cache]
    if missing:
        cache.update({query: _fetch_hit(query, email) for query in missing})
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(cache, indent=2) + "\n", encoding="utf-8")
    return [
        VariantHit(
            variant=variant,
            query=query,
            query_translation=cache[query]["query_translation"],
            count=cache[query]["count"],
            articles=tuple((pmid, title) for pmid, title in cache[query]["articles"]),
        )
        for variant, query in zip(variants, queries)
    ]


def format_variants_report(hits: list[VariantHit]) -> str:
    """Render one block per variant, flagging hits that do not concern its position."""
    blocks = [f"VCF to PubMed | Variants: {len(hits)}"]
    for hit in hits:
        v = hit.variant
        lines = [
            f"Variant: {v.chrom}:{v.pos} {v.ref}>{v.alt} (ID: {v.vid or 'none'})",
            f"Query: {hit.query}",
            f"PubMed read it as: {hit.query_translation}",
            f"Hits: {hit.count}",
        ]
        if position_dropped(v, hit.query_translation):
            lines.append(
                "WARNING: PubMed dropped the position from this query; "
                f"hits are about {_chrom_label(v)} and TP53, not this variant."
            )
        lines.extend(f"  PMID {pmid}: {title}" for pmid, title in hit.articles)
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks) + "\n"


def main(argv: list[str] | None = None) -> int:
    """Write the variant search report."""
    parser = argparse.ArgumentParser(description="Search PubMed for the variants in a VCF")
    parser.add_argument("--vcf", type=Path, default=PROJECT_ROOT / "data/tp53_variants.vcf")
    parser.add_argument("--refresh", action="store_true", help="Query PubMed instead of using the cache")
    parser.add_argument("--email", default=DEFAULT_EMAIL, help="Email sent to NCBI Entrez")
    args = parser.parse_args(argv)

    hits = load_or_fetch_hits(
        parse_vcf(args.vcf), PROJECT_ROOT / "data/variant_searches.json", refresh=args.refresh, email=args.email
    )
    out = PROJECT_ROOT / f"results/variants_{HANDLE}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(format_variants_report(hits), encoding="utf-8")
    print(f"[OK] Variant report -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
