"""Task 1: Pairwise sequence distances and Polars upper-triangular matrix."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from Bio import Align
import polars as pl

from sequence_alignment.sequences import CDS_SPANS, SequenceRecord
from sequence_alignment.task2_pairwise import format_regions, region_stats

UNAMBIGUOUS_DNA: frozenset[str] = frozenset("ACGT")


@dataclass(frozen=True)
class PairDistance:
    """Pairwise global-alignment p-distance metrics."""

    id1: str
    id2: str
    score: float
    mismatches: int
    compared_sites: int
    p_distance: float | None


@dataclass(frozen=True)
class DistanceSummary:
    """Task 1 summary containing pairwise metrics and Polars tables."""

    pairs: list[PairDistance]
    upper_matrix: pl.DataFrame
    pair_table: pl.DataFrame
    closest_pair: PairDistance
    closest_regions: dict[str, tuple[int, int, int]]


def hamming_distance(seq1: str, seq2: str) -> int:
    """Compute Hamming distance between two sequences of equal length."""
    if len(seq1) != len(seq2):
        raise ValueError(
            f"Hamming distance requires sequences of equal length ({len(seq1)} != {len(seq2)})"
        )
    return sum(a != b for a, b in zip(seq1.upper(), seq2.upper()))


def _global_alignment(
    seq1: str, seq2: str, match: float = 1.0, mismatch: float = -1.0, gap: float = -2.0
) -> Align.Alignment:
    aligner = Align.PairwiseAligner(
        mode="global",
        match_score=match,
        mismatch_score=mismatch,
        open_gap_score=gap,
        extend_gap_score=gap,
    )
    # ponytail: first of many co-optimal alignments, so counts depend on the tie-break.
    return aligner.align(seq1, seq2)[0]


def pairwise_p_distance(
    id1: str,
    seq1: str,
    id2: str,
    seq2: str,
    match: float = 1.0,
    mismatch: float = -1.0,
    gap: float = -2.0,
) -> PairDistance:
    """Align two sequences globally and compute gap/ambiguity-deleted p-distance."""
    # Tied alignments differ by argument order; sorting by ID keeps d(a, b) == d(b, a).
    if id2 < id1:
        id1, seq1, id2, seq2 = id2, seq2, id1, seq1
    a_up = seq1.upper()
    b_up = seq2.upper()
    alignment = _global_alignment(a_up, b_up, match, mismatch, gap)

    mismatches = 0
    compared_sites = 0
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        for x, y in zip(a_up[a0:a1], b_up[b0:b1]):
            if x in UNAMBIGUOUS_DNA and y in UNAMBIGUOUS_DNA:
                compared_sites += 1
                if x != y:
                    mismatches += 1

    p_dist = (mismatches / compared_sites) if compared_sites > 0 else None
    return PairDistance(
        id1=id1,
        id2=id2,
        score=float(alignment.score),
        mismatches=mismatches,
        compared_sites=compared_sites,
        p_distance=p_dist,
    )


def compute_distances(records: list[SequenceRecord]) -> DistanceSummary:
    """Compute pairwise p-distances and upper-triangular matrix across >=3 records."""
    if len(records) < 3:
        raise ValueError("compute_distances requires at least 3 sequences")

    pairs: list[PairDistance] = []
    lookup: dict[tuple[int, int], PairDistance] = {}
    for i, j in combinations(range(len(records)), 2):
        res = pairwise_p_distance(
            records[i].accession,
            records[i].seq,
            records[j].accession,
            records[j].seq,
        )
        pairs.append(res)
        lookup[(i, j)] = res

    accs = [r.accession for r in records]
    matrix_rows: list[dict[str, str]] = []
    for i, row_acc in enumerate(accs):
        row: dict[str, str] = {"accession": row_acc}
        for j, col_acc in enumerate(accs):
            if i == j:
                row[col_acc] = "0.0000"
            elif i < j:
                p = lookup[(i, j)].p_distance
                row[col_acc] = f"{p:.4f}" if p is not None else "NA"
            else:
                row[col_acc] = ""
        matrix_rows.append(row)

    upper_matrix = pl.DataFrame(matrix_rows)
    pair_table = pl.DataFrame(
        [
            {
                "id1": p.id1,
                "id2": p.id2,
                "global_score": p.score,
                "mismatches": p.mismatches,
                "compared_sites": p.compared_sites,
                "p_distance": f"{p.p_distance:.4f}" if p.p_distance is not None else "NA",
            }
            for p in pairs
        ]
    )

    valid_pairs = [p for p in pairs if p.p_distance is not None]
    closest_pair = min(valid_pairs, key=lambda p: p.p_distance) if valid_pairs else pairs[0]
    seqs = {r.accession: r.seq for r in records}
    cds = CDS_SPANS.get(closest_pair.id1)
    closest_regions = (
        region_stats(_global_alignment(seqs[closest_pair.id1], seqs[closest_pair.id2]), cds)
        if cds
        else {}
    )
    return DistanceSummary(
        pairs=pairs,
        upper_matrix=upper_matrix,
        pair_table=pair_table,
        closest_pair=closest_pair,
        closest_regions=closest_regions,
    )


def format_task1_report(summary: DistanceSummary) -> str:
    """Render Task 1 distance matrix and closest-pair summary as text."""
    with pl.Config(tbl_cols=-1, tbl_rows=-1, fmt_str_lengths=100, tbl_width_chars=180):
        mat_str = str(summary.upper_matrix)
        tbl_str = str(summary.pair_table)
    cp = summary.closest_pair
    cp_val = f"{cp.p_distance:.4f}" if cp.p_distance is not None else "NA"
    return "\n".join(
        [
            "=== TASK 1: Pairwise p-distance Matrix (global +1/-1/-2, pairwise gap deletion) ===",
            mat_str,
            "",
            "Pairwise alignment details:",
            tbl_str,
            f"Closest pair: {cp.id1} vs {cp.id2} (p-distance = {cp_val}, {cp.mismatches}/{cp.compared_sites} mismatches)",
            *([f"By {cp.id1} region: {format_regions(summary.closest_regions)}"] if summary.closest_regions else []),
        ]
    )
