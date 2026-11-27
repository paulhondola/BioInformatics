"""Task 1: pairwise p-distances, Jukes-Cantor correction, and the saved distance matrix."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import combinations
from math import log
from pathlib import Path

from Bio import Align

from phylogenetics.dataset import CdsRecord

UNAMBIGUOUS_DNA: frozenset[str] = frozenset("ACGT")
# Same scoring as Lab 02 Task 1: +1 / -1 / linear gap -2.
_ALIGNER = Align.PairwiseAligner(
    mode="global",
    match_score=1.0,
    mismatch_score=-1.0,
    open_gap_score=-2.0,
    extend_gap_score=-2.0,
)


@dataclass(frozen=True)
class PairDistance:
    """Distances for one pair, from its global alignment with pairwise gap deletion."""

    id1: str
    id2: str
    compared_sites: int
    mismatches: int
    p_distance: float
    jc_distance: float


@dataclass(frozen=True)
class DistanceResult:
    """All pairs plus the full square Jukes-Cantor matrix in ``labels`` order."""

    labels: list[str]
    matrix: list[list[float]]
    pairs: list[PairDistance]


def jukes_cantor(p: float) -> float:
    """JC69 distance ``-3/4 ln(1 - 4p/3)``; undefined (saturated) for ``p >= 0.75``."""
    if p >= 0.75:
        raise ValueError(f"p-distance {p:.4f} is saturated (>= 0.75); JC undefined")
    return -0.75 * log(1 - 4 * p / 3) if p > 0 else 0.0


def pair_distance(id1: str, seq1: str, id2: str, seq2: str) -> PairDistance:
    """Globally align two sequences and count mismatches over aligned ACGT/ACGT columns."""
    # Tied alignments differ by argument order; sorting by ID keeps d(a, b) == d(b, a).
    if id2 < id1:
        id1, seq1, id2, seq2 = id2, seq2, id1, seq1
    # ponytail: first of many co-optimal alignments, so counts depend on the tie-break.
    alignment = _ALIGNER.align(seq1, seq2)[0]
    compared = mismatches = 0
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        for x, y in zip(seq1[a0:a1], seq2[b0:b1]):
            if x in UNAMBIGUOUS_DNA and y in UNAMBIGUOUS_DNA:
                compared += 1
                mismatches += x != y
    if compared == 0:
        raise ValueError(f"No comparable sites between {id1} and {id2}")
    p = mismatches / compared
    return PairDistance(id1, id2, compared, mismatches, p, jukes_cantor(p))


def compute_distances(records: list[CdsRecord]) -> DistanceResult:
    """All-pairs distances and the symmetric JC matrix."""
    labels = [r.label for r in records]
    index = {label: i for i, label in enumerate(labels)}
    matrix = [[0.0] * len(labels) for _ in labels]
    pairs: list[PairDistance] = []
    for a, b in combinations(records, 2):
        d = pair_distance(a.label, a.seq, b.label, b.seq)
        pairs.append(d)
        i, j = index[d.id1], index[d.id2]
        matrix[i][j] = matrix[j][i] = d.jc_distance
    return DistanceResult(labels, matrix, pairs)


def write_matrix_csv(result: DistanceResult, path: Path) -> None:
    """Full square JC matrix, labels as header row and first column, 4 decimals."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["label", *result.labels])
        for label, row in zip(result.labels, result.matrix):
            w.writerow([label, *(f"{v:.4f}" for v in row)])


def write_pairs_csv(result: DistanceResult, path: Path) -> None:
    """One row per pair with sites, mismatches, p and JC distance."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["id1", "id2", "compared_sites", "mismatches", "p_distance", "jc_distance"]
        )
        for p in result.pairs:
            w.writerow(
                [
                    p.id1,
                    p.id2,
                    p.compared_sites,
                    p.mismatches,
                    f"{p.p_distance:.4f}",
                    f"{p.jc_distance:.4f}",
                ]
            )


def format_task1_report(result: DistanceResult) -> str:
    """Matrix as an aligned text table plus the closest and farthest pairs."""
    w = max(len(label) for label in result.labels)
    lines = [
        "=== TASK 1: Jukes-Cantor distance matrix (global +1/-1/-2, pairwise gap deletion) ===",
        " " * w + " " + " ".join(f"{i + 1:>6}" for i in range(len(result.labels))),
    ]
    for i, (label, row) in enumerate(zip(result.labels, result.matrix)):
        lines.append(
            f"{label:<{w}} " + " ".join(f"{v:6.4f}" for v in row) + f"  ({i + 1})"
        )
    near = min(result.pairs, key=lambda p: p.jc_distance)
    far = max(result.pairs, key=lambda p: p.jc_distance)
    for name, p in (("Closest", near), ("Farthest", far)):
        lines.append(
            f"{name} pair: {p.id1} vs {p.id2} (p = {p.p_distance:.4f}, JC = {p.jc_distance:.4f}, "
            f"{p.mismatches}/{p.compared_sites} mismatches)"
        )
    return "\n".join(lines)
