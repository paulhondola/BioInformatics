"""Task 2: Global vs. local pairwise alignment comparison using Biopython."""

from __future__ import annotations

from dataclasses import dataclass

from Bio import Align
import polars as pl

from sequence_alignment.sequences import CDS_SPANS, SequenceRecord


def region_stats(
    alignment: Align.Alignment, cds: tuple[int, int]
) -> dict[str, tuple[int, int, int]]:
    """Return ``{region: (compared sites, mismatches, gap columns)}`` keyed on seq1 coordinates."""
    stats = {region: [0, 0, 0] for region in ("5'UTR", "CDS", "3'UTR")}
    pos = int(alignment.coordinates[0][0])
    for a, b in zip(alignment[0], alignment[1]):
        region = "5'UTR" if pos < cds[0] else "CDS" if pos < cds[1] else "3'UTR"
        if a == "-" or b == "-":
            stats[region][2] += 1
        else:
            stats[region][0] += 1
            stats[region][1] += a != b
        if a != "-":
            pos += 1
    return {region: (s[0], s[1], s[2]) for region, s in stats.items()}


def format_regions(regions: dict[str, tuple[int, int, int]]) -> str:
    """Render ``region_stats`` output on one line."""
    return " | ".join(
        f"{region}: p={mism / sites if sites else 0:.3f} ({mism}/{sites}), {gaps} gap cols"
        for region, (sites, mism, gaps) in regions.items()
    )


@dataclass(frozen=True)
class AlignmentMetrics:
    """Quantitative metrics extracted from a Biopython PairwiseAligner Alignment."""

    mode: str
    score: float
    seq1_span: tuple[int, int]
    seq2_span: tuple[int, int]
    aligned_columns: int
    matches: int
    mismatches: int
    gap_openings: int
    gap_bases: int
    terminal_gap_bases: int
    internal_gap_bases: int


@dataclass(frozen=True)
class AlignmentSnippet:
    """Illustrative alignment fragment contrasting global gap-forcing vs. local match."""

    description: str
    global_block: str
    local_block: str


@dataclass(frozen=True)
class PairwiseComparison:
    """Task 2 comparison output for global vs. local pairwise alignments."""

    seq1_id: str
    seq2_id: str
    global_metrics: AlignmentMetrics
    local_metrics: AlignmentMetrics
    summary_df: pl.DataFrame
    snippet: AlignmentSnippet
    global_regions: dict[str, tuple[int, int, int]]
    local_regions: dict[str, tuple[int, int, int]]


def _count_gap_runs(gapped: str) -> int:
    runs = 0
    in_gap = False
    for ch in gapped:
        if ch == "-" and not in_gap:
            runs += 1
            in_gap = True
        elif ch != "-":
            in_gap = False
    return runs


def summarize_alignment(alignment: Align.Alignment, mode: str) -> AlignmentMetrics:
    """Compute column, match, mismatch, and terminal/internal gap counts from ``alignment``."""
    g1: str = alignment[0]
    g2: str = alignment[1]
    aligned_blocks_1, aligned_blocks_2 = alignment.aligned

    if len(aligned_blocks_1) == 0:
        return AlignmentMetrics(
            mode=mode,
            score=float(alignment.score),
            seq1_span=(0, 0),
            seq2_span=(0, 0),
            aligned_columns=len(g1),
            matches=0,
            mismatches=0,
            gap_openings=0,
            gap_bases=len(g1),
            terminal_gap_bases=len(g1),
            internal_gap_bases=0,
        )

    seq1_span = (int(aligned_blocks_1[0][0]), int(aligned_blocks_1[-1][1]))
    seq2_span = (int(aligned_blocks_2[0][0]), int(aligned_blocks_2[-1][1]))

    matches = 0
    mismatches = 0
    gap_bases = 0
    paired_indices: list[int] = []

    for idx, (a, b) in enumerate(zip(g1, g2)):
        if a == "-" or b == "-":
            gap_bases += 1
        else:
            paired_indices.append(idx)
            if a == b:
                matches += 1
            else:
                mismatches += 1

    if mode == "local" or not paired_indices:
        terminal_gaps = 0 if mode == "local" else gap_bases
    else:
        first_paired = paired_indices[0]
        last_paired = paired_indices[-1]
        terminal_gaps = first_paired + (len(g1) - 1 - last_paired)

    internal_gaps = gap_bases - terminal_gaps
    gap_openings = _count_gap_runs(g1) + _count_gap_runs(g2)

    return AlignmentMetrics(
        mode=mode,
        score=float(alignment.score),
        seq1_span=seq1_span,
        seq2_span=seq2_span,
        aligned_columns=len(g1),
        matches=matches,
        mismatches=mismatches,
        gap_openings=gap_openings,
        gap_bases=gap_bases,
        terminal_gap_bases=terminal_gaps,
        internal_gap_bases=internal_gaps,
    )


def format_alignment_block(
    id1: str,
    g1: str,
    start1_0: int,
    id2: str,
    g2: str,
    start2_0: int,
) -> str:
    """Format a gapped slice as a 3-line alignment block with 1-based coordinates."""
    bases1 = len(g1) - g1.count("-")
    bases2 = len(g2) - g2.count("-")
    s1_start = start1_0 + 1 if bases1 > 0 else start1_0
    s1_end = start1_0 + bases1
    s2_start = start2_0 + 1 if bases2 > 0 else start2_0
    s2_end = start2_0 + bases2

    label_w = max(len(id1), len(id2))
    coord_w = max(len(str(s1_start)), len(str(s2_start)))
    prefix1 = f"{id1:<{label_w}} {s1_start:>{coord_w}} "
    prefix2 = f"{id2:<{label_w}} {s2_start:>{coord_w}} "
    mid_chars = [
        "|" if (a == b and a != "-") else ("." if (a != "-" and b != "-") else " ")
        for a, b in zip(g1, g2)
    ]
    mid_line = " " * len(prefix1) + "".join(mid_chars)
    return "\n".join(
        [
            f"{prefix1}{g1} {s1_end}",
            mid_line,
            f"{prefix2}{g2} {s2_end}",
        ]
    )


def _seq1_pos_to_col(g1: str, start1_0: int) -> dict[int, int]:
    mapping: dict[int, int] = {}
    pos = start1_0
    for col, ch in enumerate(g1):
        if ch != "-":
            mapping[pos] = col
            pos += 1
    return mapping


def _col_to_seq_starts(g1: str, start1_0: int, g2: str, start2_0: int) -> list[tuple[int, int]]:
    starts: list[tuple[int, int]] = []
    p1, p2 = start1_0, start2_0
    for a, b in zip(g1, g2):
        starts.append((p1, p2))
        if a != "-":
            p1 += 1
        if b != "-":
            p2 += 1
    return starts


def _extract_divergence_snippet(
    id1: str,
    id2: str,
    global_aln: Align.Alignment,
    local_aln: Align.Alignment,
    width: int = 56,
) -> AlignmentSnippet:
    glob_g1, glob_g2 = global_aln[0], global_aln[1]
    loc_g1, loc_g2 = local_aln[0], local_aln[1]

    glob_s1 = int(global_aln.aligned[0][0][0]) if len(global_aln.aligned[0]) else 0
    glob_s2 = int(global_aln.aligned[1][0][0]) if len(global_aln.aligned[1]) else 0
    loc_s1 = int(local_aln.aligned[0][0][0]) if len(local_aln.aligned[0]) else 0
    loc_s2 = int(local_aln.aligned[1][0][0]) if len(local_aln.aligned[1]) else 0

    glob_starts = _col_to_seq_starts(glob_g1, 0, glob_g2, 0)
    loc_starts = _col_to_seq_starts(loc_g1, loc_s1, loc_g2, loc_s2)
    glob_map1 = _seq1_pos_to_col(glob_g1, 0)
    glob_map2 = {p2: c for c, (_, p2) in enumerate(glob_starts) if glob_g2[c] != "-"}
    loc_map2 = {p2: c for c, (_, p2) in enumerate(loc_starts) if loc_g2[c] != "-"}

    loc_win = min(width, len(loc_g1))
    glob_win = min(width, len(glob_g1))

    best_key: tuple[int, int] | None = None
    best_loc_col = max(0, len(loc_g1) - loc_win)
    best_glob_col = 0

    # Anchor on the shorter/second sequence (seq2) to find a homologous block preserved in Local
    # that Global shatters with >=20 extra gap bases to bridge UTR length disparity.
    for p2, lc in loc_map2.items():
        if lc + loc_win > len(loc_g1) or p2 not in glob_map2:
            continue
        gc = glob_map2[p2]
        if gc + glob_win > len(glob_g1):
            continue
        l_sub1 = loc_g1[lc : lc + loc_win]
        l_sub2 = loc_g2[lc : lc + loc_win]
        l_gaps = l_sub1.count("-") + l_sub2.count("-")
        g_gaps = (
            glob_g1[gc : gc + glob_win].count("-") + glob_g2[gc : gc + glob_win].count("-")
        )
        if g_gaps - l_gaps < 20:
            continue
        l_matches = sum(a == b and a != "-" for a, b in zip(l_sub1, l_sub2))
        cand_key = (l_matches, g_gaps - l_gaps)
        if best_key is None or cand_key > best_key:
            best_key = cand_key
            best_loc_col = lc
            best_glob_col = gc

    if best_key is None:
        # Fallback: capture the 3' (or 5') boundary of the local hit where global extends with forced gaps
        loc_end_p1 = int(local_aln.aligned[0][-1][1]) - 1 if len(local_aln.aligned[0]) else 0
        boundary_col = glob_map1.get(loc_end_p1, 0)
        best_glob_col = max(0, min(boundary_col - glob_win // 2, len(glob_g1) - glob_win))
        best_loc_col = max(0, len(loc_g1) - loc_win)
        if "-" not in glob_g1[best_glob_col : best_glob_col + glob_win] and "-" not in glob_g2[
            best_glob_col : best_glob_col + glob_win
        ]:
            best_glob_col = max(
                range(0, max(1, len(glob_g1) - glob_win + 1)),
                key=lambda c: glob_g1[c : c + glob_win].count("-")
                + glob_g2[c : c + glob_win].count("-"),
            )

    gs1, gs2 = glob_starts[best_glob_col]
    ls1, ls2 = loc_starts[best_loc_col]
    global_block = format_alignment_block(
        id1,
        glob_g1[best_glob_col : best_glob_col + glob_win],
        gs1,
        id2,
        glob_g2[best_glob_col : best_glob_col + glob_win],
        gs2,
    )
    local_block = format_alignment_block(
        id1,
        loc_g1[best_loc_col : best_loc_col + loc_win],
        ls1,
        id2,
        loc_g2[best_loc_col : best_loc_col + loc_win],
        ls2,
    )
    return AlignmentSnippet(
        description=(
            f"Illustrative {loc_win}-col window contrasting global gap inflation vs. clean local homology"
        ),
        global_block=global_block,
        local_block=local_block,
    )


def metrics_to_row(m: AlignmentMetrics) -> dict[str, object]:
    """Convert ``AlignmentMetrics`` into a Polars-ready dictionary."""
    return {
        "mode": m.mode,
        "score": m.score,
        "seq1_span": f"{m.seq1_span[0]}..{m.seq1_span[1]}",
        "seq2_span": f"{m.seq2_span[0]}..{m.seq2_span[1]}",
        "columns": m.aligned_columns,
        "matches": m.matches,
        "mismatches": m.mismatches,
        "gap_openings": m.gap_openings,
        "gap_bases": m.gap_bases,
        "internal_gaps": m.internal_gap_bases,
        "terminal_gaps": m.terminal_gap_bases,
    }


def compare_global_local(
    rec1: SequenceRecord,
    rec2: SequenceRecord,
    match: float = 2.0,
    mismatch: float = -1.0,
    open_gap: float = -2.0,
    extend_gap: float = -1.0,
) -> PairwiseComparison:
    """Run global and local Biopython alignments on two records and compare them."""
    global_aligner = Align.PairwiseAligner(
        mode="global",
        match_score=match,
        mismatch_score=mismatch,
        open_gap_score=open_gap,
        extend_gap_score=extend_gap,
    )
    local_aligner = Align.PairwiseAligner(
        mode="local",
        match_score=match,
        mismatch_score=mismatch,
        open_gap_score=open_gap,
        extend_gap_score=extend_gap,
    )

    global_aln = global_aligner.align(rec1.seq, rec2.seq)[0]
    local_aln = local_aligner.align(rec1.seq, rec2.seq)[0]

    g_metrics = summarize_alignment(global_aln, "Global")
    l_metrics = summarize_alignment(local_aln, "Local")
    summary_df = pl.DataFrame([metrics_to_row(g_metrics), metrics_to_row(l_metrics)])
    snippet = _extract_divergence_snippet(rec1.accession, rec2.accession, global_aln, local_aln)
    cds = CDS_SPANS.get(rec1.accession)

    return PairwiseComparison(
        seq1_id=rec1.accession,
        seq2_id=rec2.accession,
        global_metrics=g_metrics,
        local_metrics=l_metrics,
        summary_df=summary_df,
        snippet=snippet,
        global_regions=region_stats(global_aln, cds) if cds else {},
        local_regions=region_stats(local_aln, cds) if cds else {},
    )


def format_task2_report(comp: PairwiseComparison) -> str:
    """Render Task 2 comparison table and alignment snippets as text."""
    with pl.Config(tbl_cols=-1, tbl_rows=-1, fmt_str_lengths=100, tbl_width_chars=180):
        tbl = str(comp.summary_df)
    return "\n".join(
        [
            f"=== TASK 2: Global vs. Local Pairwise Alignment ({comp.seq1_id} vs {comp.seq2_id}, +2/-1, gap -2/-1) ===",
            tbl,
            *(
                [
                    f"By {comp.seq1_id} region (global): {format_regions(comp.global_regions)}",
                    f"By {comp.seq1_id} region (local):  {format_regions(comp.local_regions)}",
                ]
                if comp.global_regions
                else []
            ),
            "",
            f"[Global forced-gap fragment — {comp.snippet.description}]",
            comp.snippet.global_block,
            "",
            "[Local conserved match fragment]",
            comp.snippet.local_block,
        ]
    )
