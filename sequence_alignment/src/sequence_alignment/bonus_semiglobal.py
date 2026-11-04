"""Bonus: Semiglobal pairwise alignment (unpenalized terminal gaps)."""

from __future__ import annotations

from dataclasses import dataclass

from Bio import Align
import polars as pl

from sequence_alignment.sequences import SequenceRecord
from sequence_alignment.task2_pairwise import (
    AlignmentMetrics,
    _col_to_seq_starts,
    compare_global_local,
    format_alignment_block,
    metrics_to_row,
    summarize_alignment,
)


@dataclass(frozen=True)
class SemiglobalComparison:
    """Bonus semiglobal summary and 3-way Global/Local/Semiglobal comparison."""

    seq1_id: str
    seq2_id: str
    semiglobal_metrics: AlignmentMetrics
    three_way_df: pl.DataFrame
    overhang_excerpt: str


def compare_semiglobal(
    rec1: SequenceRecord,
    rec2: SequenceRecord,
    match: float = 2.0,
    mismatch: float = -1.0,
    open_gap: float = -2.0,
    extend_gap: float = -1.0,
) -> SemiglobalComparison:
    """Run semiglobal alignment (end_gap_score=0) and compare against global and local."""
    pw = compare_global_local(
        rec1,
        rec2,
        match=match,
        mismatch=mismatch,
        open_gap=open_gap,
        extend_gap=extend_gap,
    )
    semi_aligner = Align.PairwiseAligner(
        mode="global",
        match_score=match,
        mismatch_score=mismatch,
        open_gap_score=open_gap,
        extend_gap_score=extend_gap,
        end_gap_score=0.0,
    )
    semi_aln = semi_aligner.align(rec1.seq, rec2.seq)[0]
    semi_metrics = summarize_alignment(semi_aln, "Semiglobal")

    three_way_df = pl.DataFrame(
        [
            metrics_to_row(pw.global_metrics),
            metrics_to_row(pw.local_metrics),
            metrics_to_row(semi_metrics),
        ]
    )

    g1, g2 = semi_aln[0], semi_aln[1]
    starts = _col_to_seq_starts(g1, 0, g2, 0)
    win = min(56, len(g1))
    # Highlight the 3' terminal boundary where the shorter sequence ends cleanly and end-gaps begin
    paired_cols = [i for i, (a, b) in enumerate(zip(g1, g2)) if a != "-" and b != "-"]
    last_paired = paired_cols[-1] if paired_cols else 0
    col_start = max(0, min(last_paired - win // 2, len(g1) - win))
    s1_0, s2_0 = starts[col_start]
    overhang_excerpt = format_alignment_block(
        rec1.accession,
        g1[col_start : col_start + win],
        s1_0,
        rec2.accession,
        g2[col_start : col_start + win],
        s2_0,
    )

    return SemiglobalComparison(
        seq1_id=rec1.accession,
        seq2_id=rec2.accession,
        semiglobal_metrics=semi_metrics,
        three_way_df=three_way_df,
        overhang_excerpt=overhang_excerpt,
    )


def format_bonus_report(comp: SemiglobalComparison) -> str:
    """Render Bonus semiglobal 3-way table and terminal overhang excerpt as text."""
    with pl.Config(tbl_cols=-1, tbl_rows=-1, fmt_str_lengths=100, tbl_width_chars=180):
        tbl = str(comp.three_way_df)
    return "\n".join(
        [
            f"=== BONUS: Semiglobal Alignment ({comp.seq1_id} vs {comp.seq2_id}, free 5'/3' end gaps) ===",
            tbl,
            "",
            "[3' terminal overhang boundary — unpenalized end gaps]",
            comp.overhang_excerpt,
        ]
    )
