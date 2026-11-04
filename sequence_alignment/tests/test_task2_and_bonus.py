from pathlib import Path

from sequence_alignment.bonus_semiglobal import compare_semiglobal, format_bonus_report
from sequence_alignment.sequences import SequenceRecord, get_record, load_sequences
from sequence_alignment.task2_pairwise import compare_global_local, format_task2_report

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def _rec(acc: str, seq: str) -> SequenceRecord:
    return SequenceRecord(acc, "Test org", acc, seq, len(seq), Path(acc))


def test_overhang_scores_and_alignment_indexing_independent_of_str():
    long_ref = _rec("REF", "TTTTTTTTTT" + "ACGTACGTACGT" + "CCCCCCCCCC")
    short_q = _rec("QRY", "ACGTACGTACGT")
    pw = compare_global_local(long_ref, short_q)
    sg = compare_semiglobal(long_ref, short_q)
    assert pw.local_metrics.score == 24.0
    assert sg.semiglobal_metrics.score == 24.0
    assert pw.global_metrics.score < sg.semiglobal_metrics.score
    assert sg.semiglobal_metrics.internal_gap_bases == 0
    assert sg.semiglobal_metrics.terminal_gap_bases == 20


def test_human_vs_mouse_global_local_and_semiglobal():
    recs = load_sequences(DATA_DIR)
    hs = get_record(recs, "NM_000546.6")
    mm = get_record(recs, "NM_011640.3")
    pw = compare_global_local(hs, mm)
    sg = compare_semiglobal(hs, mm)
    assert pw.global_metrics.score == 1664.0
    assert pw.local_metrics.score == 2014.0
    assert sg.semiglobal_metrics.score == 2014.0
    assert pw.global_metrics.gap_bases > pw.local_metrics.gap_bases
    assert [pw.global_regions[r][2] for r in ("5'UTR", "CDS", "3'UTR")] == [64, 32, 743]
    assert [pw.local_regions[r][2] for r in ("5'UTR", "CDS", "3'UTR")] == [42, 32, 119]
    assert "3'UTR" in format_task2_report(pw)
    assert "-" in pw.snippet.global_block
    assert len(pw.snippet.local_block.splitlines()) == 3
    assert sg.three_way_df.height == 3
    assert "1664" in format_task2_report(pw)
    assert "Semiglobal" in format_bonus_report(sg)
