import gzip
from pathlib import Path

import pytest

from fastq_and_ngs.ex04_fastq_stats import (
    FastqStats,
    compute_fastq_stats,
    format_qc_report,
    write_qc_outputs,
    write_qc_plot,
)

REAL_FASTQ = Path(__file__).resolve().parents[1] / "data/work/paulhondola/lab03/your_reads.fastq.gz"
REAL_REPORT = "Reads: 81220\nMean length: 142.00\nN rate: 0.0000\nMean Phred: 37.30\n"
TWO_READS = "@r1\nACGN\n+\nIIII\n@r2\nAC\n+\n++\n"


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    if name.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("name", ["reads.fastq", "reads.fastq.gz"])
def test_stats_on_plain_and_gzipped(tmp_path, name):
    s = compute_fastq_stats(_write(tmp_path, name, TWO_READS))
    assert (s.reads, s.mean_length) == (2, 3.0)
    assert s.n_rate == pytest.approx(1 / 6)
    assert s.mean_phred == pytest.approx(180 / 6)
    assert s.read_lengths == [4, 2]
    assert s.read_mean_phreds == [40.0, 10.0]


def test_lowercase_n_counts(tmp_path):
    s = compute_fastq_stats(_write(tmp_path, "r.fastq", "@r\nacgn\n+\nIIII\n"))
    assert s.n_rate == 0.25


def test_report_has_exact_template_lines():
    stats = FastqStats(81220, 142.0, 0.0, 37.3041, [142], [37.3])
    assert format_qc_report(stats) == REAL_REPORT


def test_empty_file_gives_zeros(tmp_path):
    s = compute_fastq_stats(_write(tmp_path, "empty.fastq", ""))
    assert (s.reads, s.mean_length, s.n_rate, s.mean_phred) == (0, 0.0, 0.0, 0.0)
    assert write_qc_plot(s, tmp_path / "out" / "plot.png").stat().st_size > 0


def test_missing_file_raises_and_keeps_report(tmp_path):
    report = tmp_path / "report.txt"
    report.write_text("real report")
    with pytest.raises(FileNotFoundError):
        write_qc_outputs(tmp_path / "absent.fastq", report, tmp_path / "plot.png")
    assert report.read_text() == "real report"
    assert not (tmp_path / "plot.png").exists()


def _corrupt_gzip(tmp_path: Path) -> Path:
    whole = _write(tmp_path, "whole.fastq.gz", TWO_READS * 200)
    cut = tmp_path / "cut.fastq.gz"
    cut.write_bytes(whole.read_bytes()[:-20])
    return cut


def _corrupt_lengths(tmp_path: Path) -> Path:
    return _write(tmp_path, "bad.fastq", "@r\nACGT\n+\nII\n")


@pytest.mark.parametrize("make", [_corrupt_gzip, _corrupt_lengths])
def test_corrupt_fastq_keeps_report_and_plot(tmp_path, make):
    report, plot = tmp_path / "report.txt", tmp_path / "plot.png"
    report.write_text("real report")
    plot.write_bytes(b"real plot")
    with pytest.raises((ValueError, EOFError, OSError)):
        write_qc_outputs(make(tmp_path), report, plot)
    assert report.read_text() == "real report"
    assert plot.read_bytes() == b"real plot"


def test_write_outputs_creates_report_and_plot(tmp_path):
    fastq = _write(tmp_path, "r.fastq", TWO_READS)
    report, plot = tmp_path / "res" / "report.txt", tmp_path / "res" / "plot.png"
    stats = write_qc_outputs(fastq, report, plot)
    assert report.read_text() == format_qc_report(stats)
    assert plot.read_bytes().startswith(b"\x89PNG")


@pytest.mark.skipif(not REAL_FASTQ.exists(), reason="FASTQ not present")
def test_real_fastq_matches_committed_report():
    assert format_qc_report(compute_fastq_stats(REAL_FASTQ)) == REAL_REPORT
