"""Exercise 04: FASTQ QC on own reads (report and bonus plot).

Standalone on purpose: imports nothing from the package, so this file also runs on its own.
Pass --input, --report and --plot when running it outside this project's folders.
"""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

import matplotlib
from Bio import SeqIO

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HANDLE: str = "paulhondola"
PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class FastqStats:
    """Aggregate QC metrics plus the per-read values behind the plot."""

    reads: int
    mean_length: float
    n_rate: float
    mean_phred: float
    read_lengths: list[int]
    read_mean_phreds: list[float]


def compute_fastq_stats(path: Path) -> FastqStats:
    """Read a plain or gzipped FASTQ and compute QC metrics."""
    fastq = Path(path)
    opener = gzip.open if fastq.suffix == ".gz" else open
    total_n = 0
    total_phred = 0
    total_bases = 0
    read_lengths: list[int] = []
    read_mean_phreds: list[float] = []

    with opener(fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            phred = record.letter_annotations["phred_quality"]
            total_n += str(record.seq).upper().count("N")
            total_phred += sum(phred)
            total_bases += len(phred)
            read_lengths.append(len(record.seq))
            read_mean_phreds.append(sum(phred) / len(phred) if phred else 0.0)

    reads = len(read_lengths)
    return FastqStats(
        reads=reads,
        mean_length=total_bases / reads if reads else 0.0,
        n_rate=total_n / total_bases if total_bases else 0.0,
        mean_phred=total_phred / total_bases if total_bases else 0.0,
        read_lengths=read_lengths,
        read_mean_phreds=read_mean_phreds,
    )


def format_qc_report(stats: FastqStats) -> str:
    """Render the course template's four report lines."""
    return (
        f"Reads: {stats.reads}\n"
        f"Mean length: {stats.mean_length:.2f}\n"
        f"N rate: {stats.n_rate:.4f}\n"
        f"Mean Phred: {stats.mean_phred:.2f}\n"
    )


def write_qc_plot(stats: FastqStats, out_png: Path) -> Path:
    """Bonus: save read-length and per-read mean Phred distributions."""
    out = Path(out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.hist(stats.read_lengths, bins=25, color="#2b5c8f", edgecolor="black", alpha=0.8)
    ax1.set_title("Distribuția lungimilor citirilor")
    ax1.set_xlabel("Lungime (bp)")
    ax1.set_ylabel("Număr citiri")
    ax1.grid(True, linestyle="--", alpha=0.5)

    ax2.hist(
        stats.read_mean_phreds, bins=25, color="#2e7d32", edgecolor="black", alpha=0.8
    )
    ax2.set_title("Distribuția scorurilor Phred medii")
    ax2.set_xlabel("Scor Phred mediu per citire")
    ax2.set_ylabel("Număr citiri")
    ax2.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close(fig)
    return out


def write_qc_outputs(fastq: Path, report: Path, plot: Path) -> FastqStats:
    """Compute stats first, so an unreadable FASTQ can never overwrite existing outputs."""
    stats = compute_fastq_stats(fastq)
    if stats.reads == 0:
        raise ValueError(f"{fastq} holds no reads; existing QC outputs were kept")
    report_path = Path(report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(format_qc_report(stats), encoding="utf-8")
    write_qc_plot(stats, plot)
    return stats


def main(argv: list[str] | None = None) -> int:
    """Write the QC report and plot for one FASTQ file."""
    parser = argparse.ArgumentParser(description="FASTQ QC report and plot")
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT / f"data/work/{HANDLE}/lab03/your_reads.fastq.gz",
    )
    parser.add_argument(
        "--report", type=Path, default=PROJECT_ROOT / f"results/qc_report_{HANDLE}.txt"
    )
    parser.add_argument(
        "--plot", type=Path, default=PROJECT_ROOT / f"results/qc_plot_{HANDLE}.png"
    )
    args = parser.parse_args(argv)

    stats = write_qc_outputs(args.input, args.report, args.plot)
    print(format_qc_report(stats), end="")
    print(f"[OK] QC report -> {args.report}\n[OK] QC plot -> {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
