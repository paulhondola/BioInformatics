"""CLI entrypoint, FASTA copy, PDF compiler, and Moodle ZIP builder for Lab 03."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from fastq_and_ngs.fastq_stats import write_qc_outputs
from fastq_and_ngs.pubmed_query import (
    DEFAULT_EMAIL,
    HANDLE,
    PROJECT_ROOT,
    format_pubmed_report,
    load_or_fetch_articles,
)
from fastq_and_ngs.vcf_pubmed import (
    format_variants_report,
    load_or_fetch_hits,
    parse_vcf,
)

FASTQ_RELATIVE: str = "data/work/reads.fastq.gz"
FASTQ_SOURCE: str = (
    "SRR31709504 read 1 (https://www.ebi.ac.uk/ena/browser/view/SRR31709504)"
)
OUTPUTS: tuple[str, ...] = (
    f"pubmed_{HANDLE}.txt",
    f"fasta_{HANDLE}.fasta",
    f"qc_report_{HANDLE}.txt",
    f"qc_plot_{HANDLE}.png",
    f"variants_{HANDLE}.txt",
)


def copy_fasta(source: Path, dest: Path) -> Path:
    """Task 2: save a byte-identical copy of the Lab 1 FASTA."""
    target = Path(dest)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    return target


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _run_task(root: Path, task: str, refresh: bool, email: str) -> None:
    results = root / "results"
    if task == "1":
        articles = load_or_fetch_articles(
            root / "data/pubmed_tp53_cancer.xml", refresh=refresh, email=email
        )
        _write_text(results / f"pubmed_{HANDLE}.txt", format_pubmed_report(articles))
    elif task == "2":
        copy_fasta(root / "data/human_tp53.fa", results / f"fasta_{HANDLE}.fasta")
    elif task == "3":
        write_qc_outputs(
            root / FASTQ_RELATIVE,
            results / f"qc_report_{HANDLE}.txt",
            results / f"qc_plot_{HANDLE}.png",
        )
    elif task == "4":
        hits = load_or_fetch_hits(
            parse_vcf(root / "data/tp53_variants.vcf"),
            root / "data/variant_searches.json",
            refresh=refresh,
            email=email,
        )
        _write_text(results / f"variants_{HANDLE}.txt", format_variants_report(hits))


def run_tasks(root: Path, task: str, refresh: bool, email: str) -> int:
    """Run one task or all of them; returns a process exit code."""
    code = 0
    for current in ("1", "2", "3", "4") if task == "all" else (task,):
        if current == "3" and not (root / FASTQ_RELATIVE).exists():
            level = "Warning, skipping Task 3" if task == "all" else "Error"
            print(
                f"{level}: {FASTQ_RELATIVE} is missing. Source: {FASTQ_SOURCE}",
                file=sys.stderr,
            )
            if task != "all":
                code = 1
            continue
        try:
            _run_task(root, current, refresh, email)
            print(f"[OK] Task {current}")
        except (OSError, EOFError, ValueError) as exc:
            detail = f" (reading {FASTQ_RELATIVE})" if current == "3" else ""
            print(f"Error in Task {current}{detail}: {exc}", file=sys.stderr)
            code = 1
    return code


def _find_pdflatex() -> str:
    found = shutil.which("pdflatex")
    if found:
        return found
    mac_tex = Path("/Library/TeX/texbin/pdflatex")
    if mac_tex.is_file():
        return str(mac_tex)
    raise RuntimeError(
        "pdflatex executable not found on PATH or /Library/TeX/texbin/pdflatex"
    )


def compile_notes_pdf(tex_path: Path, output_pdf: Path) -> Path:
    """Compile ``tex_path`` via pdflatex in an isolated temporary directory."""
    src_tex = Path(tex_path).resolve()
    if not src_tex.is_file():
        raise FileNotFoundError(f"LaTeX source not found: {src_tex}")
    pdflatex = _find_pdflatex()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        proc = subprocess.run(
            [
                pdflatex,
                "-interaction=nonstopmode",
                "-halt-on-error",
                f"-output-directory={tmp_path}",
                str(src_tex),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        built_pdf = tmp_path / f"{src_tex.stem}.pdf"
        if proc.returncode != 0 or not built_pdf.is_file():
            raise RuntimeError(
                f"pdflatex failed (code {proc.returncode}):\n{proc.stdout[-800:]}\n{proc.stderr[-400:]}"
            )
        dest_pdf = Path(output_pdf)
        dest_pdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built_pdf, dest_pdf)
    return dest_pdf


def build_submission_zip(root: Path, output_zip: Path) -> Path:
    """Create the Moodle archive; the FASTQ stays out because it is not a listed deliverable."""
    base = Path(root).resolve()
    required = [
        base / "README.txt",
        base / "notes.pdf",
        base / "pyproject.toml",
        base / ".python-version",
        *(base / "results" / name for name in OUTPUTS),
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f"Missing required submission file: {path}")

    zip_dest = Path(output_zip).resolve()
    zip_dest.parent.mkdir(parents=True, exist_ok=True)
    top_files = [f for f in [*required[:4], base / "uv.lock"] if f.is_file()]

    with zipfile.ZipFile(zip_dest, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file_path in top_files:
            zf.write(file_path, arcname=file_path.relative_to(base).as_posix())
        for dir_path in (base / "data", base / "results", base / "src", base / "tests"):
            for item in sorted(dir_path.rglob("*")):
                relative = item.relative_to(base)
                if (
                    item.is_dir()
                    or relative.as_posix().startswith("data/work/")
                    or "__pycache__" in relative.parts
                    or item.name == ".DS_Store"
                    or item.suffix in {".pyc", ".pyo"}
                ):
                    continue
                zf.write(item, arcname=relative.as_posix())
    return zip_dest


def main(argv: list[str] | None = None) -> int:
    """Run the Lab 03 pipeline and optional PDF/ZIP builders."""
    parser = argparse.ArgumentParser(
        description="Lab 03 Formats and NGS Pipeline (Paul Hondola)"
    )
    parser.add_argument(
        "--task",
        choices=("all", "1", "2", "3", "4"),
        default="all",
        help="Select one task or all",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Re-query PubMed instead of using the cache in data/",
    )
    parser.add_argument(
        "--email",
        default=DEFAULT_EMAIL,
        help="Email sent to NCBI Entrez with --refresh",
    )
    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Compile src/fastq_and_ngs/notes.tex into notes.pdf",
    )
    parser.add_argument(
        "--zip", action="store_true", help=f"Package fastq_and_ngs_{HANDLE}.zip"
    )
    args = parser.parse_args(argv)
    if args.refresh and (args.pdf or args.zip):
        parser.error(
            "--refresh changes results that notes.tex hard-codes; update the notes before --pdf/--zip"
        )

    code = run_tasks(PROJECT_ROOT, args.task, args.refresh, args.email)
    if code != 0:
        return code

    if args.pdf:
        pdf_out = compile_notes_pdf(
            Path(__file__).resolve().parent / "notes.tex", PROJECT_ROOT / "notes.pdf"
        )
        print(f"Compiled report: {pdf_out}")

    if args.zip:
        zip_out = build_submission_zip(
            PROJECT_ROOT, PROJECT_ROOT / f"fastq_and_ngs_{HANDLE}.zip"
        )
        print(f"Created Moodle submission archive: {zip_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
