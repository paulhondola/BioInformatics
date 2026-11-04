"""CLI entrypoint, PDF compiler, and Moodle ZIP builder for Lab 02."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

from sequence_alignment.bonus_semiglobal import compare_semiglobal, format_bonus_report
from sequence_alignment.sequences import get_record, load_sequences
from sequence_alignment.task1_distances import compute_distances, format_task1_report
from sequence_alignment.task2_pairwise import compare_global_local, format_task2_report
from sequence_alignment.task3_msa import analyze_msa, format_task3_report, load_or_fetch_msa


def _default_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _find_pdflatex() -> str:
    found = shutil.which("pdflatex")
    if found:
        return found
    mac_tex = Path("/Library/TeX/texbin/pdflatex")
    if mac_tex.is_file():
        return str(mac_tex)
    raise RuntimeError("pdflatex executable not found on PATH or /Library/TeX/texbin/pdflatex")


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


def build_submission_zip(repo_root: Path, output_zip: Path) -> Path:
    """Create a clean Moodle submission archive ``sequence_alignment_paulhondola.zip``."""
    root = Path(repo_root).resolve()
    required_files = [
        root / "README.txt",
        root / "notes.pdf",
        root / "pyproject.toml",
        root / ".python-version",
    ]
    for req in required_files:
        if not req.is_file():
            raise FileNotFoundError(f"Missing required submission file: {req}")

    zip_dest = Path(output_zip).resolve()
    zip_dest.parent.mkdir(parents=True, exist_ok=True)

    top_files = [f for f in [*required_files, root / "uv.lock"] if f.is_file()]
    top_dirs = [root / "data", root / "src", root / "tests"]

    with zipfile.ZipFile(zip_dest, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file_path in top_files:
            zf.write(file_path, arcname=file_path.relative_to(root).as_posix())
        for dir_path in top_dirs:
            if not dir_path.is_dir():
                continue
            for item in sorted(dir_path.rglob("*")):
                if item.is_dir():
                    continue
                parts = set(item.parts)
                if "__pycache__" in parts or item.name == ".DS_Store" or item.suffix in {".pyc", ".pyo"}:
                    continue
                zf.write(item, arcname=item.relative_to(root).as_posix())
    return zip_dest


def main(argv: list[str] | None = None) -> int:
    """Run Lab 02 sequence alignment pipeline and optional PDF/ZIP builders."""
    repo_root = _default_repo_root()
    parser = argparse.ArgumentParser(description="Lab 02 Sequence Alignment Pipeline (Paul Hondola)")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=repo_root / "data",
        help="Directory containing Lab 1 FASTA files and cached MSA",
    )
    parser.add_argument(
        "--task",
        choices=("all", "1", "2", "3", "bonus"),
        default="all",
        help="Select a single task or run all tasks",
    )
    parser.add_argument(
        "--refresh-msa",
        action="store_true",
        help="Re-submit sequences to EBI online MSA instead of using data/clustalo_msa.aln",
    )
    parser.add_argument(
        "--email",
        default="paulhondola@gmail.com",
        help="Email passed to EBI Tools REST API when --refresh-msa is used",
    )
    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Compile src/sequence_alignment/notes.tex into root notes.pdf using pdflatex",
    )
    parser.add_argument(
        "--zip",
        action="store_true",
        help="Package clean sequence_alignment_paulhondola.zip at repository root",
    )
    args = parser.parse_args(argv)
    if args.refresh_msa and (args.pdf or args.zip):
        parser.error("--refresh-msa changes the MSA that notes.tex hard-codes; update the notes before --pdf/--zip")

    records = load_sequences(args.data_dir)
    hs = get_record(records, "NM_000546.6")
    mm = get_record(records, "NM_011640.3")

    sections: list[str] = [
        "Loaded Lab 1 TP53 RefSeq transcripts:\n"
        + "\n".join(f"  - {r.accession} ({r.organism}, {r.length} nt) [{r.path.name}]" for r in records)
    ]

    if args.task in {"all", "1"}:
        sections.append(format_task1_report(compute_distances(records)))

    if args.task in {"all", "2"}:
        sections.append(format_task2_report(compare_global_local(hs, mm)))

    if args.task in {"all", "3"}:
        msa_cache = Path(args.data_dir) / "clustalo_msa.aln"
        msa = load_or_fetch_msa(records, msa_cache, refresh=args.refresh_msa, email=args.email)
        sections.append(format_task3_report(analyze_msa(msa)))

    if args.task in {"all", "bonus"}:
        sections.append(format_bonus_report(compare_semiglobal(hs, mm)))

    print("\n\n".join(sections))

    if args.pdf:
        pdf_out = compile_notes_pdf(Path(__file__).resolve().parent / "notes.tex", repo_root / "notes.pdf")
        print(f"\nCompiled report: {pdf_out}")

    if args.zip:
        zip_out = build_submission_zip(repo_root, repo_root / "sequence_alignment_paulhondola.zip")
        print(f"Created Moodle submission archive: {zip_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
