"""CLI entrypoint, PDF compiler, and Moodle ZIP builder for Lab 04."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from phylogenetics.dataset import (
    DEFAULT_EMAIL,
    GROUPS,
    HANDLE,
    PROJECT_ROOT,
    load_or_fetch_dataset,
)
from phylogenetics.distances import (
    compute_distances,
    format_task1_report,
    write_matrix_csv,
    write_pairs_csv,
)
from phylogenetics.msa import (
    conserved_regions,
    format_task3_report,
    load_or_fetch_clustalo,
    robinson_foulds,
)
from phylogenetics.plot import draw_tree
from phylogenetics.tree import (
    build_nj_tree,
    format_task2_report,
    monophyly,
    write_newick,
)


def _find_pdflatex() -> str:
    found = shutil.which("pdflatex")
    if found:
        return found
    mac_tex = Path("/Library/TeX/texbin/pdflatex")
    if mac_tex.is_file():
        return str(mac_tex)
    raise RuntimeError("pdflatex not found on PATH or /Library/TeX/texbin/pdflatex")


def compile_notes_pdf(tex: Path, out: Path) -> Path:
    """Compile ``tex`` with pdflatex in a temporary directory (run from tex's folder for \\includegraphics)."""
    src = Path(tex).resolve()
    if not src.is_file():
        raise FileNotFoundError(f"LaTeX source not found: {src}")
    with tempfile.TemporaryDirectory() as tmp:
        proc = subprocess.run(
            [
                _find_pdflatex(),
                "-interaction=nonstopmode",
                "-halt-on-error",
                f"-output-directory={tmp}",
                str(src),
            ],
            cwd=src.parent,
            capture_output=True,
            text=True,
            errors="replace",  # pdflatex logs are not always valid UTF-8
            check=False,
        )
        built = Path(tmp) / f"{src.stem}.pdf"
        if proc.returncode != 0 or not built.is_file():
            raise RuntimeError(
                f"pdflatex failed (code {proc.returncode}):\n{proc.stdout[-800:]}"
            )
        dest = Path(out)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built, dest)
    return dest


def build_submission_zip(root: Path, out: Path) -> Path:
    """Create the clean Moodle archive ``phylogenetics_<handle>.zip``."""
    root = Path(root).resolve()
    required = [
        root / name
        for name in ("README.txt", "notes.pdf", "pyproject.toml", ".python-version")
    ]
    for req in required:
        if not req.is_file():
            raise FileNotFoundError(f"Missing required submission file: {req}")
    dest = Path(out).resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    top_files = [f for f in [*required, root / "uv.lock"] if f.is_file()]
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in top_files:
            zf.write(f, arcname=f.relative_to(root).as_posix())
        for d in ("data", "results", "src", "tests"):
            for item in sorted((root / d).rglob("*")):
                if (
                    item.is_dir()
                    or "__pycache__" in item.parts
                    or item.name == ".DS_Store"
                    or item.suffix in {".pyc", ".pyo"}
                ):
                    continue
                zf.write(item, arcname=item.relative_to(root).as_posix())
    return dest


def main(argv: list[str] | None = None) -> int:
    """Run the Lab 04 pipeline; each task recomputes what it needs from the cached data."""
    parser = argparse.ArgumentParser(description="Lab 04 Phylogenetics (Paul Hondola)")
    parser.add_argument(
        "--task", choices=("all", "1", "2", "3", "bonus"), default="all"
    )
    parser.add_argument(
        "--refresh-data",
        action="store_true",
        help="Re-download the TP53 CDS from NCBI (replaces data/tp53_cds.fasta)",
    )
    parser.add_argument(
        "--refresh-msa",
        action="store_true",
        help="Re-run Clustal Omega at EBI (replaces data/clustalo_*)",
    )
    parser.add_argument(
        "--email", default=DEFAULT_EMAIL, help="Contact email for NCBI/EBI refreshes"
    )
    parser.add_argument("--results-dir", type=Path, default=PROJECT_ROOT / "results")
    parser.add_argument(
        "--pdf", action="store_true", help="Compile notes.tex into notes.pdf"
    )
    parser.add_argument(
        "--zip", action="store_true", help=f"Package phylogenetics_{HANDLE}.zip"
    )
    args = parser.parse_args(argv)
    if (args.refresh_data or args.refresh_msa) and (args.pdf or args.zip):
        parser.error(
            "--refresh-* changes results that notes.tex hard-codes; update the notes before --pdf/--zip"
        )

    data, results = PROJECT_ROOT / "data", Path(args.results_dir)
    records = load_or_fetch_dataset(
        data / "tp53_cds.fasta", refresh=args.refresh_data, email=args.email
    )
    print(
        "Dataset: "
        + ", ".join(f"{r.label} ({r.accession}, {len(r.seq)} nt)" for r in records)
    )
    dist = compute_distances(records)
    sections: list[str] = []

    if args.task in {"all", "1"}:
        write_matrix_csv(dist, results / f"distance_matrix_{HANDLE}.csv")
        write_pairs_csv(dist, results / f"pairwise_{HANDLE}.csv")
        sections.append(format_task1_report(dist))

    nj = build_nj_tree(dist.labels, dist.matrix)
    if args.task in {"all", "2"}:
        write_newick(nj.tree, results / f"tree_{HANDLE}.nwk")
        sections.append(format_task2_report(nj, monophyly(nj.tree, GROUPS)))

    if args.task in {"all", "3"}:
        aln, clustal_tree = load_or_fetch_clustalo(
            records,
            data / "clustalo_msa.aln",
            data / "clustalo_phylotree.nwk",
            refresh=args.refresh_msa,
            email=args.email,
        )
        sections.append(
            format_task3_report(
                robinson_foulds(nj.tree, clustal_tree),
                monophyly(nj.tree, GROUPS),
                monophyly(clustal_tree, GROUPS),
                conserved_regions(aln),
                aln,
            )
        )

    if args.task in {"all", "bonus"}:
        png = draw_tree(nj.tree, results / f"tree_{HANDLE}.png", rooting=nj.rooting)
        sections.append(f"=== BONUS: tree plot saved to {png} ===")

    print("\n\n".join(sections))

    if args.pdf:
        pdf = compile_notes_pdf(
            Path(__file__).resolve().parent / "notes.tex", PROJECT_ROOT / "notes.pdf"
        )
        print(f"\nCompiled report: {pdf}")
    if args.zip:
        z = build_submission_zip(
            PROJECT_ROOT, PROJECT_ROOT / f"phylogenetics_{HANDLE}.zip"
        )
        print(f"Created Moodle submission archive: {z}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
