import zipfile
from pathlib import Path

import pytest

from fastq_and_ngs import __main__ as cli
from fastq_and_ngs.__main__ import build_submission_zip, copy_fasta, main, run_tasks
from fastq_and_ngs.pubmed_query import Article

EMAIL = "test@example.org"
OUTPUTS = (
    "pubmed_paulhondola.txt",
    "fasta_paulhondola.fasta",
    "qc_report_paulhondola.txt",
    "qc_plot_paulhondola.png",
    "variants_paulhondola.txt",
)
FASTQ = "data/work/reads.fastq.gz"


@pytest.fixture
def project(tmp_path, monkeypatch):
    """A fake project: real parsers, no network and no real FASTQ."""
    (tmp_path / "data").mkdir()
    (tmp_path / "data/human_tp53.fa").write_text(">NM_000546.6\nACGT\n")
    (tmp_path / "data/tp53_variants.vcf").write_text(
        "17\t1\trs1\tA\tT\n17\t2\t.\tA\tT\n"
    )
    monkeypatch.setattr(
        cli, "load_or_fetch_articles", lambda *a, **k: [Article("1", "T", (), None)]
    )
    monkeypatch.setattr(cli, "load_or_fetch_hits", lambda *a, **k: [])
    return tmp_path


def _existing_report(root: Path) -> Path:
    report = root / "results/qc_report_paulhondola.txt"
    report.parent.mkdir(exist_ok=True)
    report.write_text("real report")
    return report


def _complete_project(root: Path) -> None:
    for name in (
        "README.txt",
        "notes.pdf",
        "pyproject.toml",
        ".python-version",
        "uv.lock",
    ):
        (root / name).write_text("x")
    (root / "results").mkdir(exist_ok=True)
    for name in OUTPUTS:
        (root / "results" / name).write_text("x")
    for rel in (
        "src/fastq_and_ngs/a.py",
        "tests/test_a.py",
        FASTQ,
        "src/fastq_and_ngs/__pycache__/a.pyc",
        "data/.DS_Store",
    ):
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x")


def test_fasta_copy_is_byte_identical(tmp_path):
    source = tmp_path / "in.fa"
    source.write_bytes(b">x\r\nACGT\n\n")
    dest = copy_fasta(source, tmp_path / "res" / "out.fasta")
    assert dest.read_bytes() == source.read_bytes()


def test_task_all_without_fastq_warns_and_runs_others(project, capsys):
    report = _existing_report(project)
    assert run_tasks(project, "all", False, EMAIL) == 0
    assert "SRR31709504" in capsys.readouterr().err
    assert report.read_text() == "real report"
    for name in (
        "pubmed_paulhondola.txt",
        "fasta_paulhondola.fasta",
        "variants_paulhondola.txt",
    ):
        assert (project / "results" / name).exists()


def test_task_3_without_fastq_fails_clearly(project, capsys):
    report = _existing_report(project)
    assert run_tasks(project, "3", False, EMAIL) == 1
    err = capsys.readouterr().err
    assert "SRR31709504" in err and "ena/browser/view/SRR31709504" in err
    assert report.read_text() == "real report"


def test_corrupt_fastq_returns_error_and_keeps_report(project, capsys):
    report = _existing_report(project)
    fastq = project / FASTQ
    fastq.parent.mkdir(parents=True)
    fastq.write_bytes(b"not gzip")
    assert run_tasks(project, "3", False, EMAIL) == 1
    assert "reads" in capsys.readouterr().err
    assert report.read_text() == "real report"


def test_single_task_runs_only_that_task(project):
    assert run_tasks(project, "2", False, EMAIL) == 0
    assert [p.name for p in (project / "results").iterdir()] == [
        "fasta_paulhondola.fasta"
    ]


def test_refresh_failure_returns_error_without_traceback(project, monkeypatch, capsys):
    def down(*args, **kwargs):
        raise OSError("network down")

    monkeypatch.setattr(cli, "load_or_fetch_articles", down)
    assert run_tasks(project, "1", True, EMAIL) == 1
    assert "network down" in capsys.readouterr().err


@pytest.mark.parametrize("flag", ["--pdf", "--zip"])
def test_refresh_cannot_combine_with_pdf_or_zip(flag):
    with pytest.raises(SystemExit):
        main(["--refresh", flag])


def test_zip_layout_and_exclusions(project, tmp_path):
    _complete_project(project)
    out = project / "fastq_and_ngs_paulhondola.zip"
    names = zipfile.ZipFile(build_submission_zip(project, out)).namelist()
    for needed in (
        "README.txt",
        "notes.pdf",
        "uv.lock",
        "data/human_tp53.fa",
        "src/fastq_and_ngs/a.py",
        "tests/test_a.py",
    ):
        assert needed in names
    assert {f"results/{name}" for name in OUTPUTS} <= set(names)
    assert not any(
        n.startswith("data/work/")
        or "__pycache__" in n
        or n.endswith((".DS_Store", ".zip"))
        for n in names
    )


def test_zip_missing_output_names_file_and_leaves_no_zip(project):
    _complete_project(project)
    (project / "results/qc_plot_paulhondola.png").unlink()
    out = project / "submission.zip"
    with pytest.raises(FileNotFoundError, match="qc_plot_paulhondola.png"):
        build_submission_zip(project, out)
    assert not out.exists()
