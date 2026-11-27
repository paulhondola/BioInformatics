import zipfile

import pytest

from phylogenetics.__main__ import build_submission_zip, main
from phylogenetics.dataset import PROJECT_ROOT


def test_refresh_with_pdf_is_parser_error():
    with pytest.raises(SystemExit):
        main(["--refresh-msa", "--pdf"])


def test_task2_alone_writes_tree_from_caches(tmp_path):
    assert main(["--task", "2", "--results-dir", str(tmp_path)]) == 0
    assert (tmp_path / "tree_paulhondola.nwk").is_file()
    assert not (tmp_path / "distance_matrix_paulhondola.csv").exists()


def test_task1_writes_both_tables(tmp_path):
    main(["--task", "1", "--results-dir", str(tmp_path)])
    assert (tmp_path / "distance_matrix_paulhondola.csv").is_file()
    assert (tmp_path / "pairwise_paulhondola.csv").is_file()


def test_bonus_png_written(tmp_path):
    main(["--task", "bonus", "--results-dir", str(tmp_path)])
    png = tmp_path / "tree_paulhondola.png"
    assert png.is_file() and png.stat().st_size > 10_000


def test_zip_contents(tmp_path):
    z = build_submission_zip(PROJECT_ROOT, tmp_path / "x.zip")
    names = set(zipfile.ZipFile(z).namelist())
    assert {
        "README.txt",
        "notes.pdf",
        "pyproject.toml",
        "data/tp53_cds.fasta",
        "results/tree_paulhondola.nwk",
        "src/phylogenetics/__main__.py",
    } <= names
    assert not any("__pycache__" in n or n.endswith(".DS_Store") for n in names)


def test_plot_title_names_the_rooting(tmp_path, monkeypatch):
    from io import StringIO

    from Bio import Phylo

    from phylogenetics import plot

    titles = []
    monkeypatch.setattr(
        plot.plt, "close", lambda fig: titles.append(fig.axes[0].get_title())
    )
    tree = Phylo.read(
        StringIO("((Homo_sapiens:1,Mus_musculus:1):1,Danio_rerio:1);"), "newick"
    )
    plot.draw_tree(tree, tmp_path / "t.png", rooting="midpoint")
    assert "midpoint" in titles[0] and "fish" not in titles[0]
