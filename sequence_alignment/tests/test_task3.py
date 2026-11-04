from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch
import pytest
from Bio import AlignIO
from sequence_alignment.sequences import SequenceRecord, load_sequences
from sequence_alignment.task3_msa import (
    EBI_MUSCLE_BASE,
    analyze_msa,
    fetch_clustalo_msa,
    format_task3_report,
    load_or_fetch_msa,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

TOY_CLUSTAL = (
    "CLUSTAL O(1.2.4) multiple sequence alignment\n\n"
    "NM_000546.6      TTTACGTACGTACGTGGG\t18\n"
    "NM_011640.3      CCCACGTACGTACGTAAA\t18\n"
    "NM_131327.2      GGGACGTACGTACGTTTT\t18\n"
    "                    ************   \n"
)


def _toy_records() -> list[SequenceRecord]:
    return [
        SequenceRecord("NM_000546.6", "Homo sapiens", "", "TTTACGTACGTACGTGGG", 18, Path("a")),
        SequenceRecord("NM_011640.3", "Mus musculus", "", "CCCACGTACGTACGTAAA", 18, Path("b")),
        SequenceRecord("NM_131327.2", "Danio rerio", "", "GGGACGTACGTACGTTTT", 18, Path("c")),
    ]


def test_fetch_clustalo_msa_polls_until_finished_and_handles_failure():
    post_resp = MagicMock(status_code=200, text="job-123")
    running_resp = MagicMock(status_code=200, text="RUNNING")
    done_resp = MagicMock(status_code=200, text="FINISHED")
    aln_resp = MagicMock(status_code=200, text=TOY_CLUSTAL)
    with patch("sequence_alignment.task3_msa.requests.post", return_value=post_resp) as m_post, patch(
        "sequence_alignment.task3_msa.requests.get", side_effect=[running_resp, done_resp, aln_resp]
    ):
        text = fetch_clustalo_msa(_toy_records(), poll_interval_s=0.01, timeout_s=5.0)
        assert "CLUSTAL" in text
        assert m_post.call_args.kwargs["data"]["stype"] == "dna"

    err_resp = MagicMock(status_code=200, text="FAILURE")
    with patch("sequence_alignment.task3_msa.requests.post", return_value=post_resp), patch(
        "sequence_alignment.task3_msa.requests.get", return_value=err_resp
    ):
        with pytest.raises(RuntimeError, match="FAILURE"):
            fetch_clustalo_msa(_toy_records(), poll_interval_s=0.01, timeout_s=1.0)


def test_cache_used_offline_and_mismatched_ids_rejected(tmp_path: Path):
    cache = tmp_path / "clustalo_msa.aln"
    cache.write_text(TOY_CLUSTAL)
    with patch("sequence_alignment.task3_msa.requests.post") as m_post:
        msa = load_or_fetch_msa(_toy_records(), cache, refresh=False)
        m_post.assert_not_called()
        assert len(msa) == 3
    bad_recs = _toy_records()[:2] + [SequenceRecord("OTHER", "Org", "", "ACGT", 4, Path("d"))]
    with pytest.raises(ValueError, match="do not match"):
        load_or_fetch_msa(bad_recs, cache, refresh=False)


def test_refresh_keeps_cache_when_response_is_invalid(tmp_path: Path):
    cache = tmp_path / "clustalo_msa.aln"
    cache.write_text(TOY_CLUSTAL)
    wrong_seq = TOY_CLUSTAL.replace("TTTACGTACGTACGTGGG", "AAAACGTACGTACGTGGG")
    for bad_text in ("<html>service error</html>", wrong_seq):
        with patch("sequence_alignment.task3_msa.fetch_clustalo_msa", return_value=bad_text):
            with pytest.raises(ValueError):
                load_or_fetch_msa(_toy_records(), cache, refresh=True)
        assert cache.read_text() == TOY_CLUSTAL


def test_refresh_cannot_be_packaged_with_static_notes():
    with pytest.raises(SystemExit):
        main(["--refresh-msa", "--zip"])


def test_analyze_msa_extracts_exact_conserved_motif():
    msa = AlignIO.read(StringIO(TOY_CLUSTAL), "clustal")
    res = analyze_msa(msa, window=12)
    assert res.fully_conserved_cols == 12
    assert res.top_contiguous_motif.consensus_seq == "ACGTACGTACGT"
    assert res.top_contiguous_motif.per_seq_coords["NM_000546.6"] == (4, 15)
    rendered = format_task3_report(res)
    assert "ACGTACGTACGT" in rendered and "*" in rendered


def test_fetch_clustalo_msa_falls_back_to_ebi_muscle_on_queued_stall():
    co_post = MagicMock(status_code=200, text="clustalo-job-1")
    mu_post = MagicMock(status_code=200, text="muscle-job-2")
    queued_resp = MagicMock(status_code=200, text="QUEUED")
    done_resp = MagicMock(status_code=200, text="FINISHED")
    aln_resp = MagicMock(status_code=200, text=TOY_CLUSTAL)
    with patch("sequence_alignment.task3_msa.requests.post", side_effect=[co_post, mu_post]) as m_post, patch(
        "sequence_alignment.task3_msa.requests.get",
        side_effect=[queued_resp, queued_resp, done_resp, aln_resp],
    ):
        text = fetch_clustalo_msa(
            _toy_records(),
            poll_interval_s=0.01,
            timeout_s=1.0,
            clustalo_queue_timeout_s=0.01,
        )
        assert "CLUSTAL" in text
        assert m_post.call_args_list[1].args[0] == f"{EBI_MUSCLE_BASE}/run"


def test_real_tp53_cached_msa():
    recs = load_sequences(DATA_DIR)
    msa = load_or_fetch_msa(recs, DATA_DIR / "clustalo_msa.aln", refresh=False)
    res = analyze_msa(msa)
    assert res.alignment_length > 2500
    assert res.top_contiguous_motif.length >= 9
    assert res.best_window_motif.identity_pct > 80.0


import os
import re
import subprocess
import sys
import zipfile
from sequence_alignment.__main__ import build_submission_zip, main

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_cli_runs_and_deliverables_meet_moodle_constraints(tmp_path: Path, capsys):
    assert main(["--data-dir", str(REPO_ROOT / "data")]) == 0
    out = capsys.readouterr().out
    assert "0.1641" in out and "2014" in out and "NM_131327.2" in out and "*" in out

    readme_lines = (REPO_ROOT / "README.txt").read_text().strip().splitlines()
    assert 1 <= len(readme_lines) <= 10

    pdf_bytes = (REPO_ROOT / "notes.pdf").read_bytes()
    assert pdf_bytes.startswith(b"%PDF-")
    assert len(re.findall(rb"/Type\s*/Page\b", pdf_bytes)) == 1

    with pytest.raises(FileNotFoundError):
        build_submission_zip(tmp_path, tmp_path / "missing.zip")

    zip_path = build_submission_zip(REPO_ROOT, tmp_path / "sequence_alignment_paulhondola.zip")
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        assert {"README.txt", "notes.pdf", "pyproject.toml", "data/clustalo_msa.aln"} <= names
        assert not any(
            "__pycache__" in n or n.startswith("docs/") or ".DS_Store" in n for n in names
        )
        extract_dir = tmp_path / "unzipped"
        zf.extractall(extract_dir)
    proc = subprocess.run(
        [sys.executable, "-W", "error", "-m", "sequence_alignment"],
        cwd=extract_dir,
        env={**os.environ, "PYTHONPATH": str(extract_dir / "src")},
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stderr == ""
    assert "0.1641" in proc.stdout and "2014" in proc.stdout
