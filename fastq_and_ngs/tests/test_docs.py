from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QC_REPORT = ROOT / "results/qc_report_paulhondola.txt"
VARIANTS = ROOT / "results/variants_paulhondola.txt"


def test_notes_numbers_match_generated_report():
    notes = (ROOT / "src/fastq_and_ngs/notes.tex").read_text(encoding="utf-8")
    report = QC_REPORT.read_text(encoding="utf-8")
    for value in ("81220", "142", "37.30"):
        assert value in report and value in notes


def test_notes_hit_counts_match_generated_report():
    notes = (ROOT / "src/fastq_and_ngs/notes.tex").read_text(encoding="utf-8")
    for line in VARIANTS.read_text(encoding="utf-8").splitlines():
        if line.startswith("Hits: "):
            assert line.removeprefix("Hits: ") in notes
