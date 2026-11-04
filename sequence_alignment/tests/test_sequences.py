from pathlib import Path
import pytest
from Bio.Seq import Seq
from sequence_alignment.sequences import CDS_SPANS, get_record, load_sequences

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_human_cds_span_translates_to_full_length_p53():
    start, end = CDS_SPANS["NM_000546.6"]
    protein = str(Seq(get_record(load_sequences(DATA_DIR), "NM_000546.6").seq[start:end]).translate())
    assert protein.startswith("MEEPQSDP") and protein.endswith("*") and len(protein) == 394
    assert protein[241:256] == "CMGGMNRRPILTIIT"  # residues Cys242..Thr256 cited in notes.tex


def test_load_tp53_sequences_canonical_order_and_lengths():
    recs = load_sequences(DATA_DIR)
    assert [r.accession for r in recs] == ["NM_000546.6", "NM_011640.3", "NM_131327.2"]
    assert [r.length for r in recs] == [2512, 1781, 2233]
    assert get_record(recs, "NM_000546.6").organism == "Homo sapiens"
    assert get_record(recs, "NM_011640.3").organism == "Mus musculus"
    assert get_record(recs, "NM_131327.2").organism == "Danio rerio"


def test_lowercase_whitespace_normalized_and_duplicate_ids_rejected(tmp_path: Path):
    (tmp_path / "a.fa").write_text(
        ">ID1 Homo sapiens test\n acgt \n>ID2 Mus musculus test\nacgt\n>ID3 Danio rerio test\nacgt\n"
    )
    recs = load_sequences(tmp_path)
    assert recs[0].seq == "ACGT"
    (tmp_path / "dup.fa").write_text(">ID1 Duplicate\nACGT\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_sequences(tmp_path)


def test_invalid_inputs_rejected(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_sequences(tmp_path / "missing")
    (tmp_path / "few.fa").write_text(">A Org a\nACGT\n>B Org b\nACGT\n")
    with pytest.raises(ValueError, match="at least 3"):
        load_sequences(tmp_path)
    (tmp_path / "bad.fa").write_text(">C Org c\nACGTZ\n")
    with pytest.raises(ValueError, match="invalid"):
        load_sequences(tmp_path)
