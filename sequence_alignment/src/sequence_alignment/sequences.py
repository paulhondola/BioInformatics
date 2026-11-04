"""FASTA sequence loader and metadata validation for Lab 02."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

CANONICAL_ORDER: tuple[str, ...] = ("NM_000546.6", "NM_011640.3", "NM_131327.2")
ALLOWED_IUPAC_DNA: frozenset[str] = frozenset("ACGTRYSWKMBDHVN")
# 0-based half-open CDS including the stop codon (RefSeq CDS 143..1324).
CDS_SPANS: dict[str, tuple[int, int]] = {"NM_000546.6": (142, 1324)}


@dataclass(frozen=True)
class SequenceRecord:
    """Validated nucleotide sequence record."""

    accession: str
    organism: str
    description: str
    seq: str
    length: int
    path: Path


def _extract_organism(accession: str, description: str) -> str:
    remainder = (
        description[len(accession) :].strip()
        if description.startswith(accession)
        else description.strip()
    )
    tokens = remainder.split()
    if len(tokens) >= 2:
        return f"{tokens[0]} {tokens[1]}"
    if len(tokens) == 1:
        return tokens[0]
    return "Unknown"


def load_sequences(data_dir: Path) -> list[SequenceRecord]:
    """Load, normalize, validate, and canonically order FASTA records from ``data_dir``."""
    data_path = Path(data_dir)
    if not data_path.is_dir():
        raise FileNotFoundError(f"Data directory does not exist: {data_path}")

    fasta_files = sorted({*data_path.glob("*.fa"), *data_path.glob("*.fasta")})
    records: list[SequenceRecord] = []
    seen_ids: set[str] = set()

    for fasta_file in fasta_files:
        with fasta_file.open(encoding="utf-8") as handle:
            for bio_rec in SeqIO.parse(handle, "fasta"):
                acc = bio_rec.id.strip()
                if acc in seen_ids:
                    raise ValueError(
                        f"Found duplicate accession '{acc}' in {fasta_file.name}"
                    )
                seen_ids.add(acc)

                clean_seq = "".join(str(bio_rec.seq).split()).upper()
                if not clean_seq:
                    raise ValueError(f"Sequence '{acc}' in {fasta_file.name} is empty")

                invalid_chars = sorted(set(clean_seq) - ALLOWED_IUPAC_DNA)
                if invalid_chars:
                    raise ValueError(
                        f"Sequence '{acc}' contains invalid nucleotide characters: {''.join(invalid_chars)}"
                    )

                records.append(
                    SequenceRecord(
                        accession=acc,
                        organism=_extract_organism(acc, bio_rec.description),
                        description=bio_rec.description.strip(),
                        seq=clean_seq,
                        length=len(clean_seq),
                        path=fasta_file,
                    )
                )

    if len(records) < 3:
        raise ValueError(
            f"Expected at least 3 sequences in {data_path}, found {len(records)}"
        )

    order_map = {acc: idx for idx, acc in enumerate(CANONICAL_ORDER)}
    records.sort(
        key=lambda r: (order_map.get(r.accession, len(CANONICAL_ORDER)), r.accession)
    )
    return records


def get_record(records: list[SequenceRecord], accession: str) -> SequenceRecord:
    """Return the record matching ``accession`` or raise ``KeyError``."""
    for rec in records:
        if rec.accession == accession:
            return rec
    raise KeyError(f"Sequence accession '{accession}' not found")
