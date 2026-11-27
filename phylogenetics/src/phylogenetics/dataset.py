"""TP53 coding sequences (CDS) for 12 vertebrates: species table, NCBI fetch, validation and cache."""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO

HANDLE: str = "paulhondola"
DEFAULT_EMAIL: str = "paulhondola@gmail.com"
PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]

STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})
ALLOWED: frozenset[str] = frozenset("ACGTRYSWKMBDHVN")  # IUPAC DNA; medaka has a Y


@dataclass(frozen=True)
class Species:
    """One dataset entry: Newick-safe label, RefSeq mRNA accession, expected clade."""

    label: str
    accession: str
    group: str


SPECIES: tuple[Species, ...] = (
    Species("Homo_sapiens", "NM_000546.6", "Primates"),
    Species("Pan_troglodytes", "XM_001172077.6", "Primates"),  # no curated NM_
    Species("Macaca_mulatta", "NM_001047151.2", "Primates"),
    Species("Mus_musculus", "NM_011640.3", "Rodents"),
    Species("Rattus_norvegicus", "NM_030989.4", "Rodents"),  # transcript variant 1
    Species("Canis_lupus_familiaris", "NM_001389218.1", "Laurasiatheria"),
    Species("Bos_taurus", "NM_174201.2", "Laurasiatheria"),
    Species("Sus_scrofa", "NM_213824.3", "Laurasiatheria"),
    Species("Gallus_gallus", "NM_205264.1", "Bird"),
    Species("Xenopus_tropicalis", "NM_001001903.1", "Amphibian"),
    Species("Danio_rerio", "NM_131327.2", "Fish"),
    Species("Oryzias_latipes", "NM_001104742.1", "Fish"),
)


def _group(*names: str) -> frozenset[str]:
    return frozenset(s.label for s in SPECIES if s.group in names)


GROUPS: dict[str, frozenset[str]] = {
    **{g: _group(g) for g in dict.fromkeys(s.group for s in SPECIES)},
    "Mammals": _group("Primates", "Rodents", "Laurasiatheria"),
    "Tetrapods": _group("Primates", "Rodents", "Laurasiatheria", "Bird", "Amphibian"),
}


@dataclass(frozen=True)
class CdsRecord:
    """One validated coding sequence."""

    label: str
    accession: str
    seq: str


def validate_cds(label: str, seq: str) -> None:
    """Raise ``ValueError`` unless ``seq`` is a complete CDS: ATG ... stop, whole codons, IUPAC DNA only."""
    if not seq:
        raise ValueError(f"{label}: empty CDS")
    bad = sorted(set(seq) - ALLOWED)
    if bad:
        raise ValueError(f"{label}: invalid characters {''.join(bad)}")
    if not seq.startswith("ATG"):
        raise ValueError(f"{label}: CDS does not start with ATG")
    if len(seq) % 3:
        raise ValueError(f"{label}: CDS length {len(seq)} is not a multiple of 3")
    if seq[-3:] not in STOP_CODONS:
        raise ValueError(f"{label}: CDS does not end in a stop codon")


def parse_cds_fasta(text: str) -> list[CdsRecord]:
    """Parse the cache format (``>Label accession``); require exactly the 12 species, in SPECIES order."""
    by_label: dict[str, CdsRecord] = {}
    for rec in SeqIO.parse(StringIO(text), "fasta"):
        parts = rec.description.split()
        if len(parts) != 2:
            raise ValueError(f"Bad header (want '>Label accession'): {rec.description}")
        label, accession = parts
        if label in by_label:
            raise ValueError(f"Duplicate label {label}")
        seq = str(rec.seq).upper()
        validate_cds(label, seq)
        by_label[label] = CdsRecord(label, accession, seq)
    expected = [s.label for s in SPECIES]
    if set(by_label) != set(expected):
        raise ValueError(
            f"Dataset labels {sorted(by_label)} do not match SPECIES {sorted(expected)}"
        )
    return [by_label[label] for label in expected]


def fetch_cds_fasta(email: str) -> str:
    """Fetch each species' CDS from NCBI and return it in cache format."""
    Entrez.email = email
    accessions = [s.accession for s in SPECIES]
    with Entrez.efetch(
        db="nuccore", id=",".join(accessions), rettype="fasta_cds_na", retmode="text"
    ) as handle:
        reply = handle.read()
    found: dict[str, list[str]] = {acc: [] for acc in accessions}
    for rec in SeqIO.parse(StringIO(reply), "fasta"):
        # ids look like lcl|NM_000546.6_cds_NP_000537.3_1; match by accession, never by reply order.
        for acc in accessions:
            if f"|{acc}_cds" in rec.id:
                found[acc].append(str(rec.seq).upper())
    lines: list[str] = []
    for s in SPECIES:
        seqs = found[s.accession]
        if len(seqs) != 1:
            raise ValueError(
                f"{s.accession} ({s.label}): expected 1 CDS, got {len(seqs)}"
            )
        lines += [f">{s.label} {s.accession}", seqs[0]]
    return "\n".join(lines) + "\n"


def load_or_fetch_dataset(
    cache: Path, refresh: bool = False, email: str = DEFAULT_EMAIL
) -> list[CdsRecord]:
    """Read the cached dataset, or fetch it; a fetched reply is validated before it replaces the cache."""
    cache = Path(cache)
    if refresh or not cache.exists():
        text = fetch_cds_fasta(email)
        records = parse_cds_fasta(text)
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(text, encoding="utf-8")
        return records
    return parse_cds_fasta(cache.read_text(encoding="utf-8"))
