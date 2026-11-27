"""Task 3: EBI Clustal Omega MSA and tree (fetch, validate, cache), tree comparison, conserved regions."""

from __future__ import annotations

import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import requests
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree

from phylogenetics.dataset import DEFAULT_EMAIL, CdsRecord
from phylogenetics.tree import splits

EBI_BASE: str = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"


def _run_job(payload: dict[str, str], poll_s: float, timeout_s: float) -> str:
    """Submit a Clustal Omega job and wait for FINISHED; return the job id."""
    resp = requests.post(f"{EBI_BASE}/run", data=payload, timeout=30.0)
    if resp.status_code != 200:
        raise RuntimeError(
            f"EBI submission failed ({resp.status_code}): {resp.text[:200]}"
        )
    job_id = resp.text.strip()
    deadline = time.monotonic() + timeout_s
    while True:
        state = requests.get(f"{EBI_BASE}/status/{job_id}", timeout=30.0).text.strip()
        if state == "FINISHED":
            return job_id
        if state in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"EBI job {job_id} ended with status {state}")
        if time.monotonic() >= deadline:
            raise RuntimeError(f"Timed out waiting for EBI job {job_id} (last={state})")
        time.sleep(poll_s)


def _result(job_id: str, result_type: str) -> str:
    resp = requests.get(f"{EBI_BASE}/result/{job_id}/{result_type}", timeout=30.0)
    if resp.status_code != 200:
        raise RuntimeError(
            f"EBI result {result_type} for {job_id} failed ({resp.status_code})"
        )
    return resp.text


def fetch_clustalo(
    records: list[CdsRecord],
    email: str,
    poll_s: float = 3.0,
    timeout_s: float = 300.0,
) -> tuple[str, str]:
    """Run Clustal Omega at EBI; return ``(clustal_text, newick_text)``. No other aligner as fallback."""
    payload = {
        "email": email,
        "sequence": "".join(f">{r.label}\n{r.seq}\n" for r in records),
        "stype": "dna",
        "outfmt": "clustal_num",
    }
    job_id = _run_job(payload, poll_s, timeout_s)
    return _result(job_id, "aln-clustal_num"), _result(job_id, "phylotree")


def parse_validated_msa(text: str, records: list[CdsRecord]) -> MultipleSeqAlignment:
    """Parse Clustal text; it must hold exactly ``records`` with unchanged sequences. Rows in record order."""
    by_id = {rec.id: rec for rec in AlignIO.read(StringIO(text), "clustal")}
    expected = [r.label for r in records]
    if set(by_id) != set(expected):
        raise ValueError(
            f"MSA ids {sorted(by_id)} do not match dataset {sorted(expected)}"
        )
    for r in records:
        if str(by_id[r.label].seq).replace("-", "").upper() != r.seq:
            raise ValueError(
                f"Ungapped MSA row for {r.label} differs from its dataset CDS"
            )
    return MultipleSeqAlignment([by_id[label] for label in expected])


def parse_validated_tree(text: str, records: list[CdsRecord]) -> Tree:
    """Parse Newick; its leaf names must equal the dataset labels exactly."""
    tree = Phylo.read(StringIO(text), "newick")
    leaves = sorted(t.name for t in tree.get_terminals())
    expected = sorted(r.label for r in records)
    if leaves != expected:
        raise ValueError(f"Tree leaves {leaves} do not match dataset {expected}")
    return tree


def load_or_fetch_clustalo(
    records: list[CdsRecord],
    msa_cache: Path,
    tree_cache: Path,
    refresh: bool = False,
    email: str = DEFAULT_EMAIL,
) -> tuple[MultipleSeqAlignment, Tree]:
    """Read cached MSA + tree, or fetch both; a fetched pair replaces the caches only if both validate."""
    msa_cache, tree_cache = Path(msa_cache), Path(tree_cache)
    if refresh or not (msa_cache.exists() and tree_cache.exists()):
        aln_text, tree_text = fetch_clustalo(records, email)
        aln = parse_validated_msa(aln_text, records)
        tree = parse_validated_tree(tree_text, records)
        msa_cache.parent.mkdir(parents=True, exist_ok=True)
        msa_cache.write_text(aln_text, encoding="utf-8")
        tree_cache.write_text(tree_text, encoding="utf-8")
        return aln, tree
    return (
        parse_validated_msa(msa_cache.read_text(encoding="utf-8"), records),
        parse_validated_tree(tree_cache.read_text(encoding="utf-8"), records),
    )


# --- Tree comparison --------------------------------------------------------------------------

WINDOW: int = 30
# 80% found no window on the 12-species MSA; 60% is ~3x the 21.5% background identity.
MIN_IDENTITY: float = 0.6
REFERENCE: str = "Homo_sapiens"


@dataclass(frozen=True)
class RFResult:
    """Robinson-Foulds comparison of two trees' unrooted splits."""

    shared: set[frozenset[str]]
    only_ours: set[frozenset[str]]
    only_clustal: set[frozenset[str]]
    rf: int
    normalized: float


def robinson_foulds(ours: Tree, theirs: Tree) -> RFResult:
    """RF = splits in one tree but not the other; normalized by the maximum ``2(n - 3)``."""
    leaves = {t.name for t in ours.get_terminals()}
    if leaves != {t.name for t in theirs.get_terminals()}:
        raise ValueError("Trees have different leaf sets")
    a, b = splits(ours), splits(theirs)
    rf = len(a - b) + len(b - a)
    return RFResult(a & b, a - b, b - a, rf, rf / (2 * (len(leaves) - 3)))


# --- Conserved regions ------------------------------------------------------------------------


@dataclass(frozen=True)
class ConservedRegion:
    """MSA columns [aln_start, aln_end) and the matching 1-based nt/codon span in the reference CDS."""

    aln_start: int
    aln_end: int
    identity: float
    ref_nt: tuple[int, int]
    ref_codons: tuple[int, int]


@dataclass(frozen=True)
class ConservationSummary:
    columns: int
    identical_columns: int
    longest_run: int
    regions: list[ConservedRegion]
    best_window: ConservedRegion


def _region(ref_row: str, flags: list[bool], start: int, end: int) -> ConservedRegion:
    before = sum(ch != "-" for ch in ref_row[:start])
    upto = sum(ch != "-" for ch in ref_row[:end])
    nt = (before + 1, upto)
    return ConservedRegion(
        start,
        end,
        sum(flags[start:end]) / (end - start),
        nt,
        ((nt[0] - 1) // 3 + 1, (nt[1] - 1) // 3 + 1),
    )


def conserved_regions(
    aln: MultipleSeqAlignment,
    reference: str = REFERENCE,
    window: int = WINDOW,
    min_identity: float = MIN_IDENTITY,
) -> ConservationSummary:
    """Windows of ``window`` columns with >= ``min_identity`` fully identical, ungapped columns,
    merged when they overlap or touch, then trimmed to their first..last identical column."""
    rows = [str(r.seq).upper() for r in aln]
    ref_row = next(str(r.seq).upper() for r in aln if r.id == reference)
    length = aln.get_alignment_length()
    flags = [("-" not in col) and len(set(col)) == 1 for col in zip(*rows)]

    longest = run = 0
    for f in flags:
        run = run + 1 if f else 0
        longest = max(longest, run)

    win = min(window, length)
    counts = [sum(flags[s : s + win]) for s in range(length - win + 1)]
    merged: list[list[int]] = []
    for s, c in enumerate(counts):
        if c < min_identity * win:
            continue
        if merged and s <= merged[-1][1]:
            merged[-1][1] = s + win
        else:
            merged.append([s, s + win])
    regions = []
    for start, end in merged:
        idx = [c for c in range(start, end) if flags[c]]
        regions.append(_region(ref_row, flags, idx[0], idx[-1] + 1))

    best = max(range(len(counts)), key=lambda s: (counts[s], -s))
    return ConservationSummary(
        length,
        sum(flags),
        longest,
        regions,
        _region(ref_row, flags, best, best + win),
    )


def _excerpt(aln: MultipleSeqAlignment, start: int, end: int) -> str:
    w = max(len(r.id) for r in aln)
    lines = [f"{r.id:<{w}}  {str(r.seq)[start:end]}" for r in aln]
    rows = [str(r.seq).upper() for r in aln]
    stars = "".join(
        "*" if ("-" not in col and len(set(col)) == 1) else " "
        for col in zip(*(row[start:end] for row in rows))
    )
    return "\n".join([*lines, " " * (w + 2) + stars])


def _fmt_region(r: ConservedRegion) -> str:
    return (
        f"MSA cols {r.aln_start + 1}..{r.aln_end} ({r.aln_end - r.aln_start} cols, "
        f"{100 * r.identity:.1f}% identical) -> {REFERENCE} nt {r.ref_nt[0]}..{r.ref_nt[1]}, "
        f"codons {r.ref_codons[0]}..{r.ref_codons[1]}"
    )


def format_task3_report(
    rf: RFResult,
    ours: dict[str, bool],
    theirs: dict[str, bool],
    cons: ConservationSummary,
    aln: MultipleSeqAlignment,
) -> str:
    """RF comparison, side-by-side monophyly, conservation stats and starred excerpts."""

    def fmt(splits_: set[frozenset[str]]) -> list[str]:
        return [
            "    {" + ", ".join(sorted(s)) + "}"
            for s in sorted(splits_, key=lambda s: (len(s), sorted(s)))
        ] or ["    (none)"]

    lines = [
        "=== TASK 3: NJ tree vs Clustal Omega (EBI) ===",
        (
            f"Robinson-Foulds distance: {rf.rf} (normalized {rf.normalized:.2f}); "
            f"shared splits: {len(rf.shared)}"
        ),
        "  Splits only in our NJ tree:",
        *fmt(rf.only_ours),
        "  Splits only in the Clustal Omega tree:",
        *fmt(rf.only_clustal),
        "  Group            NJ    Clustal",
        *(
            f"  {g:<15}  {'yes' if ok else 'no':<4}  {'yes' if theirs[g] else 'no'}"
            for g, ok in ours.items()
        ),
        "",
        (
            f"MSA: {cons.columns} columns, {cons.identical_columns} fully identical "
            f"({100 * cons.identical_columns / cons.columns:.1f}%), longest identical run "
            f"{cons.longest_run} nt"
        ),
        f"Conserved regions (>= {MIN_IDENTITY:.0%} identical columns in a {WINDOW}-column window):",
    ]
    if not cons.regions:
        lines += [
            "  none qualify; best window:",
            "  " + _fmt_region(cons.best_window),
            _excerpt(aln, cons.best_window.aln_start, cons.best_window.aln_end),
        ]
    for r in cons.regions:
        lines += ["  " + _fmt_region(r), _excerpt(aln, r.aln_start, r.aln_end)]
    return "\n".join(lines)
