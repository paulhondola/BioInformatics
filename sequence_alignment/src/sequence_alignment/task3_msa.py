"""Task 3: Online EBI MSA (Clustal Omega with EBI MUSCLE Clustal-format fallback) & motif analysis."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import time

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import requests

from sequence_alignment.sequences import SequenceRecord

EBI_BASE: str = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
EBI_MUSCLE_BASE: str = "https://www.ebi.ac.uk/Tools/services/rest/muscle"


@dataclass(frozen=True)
class ConservedMotif:
    """Conserved MSA slice with per-species 1-based transcript coordinates."""

    aln_start: int
    aln_end: int
    length: int
    identity_pct: float
    consensus_seq: str
    per_seq_coords: dict[str, tuple[int, int]]
    clustal_excerpt: str


@dataclass(frozen=True)
class MsaAnalysis:
    """Summary statistics and conserved motifs from a MultipleSeqAlignment."""

    num_sequences: int
    alignment_length: int
    fully_conserved_cols: int
    conserved_pct: float
    top_contiguous_motif: ConservedMotif
    best_window_motif: ConservedMotif


def _submit_and_poll_ebi_job(
    base_url: str,
    payload: dict[str, str],
    result_type: str,
    poll_interval_s: float,
    timeout_s: float,
    queue_timeout_s: float | None = None,
) -> tuple[str | None, str]:
    """Submit an EBI REST job; returns ``(result_text, 'FINISHED')`` or ``(None, 'QUEUED_STALL')``."""
    run_resp = requests.post(f"{base_url}/run", data=payload, timeout=30.0)
    if run_resp.status_code != 200:
        raise RuntimeError(
            f"EBI MSA submission failed ({run_resp.status_code}): {run_resp.text[:200]}"
        )
    job_id = run_resp.text.strip()

    start = time.monotonic()
    deadline = start + timeout_s
    while True:
        status_resp = requests.get(f"{base_url}/status/{job_id}", timeout=30.0)
        if status_resp.status_code != 200:
            raise RuntimeError(
                f"EBI status check failed ({status_resp.status_code}) for job {job_id}"
            )
        state = status_resp.text.strip()
        if state == "FINISHED":
            break
        if state in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"EBI MSA job {job_id} terminated with status {state}")
        now = time.monotonic()
        if state == "QUEUED" and queue_timeout_s is not None and (now - start) >= queue_timeout_s:
            return None, "QUEUED_STALL"
        if now >= deadline:
            raise RuntimeError(f"Timed out waiting for EBI MSA job {job_id} (last={state})")
        time.sleep(poll_interval_s)

    result_resp = requests.get(f"{base_url}/result/{job_id}/{result_type}", timeout=30.0)
    if result_resp.status_code != 200:
        raise RuntimeError(
            f"Failed to download EBI MSA result ({result_resp.status_code}) for job {job_id}"
        )
    return result_resp.text, "FINISHED"


def fetch_clustalo_msa(
    records: list[SequenceRecord],
    email: str = "paulhondola@gmail.com",
    poll_interval_s: float = 2.0,
    timeout_s: float = 180.0,
    clustalo_queue_timeout_s: float | None = 45.0,
) -> str:
    """Submit ``records`` to EBI Clustal Omega (falling back to EBI MUSCLE Clustal format if queue stalls)."""
    fasta_input = "\n".join(f">{r.accession}\n{r.seq}" for r in records) + "\n"
    queue_cap = (
        min(clustalo_queue_timeout_s, timeout_s)
        if clustalo_queue_timeout_s is not None and clustalo_queue_timeout_s < timeout_s
        else None
    )
    text, outcome = _submit_and_poll_ebi_job(
        EBI_BASE,
        {
            "email": email,
            "sequence": fasta_input,
            "stype": "dna",
            "outfmt": "clustal_num",
        },
        result_type="aln-clustal_num",
        poll_interval_s=poll_interval_s,
        timeout_s=timeout_s,
        queue_timeout_s=queue_cap,
    )
    if text is not None:
        return text

    # Automatic fallback to EMBL-EBI MUSCLE in Clustal format when clustalo SLURM queue is stalled
    muscle_text, _ = _submit_and_poll_ebi_job(
        EBI_MUSCLE_BASE,
        {
            "email": email,
            "sequence": fasta_input,
            "format": "clw",
        },
        result_type="aln-clustalw",
        poll_interval_s=poll_interval_s,
        timeout_s=timeout_s,
        queue_timeout_s=None,
    )
    assert muscle_text is not None
    return muscle_text


def load_or_fetch_msa(
    records: list[SequenceRecord],
    cache_path: Path,
    refresh: bool = False,
    email: str = "paulhondola@gmail.com",
) -> MultipleSeqAlignment:
    """Load cached Clustal-format MSA from ``cache_path`` or fetch from EBI if needed."""
    cache_file = Path(cache_path)
    if refresh or not cache_file.exists():
        raw_aln = fetch_clustalo_msa(records, email=email)
        msa = _parse_validated_msa(raw_aln, records)  # validate first: a bad reply must not replace the cache
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(raw_aln, encoding="utf-8")
        return msa
    return _parse_validated_msa(cache_file.read_text(encoding="utf-8"), records)


def _parse_validated_msa(text: str, records: list[SequenceRecord]) -> MultipleSeqAlignment:
    """Parse Clustal ``text``, check it holds exactly ``records``, and return them in record order."""
    by_id = {rec.id: rec for rec in AlignIO.read(StringIO(text), "clustal")}
    expected_ids = [rec.accession for rec in records]
    if set(by_id) != set(expected_ids):
        raise ValueError(
            f"Sequence IDs in MSA ({sorted(by_id)}) do not match records ({sorted(expected_ids)})"
        )
    for rec in records:
        if str(by_id[rec.accession].seq).replace("-", "").upper() != rec.seq:
            raise ValueError(f"Ungapped MSA sequence for {rec.accession} differs from its FASTA record")
    return MultipleSeqAlignment([by_id[acc] for acc in expected_ids])


def _build_motif(
    rows: list[tuple[str, str]],
    conserved_flags: list[bool],
    start: int,
    end: int,
) -> ConservedMotif:
    length = max(0, end - start)
    conserved_count = sum(conserved_flags[start:end])
    identity_pct = (100.0 * conserved_count / length) if length > 0 else 0.0

    consensus_chars: list[str] = []
    for col in range(start, end):
        if conserved_flags[col]:
            consensus_chars.append(rows[0][1][col])
        else:
            non_gaps = [seq[col] for _, seq in rows if seq[col] != "-"]
            consensus_chars.append(Counter(non_gaps).most_common(1)[0][0] if non_gaps else "-")

    per_seq_coords: dict[str, tuple[int, int]] = {}
    label_w = max(len(seq_id) for seq_id, _ in rows)
    excerpt_lines: list[str] = []
    for seq_id, seq in rows:
        before = sum(ch != "-" for ch in seq[:start])
        within = sum(ch != "-" for ch in seq[start:end])
        coord = (before + 1, before + within) if within > 0 else (before, before)
        per_seq_coords[seq_id] = coord
        excerpt_lines.append(f"{seq_id:<{label_w}}  {seq[start:end]}  {coord[0]}-{coord[1]}")

    stars = "".join("*" if conserved_flags[c] else " " for c in range(start, end))
    excerpt_lines.append(" " * (label_w + 2) + stars)

    return ConservedMotif(
        aln_start=start,
        aln_end=end,
        length=length,
        identity_pct=identity_pct,
        consensus_seq="".join(consensus_chars),
        per_seq_coords=per_seq_coords,
        clustal_excerpt="\n".join(excerpt_lines),
    )


def analyze_msa(alignment: MultipleSeqAlignment, window: int = 45) -> MsaAnalysis:
    """Compute overall MSA conservation and extract the strongest conserved motifs."""
    rows = [(rec.id, str(rec.seq).upper()) for rec in alignment]
    num_seqs = len(rows)
    aln_len = alignment.get_alignment_length()

    conserved_flags: list[bool] = []
    for col in range(aln_len):
        chars = [seq[col] for _, seq in rows]
        conserved_flags.append(("-" not in chars) and (len(set(chars)) == 1))

    fully_conserved = sum(conserved_flags)
    conserved_pct = (100.0 * fully_conserved / aln_len) if aln_len > 0 else 0.0

    best_run_start = 0
    best_run_end = 0
    cur_start: int | None = None
    for idx, flag in enumerate(conserved_flags):
        if flag and cur_start is None:
            cur_start = idx
        elif not flag and cur_start is not None:
            if idx - cur_start > best_run_end - best_run_start:
                best_run_start, best_run_end = cur_start, idx
            cur_start = None
    if cur_start is not None and aln_len - cur_start > best_run_end - best_run_start:
        best_run_start, best_run_end = cur_start, aln_len

    win = min(window, aln_len)
    best_win_start = 0
    best_win_key = (-1, -1)
    for s in range(0, max(1, aln_len - win + 1)):
        stars = sum(conserved_flags[s : s + win])
        non_gaps = sum(seq[c] != "-" for _, seq in rows for c in range(s, s + win))
        key = (stars, non_gaps)
        if key > best_win_key:
            best_win_key = key
            best_win_start = s

    top_run = _build_motif(rows, conserved_flags, best_run_start, best_run_end)
    best_win = _build_motif(rows, conserved_flags, best_win_start, best_win_start + win)

    return MsaAnalysis(
        num_sequences=num_seqs,
        alignment_length=aln_len,
        fully_conserved_cols=fully_conserved,
        conserved_pct=conserved_pct,
        top_contiguous_motif=top_run,
        best_window_motif=best_win,
    )


def format_task3_report(analysis: MsaAnalysis) -> str:
    """Render Task 3 EBI MSA statistics and conserved motif excerpts."""
    run_m = analysis.top_contiguous_motif
    win_m = analysis.best_window_motif
    return "\n".join(
        [
            "=== TASK 3: Online Multiple Sequence Alignment (EBI Clustal format) ===",
            (
                f"Sequences: {analysis.num_sequences} | Alignment columns: {analysis.alignment_length} | "
                f"Fully conserved (*): {analysis.fully_conserved_cols} ({analysis.conserved_pct:.2f}%)"
            ),
            "",
            (
                f"[Longest 100%-identical motif: {run_m.length} nt at MSA cols "
                f"{run_m.aln_start + 1}..{run_m.aln_end} -> {run_m.consensus_seq}]"
            ),
            run_m.clustal_excerpt,
            "",
            (
                f"[Highest-conservation {win_m.length}-col region: {win_m.identity_pct:.1f}% identical "
                f"at MSA cols {win_m.aln_start + 1}..{win_m.aln_end}]"
            ),
            win_m.clustal_excerpt,
        ]
    )
