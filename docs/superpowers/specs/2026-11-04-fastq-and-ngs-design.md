# Lab 03 (Formats & NGS): design

Date: 2026-11-04. Brief: `fastq_and_ngs/docs/task.md` (Romanian). Reference for style: `sequence_alignment/`.

## Intent

Paul wants this lab completed for full marks (10 points + 1 bonus), in the same style as the sequence alignment lab, using his own data wherever the brief caps marks for sample data. One command should regenerate every deliverable, and the submission must be honest about where its data came from and what its searches really found.

## Goals

- Produce every deliverable named in the brief, with the handle `paulhondola`:
  - scripts: `pubmed_query.py`, `vcf_pubmed.py`, `ex04_fastq_stats.py`
  - outputs: `pubmed_paulhondola.txt`, `fasta_paulhondola.fasta`, `qc_report_paulhondola.txt`, `variants_paulhondola.txt`, `qc_plot_paulhondola.png`
  - `README.txt`, `notes.pdf`, and a ZIP
- Use the real FASTQ and a VCF that is not the course sample, so neither mark cap applies.
- Tests run offline.

## Non-goals

- Calling variants from the reads. Paul chose public variants for the VCF.
- Read trimming or alignment.
- Clinical interpretation of the variants.
- Changing the course repo or its PRs. Paul decided this assignment is submitted only as a Moodle ZIP; the course repo's PRs are for the separate weekly exercises. No export into `submissions/` is needed.
- Uploading the ZIP. Paul does that himself; this work only builds it.

## Anti-goals

Any of these is a failure even if the code runs:

- An invented or unverifiable variant in the VCF.
- Presenting no-ID search hits as literature about that specific position (see Task 4).
- Numbers in `notes.pdf` that disagree with the generated reports.

## Verified facts (checked 2026-11-04)

**FASTQ.** `/Users/paulhondola/Faculta/bioinf-y4-2026-2027-paulhondola/data/work/paulhondola/lab03/your_reads.fastq.gz`

- It is read 1 of paired run SRR31709504: its MD5 `1ba7d613baf286e55523f25bc4029bb1` equals ENA's `SRR31709504_1.fastq.gz`.
- Study PRJNA1198040, "cDNA sequencing of TP53-mRNA transcripts". Sample: MOLM13 isogenic AML cell line with one TP53 WT allele. Illumina MiSeq, amplicon.
- 81,220 reads, every read 142 nt, N rate 0, mean Phred 37.30.

**Paul's existing Lab 3 work.** Branch `feat/lab03-paulhondola` of the course repo (open PR #3) already holds a finished `ex04_fastq_stats.py`, `ex01_fetch_fastq.py`, `qc_report_paulhondola.txt`, `qc_plot_paulhondola.png` and `paulhondola_notes.md`. It holds no `pubmed_query.py`, `vcf_pubmed.py`, or their outputs. The committed report is:

```
Reads: 81220
Mean length: 142.00
N rate: 0.0000
Mean Phred: 37.30
```

This matches an independent recount of the FASTQ. The script's paths assume the course repo layout, and it runs everything at module level, so it can't be imported or tested as it stands. Its notes credit Google Antigravity as the AI tool used.

**Variants.** All from ClinVar, GRCh38. POS is the SPDI position + 1. REF and ALT come from SPDI, so they are on the forward strand (TP53 is on the reverse strand, so they differ from the c. notation).

| CHROM | POS | ID | REF | ALT | ClinVar | Change |
|---|---|---|---|---|---|---|
| 17 | 7676154 | rs1042522 | G | C | 12351 | c.215C>G, p.Pro72Arg |
| 17 | 7675088 | rs28934578 | C | T | 12374 | c.524G>A, p.Arg175His |
| 17 | 7674220 | rs11540652 | C | T | 12356 | c.743G>A, p.Arg248Gln |
| 17 | 7675091 | . | C | G | 5031076 | c.521G>C, p.Arg174Thr |
| 17 | 7674224 | . | T | G | 5031067 | c.739A>C, p.Asn247His |

PubMed hit counts on the day: rs1042522 351, rs28934578 1, rs11540652 3.

The course sample `tp53_demo.vcf` lists rs1042522 as REF C, ALT G, which is the transcript strand. ClinVar's forward-strand form is G>C. Do not copy from the sample.

**PubMed ignores the position in the brief's no-ID query.** `chr17:7675091 AND TP53` is translated by PubMed to `"chr17"[All Fields] AND (tp53 terms)`. It returns the same 13 papers as `chr17 AND TP53`, and the same as any other position. Quoting the term does not help.

## What gets built

The project is `fastq_and_ngs/`, in the same shape as `sequence_alignment/`: `uv` project, `src/fastq_and_ngs/` package, `__main__.py` CLI, frozen dataclasses, `format_*_report` functions, pytest, `notes.tex` compiled by `pdflatex`.

### Task 1: `pubmed_query.py` (2p)

Search PubMed for exactly `TP53 AND cancer`, at most 5 articles, using Biopython `Entrez`. Write title, authors and abstract per article to `pubmed_paulhondola.txt`. Articles with no abstract are marked as such and do not crash the run.

### Task 2: FASTA (2p)

Copy `sequence_alignment/data/human_tp53.fa` into `fastq_and_ngs/data/`, then write it as `fasta_paulhondola.fasta`, byte-identical. It holds one record, NM_000546.6.

### Task 3: `ex04_fastq_stats.py` (3p)

Start from Paul's finished script on `feat/lab03-paulhondola`, keeping its calculations and report format, and do not rewrite them from scratch. Only reshape it: functions instead of module-level code, and project-relative paths. Default input is `fastq_and_ngs/data/work/paulhondola/lab03/your_reads.fastq.gz`, the path the brief names, relative to the project. Output is `qc_report_paulhondola.txt` and must stay identical to the committed report. Division by zero on an empty file gives 0, as the existing script does.

The implementation copies the FASTQ from the course repo to that path once, as a setup step, and checks its MD5 against the value above. The file stays out of git and the ZIP, so it can be absent on a fresh clone or an unzipped submission:

- `--task 3` fails with an error naming the accession and the ENA source.
- `--task all` warns, skips Task 3 and still runs the others.
- In both cases, the existing `qc_report` and `qc_plot` are left untouched. A missing file must never be treated like an empty one and overwrite the real report.

The script imports nothing from the package, so it also runs as a single file, and it takes optional input and output paths so it can run from another folder.

### Task 4: `vcf_pubmed.py` (3p)

Reads `fastq_and_ngs/data/tp53_variants.vcf`: a tab-separated VCFv4.2 file holding the five records above, with header lines stating GRCh38 and ClinVar as the source and each record's ClinVar ID in INFO.

For each variant, search PubMed as the brief says: the rsID if present, otherwise `chr<CHROM>:<POS> AND TP53`. A CHROM already starting with `chr` must not become `chrchr`.

`variants_paulhondola.txt` records, per variant: the variant, the query sent, PubMed's `QueryTranslation`, the hit count, and up to 3 PMIDs with titles. Zero hits is a normal result.

No-ID variants carry a visible warning when the translation dropped the position, saying that the hits are about the chromosome and TP53 rather than this variant. The detection compares the translation with the position, so it stays correct if PubMed changes behaviour.

### Bonus: `qc_plot_paulhondola.png` (+1p)

Written by `ex04_fastq_stats.py`, keeping Paul's existing two panels and Romanian labels: read-length distribution and per-read mean Phred distribution. Use matplotlib's non-interactive backend. The length panel is a single bar at 142, which is the real result of this amplicon run.

### CLI: `uv run fastq-and-ngs`

- `--task {all,1,2,3,4}`, default `all`.
- `--refresh` re-queries PubMed. Otherwise raw PubMed responses cached in `data/` are used, the same idea as `clustalo_msa.aln`. A bad reply must not overwrite the cache.
- `--email`, default `paulhondola@gmail.com`, as the alignment lab does for EBI. NCBI requires an email for Entrez.
- `--pdf` and `--zip`, as in the alignment lab. `--refresh` cannot be combined with them, since the notes hard-code results.

Generated outputs go in `fastq_and_ngs/results/`. `README.txt` and `notes.pdf` sit at the project root.

### ZIP: `fastq_and_ngs_paulhondola.zip`

Same layout and exclusions as the alignment lab's ZIP: `README.txt`, `notes.pdf`, `pyproject.toml`, `.python-version`, `uv.lock`, `data/`, `results/`, `src/`, `tests/`. It leaves out the FASTQ, which is not a listed deliverable.

### `README.txt`

Author and handle (Paul Hondola, `paulhondola`), run steps, AI disclosure, FASTQ source with the SRR accession and link, a statement that it is read 1 of a pair, and the VCF's ClinVar source. The AI disclosure must cover both Google Antigravity (the earlier QC script) and Claude (this work). Paul must review it before submitting, because it is his statement.

### `notes.pdf` (Romanian, max 1 page)

- 3–5 sentences: why QC matters before variant calling, using the real numbers. Condense the reflection already in Paul's `paulhondola_notes.md`.
- 2–3 sentences: how the PubMed searches were formulated, including the finding that PubMed drops the position from the no-ID query.

Claude drafts both; Paul edits.

### Tests

Offline and self-contained:

- QC stats on a tiny generated FASTQ, plain and gzipped, plus an empty file.
- VCF parsing, including an ID of `.`, and query building for both paths.
- The position-dropped warning, for a translation with and without the position.
- Report formatting with stubbed PubMed data.
- FASTA copy is byte-identical.
- The committed VCF matches the table above.
- On the real FASTQ when present (skipped otherwise): 81,220 reads, 142.00, 0.0000, 37.30.

### Project setup changes

- `.python-version` to 3.12 and `requires-python = ">=3.12"`, matching the alignment lab.
- Dependencies: `biopython`, `matplotlib`; dev: `pytest`.
- `readme = "README.txt"`, removing the empty `README.md`.
- `.gitignore` copied from the alignment lab, plus `data/work/` so the 4.5 MB FASTQ stays out of git.
- Update the root `README.md` line saying `sequence_alignment/README.txt` is the only README filled in.

## Open items for Paul

- **Whether the ClinVar-derived VCF counts as "not the sample".** It is not the course file, but only the instructor can confirm that it earns the full 3 points.

## Calls Claude made

- Copying the five variants above: three well-studied hotspots or polymorphisms with PubMed hits, and two real no-ID records.
- Keeping the brief's no-ID query as written and adding the warning, rather than replacing the query.
- Up to 3 PMIDs per variant.
- `results/` folder, ZIP keeping the project layout rather than flat, FASTQ left out of the ZIP and git.
- Skipping Task 3 with a warning under `--task all` when the FASTQ is absent, rather than downloading it or crashing.
- The named scripts live in the package and run with `python -m fastq_and_ngs.<name>`. `vcf_pubmed.py` reuses the Entrez and cache code in `pubmed_query.py`; only `ex04_fastq_stats.py` is fully standalone.
- Reports written in English, as in the alignment lab's console output; notes in Romanian.

## Left to the builder

Dataclass shapes, cache file format, plot styling, function names other than the required script names.
