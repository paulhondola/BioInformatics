Lab 03 — Formats & NGS
Author: Paul Hondola (GitHub: paulhondola)
Python: 3.12+. Dependencies: biopython, matplotlib; pytest for tests; pdflatex only for --pdf.

Running
1. Full pipeline: uv sync && uv run fastq-and-ngs
2. One task only: uv run fastq-and-ngs --task {1|2|3|4}
3. Re-query PubMed (replaces the cache in data/; update notes.tex afterwards): uv run fastq-and-ngs --refresh
4. Compile PDF and build ZIP: uv run fastq-and-ngs --pdf --zip
5. Tests (no network needed): uv run pytest -v
Single scripts: uv run python -m fastq_and_ngs.pubmed_query (or .vcf_pubmed, .ex04_fastq_stats)
Outputs are in results/.

FASTQ source
SRR31709504, read 1 of a pair (SRR31709504_1.fastq.gz), study PRJNA1198040: cDNA sequencing of TP53 mRNA
transcripts from a CRISPR-edited MOLM13 AML cell line, Illumina MiSeq, amplicon.
https://www.ebi.ac.uk/ena/browser/view/SRR31709504
The file is not in this ZIP. To rerun Task 3, place it at data/work/paulhondola/lab03/your_reads.fastq.gz.
Without it, the full pipeline warns, skips Task 3 and keeps the existing QC report and plot.

VCF source
data/tp53_variants.vcf is built from ClinVar records (GRCh38, accessed 2026-11-04), not the course sample.
ClinVar IDs: 12351, 12374, 12356, 5031076, 5031067. The last two have no rsID. REF and ALT are forward-strand.

AI disclosure
- Google Antigravity: earlier help with the ENA API, the FASTQ download and the QC plotting script.
- Claude (Claude Code): design, implementation of this project and its tests, and a draft of the notes.
- Verification: automated tests; the FASTQ MD5 was compared with ENA; the QC numbers were recounted separately;
  variant coordinates were read from ClinVar; PubMed's own query translations were inspected.
