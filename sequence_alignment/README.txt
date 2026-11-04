Lab 02 — Sequence Alignment (Paul Hondola)
Python: 3.12+
Dependencies: biopython>=1.84, polars>=1.0, requests>=2.31; pdflatex only for --pdf
1. Run full pipeline with uv: uv sync && uv run sequence-alignment
2. Run in active venv without uv: pip install biopython polars requests && PYTHONPATH=src python3 -m sequence_alignment
3. Run single task only: uv run sequence-alignment --task {1|2|3|bonus}
4. Compile PDF and build ZIP: uv run sequence-alignment --pdf --zip
5. Re-query EBI MSA (replaces data/clustalo_msa.aln; update notes.tex afterwards): uv run sequence-alignment --task 3 --refresh-msa
6. Run automated tests: uv run pytest -v
