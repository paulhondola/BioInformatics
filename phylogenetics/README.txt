Lab 04 — Phylogenetics (Paul Hondola). Python 3.12+; deps: biopython, matplotlib, requests; pytest for tests; pdflatex only for --pdf.
1. Full pipeline (offline, uses cached data/): uv sync && uv run phylogenetics
2. One task only: uv run phylogenetics --task {1|2|3|bonus}
3. Compile notes.pdf and build the ZIP: uv run phylogenetics --pdf --zip
4. Tests (no network): uv run pytest -v
Outputs in results/: distance_matrix_paulhondola.csv, pairwise_paulhondola.csv, tree_paulhondola.nwk (reused by the clustering lab), tree_paulhondola.png.
Data: data/tp53_cds.fasta = TP53 CDS of 12 vertebrates from NCBI RefSeq (accessions in notes.pdf; chimp XM_001172077.6 is a predicted model, no NM_ exists).
Refresh: --refresh-data re-downloads it; --refresh-msa reruns Clustal Omega at EBI (update notes.tex afterwards).
Clustal Omega: run on the EBI web form, job clustalo-I20261127-003011-0299-43071492-p1m -> data/clustalo_msa.aln + data/clustalo_phylotree.nwk.
