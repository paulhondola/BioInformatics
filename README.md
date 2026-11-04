# BioInformatics — Lab Projects

Python lab assignments (managed with `uv`, one project per folder). Each has its brief in `docs/task.md` (Romanian).

| Folder | Topic | Tasks (points) |
|---|---|---|
| [`sequence_alignment/`](sequence_alignment/) | Sequence alignment | Pairwise distances (3), global vs. local alignment (4), online MSA with Clustal Omega (3), bonus semiglobal (+1) |
| [`fastq_and_ngs/`](fastq_and_ngs/) | Formats & NGS | PubMed parsing (2), FASTA copy (2), FASTQ QC (3), VCF → PubMed (3), bonus QC plot (+1) |
| [`phylogenetics/`](phylogenetics/) | Phylogenetics | Distance matrix for ≥10 sequences (3), Neighbor-Joining tree + cluster interpretation (4), compare with Clustal Omega MSA and mark conserved regions (3), bonus tree plot (+1) |
| [`clustering/`](clustering/) | Clustering + phylogeny (extends the phylogenetics tree) | Reuse dataset and tree (1), two clustering methods with Silhouette (3), dendrogram/PCA/annotated tree (3), clusters vs. clades (2), report (1), bonus third method (+1) |

## Running a project

```bash
cd <folder>
uv sync
uv run <package-name>      # e.g. uv run sequence-alignment
uv run pytest -v
```

Python 3.12+. See each project's README for details (`sequence_alignment/README.txt` is the only one filled in).
