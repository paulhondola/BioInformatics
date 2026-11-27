from phylogenetics.dataset import GROUPS, PROJECT_ROOT, load_or_fetch_dataset
from phylogenetics.distances import compute_distances
from phylogenetics.msa import conserved_regions, load_or_fetch_clustalo, robinson_foulds
from phylogenetics.tree import build_nj_tree, monophyly

NOTES = (PROJECT_ROOT / "src/phylogenetics/notes.tex").read_text(encoding="utf-8")


def test_readme_is_short_and_names_the_cli():
    lines = (PROJECT_ROOT / "README.txt").read_text(encoding="utf-8").splitlines()
    assert len(lines) <= 10
    assert any("uv run phylogenetics" in line for line in lines)


def test_notes_numbers_match_the_pipeline():
    recs = load_or_fetch_dataset(PROJECT_ROOT / "data/tp53_cds.fasta")
    dist = compute_distances(recs)
    nj = build_nj_tree(dist.labels, dist.matrix)
    aln, clustal = load_or_fetch_clustalo(
        recs,
        PROJECT_ROOT / "data/clustalo_msa.aln",
        PROJECT_ROOT / "data/clustalo_phylotree.nwk",
    )
    rf = robinson_foulds(nj.tree, clustal)
    assert f"{rf.normalized:.2f}" in NOTES and f"RF={rf.rf}" in NOTES
    for r in conserved_regions(aln).regions:
        assert f"{r.ref_codons[0]}–{r.ref_codons[1]}" in NOTES
    assert nj.rooting == "outgroup"  # notes say the tree is rooted on the fish
    assert not monophyly(nj.tree, GROUPS)["Laurasiatheria"]  # notes discuss the dog


def test_notes_name_the_saved_tree():
    assert r"tree\_paulhondola.nwk" in NOTES


ABBR = {
    "Homo_sapiens": "Hsa",
    "Pan_troglodytes": "Ptr",
    "Macaca_mulatta": "Mml",
    "Mus_musculus": "Mmu",
    "Rattus_norvegicus": "Rno",
    "Canis_lupus_familiaris": "Cfa",
    "Bos_taurus": "Bta",
    "Sus_scrofa": "Ssc",
    "Gallus_gallus": "Gga",
    "Xenopus_tropicalis": "Xtr",
    "Danio_rerio": "Dre",
    "Oryzias_latipes": "Ola",
}


def test_notes_matrix_cells_match_distances():
    dist = compute_distances(
        load_or_fetch_dataset(PROJECT_ROOT / "data/tp53_cds.fasta")
    )
    for i, label in enumerate(dist.labels[1:], start=1):
        expected = [f"{v:.3f}"[1:] for v in dist.matrix[i][:i]]
        row = next(
            line
            for line in NOTES.splitlines()
            if line.startswith(rf"\textbf{{{ABBR[label]}}} &")
        )
        cells = [c.strip() for c in row.removesuffix(r"\\").split("&")[1:]]
        assert cells[:i] == expected, label


def test_notes_r273_codon_claim_matches_msa():
    recs = load_or_fetch_dataset(PROJECT_ROOT / "data/tp53_cds.fasta")
    aln, _ = load_or_fetch_clustalo(
        recs,
        PROJECT_ROOT / "data/clustalo_msa.aln",
        PROJECT_ROOT / "data/clustalo_phylotree.nwk",
    )
    human = next(str(r.seq) for r in aln if r.id == "Homo_sapiens")
    col = [i for i, ch in enumerate(human) if ch != "-"][(273 - 1) * 3]
    codon = {r.id: str(r.seq)[col : col + 3] for r in aln}
    assert codon["Homo_sapiens"] == "CGT" and codon["Oryzias_latipes"] == "AGA"
    if {codon[m] for m in GROUPS["Mammals"]} != {"CGT"}:
        assert r"\texttt{CGT} la mamifere" not in NOTES
        assert r"\texttt{CGT} la om" in NOTES
