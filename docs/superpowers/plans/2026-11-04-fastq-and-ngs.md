# Lab 03 (Formats & NGS) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `fastq_and_ngs/` so one command regenerates every Lab 03 deliverable and packages the Moodle ZIP.

**Architecture:** A `uv` package shaped like `sequence_alignment/`. Three modules carry the names the brief requires, and `__main__.py` holds the CLI, FASTA copy, PDF build and ZIP build. Raw PubMed replies are cached in `data/`, so reruns and tests are offline.

**Tech Stack:** Python 3.12+, Biopython (`SeqIO`, `Entrez`), matplotlib, stdlib `xml.etree.ElementTree`, pytest, `pdflatex`.

**Spec:** `docs/superpowers/specs/2026-11-04-fastq-and-ngs-design.md`

## Global Constraints

- Handle `paulhondola` in every output name. Default email `paulhondola@gmail.com`.
- PubMed query for Task 1 is exactly `TP53 AND cancer`, at most 5 articles.
- No-ID variant query is exactly `chr<CHROM>:<POS> AND TP53`.
- Report text is English; `notes.pdf` is Romanian and at most 1 page; plot labels stay Romanian.
- `requires-python = ">=3.12"`; `.python-version` is `3.12`; dependencies only `biopython`, `matplotlib`; dev only `pytest`.
- The FASTQ never enters git or the ZIP (`data/work/` is excluded from both).
- Tests never touch the network.
- Do not change the course repo. Read from it only.
- All work happens on branch `feat/fastq-and-ngs`. Stage named paths only; `clustering/` and `phylogenetics/` are not part of this work.
- Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. Truncated or malformed FASTQ: the run errors and the existing `qc_report` and `qc_plot` stay untouched. Test in Task 1.
2. Empty or malformed PubMed reply during `--refresh`: an error, and the cache stays as it was. Test in Task 2.
3. VCF edited after caching, without `--refresh`: new queries are fetched rather than crashing or reusing stale results. Test in Task 3.
4. PubMed articles with no abstract, no authors, a collective author, a structured abstract, or inline markup in the title: all render without crashing or losing text. Test in Task 2.
5. `--zip` before outputs exist: a clear error naming the missing file, and no partial ZIP left behind. Test in Task 4.

---

### Task 1: Project setup and FASTQ QC

**Files:**
- Modify: `fastq_and_ngs/pyproject.toml`, `fastq_and_ngs/.python-version`, `fastq_and_ngs/src/fastq_and_ngs/__init__.py`
- Delete: `fastq_and_ngs/README.md`
- Create: `fastq_and_ngs/README.txt` (one placeholder line; Task 5 fills it), `fastq_and_ngs/.gitignore`
- Create: `fastq_and_ngs/src/fastq_and_ngs/ex04_fastq_stats.py`
- Create: `fastq_and_ngs/data/work/paulhondola/lab03/your_reads.fastq.gz` (ignored by git)
- Test: `fastq_and_ngs/tests/test_ex04.py`

**Interfaces:**
- Produces, in `fastq_and_ngs.ex04_fastq_stats`:
  - `FastqStats`, a frozen dataclass: `reads: int, mean_length: float, n_rate: float, mean_phred: float, read_lengths: list[int], read_mean_phreds: list[float]`
  - `compute_fastq_stats(path: Path) -> FastqStats`
  - `format_qc_report(stats: FastqStats) -> str`
  - `write_qc_plot(stats: FastqStats, out_png: Path) -> Path`
  - `write_qc_outputs(fastq: Path, report: Path, plot: Path) -> FastqStats`
  - `main(argv: list[str] | None = None) -> int`, with optional args `--input`, `--report`, `--plot`

- [ ] **Step 1: Create the branch and set up the project**

Run `git switch -c feat/fastq-and-ngs`. Copy `sequence_alignment/.gitignore` and append `data/work/`. In `pyproject.toml`: description, `readme = "README.txt"`, `requires-python = ">=3.12"`, the dependencies above, script `fastq-and-ngs = "fastq_and_ngs.__main__:main"`, `[tool.uv.build-backend] module-name = "fastq_and_ngs"`, `uv_build>=0.13.6,<0.14.0`. Set `__init__.py` to a docstring and `__all__: list[str] = []`, as in the alignment lab. Then:

```bash
mkdir -p fastq_and_ngs/data/work/paulhondola/lab03
cp /Users/paulhondola/Faculta/bioinf-y4-2026-2027-paulhondola/data/work/paulhondola/lab03/your_reads.fastq.gz fastq_and_ngs/data/work/paulhondola/lab03/
md5 -q fastq_and_ngs/data/work/paulhondola/lab03/your_reads.fastq.gz
cd fastq_and_ngs && uv sync
```

Expected: MD5 `1ba7d613baf286e55523f25bc4029bb1`, and `uv sync` succeeds.

- [ ] **Step 2: Write the failing tests**

`test_ex04.py` writes tiny FASTQ files into `tmp_path` (use Phred `I` = 40 and `+` = 10):

```python
def test_stats_on_plain_and_gzipped(tmp_path): ...
    # reads "ACGN"/"IIII" and "AC"/"++": both file kinds give
    assert (s.reads, s.mean_length, s.n_rate) == (2, 3.0, 1 / 6)
    assert s.mean_phred == pytest.approx(180 / 6)
    assert s.read_lengths == [4, 2] and s.read_mean_phreds == [40.0, 10.0]

def test_lowercase_n_counts(tmp_path): ...       # "acgn" gives n_rate 0.25

def test_report_has_exact_template_lines():
    assert format_qc_report(stats) == "Reads: 81220\nMean length: 142.00\nN rate: 0.0000\nMean Phred: 37.30\n"

def test_empty_file_gives_zeros(tmp_path): ...   # all four values 0, and write_qc_plot still writes a PNG

def test_missing_file_raises_and_keeps_report(tmp_path): ...  # FileNotFoundError; existing report text unchanged

def test_corrupt_fastq_keeps_report_and_plot(tmp_path): ...
    # truncated gzip and a seq/quality length mismatch each raise,
    # and write_qc_outputs leaves pre-existing report and plot bytes identical

@pytest.mark.skipif(not REAL_FASTQ.exists(), reason="FASTQ not present")
def test_real_fastq_matches_committed_report():
    assert format_qc_report(compute_fastq_stats(REAL_FASTQ)) == "Reads: 81220\nMean length: 142.00\nN rate: 0.0000\nMean Phred: 37.30\n"
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `cd fastq_and_ngs && uv run pytest tests/test_ex04.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 4: Port Paul's script**

Get the source with `git -C /Users/paulhondola/Faculta/bioinf-y4-2026-2027-paulhondola cat-file -p "feat/lab03-paulhondola:labs/03_formats&NGS/submissions/paulhondola/ex04_fastq_stats.py"`. Keep its calculations, four-line report and two plot panels with their titles and labels. Replace module-level code and globals with the functions above, and set `matplotlib.use("Agg")` before importing `pyplot`. The module imports nothing from the package. `write_qc_outputs` computes stats fully before opening any output. Defaults derive from `Path(__file__).resolve().parents[2]`: input `data/work/paulhondola/lab03/your_reads.fastq.gz`, outputs `results/qc_report_paulhondola.txt` and `results/qc_plot_paulhondola.png`. Create parent folders.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest tests/test_ex04.py -v`
Expected: all PASS, including the real-FASTQ test.

- [ ] **Step 6: Commit**

```bash
git add docs/superpowers/specs/2026-11-04-fastq-and-ngs-design.md docs/superpowers/plans/2026-11-04-fastq-and-ngs.md fastq_and_ngs/pyproject.toml fastq_and_ngs/uv.lock fastq_and_ngs/.python-version fastq_and_ngs/.gitignore fastq_and_ngs/README.txt fastq_and_ngs/docs fastq_and_ngs/src fastq_and_ngs/tests
git commit -m "Add Lab 03 project setup and FASTQ QC"
```

---

### Task 2: PubMed query

**Files:**
- Create: `fastq_and_ngs/src/fastq_and_ngs/pubmed_query.py`
- Create: `fastq_and_ngs/data/pubmed_tp53_cancer.xml` (cache, generated in Step 5)
- Test: `fastq_and_ngs/tests/test_pubmed_query.py`

**Interfaces:**
- Produces, in `fastq_and_ngs.pubmed_query`:
  - Constants `QUERY = "TP53 AND cancer"`, `MAX_ARTICLES = 5`, `DEFAULT_EMAIL = "paulhondola@gmail.com"`
  - `Article`, frozen: `pmid: str, title: str, authors: tuple[str, ...], abstract: str | None`
  - `SearchResult`, frozen: `count: int, pmids: list[str], query_translation: str`
  - `search_pubmed(term: str, retmax: int, email: str) -> SearchResult`
  - `fetch_articles_xml(pmids: list[str], email: str) -> str`
  - `parse_articles(xml: str) -> list[Article]`
  - `load_or_fetch_articles(cache_path: Path, refresh: bool = False, email: str = DEFAULT_EMAIL) -> list[Article]`
  - `format_pubmed_report(articles: list[Article]) -> str`
  - `main(argv: list[str] | None = None) -> int`, with `--refresh`, `--email`

- [ ] **Step 1: Write the failing tests**

Tests hold a small efetch-style XML string and monkeypatch `search_pubmed` and `fetch_articles_xml`:

```python
def test_parse_handles_real_world_shapes():
    # article 1: LastName+Initials authors, CollectiveName, title with <i>,
    #            AbstractText Label="BACKGROUND" and Label="RESULTS"
    # article 2: no AuthorList, no Abstract
    assert a1.title == "TP53 in cancer"          # <i>TP53</i> keeps its text
    assert a1.authors == ("Smith J", "TP53 Consortium")
    assert a1.abstract == "BACKGROUND: first. RESULTS: second."
    assert (a2.authors, a2.abstract) == ((), None)

def test_report_lists_title_authors_abstract():
    text = format_pubmed_report([a1, a2])
    for needed in ("PMID: 1", "Title: TP53 in cancer", "Authors: Smith J, TP53 Consortium", "Abstract: BACKGROUND"):
        assert needed in text
    assert "Authors: (none listed)" in text and "Abstract: (none available)" in text

def test_uses_exact_query_and_cap(tmp_path, monkeypatch): ...
    # first run asserts search_pubmed got ("TP53 AND cancer", 5, email) and writes the cache
    # second run with both network functions raising still returns the articles

def test_refresh_refetches(tmp_path, monkeypatch): ...

@pytest.mark.parametrize("reply", ["", "<not xml", "<PubmedArticleSet></PubmedArticleSet>"])
def test_bad_reply_keeps_cache(tmp_path, monkeypatch, reply):
    # ValueError on refresh; cache file bytes unchanged
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_pubmed_query.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement the module**

Parse with `ElementTree` rather than `Entrez.read`, which can reach for DTDs online; malformed XML becomes `ValueError`. Read text with `itertext()`. Iterate `PubmedArticle` elements only. `search_pubmed` uses `Entrez.read(Entrez.esearch(...))` and reads `Count`, `IdList`, `QueryTranslation`. `load_or_fetch_articles` raises `ValueError` on zero parsed articles, and writes the cache only after a successful parse. `main` writes `results/pubmed_paulhondola.txt` using the default cache `data/pubmed_tp53_cancer.xml`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_pubmed_query.py -v`
Expected: all PASS

- [ ] **Step 5: Generate the real cache and check it**

Run: `uv run python -m fastq_and_ngs.pubmed_query --refresh`
Expected: `data/pubmed_tp53_cancer.xml` and `results/pubmed_paulhondola.txt` exist, with 5 articles. Open the report and confirm it reads correctly.

- [ ] **Step 6: Commit**

```bash
git add fastq_and_ngs/src/fastq_and_ngs/pubmed_query.py fastq_and_ngs/tests/test_pubmed_query.py fastq_and_ngs/data/pubmed_tp53_cancer.xml fastq_and_ngs/results/pubmed_paulhondola.txt
git commit -m "Add PubMed query for TP53 and cancer"
```

---

### Task 3: VCF and variant searches

**Files:**
- Create: `fastq_and_ngs/data/tp53_variants.vcf`, `fastq_and_ngs/data/variant_searches.json` (cache, generated in Step 5)
- Create: `fastq_and_ngs/src/fastq_and_ngs/vcf_pubmed.py`
- Test: `fastq_and_ngs/tests/test_vcf_pubmed.py`

**Interfaces:**
- Consumes: `search_pubmed`, `fetch_articles_xml`, `parse_articles`, `DEFAULT_EMAIL` from `pubmed_query`
- Produces, in `fastq_and_ngs.vcf_pubmed`:
  - `Variant`, frozen: `chrom: str, pos: int, vid: str | None, ref: str, alt: str`
  - `VariantHit`, frozen: `variant: Variant, query: str, query_translation: str, count: int, articles: tuple[tuple[str, str], ...]`
  - `parse_vcf(path: Path) -> list[Variant]`
  - `build_query(variant: Variant) -> str`
  - `position_dropped(variant: Variant, translation: str) -> bool`
  - `load_or_fetch_hits(variants: list[Variant], cache_path: Path, refresh: bool = False, email: str = DEFAULT_EMAIL) -> list[VariantHit]`
  - `format_variants_report(hits: list[VariantHit]) -> str`
  - `main(argv: list[str] | None = None) -> int`, with `--vcf`, `--refresh`, `--email`

- [ ] **Step 1: Write the failing tests**

```python
EXPECTED = [
    ("17", 7676154, "rs1042522", "G", "C"),
    ("17", 7675088, "rs28934578", "C", "T"),
    ("17", 7674220, "rs11540652", "C", "T"),
    ("17", 7675091, None, "C", "G"),
    ("17", 7674224, None, "T", "G"),
]

def test_committed_vcf_matches_spec():
    assert [(v.chrom, v.pos, v.vid, v.ref, v.alt) for v in parse_vcf(VCF)] == EXPECTED
    header = VCF.read_text()
    assert "##reference=GRCh38" in header and "ClinVar" in header
    for line in (l for l in header.splitlines() if not l.startswith("#")):
        assert line.count("\t") >= 7

def test_parse_accepts_space_aligned_vcf(tmp_path): ...     # like the course sample
def test_parse_rejects_fewer_than_two_variants(tmp_path): ...   # ValueError mentioning "at least 2"
def test_parse_rejects_short_line(tmp_path): ...            # ValueError naming the line number

def test_build_query():
    assert build_query(Variant("17", 1, "rs5", "A", "T")) == "rs5"
    assert build_query(Variant("17", 7675091, None, "C", "G")) == "chr17:7675091 AND TP53"
    assert build_query(Variant("chr17", 7675091, None, "C", "G")) == "chr17:7675091 AND TP53"
    assert build_query(Variant("17", 1, "rs5;rs6", "A", "T")) == "rs5 OR rs6"

def test_position_dropped():
    no_id = Variant("17", 7675091, None, "C", "G")
    assert position_dropped(no_id, '"chr17"[All Fields] AND "tp53"[All Fields]') is True
    assert position_dropped(no_id, '"chr17:7675091"[All Fields] AND "tp53"[All Fields]') is False
    assert position_dropped(Variant("17", 1, "rs5", "A", "T"), '"rs5"[All Fields]') is False

def test_report_shows_translation_zero_hits_and_warning(): ...
    # every hit shows variant, query, translation, count
    # count 0 prints "Hits: 0" with no PMID lines
    # a dropped-position hit contains "hits are about chr17 and TP53, not this variant"
    # a kept-position hit and an rsID hit do not contain it

def test_cache_avoids_network_and_caps_titles(tmp_path, monkeypatch): ...
    # fetches at most 3 articles per variant; a second run with network stubs raising returns the same hits

def test_new_vcf_record_fetches_only_new_query(tmp_path, monkeypatch): ...
def test_refresh_failure_keeps_cache(tmp_path, monkeypatch): ...   # a raise midway leaves cache bytes unchanged
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_vcf_pubmed.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `data/tp53_variants.vcf`**

Tab-separated VCFv4.2 with `##reference=GRCh38`, `##source=ClinVar (accessed 2026-11-04)`, and `##INFO` lines defining `CLNID` and `CLNHGVS`. Use the spec's five records, with QUAL and FILTER `.` and INFO such as `CLNID=12351;CLNHGVS=NM_000546.6:c.215C>G`.

- [ ] **Step 4: Implement the module**

`parse_vcf` splits on any whitespace and maps ID `.` to `None`. The cache is JSON keyed by query: `{query: {"count": int, "query_translation": str, "articles": [[pmid, title], ...]}}`. With zero PMIDs, skip `fetch_articles_xml`. Write the cache only after all variants succeed. `main` writes `results/variants_paulhondola.txt`.

Warning text, on its own line: `WARNING: PubMed dropped the position from this query; hits are about chr<CHROM> and TP53, not this variant.`

- [ ] **Step 5: Run the tests, then generate the real cache**

Run: `uv run pytest tests/test_vcf_pubmed.py -v`, then `uv run python -m fastq_and_ngs.vcf_pubmed --refresh`
Expected: tests PASS. The report has 5 variants, three rsID hit counts, and the warning on both no-ID variants. Counts may differ from the spec's 2026-11-04 values if PubMed has changed.

- [ ] **Step 6: Commit**

```bash
git add fastq_and_ngs/src/fastq_and_ngs/vcf_pubmed.py fastq_and_ngs/tests/test_vcf_pubmed.py fastq_and_ngs/data/tp53_variants.vcf fastq_and_ngs/data/variant_searches.json fastq_and_ngs/results/variants_paulhondola.txt
git commit -m "Add VCF to PubMed variant searches"
```

---

### Task 4: CLI, FASTA copy, PDF and ZIP

**Files:**
- Create: `fastq_and_ngs/src/fastq_and_ngs/__main__.py`
- Create: `fastq_and_ngs/data/human_tp53.fa` (copied from `sequence_alignment/data/human_tp53.fa`)
- Test: `fastq_and_ngs/tests/test_main.py`

**Interfaces:**
- Consumes: `write_qc_outputs` from `ex04_fastq_stats`; `load_or_fetch_articles`, `format_pubmed_report` from `pubmed_query`; `parse_vcf`, `load_or_fetch_hits`, `format_variants_report` from `vcf_pubmed`
- Produces, in `fastq_and_ngs.__main__`:
  - `copy_fasta(source: Path, dest: Path) -> Path`
  - `run_tasks(root: Path, task: str, refresh: bool, email: str) -> int`
  - `compile_notes_pdf(tex_path: Path, output_pdf: Path) -> Path`
  - `build_submission_zip(root: Path, output_zip: Path) -> Path`
  - `main(argv: list[str] | None = None) -> int`

- [ ] **Step 1: Write the failing tests**

Tests build a fake project in `tmp_path` and monkeypatch the network functions:

```python
def test_fasta_copy_is_byte_identical(tmp_path): ...

def test_task_all_without_fastq_warns_and_runs_others(tmp_path, monkeypatch, capsys):
    # with a pre-existing results/qc_report_paulhondola.txt
    assert run_tasks(root, "all", False, EMAIL) == 0
    assert "SRR31709504" in capsys.readouterr().err
    assert report.read_text() == before
    # pubmed, fasta and variants outputs now exist

def test_task_3_without_fastq_fails_clearly(tmp_path, capsys):
    assert run_tasks(root, "3", False, EMAIL) == 1
    assert "SRR31709504" in capsys.readouterr().err   # and no traceback; report unchanged

def test_refresh_cannot_combine_with_pdf_or_zip():
    with pytest.raises(SystemExit):
        main(["--refresh", "--zip"])

def test_zip_layout_and_exclusions(tmp_path):
    names = zipfile.ZipFile(build_submission_zip(root, out)).namelist()
    assert "results/variants_paulhondola.txt" in names and "notes.pdf" in names
    assert not any(n.startswith("data/work/") or "__pycache__" in n or n.endswith(".DS_Store") for n in names)

def test_zip_missing_output_names_file_and_leaves_no_zip(tmp_path):
    with pytest.raises(FileNotFoundError, match="qc_plot_paulhondola.png"):
        build_submission_zip(root, out)
    assert not out.exists()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_main.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement the module**

Reuse `_find_pdflatex` and `compile_notes_pdf` from `sequence_alignment/src/sequence_alignment/__main__.py` unchanged. Base `build_submission_zip` on the same file. Required before opening the ZIP: `README.txt`, `notes.pdf`, `pyproject.toml`, `.python-version`, and all five `results/` deliverables. Also include `uv.lock`, `data/` (except `data/work/`), `src/`, `tests/`, and `results/`. Args are `--task {all,1,2,3,4}`, `--refresh`, `--email`, `--pdf`, `--zip`; `--refresh` with `--pdf` or `--zip` is a `parser.error`. The ZIP is `fastq_and_ngs_paulhondola.zip` and the PDF source is `src/fastq_and_ngs/notes.tex`. Missing-FASTQ messages name `SRR31709504` and `https://www.ebi.ac.uk/ena/browser/view/SRR31709504`, and go to stderr. Corrupt-FASTQ errors in `run_tasks` return 1 the same way.

- [ ] **Step 4: Run the tests and the real pipeline**

Run: `uv run pytest -v`, then `uv run fastq-and-ngs`, then `cmp data/human_tp53.fa results/fasta_paulhondola.fasta`
Expected: all PASS; all five outputs in `results/`; `cmp` prints nothing. Compare the plot visually with Paul's: two panels, Romanian labels, single bar at 142.

- [ ] **Step 5: Commit**

```bash
git add fastq_and_ngs/src/fastq_and_ngs/__main__.py fastq_and_ngs/tests/test_main.py fastq_and_ngs/data/human_tp53.fa fastq_and_ngs/results
git commit -m "Add Lab 03 CLI, FASTA copy and ZIP builder"
```

---

### Task 5: Notes, README and final build

**Files:**
- Create: `fastq_and_ngs/src/fastq_and_ngs/notes.tex`, `fastq_and_ngs/notes.pdf` (generated)
- Modify: `fastq_and_ngs/README.txt`, `README.md` (root)
- Test: `fastq_and_ngs/tests/test_docs.py`

**Interfaces:**
- Consumes: `main` from Task 4; generated `results/` files

- [ ] **Step 1: Write the failing tests**

```python
def test_notes_numbers_match_generated_report():
    notes = (SRC / "notes.tex").read_text()
    for value in ("81220", "142", "37.30"):
        assert value in QC_REPORT.read_text() and value in notes

def test_readme_has_required_content():
    text = README.read_text()
    for needed in ("Paul Hondola", "paulhondola", "SRR31709504", "PRJNA1198040", "read 1", "ClinVar", "Antigravity", "Claude", "uv run fastq-and-ngs"):
        assert needed in text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_docs.py -v`
Expected: FAIL on the missing `notes.tex` and placeholder README

- [ ] **Step 3: Write `notes.tex`**

Reuse the preamble and heading style of `sequence_alignment/src/sequence_alignment/notes.tex`, including `\pdfobjcompresslevel=0`. In Romanian:

- **QC:** 3–5 sentences condensed from section 2 of Paul's `paulhondola_notes.md` (false positives, artefacts, Phred-based variant callers), using the real numbers. Note that the uniform 142 length fits an amplicon run. Do not repeat that notes file's claims about "panel" optimality.
- **Searches:** 2–3 sentences: rsID searched directly, no-ID queried as `chr<CHROM>:<POS> AND TP53`, and PubMed's `QueryTranslation` shows it drops the position, so those hits are not specific to the variant. Take counts from `results/variants_paulhondola.txt`.

- [ ] **Step 4: Write `README.txt`**

Author and handle, Python version and dependencies, run steps (`uv sync`, `uv run fastq-and-ngs`, `--task`, `--refresh`, `--pdf --zip`, `uv run pytest -v`). FASTQ source: SRR31709504, study PRJNA1198040, read 1 of a pair, ENA link, and how to place the file, since it's not in the ZIP. VCF source: ClinVar IDs, accessed 2026-11-04. AI disclosure: Google Antigravity for the earlier QC script and ENA download; Claude for this project; how it was verified. Mark the disclosure for Paul to review.

- [ ] **Step 5: Update the root `README.md`**

Replace `(`sequence_alignment/README.txt` is the only one filled in)` so it also names `fastq_and_ngs/README.txt`.

- [ ] **Step 6: Build and verify everything**

```bash
uv run pytest -v
uv run fastq-and-ngs --pdf --zip
grep -ac '^/Type /Page$' notes.pdf
unzip -l fastq_and_ngs_paulhondola.zip
```

Expected: all tests PASS; page count `1`; the ZIP lists all deliverables, no `data/work/`, and is well under 4 MB. Read the rendered `notes.pdf` for layout and correct Romanian diacritics.

- [ ] **Step 7: Commit**

```bash
git add fastq_and_ngs/src/fastq_and_ngs/notes.tex fastq_and_ngs/notes.pdf fastq_and_ngs/README.txt fastq_and_ngs/tests/test_docs.py fastq_and_ngs/results README.md
git commit -m "Add Lab 03 notes and README"
```

Paul then reviews the notes and AI disclosure, and uploads the ZIP to Moodle himself.
