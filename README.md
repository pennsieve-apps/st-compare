# st-compare

Which tissue regions or cell types become more or less common between conditions, and which
genes change within each of them? This step of the **Spatial Domain Mapping & Condition
Comparison** workflow answers both questions for spatial transcriptomics samples, per
**label**: a spatial domain from `st-domains`, a cell type or cluster from `st-celltype`, or
the dominant cell type from `cell2location`.

- **Composition (differential abundance).** Does a label make up a larger or smaller share
  of the tissue in one condition than in the reference? E.g. tertiary lymphoid structures
  covering more of the section in lung than in kidney tumours.
- **Pseudobulk differential expression (DE).** Within a label, which genes are expressed at a
  different level between conditions? E.g. the tumour regions of lung vs kidney cancers, or
  TLS in one condition vs another. Each sample's spots (or bins, or cells) of that label are
  summed into one **pseudobulk** profile, and the profiles are compared with DESeq2's
  negative binomial model (PyDESeq2).

Both are tests across **samples**, not spots or cells: each section is one replicate, the
variability between sections within a condition is what a change is judged against, and a
donor or batch can be accounted for. That is the point of pseudobulk over per-spot tests,
which treat thousands of spots from one section as independent replicates and call many
false discoveries (Squair et al. 2021). It also means the statistics are only as good as the
number of samples: with 2-3 sections per condition only large, consistent changes can be
found.

`workflow A: st-prep -> st-qc -> st-domains -> {st-svg, st-compare} -> st-report`

The app is light (tables of samples x labels and pseudobulks x genes; the counts are read a
block of rows at a time), so it runs as a container or a Lambda function.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `spatial` | input (`placement: subfolder`, required) | `<sample_id>.h5ad` | From `st-qc`: raw counts in `layers["counts"]` (used when present, else `X`, which must then hold integer counts, as `st-prep`'s output does). Only `obs_names`, `var_names` and the count matrix are read. |
| `labels` | input (`placement: subfolder`, required) | `labels.csv` (`sample_id`, `obs_id`, label columns) | From `st-domains` (`domain`), `st-celltype` (`cluster`, `cell_type`) or `cell2location` (`dominant_type`). `obs_id` matches the h5ad's `obs_names`. |
| `metadata` | input (`placement: subfolder`, required) | `samples.csv` | From `st-prep`: `sample_id` plus the columns of your sample sheet (`condition`, `tissue`, `donor`, ...). |
| `compare` | output | `<LABEL>/...` (below) | Everything under a folder named after `LABEL`, so several runs (domains, cell types) can be linked into one port of `st-report`. |

Linked from a workflow, the inputs arrive in `INPUT_DIR/spatial/`, `INPUT_DIR/labels/` and
`INPUT_DIR/metadata/`. Without links (everything side by side under `INPUT_DIR`, possibly
nested, entries symlinked) the files are found by name anywhere under `INPUT_DIR`.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `LABEL` | string | *(empty)* | Output folder name (`[a-z0-9-]+`). Empty: `LABEL_COLUMN` with `_` as `-` (`domain`, `cell-type`). |
| `LABEL_COLUMN` | string | `domain` | Column of `labels.csv` to compare. |
| `CONDITION_COLUMN` | string | `condition` | Column of `samples.csv` with each sample's condition. For the TLS demo: `tissue`. |
| `REFERENCE_LEVEL` | string | *(empty)* | The level the others are compared with: log2 fold change > 0 and CLR difference > 0 mean higher than in the reference. Empty: the first level alphabetically. |
| `BLOCK_COLUMN` | string | *(empty)* | Column of `samples.csv` with a blocking factor (donor, patient, slide, batch). Added to the DE design as `~ block + condition` (a fixed effect). Dropped for a label, with a note in `de_summary.csv`, when it has one level among that label's samples or is confounded with the condition (design not of full rank, or no residual degrees of freedom). The composition tests do not use it (the boxplots show it as point shape). |
| `MIN_OBS_PER_PSEUDOBULK` | integer | `20` | A sample contributes a pseudobulk to a label only with at least this many observations of it. |
| `MIN_SAMPLES_PER_GROUP` | integer | `2` | A label is tested for DE only when the reference and at least one other level each have this many pseudobulks; levels with fewer are left out for that label (noted). At least 2. |
| `MIN_COUNT_GENE` | integer | `10` | A gene is tested in a label when it has at least this many counts in at least as many of the label's pseudobulks as the smallest condition group has (edgeR's `filterByExpr` idea). |
| `FDR` | number | `0.05` | Benjamini-Hochberg adjusted p-value below which a result is `significant`. |
| `MAX_PLOTS` | integer | `8` | Volcano plots for at most this many label x comparison tests (most hits first). |
| `SEED` | integer | `42` | Seeds numpy; the analysis has no random step, so results are deterministic regardless. |
| `THREADS` | integer | `0` | Threads for PyDESeq2 (joblib) and BLAS. `0`: from the cgroup CPU quota; 2 when there is none (Lambda). |

**Samples.** A sample is used when it is in `samples.csv` with a condition (empty and `NA`
count as none), a block value when `BLOCK_COLUMN` is set, and labelled observations in
`labels.csv`; it enters the pseudobulk analysis when its `<sample_id>.h5ad` is found and
readable. Each sample's status and the reason it was left out are in `provenance.json` and the
log (one line per sample).

**More than two condition levels.** The composition tests become one-way ANOVA and
Kruskal-Wallis (is there any difference between the levels?). For DE, one model is fitted per
label to all levels with enough pseudobulks and each non-reference level is tested against the
reference with its own Wald contrast (`<level> vs <reference>`); the DE tables then have one
block of rows per comparison. Comparisons between two non-reference levels are not run: change
`REFERENCE_LEVEL` for those.

**One condition level or too few samples.** The composition table is still written
(descriptive), the tests are left empty, the DE tables have only a header, and
`summary.json` (`notes`) and the log say why. The run succeeds (exit 0).

Exit codes: `0` success; `1` no usable input (no `labels.csv`, `samples.csv` or `.h5ad`, no
sample with a condition and labels, or every `.h5ad` failed); `2` bad parameters
(`LABEL_COLUMN`, `CONDITION_COLUMN`, `BLOCK_COLUMN` not a column, `REFERENCE_LEVEL` not a
level, a value out of range).

## Outputs (`OUTPUT_DIR/<LABEL>/`)

`composition.csv`: one row per sample and label.

| Column | |
|---|---|
| `sample_id`, `condition`, `block` | `block` only with `BLOCK_COLUMN` |
| `label` | the label value |
| `n_obs`, `n_total` | the label's observations in the sample; the sample's labelled observations |
| `fraction` | `n_obs / n_total` |
| `clr` | centred log-ratio: `log(n_obs + 0.5)` minus the mean of that over the sample's labels |

`composition_tests.csv`: one row per label, sorted by p-value.

| Column | |
|---|---|
| `label`, `comparison` | `<level> vs <reference>` (two levels) or `any difference between ...` |
| `n_samples`, `n_<level>` | samples in total and per level |
| `mean_fraction_<level>`, `mean_clr_<level>` | per level |
| `clr_difference` | two levels: mean CLR of the level minus the reference (a log-ratio; > 0 = more common than in the reference, relative to the other labels) |
| `parametric_test`, `parametric_stat`, `parametric_p`, `parametric_p_adj` | Welch's t (two levels) or one-way ANOVA F, on CLR; BH across labels |
| `rank_test`, `rank_stat`, `rank_p`, `rank_p_adj` | Mann-Whitney U or Kruskal-Wallis H, on CLR; BH across labels |
| `tested`, `significant` | every level had >= 2 samples and the values were not all equal; `parametric_p_adj < FDR` |

`de.csv`: one row per label, comparison and tested gene.

| Column | |
|---|---|
| `label`, `comparison`, `gene` | |
| `base_mean` | mean of the size-factor-normalized counts over the label's pseudobulks |
| `log2_fold_change`, `lfc_se`, `stat` | DESeq2's MLE coefficient (not shrunken), its standard error, the Wald statistic |
| `p_val`, `p_adj` | Wald p-value; BH-adjusted within the label x comparison after DESeq2's independent filtering (empty for genes filtered out or flagged by Cook's distance) |
| `significant` | `p_adj < FDR` |
| `n_samples_level`, `n_samples_reference` | pseudobulks in the comparison |

`de_summary.csv`: one row per label and comparison: `status` (`tested`, `skipped`, `failed`),
`reason`, `design` (`~ condition` or `~ block + condition`), `note` (block dropped, levels
left out), `n_samples_<level>`, `n_genes_tested`, `n_significant`, `n_up`, `n_down`,
`top_up`, `top_down` (up to five genes by p-value).

`pseudobulk_counts.csv.gz`: genes x pseudobulks (`<sample_id>:<label>`), the summed raw
counts of every pseudobulk with at least `MIN_OBS_PER_PSEUDOBULK` observations (genes present
in every sample). `pseudobulk_samples.csv`: one row per sample and label: `pseudobulk`,
`sample_id`, `label`, `condition`, (`block`), `n_obs`, `total_counts`, `used`.

Figures: `composition.png` (a panel per label, % of each sample's observations per level,
one point per sample, point shape = block, the adjusted p-values in the title);
`volcano_<label>.png` (`volcano_<label>_<level>.png` with more than two levels; log2 fold
change vs -log10 p, significant genes coloured and the top ten named) for at most
`MAX_PLOTS` tests; `de_heatmap.png` (the top five genes of every test, up to 50, x tests,
coloured by log2 fold change, dots where significant).

`summary.json`: label column, condition, block, reference, levels and samples per level, the
labels; composition (tests, number tested and significant, the significant labels with CLR
difference, adjusted p-values and mean fractions); DE (run or not, labels tested / skipped,
genes, significant genes in total, per label x comparison status and hits); `notes`;
figures. `provenance.json`: app name and version, tool versions, parameters, resolved
reference / levels / threads, method, input files (relative to `INPUT_DIR`), timings, peak
memory, per-sample status (counts read from `layers/counts` or `X`, observations matched to
labels, pseudobulks).

## Methods in brief

**Composition.** Counts per sample x label come from `labels.csv`. Proportions are
compositional (one label growing makes every other label's share shrink), so the tests are on
the centred log-ratio, CLR = log(count + 0.5) - mean over labels of log(count + 0.5)
(Aitchison 1982), one value per sample. With two levels: Welch's t-test and the Mann-Whitney U
test; with more: one-way ANOVA and Kruskal-Wallis; BH across labels per test (Benjamini &
Hochberg 1995). `significant` follows the parametric test; the rank test is there as a check
that does not assume normality, but with 3 vs 3 samples its smallest possible two-sided p is
0.1, and with 3 vs 5 it is 0.036, so it cannot reach significance in small designs.
**scCODA** (Büttner et al. 2021), a Bayesian multinomial-Dirichlet model with a reference
label, is the fuller treatment of compositional changes and handles few samples better; this
app does not run it (it needs TensorFlow, which is too heavy for a Lambda-sized step), and
its use is recommended for a definitive differential-abundance analysis.

**Pseudobulk DE.** For each sample and label, the raw counts of the label's observations are
summed with h5py, a block of 20,000 rows at a time (so a Visium HD sample never sits in memory
whole). Pseudobulks with fewer than `MIN_OBS_PER_PSEUDOBULK` observations are dropped; genes
are restricted to those present in every sample and then filtered per label
(`MIN_COUNT_GENE`). Per label, PyDESeq2 (Muzellec et al. 2023), a Python reimplementation of
DESeq2 (Love et al. 2014), estimates size factors (median of ratios), gene-wise dispersions
with a parametric trend and MAP shrinkage, fits the negative binomial GLM `~ condition` or
`~ block + condition`, refits Cook's-distance outliers when a level has >= 7 replicates, and
tests each level against the reference with a Wald test; DESeq2's independent filtering and
BH adjustment give `p_adj` within each label x comparison (not across labels). Levels are
encoded internally (`L00` = reference, `L01`, ...), so any level name works.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python   # pytest on the synthetic fixture, both input layouts
make run PYTHON=.venv/bin/python    # builds the image (amd64), runs it over the fixture into data/output/
make run INPUT=/path/to/folder      # or over a folder with spatial/, labels/, metadata/ (see sample-input/README.md)
```

The tests check that the planted composition change (one label 20 % -> 35 %) and the eight
planted DE genes are found with no other hits, that a pseudobulk equals the sum of the counts
layer, the block design, a block confounded with the condition, `REFERENCE_LEVEL`, three
levels, one level, too few samples, `MIN_OBS_PER_PSEUDOBULK` skipping a label, bad parameters,
missing inputs, a log-normalized `X` without counts, and the Lambda handler.
`tests/labels_from_obs.py` writes a `labels.csv` from an `obs` column (e.g. the pathologist
annotations `st-prep` joins onto `obs["annotation"]`), to run the app without a labelling step.

## Runtime

`computeTypes: [standard, lambda]`, 2 vCPU / 3008 MB, 15 min timeout. One image serves both:
`processor/entry.py` fixes the thread count and runs `main()` as a container, or, when
`AWS_LAMBDA_RUNTIME_API` is set, starts the Lambda runtime interface client with
`processor/handler.py`, which turns the payload's `inputDir`, `outputDir` and parameters into
the same environment variables and runs the same `main()` (non-zero exit -> the invocation
raises, so the step fails). Writes go only to `OUTPUT_DIR` and `/tmp`; checked with
`docker run --read-only --tmpfs /tmp` through the handler.

Measured with `docker run --platform linux/amd64 --memory 3g --cpus 2 --network none` (amd64
emulation on an M-series Mac), peak RSS from `getrusage` and cgroup `memory.peak` (which
includes page cache from reading the `.h5ad` files):

| Data | labels tested | wall time | peak RSS | cgroup peak |
|---|---|---|---|---|
| TLS Visium demo, st-qc + st-domains output: 8 samples, 27,497 spots, 15,392 genes in every sample, 30 domains | 10 of 30 | 134 s (PyDESeq2 127 s) | 681 MB | 1.07 GB |
| TLS Visium demo, st-qc output with the pathologist annotations as labels (6 labels) | 5 of 6 | 88 s | 680 MB | 1.06 GB |
| Synthetic Visium HD scale: 4 samples x 500,000 bins x 18,000 genes (50 M counts each), 10 labels | 10 | 317 s (reading 30 s, PyDESeq2 275 s) | 686 MB | 2.88 GB (page cache of 3.1 GB of input) |
| Synthetic fixture: 6 samples x 300 spots x 200 genes | 3 | 8 s | 295 MB | |

So 2 vCPU / 3008 MB (Lambda's ceiling) has ample room: the process itself stays under 0.7 GB
because the counts are streamed and only labels x genes sums are kept.

Time is dominated by PyDESeq2: about 10-40 s per label under emulation, growing with the
number of genes that pass the filter and the number of pseudobulks; memory barely depends on
the number of observations. Many labels (e.g. 20+ cell types) on a large gene set can approach
Lambda's 15-minute limit; run those on `standard`. Image: 161 MB compressed, 687 MB unpacked (no scanpy, numba or torch).

## Results on the demo data

TLS Visium demo (Dawo et al. 2025; 3 kidney and 5 lung tumour sections), `st-prep -> st-qc ->
st-domains`, then `LABEL_COLUMN=domain CONDITION_COLUMN=tissue` (reference `kidney`, the
first alphabetically), defaults otherwise. Output in the workflow repository under
`data/work/tls/compare/domain/`.

- **Composition: 0 of 30 domains change significantly** (FDR 0.05). The strongest
  candidates are kidney-specific domains (D13: 7.4 % of kidney spots vs ~0.2 % of lung,
  p = 0.006, q = 0.19; D08 and D14 also near-absent in lung) and D05, more common in lung
  (8.7 % vs 5.2 %, p = 0.037, q = 0.55). Many of st-domains' 30 domains occur in only one or
  two sections (sample- or organ-specific), which with 3 vs 5 samples cannot reach
  significance after correction across 30 labels.
- **Pseudobulk DE: 10 of 30 domains tested** (20 skipped: fewer than 2 sections with >= 20
  spots of the domain in kidney or in lung), **1,663 significant gene x domain results**.
  Every tested domain shows the organ: lung surfactant genes (SFTPB, SFTPC, SLC34A2) and
  immunoglobulins (IGHG1/3/4, IGHD) up in lung; kidney tubule genes (FXYD2, HAVCR1, CDH16,
  KCNJ16, DEFB1, BHMT2) and CRYAB down. D02, the domain that best matches the pathologists'
  TLS annotation (F1 0.52 in st-domains' `annotation_match.csv`), has 125 hits, with CCL18
  and FASN up in lung TLS.
- With the **pathologist annotations** as labels instead (`tests/labels_from_obs.py` on
  st-qc's output; `data/work/tls/compare/annotation/`): no composition change (TLS 2.6 % of
  kidney vs 4.2 % of lung spots, q = 0.88); DE in 5 of 6 annotations (Lymph_node only in one
  lung section), 2,124 significant results: Tumor 1,137, Unassigned 558, Normal 383, TLS 31,
  Immune 15 (only 2 kidney sections with >= 20 Immune spots).

As the limitations below say, kidney vs lung here is organ, tumour type and patient at once;
the DE mostly recovers organ identity, which is a sanity check rather than a finding.

## Limitations

- **Few samples, little power.** The sample is the replicate. With 3 vs 5 sections, only large
  and consistent changes are detectable, and the rank tests cannot reach significance at 3 vs 3.
  A "not significant" is not evidence of no change.
- **Condition vs everything that goes with it.** In the demo, kidney and lung differ in organ,
  tumour type and patient; DE between them mixes all three. A block (donor, slide) is a fixed
  effect and only helps when it is not confounded with the condition.
- **Labels are inputs, not truth.** Domains or cell types found on all samples together can
  absorb condition differences (a domain only present in one condition); look at
  `composition.csv` before reading DE for a label that is rare in one level. Spots on Visium
  hold several cells, so a domain's pseudobulk mixes cell types and a DE gene can reflect a
  change in cell-type mix within the domain rather than in expression.
- **Composition tests are simple.** CLR + t / rank tests per label, BH across labels; no
  model of the multinomial sampling (scCODA does that), no block in the composition tests.
- **Log fold changes are not shrunken** (PyDESeq2's `lfc_shrink` is not applied); for ranking
  genes by effect size, weigh `lfc_se`. BH is within each label x comparison, not across labels.
- Genes missing from some samples' `.h5ad` are dropped from DE (only genes in every sample are
  used); samples processed with different gene filters lose genes.
- Only each-level-vs-reference contrasts; no interactions or continuous covariates.

## Licences

Our code: Apache-2.0 (`LICENSE`). Everything in the image is under permissive licences
(PyDESeq2 MIT; anndata, h5py, NumPy, SciPy, pandas, scikit-learn BSD-3-Clause; Matplotlib's
PSF-based licence); nothing is GPL. See `THIRD_PARTY.md`.

## Citations

Checked against Crossref (title, authors, journal, year, volume, pages) and PubMed
(E-utilities `esearch` on the DOI, then `esummary`: title, year, first author, author count);
the demo data via DataCite. Full author lists are in `CITATION.cff`.

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq
  data with DESeq2. *Genome Biol* 15(12):550 (2014). doi:10.1186/s13059-014-0550-8.
  PMID 25516281.
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq
  differential expression analysis. *Bioinformatics* 39(9):btad547 (2023).
  doi:10.1093/bioinformatics/btad547. PMID 37669147.
- Squair JW, Gautier M, Kathe C, Anderson MA, James ND, Hutson TH, Hudelle R, Qaiser T,
  Matson KJE, Barraud Q, Levine AJ, La Manno G, Skinnider MA, Courtine G. Confronting false
  discoveries in single-cell differential expression. *Nat Commun* 12(1):5692 (2021).
  doi:10.1038/s41467-021-25960-2. PMID 34584091. (why pseudobulk)
- Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B. scCODA is a Bayesian model for
  compositional single-cell data analysis. *Nat Commun* 12(1):6876 (2021).
  doi:10.1038/s41467-021-27150-6. PMID 34824236. (the fuller composition model; not run)
- Aitchison J. The Statistical Analysis of Compositional Data. *J R Stat Soc Series B*
  44(2):139-160 (1982). doi:10.1111/j.2517-6161.1982.tb01195.x. (CLR; not in PubMed, Crossref only)
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful
  Approach to Multiple Testing. *J R Stat Soc Series B* 57(1):289-300 (1995).
  doi:10.1111/j.2517-6161.1995.tb02031.x. (not in PubMed, Crossref only)
- Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific
  computing in Python. *Nat Methods* 17(3):261-272 (2020). doi:10.1038/s41592-019-0686-2.
  PMID 32015543. (the t, U, F and H tests)
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated
  data matrices. *J Open Source Softw* 9(101):4371 (2024). doi:10.21105/joss.04371. (not in
  PubMed, Crossref only)

Demo data: Dawo S, Nonchev K, Silina K. 10x Visium Spatial Transcriptomics Dataset: Kidney (3)
and Lung (5) Cancer with Tertiary Lymphoid Structures. Zenodo (2025).
doi:10.5281/zenodo.14620362. CC BY 4.0.
