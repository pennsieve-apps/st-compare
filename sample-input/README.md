# sample-input

This app reads three inputs (the contract is in
`workflow-spatial-transcriptomics/docs/conventions.md`):

- `spatial/<sample_id>.h5ad` from `st-qc` (raw counts in `layers["counts"]`; `st-prep`'s output,
  with counts in `X`, also works),
- `labels/labels.csv` from `st-domains` (`domain`), `st-celltype` (`cluster`, `cell_type`) or
  `cell2location` (`dominant_type`),
- `metadata/samples.csv` from `st-prep`, with a condition column (and optionally a block column
  such as donor or slide) carried over from your sample sheet.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/domain/
```

Six Visium-like samples (`a1`-`a3` in condition `a`, `b1`-`b3` in `b`; `donor` pairs one of
each), 300 spots and 200 genes each, labelled `D01`-`D03` as st-domains would. Planted: `D03`
holds 20 % of the spots in `a` and 35 % in `b` (`D01` 50 % -> 35 %), and in `D02` five genes
(`G100`-`G104`) are 4x higher and three (`G105`-`G107`) 4x lower in `b`. The app finds both
composition changes and all eight genes (and nothing else) at FDR 0.05.

## Real data

The TLS Visium demo of the workflow repository (`scripts/fetch-sample-data.sh`): Dawo S,
Nonchev K, Silina K. *10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung (5)
Cancer with Tertiary Lymphoid Structures.* Zenodo 2025, https://doi.org/10.5281/zenodo.14620362
(CC BY 4.0): 8 FFPE tumour sections, `samples.csv` with `tissue` (kidney / lung).

Run `st-prep -> st-qc -> st-domains` on it, then put the three outputs side by side:

```sh
mkdir -p data/tls && cd data/tls
ln -s /path/to/st-qc/output/spatial  spatial
ln -s /path/to/st-domains/output     labels      # holds labels.csv
ln -s /path/to/st-prep/output/metadata metadata  # holds samples.csv
cd ../..
docker run --rm --platform linux/amd64 --memory 3g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -e LABEL_COLUMN=domain -e CONDITION_COLUMN=tissue \
  -v "$PWD/data/tls:/data/input:ro" -v "$PWD/data/output:/data/output" st-compare
```

(With symlinks, mount the real folders instead, or copy them; Docker does not follow links
that point outside the mount.) Or `make run INPUT=/path/to/folder` with the three subfolders.
