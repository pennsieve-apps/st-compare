"""
A synthetic fixture that follows the spatial transcriptomics conventions, with planted effects:

  6 samples (a1-a3 in condition `a`, b1-b3 in condition `b`; donor d1-d3 pairs one of each),
  3 labels (D01-D03, as st-domains writes them), 200 genes, 300 observations per sample.

  Planted: in D02 the genes DE_UP (x4) are higher and DE_DOWN (x0.25) lower in `b` than in `a`;
  D03 holds 20 % of the observations in `a` and 35 % in `b` (D01 the other way round).
  Each label also has its own marker genes, so pseudobulks of different labels differ, as they do
  in real data.

Writes INPUT_DIR layouts: linked (spatial/, labels/, metadata/) and flat (everything side by side).

    python -m tests.make_fixture data/fixture --layout both
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

SAMPLES = [("a1", "a", "d1"), ("a2", "a", "d2"), ("a3", "a", "d3"),
           ("b1", "b", "d1"), ("b2", "b", "d2"), ("b3", "b", "d3")]
LABELS = ["D01", "D02", "D03"]
COMPOSITION = {"a": [0.50, 0.30, 0.20], "b": [0.35, 0.30, 0.35]}
N_GENES = 200
N_OBS = 300
DE_LABEL = "D02"
DE_UP = ["G100", "G101", "G102", "G103", "G104"]
DE_DOWN = ["G105", "G106", "G107"]
MARKERS = {"D01": range(0, 10), "D02": range(10, 20), "D03": range(20, 30)}


# Gene baselines are shared by every sample (only the sample-level noise below differs).
BASE = np.exp(np.random.default_rng(0).normal(np.log(2.0), 1.0, N_GENES))


def gene_names() -> list[str]:
    return [f"G{i:03d}" for i in range(N_GENES)]


def make_sample(sample_id: str, condition: str, donor: str, rng: np.random.Generator):
    import anndata as ad

    genes = gene_names()
    base = BASE                                                         # mean counts per observation (shared)
    sample_factor = np.exp(rng.normal(0.0, 0.2, N_GENES))               # sample-level biological noise
    labels = rng.choice(LABELS, size=N_OBS, p=COMPOSITION[condition])
    counts = np.zeros((N_OBS, N_GENES), dtype=np.int32)
    for lab in LABELS:
        rows = np.flatnonzero(labels == lab)
        mu = base * sample_factor
        mu = mu.copy()
        mu[list(MARKERS[lab])] *= 3.0
        if lab == DE_LABEL and condition == "b":
            mu[[genes.index(g) for g in DE_UP]] *= 4.0
            mu[[genes.index(g) for g in DE_DOWN]] *= 0.25
        size = rng.gamma(0.7, 1.0, size=(len(rows), 1))                  # per-observation depth
        rate = mu[None, :] * size * rng.gamma(1 / 0.3, 0.3, size=(len(rows), N_GENES))
        counts[rows] = rng.poisson(rate).astype(np.int32)
    # Hex-like grid coordinates in "pixels".
    n_side = int(np.ceil(np.sqrt(N_OBS)))
    ij = np.array([(i, j) for i in range(n_side) for j in range(n_side)])[:N_OBS]
    xy = np.column_stack([ij[:, 1] * 100.0 + (ij[:, 0] % 2) * 50.0, ij[:, 0] * 87.0]) + 500.0
    barcodes = [f"BC{i:05d}-1" for i in range(N_OBS)]
    obs_names = [f"{sample_id}:{b}" for b in barcodes]

    X = sp.csr_matrix(counts, dtype=np.float32)
    totals = np.asarray(X.sum(axis=1)).ravel()
    norm = sp.diags(np.where(totals > 0, 1e4 / np.maximum(totals, 1), 0.0)) @ X
    norm.data = np.log1p(norm.data).astype(np.float32)
    obs = pd.DataFrame({
        "sample_id": pd.Categorical([sample_id] * N_OBS),
        "x": xy[:, 0], "y": xy[:, 1], "in_tissue": 1,
        "array_row": ij[:, 0], "array_col": ij[:, 1],
        "total_counts": totals.astype(np.float32),
        "n_genes_by_counts": np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int32),
        "pct_counts_mt": np.float32(0.0),
        "condition": condition, "donor": donor,
    }, index=obs_names)
    var = pd.DataFrame({"gene_ids": [f"ENSG{i:011d}" for i in range(N_GENES)],
                        "feature_types": "Gene Expression", "mt": False}, index=genes)
    adata = ad.AnnData(X=sp.csr_matrix(norm, dtype=np.float32), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(counts, dtype=np.int32)
    adata.obsm["spatial"] = xy.astype(np.float64)
    adata.uns["pennsieve"] = {"platform": "visium", "kind": "spot", "coordinate_unit": "pixel",
                              "spot_diameter": 90.0, "grid": "hex", "source": f"fixture/{sample_id}",
                              "prep_version": "fixture", "normalized": True, "target_sum": 10000.0}
    labels_df = pd.DataFrame({"sample_id": sample_id, "obs_id": obs_names, "domain": labels})
    return adata, labels_df


def make_fixture(root: Path, layout: str = "both", seed: int = 7) -> dict[str, Path]:
    root = Path(root)
    rng = np.random.default_rng(seed)
    adatas, label_frames, rows = [], [], []
    for sample_id, condition, donor in SAMPLES:
        adata, labels_df = make_sample(sample_id, condition, donor, rng)
        adatas.append((sample_id, adata))
        label_frames.append(labels_df)
        rows.append({"sample_id": sample_id, "path": f"fixture/{sample_id}", "platform": "visium", "kind": "spot",
                     "n_obs": adata.n_obs, "n_genes": adata.n_vars, "n_obs_in_tissue": adata.n_obs,
                     "condition": condition, "donor": donor})
    labels = pd.concat(label_frames, ignore_index=True)
    samples = pd.DataFrame(rows)
    out: dict[str, Path] = {}
    layouts = ["linked", "flat"] if layout == "both" else [layout]
    for name in layouts:
        base = root / name
        if base.exists():
            shutil.rmtree(base)
        spatial = base / "spatial" if name == "linked" else base
        labels_dir = base / "labels" if name == "linked" else base
        meta = base / "metadata" if name == "linked" else base
        for d in (spatial, labels_dir, meta):
            d.mkdir(parents=True, exist_ok=True)
        for sample_id, adata in adatas:
            adata.write_h5ad(spatial / f"{sample_id}.h5ad")
        labels.to_csv(labels_dir / "labels.csv", index=False)
        samples.to_csv(meta / "samples.csv", index=False)
        out[name] = base
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)
    paths = make_fixture(args.out, args.layout, args.seed)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
