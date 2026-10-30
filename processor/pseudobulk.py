"""
Pseudobulk differential expression between condition levels, per label.

For every sample and label, the raw counts (layers["counts"], else X when it
holds integer counts) of that label's observations are summed into one
pseudobulk profile; a pseudobulk needs at least MIN_OBS_PER_PSEUDOBULK
observations. Per label, PyDESeq2 (the Python implementation of DESeq2) fits a
negative binomial GLM across the samples' pseudobulks with `~ condition` or,
with a block, `~ block + condition`, and each non-reference level is tested
against the reference with a Wald test. The sample is the replicate, as in
Squair et al. 2021: tests on pseudobulks avoid the false discoveries that
per-spot or per-cell tests make by treating thousands of observations from one
section as independent.

The counts are read with h5py, a block of rows at a time, so a large sample
(Visium HD, a Xenium slide) never sits in memory as a whole.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("st-compare")
CHUNK_ROWS = 20000


@dataclass
class SamplePseudobulk:
    sample_id: str
    genes: list[str]
    labels: list[str]
    sums: np.ndarray            # labels x genes, float64 (integer-valued)
    n_obs: dict[str, int]       # label -> observations summed
    source: str                 # "layers/counts" or "X"
    n_obs_file: int = 0
    n_obs_labelled: int = 0


def _read_index(group) -> np.ndarray:
    from anndata.io import read_elem

    key = group.attrs.get("_index", "_index")
    return np.asarray(read_elem(group[key])).astype(str)


def _matrix_elem(f):
    if "layers" in f and "counts" in f["layers"]:
        return f["layers"]["counts"], "layers/counts"
    return f["X"], "X"


def _iter_rows(elem, n_obs: int, n_vars: int, chunk: int):
    """Yield (start, stop, csr block) over the rows of a CSR group, a dense dataset or (read whole) anything else."""
    import h5py

    enc = elem.attrs.get("encoding-type", "")
    if isinstance(elem, h5py.Group) and enc == "csr_matrix":
        indptr = elem["indptr"][:]
        for start in range(0, n_obs, chunk):
            stop = min(start + chunk, n_obs)
            a, b = int(indptr[start]), int(indptr[stop])
            block = sp.csr_matrix((elem["data"][a:b], elem["indices"][a:b], indptr[start:stop + 1] - a),
                                  shape=(stop - start, n_vars))
            yield start, stop, block
    elif isinstance(elem, h5py.Dataset):
        for start in range(0, n_obs, chunk):
            stop = min(start + chunk, n_obs)
            yield start, stop, sp.csr_matrix(elem[start:stop])
    else:
        from anndata.io import read_elem

        m = sp.csr_matrix(read_elem(elem))
        for start in range(0, n_obs, chunk):
            stop = min(start + chunk, n_obs)
            yield start, stop, m[start:stop]


def sample_pseudobulk(path: Path, sample_id: str, obs_labels: pd.Series, chunk: int = CHUNK_ROWS) -> SamplePseudobulk:
    """Sum the counts of each label's observations. obs_labels: obs_id -> label for this sample."""
    import h5py

    with h5py.File(path, "r") as f:
        obs_names = _read_index(f["obs"])
        genes = _read_index(f["var"]).tolist()
        elem, source = _matrix_elem(f)
        lab = obs_labels.reindex(obs_names)
        labels = sorted(pd.unique(lab.dropna()).tolist())
        code = {v: i for i, v in enumerate(labels)}
        codes = np.array([code.get(v, -1) if isinstance(v, str) else -1 for v in lab.to_numpy()], dtype=np.int64)
        sums = np.zeros((len(labels), len(genes)), dtype=np.float64)
        for start, stop, block in _iter_rows(elem, len(obs_names), len(genes), chunk):
            data = block.data
            if data.size and not np.all(np.equal(np.mod(data, 1), 0)):
                raise ValueError(f"{source} does not hold integer counts (log-normalized X and no layers['counts']?)")
            if data.size and data.min() < 0:
                raise ValueError(f"{source} has negative values; not counts")
            c = codes[start:stop]
            keep = np.flatnonzero(c >= 0)
            if keep.size == 0:
                continue
            onehot = sp.csr_matrix((np.ones(keep.size), (c[keep], keep)), shape=(len(labels), stop - start))
            sums += np.asarray((onehot @ block.astype(np.float64)).todense())
    n_obs = {v: int((codes == i).sum()) for i, v in enumerate(labels)}
    return SamplePseudobulk(sample_id, genes, labels, sums, n_obs, source, len(obs_names), int((codes >= 0).sum()))


@dataclass
class LabelResult:
    label: str
    status: str                                        # tested | skipped | failed
    reason: str = ""
    design: str = ""
    block_note: str = ""
    levels: list[str] = field(default_factory=list)
    n_per_level: dict[str, int] = field(default_factory=dict)
    n_genes_tested: int = 0
    results: list[pd.DataFrame] = field(default_factory=list)   # one per comparison


def _full_rank(meta: pd.DataFrame, with_block: bool) -> tuple[bool, str]:
    cols = [pd.Series(1.0, index=meta.index, name="intercept")]
    factors = (["block"] if with_block else []) + ["condition"]
    for fct in factors:
        cols.append(pd.get_dummies(meta[fct], drop_first=True, dtype=float))
    X = pd.concat(cols, axis=1).to_numpy()
    rank = np.linalg.matrix_rank(X)
    if rank < X.shape[1]:
        return False, "confounded with the condition (design not of full rank)"
    if X.shape[0] - rank < 1:
        return False, "leaves no residual degrees of freedom"
    return True, ""


def filter_genes(counts: pd.DataFrame, min_count: int, min_samples: int) -> pd.DataFrame:
    """counts: pseudobulks x genes. Keep genes with >= min_count counts in >= min_samples pseudobulks
    (edgeR filterByExpr's idea, with min_samples = the smallest group's size)."""
    keep = (counts >= min_count).sum(axis=0) >= max(1, min_samples)
    return counts.loc[:, keep & (counts.sum(axis=0) > 0)]


def de_for_label(label: str, counts: pd.DataFrame, meta: pd.DataFrame, reference: str, level_order: list[str],
                 use_block: bool, params: dict, n_cpus: int) -> LabelResult:
    """counts: pseudobulks x genes (index = sample_id); meta: index sample_id, columns condition [, block].
    Levels with fewer than MIN_SAMPLES_PER_GROUP pseudobulks are left out for this label."""
    min_s = params["MIN_SAMPLES_PER_GROUP"]
    n_per = meta["condition"].value_counts().to_dict()
    res = LabelResult(label, "skipped", n_per_level={lv: int(n_per.get(lv, 0)) for lv in level_order})
    levels = [lv for lv in level_order if n_per.get(lv, 0) >= min_s]
    if reference not in levels:
        res.reason = (f"reference level {reference!r} has {n_per.get(reference, 0)} pseudobulks "
                      f"(< MIN_SAMPLES_PER_GROUP={min_s})")
        return res
    if len(levels) < 2:
        res.reason = f"no other level has >= MIN_SAMPLES_PER_GROUP={min_s} pseudobulks"
        return res
    dropped = [lv for lv in level_order if lv not in levels and n_per.get(lv, 0) > 0]
    meta = meta[meta["condition"].isin(levels)].copy()
    counts = counts.loc[meta.index]
    res.levels = levels

    # Internal tokens keep formulaic away from odd level names; reference first so it is the base level.
    cond_tok = {lv: f"L{i:02d}" for i, lv in enumerate(levels)}
    md = pd.DataFrame({"condition": meta["condition"].map(cond_tok)}, index=meta.index)
    design = "~ condition"
    notes = []
    if use_block:
        blocks = sorted(meta["block"].unique())
        if len(blocks) < 2:
            notes.append(f"block dropped for this label: one level ({blocks[0]})")
        else:
            md["block"] = meta["block"].map({b: f"B{i:02d}" for i, b in enumerate(blocks)})
            ok, why = _full_rank(md, True)
            if ok:
                design = "~ block + condition"
            else:
                md = md.drop(columns="block")
                notes.append(f"block dropped for this label: {why}")
    ok, why = _full_rank(md, False)
    if not ok:
        res.reason = f"condition design {why}"
        return res
    if dropped:
        notes.append("levels left out (too few pseudobulks): " + ", ".join(dropped))
    res.design = design
    res.block_note = "; ".join(notes)

    smallest = int(md["condition"].value_counts().min())
    filt = filter_genes(counts, params["MIN_COUNT_GENE"], smallest)
    res.n_genes_tested = int(filt.shape[1])
    if filt.shape[1] < 2:
        res.reason = f"{filt.shape[1]} genes pass the MIN_COUNT_GENE filter"
        return res

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dds = DeseqDataSet(counts=filt.round().astype(np.int64), metadata=md, design=design,
                               refit_cooks=True, n_cpus=n_cpus, quiet=True)
            dds.deseq2()
            for lv in levels[1:]:
                ds = DeseqStats(dds, contrast=["condition", cond_tok[lv], cond_tok[reference]],
                                alpha=params["FDR"], quiet=True, n_cpus=n_cpus)
                ds.summary()
                r = ds.results_df.copy()
                r.index.name = "gene"
                r = r.reset_index().rename(columns={"baseMean": "base_mean", "log2FoldChange": "log2_fold_change",
                                                    "lfcSE": "lfc_se", "pvalue": "p_val", "padj": "p_adj"})
                r.insert(0, "label", label)
                r.insert(1, "comparison", f"{lv} vs {reference}")
                r["n_samples_level"] = int((md["condition"] == cond_tok[lv]).sum())
                r["n_samples_reference"] = int((md["condition"] == cond_tok[reference]).sum())
                res.results.append(r)
    except Exception as exc:  # one label failing is reported, not fatal
        res.status, res.reason = "failed", f"{type(exc).__name__}: {exc}"
        return res
    res.status = "tested"
    return res
