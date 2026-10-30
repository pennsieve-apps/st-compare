"""
st-compare: which labels (spatial domains, cell types) change in abundance
between conditions, and which genes change within each label?

Reads st-qc's <sample_id>.h5ad files (port `spatial`), a labels.csv (port
`labels`, from st-domains, st-celltype or cell2location) and st-prep's
samples.csv (port `metadata`), and writes, all under OUTPUT_DIR/<LABEL>/ so
several runs (domains, cell types) can be linked into one report folder:

  composition.csv          sample x label: n_obs, fraction, CLR
  composition_tests.csv    per label: Welch t / Mann-Whitney (2 levels) or ANOVA / Kruskal-Wallis, BH
  pseudobulk_counts.csv.gz genes x pseudobulks (sample:label), summed raw counts
  pseudobulk_samples.csv   one row per pseudobulk: sample, label, condition, n_obs, used
  de.csv                   per label x comparison x gene: PyDESeq2 Wald test
  de_summary.csv           per label x comparison: status, design, samples, genes, hits
  composition.png, volcano_*.png, de_heatmap.png, summary.json, provenance.json

See processor/composition.py and processor/pseudobulk.py for the methods.
"""
from __future__ import annotations

import json
import logging
import math
import os
import platform
import re
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("st-compare")
log.setLevel(logging.INFO)   # the Lambda runtime sets the root logger to WARNING first
for noisy in ("matplotlib", "PIL", "fontTools", "h5py", "numba"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-compare"
REPO_ROOT = Path(__file__).resolve().parent.parent
LABEL_RE = re.compile(r"^[a-z0-9-]+$")

DEFAULTS = {
    "LABEL": "",
    "LABEL_COLUMN": "domain",
    "CONDITION_COLUMN": "condition",
    "REFERENCE_LEVEL": "",
    "BLOCK_COLUMN": "",
    "MIN_OBS_PER_PSEUDOBULK": "20",
    "MIN_SAMPLES_PER_GROUP": "2",
    "MIN_COUNT_GENE": "10",
    "FDR": "0.05",
    "MAX_PLOTS": "8",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (str(env.get(k)).strip() if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {k: raw[k] for k in ("LABEL_COLUMN", "CONDITION_COLUMN", "REFERENCE_LEVEL", "BLOCK_COLUMN")}
    if p["BLOCK_COLUMN"].lower() == "none":
        p["BLOCK_COLUMN"] = ""
    if not p["LABEL_COLUMN"] or not p["CONDITION_COLUMN"]:
        raise ParameterError("LABEL_COLUMN and CONDITION_COLUMN must not be empty")
    label = raw["LABEL"] or re.sub(r"[^a-z0-9-]+", "-", p["LABEL_COLUMN"].lower().replace("_", "-")).strip("-")
    if not LABEL_RE.match(label):
        raise ParameterError(f"LABEL must match [a-z0-9-]+, got {label!r}")
    p["LABEL"] = label
    try:
        for k in ("MIN_OBS_PER_PSEUDOBULK", "MIN_SAMPLES_PER_GROUP", "MIN_COUNT_GENE", "MAX_PLOTS", "SEED", "THREADS"):
            p[k] = int(raw[k])
        p["FDR"] = float(raw["FDR"])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if min(p["MIN_OBS_PER_PSEUDOBULK"], p["MIN_COUNT_GENE"], p["MAX_PLOTS"], p["THREADS"]) < 0:
        raise ParameterError("MIN_OBS_PER_PSEUDOBULK, MIN_COUNT_GENE, MAX_PLOTS and THREADS must be >= 0")
    if p["MIN_SAMPLES_PER_GROUP"] < 2:
        raise ParameterError("MIN_SAMPLES_PER_GROUP must be >= 2 (a dispersion needs replicates)")
    if not 0 < p["FDR"] < 1:
        raise ParameterError("FDR must be in (0, 1)")
    return p


def _rel(path, root: Path):
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _clean_json(obj):
    import numpy as np

    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        obj = float(obj)
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def write_csv(df, path: Path) -> None:
    df.to_csv(path, index=False, float_format="%.6g", lineterminator="\n")


def tool_versions() -> dict:
    out = {"python": platform.python_version()}
    for pkg in ("pydeseq2", "anndata", "h5py", "numpy", "scipy", "pandas", "scikit-learn", "formulaic",
                "formulaic-contrasts", "matplotlib"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def peak_rss_mb():
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def cgroup_peak_mb():
    try:
        return round(int(open("/sys/fs/cgroup/memory.peak").read()) / 1024 / 1024, 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def _text(v) -> str:
    import pandas as pd

    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    s = str(v).strip()
    return "" if s.upper() in ("NA", "NAN", "NONE") else s


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    n_cpus = int(os.environ.get("ST_COMPARE_THREADS") or params["THREADS"] or 1)

    import numpy as np
    import pandas as pd

    from processor import composition as C
    from processor import figures as F
    from processor import pseudobulk as P
    from processor.inputs import discover

    np.random.seed(params["SEED"])
    found = discover(input_dir)
    if found.labels_csv is None or found.samples_csv is None:
        log.error("need labels.csv (port labels) and samples.csv (port metadata) under %s; found labels=%s samples=%s",
                  input_dir, found.labels_csv, found.samples_csv)
        return 1
    if not found.h5ad:
        log.error("no <sample_id>.h5ad files (port spatial) under %s", input_dir)
        return 1
    labels = pd.read_csv(found.labels_csv, dtype=str, keep_default_na=False)
    lc, cc, bc = params["LABEL_COLUMN"], params["CONDITION_COLUMN"], params["BLOCK_COLUMN"]
    missing = [c for c in ("sample_id", "obs_id") if c not in labels.columns]
    if missing:
        log.error("labels.csv lacks %s", ", ".join(missing))
        return 1
    if lc not in labels.columns:
        log.error("LABEL_COLUMN %r is not a column of labels.csv (columns: %s)", lc, ", ".join(labels.columns))
        return 2
    samples = pd.read_csv(found.samples_csv, dtype=str, keep_default_na=False)
    if "sample_id" not in samples.columns:
        log.error("samples.csv has no sample_id column")
        return 1
    for col, what in ((cc, "CONDITION_COLUMN"), (bc, "BLOCK_COLUMN")):
        if col and col not in samples.columns:
            log.error("%s %r is not a column of samples.csv (columns: %s)", what, col, ", ".join(samples.columns))
            return 2

    labels[lc] = labels[lc].map(_text)
    labels = labels[labels[lc] != ""]
    label_samples = set(labels["sample_id"])
    statuses: dict[str, dict] = {}
    meta_rows = []
    for _, r in samples.drop_duplicates("sample_id").iterrows():
        sid = _text(r["sample_id"])
        if not sid:
            continue
        st = {"sample_id": sid, "condition": _text(r[cc]), "block": _text(r[bc]) if bc else None,
              "h5ad": _rel(found.h5ad.get(sid), input_dir)}
        if not st["condition"]:
            st.update(status="excluded", reason=f"no {cc} value")
        elif bc and not st["block"]:
            st.update(status="excluded", reason=f"no {bc} value")
        elif sid not in label_samples:
            st.update(status="excluded", reason=f"no observation with a {lc} in labels.csv")
        else:
            st["status"] = "ok"
            meta_rows.append({"sample_id": sid, "condition": st["condition"], **({"block": st["block"]} if bc else {})})
        statuses[sid] = st
    for sid in sorted(label_samples - set(statuses)):
        statuses[sid] = {"sample_id": sid, "status": "excluded", "reason": "in labels.csv but not in samples.csv"}
    for sid, st in statuses.items():
        if st["status"] != "ok":
            log.warning("%s: excluded, %s", sid, st["reason"])
    if not meta_rows:
        log.error("no sample has both a %s in samples.csv and labels in labels.csv", cc)
        return 1
    meta = pd.DataFrame(meta_rows).set_index("sample_id")
    levels_all = sorted(meta["condition"].unique())
    reference = params["REFERENCE_LEVEL"] or levels_all[0]
    if reference not in levels_all:
        log.error("REFERENCE_LEVEL %r is not a level of %s (levels: %s)", reference, cc, ", ".join(levels_all))
        return 2
    levels = [reference] + [lv for lv in levels_all if lv != reference]
    n_per_level = meta["condition"].value_counts().reindex(levels).astype(int).to_dict()
    log.info("%d samples, %s levels %s (reference %s), label column %s, block %s, LABEL %s", len(meta), cc,
             ", ".join(f"{lv} ({n_per_level[lv]})" for lv in levels), reference, lc, bc or "none", params["LABEL"])
    notes: list[str] = []

    out = output_dir / params["LABEL"]
    out.mkdir(parents=True, exist_ok=True)

    # ---- 1. composition ------------------------------------------------------------------------
    t_comp = time.monotonic()
    lab_use = labels[labels["sample_id"].isin(meta.index)]
    label_values = sorted(lab_use[lc].unique())
    counts = (lab_use.groupby(["sample_id", lc]).size().unstack(fill_value=0)
              .reindex(index=meta.index, columns=label_values, fill_value=0))
    comp = C.composition_table(counts, meta, bc or None)
    tests = C.composition_tests(comp, levels, reference, params["FDR"])
    write_csv(comp, out / "composition.csv")
    write_csv(tests, out / "composition_tests.csv")
    if len(levels) < 2:
        notes.append(f"Only one {cc} level ({levels[0]}): composition is descriptive, nothing was tested.")
    elif not tests["tested"].any():
        notes.append(f"Composition not tested: every {cc} level needs at least 2 samples.")
    for sid in meta.index:
        statuses[sid]["n_obs_labelled"] = int(counts.loc[sid].sum())
        log.info("%s: %s=%s, %d labelled observations in %d labels", sid, cc, meta.at[sid, "condition"],
                 counts.loc[sid].sum(), int((counts.loc[sid] > 0).sum()))
    t_comp = time.monotonic() - t_comp

    # ---- 2. pseudobulk DE ----------------------------------------------------------------------
    t_pb = time.monotonic()
    min_s = params["MIN_SAMPLES_PER_GROUP"]
    with_h5ad = [s for s in meta.index if s in found.h5ad]
    for s in meta.index:
        if s not in found.h5ad:
            statuses[s]["pseudobulk"] = "no h5ad"
            log.warning("%s: no %s.h5ad, left out of the pseudobulk analysis", s, s)
    per_level_h5 = meta.loc[with_h5ad, "condition"].value_counts()
    de_possible = (len(levels) >= 2 and per_level_h5.get(reference, 0) >= min_s
                   and any(per_level_h5.get(lv, 0) >= min_s for lv in levels[1:]))
    pbs: list = []
    de_results: list = []
    pb_meta = pd.DataFrame()
    genes: list[str] = []
    if not de_possible:
        if len(levels) >= 2:
            notes.append(f"Pseudobulk DE not run: needs >= MIN_SAMPLES_PER_GROUP={min_s} samples (with an .h5ad) in "
                         f"the reference level and in at least one other level; have "
                         + ", ".join(f"{lv}: {int(per_level_h5.get(lv, 0))}" for lv in levels) + ".")
        else:
            notes.append("Pseudobulk DE not run: one condition level.")
        log.warning("%s", notes[-1])
    else:
        lab_by_sample = {s: g.set_index("obs_id")[lc] for s, g in lab_use.groupby("sample_id")}
        for s in with_h5ad:
            t = time.monotonic()
            try:
                pb = P.sample_pseudobulk(found.h5ad[s], s, lab_by_sample[s])
                pbs.append(pb)
                kept = sum(1 for v in pb.n_obs.values() if v >= params["MIN_OBS_PER_PSEUDOBULK"])
                statuses[s].update(pseudobulk="ok", counts_from=pb.source, n_obs_file=pb.n_obs_file,
                                   n_obs_matched=pb.n_obs_labelled, n_genes=len(pb.genes),
                                   n_pseudobulks=kept, seconds=round(time.monotonic() - t, 2))
                if pb.n_obs_labelled == 0:
                    log.warning("%s: none of the h5ad's observations are in labels.csv (obs_id = obs_names?)", s)
                log.info("%s: %d obs, %d with a label, %d genes, counts from %s, %d pseudobulks >= %d obs, %.2fs", s,
                         pb.n_obs_file, pb.n_obs_labelled, len(pb.genes), pb.source, kept,
                         params["MIN_OBS_PER_PSEUDOBULK"], time.monotonic() - t)
            except Exception as exc:   # one bad file is skipped, not fatal
                statuses[s].update(pseudobulk="failed", error=f"{type(exc).__name__}: {exc}")
                log.error("%s: failed: %s", s, exc)
        if not pbs:
            notes.append("Pseudobulk DE not run: every sample's .h5ad failed to read.")
            log.error("every sample's .h5ad failed")
    if pbs:
        genes = list(pbs[0].genes)
        common = set(genes)
        for pb in pbs[1:]:
            common &= set(pb.genes)
        if len(common) < len(set().union(*[set(pb.genes) for pb in pbs])):
            notes.append(f"Samples have different genes; the {len(common)} genes in every sample were used.")
        genes = [g for g in genes if g in common]
        cols, rows = {}, []
        for pb in pbs:
            idx = pd.Index(pb.genes)
            pos = idx.get_indexer(genes)
            for i, lab in enumerate(pb.labels):
                n = pb.n_obs[lab]
                used = n >= params["MIN_OBS_PER_PSEUDOBULK"]
                rows.append({"pseudobulk": f"{pb.sample_id}:{lab}", "sample_id": pb.sample_id, "label": lab,
                             "condition": meta.at[pb.sample_id, "condition"],
                             **({"block": meta.at[pb.sample_id, "block"]} if bc else {}),
                             "n_obs": n, "total_counts": int(pb.sums[i, pos].sum()), "used": used})
                if used:
                    cols[f"{pb.sample_id}:{lab}"] = pb.sums[i, pos]
        pb_meta = pd.DataFrame(rows)
        pb_counts = pd.DataFrame(cols, index=pd.Index(genes, name="gene"))
        pb_counts.round().astype(np.int64).to_csv(out / "pseudobulk_counts.csv.gz", compression="gzip")
        write_csv(pb_meta, out / "pseudobulk_samples.csv")
        pb_counts = pb_counts.T   # pseudobulks x genes
        del pbs
        for lab in label_values:
            sel = pb_meta[(pb_meta["label"] == lab) & pb_meta["used"]]
            t = time.monotonic()
            if sel.empty:
                r = P.LabelResult(lab, "skipped", reason=f"no pseudobulk with >= {params['MIN_OBS_PER_PSEUDOBULK']} obs",
                                  n_per_level={lv: 0 for lv in levels})
            else:
                cnt = pb_counts.loc[sel["pseudobulk"]].set_axis(sel["sample_id"].to_numpy())
                md = meta.loc[sel["sample_id"]]
                r = P.de_for_label(lab, cnt, md, reference, levels, bool(bc), params, n_cpus)
            de_results.append(r)
            n_sig = sum(int((x["p_adj"] < params["FDR"]).sum()) for x in r.results)
            log.info("label %s: %s%s, %s, %d genes, %d significant, %.1fs", lab, r.status,
                     f" ({r.reason})" if r.reason else "",
                     ", ".join(f"{lv} {n}" for lv, n in r.n_per_level.items()), r.n_genes_tested, n_sig,
                     time.monotonic() - t)
    t_pb = time.monotonic() - t_pb

    # ---- DE tables -----------------------------------------------------------------------------
    de_cols = ["label", "comparison", "gene", "base_mean", "log2_fold_change", "lfc_se", "stat", "p_val", "p_adj",
               "significant", "n_samples_level", "n_samples_reference"]
    frames = [x for r in de_results for x in r.results]
    de = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=de_cols)
    de["significant"] = de["p_adj"].astype(float) < params["FDR"] if len(de) else []
    de = de[de_cols]
    if len(de):
        de = de.sort_values(["label", "comparison", "p_val"], na_position="last")
    write_csv(de, out / "de.csv")
    summ_rows = []
    for r in de_results:
        comps = [f"{lv} vs {reference}" for lv in (r.levels[1:] if r.levels else levels[1:])]
        for comp_name in comps:
            x = next((f for f in r.results if f["comparison"].iloc[0] == comp_name), None) if r.results else None
            row = {"label": r.label, "comparison": comp_name, "status": r.status if x is not None or r.status != "tested"
                   else "skipped", "reason": r.reason, "design": r.design, "note": r.block_note,
                   **{f"n_samples_{lv}": r.n_per_level.get(lv, 0) for lv in levels},
                   "n_genes_tested": r.n_genes_tested}
            if x is not None:
                sig = x[x["p_adj"] < params["FDR"]]
                row.update(n_significant=len(sig), n_up=int((sig["log2_fold_change"] > 0).sum()),
                           n_down=int((sig["log2_fold_change"] < 0).sum()),
                           top_up=";".join(sig[sig["log2_fold_change"] > 0].nsmallest(5, "p_val")["gene"]),
                           top_down=";".join(sig[sig["log2_fold_change"] < 0].nsmallest(5, "p_val")["gene"]))
            else:
                row.update(n_significant=None, n_up=None, n_down=None, top_up="", top_down="")
            summ_rows.append(row)
    de_summary = pd.DataFrame(summ_rows, columns=["label", "comparison", "status", "reason", "design", "note"]
                              + [f"n_samples_{lv}" for lv in levels]
                              + ["n_genes_tested", "n_significant", "n_up", "n_down", "top_up", "top_down"])
    write_csv(de_summary, out / "de_summary.csv")

    # ---- figures -------------------------------------------------------------------------------
    t_fig = time.monotonic()
    figures = []
    try:
        figures.append(F.composition_boxplots(comp, tests, levels, f"{lc} composition by {cc} · {params['LABEL']}",
                                              out / "composition.png"))
        tested = de_summary[de_summary["status"] == "tested"].sort_values("n_significant", ascending=False)
        multi = len(levels) > 2
        for _, row in tested.head(params["MAX_PLOTS"]).iterrows():
            sub = de[(de["label"] == row["label"]) & (de["comparison"] == row["comparison"])]
            level = row["comparison"].split(" vs ")[0]
            name = f"volcano_{F.safe_name(row['label'])}" + (f"_{F.safe_name(level)}" if multi else "") + ".png"
            figures.append(F.volcano(sub, params["FDR"], f"{lc} {row['label']}: {row['comparison']} "
                                     f"({row['n_significant']} at FDR {params['FDR']:g})", out / name))
        if len(de):
            figures.append(F.de_heatmap(de, params["FDR"], f"Top DE genes per {lc} · {params['LABEL']}",
                                        out / "de_heatmap.png"))
    except Exception as exc:   # figures are a convenience; the tables stand without them
        log.exception("figures failed: %s", exc)
    t_fig = time.monotonic() - t_fig

    # ---- summary / provenance ------------------------------------------------------------------
    sig_comp = tests[tests["significant"]] if len(tests) else tests
    tested_de = de_summary[de_summary["status"] == "tested"]
    summary = {
        "label": params["LABEL"], "label_column": lc, "condition_column": cc, "block_column": bc or None,
        "reference": reference, "levels": levels, "n_samples": len(meta), "n_samples_per_level": n_per_level,
        "labels": label_values,
        "composition": {
            "tests": ("welch_t + mann_whitney_u" if len(levels) == 2 else "anova + kruskal_wallis")
            if len(levels) >= 2 else None,
            "n_labels": len(label_values), "n_tested": int(tests["tested"].sum()) if len(tests) else 0,
            "n_significant": int(len(sig_comp)),
            "significant": [{"label": r["label"], "comparison": r["comparison"],
                             "clr_difference": r["clr_difference"], "p_adj": r["parametric_p_adj"],
                             "rank_p_adj": r["rank_p_adj"],
                             **{f"mean_fraction_{lv}": r[f"mean_fraction_{lv}"] for lv in levels}}
                            for _, r in sig_comp.iterrows()],
        },
        "de": {
            "run": bool(de_possible and len(de_results) > 0),
            "n_labels_tested": int(tested_de["label"].nunique()),
            "n_labels_skipped": int(de_summary.loc[de_summary["status"] != "tested", "label"].nunique()),
            "n_genes": len(genes),
            "n_significant_total": int(tested_de["n_significant"].fillna(0).sum()) if len(tested_de) else 0,
            "per_label": [{k: r[k] for k in ("label", "comparison", "status", "reason", "n_genes_tested",
                                             "n_significant", "n_up", "n_down", "top_up", "top_down")}
                          for _, r in de_summary.iterrows()],
        },
        "fdr": params["FDR"], "notes": notes, "figures": [p.name for p in figures],
    }
    write_json(out / "summary.json", summary)
    n_failed = sum(1 for s in statuses.values() if s.get("pseudobulk") == "failed")
    write_json(out / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(),
        "parameters": params, "resolved": {"reference": reference, "levels": levels,
                                           "threads": {"n": n_cpus, "source": os.environ.get("ST_COMPARE_THREADS_SOURCE", "THREADS")}},
        "method": {
            "composition": "per sample x label counts from labels.csv; CLR with pseudocount 0.5; per label Welch t + "
                           "Mann-Whitney U (2 levels) or one-way ANOVA + Kruskal-Wallis (> 2 levels) on CLR; "
                           "Benjamini-Hochberg across labels per test",
            "pseudobulk": "sum of layers['counts'] (else integer X) per sample x label, >= MIN_OBS_PER_PSEUDOBULK obs",
            "gene_filter": ">= MIN_COUNT_GENE counts in >= (smallest group size) pseudobulks of the label",
            "de": "PyDESeq2 DeseqDataSet(design '~ condition' or '~ block + condition', refit_cooks) + DeseqStats Wald "
                  "test per level vs reference, Cook's and independent filtering, BH per label x comparison",
        },
        "inputs": {"h5ad": [_rel(p, input_dir) for p in found.h5ad.values()],
                   "labels_csv": _rel(found.labels_csv, input_dir), "samples_csv": _rel(found.samples_csv, input_dir),
                   "skipped": found.skipped},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"composition": round(t_comp, 2), "pseudobulk_de": round(t_pb, 2),
                            "figures": round(t_fig, 2), "total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": peak_rss_mb(), "cgroup_memory_peak_mb": cgroup_peak_mb(),
        "samples": list(statuses.values()), "notes": notes,
        "outputs": sorted(p.name for p in out.iterdir()) + ["provenance.json"],
    })
    log.info("done: %d samples, %d labels, composition %d/%d significant, DE %d labels tested, %d genes significant, "
             "%.1fs, peak RSS %s MB -> %s", len(meta), len(label_values), summary["composition"]["n_significant"],
             summary["composition"]["n_tested"], summary["de"]["n_labels_tested"],
             summary["de"]["n_significant_total"], time.monotonic() - t0, peak_rss_mb(), out)
    for n in notes:
        log.info("note: %s", n)
    if de_possible and with_h5ad and n_failed == len(with_h5ad):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
