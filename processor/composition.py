"""
Label composition per sample and per-label tests between condition levels.

Each sample contributes one value per label: its fraction of the sample's
observations, and the centred log-ratio (CLR) of the counts with a 0.5
pseudocount. CLR removes the unit-sum constraint of proportions (a label
growing makes every other label's proportion shrink); the tests compare the
CLR values of the samples between condition levels:

  two levels:   Welch's t-test and the Mann-Whitney U test (two-sided)
  more levels:  one-way ANOVA and the Kruskal-Wallis H test

Benjamini-Hochberg is applied across labels, separately for each test. The
sample, not the spot or cell, is the unit: n is the number of samples.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

PSEUDOCOUNT = 0.5


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def composition_table(counts: pd.DataFrame, meta: pd.DataFrame, block: str | None) -> pd.DataFrame:
    """counts: samples x labels (n_obs). meta: index sample_id with `condition` (and `block`).
    Long table: sample_id, condition, [block], label, n_obs, n_total, fraction, clr."""
    counts = counts.loc[meta.index]
    totals = counts.sum(axis=1)
    frac = counts.div(totals.where(totals > 0), axis=0)
    logc = np.log(counts.astype(float) + PSEUDOCOUNT)
    clr = logc.sub(logc.mean(axis=1), axis=0)
    rows = []
    for sid in counts.index:
        for lab in counts.columns:
            r = {"sample_id": sid, "condition": meta.at[sid, "condition"]}
            if block:
                r["block"] = meta.at[sid, "block"]
            r.update(label=lab, n_obs=int(counts.at[sid, lab]), n_total=int(totals[sid]),
                     fraction=float(frac.at[sid, lab]), clr=float(clr.at[sid, lab]))
            rows.append(r)
    return pd.DataFrame(rows)


def composition_tests(comp: pd.DataFrame, levels: list[str], reference: str, fdr: float) -> pd.DataFrame:
    """One row per label. `levels` are the condition levels with samples, reference first."""
    rows = []
    two = len(levels) == 2
    for lab, sub in comp.groupby("label", sort=False):
        r: dict = {"label": lab, "n_samples": int(len(sub))}
        groups = [sub.loc[sub["condition"] == lv, "clr"].to_numpy() for lv in levels]
        for lv, g in zip(levels, groups):
            r[f"n_{lv}"] = int(len(g))
            r[f"mean_fraction_{lv}"] = float(sub.loc[sub["condition"] == lv, "fraction"].mean()) if len(g) else np.nan
            r[f"mean_clr_{lv}"] = float(g.mean()) if len(g) else np.nan
        testable = len(levels) >= 2 and all(len(g) >= 2 for g in groups)
        if two:
            r["comparison"] = f"{levels[1]} vs {reference}"
            r["clr_difference"] = r[f"mean_clr_{levels[1]}"] - r[f"mean_clr_{reference}"]
            r["parametric_test"], r["rank_test"] = "welch_t", "mann_whitney_u"
        else:
            r["comparison"] = "any difference between " + ", ".join(levels)
            r["clr_difference"] = np.nan
            r["parametric_test"], r["rank_test"] = "anova", "kruskal_wallis"
        r.update(parametric_stat=np.nan, parametric_p=np.nan, rank_stat=np.nan, rank_p=np.nan, tested=False)
        constant = np.allclose(np.concatenate(groups), np.concatenate(groups)[0]) if testable else True
        if testable and not constant:
            with np.errstate(all="ignore"):
                if two:
                    t = stats.ttest_ind(groups[1], groups[0], equal_var=False)
                    u = stats.mannwhitneyu(groups[1], groups[0], alternative="two-sided")
                    r.update(parametric_stat=float(t.statistic), parametric_p=float(t.pvalue),
                             rank_stat=float(u.statistic), rank_p=float(u.pvalue))
                else:
                    f = stats.f_oneway(*groups)
                    h = stats.kruskal(*groups)
                    r.update(parametric_stat=float(f.statistic), parametric_p=float(f.pvalue),
                             rank_stat=float(h.statistic), rank_p=float(h.pvalue))
            r["tested"] = True
        rows.append(r)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["parametric_p_adj"] = bh(out["parametric_p"].to_numpy())
    out["rank_p_adj"] = bh(out["rank_p"].to_numpy())
    out["significant"] = out["parametric_p_adj"] < fdr
    lead = ["label", "comparison", "n_samples"] + [f"n_{lv}" for lv in levels]
    means = [f"mean_fraction_{lv}" for lv in levels] + [f"mean_clr_{lv}" for lv in levels]
    tail = ["clr_difference", "parametric_test", "parametric_stat", "parametric_p", "parametric_p_adj",
            "rank_test", "rank_stat", "rank_p", "rank_p_adj", "tested", "significant"]
    return out[lead + means + tail].sort_values(["parametric_p", "label"], na_position="last").reset_index(drop=True)
