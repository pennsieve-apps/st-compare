"""Figures: composition boxplots, volcano plots, DE heatmap (PNG, 150 dpi)."""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

DPI = 150
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
UP, DOWN, NS = "#e34948", "#2a78d6", "#b9b8b2"
TEXT, TEXT2, GRID, MID = "#0b0b0b", "#52514e", "#e4e3df", "#f0efec"


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "x"


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                         "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white"})
    return plt


def composition_boxplots(comp: pd.DataFrame, tests: pd.DataFrame, levels: list[str], title: str,
                         path: Path) -> Path:
    """One panel per label: fraction of each sample's observations, a box per level with the samples as points."""
    plt = _plt()
    labels = list(tests["label"]) if not tests.empty else sorted(comp["label"].unique())
    n = len(labels)
    ncol = min(4, n)
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol + 0.6, 2.5 * nrow + 0.9), squeeze=False)
    blocks = sorted(comp["block"].unique()) if "block" in comp.columns else []
    rng = np.random.default_rng(0)
    for k, lab in enumerate(labels):
        ax = axes[k // ncol][k % ncol]
        sub = comp[comp["label"] == lab]
        data = [sub.loc[sub["condition"] == lv, "fraction"].to_numpy() * 100 for lv in levels]
        bp = ax.boxplot(data, positions=range(len(levels)), widths=0.55, showfliers=False, patch_artist=True,
                        medianprops={"color": TEXT, "linewidth": 1.5}, whiskerprops={"color": TEXT2},
                        capprops={"color": TEXT2})
        for i, box in enumerate(bp["boxes"]):
            box.set(facecolor="none", edgecolor=TEXT2)
        for i, lv in enumerate(levels):
            s = sub[sub["condition"] == lv]
            x = i + rng.uniform(-0.12, 0.12, len(s))
            if blocks:
                for b, (xx, yy) in zip(s["block"], zip(x, s["fraction"] * 100)):
                    ax.scatter(xx, yy, s=28, color=SERIES[i % len(SERIES)], marker=MARKERS[blocks.index(b) % 8],
                               edgecolor="white", linewidth=0.8, zorder=3)
            else:
                ax.scatter(x, s["fraction"] * 100, s=28, color=SERIES[i % len(SERIES)], edgecolor="white",
                           linewidth=0.8, zorder=3)
        ax.set_xticks(range(len(levels)), levels, rotation=30 if len(levels) > 2 else 0)
        row = tests[tests["label"] == lab]
        sub_t = ""
        if not row.empty and bool(row["tested"].iloc[0]):
            sub_t = f"\nq={row['parametric_p_adj'].iloc[0]:.2g} (t/F), q={row['rank_p_adj'].iloc[0]:.2g} (rank)"
        ax.set_title(f"{lab}{sub_t}", fontsize=9, color=TEXT)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        if k % ncol == 0:
            ax.set_ylabel("% of the sample's observations")
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    if blocks:
        from matplotlib.lines import Line2D

        handles = [Line2D([], [], marker=MARKERS[i % 8], linestyle="", color=TEXT2, label=str(b))
                   for i, b in enumerate(blocks)]
        fig.legend(handles=handles, title="block", loc="upper right", frameon=False, fontsize=8)
    fig.suptitle(title, fontsize=10, color=TEXT)
    fig.tight_layout(rect=(0, 0, 0.93 if blocks else 1, 0.97))
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def volcano(res: pd.DataFrame, fdr: float, title: str, path: Path, n_labels: int = 10) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    r = res.dropna(subset=["p_val", "log2_fold_change"]).copy()
    r["nlp"] = -np.log10(np.clip(r["p_val"].to_numpy(dtype=float), 1e-300, 1))
    sig = r["p_adj"].fillna(1) < fdr
    up, down = sig & (r["log2_fold_change"] > 0), sig & (r["log2_fold_change"] < 0)
    ax.scatter(r.loc[~sig, "log2_fold_change"], r.loc[~sig, "nlp"], s=6, color=NS, linewidth=0, label="not significant")
    ax.scatter(r.loc[up, "log2_fold_change"], r.loc[up, "nlp"], s=12, color=UP, linewidth=0, label=f"up ({int(up.sum())})")
    ax.scatter(r.loc[down, "log2_fold_change"], r.loc[down, "nlp"], s=12, color=DOWN, linewidth=0,
               label=f"down ({int(down.sum())})")
    for _, g in r[sig].nsmallest(n_labels, "p_val").iterrows():
        ax.annotate(str(g["gene"]), (g["log2_fold_change"], g["nlp"]), fontsize=7, color=TEXT,
                    xytext=(3, 2), textcoords="offset points")
    ax.axvline(0, color=GRID, linewidth=0.8, zorder=0)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p")
    ax.set_title(title, fontsize=9, color=TEXT)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def de_heatmap(de: pd.DataFrame, fdr: float, title: str, path: Path, per_test: int = 5, max_genes: int = 50) -> Path:
    """Top genes (by adjusted p) of each label x comparison, coloured by log2FC in every label x comparison."""
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap

    de = de.copy()
    de["test"] = de["label"].astype(str) + (": " + de["comparison"] if de["comparison"].nunique() > 1 else "")
    sig = de[de["p_adj"].fillna(1) < fdr]
    genes: list[str] = []
    for _, g in sig.sort_values("p_adj").groupby("test", sort=False):
        for gene in g["gene"].head(per_test):
            if gene not in genes:
                genes.append(gene)
    genes = genes[:max_genes]
    tests = list(dict.fromkeys(de["test"]))
    fig_h = 1.4 + 0.18 * max(len(genes), 4)
    fig, ax = plt.subplots(figsize=(1.6 + 0.55 * len(tests), fig_h))
    if not genes:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No gene reached FDR < {fdr:g} in any label", ha="center", va="center", color=TEXT2)
    else:
        mat = de.pivot_table(index="gene", columns="test", values="log2_fold_change", aggfunc="first")
        mat = mat.reindex(index=genes, columns=tests)
        q = de.pivot_table(index="gene", columns="test", values="p_adj", aggfunc="first").reindex(index=genes, columns=tests)
        vmax = float(np.nanpercentile(np.abs(mat.to_numpy()), 98)) if np.isfinite(mat.to_numpy()).any() else 1.0
        vmax = max(vmax, 0.5)
        cmap = LinearSegmentedColormap.from_list("div", [DOWN, MID, UP]).with_extremes(bad="#ffffff")
        im = ax.imshow(np.ma.masked_invalid(mat.to_numpy()), cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
        ys, xs = np.nonzero(q.fillna(1).to_numpy() < fdr)
        ax.scatter(xs, ys, s=6, color=TEXT, marker="o", linewidth=0)
        ax.set_xticks(range(len(tests)), tests, rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(len(genes)), genes, fontsize=7)
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label("log2 fold change", fontsize=8)
        ax.spines[:].set_visible(False)
    ax.set_title(title + ("\n(dot: FDR < %g; blank: not tested)" % fdr if genes else ""), fontsize=9, color=TEXT)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path
