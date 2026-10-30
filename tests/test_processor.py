"""End-to-end tests on the synthetic fixture (tests/make_fixture.py), no Docker."""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from tests.make_fixture import DE_DOWN, DE_LABEL, DE_UP, make_fixture

PARAMS = ("LABEL", "LABEL_COLUMN", "CONDITION_COLUMN", "REFERENCE_LEVEL", "BLOCK_COLUMN", "MIN_OBS_PER_PSEUDOBULK",
          "MIN_SAMPLES_PER_GROUP", "MIN_COUNT_GENE", "FDR", "MAX_PLOTS", "SEED", "THREADS")


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("fixture")
    make_fixture(root, "both")
    return root


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("THREADS", "1")
    monkeypatch.setenv("ST_COMPARE_THREADS", "1")
    for k, v in params.items():
        monkeypatch.setenv(k, str(v))
    from processor.main import main

    return main()


def copy_with_samples(src: Path, dst: Path, edit) -> Path:
    shutil.copytree(src, dst)
    p = next(dst.rglob("samples.csv"))
    df = pd.read_csv(p, dtype=str)
    edit(df).to_csv(p, index=False)
    return dst


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_effects(fixture_root, tmp_path, monkeypatch, layout):
    assert run(fixture_root / layout, tmp_path, monkeypatch) == 0
    out = tmp_path / "domain"
    for name in ("composition.csv", "composition_tests.csv", "de.csv", "de_summary.csv", "pseudobulk_counts.csv.gz",
                 "pseudobulk_samples.csv", "summary.json", "provenance.json", "composition.png", "de_heatmap.png"):
        assert (out / name).is_file(), name
    comp = pd.read_csv(out / "composition.csv")
    assert set(comp.columns) >= {"sample_id", "condition", "label", "n_obs", "n_total", "fraction", "clr"}
    assert np.allclose(comp.groupby("sample_id")["fraction"].sum(), 1.0)
    assert np.allclose(comp.groupby("sample_id")["clr"].sum(), 0.0, atol=1e-5)
    tests = pd.read_csv(out / "composition_tests.csv").set_index("label")
    # Planted: D03 20 % -> 35 %, D01 50 % -> 35 %; D02 unchanged.
    assert tests.at["D03", "clr_difference"] > 0 and tests.at["D01", "clr_difference"] < 0
    assert tests.at["D03", "parametric_p_adj"] < 0.05 and tests.at["D01", "parametric_p_adj"] < 0.05
    assert tests.at["D02", "parametric_p_adj"] > 0.05
    assert tests.at["D03", "parametric_test"] == "welch_t" and tests.at["D03", "rank_test"] == "mann_whitney_u"

    de = pd.read_csv(out / "de.csv")
    sig = de[de["significant"]]
    hits = set(sig.loc[sig["label"] == DE_LABEL, "gene"])
    assert len(hits & set(DE_UP + DE_DOWN)) >= 6
    assert (sig.loc[sig["gene"].isin(DE_UP) & (sig["label"] == DE_LABEL), "log2_fold_change"] > 0).all()
    assert (sig.loc[sig["gene"].isin(DE_DOWN) & (sig["label"] == DE_LABEL), "log2_fold_change"] < 0).all()
    assert len(sig[sig["label"] != DE_LABEL]) <= 2          # few false positives elsewhere
    assert set(de["comparison"]) == {"b vs a"}

    pbc = pd.read_csv(out / "pseudobulk_counts.csv.gz", index_col=0)
    assert pbc.shape == (200, 18)
    # The pseudobulk equals the sum of the counts layer over the label's observations.
    labels = pd.read_csv(next((fixture_root / layout).rglob("labels.csv")))
    with h5py.File(next((fixture_root / layout).rglob("a1.h5ad")), "r") as f:
        g = f["layers"]["counts"]
        import scipy.sparse as sp

        m = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
    rows = (labels["sample_id"] == "a1") & (labels["domain"] == "D02")
    idx = np.flatnonzero(rows.to_numpy()[labels["sample_id"].to_numpy() == "a1"])
    assert np.array_equal(np.asarray(m[idx].sum(axis=0)).ravel(), pbc["a1:D02"].to_numpy())

    summary = json.loads((out / "summary.json").read_text())
    assert summary["reference"] == "a" and summary["levels"] == ["a", "b"]
    assert summary["de"]["n_labels_tested"] == 3
    assert summary["composition"]["n_significant"] == 2
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["app"] == "st-compare" and prov["version"] == "0.1.0" and prov["tools"]["pydeseq2"]
    assert len(list(out.glob("volcano_*.png"))) == 3


def test_block_and_reference_level(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root / "linked", tmp_path, monkeypatch, BLOCK_COLUMN="donor", REFERENCE_LEVEL="b",
               LABEL="by-donor", MAX_PLOTS=1) == 0
    out = tmp_path / "by-donor"
    s = pd.read_csv(out / "de_summary.csv")
    assert (s["design"] == "~ block + condition").all()
    assert set(s["comparison"]) == {"a vs b"}
    de = pd.read_csv(out / "de.csv")
    up = de[(de["label"] == DE_LABEL) & de["gene"].isin(DE_UP) & de["significant"]]
    assert len(up) >= 3 and (up["log2_fold_change"] < 0).all()      # reversed: a vs b
    assert len(list(out.glob("volcano_*.png"))) == 1
    assert "block" in pd.read_csv(out / "composition.csv").columns


def test_confounded_block_is_dropped(fixture_root, tmp_path, monkeypatch):
    def edit(df):
        df["batch"] = df["condition"]
        return df
    inp = copy_with_samples(fixture_root / "linked", tmp_path / "in", edit)
    assert run(inp, tmp_path / "out", monkeypatch, BLOCK_COLUMN="batch") == 0
    s = pd.read_csv(tmp_path / "out" / "domain" / "de_summary.csv")
    assert (s["design"] == "~ condition").all()
    assert s["note"].str.contains("confounded").all()


def test_three_levels(fixture_root, tmp_path, monkeypatch):
    def edit(df):
        df.loc[df["sample_id"].isin(["b2", "b3"]), "condition"] = "c"
        df.loc[df["sample_id"] == "a3", "condition"] = "b"
        return df
    # a: a1, a2; b: a3, b1; c: b2, b3
    inp = copy_with_samples(fixture_root / "linked", tmp_path / "in", edit)
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    out = tmp_path / "out" / "domain"
    tests = pd.read_csv(out / "composition_tests.csv")
    assert set(tests["parametric_test"]) == {"anova"} and set(tests["rank_test"]) == {"kruskal_wallis"}
    de = pd.read_csv(out / "de.csv")
    assert set(de["comparison"]) == {"b vs a", "c vs a"}
    assert (out / f"volcano_{DE_LABEL}_c.png").is_file()


def test_one_condition_composition_only(fixture_root, tmp_path, monkeypatch):
    def edit(df):
        df["condition"] = "all"
        return df
    inp = copy_with_samples(fixture_root / "linked", tmp_path / "in", edit)
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    out = tmp_path / "out" / "domain"
    summary = json.loads((out / "summary.json").read_text())
    assert any("one" in n.lower() for n in summary["notes"])
    assert not summary["de"]["run"]
    assert len(pd.read_csv(out / "composition.csv")) == 18
    assert not pd.read_csv(out / "composition_tests.csv")["tested"].any()
    assert len(pd.read_csv(out / "de.csv")) == 0
    assert (out / "composition.png").is_file()


def test_too_few_samples(fixture_root, tmp_path, monkeypatch):
    def edit(df):
        return df[df["sample_id"].isin(["a1", "a2", "b1"])]
    inp = copy_with_samples(fixture_root / "linked", tmp_path / "in", edit)
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    summary = json.loads((tmp_path / "out" / "domain" / "summary.json").read_text())
    assert not summary["de"]["run"] and any("MIN_SAMPLES_PER_GROUP" in n for n in summary["notes"])
    assert summary["composition"]["n_tested"] == 0


def test_min_obs_skips_labels(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root / "linked", tmp_path, monkeypatch, MIN_OBS_PER_PSEUDOBULK=95) == 0
    s = pd.read_csv(tmp_path / "domain" / "de_summary.csv").set_index("label")
    # D03 holds ~60 obs in `a` samples: too few pseudobulks in the reference.
    assert s.at["D03", "status"] == "skipped" and "reference" in s.at["D03", "reason"]
    assert s.at["D01", "status"] == "tested"


def test_bad_parameters(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root / "linked", tmp_path, monkeypatch, LABEL_COLUMN="nope") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, CONDITION_COLUMN="nope") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, REFERENCE_LEVEL="zzz") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, LABEL="Bad_Label") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, FDR="2") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, MIN_SAMPLES_PER_GROUP="1") == 2


def test_missing_inputs(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_log_normalized_x_without_counts_fails(fixture_root, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    shutil.copytree(fixture_root / "linked", inp)
    for p in (inp / "spatial").glob("*.h5ad"):
        with h5py.File(p, "a") as f:
            del f["layers"]["counts"]
    assert run(inp, tmp_path / "out", monkeypatch) == 1
    prov = json.loads((tmp_path / "out" / "domain" / "provenance.json").read_text())
    assert all(s.get("pseudobulk") == "failed" for s in prov["samples"])
    assert (tmp_path / "out" / "domain" / "composition_tests.csv").is_file()


def test_handler(fixture_root, tmp_path, monkeypatch):
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("ST_COMPARE_THREADS", "1")
    monkeypatch.setenv("FDR", "0.5")          # must be cleared by the handler
    from processor.handler import handler

    r = handler({"inputDir": str(fixture_root / "linked"), "outputDir": str(tmp_path), "LABEL": "lam",
                 "MAX_PLOTS": 0}, None)
    assert r["status"] == "ok"
    summary = json.loads((tmp_path / "lam" / "summary.json").read_text())
    assert summary["fdr"] == 0.05
    assert not list((tmp_path / "lam").glob("volcano_*.png"))
    with pytest.raises(RuntimeError):
        handler({"inputDir": str(fixture_root / "linked"), "outputDir": str(tmp_path), "LABEL_COLUMN": "x"}, None)
