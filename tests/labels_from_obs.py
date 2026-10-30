"""
Test harness: write a labels.csv (sample_id, obs_id, <column>) from an obs column of
<sample_id>.h5ad files, e.g. the pathologist annotations st-prep joins onto obs["annotation"],
so st-compare can be run without st-domains / st-celltype.

    python -m tests.labels_from_obs <folder with .h5ad> <out>/labels.csv --column annotation
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spatial", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--column", default="annotation")
    args = ap.parse_args(argv)
    frames = []
    for p in sorted(args.spatial.glob("*.h5ad")):
        with h5py.File(p, "r") as f:
            obs = f["obs"]
            if args.column not in obs:
                print(f"{p.name}: no obs[{args.column!r}], skipped")
                continue
            idx = read_elem(obs[obs.attrs["_index"]])
            vals = pd.Series(read_elem(obs[args.column])).astype(str)
        frames.append(pd.DataFrame({"sample_id": p.stem, "obs_id": idx, args.column: vals.to_numpy()}))
        print(f"{p.name}: {len(idx)} observations, {vals.value_counts().to_dict()}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(args.out, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
