"""
Finding the inputs: <sample_id>.h5ad files (port `spatial`, from st-qc),
labels.csv (port `labels`, from st-domains, st-celltype or cell2location) and
samples.csv (port `metadata`, from st-prep). Linked layout: INPUT_DIR/spatial/,
INPUT_DIR/labels/, INPUT_DIR/metadata/; without links the files are found by
name anywhere under INPUT_DIR (symlinks followed).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("st-compare")


@dataclass
class Inputs:
    h5ad: dict[str, Path]                 # sample_id (file stem) -> path
    labels_csv: Path | None
    samples_csv: Path | None
    files_seen: list[Path]
    skipped: dict[str, str] = field(default_factory=dict)


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinks (INPUT_DIR entries are symlinks), sorted."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return found


def _find(root: Path, port: str, predicate) -> list[Path]:
    """Files of a port: its own folder when it exists, else anywhere under INPUT_DIR."""
    folder = root / port
    base = folder if folder.is_dir() else root
    return [p for p in walk(base) if predicate(p)]


def _one(paths: list[Path], what: str) -> Path | None:
    if len(paths) > 1:
        log.warning("%d %s files; using %s", len(paths), what, paths[0])
    return paths[0] if paths else None


def discover(input_dir: Path) -> Inputs:
    h5ads = _find(input_dir, "spatial", lambda p: p.suffix.lower() == ".h5ad")
    labels = _find(input_dir, "labels", lambda p: p.name.lower() == "labels.csv")
    samples = _find(input_dir, "metadata", lambda p: p.name.lower() == "samples.csv")
    skipped: dict[str, str] = {}
    by_stem: dict[str, Path] = {}
    for p in h5ads:
        if p.stem in by_stem:
            skipped[str(p)] = f"two files for sample {p.stem}; using {by_stem[p.stem]}"
            continue
        by_stem[p.stem] = p
    return Inputs(by_stem, _one(labels, "labels.csv"), _one(samples, "samples.csv"),
                  sorted(h5ads + labels + samples), skipped)
