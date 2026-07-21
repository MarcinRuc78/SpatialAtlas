#!/usr/bin/env python3
"""Check H5AD, tissue-image, and marker files used by SpatialAtlas."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--image", required=True, type=Path)
    parser.add_argument("--markers", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"name": name, "status": "PASS" if condition else "FAIL", "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    dataset = ad.read_h5ad(args.h5ad)
    check("non_empty", dataset.n_obs > 0 and dataset.n_vars > 0,
          f"shape={dataset.n_obs}x{dataset.n_vars}")
    check("unique_barcodes", dataset.obs_names.is_unique, "obs_names are unique")
    check("unique_genes", dataset.var_names.is_unique, "var_names are unique")
    check("cluster", "cluster" in dataset.obs, "obs['cluster'] is present")
    check("spatial", "spatial" in dataset.obsm, "obsm['spatial'] is present")
    spatial = np.asarray(dataset.obsm["spatial"])
    check("spatial_shape", spatial.shape == (dataset.n_obs, 2), f"shape={spatial.shape}")
    check("spatial_finite", bool(np.isfinite(spatial).all()), "all coordinates are finite")
    has_umap = "UMAP" in dataset.obsm
    if has_umap:
        umap = np.asarray(dataset.obsm["UMAP"])
        check("umap_shape", umap.shape == (dataset.n_obs, 2), f"shape={umap.shape}")
        check("umap_finite", bool(np.isfinite(umap).all()), "all UMAP values are finite")

    with Image.open(args.image) as image:
        image.verify()
    with Image.open(args.image) as image:
        image_size = image.size
    check("image", image_size[0] > 0 and image_size[1] > 0, f"size={image_size}")

    marker_overlap = None
    if args.markers:
        markers = pd.read_csv(args.markers)
        check("markers_columns", {"gene", "cluster"}.issubset(markers.columns),
              f"columns={list(markers.columns)}")
        marker_overlap = int(markers["gene"].isin(dataset.var_names).sum())
        check("markers_overlap", marker_overlap > 0,
              f"overlapping rows={marker_overlap}/{len(markers)}")

    report = {
        "status": "PASS",
        "artifact": str(args.h5ad),
        "sha256": sha256(args.h5ad),
        "image_sha256": sha256(args.image),
        "n_obs": dataset.n_obs,
        "n_vars": dataset.n_vars,
        "clusters": dataset.obs["cluster"].astype(str).value_counts().sort_index().to_dict(),
        "matrix_type": type(dataset.X).__name__,
        "has_umap": has_umap,
        "image_size": list(image_size),
        "marker_overlap_rows": marker_overlap,
        "python": sys.version,
        "platform": platform.platform(),
        "anndata": ad.__version__,
        "checks": checks,
    }
    (args.output_dir / f"{args.prefix}_artifact_validation.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
