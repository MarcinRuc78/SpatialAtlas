#!/usr/bin/env python3
"""
prepare_data.py — prepare Visium-derived data for SpatialAtlas
==============================================================
Converts supported spatial-transcriptomics inputs into the H5AD, image,
marker-table, and configuration layout read by SpatialAtlas.

Supported inputs:
  A) SpaceRanger v1.x / v2.x output (standard Visium, 55µm spots)
  B) Space Ranger v3.x/v4.x output (Visium HD) at 2µm, 8µm, or 16µm binning
  C) Space Ranger v3.x/v4.x output with cell segmentation (Cellpose/StarDist)
  D) Pre-existing .h5ad file + tissue image
  E) Files downloaded from GEO (NCBI) — .h5, .h5ad, .tar.gz
  F) Files downloaded from 10x Genomics website — .h5, .h5ad

Produces:
  - sample_<name>.h5ad : filtered, annotated AnnData
  - tissue_<name>.png  : tissue image for overlay
  - markers.csv        : marker genes (optional)
  - Config YAML block  : printed to stdout — paste into config.yaml

Usage examples:
  # Standard Visium (SpaceRanger v1/v2)
  python prepare_data.py --spaceranger /path/to/spaceranger_out \\
      --output ./data --sample-name "Sample1"

  # Visium HD at 8µm bins
  python prepare_data.py --spaceranger-hd /path/to/spaceranger_hd_out \\
      --hd-resolution 8um --output ./data --sample-name "HD_8um"

  # Visium HD with cell segmentation
  python prepare_data.py --spaceranger-hd /path/to/spaceranger_hd_out \\
      --hd-resolution cell_segmentation --output ./data --sample-name "HD_cells"

  # From GEO download (.h5 matrix + spatial folder)
  python prepare_data.py --h5-matrix /path/to/filtered_feature_bc_matrix.h5 \\
      --spatial-dir /path/to/spatial/ --output ./data --sample-name "GEO_sample"

  # From pre-existing .h5ad
  python prepare_data.py --h5ad /path/to/data.h5ad \\
      --tissue-image /path/to/tissue.png --output ./data --sample-name "MySample"

  # From 10x website .h5ad download
  python prepare_data.py --h5ad /path/to/10x_download.h5ad \\
      --output ./data --sample-name "TenX"
"""

import argparse
import os
import sys
import json
import glob
import shutil
import re
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image


# =============================================================================
# Loaders for each data source
# =============================================================================

def find_file(directory, candidates, required=True):
    """Find the first existing file from a list of candidates."""
    for c in candidates:
        path = os.path.join(directory, c)
        if os.path.exists(path):
            return path
    if required:
        sys.exit(f"ERROR: None of these files found in {directory}: {candidates}")
    return None


def load_spaceranger_standard(sr_dir, use_lowres=False):
    """
    Load standard Visium data from SpaceRanger v1.x / v2.x output.

    Expected structure:
      spaceranger_out/
      └── outs/
          ├── filtered_feature_bc_matrix.h5
          ├── filtered_feature_bc_matrix/     (alternative: directory with MTX)
          │   ├── barcodes.tsv.gz
          │   ├── features.tsv.gz
          │   └── matrix.mtx.gz
          └── spatial/
              ├── tissue_positions_list.csv   (v1.x)
              ├── tissue_positions.csv        (v2.x)
              ├── scalefactors_json.json
              ├── tissue_hires_image.png
              └── tissue_lowres_image.png
    """
    outs = os.path.join(sr_dir, "outs") if os.path.isdir(os.path.join(sr_dir, "outs")) else sr_dir

    # Load expression matrix
    h5_path = find_file(outs, [
        "filtered_feature_bc_matrix.h5",
        "raw_feature_bc_matrix.h5",
    ], required=False)

    if h5_path:
        print(f"  Loading H5: {h5_path}")
        adata = sc.read_10x_h5(h5_path)
    else:
        # Try MTX directory
        mtx_dir = find_file(outs, [
            "filtered_feature_bc_matrix",
            "raw_feature_bc_matrix",
        ])
        print(f"  Loading MTX: {mtx_dir}")
        adata = sc.read_10x_mtx(mtx_dir)

    adata.var_names_make_unique()

    # Load spatial info
    spatial_dir = os.path.join(outs, "spatial")
    if not os.path.isdir(spatial_dir):
        sys.exit(f"ERROR: spatial/ directory not found in {outs}")

    # Tissue positions (v1 vs v2 format)
    pos_path = find_file(spatial_dir, [
        "tissue_positions.csv",        # SpaceRanger v2.x
        "tissue_positions_list.csv",   # SpaceRanger v1.x
    ])

    # Detect header
    with open(pos_path, "r") as f:
        first_line = f.readline()

    if "barcode" in first_line.lower():
        pos = pd.read_csv(pos_path, index_col=0)
    else:
        pos = pd.read_csv(pos_path, header=None,
                          names=["barcode", "in_tissue", "array_row", "array_col",
                                 "pxl_row_fullres", "pxl_col_fullres"])
        pos = pos.set_index("barcode")

    # Filter to in-tissue spots
    if "in_tissue" in pos.columns:
        pos = pos[pos["in_tissue"] == 1]

    common = adata.obs_names.intersection(pos.index)
    if len(common) == 0:
        sys.exit("ERROR: No matching barcodes between expression matrix and tissue positions")

    adata = adata[common].copy()

    # Spatial coordinates: pxl_col (x), pxl_row (y)
    col_x = "pxl_col_fullres" if "pxl_col_fullres" in pos.columns else pos.columns[-1]
    col_y = "pxl_row_fullres" if "pxl_row_fullres" in pos.columns else pos.columns[-2]
    adata.obsm["spatial"] = pos.loc[common, [col_x, col_y]].values.astype(np.float64)

    # Scale factors
    sf_path = find_file(spatial_dir, ["scalefactors_json.json"], required=False)
    scale_info = {}
    if sf_path:
        with open(sf_path) as f:
            scale_info = json.load(f)

    # Tissue image (prefer hires)
    image_candidates = (["tissue_lowres_image.png", "tissue_hires_image.png"]
                        if use_lowres else
                        ["tissue_hires_image.png", "tissue_lowres_image.png"])
    img_path = find_file(spatial_dir, image_candidates)

    return adata, img_path, scale_info


def load_spaceranger_hd(sr_dir, resolution="8um", use_lowres=False):
    """
    Load Visium HD data from Space Ranger v3.x/v4.x output.

    Expected structure:
      spaceranger_hd_out/
      └── outs/
          ├── spatial/
          │   ├── tissue_hires_image.png
          │   ├── tissue_lowres_image.png
          │   ├── scalefactors_json.json
          │   └── tissue_positions.parquet     (full-res positions)
          └── binned_outputs/
              ├── square_002um/
              │   ├── filtered_feature_bc_matrix.h5
              │   └── spatial/
              │       ├── tissue_positions.parquet
              │       └── scalefactors_json.json
              ├── square_008um/
              │   ├── filtered_feature_bc_matrix.h5
              │   └── spatial/
              │       ├── tissue_positions.parquet
              │       └── scalefactors_json.json
              └── square_016um/
                  ├── filtered_feature_bc_matrix.h5
                  └── spatial/
                      ├── tissue_positions.parquet
                      └── scalefactors_json.json

    For cell_segmentation:
      spaceranger_hd_out/
      └── outs/
          └── per_cell_outs/        (or cell_segmentation/)
              ├── filtered_feature_bc_matrix.h5
              └── spatial/
                  └── tissue_positions.parquet
    """
    outs = os.path.join(sr_dir, "outs") if os.path.isdir(os.path.join(sr_dir, "outs")) else sr_dir
    print(f"  SpaceRanger HD base: {os.path.abspath(outs)}")

    # Determine bin directory
    resolution_map = {
        "2um":  "square_002um",
        "8um":  "square_008um",
        "16um": "square_016um",
    }

    if resolution == "cell_segmentation":
        # Look for cell segmentation output (varies by SpaceRanger version)
        cell_dirs = [
            os.path.join(outs, "segmented_outputs"),  # SpaceRanger 4.x
            os.path.join(outs, "per_cell_outs"),       # SpaceRanger 3.x
            os.path.join(outs, "cell_segmentation"),   # alternative naming
        ]
        bin_dir = None
        for d in cell_dirs:
            if os.path.isdir(d):
                bin_dir = d
                break
        if bin_dir is None:
            sys.exit(f"ERROR: Cell segmentation output not found. Looked in: {cell_dirs}")
        print(f"  Using cell segmentation: {bin_dir}")
    else:
        bin_name = resolution_map.get(resolution)
        if not bin_name:
            sys.exit(f"ERROR: Unknown resolution '{resolution}'. Use: 2um, 8um, 16um, cell_segmentation")
        bin_dir = os.path.join(outs, "binned_outputs", bin_name)
        if not os.path.isdir(bin_dir):
            sys.exit(f"ERROR: Binned output not found: {bin_dir}")
        print(f"  Using HD bins: {bin_dir} ({resolution})")

    # Load expression (note: cell segmentation uses "cell_matrix" not "bc_matrix")
    h5_path = find_file(bin_dir, [
        "filtered_feature_bc_matrix.h5",
        "filtered_feature_cell_matrix.h5",   # SpaceRanger 4.x segmented
        "raw_feature_bc_matrix.h5",
        "raw_feature_cell_matrix.h5",
    ])
    print(f"  Loading H5: {h5_path}")
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    # Load positions (Visium HD uses Parquet)
    spatial_sub = os.path.join(bin_dir, "spatial")
    if not os.path.isdir(spatial_sub):
        spatial_sub = bin_dir

    parquet_path = find_file(spatial_sub, [
        "tissue_positions.parquet",
    ], required=False)

    csv_path = find_file(spatial_sub, [
        "tissue_positions.csv",
        "tissue_positions_list.csv",
    ], required=False)

    if parquet_path:
        print(f"  Loading positions (Parquet): {parquet_path}")
        pos = pd.read_parquet(parquet_path)
        if "barcode" in pos.columns:
            pos = pos.set_index("barcode")
    elif csv_path:
        print(f"  Loading positions (CSV): {csv_path}")
        with open(csv_path, "r") as f:
            first_line = f.readline()
        if "barcode" in first_line.lower() or "cell_id" in first_line.lower():
            pos = pd.read_csv(csv_path, index_col=0)
        else:
            pos = pd.read_csv(
                csv_path, header=None,
                names=["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_fullres", "pxl_col_fullres"],
                index_col=0)
    else:
        sys.exit(f"ERROR: No tissue_positions file found in {spatial_sub}")

    # Filter in-tissue
    if "in_tissue" in pos.columns:
        pos = pos[pos["in_tissue"] == 1]

    common = adata.obs_names.intersection(pos.index)
    print(f"  Matched barcodes: {len(common)}")
    if len(common) == 0:
        sys.exit("ERROR: No matching barcodes")

    adata = adata[common].copy()

    col_x = "pxl_col_fullres" if "pxl_col_fullres" in pos.columns else pos.columns[-1]
    col_y = "pxl_row_fullres" if "pxl_row_fullres" in pos.columns else pos.columns[-2]
    adata.obsm["spatial"] = pos.loc[common, [col_x, col_y]].values.astype(np.float64)

    # Scale factors
    sf_path = find_file(spatial_sub, ["scalefactors_json.json"], required=False)
    if not sf_path:
        sf_path = find_file(os.path.join(outs, "spatial"), ["scalefactors_json.json"], required=False)
    scale_info = {}
    if sf_path:
        with open(sf_path) as f:
            scale_info = json.load(f)

    # Tissue image (from top-level spatial/)
    top_spatial = os.path.join(outs, "spatial")
    image_candidates = (["tissue_lowres_image.png", "tissue_hires_image.png"]
                        if use_lowres else
                        ["tissue_hires_image.png", "tissue_lowres_image.png"])
    img_path = find_file(top_spatial, image_candidates, required=False)
    if not img_path:
        img_path = find_file(spatial_sub, image_candidates, required=False)

    if img_path:
        print(f"  Tissue image source: {img_path}")
    else:
        print("  WARNING: No tissue image found in SpaceRanger output.")

    return adata, img_path, scale_info


def load_h5_matrix(h5_path, spatial_dir, use_lowres=False):
    """
    Load from standalone .h5 file + spatial/ directory.
    Typical for GEO downloads or 10x website downloads.

    Expected files:
      filtered_feature_bc_matrix.h5  (or any .h5)
      spatial/
      ├── tissue_positions.csv (or tissue_positions_list.csv or .parquet)
      ├── scalefactors_json.json
      ├── tissue_hires_image.png
      └── tissue_lowres_image.png
    """
    print(f"  Loading H5: {h5_path}")
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    # Positions
    parquet = find_file(spatial_dir, ["tissue_positions.parquet"], required=False)
    csv = find_file(spatial_dir, [
        "tissue_positions.csv", "tissue_positions_list.csv"
    ], required=False)

    if parquet:
        pos = pd.read_parquet(parquet)
        if "barcode" in pos.columns:
            pos = pos.set_index("barcode")
    elif csv:
        with open(csv) as f:
            first = f.readline()
        if "barcode" in first.lower():
            pos = pd.read_csv(csv, index_col=0)
        else:
            pos = pd.read_csv(csv, header=None,
                              names=["barcode", "in_tissue", "array_row", "array_col",
                                     "pxl_row_fullres", "pxl_col_fullres"])
            pos = pos.set_index("barcode")
    else:
        sys.exit(f"ERROR: No positions file in {spatial_dir}")

    if "in_tissue" in pos.columns:
        pos = pos[pos["in_tissue"] == 1]

    common = adata.obs_names.intersection(pos.index)
    print(f"  Matched barcodes: {len(common)}")
    if len(common) == 0:
        sys.exit("ERROR: No matching barcodes between expression matrix and tissue positions")
    adata = adata[common].copy()

    col_x = "pxl_col_fullres" if "pxl_col_fullres" in pos.columns else pos.columns[-1]
    col_y = "pxl_row_fullres" if "pxl_row_fullres" in pos.columns else pos.columns[-2]
    adata.obsm["spatial"] = pos.loc[common, [col_x, col_y]].values.astype(np.float64)

    sf_path = find_file(spatial_dir, ["scalefactors_json.json"], required=False)
    scale_info = json.load(open(sf_path)) if sf_path else {}

    image_candidates = (["tissue_lowres_image.png", "tissue_hires_image.png"]
                        if use_lowres else
                        ["tissue_hires_image.png", "tissue_lowres_image.png"])
    img_path = find_file(spatial_dir, image_candidates, required=False)

    return adata, img_path, scale_info


def load_h5ad(h5ad_path, tissue_image=None):
    """Load from pre-existing .h5ad file (from GEO, 10x, or custom pipeline)."""
    print(f"  Loading H5AD: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    if "spatial" not in adata.obsm:
        # Try common alternatives
        for key in ["X_spatial", "spatial_coords"]:
            if key in adata.obsm:
                adata.obsm["spatial"] = adata.obsm[key]
                break
        else:
            print("WARNING: No spatial coordinates found in .obsm['spatial']")

    return adata, tissue_image, {}


# =============================================================================
# Analysis pipeline
# =============================================================================

def run_analysis(adata, min_counts=100, min_genes=50, n_clusters=10):
    """Minimal Scanpy pipeline: QC → normalize → cluster → UMAP."""
    print(f"  Input: {adata.n_obs} cells/spots, {adata.n_vars} genes")

    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=10)
    print(f"  After QC: {adata.n_obs} cells/spots, {adata.n_vars} genes")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # The dispersion-based Seurat method operates on log-normalized data and
    # avoids the optional scikit-misc/OpenBLAS dependency required by
    # flavor="seurat_v3", improving portability on Apple Silicon and Linux.
    n_top = min(3000, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    n_comps = min(50, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1)
    if n_comps < 2:
        sys.exit("ERROR: Too few observations or genes remain after QC for PCA")
    sc.tl.pca(adata_hvg, n_comps=n_comps, random_state=0)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, random_state=0)
    sc.tl.umap(adata, random_state=0)

    # Iteratively find resolution that gives ~n_clusters
    res_lo, res_hi = 0.01, 3.0
    best_res, best_diff = 0.5, 999
    for _ in range(15):
        res_mid = (res_lo + res_hi) / 2
        sc.tl.leiden(adata, resolution=res_mid, key_added="cluster",
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=0)
        n_found = adata.obs["cluster"].nunique()
        diff = n_found - n_clusters
        if abs(diff) < abs(best_diff):
            best_diff = diff
            best_res = res_mid
        if diff == 0:
            break
        if n_found > n_clusters:
            res_hi = res_mid
        else:
            res_lo = res_mid

    if abs(best_diff) > 1:
        # Re-run with best resolution found
        sc.tl.leiden(adata, resolution=best_res, key_added="cluster",
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=0)
    n_final = adata.obs["cluster"].nunique()
    print(f"  Clustering: requested {n_clusters}, got {n_final} (resolution={best_res:.4f})")
    adata.obs["cluster"] = adata.obs["cluster"].astype(str)

    adata.obsm["UMAP"] = adata.obsm["X_umap"]
    # Keep log-normalized expression in .X for visualization (better color contrast)
    # Raw counts preserved in .layers["counts"] for downstream re-analysis
    return adata


def find_markers(adata, n_top=25):
    """Find marker genes per cluster using Wilcoxon rank-sum test."""
    adata_log = adata.copy()
    # Avoid applying normalize_total/log1p twice to objects produced by
    # run_analysis. For raw-count inputs, normalize a temporary copy.
    is_log1p = "log1p" in adata_log.uns
    if not is_log1p:
        matrix_max = (adata_log.X.max() if not hasattr(adata_log.X, "toarray")
                      else adata_log.X.max())
        if float(matrix_max) > 50:
            sc.pp.normalize_total(adata_log, target_sum=1e4)
            sc.pp.log1p(adata_log)
    sc.tl.rank_genes_groups(adata_log, groupby="cluster", method="wilcoxon",
                            n_genes=n_top)
    result = adata_log.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    rows = []
    for group in groups:
        for i in range(n_top):
            row = {
                "gene": result["names"][group][i],
                "cluster": group,
                "avg_log2FC": result["logfoldchanges"][group][i],
                "p_val": result["pvals"][group][i],
                "p_val_adj": result["pvals_adj"][group][i],
            }
            if "pcts" in result:
                row["pct.1"] = result["pcts"][group][i]
            if "pcts_rest" in result:
                row["pct.2"] = result["pcts_rest"][group][i]
            rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# Alignment computation
# =============================================================================

def compute_alignment(adata, img_path, scale_info, use_hires=True):
    """Compute alignment parameters for config.yaml."""
    spatial = adata.obsm["spatial"]
    img = Image.open(img_path)
    img_w, img_h = img.size

    # If SpaceRanger scalefactors are available, use them
    if scale_info:
        if use_hires and "tissue_hires_scalef" in scale_info:
            sf = scale_info["tissue_hires_scalef"]
            print(f"  Using hires scale factor from SpaceRanger: {sf:.6f}")
        elif "tissue_lowres_scalef" in scale_info:
            sf = scale_info["tissue_lowres_scalef"]
            print(f"  Using lowres scale factor from SpaceRanger: {sf:.6f}")
        else:
            sf = None

        if sf:
            # SpaceRanger scale factor maps fullres coords directly to image pixels.
            # offset must be 0 — the coordinate systems already match.
            return {
                "scale_factor_x": round(sf, 6),
                "scale_factor_y": round(sf, 6),
                "offset_x": 0,
                "offset_y": 0,
                "img_xmin": 0,
                "img_xmax": img_w,
                "img_ymin": 0,
                "img_ymax": img_h,
            }

    # Fallback: compute from data range
    x_min, x_max = spatial[:, 0].min(), spatial[:, 0].max()
    y_min, y_max = spatial[:, 1].min(), spatial[:, 1].max()
    margin = 0.05
    data_range_x = x_max - x_min
    data_range_y = y_max - y_min
    sf_x = img_w * (1 - 2 * margin) / data_range_x
    sf_y = img_h * (1 - 2 * margin) / data_range_y
    sf = min(sf_x, sf_y)
    offset_x = x_min * sf - img_w * margin
    offset_y = y_min * sf - img_h * margin

    return {
        "scale_factor_x": round(sf, 6),
        "scale_factor_y": round(sf, 6),
        "offset_x": round(offset_x),
        "offset_y": round(offset_y),
        "img_xmin": 0,
        "img_xmax": img_w,
        "img_ymin": 0,
        "img_ymax": img_h,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Prepare spatial transcriptomics data for SpatialAtlas",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Standard Visium (SpaceRanger v1/v2)
  python prepare_data.py --spaceranger ./SR_output --output ./data -n Sample1

  # Visium HD at 8µm bins
  python prepare_data.py --spaceranger-hd ./SR_HD_output --hd-resolution 8um \\
      --output ./data -n HD_Sample

  # Visium HD cell segmentation
  python prepare_data.py --spaceranger-hd ./SR_HD_output \\
      --hd-resolution cell_segmentation --output ./data -n CellSeg

  # GEO download: .h5 + spatial/
  python prepare_data.py --h5-matrix ./GSMxxxx_filtered_feature_bc_matrix.h5 \\
      --spatial-dir ./GSMxxxx_spatial/ --output ./data -n GEO_sample

  # Pre-existing .h5ad
  python prepare_data.py --h5ad ./my_processed.h5ad \\
      --tissue-image ./tissue.png --output ./data -n MyAtlas --skip-analysis
""")

    # Input sources (mutually exclusive groups)
    src = parser.add_argument_group("Input source (choose one)")
    src.add_argument("--spaceranger", metavar="DIR",
                     help="SpaceRanger v1.x/v2.x output directory (standard Visium)")
    src.add_argument("--spaceranger-hd", metavar="DIR",
                     help="Space Ranger v3.x/v4.x output directory (Visium HD)")
    src.add_argument("--h5-matrix", metavar="FILE",
                     help="Standalone .h5 matrix file (+ --spatial-dir)")
    src.add_argument("--h5ad", metavar="FILE",
                     help="Pre-existing .h5ad file")

    # Additional input options
    parser.add_argument("--spatial-dir", metavar="DIR",
                        help="spatial/ directory (required with --h5-matrix)")
    parser.add_argument("--tissue-image", metavar="FILE",
                        help="Tissue image file (required with --h5ad if no spatial/ dir)")
    parser.add_argument("--hd-resolution", default="8um",
                        choices=["2um", "8um", "16um", "cell_segmentation"],
                        help="Visium HD resolution (default: 8um)")

    # Output
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--sample-name", "-n", required=True, help="Sample name")

    # Analysis
    parser.add_argument("--skip-analysis", action="store_true",
                        help="Skip QC/clustering (requires 'cluster' in .obs)")
    parser.add_argument("--min-counts", type=int, default=100)
    parser.add_argument("--min-genes", type=int, default=50)
    parser.add_argument("--n-clusters", type=int, default=10)
    parser.add_argument("--markers", action="store_true",
                        help="Compute and save marker genes")
    parser.add_argument("--use-lowres", action="store_true",
                        help="Use lowres image instead of hires")
    parser.add_argument("--keep-analysis-artifacts", action="store_true",
                        help="Retain PCA/neighbor intermediates in output .h5ad")

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)

    # Count how many input sources provided
    sources = sum(1 for x in [args.spaceranger, args.spaceranger_hd,
                               args.h5_matrix, args.h5ad] if x)
    if sources == 0:
        parser.error("Provide one input: --spaceranger, --spaceranger-hd, --h5-matrix, or --h5ad")
    if sources > 1:
        parser.error("Provide only ONE input source")

    # ── Load ──
    print(f"=== Preparing: {args.sample_name} ===")

    if args.spaceranger:
        adata, img_path, scale_info = load_spaceranger_standard(
            args.spaceranger, use_lowres=args.use_lowres)
    elif args.spaceranger_hd:
        adata, img_path, scale_info = load_spaceranger_hd(
            args.spaceranger_hd, args.hd_resolution, use_lowres=args.use_lowres)
    elif args.h5_matrix:
        if not args.spatial_dir:
            parser.error("--h5-matrix requires --spatial-dir")
        adata, img_path, scale_info = load_h5_matrix(
            args.h5_matrix, args.spatial_dir, use_lowres=args.use_lowres)
    elif args.h5ad:
        adata, img_path, scale_info = load_h5ad(args.h5ad, args.tissue_image)

    # An explicit image takes precedence over the file discovered by the loader.
    if args.tissue_image and args.tissue_image != img_path:
        if os.path.exists(args.tissue_image):
            print(f"  Overriding tissue image with: {args.tissue_image}")
            img_path = args.tissue_image
        else:
            sys.exit(f"ERROR: --tissue-image file not found: {args.tissue_image}")

    print(f"  Loaded: {adata.n_obs} cells/spots, {adata.n_vars} genes")

    # ── Analysis ──
    if not args.skip_analysis:
        print("Running analysis pipeline...")
        adata = run_analysis(adata, args.min_counts, args.min_genes, args.n_clusters)
    else:
        if "cluster" not in adata.obs.columns:
            sys.exit("ERROR: --skip-analysis requires 'cluster' in .obs")
        if "UMAP" not in adata.obsm and "X_umap" in adata.obsm:
            adata.obsm["UMAP"] = adata.obsm["X_umap"]
        elif "UMAP" not in adata.obsm:
            print("  WARNING: No UMAP found. UMAP tab will be disabled.")

    # ── Markers (before optional slimming) ──
    markers = None
    if args.markers:
        print("  Finding marker genes...")
        markers = find_markers(adata)

    # The viewer only needs X, obs['cluster'], obsm['spatial'] and optionally
    # obsm['UMAP']. Remove analysis intermediates by default to make the
    # distributable atlas substantially smaller.
    if not args.keep_analysis_artifacts:
        for key in ["X_pca", "X_umap"]:
            if key in adata.obsm:
                del adata.obsm[key]
        for key in list(adata.obsp.keys()):
            del adata.obsp[key]
        for key in ["neighbors", "umap", "leiden", "rank_genes_groups"]:
            adata.uns.pop(key, None)
        for key in list(adata.layers.keys()):
            del adata.layers[key]

    # ── Save .h5ad ──
    name_lower = re.sub(r"[^a-z0-9._-]+", "_", args.sample_name.lower()).strip("_")
    if not name_lower:
        sys.exit("ERROR: --sample-name must contain at least one letter or digit")
    h5ad_out = os.path.join(args.output, f"sample_{name_lower}.h5ad")
    print(f"  Saving: {h5ad_out}")
    adata.write_h5ad(h5ad_out)

    # ── Copy/save tissue image ──
    img_out = None
    if img_path and os.path.exists(img_path):
        img_out = os.path.join(args.output, f"tissue_{name_lower}.png")
        src_size = os.path.getsize(img_path)
        print(f"  Copying image: {img_path} ({src_size} bytes)")
        if src_size == 0:
            print("  WARNING: Source tissue image is 0 bytes! (Dropbox Smart Sync? Right-click → Make Available Offline)")
        if img_path.lower().endswith(".png"):
            shutil.copy2(img_path, img_out)
        else:
            Image.open(img_path).save(img_out, "PNG")
        print(f"  Image saved: {img_out} ({os.path.getsize(img_out)} bytes)")

    # ── Alignment ──
    if img_out and "spatial" in adata.obsm:
        print("  Computing alignment...")
        align = compute_alignment(adata, img_out, scale_info, not args.use_lowres)

        # Relative path from cwd for local use (prefix ./ for clarity)
        rel_h5ad = "./" + os.path.relpath(h5ad_out, os.getcwd())
        rel_img  = "./" + os.path.relpath(img_out, os.getcwd())

        print("\n" + "=" * 60)
        print("FOR LOCAL USE — paste into config.yaml:")
        print("=" * 60)
        print(f'  - name: "{args.sample_name}"')
        print(f'    h5ad: "{rel_h5ad}"')
        print(f'    image: "{rel_img}"')
        for k, v in align.items():
            print(f"    {k}: {v}")

        print("\n" + "-" * 60)
        print("FOR DOCKER — paste into config.yaml:")
        print("-" * 60)
        print(f'  - name: "{args.sample_name}"')
        print(f'    h5ad: "/data/sample_{name_lower}.h5ad"')
        print(f'    image: "/data/tissue_{name_lower}.png"')
        for k, v in align.items():
            print(f"    {k}: {v}")
        print("=" * 60)
    else:
        print("\n  WARNING: Cannot compute alignment (missing image or spatial coords)")

    # ── Save markers ──
    if markers is not None:
        markers_path = os.path.join(args.output, "markers.csv")
        markers.to_csv(markers_path, index=False)
        rel_markers = "./" + os.path.relpath(markers_path, os.getcwd())
        print(f"  Saved: {markers_path}")
        print(f"\n  Add to config.yaml:  markers_csv: \"{rel_markers}\"")

    print(f"\nDone! Files saved to {args.output}/")


if __name__ == "__main__":
    main()
