#!/usr/bin/env python3
"""
integrate_samples.py — joint clustering of SpatialAtlas samples
================================================================
Integrates prepared samples so that they share clusters and UMAP coordinates,
while retaining the full expression matrix of each input object.

Usage:
  python integrate_samples.py \
    --samples ./data/sample1/sample_brain_cells.h5ad \
              ./data/sample2/sample_brain_cells2.h5ad \
    --n-clusters 8 \
    --markers

Steps:
  1. Load all .h5ad files
  2. Find common genes across all samples
  3. Normalize raw-count inputs on the common gene set; retain log-normalized inputs
  4. Run PCA and Harmony batch correction (or the requested fallback)
  5. Joint Leiden clustering to target --n-clusters
  6. Joint UMAP embedding
  7. Write integrated copies to a new directory (or overwrite only with --in-place)
  8. Optionally compute shared marker genes
"""

import argparse
import os
import sys
import tempfile
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp


def main():
    parser = argparse.ArgumentParser(
        description="Joint integration and clustering of SpatialAtlas samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--samples", "-s", nargs="+", required=True,
                        help="Paths to .h5ad files from prepare_data.py")
    parser.add_argument("--n-clusters", "-n", type=int, default=10,
                        help="Target number of clusters (default: 10)")
    parser.add_argument("--markers", action="store_true",
                        help="Compute and save shared markers.csv")
    parser.add_argument("--output-markers", default=None,
                        help="Path for markers.csv (default: same dir as first sample)")
    parser.add_argument("--output-dir", default=None,
                        help="Directory for integrated copies (default: integrated/ beside first sample)")
    parser.add_argument("--in-place", action="store_true",
                        help="Overwrite inputs atomically instead of creating integrated copies")
    args = parser.parse_args()
    if args.in_place and args.output_dir:
        parser.error("Use either --in-place or --output-dir, not both")

    # ── Load all samples ──
    adatas = []
    sample_names = []
    sample_n_cells = []  # track cell count per sample for write-back
    for path in args.samples:
        if not os.path.exists(path):
            sys.exit(f"ERROR: File not found: {path}")
        print(f"Loading: {path}")
        adata = sc.read_h5ad(path)
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        # Extract sample name from filename
        name = os.path.basename(path).replace("sample_", "").replace(".h5ad", "")
        adata.obs["_sample"] = name
        sample_names.append(name)
        sample_n_cells.append(adata.n_obs)
        adatas.append(adata)
        print(f"  {name}: {adata.n_obs} cells, {adata.n_vars} genes")

    # ── Find common genes ──
    common_genes = set(adatas[0].var_names)
    for a in adatas[1:]:
        common_genes &= set(a.var_names)
    common_genes = sorted(common_genes)
    print(f"\nCommon genes across {len(adatas)} samples: {len(common_genes)}")

    if len(common_genes) < 100:
        sys.exit("ERROR: Too few common genes (<100). Check that samples are from the same species/reference.")

    # Subset to common genes
    for i in range(len(adatas)):
        adatas[i] = adatas[i][:, common_genes].copy()

    # ── Concatenate ──
    adata = ad.concat(adatas, label="_sample", keys=sample_names, join="inner")
    adata.obs_names_make_unique()
    print(f"Combined: {adata.n_obs} cells, {adata.n_vars} genes")

    # ── Re-normalize jointly ──
    print("\nRunning joint analysis pipeline...")

    # Check if data is already log-transformed
    max_val = adata.X.max()
    if max_val > 50:
        # Raw counts — normalize
        print("  Normalizing (target sum 10,000 + log1p)...")
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("  Data appears log-normalized, using as-is.")
        adata.layers["counts"] = adata.X.copy()

    # HVG + PCA. Use the log-expression-compatible Seurat dispersion method
    # to avoid the optional scikit-misc/OpenBLAS dependency.
    print("  Finding HVGs...")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=min(3000, adata.n_vars), flavor="seurat",
        batch_key="_sample")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    n_comps = min(50, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1)
    if n_comps < 2:
        sys.exit("ERROR: Too few observations or genes remain for PCA")
    sc.tl.pca(adata_hvg, n_comps=n_comps, random_state=0)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    # ── Batch correction ──
    try:
        import harmonypy as hm
        print("  Batch correction: Harmony...")
        # A contiguous copy is required: PCA output can carry negative strides,
        # which the Torch backend of recent harmonypy releases rejects.
        pca_data = np.ascontiguousarray(adata.obsm["X_pca"], dtype=np.float32)
        meta = adata.obs[["_sample"]].copy()
        ho = hm.run_harmony(pca_data, meta, "_sample", max_iter_harmony=20)
        # Get corrected PCs — shape must be (n_cells, n_pcs)
        Z = ho.Z_corr
        if hasattr(Z, 'numpy'):
            Z = Z.numpy()  # PyTorch tensor → numpy
        Z = np.asarray(Z, dtype=np.float32)
        if Z.shape[0] != adata.n_obs:
            Z = Z.T
        adata.obsm["X_pca_harmony"] = Z
        use_rep = "X_pca_harmony"
        print(f"  Harmony done: {Z.shape}")
    except ImportError:
        try:
            import bbknn
            print("  Batch correction: BBKNN (harmonypy not installed)...")
            bbknn.bbknn(adata, batch_key="_sample", use_rep="X_pca")
            use_rep = None  # BBKNN already computed neighbors
        except ImportError:
            print("  WARNING: No batch correction (install harmonypy or bbknn for better results)")
            use_rep = "X_pca"

    # ── Neighbors + UMAP ──
    if use_rep is not None:
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=30, random_state=0)
    sc.tl.umap(adata, random_state=0)

    # ── Iterative Leiden clustering ──
    n_clusters = args.n_clusters
    print(f"  Clustering (target: {n_clusters})...")
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
        sc.tl.leiden(adata, resolution=best_res, key_added="cluster",
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=0)
    n_final = adata.obs["cluster"].nunique()
    print(f"  Clustering: requested {n_clusters}, got {n_final} (resolution={best_res:.4f})")
    adata.obs["cluster"] = adata.obs["cluster"].astype(str)

    # ── Write back to individual .h5ad files ──
    # Use positional slicing — ad.concat preserves within-sample order
    default_output_dir = os.path.join(os.path.dirname(args.samples[0]), "integrated")
    output_dir = args.output_dir or default_output_dir
    if not args.in_place:
        os.makedirs(output_dir, exist_ok=True)
        print(f"\nSaving non-destructive integrated copies to: {output_dir}")
    else:
        print("\nSaving integrated results in place (atomic replacement)...")
    offset = 0
    for i, (path, name, n_orig) in enumerate(zip(args.samples, sample_names, sample_n_cells)):
        # Load original to preserve spatial coords, image alignment and full gene set
        orig = sc.read_h5ad(path)
        orig.var_names_make_unique()
        orig.obs_names_make_unique()

        # Sanity check
        if orig.n_obs != n_orig:
            print(f"  WARNING: {name} cell count changed ({n_orig} → {orig.n_obs}). Skipping.")
            offset += n_orig
            continue

        # Extract cluster labels and UMAP by positional slice
        cluster_vals = adata.obs["cluster"].values[offset:offset + n_orig]
        umap_vals = np.asarray(adata.obsm["X_umap"][offset:offset + n_orig])

        orig.obs["cluster"] = cluster_vals
        orig.obsm["UMAP"] = umap_vals

        # Show cluster distribution for this sample
        from collections import Counter
        dist = Counter(cluster_vals)
        dist_str = ", ".join(f"cl{k}:{v}" for k, v in sorted(dist.items()))
        print(f"  {name}: {dist_str}")

        target_path = path if args.in_place else os.path.join(output_dir, os.path.basename(path))
        target_parent = os.path.dirname(os.path.abspath(target_path))
        fd, tmp_path = tempfile.mkstemp(prefix=".integrating-", suffix=".h5ad", dir=target_parent)
        os.close(fd)
        try:
            orig.write_h5ad(tmp_path)
            os.replace(tmp_path, target_path)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        print(f"  Saved: {target_path} ({n_orig} cells)")
        offset += n_orig

    # ── Shared markers ──
    if args.markers:
        print("\nFinding shared marker genes...")
        sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon",
                                use_raw=False, n_genes=25)
        result = adata.uns["rank_genes_groups"]
        rows = []
        for cluster in result["names"].dtype.names:
            for j in range(25):
                rows.append({
                    "cluster": cluster,
                    "gene": result["names"][cluster][j],
                    "score": result["scores"][cluster][j],
                    "pval_adj": result["pvals_adj"][cluster][j],
                    "logfoldchange": result["logfoldchanges"][cluster][j],
                })
        markers_df = pd.DataFrame(rows)
        markers_path = args.output_markers or os.path.join(
            os.path.dirname(args.samples[0]) if args.in_place else output_dir,
            "markers.csv")
        markers_df.to_csv(markers_path, index=False)
        rel_markers = "./" + os.path.relpath(markers_path, os.getcwd())
        print(f"  Saved: {markers_path}")
        print(f"\n  Add to config.yaml:  markers_csv: \"{rel_markers}\"")

    print(f"\nDone! All {len(adatas)} output samples share {n_final} joint clusters.")
    print("Point config.yaml to the integrated files, then restart app.py.")


if __name__ == "__main__":
    main()
