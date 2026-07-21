#!/usr/bin/env python3
"""Build the SpatialAtlas input from the normalized Seurat adrenal object."""

from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE = ROOT / "prepared" / "intermediate"
OUTPUT = ROOT / "prepared"

matrix = mmread(INTERMEDIATE / "expression_log1p.mtx").tocsr().astype(np.float32)
genes = pd.read_csv(INTERMEDIATE / "features.tsv", header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(INTERMEDIATE / "barcodes.tsv", header=None)[0].astype(str).tolist()
coordinates = pd.read_csv(INTERMEDIATE / "coordinates_lowres.csv").set_index("barcode").loc[barcodes]
identities = pd.read_csv(INTERMEDIATE / "identities.csv").set_index("barcode").loc[barcodes]

atlas = ad.AnnData(
    X=matrix,
    obs=pd.DataFrame(index=pd.Index(barcodes, name="barcode")),
    var=pd.DataFrame(index=pd.Index(genes, name="gene")),
)
atlas.obs["cluster"] = pd.Categorical(identities["cluster"].astype(str))
atlas.obs["sample"] = pd.Categorical(["adrenal_visium_classic"] * atlas.n_obs)
atlas.obsm["spatial"] = coordinates[["imagecol", "imagerow"]].to_numpy(dtype=np.float32)
atlas.uns["normalization"] = "Seurat data layer supplied by the first adrenal-atlas application"
atlas.write_h5ad(OUTPUT / "sample_adrenal_visium_classic.h5ad", compression="gzip")
shutil.copy2(
    INTERMEDIATE / "tissue_adrenal_visium.png",
    OUTPUT / "tissue_adrenal_visium_classic.png",
)

print(f"profiles={atlas.n_obs}")
print(f"genes={atlas.n_vars}")
print(f"matrix_nnz={atlas.X.nnz}")
print(f"image_width=600")
print(f"image_height=546")
