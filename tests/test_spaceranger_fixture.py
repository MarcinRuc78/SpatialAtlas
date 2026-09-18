"""Space Ranger loaders exercised on a synthetic miniature Visium run.

The fixture writes a minimal but structurally faithful Space Ranger directory —
a 10x-format HDF5 matrix, tissue positions, scale factors and a tissue image —
so that the standard-Visium and MEX code paths are covered without shipping a
real Space Ranger output in the test suite.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import h5py
import numpy as np
import pytest
import scipy.sparse as sp
from PIL import Image

from _cli import load_script

N_SPOTS = 24
N_GENES = 16


@pytest.fixture(scope="module")
def prepare_module():
    return load_script("prepare_data")


def _write_10x_h5(path: Path, matrix: sp.csc_matrix, barcodes, genes):
    """Write a CellRanger v3 style filtered_feature_bc_matrix.h5."""
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=matrix.data.astype(np.int32))
        group.create_dataset("indices", data=matrix.indices.astype(np.int64))
        group.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        group.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int32))
        group.create_dataset("barcodes",
                             data=np.array(barcodes, dtype="S"))
        features = group.create_group("features")
        features.create_dataset("id", data=np.array(
            [f"ENSTEST{i:08d}" for i in range(len(genes))], dtype="S"))
        features.create_dataset("name", data=np.array(genes, dtype="S"))
        features.create_dataset("feature_type", data=np.array(
            ["Gene Expression"] * len(genes), dtype="S"))
        features.create_dataset("genome", data=np.array(
            ["synthetic"] * len(genes), dtype="S"))


@pytest.fixture(scope="module")
def spaceranger_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("spaceranger")
    outs = root / "outs"
    spatial = outs / "spatial"
    spatial.mkdir(parents=True)

    rng = np.random.default_rng(101)
    barcodes = [f"AAACAAGTATCTCCCA-{i}" for i in range(N_SPOTS)]
    genes = [f"Gene{i:02d}" for i in range(N_GENES)]
    counts = rng.poisson(5.0, size=(N_SPOTS, N_GENES)).astype(np.int32)
    # 10x stores genes × barcodes in CSC form.
    _write_10x_h5(outs / "filtered_feature_bc_matrix.h5",
                  sp.csc_matrix(counts.T), barcodes, genes)

    # Space Ranger v2 positions, with a header and two out-of-tissue spots.
    rows = ["barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"]
    for index, barcode in enumerate(barcodes):
        in_tissue = 0 if index >= N_SPOTS - 2 else 1
        rows.append(f"{barcode},{in_tissue},{index},{index},"
                    f"{100 + index * 7},{200 + index * 5}")
    (spatial / "tissue_positions.csv").write_text("\n".join(rows) + "\n")

    (spatial / "scalefactors_json.json").write_text(json.dumps({
        "tissue_hires_scalef": 0.2,
        "tissue_lowres_scalef": 0.05,
        "spot_diameter_fullres": 60.0,
    }))
    Image.new("RGB", (300, 200), color=(200, 180, 180)).save(
        spatial / "tissue_hires_image.png")
    Image.new("RGB", (75, 50), color=(200, 180, 180)).save(
        spatial / "tissue_lowres_image.png")
    return root


def test_standard_loader_reads_matrix_positions_and_image(
        prepare_module, spaceranger_run):
    adata, img_path, scale_info = prepare_module.load_spaceranger_standard(
        spaceranger_run)
    # Out-of-tissue spots are dropped.
    assert adata.n_obs == N_SPOTS - 2
    assert adata.n_vars == N_GENES
    assert adata.obsm["spatial"].shape == (N_SPOTS - 2, 2)
    assert np.isfinite(adata.obsm["spatial"]).all()
    assert img_path.endswith("tissue_hires_image.png")
    assert scale_info["tissue_hires_scalef"] == pytest.approx(0.2)


def test_standard_loader_can_prefer_the_lowres_image(
        prepare_module, spaceranger_run):
    _, img_path, _ = prepare_module.load_spaceranger_standard(
        spaceranger_run, use_lowres=True)
    assert img_path.endswith("tissue_lowres_image.png")


def test_headerless_v1_positions_are_accepted(prepare_module, spaceranger_run,
                                              tmp_path):
    """Space Ranger 1.x wrote tissue_positions_list.csv without a header."""
    import shutil

    root = tmp_path / "v1"
    shutil.copytree(spaceranger_run, root)
    spatial = root / "outs" / "spatial"
    source = (spatial / "tissue_positions.csv").read_text().splitlines()[1:]
    (spatial / "tissue_positions_list.csv").write_text("\n".join(source) + "\n")
    (spatial / "tissue_positions.csv").unlink()

    adata, _, _ = prepare_module.load_spaceranger_standard(root)
    assert adata.n_obs == N_SPOTS - 2


def test_missing_spatial_directory_is_reported(prepare_module, spaceranger_run,
                                               tmp_path):
    import shutil

    root = tmp_path / "no_spatial"
    shutil.copytree(spaceranger_run, root)
    shutil.rmtree(root / "outs" / "spatial")
    with pytest.raises(SystemExit, match="spatial"):
        prepare_module.load_spaceranger_standard(root)


def test_full_standard_workflow_produces_a_viewable_atlas(
        prepare_module, spaceranger_run, tmp_path, capsys):
    """Loader → QC/clustering → alignment, on the synthetic run."""
    adata, img_path, scale_info = prepare_module.load_spaceranger_standard(
        spaceranger_run)
    result = prepare_module.run_analysis(adata, min_counts=5, min_genes=2,
                                         n_clusters=2)
    assert "cluster" in result.obs and "UMAP" in result.obsm
    alignment = prepare_module.compute_alignment(result, img_path, scale_info)
    assert alignment["scale_factor_x"] == pytest.approx(0.2)
    assert alignment["img_xmax"] == 300 and alignment["img_ymax"] == 200
