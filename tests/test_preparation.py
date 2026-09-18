"""Preparation workflow: QC filters, alignment, and the written data contract.

These tests use synthetic count matrices small enough to run in seconds. The
full Space Ranger and Visium HD workflows remain covered by the release-level
integration protocols in ``examples/`` and ``validation/``.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from PIL import Image

from _cli import load_script, run_main


@pytest.fixture(scope="module")
def prepare_module():
    return load_script("prepare_data")


def _counts_adata(n_obs=180, n_genes=60, seed=7):
    """Synthetic counts with a block of genes seen in fewer than ten spots."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(3.0, size=(n_obs, n_genes)).astype(np.float32)
    # The last five genes appear in exactly three observations each, i.e. below
    # the minimum-observation threshold used by the preparation workflow.
    counts[:, -5:] = 0
    counts[:3, -5:] = 4
    adata = ad.AnnData(X=sp.csr_matrix(counts))
    adata.var_names = [f"Gene{i:03d}" for i in range(n_genes)]
    adata.obs_names = [f"spot{i:04d}" for i in range(n_obs)]
    adata.obsm["spatial"] = rng.uniform(0, 1000, size=(n_obs, 2)).astype(np.float32)
    return adata


def test_rarely_detected_genes_are_removed(prepare_module):
    adata = _counts_adata()
    rare = list(adata.var_names[-5:])
    result = prepare_module.run_analysis(adata.copy(), min_counts=10,
                                         min_genes=5, n_clusters=3)
    for gene in rare:
        assert gene not in set(result.var_names), (
            "genes detected in fewer than ten observations must be dropped; "
            "their later absence is not a measured zero")
    assert result.n_vars < adata.n_vars


def test_analysis_writes_the_viewer_contract(prepare_module):
    result = prepare_module.run_analysis(_counts_adata(), min_counts=10,
                                         min_genes=5, n_clusters=3)
    assert "cluster" in result.obs
    assert result.obs["cluster"].dtype == object
    assert "UMAP" in result.obsm
    assert result.obsm["UMAP"].shape == (result.n_obs, 2)
    assert result.obsm["spatial"].shape == (result.n_obs, 2)
    assert np.isfinite(np.asarray(result.obsm["spatial"])).all()


def test_expression_is_log_normalized_not_raw(prepare_module):
    result = prepare_module.run_analysis(_counts_adata(), min_counts=10,
                                         min_genes=5, n_clusters=3)
    values = result.X.toarray() if sp.issparse(result.X) else np.asarray(result.X)
    assert values.max() < 15.0
    assert not np.allclose(values, np.round(values))


def test_markers_are_restricted_to_retained_genes(prepare_module):
    result = prepare_module.run_analysis(_counts_adata(), min_counts=10,
                                         min_genes=5, n_clusters=3)
    markers = prepare_module.find_markers(result, n_top=5)
    assert set(markers["gene"]).issubset(set(result.var_names))
    assert {"gene", "cluster"}.issubset(set(markers.columns))


def test_alignment_prefers_spaceranger_scale_factors(prepare_module, tmp_path):
    adata = _counts_adata(n_obs=20, n_genes=10)
    image = tmp_path / "tissue.png"
    Image.new("RGB", (500, 300)).save(image)
    alignment = prepare_module.compute_alignment(
        adata, image, {"tissue_hires_scalef": 0.25})
    assert alignment["scale_factor_x"] == pytest.approx(0.25)
    assert alignment["offset_x"] == 0 and alignment["offset_y"] == 0
    assert alignment["img_xmax"] == 500 and alignment["img_ymax"] == 300


def test_alignment_falls_back_to_the_data_range(prepare_module, tmp_path):
    adata = _counts_adata(n_obs=20, n_genes=10)
    image = tmp_path / "tissue.png"
    Image.new("RGB", (400, 400)).save(image)
    alignment = prepare_module.compute_alignment(adata, image, None)
    coords = np.asarray(adata.obsm["spatial"])
    x = coords[:, 0] * alignment["scale_factor_x"] - alignment["offset_x"]
    y = coords[:, 1] * alignment["scale_factor_y"] - alignment["offset_y"]
    assert x.min() >= -1 and x.max() <= 401
    assert y.min() >= -1 and y.max() <= 401


def test_h5ad_input_path_end_to_end(tmp_path, capsys):
    """A preprocessed H5AD passes through the CLI and yields a usable atlas."""
    rng = np.random.default_rng(11)
    n_obs, n_genes = 40, 12
    adata = ad.AnnData(X=sp.csr_matrix(rng.random((n_obs, n_genes), dtype=np.float32)))
    adata.var_names = [f"G{i}" for i in range(n_genes)]
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    adata.obs["cluster"] = [str(i % 2) for i in range(n_obs)]
    adata.obsm["spatial"] = rng.uniform(0, 200, size=(n_obs, 2)).astype(np.float32)
    source = tmp_path / "input.h5ad"
    adata.write_h5ad(source)

    image = tmp_path / "tissue.png"
    Image.new("RGB", (256, 256), color=(120, 90, 90)).save(image)

    output = tmp_path / "atlas"
    code, _ = run_main("prepare_data", [
        "--h5ad", str(source), "--tissue-image", str(image),
        "--output", str(output), "--sample-name", "synthetic",
        "--skip-analysis",
    ])
    assert code == 0
    stdout = capsys.readouterr().out

    written = list(output.rglob("*.h5ad"))
    assert written, f"no atlas H5AD written: {stdout}"
    result = ad.read_h5ad(written[0])
    assert "cluster" in result.obs
    assert "spatial" in result.obsm

    # The tissue image is copied next to the atlas and the alignment block is
    # emitted for config.yaml rather than written silently.
    assert (output / "tissue_synthetic.png").is_file()
    for key in ("scale_factor_x", "scale_factor_y", "offset_x", "img_xmax"):
        assert key in stdout
    assert "FOR DOCKER" in stdout
