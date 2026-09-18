"""Multi-sample Harmony integration on a compact synthetic pair."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from _cli import run_main


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def sample_pair(tmp_path_factory):
    root = tmp_path_factory.mktemp("integration")
    paths = []
    # The two synthetic samples deliberately carry unequal gene sets, as
    # independently prepared sections do.
    for name, seed, first_gene in (("a", 21, 0), ("b", 22, 20)):
        rng = np.random.default_rng(seed)
        n_obs, n_genes = 120, 300
        counts = rng.poisson(3.0, size=(n_obs, n_genes)).astype(np.float32)
        adata = ad.AnnData(X=sp.csr_matrix(counts))
        adata.var_names = [f"Gene{first_gene + i:03d}" for i in range(n_genes)]
        adata.obs_names = [f"{name}{i:04d}" for i in range(n_obs)]
        adata.obs["cluster"] = [str(i % 2) for i in range(n_obs)]
        adata.obsm["spatial"] = rng.uniform(0, 500, size=(n_obs, 2)).astype(np.float32)
        path = root / f"sample_{name}.h5ad"
        adata.write_h5ad(path)
        paths.append(path)
    return root, paths


@pytest.fixture(scope="module")
def integrated(sample_pair):
    root, paths = sample_pair
    before = {p: _digest(p) for p in paths}
    output = root / "integrated"
    code, _ = run_main("integrate_samples", [
        "--samples", *map(str, paths),
        "--n-clusters", "3",
        "--output-dir", str(output),
    ])
    assert code == 0
    return output, paths, before


def test_inputs_are_not_modified_without_in_place(integrated):
    _, paths, before = integrated
    for path in paths:
        assert _digest(path) == before[path], "integration must not overwrite inputs"


def test_outputs_are_written_to_a_new_directory(integrated):
    output, paths, _ = integrated
    written = sorted(output.glob("*.h5ad"))
    assert len(written) == len(paths)


def test_integrated_objects_share_clusters_and_umap(integrated):
    output, _, _ = integrated
    written = sorted(output.glob("*.h5ad"))
    cluster_sets = []
    for path in written:
        adata = ad.read_h5ad(path)
        assert "cluster" in adata.obs
        assert "UMAP" in adata.obsm
        assert adata.obsm["UMAP"].shape == (adata.n_obs, 2)
        assert np.isfinite(np.asarray(adata.obsm["UMAP"])).all()
        assert "spatial" in adata.obsm
        cluster_sets.append(set(adata.obs["cluster"].astype(str)))
    assert cluster_sets[0] & cluster_sets[1], "samples must share cluster labels"


def test_integration_keeps_each_sample_gene_set(integrated):
    """Genes are intersected only for the shared embedding.

    The written objects retain the gene set each sample was prepared with, so
    an integrated atlas legitimately contains samples with different genes.
    The viewer therefore has to report per-sample availability rather than
    assume a common gene space.
    """
    output, paths, _ = integrated
    written = sorted(output.glob("*.h5ad"))
    written_sets = [set(ad.read_h5ad(p).var_names.astype(str)) for p in written]
    source_sets = [set(ad.read_h5ad(p).var_names.astype(str)) for p in paths]
    assert written_sets == source_sets
    assert written_sets[0] != written_sets[1]
    assert written_sets[0] & written_sets[1]


def test_in_place_and_output_dir_are_mutually_exclusive(sample_pair, tmp_path):
    _, paths = sample_pair
    code, _ = run_main("integrate_samples", [
        "--samples", *map(str, paths), "--in-place",
        "--output-dir", str(tmp_path / "out"),
    ])
    assert code != 0
