"""Shared fixtures for the SpatialAtlas test suite.

Every fixture here builds a small synthetic atlas on disk so that the unit and
functional tests run in seconds and need no biological data. The large example
workflows remain covered by the release-level integration tests in
``validation/`` and ``examples/``.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
import yaml
from PIL import Image

SOFTWARE_DIR = Path(__file__).resolve().parents[1] / "software"

# Genes shared by both synthetic samples, plus one gene that only sample A
# retains. The exclusive gene reproduces the situation the viewer has to
# communicate: absent from a matrix is not the same as measured zero.
SHARED_GENES = ["Aaa1", "Bbb2", "Ccc3", "Ddd4", "Eee5", "Fff6"]
SAMPLE_A_ONLY = "OnlyInA"
SAMPLE_B_ONLY = "OnlyInB"


def _write_sample(folder: Path, name: str, genes: list[str], n_obs: int,
                  seed: int, with_umap: bool = True) -> Path:
    """Write one synthetic H5AD plus its tissue image."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    dense = rng.random((n_obs, len(genes))).astype(np.float32)
    dense[dense < 0.55] = 0.0  # keep the matrix genuinely sparse
    matrix = sp.csr_matrix(dense)

    obs_names = [f"{name}_{i:04d}" for i in range(n_obs)]
    adata = ad.AnnData(X=matrix)
    adata.obs_names = obs_names
    adata.var_names = genes
    adata.obs["cluster"] = [str(i % 3) for i in range(n_obs)]
    adata.obsm["spatial"] = rng.uniform(0, 400, size=(n_obs, 2)).astype(np.float32)
    if with_umap:
        adata.obsm["UMAP"] = rng.normal(0, 3, size=(n_obs, 2)).astype(np.float32)

    h5ad_path = folder / f"{name}.h5ad"
    adata.write_h5ad(h5ad_path)

    image_path = folder / f"{name}.png"
    Image.fromarray(
        rng.integers(0, 255, size=(64, 64, 3), dtype=np.uint8)
    ).save(image_path)
    return h5ad_path


def _write_markers(path: Path, genes: list[str]) -> None:
    lines = ["cluster,gene,avg_log2FC"]
    for index, gene in enumerate(genes):
        lines.append(f"{index % 3},{gene},{2.0 - index * 0.1:.2f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture(scope="session")
def atlas_dir(tmp_path_factory) -> Path:
    """A two-sample synthetic atlas with deliberately unequal gene sets."""
    root = tmp_path_factory.mktemp("atlas")
    _write_sample(root / "sample_a", "sample_a",
                  SHARED_GENES + [SAMPLE_A_ONLY], n_obs=90, seed=1)
    _write_sample(root / "sample_b", "sample_b",
                  SHARED_GENES + [SAMPLE_B_ONLY], n_obs=110, seed=2)
    _write_markers(root / "markers.csv", SHARED_GENES + [SAMPLE_A_ONLY])

    config = {
        "atlas": {"title": "Synthetic test atlas", "subtitle": "unit tests"},
        "markers_csv": "markers.csv",
        "samples": [
            {"name": "Sample A", "h5ad": "sample_a/sample_a.h5ad",
             "image": "sample_a/sample_a.png"},
            {"name": "Sample B", "h5ad": "sample_b/sample_b.h5ad",
             "image": "sample_b/sample_b.png"},
        ],
        "display": {"max_gene_dropdown": 50, "spatial_plot_height": 400},
        "multi_gene": {"gene_a": "Aaa1", "gene_b": "Bbb2",
                       "gene_c": "Ccc3", "gene_d": SAMPLE_A_ONLY},
        "citation": {"text": "Synthetic atlas used by the test suite."},
        "server": {"host": "127.0.0.1", "port": 8050},
    }
    (root / "config.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    return root


@pytest.fixture(scope="session")
def single_sample_dir(tmp_path_factory) -> Path:
    """A one-sample atlas without UMAP, to cover the optional-field paths."""
    root = tmp_path_factory.mktemp("atlas_single")
    _write_sample(root / "only", "only", SHARED_GENES, n_obs=60, seed=3,
                  with_umap=False)
    config = {
        "atlas": {"title": "Single sample"},
        "samples": [{"name": "Only", "h5ad": "only/only.h5ad",
                     "image": "only/only.png"}],
        "display": {"spatial_plot_height": 300},
    }
    (root / "config.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    return root


def load_app(config_path: Path, module_name: str):
    """Import app.py against a given configuration file."""
    import os

    os.environ["SPATIALATLAS_CONFIG"] = str(config_path.resolve())
    sys.modules.pop(module_name, None)
    spec = importlib.util.spec_from_file_location(
        module_name, SOFTWARE_DIR / "app.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def app_module(atlas_dir: Path):
    """The two-sample application, imported once for the whole session."""
    return load_app(atlas_dir / "config.yaml", "spatialatlas_under_test")


@pytest.fixture(scope="session")
def single_app_module(single_sample_dir: Path):
    return load_app(single_sample_dir / "config.yaml",
                    "spatialatlas_single_under_test")
