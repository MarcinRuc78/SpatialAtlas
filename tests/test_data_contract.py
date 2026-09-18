"""Configuration handling, the input data contract, and HTTP routes."""

from __future__ import annotations

import numpy as np
import pytest
import yaml

from conftest import SHARED_GENES, load_app


def test_relative_paths_resolve_from_the_config_file(app_module, atlas_dir):
    assert app_module.CONFIG_DIR == str(atlas_dir)
    resolved = app_module.resolve_config_path("markers.csv")
    assert resolved == str(atlas_dir / "markers.csv")


def test_absolute_paths_are_left_untouched(app_module):
    assert app_module.resolve_config_path("/tmp/x.h5ad") == "/tmp/x.h5ad"


def test_samples_expose_the_required_contract(app_module):
    for sample in app_module.samples:
        assert sample["coords"].shape[1] == 2
        assert np.isfinite(sample["coords"]).all()
        assert sample["clusters"].shape[0] == sample["n_obs"]
        assert sample["X"].shape[0] == sample["n_obs"]
        assert len(sample["gene_index"]) == sample["X"].shape[1]


def test_expression_matrices_stay_sparse_and_column_oriented(app_module):
    import scipy.sparse as sp
    for sample in app_module.samples:
        assert sp.issparse(sample["X"])
        assert sample["X"].format == "csc"


def test_combined_metadata_is_consistent(app_module):
    total = sum(s["n_obs"] for s in app_module.samples)
    assert app_module.n_obs_combined == total
    assert app_module.clusters_combined.shape[0] == total
    assert app_module.sample_names_combined.shape[0] == total


def test_sample_slice_partitions_the_combined_vector(app_module):
    covered = 0
    for index in range(len(app_module.samples)):
        start, end = app_module.sample_slice(index)
        assert start == covered
        covered = end
    assert covered == app_module.n_obs_combined


def test_image_mimetypes(app_module):
    assert app_module.image_mimetype("a.PNG") == "image/png"
    assert app_module.image_mimetype("a.webp") == "image/webp"
    assert app_module.image_mimetype("a.tif") == "image/tiff"
    assert app_module.image_mimetype("a.unknown") is None


def test_missing_cluster_field_is_rejected(tmp_path, atlas_dir):
    import anndata as ad

    broken = tmp_path / "broken.h5ad"
    source = ad.read_h5ad(atlas_dir / "sample_a" / "sample_a.h5ad")
    del source.obs["cluster"]
    source.write_h5ad(broken)

    module = load_app(atlas_dir / "config.yaml", "spatialatlas_contract_probe")
    with pytest.raises(ValueError, match=r"cluster"):
        module.load_sample({
            "name": "Broken", "h5ad": str(broken),
            "image": str(atlas_dir / "sample_a" / "sample_a.png"),
        })


def test_missing_file_is_reported(app_module, atlas_dir):
    with pytest.raises(FileNotFoundError):
        app_module.load_sample({
            "name": "Absent", "h5ad": str(atlas_dir / "nope.h5ad"),
            "image": str(atlas_dir / "sample_a" / "sample_a.png"),
        })


def test_more_than_four_samples_is_rejected(tmp_path, atlas_dir):
    config = yaml.safe_load((atlas_dir / "config.yaml").read_text())
    config["samples"] = config["samples"] * 3  # six entries
    for sample in config["samples"]:
        sample["h5ad"] = str(atlas_dir / sample["h5ad"])
        sample["image"] = str(atlas_dir / sample["image"])
    config["markers_csv"] = str(atlas_dir / "markers.csv")
    path = tmp_path / "too_many.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ValueError, match="between 1 and 4 samples"):
        load_app(path, "spatialatlas_too_many")


def test_tissue_routes(app_module):
    client = app_module.server.test_client()
    response = client.get("/_spatialatlas/tissue/0")
    assert response.status_code == 200
    assert response.headers["Content-Type"].startswith("image/")
    assert "max-age=3600" in response.headers.get("Cache-Control", "")

    archival = client.get("/_spatialatlas/tissue-full/1")
    assert archival.status_code == 200

    assert client.get("/_spatialatlas/tissue/99").status_code == 404


def test_application_root_serves_the_layout(app_module):
    client = app_module.server.test_client()
    response = client.get("/")
    assert response.status_code == 200
    assert b"react-entry-point" in response.data


def test_gene_search_is_prefix_first(app_module):
    results = app_module.gene_search_results("aa", None)
    assert results[0]["value"] == SHARED_GENES[0]


def test_gene_search_keeps_the_current_selection(app_module):
    results = app_module.gene_search_results("zzz", SHARED_GENES[0])
    assert results[0]["value"] == SHARED_GENES[0]
