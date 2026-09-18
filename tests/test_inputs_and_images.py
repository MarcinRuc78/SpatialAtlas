"""Input discovery, image handling, and the optional display-image contract."""

from __future__ import annotations

import numpy as np
import pytest
from PIL import Image

from _cli import load_script


@pytest.fixture(scope="module")
def prepare_module():
    return load_script("prepare_data")


def test_find_file_returns_the_first_candidate(prepare_module, tmp_path):
    (tmp_path / "second.h5").write_bytes(b"x")
    found = prepare_module.find_file(
        tmp_path, ["first.h5", "second.h5"], required=True)
    assert found.endswith("second.h5")


def test_find_file_optional_returns_none(prepare_module, tmp_path):
    assert prepare_module.find_file(tmp_path, ["absent.h5"], required=False) is None


def test_find_file_required_exits(prepare_module, tmp_path):
    with pytest.raises(SystemExit):
        prepare_module.find_file(tmp_path, ["absent.h5"], required=True)


def test_load_h5ad_warns_when_spatial_coordinates_are_absent(
        prepare_module, tmp_path, capsys):
    import anndata as ad
    import scipy.sparse as sp

    adata = ad.AnnData(X=sp.csr_matrix(np.ones((5, 4), dtype=np.float32)))
    adata.var_names = [f"G{i}" for i in range(4)]
    path = tmp_path / "no_spatial.h5ad"
    adata.write_h5ad(path)
    prepare_module.load_h5ad(path)
    assert "No spatial coordinates" in capsys.readouterr().out


def test_display_image_must_preserve_aspect_ratio(app_module, atlas_dir, tmp_path):
    wrong = tmp_path / "wrong_ratio.png"
    Image.new("RGB", (200, 50)).save(wrong)
    with pytest.raises(ValueError, match="aspect ratio"):
        app_module.load_sample({
            "name": "Skewed",
            "h5ad": str(atlas_dir / "sample_a" / "sample_a.h5ad"),
            "image": str(atlas_dir / "sample_a" / "sample_a.png"),
            "display_image": str(wrong),
        })


def test_display_image_of_matching_ratio_is_accepted(app_module, atlas_dir, tmp_path):
    smaller = tmp_path / "small.png"
    Image.new("RGB", (32, 32)).save(smaller)
    sample = app_module.load_sample({
        "name": "Scaled",
        "h5ad": str(atlas_dir / "sample_a" / "sample_a.h5ad"),
        "image": str(atlas_dir / "sample_a" / "sample_a.png"),
        "display_image": str(smaller),
    })
    # Coordinate bounds follow the archival image, not the display copy.
    assert sample["img_bounds"]["xmax"] == 64
    assert sample["display_image_path"] == str(smaller)


def test_coordinates_use_the_configured_scale_and_offset(app_module, atlas_dir):
    import anndata as ad

    source = ad.read_h5ad(atlas_dir / "sample_a" / "sample_a.h5ad")
    sample = app_module.load_sample({
        "name": "Aligned",
        "h5ad": str(atlas_dir / "sample_a" / "sample_a.h5ad"),
        "image": str(atlas_dir / "sample_a" / "sample_a.png"),
        "scale_factor_x": 0.5, "scale_factor_y": 0.25,
        "offset_x": 10, "offset_y": 4,
    })
    spatial = np.asarray(source.obsm["spatial"])
    expected_x = spatial[:, 0] * 0.5 - 10
    expected_y = spatial[:, 1] * 0.25 - 4
    assert np.allclose(sample["coords"][:, 0], expected_x, atol=1e-3)
    assert np.allclose(sample["coords"][:, 1], expected_y, atol=1e-3)
