"""Every visualization callback must run and must respect gene availability."""

from __future__ import annotations

import json

import numpy as np
import pytest

from conftest import SAMPLE_A_ONLY, SHARED_GENES

SHARED = SHARED_GENES[0]


def _annotation_texts(figure):
    return [a.text or "" for a in figure.layout.annotations]


def test_spatial_cluster_mode_returns_traces(app_module):
    figure = app_module.update_spatial(None, 3, 0.8, "Viridis", "tab-spatial")
    assert len(figure.data) > 0
    assert len(figure.layout.images) == 2


def test_spatial_expression_mode_returns_traces(app_module):
    figure = app_module.update_spatial(SHARED, 3, 0.8, "Viridis", "tab-spatial")
    assert len(figure.data) == 2
    assert all(np.isfinite(np.asarray(trace.marker.color)).all()
               for trace in figure.data)


def test_spatial_marks_sample_without_the_gene(app_module):
    figure = app_module.update_spatial(SAMPLE_A_ONLY, 3, 0.8, "Viridis",
                                       "tab-spatial")
    # Only the sample that retains the gene contributes a point cloud; the
    # other subplot gets an invisible extent placeholder so its histology
    # stays framed.
    visible = [trace for trace in figure.data if trace.marker.opacity != 0]
    assert len(visible) == 1
    placeholder = [trace for trace in figure.data if trace.marker.opacity == 0]
    assert len(placeholder) == 1
    texts = " ".join(_annotation_texts(figure))
    assert "no data for this gene in this sample" in texts
    assert "no data for" in texts  # subplot title of the absent sample


def test_spatial_colorbar_follows_the_available_sample(app_module):
    figure = app_module.update_spatial(SAMPLE_A_ONLY, 3, 0.8, "Viridis",
                                       "tab-spatial")
    assert any(trace.marker.showscale for trace in figure.data)


def test_unavailable_subplot_keeps_the_histology_framed(app_module):
    """A subplot without data must still show its tissue image."""
    figure = app_module.update_spatial(SAMPLE_A_ONLY, 3, 0.8, "Viridis",
                                       "tab-spatial")
    assert len(figure.layout.images) == 2
    placeholder = next(t for t in figure.data if t.marker.opacity == 0)
    bounds = app_module.samples[1]["img_bounds"]
    assert list(placeholder.x) == [bounds["xmin"], bounds["xmax"]]
    assert list(placeholder.y) == [bounds["ymin"], bounds["ymax"]]


def test_spatial_ignores_inactive_tab(app_module):
    from dash import no_update
    assert app_module.update_spatial(SHARED, 3, 0.8, "Viridis", "tab-umap") is no_update


def test_umap_separates_unavailable_observations(app_module):
    figure = app_module.update_umap(SAMPLE_A_ONLY, "tab-umap")
    names = [trace.name or "" for trace in figure.data]
    assert any("no data for" in name for name in names)
    grey = figure.data[0]
    assert len(grey.x) == app_module.samples[1]["n_obs"]


def test_umap_for_shared_gene_has_no_unavailable_trace(app_module):
    figure = app_module.update_umap(SHARED, "tab-umap")
    names = [trace.name or "" for trace in figure.data]
    assert not any("no data for" in name for name in names)


def test_violin_omits_sample_without_the_gene(app_module):
    figure = app_module.update_violin(SAMPLE_A_ONLY, "tab-violin")
    plotted = {trace.name for trace in figure.data}
    assert plotted == {"Sample A"}
    texts = " ".join(_annotation_texts(figure))
    assert f"is missing from Sample B" in texts


def test_violin_uses_both_samples_for_a_shared_gene(app_module):
    figure = app_module.update_violin(SHARED, "tab-violin")
    assert {trace.name for trace in figure.data} == {"Sample A", "Sample B"}
    assert not _annotation_texts(figure)


def test_violin_values_are_finite(app_module):
    figure = app_module.update_violin(SAMPLE_A_ONLY, "tab-violin")
    for trace in figure.data:
        assert np.isfinite(np.asarray(trace.y, dtype=float)).all()


def test_heatmap_flags_partial_genes(app_module):
    figure = app_module.update_heatmap(3, "tab-heatmap")
    labels = list(figure.data[0].x)
    assert any(label.endswith(" *") for label in labels)
    assert "not present in every sample" in figure.layout.title.text
    assert np.isfinite(np.asarray(figure.data[0].z, dtype=float)).all()


def test_multi_gene_overlay_marks_absent_channels(app_module):
    enables = [["on"]] * 4
    opacities = [0.8] * 4
    figure = app_module.update_multi(
        SHARED_GENES[0], SHARED_GENES[1], SHARED_GENES[2], SAMPLE_A_ONLY,
        *enables, *opacities, 3, "tab-multiple")
    texts = " ".join(_annotation_texts(figure))
    assert f"no data here for: {SAMPLE_A_ONLY}" in texts


def test_multi_gene_overlay_without_partial_genes_is_unannotated(app_module):
    enables = [["on"]] * 4
    opacities = [0.8] * 4
    figure = app_module.update_multi(
        SHARED_GENES[0], SHARED_GENES[1], SHARED_GENES[2], SHARED_GENES[3],
        *enables, *opacities, 3, "tab-multiple")
    texts = " ".join(_annotation_texts(figure))
    assert "no data here for" not in texts


@pytest.mark.parametrize("gene", [SHARED, SAMPLE_A_ONLY])
def test_every_figure_serializes(app_module, gene):
    figures = [
        app_module.update_spatial(gene, 3, 0.8, "Viridis", "tab-spatial"),
        app_module.update_umap(gene, "tab-umap"),
        app_module.update_violin(gene, "tab-violin"),
        app_module.update_heatmap(2, "tab-heatmap"),
    ]
    for figure in figures:
        payload = figure.to_json()
        assert json.loads(payload)


def test_single_sample_atlas_has_no_umap_tab(single_app_module):
    assert single_app_module.has_umap is False
    figure = single_app_module.update_spatial(SHARED, 3, 0.8, "Viridis",
                                              "tab-spatial")
    assert len(figure.data) == 1
