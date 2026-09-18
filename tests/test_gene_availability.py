"""A gene absent from a sample matrix must never be shown as measured zero."""

from __future__ import annotations

import numpy as np

from conftest import SAMPLE_A_ONLY, SAMPLE_B_ONLY, SHARED_GENES


def test_union_gene_set_is_exposed(app_module):
    assert set(app_module.gene_set) == set(
        SHARED_GENES + [SAMPLE_A_ONLY, SAMPLE_B_ONLY])


def test_partial_genes_are_identified(app_module):
    assert app_module.partial_genes == {SAMPLE_A_ONLY, SAMPLE_B_ONLY}
    for gene in SHARED_GENES:
        assert gene not in app_module.partial_genes


def test_presence_map_matches_sample_matrices(app_module):
    assert app_module.gene_presence[SAMPLE_A_ONLY] == (True, False)
    assert app_module.gene_presence[SAMPLE_B_ONLY] == (False, True)
    assert app_module.gene_presence[SHARED_GENES[0]] == (True, True)


def test_missing_sample_names(app_module):
    assert app_module.missing_samples_for(SAMPLE_A_ONLY) == ["Sample B"]
    assert app_module.missing_samples_for(SAMPLE_B_ONLY) == ["Sample A"]
    assert app_module.missing_samples_for(SHARED_GENES[0]) == []


def test_absent_gene_returns_nan_not_zero(app_module):
    expr = app_module.get_expr(SAMPLE_A_ONLY)
    start, end = app_module.sample_slice(1)
    absent_block = expr[start:end]
    assert np.isnan(absent_block).all(), "absent sample must not be filled with 0"
    present_block = expr[:start]
    assert np.isfinite(present_block).all()
    assert (present_block >= 0).all()


def test_zero_fill_remains_available_on_request(app_module):
    expr = app_module.get_expr(SAMPLE_A_ONLY, missing=0.0)
    start, end = app_module.sample_slice(1)
    assert (expr[start:end] == 0).all()


def test_unknown_gene_is_not_zero_expression(app_module):
    expr = app_module.get_expr("NotAGeneAtAll")
    assert np.isnan(expr).all()


def test_vmax_ignores_unavailable_observations(app_module):
    expr = app_module.get_expr(SAMPLE_A_ONLY)
    vmax = app_module.expression_vmax(expr)
    measured = expr[np.isfinite(expr)]
    expected = max(float(np.percentile(measured[measured > 0], 99)), 0.1)
    assert np.isclose(vmax, expected)


def test_dropdown_label_states_availability(app_module):
    assert app_module.gene_option_label(SHARED_GENES[0]) == SHARED_GENES[0]
    label = app_module.gene_option_label(SAMPLE_A_ONLY)
    assert SAMPLE_A_ONLY in label and "1 of 2 samples" in label


def test_availability_banner_only_for_partial_genes(app_module):
    assert app_module.gene_availability_note(SHARED_GENES[0]) is None
    assert app_module.gene_availability_note(None) is None
    banner = app_module.gene_availability_note(SAMPLE_A_ONLY)
    assert banner is not None
    text = str(banner)
    assert "Sample A only" not in text
    assert "Sample B" in text          # the sample that lacks the gene
    assert "is missing from" in text


def test_natural_list_reads_as_a_sentence(app_module):
    assert app_module.natural_list([]) == ""
    assert app_module.natural_list(["A"]) == "A"
    assert app_module.natural_list(["A", "B"]) == "A and B"
    assert app_module.natural_list(["A", "B", "C"]) == "A, B and C"


def test_single_sample_atlas_has_no_partial_genes(single_app_module):
    assert single_app_module.partial_genes == set()
