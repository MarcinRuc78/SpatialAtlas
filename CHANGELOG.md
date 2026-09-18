# Changelog

All notable changes to SpatialAtlas are documented in this file.

## 1.1.0 — 2026-09-10

Maintenance release: explicit reporting of genes absent from a sample, and an automated test suite with continuous integration.

### Changed

- A gene that is absent from a sample's expression matrix is now reported as
  unavailable for that sample instead of being represented by zeros. Absence
  from a matrix means the gene was undetected or removed by the
  minimum-observation filter during preparation, which is not a measured
  expression of zero. The distinction is carried through every view:
  - the gene selector labels genes retained in only some samples and shows an
    explanatory banner naming the samples that lack the selected gene;
  - spatial and multi-gene panels for an unavailable sample are labelled rather
    than drawn as an all-zero profile;
  - UMAP draws such observations in a separate, uncolored "not available"
    layer;
  - violin plots omit the affected sample instead of adding a distribution at
    zero;
  - marker-heatmap cluster means are computed only from samples that measured
    the gene, and those genes are flagged on the axis.
- `get_expr()` now fills unavailable samples with NaN by default;
  `get_expr(gene, missing=0.0)` preserves the previous zero-fill for callers
  that need it.

### Added

- Automated test suite (`tests/`) with synthetic fixtures, covering gene
  availability, all visualization callbacks, the data contract and HTTP routes,
  the preparation workflow, a miniature synthetic Space Ranger run, and
  multi-sample integration.
- Continuous-integration workflow `.github/workflows/tests.yml` running the
  suite with coverage reporting, artifact upload, a job-summary coverage table,
  and a minimum coverage threshold.
- `TESTING.md` describing the two testing levels and current coverage.

### Fixed

- Multi-sample integration passed a possibly non-contiguous PCA array to
  Harmony; recent harmonypy releases with the Torch backend rejected it.

## 1.0.0 — 2026-07-21

- Added self-hosted visualization of processed Visium and Visium HD H5AD objects.
- Added spatial, UMAP, violin, marker-heatmap, and four-gene composite views.
- Added YAML configuration, Docker Compose deployment, and optional Caddy proxying.
- Added preparation workflows for standard Visium, binned Visium HD, cell-segmented Visium HD, standalone 10x matrices, and processed AnnData.
- Added multi-sample Harmony integration and coordinated one-to-four-sample layouts.
- Added browser-efficient tissue-image delivery, typed-array serialization, sparse column access, and deferred hidden-tab callbacks.
- Added five worked examples, validation scripts, matched Shiny benchmarks, an independent TissUUmaps benchmark, and concurrent-user measurements.
- Added reproducibility protocols, machine-readable reports, source-data inventories, and SHA-256 manifests.
