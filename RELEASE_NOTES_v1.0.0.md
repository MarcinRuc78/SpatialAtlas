# SpatialAtlas v1.0.0

SpatialAtlas v1.0.0 provides a reproducible route from processed Visium or Visium HD objects to a self-hosted interactive atlas.

Version-specific DOI: `10.5281/zenodo.21476516` (https://doi.org/10.5281/zenodo.21476516)

## Included

- Python/Dash application with spatial, UMAP, violin, heatmap, and four-gene views.
- Sparse H5AD access with one-gene-at-a-time materialization.
- YAML configuration and Docker Compose deployment.
- Standard Visium, Visium HD, cell-segmentation, and multi-sample preparation workflows.
- Five worked examples with dataset-specific step-by-step protocols.
- Deployment and data-access guidance in `docs/deployment.md`.
- Artifact, viewer, integration, Docker, and performance validation code with compact recorded reports.
- Matched Shiny, independent TissUUmaps, browser-readiness, and concurrent-session benchmarks.

## Release assets

Large biological inputs, prepared H5AD objects, retained Shiny objects, and TissUUmaps files are distributed as five versioned data archives rather than tracked in Git. Verify them with `release-assets.sha256` before use.
