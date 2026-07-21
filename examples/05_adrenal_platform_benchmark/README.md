# Example 5 — adrenal platform and resolution benchmark

This example uses data from the authors' adrenal-gland studies and contains matched and independent application comparisons:

1. a matched implementation comparison of the same classic Visium adrenal object in the original Seurat/Shiny application and SpatialAtlas;
2. the same 143,112-bin, 8 µm Visium HD adrenal object in a purpose-built R/Shiny reference application and SpatialAtlas;
3. the same HD object in the official TissUUmaps Docker image, including browser-readiness and conversion-equivalence checks;
4. SpatialAtlas first-visit bundles at 1, 5, 10 and 20 concurrent users.

Within each pair, expression values, coordinates, histology and marker genes are identical. The HD reference is a reproducible server-rasterized R/Shiny implementation constructed directly from the validated H5AD object.

## Biological material

The classic Visium object contains four mouse adrenal sections and 992 tissue-covered spots across 19,465 genes. It accompanies:

> Blatkiewicz M, Hryhorowicz S, Szyszka M, et al. Molecular landscape of the mouse adrenal gland and adjacent adipose tissue by spatial transcriptomics. DOI: 10.5603/fhc.108988. GEO: GSE283302.

The Visium HD example contains 143,112 displayed 8 µm bins across 15,109 genes and is shared with `examples/02_visium_hd_8um`. It accompanies:

> Blatkiewicz M, Hryhorowicz S, Szyszka M, et al. Single-cell spatial transcriptomics reveals sex-dependent gene expression and intercellular signalling in mouse adrenal cortex. DOI: 10.1038/s42003-026-10697-9.

## Directory contents

```text
source_app/app.R                         application retained from the first adrenal study
source_data/adrenal_visium_seurat.rds   normalized Seurat object used by that application
source_data/adrenal_visium_hd_shiny.rds matched sparse R object built from the HD H5AD
prepared/                               matched classic H5AD/histology and retained viewer assets
prepared/spatialatlas/                  full-dimension WebP display image used in browser tests
prepared/tissuumaps/                    adapter H5AD, converted CSC H5AD and image pyramids
scripts/                                deterministic classic and HD conversion scripts
benchmark/reference_app/                Seurat 5 compatibility-adjusted application
benchmark/reference_app_hd/             purpose-built matched Visium HD Shiny application
benchmark/reference_app_hd/preview_*    rendered Cyp11b2 application preview
benchmark/*.csv                         32 per-gene measurements per implementation
benchmark/hd_object_equivalence.json    full sparse-array and coordinate identity check
benchmark/startup_memory_results.json   six container starts per scenario
benchmark/benchmark_summary.*           aggregate benchmark table
benchmark/prepare_tissuumaps_h5ad.py    deterministic image-metadata adapter
benchmark/prepare_display_image.py      deterministic WebP preparation and validation
benchmark/display_image_validation.json dimensions, PSNR, sizes and SHA-256 digests
benchmark/benchmark_data_inventory.csv retained data paths, sizes and SHA-256 digests
benchmark/build_data_inventory.py       deterministic inventory generator
benchmark/tissuumaps_object_equivalence.json complete CSR/CSC identity check
benchmark/tissuumaps_startup_memory_results.json independent-viewer measurements
benchmark/browser_ready_results.json    five browser runs per prepared viewer
benchmark/concurrent_users_results.json 1/5/10/20-user raw measurements
benchmark/independent_browser_summary.csv concise benchmark values
benchmark/summarize_independent_viewer.py deterministic independent-viewer aggregation
benchmark/viewer_feature_comparison.csv six-viewer feature matrix
benchmark/viewer_feature_sources.md     official source URLs and classification note
config.classic.docker.yaml              classic Visium SpatialAtlas configuration
config.hd.docker.yaml                   Visium HD SpatialAtlas configuration
protocol.md                             complete reproduction commands
provenance.json                         source identifiers and SHA-256 digests
```

The classic compatibility copy changes only the object path and Seurat 5 argument name (`slot` to `layer`) and loads Cairo for deterministic PNG response generation. The HD Shiny object preserves the H5AD sparse values without normalization or dense conversion. Neither Shiny path uses response caching. All TissUUmaps inputs and generated outputs used in the reported measurements are physically retained in `prepared/tissuumaps`; the equivalence report verifies the biological arrays after CSR-to-CSC conversion.

Biological data reuse remains subject to the originating studies, GEO record and institutional approvals. The software license does not redefine those terms.
