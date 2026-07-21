# Adrenal Visium and Visium HD benchmark protocol

Run commands from the package root. The benchmark was executed on an Apple M4 Pro host with all four applications forced to `linux/amd64` containers.

## 1. Verify retained inputs

```bash
test -s examples/05_adrenal_platform_benchmark/source_data/adrenal_visium_seurat.rds
test -s examples/05_adrenal_platform_benchmark/source_app/app.R
test -s examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad
shasum -a 256 examples/05_adrenal_platform_benchmark/source_data/adrenal_visium_seurat.rds
```

The expected RDS digest is `cd338178ff885d5d6eb7080b66bdaf22356b5b619324d51422cc88cb6dbbcb94`.

## 2. Export the normalized Seurat data layer and image

```bash
docker run --rm --platform linux/amd64 \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e \
  Rscript /work/scripts/export_adrenal_visium.R \
  /work/source_data/adrenal_visium_seurat.rds \
  /work/prepared/intermediate
```

Expected output: 992 profiles, 19,465 genes and a 600 × 546 low-resolution histology image.

## 3. Build the matched SpatialAtlas object

Use an environment containing `anndata`, `numpy`, `pandas` and `scipy`:

```bash
python examples/05_adrenal_platform_benchmark/scripts/build_adrenal_visium_h5ad.py
```

The script retains the normalized Seurat `data` layer; it does not normalize the matrix a second time.

## 4. Build the matched Visium HD Shiny object

The conversion uses the H5AD CSR arrays directly as a gene-by-bin `dgCMatrix`. It does not normalize, filter or densify the values.

```bash
docker run --rm --platform linux/amd64 \
  -e ADRENAL_HD_H5AD=/package/examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad \
  -e ADRENAL_HD_SHINY_RDS=/work/source_data/adrenal_visium_hd_shiny.rds \
  -v "$PWD:/package:ro" \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e \
  Rscript /work/scripts/build_adrenal_hd_shiny_object.R
```

Expected output: 143,112 profiles, 15,109 genes and a compressed sparse RDS of approximately 51 MiB.

Verify that the complete sparse arrays, identifiers and coordinates are identical:

```bash
docker run --rm --platform linux/amd64 \
  -e ADRENAL_HD_H5AD=/package/examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad \
  -e ADRENAL_HD_SHINY_RDS=/work/source_data/adrenal_visium_hd_shiny.rds \
  -v "$PWD:/package:ro" \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e \
  Rscript /work/scripts/validate_adrenal_hd_shiny_object.R
```

The machine-readable result is written to `benchmark/hd_object_equivalence.json`.

## 5. Validate both SpatialAtlas inputs

```bash
python validation/validate_artifacts.py \
  --config examples/05_adrenal_platform_benchmark/config.classic.docker.yaml

python validation/validate_artifacts.py \
  --config examples/05_adrenal_platform_benchmark/config.hd.docker.yaml
```

When running the YAML files directly, mount the classic `prepared` directory or the Visium HD `output` directory at `/data`.

## 6. Measure the original classic Seurat/Shiny response path

```bash
docker run --rm --platform linux/amd64 \
  -e SHINY_SEURAT_RDS=/work/source_data/adrenal_visium_seurat.rds \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  -w /work/benchmark/reference_app \
  satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e \
  Rscript /work/benchmark/benchmark_shiny_core.R
```

The script warms each marker once and then records eight measurements for each of `Npy`, `Cyp11b1`, `Cyp11b2` and `Th`. Each response is an 800 × 800 Cairo PNG.

## 7. Measure the SpatialAtlas classic Visium response path

```bash
docker run --rm --platform linux/amd64 \
  -e SPATIALATLAS_CONFIG=/work/config.classic.docker.yaml \
  -e BENCHMARK_LABEL=spatialatlas_visium \
  -e BENCHMARK_OUTPUT=/work/benchmark/spatialatlas_classic_core_results.csv \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  -v "$PWD/examples/05_adrenal_platform_benchmark/prepared:/data:ro" \
  -w /app spatialatlas-benchmark:amd64 \
  python /work/benchmark/benchmark_spatialatlas_core.py
```

## 8. Measure the matched Visium HD Shiny response path

```bash
docker run --rm --platform linux/amd64 \
  -e SHINY_HD_RDS=/work/source_data/adrenal_visium_hd_shiny.rds \
  -e SHINY_HD_IMAGE=/data/tissue_adrenal_female_hd_8um.png \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  -v "$PWD/examples/02_visium_hd_8um/output:/data:ro" \
  -w /work/benchmark/reference_app_hd \
  satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e \
  Rscript /work/benchmark/benchmark_shiny_hd_core.R
```

The reference application uses the same 143,112 expression profiles, coordinates and histology as SpatialAtlas and renders a 900 × 900 Cairo PNG on the server. Its complete source and deterministic H5AD-to-RDS conversion are included with the benchmark.

## 9. Measure the SpatialAtlas Visium HD response path

```bash
docker run --rm --platform linux/amd64 \
  -e SPATIALATLAS_CONFIG=/work/config.hd.docker.yaml \
  -e BENCHMARK_LABEL=spatialatlas_visium_hd \
  -e BENCHMARK_OUTPUT=/work/benchmark/spatialatlas_hd_core_results.csv \
  -v "$PWD/examples/05_adrenal_platform_benchmark:/work" \
  -v "$PWD/examples/02_visium_hd_8um/output:/data:ro" \
  -w /app spatialatlas-benchmark:amd64 \
  python /work/benchmark/benchmark_spatialatlas_core.py
```

SpatialAtlas serializes Plotly JSON for client-side WebGL rendering. Shiny renders a PNG on the server. Response time and payload therefore characterize the actual application paths and are not interchangeable graphics-engine benchmarks.

## 10. Repeat startup and memory measurements

```bash
python examples/05_adrenal_platform_benchmark/benchmark/benchmark_startup_memory.py
```

Six starts are recorded per scenario. The first is retained as the cold observation; the reported summary uses the median and minimum–maximum range of the next five. Memory is read from `memory.current` two seconds after HTTP readiness.

## 11. Aggregate results and regenerate benchmark plots

```bash
python examples/05_adrenal_platform_benchmark/benchmark/summarize_results.py
python examples/05_adrenal_platform_benchmark/benchmark/plot_operational_benchmark.py
```

The plotting script writes PNG and SVG versions of the operational benchmark and the classic-Visium-versus-Visium-HD adrenal comparison to `benchmark/plots/`.

## 12. Prepare the SpatialAtlas display image

The archival PNG remains unchanged. The browser path uses a full-dimension WebP at quality 90; the preparation script verifies dimensions, measures PSNR, and records both SHA-256 digests.

```bash
python examples/05_adrenal_platform_benchmark/benchmark/prepare_display_image.py \
  examples/02_visium_hd_8um/output/tissue_adrenal_female_hd_8um.png \
  examples/05_adrenal_platform_benchmark/prepared/spatialatlas/tissue_adrenal_female_hd_8um_display.webp \
  examples/05_adrenal_platform_benchmark/benchmark/display_image_validation.json
```

The supplied HD configuration points `display_image` to this WebP and retains `image` as the archival PNG.

## 13. Prepare the independent TissUUmaps object

Use the retained writable benchmark directory because TissUUmaps creates a CSC copy and URL-specific image pyramids. The adapter adds only `library_id`, the existing histology and its scale factor; it does not change expression values, identifiers or coordinates.

```bash
BENCH_DIR="$PWD/examples/05_adrenal_platform_benchmark/prepared/tissuumaps"
mkdir -p "$BENCH_DIR"

docker run --rm --platform linux/amd64 \
  -e HDF5_DISABLE_VERSION_CHECK=1 \
  -v "$PWD:/package:ro" -v "$BENCH_DIR:/bench" \
  cavenel/tissuumaps@sha256:e2f7898c700cd550432cfe546c14546681a4235c27b40c49cbb019e97ab86c06 \
  python /package/examples/05_adrenal_platform_benchmark/benchmark/prepare_tissuumaps_h5ad.py \
  /package/examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad \
  /package/examples/02_visium_hd_8um/output/tissue_adrenal_female_hd_8um.png \
  /bench/adrenal_hd_tissuumaps.h5ad
```

Start the official image and request the adapter URL once:

```bash
docker run -d --platform linux/amd64 --name spatialatlas-tissuumaps \
  -p 127.0.0.1:18051:80 -e HDF5_DISABLE_VERSION_CHECK=1 \
  -v "$BENCH_DIR:/mnt/data/shared" \
  cavenel/tissuumaps@sha256:e2f7898c700cd550432cfe546c14546681a4235c27b40c49cbb019e97ab86c06

curl -f http://127.0.0.1:18051/adrenal_hd_tissuumaps.h5ad >/dev/null
```

For deployed and repeated tests, open `http://127.0.0.1:18051/adrenal_hd_tissuumaps_tmap.h5ad`; requesting the original CSR URL again repeats conversion in version 3.2.1.14.

## 14. Validate the TissUUmaps conversion

```bash
docker run --rm --platform linux/amd64 \
  -e HDF5_DISABLE_VERSION_CHECK=1 \
  -v "$PWD:/package" -v "$BENCH_DIR:/bench" \
  cavenel/tissuumaps@sha256:e2f7898c700cd550432cfe546c14546681a4235c27b40c49cbb019e97ab86c06 \
  python /package/examples/05_adrenal_platform_benchmark/benchmark/validate_tissuumaps_equivalence.py \
  /package/examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad \
  /bench/adrenal_hd_tissuumaps.h5ad /bench/adrenal_hd_tissuumaps_tmap.h5ad \
  /package/examples/05_adrenal_platform_benchmark/benchmark/tissuumaps_object_equivalence.json
```

Expected: all adapter digests pass; the reconstructed sparse-matrix difference has zero nonzero values; identifiers and spatial coordinates are identical.

## 15. Measure independent-viewer startup and memory

```bash
python examples/05_adrenal_platform_benchmark/benchmark/benchmark_tissuumaps_startup.py \
  --prepared-dir "$BENCH_DIR" --runs 6 \
  --output examples/05_adrenal_platform_benchmark/benchmark/tissuumaps_startup_memory_results.json
```

The first start is retained but excluded; the summary uses the following five runs. Readiness is HTTP 200 for the prepared HD viewer document, not the empty TissUUmaps landing page.

## 16. Browser readiness

Start SpatialAtlas on port 18052 with the HD configuration and retain TissUUmaps on port 18051. In Chromium, reload each prepared viewer five times. Record wall time until:

1. SpatialAtlas displays the Plotly WebGL canvas containing the 143,112-profile cluster map after its initial callback;
2. TissUUmaps displays the OpenSeadragon view and the 15,109-gene expression selector.

The recorded values and complete definition are in `benchmark/browser_ready_results.json`. They include loopback transfer, JavaScript parsing, layout and rendering on the tested browser/GPU path. They do not isolate GPU time.

## 17. Concurrent first-visit bundles

With the HD SpatialAtlas container listening on port 18052:

```bash
python examples/05_adrenal_platform_benchmark/benchmark/benchmark_concurrent_users.py \
  --base-url http://127.0.0.1:18052 --rounds 5 \
  --output examples/05_adrenal_platform_benchmark/benchmark/concurrent_users_results.json
```

Each uncached virtual session requests the HTML shell, Dash layout, dependency graph and configured display image. The test uses 1, 5, 10 and 20 concurrent sessions and reports errors, median, 95th percentile, throughput and bytes transferred.

## 18. Generate the independent-viewer plot

```bash
python examples/05_adrenal_platform_benchmark/benchmark/summarize_independent_viewer.py
python examples/05_adrenal_platform_benchmark/benchmark/plot_independent_viewer_benchmark.py
```

The script produces GraphPad-style PNG and SVG files in `benchmark/plots/`. Machine-readable feature classifications and source URLs are in `benchmark/viewer_feature_comparison.csv` and `benchmark/viewer_feature_sources.md`.
