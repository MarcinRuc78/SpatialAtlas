# SpatialAtlas

SpatialAtlas is a self-hosted Dash/Plotly application for publishing processed Visium and Visium HD data as an interactive web atlas. It accepts sparse AnnData objects, keeps expression matrices on the server, and materializes only the selected gene vector for plotting. Dataset preparation, deployment, validation and benchmark procedures are documented alongside the code.

## Main capabilities

- standard Visium, Visium HD bins, Visium HD cell segmentation, standalone 10x matrices, and processed H5AD input;
- one to four coordinated samples configured in YAML;
- spatial cluster and expression maps, UMAP, violin plots, marker heatmaps, and a four-gene overlay;
- sparse expression access and WebGL point-cloud rendering;
- local Python, Gunicorn, Docker Compose, and optional Caddy deployment;
- explicit per-sample gene availability: a gene absent from one sample's matrix is reported as unavailable, never drawn or summarized as a measured zero;
- an automated test suite on synthetic fixtures, run in continuous integration with coverage reporting;
- reproducible artifact, viewer, implementation, browser-readiness, and concurrent-session validation.

## Repository layout

```text
software/     application, preparation and integration scripts, dependencies, and Docker files
examples/     configurations, protocols, benchmark code, and compact recorded reports
tests/        automated unit and functional tests on synthetic fixtures
validation/   artifact, viewer, and performance-validation scripts
docs/         deployment and data-access guidance
```

Large biological inputs and prepared objects are available as versioned release assets with SHA-256 digests.

## Installation

Use Python 3.11 in a clean virtual environment:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r software/requirements.txt
```

For application-only deployment, use `software/requirements-app.txt`.

## Prepare data

Standard Visium:

```bash
python software/prepare_data.py \
  --spaceranger /path/to/spaceranger/outs \
  --output /path/to/output \
  --sample-name sample_1 \
  --n-clusters 8 \
  --markers
```

Visium HD at 8 µm:

```bash
python software/prepare_data.py \
  --spaceranger-hd /path/to/spaceranger/outs \
  --hd-resolution 8um \
  --output /path/to/output \
  --sample-name sample_hd_8um \
  --n-clusters 8 \
  --markers
```

Additional starting formats and exact commands are documented in the numbered `protocol.md` file in each example directory.

## Run an atlas

After downloading and extracting the selected example release asset:

```bash
export SPATIALATLAS_CONFIG="$PWD/examples/01_public_visium/config.local.yaml"
python software/app.py
```

Open `http://127.0.0.1:8050` and stop the server with `Ctrl-C`.

For Docker deployment:

```bash
cd software
docker compose config --quiet
docker compose up --build
```

Production paths, HTTPS, access control, monitoring, and updates are covered in [docs/deployment.md](docs/deployment.md).

## Examples and data assets

| Release asset | Validated content |
|---|---|
| `SpatialAtlas-v1.0.0-example01-public-visium.tar.gz` | 4,037 post-QC lymph-node spots and 18,500 genes |
| `SpatialAtlas-v1.0.0-example02-visium-hd-8um.tar.gz` | 143,112 adrenal 8 µm profiles and 15,109 genes |
| `SpatialAtlas-v1.0.0-example03-multisample-adrenal.tar.gz` | female and male adrenal atlas with 110,823 profiles |
| `SpatialAtlas-v1.0.0-example04-harmony-integration.tar.gz` | compact 2,000-profile integration example |
| `SpatialAtlas-v1.0.0-example05-platform-benchmark.tar` | matched Shiny data and retained TissUUmaps benchmark objects |

Verify downloaded assets from the release directory:

```bash
shasum -a 256 -c release-assets.sha256
```

After extracting the archives into the repository root, verify the individual retained data files with:

```bash
shasum -a 256 -c example-data.sha256
```

`checksums.sha256` covers the compact files stored directly in the repository.

Data provenance and reuse conditions are listed in [DATASETS.md](DATASETS.md). The classic adrenal object is associated with DOI `10.5603/fhc.108988` and GEO `GSE283302`; the Visium HD adrenal data are associated with DOI `10.1038/s42003-026-10697-9` and GEO `GSE312015`.

## Automated tests

```bash
pip install -r software/requirements.txt -r requirements-dev.txt
pytest --cov --cov-report=term-missing
```

The suite in `tests/` builds every fixture synthetically in a temporary
directory, runs in roughly fifteen seconds, and needs no biological data. It
runs on each push through `.github/workflows/tests.yml`, which publishes the
coverage and JUnit reports as build artifacts and fails below 70 % line
coverage. See [TESTING.md](TESTING.md) for what is covered at each level.

## Validation and benchmarks

Artifact and viewer checks are in `validation/`, with recorded JSON results in `validation/results/`. The complete matched Shiny, independent TissUUmaps, browser-readiness, and concurrent-session workflow is in `examples/05_adrenal_platform_benchmark/protocol.md`.

Recorded timings describe the tested Apple M4 Pro host and pinned container configurations; they are not hardware-independent performance guarantees. Response-generation measurements also preserve each application's actual rendering path: Shiny generates PNG output on the server, whereas SpatialAtlas returns Plotly JSON for client-side WebGL rendering.

## Scope and data access

SpatialAtlas is a dissemination layer, not a replacement for spatial statistics, trajectory inference, deconvolution, or cell-cell communication analysis. Do not expose restricted biological data without confirming consent, licensing, institutional policy, and suitable authentication or network controls.

Software is released under the MIT License. Biological data retain the conditions of their originating studies and repositories. Citation metadata are provided in `CITATION.cff`.

## Citation

The evaluated SpatialAtlas v1.0.0 source release is archived in Zenodo: https://doi.org/10.5281/zenodo.21476516. The five versioned example and benchmark data archives are distributed with GitHub Release v1.0.0 and verified by `release-assets.sha256`.
