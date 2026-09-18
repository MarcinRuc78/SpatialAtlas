# Testing and quality assurance

SpatialAtlas is verified at two levels.

| Level | Location | Runtime | Data |
| --- | --- | --- | --- |
| Automated unit and functional tests | `tests/` | ~15 s | Synthetic fixtures generated in a temporary directory |
| Release-level integration and benchmark workflows | `validation/`, `examples/` | Minutes to hours | Real Visium, Visium HD and multi-sample atlases |

The fast suite runs on every push and pull request through
`.github/workflows/tests.yml`. The release-level workflows remain the evidence
for the published measurements and are executed before a release, as recorded
in Supplementary Material S3 and S4.

## Running the fast suite

```bash
pip install -r software/requirements.txt
pip install -r requirements-dev.txt
pytest
```

With coverage:

```bash
pytest --cov --cov-report=term-missing --cov-report=html
```

The HTML report is written to `htmlcov/index.html`. Continuous integration
publishes `coverage.xml`, the HTML report and the JUnit results as build
artifacts, prints a per-file coverage table in the job summary, and fails the
build below 70 % line coverage.

## What the fast suite covers

No biological data is required: every fixture is built in a temporary
directory by `tests/conftest.py` and `tests/test_spaceranger_fixture.py`.

- **Gene availability across samples** (`test_gene_availability.py`) — a gene
  that is absent from a sample matrix is reported as unavailable and is never
  substituted by a measured expression of zero. The synthetic atlas contains
  genes deliberately present in only one of its two samples.
- **Visualization callbacks** (`test_callbacks.py`) — spatial, UMAP, violin,
  marker-heatmap and multi-gene views are executed and serialized, including
  the unavailable-gene branch of each view.
- **Data contract and HTTP routes** (`test_data_contract.py`) — configuration
  path resolution, required H5AD fields, the 1–4 sample limit, sparse
  column-oriented storage, cache headers on the tissue-image routes, and the
  gene search behaviour.
- **Input handling** (`test_inputs_and_images.py`) — file discovery, the
  optional display image and its aspect-ratio check, and coordinate alignment.
- **Preparation workflow** (`test_preparation.py`) — QC thresholds including
  removal of genes detected in fewer than ten observations, log-normalization,
  marker calling, and both alignment paths.
- **Space Ranger loaders** (`test_spaceranger_fixture.py`) — a miniature but
  structurally faithful Space Ranger run, covering the 10x HDF5 matrix, both
  the v1.x and v2.x tissue-position layouts, in-tissue filtering, and scale
  factors.
- **Multi-sample integration** (`test_integration_workflow.py`) — Harmony
  integration writes to a new directory, leaves inputs byte-identical, assigns
  shared clusters and joint UMAP coordinates, and preserves each sample's own
  gene set.

## Current coverage

Measured with the pinned Python 3.11 environment:

| File | Line coverage |
| --- | ---: |
| `software/app.py` | 88 % |
| `software/integrate_samples.py` | 78 % |
| `software/prepare_data.py` | 60 % |
| **Total** | **76 %** |

The uncovered remainder of `prepare_data.py` is dominated by the Visium HD
binned and cell-segmentation loaders and the standalone 10x/GEO matrix path,
which require real Space Ranger outputs; these are exercised by the
release-level protocols in `examples/02_visium_hd_8um` and
`examples/03_multisample_adrenal`.
