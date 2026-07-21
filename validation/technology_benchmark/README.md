# Technology benchmark

This directory retains the general three-example SpatialAtlas measurements. The matched adrenal application benchmark and its raw data are in `examples/05_adrenal_platform_benchmark`.

## General SpatialAtlas workloads

`benchmark_three_examples.py` measures container startup and steady memory for the public Visium, Visium HD, and multi-sample adrenal examples. `benchmark_dash_core.py` measures expression-vector extraction, Plotly figure construction, JSON serialization, and payload size for four prespecified genes.

Run from the package root:

```bash
python validation/technology_benchmark/benchmark_three_examples.py
```

The committed JSON files are direct outputs of those runs.

## Adrenal platform and resolution comparison

The controlled comparison uses the 992-spot classic Visium adrenal object and Seurat/Shiny application from Blatkiewicz et al. (doi:10.5603/fhc.108988; GEO GSE283302), together with a matched SpatialAtlas H5AD. The same 143,112-bin Visium HD adrenal object is evaluated in a purpose-built R/Shiny reference application and SpatialAtlas.

All four applications ran as `linux/amd64` containers on the same Apple M4 Pro host. Six starts were recorded per scenario; the first was treated as cold and the following five were summarized. `Npy`, `Cyp11b1`, `Cyp11b2`, and `Th` were each measured eight times after one warm-up.

The complete scripts, raw CSV/JSON measurements, source objects, compatibility notes, digests, and commands are in:

```text
examples/05_adrenal_platform_benchmark/
```

The aggregate result is copied here as `adrenal_platform_summary.json` for machine-readable validation.

## Interpretation constraints

- The classic Seurat/Shiny and SpatialAtlas rows are matched on the same object.
- The Visium HD Shiny and SpatialAtlas rows use the same expression values, coordinates, histology and markers.
- The HD Shiny application is a fully supplied reference implementation constructed from the validated H5AD object.
- The Shiny path returns a server-rendered PNG; SpatialAtlas returns Plotly JSON for client-side WebGL.
- Container memory includes the complete application process and loaded data.
- Network transfer and browser rendering are outside the server-response measurements.
