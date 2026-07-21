# Example 2 — mouse adrenal gland (Visium HD, 8 µm bins)

This example prepares a female mouse adrenal section from the Visium HD study by Blatkiewicz et al. (doi:10.1038/s42003-026-10697-9), using the minimal Space Ranger 8 µm output required by SpatialAtlas. The complete run retained 143,112 profiles and reached 7.91 GB resident memory, so 12 GB of available RAM is a practical minimum. Data reuse remains subject to the originating adrenal study described in `DATASETS.md`.

## 1. Check the Space Ranger output

```bash
HD=examples/02_visium_hd_8um/input/binned_outputs/square_008um
test -s "$HD/filtered_feature_bc_matrix.h5"
test -s "$HD/spatial/tissue_positions.parquet"
test -s "$HD/spatial/scalefactors_json.json"
test -s "$HD/spatial/tissue_hires_image.png"
```

## 2. Prepare the 8 µm atlas

```bash
python software/prepare_data.py \
  --spaceranger-hd examples/02_visium_hd_8um/input \
  --hd-resolution 8um \
  --output /tmp/spatialatlas_hd_8um \
  --sample-name adrenal_female_hd_8um \
  --n-clusters 8 \
  --markers
```

The reference run matched 250,669 spatial barcodes. Default filtering retained 143,112 profiles and 15,109 genes; the eight-cluster solution was obtained near Leiden resolution 0.5706. Preparation took 144.5 s on the reference machine.

## 3. Check the generated artifact

```bash
python validation/validate_artifacts.py \
  --h5ad /tmp/spatialatlas_hd_8um/sample_adrenal_female_hd_8um.h5ad \
  --image /tmp/spatialatlas_hd_8um/tissue_adrenal_female_hd_8um.png \
  --markers /tmp/spatialatlas_hd_8um/markers.csv \
  --output-dir /tmp/spatialatlas_hd_validation \
  --prefix visium_hd_8um
```

The JSON report should contain the dimensions above, finite spatial and UMAP coordinates, eight cluster labels, and non-zero marker overlap.

## 4. Inspect the dense spatial view

```bash
export SPATIALATLAS_CONFIG="$PWD/examples/02_visium_hd_8um/config.local.yaml"
python software/app.py
```

Open `http://127.0.0.1:8051`. Use a point size close to 0.5 and inspect `Cyp11b2`, `Gpc3`, `Npy`, and `Akr1c18`. Switching genes should update the point layer without reloading the tissue image.

For a repeatable server-side check, run:

```bash
python validation/validate_viewer.py \
  --app software/app.py \
  --config examples/02_visium_hd_8um/config.local.yaml \
  --output validation/results/visium_hd_8um_viewer_validation.json
```

The final recorded import time was 2.16 s and the sparse expression matrix occupied 0.291 GiB. Callback timings in the report exclude network transfer and browser rendering.
