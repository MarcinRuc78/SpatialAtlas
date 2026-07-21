# Example 1 — human lymph node (standard Visium)

This example starts from the 10x Genomics Human Lymph Node dataset processed with Space Ranger 1.0.0 and produces the files used by the public SpatialAtlas demonstration. The source dataset is distributed by 10x Genomics under CC BY 4.0; the files used here were downloaded on 20 July 2026.

Run all commands from the package root in a Python 3.11 environment containing `software/requirements.txt`. Allow about 4 GB of memory and 1 GB of free disk space.

## 1. Check the source files

```bash
test -s examples/01_public_visium/input/filtered_feature_bc_matrix.h5
test -s examples/01_public_visium/input/spatial/tissue_positions_list.csv
test -s examples/01_public_visium/input/spatial/scalefactors_json.json
test -s examples/01_public_visium/input/spatial/tissue_hires_image.png
```

The `input/` directory is ready to use. To reconstruct it from the downloaded files instead:

```bash
mkdir -p /tmp/spatialatlas_public/input
cp examples/01_public_visium/downloads/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5 \
  /tmp/spatialatlas_public/input/filtered_feature_bc_matrix.h5
tar -xzf examples/01_public_visium/downloads/V1_Human_Lymph_Node_spatial.tar.gz \
  -C /tmp/spatialatlas_public/input
```

## 2. Prepare the atlas

```bash
python software/prepare_data.py \
  --h5-matrix examples/01_public_visium/input/filtered_feature_bc_matrix.h5 \
  --spatial-dir examples/01_public_visium/input/spatial \
  --output /tmp/spatialatlas_public/output \
  --sample-name human_lymph_node \
  --n-clusters 8 \
  --markers
```

The bundled run matched 4,039 tissue barcodes. After the default filters, the output contained 4,037 spots, 18,500 genes, eight clusters, and 200 marker-table rows. The first run on the reference machine took 69.2 s and reached 2.46 GB resident memory; these values are useful for planning but will vary by computer.

## 3. Inspect the output files

```bash
python validation/validate_artifacts.py \
  --h5ad /tmp/spatialatlas_public/output/sample_human_lymph_node.h5ad \
  --image /tmp/spatialatlas_public/output/tissue_human_lymph_node.png \
  --markers /tmp/spatialatlas_public/output/markers.csv \
  --output-dir /tmp/spatialatlas_public/validation \
  --prefix public_visium
```

The command checks identifiers, coordinate arrays, cluster annotations, UMAP values, image decoding, and marker overlap. Its JSON output should report the dimensions above and finish with `"status": "PASS"`.

## 4. Open the atlas

```bash
export SPATIALATLAS_CONFIG="$PWD/examples/01_public_visium/config.local.yaml"
python software/app.py
```

Open `http://127.0.0.1:8050`. Search for `CD74`, display it in the Spatial and UMAP panels, and change point size and opacity to confirm that the coordinates remain aligned with the lymph-node image. The Violin, Marker heatmap, Multi-gene, and Citation tabs should also render.

The same views can be exercised without a browser:

```bash
python validation/validate_viewer.py \
  --app software/app.py \
  --config examples/01_public_visium/config.local.yaml \
  --output validation/results/public_visium_viewer_validation.json
```

## 5. Optional container run

```bash
cd software
docker build -t spatialatlas-public:test .
docker run --rm -p 8050:8050 \
  -v "$PWD/../examples/01_public_visium/output:/data:ro" \
  -v "$PWD/../examples/01_public_visium/config.docker.yaml:/app/config.yaml:ro" \
  spatialatlas-public:test
```

From a second terminal, request `http://127.0.0.1:8050/`. Stop the container with `Ctrl-C` when the page and tissue image have loaded.
