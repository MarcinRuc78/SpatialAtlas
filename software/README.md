# SpatialAtlas software

SpatialAtlas is a self-hosted Dash/Plotly viewer for one to four processed Visium or Visium HD samples. It reads sparse AnnData matrices, materializes one requested gene at a time, and serves tissue images through a cacheable endpoint.

## Minimum input contract

Each configured sample needs:

- an `.h5ad` file with expression in `.X`, `obs['cluster']`, and an `N x 2` `obsm['spatial']` matrix;
- an optional `N x 2` `obsm['UMAP']` matrix;
- a browser-readable tissue image and, optionally, a full-dimension `display_image` in a more compact encoding;
- optional markers with at least `gene` and `cluster` columns.

All paths in a YAML file are resolved relative to that YAML file. The example configurations are therefore portable.

## Clean local installation

Run these commands from the package root:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r software/requirements.txt
```

Prepare standard Visium input:

```bash
python software/prepare_data.py \
  --spaceranger /path/to/spaceranger/outs \
  --output /path/to/output \
  --sample-name sample_1 \
  --n-clusters 8 \
  --markers
```

Prepare Visium HD 8 µm input:

```bash
python software/prepare_data.py \
  --spaceranger-hd /path/to/spaceranger/outs \
  --hd-resolution 8um \
  --output /path/to/output \
  --sample-name sample_hd_8um \
  --n-clusters 8 \
  --markers
```

Launch a validated local example:

```bash
export SPATIALATLAS_CONFIG="$PWD/examples/01_public_visium/config.local.yaml"
python software/app.py
```

Open `http://127.0.0.1:8050`. Stop the server with `Ctrl-C`.

## Multi-sample integration

The default is non-destructive: integrated copies are written to the selected output directory.

```bash
python software/integrate_samples.py \
  --samples sample_a.h5ad sample_b.h5ad \
  --n-clusters 8 \
  --markers \
  --output-dir integrated
```

Use `--in-place` only when intentional replacement of the source files is desired.

## Docker

The supplied Compose file defaults to the public lymph-node example:

```bash
cd software
docker compose config --quiet
docker compose up --build
```

For a different dataset, set `SPATIALATLAS_DATA_DIR`, `SPATIALATLAS_CONFIG_FILE`, and optionally `SPATIALATLAS_MEMORY_LIMIT`. Public HTTPS additionally requires replacing `YOUR_DOMAIN` in `Caddyfile`, valid DNS, firewall/ingress configuration, and a review of data-access policy.

## Detailed protocols

See the package [README](../README.md), the per-example `protocol.md` files, and [deployment guide](../docs/deployment.md). Exact validation reports are stored in `validation/results/`.
