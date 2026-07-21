# Example 4 — non-destructive multi-sample integration

Two compact AnnData objects are included to exercise the Harmony integration workflow without rerunning a large biological study. The example also demonstrates that source files are left unchanged unless `--in-place` is supplied explicitly.

## 1. Record the source checksums

```bash
shasum -a 256 examples/04_multisample_integration/input/*.h5ad \
  > /tmp/spatialatlas_integration_before.sha256
```

## 2. Integrate the samples

```bash
python software/integrate_samples.py \
  --samples \
    examples/04_multisample_integration/input/sample_replicate_a.h5ad \
    examples/04_multisample_integration/input/sample_replicate_b.h5ad \
  --n-clusters 6 \
  --markers \
  --output-dir /tmp/spatialatlas_integrated
```

The combined input contains 2,000 profiles and 2,000 common genes. The expected output is two H5AD files with shared cluster and UMAP coordinates plus `markers.csv`. The final reference run used Harmony, completed in 11.7 s, and reached 0.98 GB resident memory.

## 3. Confirm that the inputs were preserved

```bash
shasum -a 256 -c /tmp/spatialatlas_integration_before.sha256
test -s /tmp/spatialatlas_integrated/sample_replicate_a.h5ad
test -s /tmp/spatialatlas_integrated/sample_replicate_b.h5ad
test -s /tmp/spatialatlas_integrated/markers.csv
```

Both checksum lines should report `OK`. Use `--in-place` only when replacing the input objects is intentional; it cannot be combined with `--output-dir`.

## 4. Inspect the integrated fields

```bash
python - <<'PY'
import anndata as ad
from pathlib import Path

for path in sorted(Path('/tmp/spatialatlas_integrated').glob('*.h5ad')):
    atlas = ad.read_h5ad(path)
    assert atlas.obs['cluster'].nunique() == 6
    assert atlas.obsm['UMAP'].shape == (atlas.n_obs, 2)
    print(path.name, atlas.shape)
PY
```
