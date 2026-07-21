# Example 3 — female and male mouse adrenal atlas

This example opens two processed adrenal sections to verify the atlas data contract and reproduce expected biological landmarks; it does not repeat the upstream differential-expression analysis.

## 1. Inspect the AnnData objects

```bash
python - <<'PY'
import anndata as ad
from pathlib import Path

for path in sorted(Path('examples/03_multisample_adrenal/input').glob('*/*.h5ad')):
    atlas = ad.read_h5ad(path, backed='r')
    print(path.name, atlas.shape,
          'cluster' in atlas.obs,
          'spatial' in atlas.obsm,
          'UMAP' in atlas.obsm)
PY
```

The female object contains 49,447 profiles and 15,009 genes; the male object contains 61,376 profiles and 17,141 genes. Together they represent 110,823 profiles and 17,173 distinct genes. Both use the same eight cluster labels and share a joint UMAP coordinate system.

## 2. Exercise the viewer

```bash
python validation/validate_viewer.py \
  --app software/app.py \
  --config examples/03_multisample_adrenal/config.local.yaml \
  --output validation/results/multisample_adrenal_viewer_validation.json
```

On the reference machine the two samples were imported in 2.66 s and their sparse expression matrices occupied 0.624 GiB. The script also requests the atlas page and both tissue-image routes.

## 3. Reproduce the views used in Figure 3

```bash
export SPATIALATLAS_CONFIG="$PWD/examples/03_multisample_adrenal/config.local.yaml"
python software/app.py
```

Open `http://127.0.0.1:8052` and inspect four anatomical markers:

1. `Cyp11b2`, enriched in the zona glomerulosa;
2. `Gpc3`, prominent in capsular and stromal regions;
3. `Npy`, marking medullary expression;
4. `Akr1c18`, showing the female adrenal X-zone.

Use the same genes in the four-channel view (green, red, cyan, and yellow). A gene that is absent from one sample should appear as zero in that sample without interrupting the other panels. These visual checks establish orientation and expected localization; they are not a substitute for statistical inference.
