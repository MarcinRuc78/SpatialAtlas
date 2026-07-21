#!/usr/bin/env python3
"""Write sizes and SHA-256 digests for the retained benchmark data files."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).resolve().parent / "benchmark_data_inventory.csv"

FILES = {
    "classic Seurat reference object": "examples/05_adrenal_platform_benchmark/source_data/adrenal_visium_seurat.rds",
    "classic matched H5AD": "examples/05_adrenal_platform_benchmark/prepared/sample_adrenal_visium_classic.h5ad",
    "classic histology": "examples/05_adrenal_platform_benchmark/prepared/tissue_adrenal_visium_classic.png",
    "HD Shiny reference object": "examples/05_adrenal_platform_benchmark/source_data/adrenal_visium_hd_shiny.rds",
    "HD source H5AD": "examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad",
    "HD archival histology": "examples/02_visium_hd_8um/output/tissue_adrenal_female_hd_8um.png",
    "HD SpatialAtlas display image": "examples/05_adrenal_platform_benchmark/prepared/spatialatlas/tissue_adrenal_female_hd_8um_display.webp",
    "TissUUmaps adapter H5AD": "examples/05_adrenal_platform_benchmark/prepared/tissuumaps/adrenal_hd_tissuumaps.h5ad",
    "TissUUmaps prepared CSC H5AD": "examples/05_adrenal_platform_benchmark/prepared/tissuumaps/adrenal_hd_tissuumaps_tmap.h5ad",
    "TissUUmaps adapter image pyramid": "examples/05_adrenal_platform_benchmark/prepared/tissuumaps/adrenal_hd_tissuumaps.h5ad_files/adrenal_hd_8um/img/tissue.tif",
    "TissUUmaps deployed image pyramid": "examples/05_adrenal_platform_benchmark/prepared/tissuumaps/adrenal_hd_tissuumaps_tmap.h5ad_files/adrenal_hd_8um/img/tissue.tif",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


rows = []
for role, relative in FILES.items():
    path = PACKAGE / relative
    if not path.is_file():
        raise FileNotFoundError(path)
    rows.append({
        "role": role,
        "path": relative,
        "bytes": path.stat().st_size,
        "mib": f"{path.stat().st_size / 2**20:.3f}",
        "sha256": sha256(path),
    })

with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {OUTPUT} with {len(rows)} files")
