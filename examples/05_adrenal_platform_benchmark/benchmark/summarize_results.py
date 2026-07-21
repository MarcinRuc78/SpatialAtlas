#!/usr/bin/env python3
"""Create aggregate tables from raw startup, memory and response measurements."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmark"
PACKAGE = ROOT.parents[1]

LABELS = {
    "seurat_shiny_visium": "Seurat/Shiny — classic Visium adrenal",
    "spatialatlas_visium": "SpatialAtlas — classic Visium adrenal",
    "seurat_shiny_visium_hd": "Shiny — Visium HD adrenal, 8 µm",
    "spatialatlas_visium_hd": "SpatialAtlas — Visium HD adrenal, 8 µm",
}
PROFILES = {
    "seurat_shiny_visium": 992,
    "spatialatlas_visium": 992,
    "seurat_shiny_visium_hd": 143_112,
    "spatialatlas_visium_hd": 143_112,
}
GENES = {
    "seurat_shiny_visium": 19_465,
    "spatialatlas_visium": 19_465,
    "seurat_shiny_visium_hd": 15_109,
    "spatialatlas_visium_hd": 15_109,
}
OBJECTS = {
    "seurat_shiny_visium": ROOT / "source_data" / "adrenal_visium_seurat.rds",
    "spatialatlas_visium": ROOT / "prepared" / "sample_adrenal_visium_classic.h5ad",
    "seurat_shiny_visium_hd": ROOT / "source_data" / "adrenal_visium_hd_shiny.rds",
    "spatialatlas_visium_hd": PACKAGE / "examples" / "02_visium_hd_8um" / "output" / "sample_adrenal_female_hd_8um.h5ad",
}
CORE_FILES = {
    "seurat_shiny_visium": BENCH / "shiny_core_results.csv",
    "spatialatlas_visium": BENCH / "spatialatlas_classic_core_results.csv",
    "seurat_shiny_visium_hd": BENCH / "shiny_hd_core_results.csv",
    "spatialatlas_visium_hd": BENCH / "spatialatlas_hd_core_results.csv",
}


startup = json.loads((BENCH / "startup_memory_results.json").read_text(encoding="utf-8"))
rows = []
for key in LABELS:
    with CORE_FILES[key].open(newline="", encoding="utf-8") as handle:
        core = list(csv.DictReader(handle))
    summary = startup["summary"][key]
    rows.append({
        "scenario": key,
        "label": LABELS[key],
        "profiles": PROFILES[key],
        "genes": GENES[key],
        "object_mib": OBJECTS[key].stat().st_size / 2**20,
        "warm_start_median_s": summary["warm_start_median_s"],
        "warm_start_min_s": summary["warm_start_min_s"],
        "warm_start_max_s": summary["warm_start_max_s"],
        "steady_memory_median_mib": summary["warm_memory_median_bytes"] / 2**20,
        "steady_memory_min_mib": summary["warm_memory_min_bytes"] / 2**20,
        "steady_memory_max_mib": summary["warm_memory_max_bytes"] / 2**20,
        "gene_response_median_ms": statistics.median(float(row["total_server_s"]) for row in core) * 1000,
        "gene_response_min_ms": min(float(row["total_server_s"]) for row in core) * 1000,
        "gene_response_max_ms": max(float(row["total_server_s"]) for row in core) * 1000,
        "payload_median_kib": statistics.median(float(row["response_bytes"]) for row in core) / 1024,
        "payload_min_kib": min(float(row["response_bytes"]) for row in core) / 1024,
        "payload_max_kib": max(float(row["response_bytes"]) for row in core) / 1024,
    })

fields = list(rows[0])
with (BENCH / "benchmark_summary.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

by_key = {row["scenario"]: row for row in rows}
shiny = by_key["seurat_shiny_visium"]
classic = by_key["spatialatlas_visium"]
hd = by_key["spatialatlas_visium_hd"]
shiny_hd = by_key["seurat_shiny_visium_hd"]
ratios = {
    "classic_spatialatlas_vs_shiny": {
        "startup_fold_faster": shiny["warm_start_median_s"] / classic["warm_start_median_s"],
        "memory_fold_lower": shiny["steady_memory_median_mib"] / classic["steady_memory_median_mib"],
        "gene_response_fold_faster": shiny["gene_response_median_ms"] / classic["gene_response_median_ms"],
        "payload_fold_smaller": shiny["payload_median_kib"] / classic["payload_median_kib"],
    },
    "hd_spatialatlas_vs_shiny": {
        "startup_fold_faster": shiny_hd["warm_start_median_s"] / hd["warm_start_median_s"],
        "spatialatlas_memory_vs_shiny": hd["steady_memory_median_mib"] / shiny_hd["steady_memory_median_mib"],
        "gene_response_fold_faster": shiny_hd["gene_response_median_ms"] / hd["gene_response_median_ms"],
        "payload_fold_smaller": shiny_hd["payload_median_kib"] / hd["payload_median_kib"],
    },
    "visium_hd_vs_classic_spatialatlas": {
        "profile_count_fold": hd["profiles"] / classic["profiles"],
        "startup_fold": hd["warm_start_median_s"] / classic["warm_start_median_s"],
        "memory_fold": hd["steady_memory_median_mib"] / classic["steady_memory_median_mib"],
        "gene_response_fold": hd["gene_response_median_ms"] / classic["gene_response_median_ms"],
        "payload_fold": hd["payload_median_kib"] / classic["payload_median_kib"],
    },
}
payload = {
    "protocol": startup["protocol"],
    "summary": rows,
    "ratios": ratios,
    "interpretation": (
        "The classic Visium pair and the Visium HD pair are matched implementation comparisons. "
        "Within each pair, Shiny and SpatialAtlas use the same expression values, coordinates, histology and markers."
    ),
}
(BENCH / "benchmark_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
