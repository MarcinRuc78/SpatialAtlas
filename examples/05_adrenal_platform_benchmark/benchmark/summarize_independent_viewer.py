#!/usr/bin/env python3
"""Create the independent-viewer and concurrency summary table."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path


BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
PACKAGE = ROOT.parents[1]


def read_json(name: str) -> dict:
    return json.loads((BENCH / name).read_text(encoding="utf-8"))


startup = read_json("startup_memory_results.json")["summary"]["spatialatlas_visium_hd"]
tissuumaps = read_json("tissuumaps_startup_memory_results.json")["summary"]
browser = read_json("browser_ready_results.json")
concurrent = read_json("concurrent_users_results.json")["summary"]
display = read_json("display_image_validation.json")

spatial_ready = statistics.median(browser["spatialatlas"]["ready_ms"]) / 1000
tissuumaps_ready = statistics.median(
    browser["tissuumaps_prepared_csc"]["ready_ms"]
) / 1000
tmap = ROOT / "prepared" / "tissuumaps" / "adrenal_hd_tissuumaps_tmap.h5ad"
deployed_pyramid = (
    ROOT / "prepared" / "tissuumaps"
    / "adrenal_hd_tissuumaps_tmap.h5ad_files" / "adrenal_hd_8um"
    / "img" / "tissue.tif"
)
source_h5ad = (
    PACKAGE / "examples" / "02_visium_hd_8um" / "output"
    / "sample_adrenal_female_hd_8um.h5ad"
)

rows = [
    {
        "metric": "warm_server_startup_median",
        "spatialatlas": startup["warm_start_median_s"],
        "tissuumaps": tissuumaps["warm_start_median_s"],
        "unit": "s",
        "n_per_viewer_or_level": 5,
        "note": "SpatialAtlas root HTTP 200; prepared TissUUmaps HD viewer HTTP 200",
    },
    {
        "metric": "steady_container_memory_median",
        "spatialatlas": startup["warm_memory_median_bytes"] / 2**20,
        "tissuumaps": tissuumaps["warm_memory_median_bytes"] / 2**20,
        "unit": "MiB",
        "n_per_viewer_or_level": 5,
        "note": "cgroup memory.current two seconds after readiness",
    },
    {
        "metric": "browser_visual_readiness_median",
        "spatialatlas": spatial_ready,
        "tissuumaps": tissuumaps_ready,
        "unit": "s",
        "n_per_viewer_or_level": 5,
        "note": "rendered 143112-profile WebGL map versus prepared OpenSeadragon viewer and gene selector",
    },
    {
        "metric": "prepared_primary_object",
        "spatialatlas": source_h5ad.stat().st_size / 2**20,
        "tissuumaps": tmap.stat().st_size / 2**20,
        "unit": "MiB",
        "n_per_viewer_or_level": 1,
        "note": f"deployed TissUUmaps viewer also uses an {deployed_pyramid.stat().st_size / 2**20:.1f} MiB pyramid; both URL-specific copies are retained",
    },
    {
        "metric": "display_image_transfer",
        "spatialatlas": display["display_bytes"] / 2**20,
        "tissuumaps": "",
        "unit": "MiB",
        "n_per_viewer_or_level": 1,
        "note": f"full-dimension WebP, PSNR {display['psnr_db']:.2f} dB; archival PNG retained",
    },
    {
        "metric": "first_unprepared_request",
        "spatialatlas": "",
        "tissuumaps": browser["tissuumaps_unprepared_csr"]["first_server_response_s"],
        "unit": "s",
        "n_per_viewer_or_level": 1,
        "note": "CSR-to-CSC conversion and image-pyramid creation",
    },
]

for users in (1, 5, 10, 20):
    result = concurrent[str(users)]
    rows.append({
        "metric": f"concurrent_{users}_user_p95" if users == 1 else f"concurrent_{users}_users_p95",
        "spatialatlas": result["session_p95_s"],
        "tissuumaps": "",
        "unit": "s",
        "n_per_viewer_or_level": result["sessions"],
        "note": f"{result['bytes_per_successful_session'] / 2**20:.2f} MiB uncached first-visit bundle; {result['errors']} errors",
    })

with (BENCH / "independent_browser_summary.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print(json.dumps(rows, indent=2))
