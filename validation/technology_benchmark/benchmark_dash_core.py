#!/usr/bin/env python3
"""Measure the server-side gene-view path for one configured atlas."""

import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, "/app")
import app as atlas

genes = [gene for gene in os.environ["BENCH_GENES"].split(",") if gene]
repetitions = int(os.environ.get("BENCH_REPETITIONS", "8"))
output = Path(os.environ["BENCH_OUTPUT"])


def cgroup_value(name):
    path = Path("/sys/fs/cgroup") / name
    try:
        value = path.read_text().strip()
        return None if value == "max" else int(value)
    except (FileNotFoundError, ValueError):
        return None


def measure(gene, iteration):
    start = time.perf_counter()
    atlas.get_expr(gene)
    extraction = time.perf_counter() - start

    start = time.perf_counter()
    point_size = atlas.display_cfg.get("default_point_size", atlas._default_pt_size)
    figure = atlas.update_spatial(gene, point_size, 0.7, "Transparent-Green")
    figure_build = time.perf_counter() - start

    start = time.perf_counter()
    payload = figure.to_json().encode("utf-8")
    serialization = time.perf_counter() - start
    return {
        "gene": gene,
        "iteration": iteration,
        "extraction_s": extraction,
        "figure_build_s": figure_build,
        "serialization_s": serialization,
        "server_total_s": figure_build + serialization,
        "response_bytes": len(payload),
        "memory_after_bytes": cgroup_value("memory.current"),
    }


for gene in genes:
    measure(gene, 0)

rows = [measure(gene, iteration)
        for iteration in range(1, repetitions + 1)
        for gene in genes]
numeric = ("extraction_s", "figure_build_s", "serialization_s", "server_total_s", "response_bytes")
summary = {name: statistics.median(row[name] for row in rows) for name in numeric}
summary["memory_peak_bytes"] = cgroup_value("memory.peak")
summary["n_measurements"] = len(rows)
output.write_text(json.dumps({"genes": genes, "runs": rows, "summary": summary}, indent=2) + "\n")
print(json.dumps(summary, indent=2))
