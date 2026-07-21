#!/usr/bin/env python3
"""Measure expression lookup, figure creation and payload serialization."""

from __future__ import annotations

import csv
import os
import statistics
import sys
import time

sys.path.insert(0, "/app")
import app as spatialatlas  # noqa: E402


GENES = ("Npy", "Cyp11b1", "Cyp11b2", "Th")
REPETITIONS = 8
LABEL = os.environ["BENCHMARK_LABEL"]
OUTPUT = os.environ["BENCHMARK_OUTPUT"]


def measure(gene: str, iteration: int) -> dict[str, object]:
    start = time.perf_counter()
    spatialatlas.get_expr(gene)
    extraction_s = time.perf_counter() - start

    start = time.perf_counter()
    figure = spatialatlas.update_spatial(gene, 3.0, 0.7, "Transparent-Green")
    figure_build_s = time.perf_counter() - start

    start = time.perf_counter()
    payload = figure.to_json().encode("utf-8")
    response_render_s = time.perf_counter() - start
    return {
        "implementation": LABEL,
        "gene": gene,
        "iteration": iteration,
        "extraction_s": extraction_s,
        "figure_build_s": figure_build_s,
        "response_render_s": response_render_s,
        "total_server_s": extraction_s + figure_build_s + response_render_s,
        "response_bytes": len(payload),
    }


for gene in GENES:
    measure(gene, 0)

rows = [
    measure(gene, iteration)
    for iteration in range(1, REPETITIONS + 1)
    for gene in GENES
]

with open(OUTPUT, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

for field in ("extraction_s", "figure_build_s", "response_render_s", "total_server_s", "response_bytes"):
    value = statistics.median(float(row[field]) for row in rows)
    print(f"{field}: {value:.9f}")
