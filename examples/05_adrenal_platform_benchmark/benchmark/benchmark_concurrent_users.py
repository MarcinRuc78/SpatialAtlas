#!/usr/bin/env python3
"""Measure concurrent first-load request bundles for the Visium HD atlas.

Each virtual user requests the HTML shell, Dash layout, callback dependency
graph and the configured display histology image. urllib has no browser cache,
so this is deliberately a conservative first-visit network test. Browser
parsing and WebGL readiness are measured separately with the recorded browser
protocol.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import statistics
import time
import urllib.request
from pathlib import Path


PATHS = ("/", "/_dash-layout", "/_dash-dependencies", "/_spatialatlas/tissue/0")


def percentile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def one_user(base_url: str, user_id: int) -> dict[str, object]:
    started = time.perf_counter()
    total_bytes = 0
    statuses: list[int] = []
    error = None
    try:
        for path in PATHS:
            request = urllib.request.Request(
                base_url.rstrip("/") + path,
                headers={"User-Agent": f"SpatialAtlas-benchmark/{user_id}", "Cache-Control": "no-cache"},
            )
            with urllib.request.urlopen(request, timeout=180) as response:
                body = response.read()
                statuses.append(response.status)
                total_bytes += len(body)
    except Exception as exc:  # retained in the machine-readable report
        error = repr(exc)
    return {
        "user_id": user_id,
        "elapsed_s": time.perf_counter() - started,
        "bytes": total_bytes,
        "statuses": statuses,
        "error": error,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:18052")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rounds", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []
    summary: dict[str, object] = {}
    for concurrency in (1, 5, 10, 20):
        level_rows: list[dict[str, object]] = []
        batch_durations: list[float] = []
        for round_index in range(args.rounds):
            batch_start = time.perf_counter()
            with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as pool:
                futures = [
                    pool.submit(one_user, args.base_url, round_index * 1000 + index)
                    for index in range(concurrency)
                ]
                batch = [future.result() for future in futures]
            batch_elapsed = time.perf_counter() - batch_start
            batch_durations.append(batch_elapsed)
            for row in batch:
                row.update({"concurrent_users": concurrency, "round": round_index})
                rows.append(row)
                level_rows.append(row)
        times = [float(row["elapsed_s"]) for row in level_rows]
        errors = sum(row["error"] is not None for row in level_rows)
        transferred = sum(int(row["bytes"]) for row in level_rows)
        summary[str(concurrency)] = {
            "sessions": len(level_rows),
            "errors": errors,
            "session_median_s": statistics.median(times),
            "session_p95_s": percentile(times, 0.95),
            "session_min_s": min(times),
            "session_max_s": max(times),
            "median_batch_s": statistics.median(batch_durations),
            "aggregate_sessions_per_s": len(level_rows) / sum(batch_durations),
            "aggregate_mib_per_s": transferred / sum(batch_durations) / 1024**2,
            "bytes_per_successful_session": next((row["bytes"] for row in level_rows if row["error"] is None), None),
        }
        print(concurrency, json.dumps(summary[str(concurrency)]), flush=True)

    payload = {
        "protocol": {
            "application": "SpatialAtlas Visium HD adrenal 8 um",
            "base_url": args.base_url,
            "concurrency_levels": [1, 5, 10, 20],
            "rounds_per_level": args.rounds,
            "request_bundle": list(PATHS),
            "cache": "disabled in urllib; configured display image transferred for every virtual user's first visit",
            "scope": "HTTP bootstrap and image transfer; browser parsing and WebGL are reported separately",
        },
        "summary": summary,
        "runs": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
