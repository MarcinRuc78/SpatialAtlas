#!/usr/bin/env python3
"""Measure startup-to-HD-viewer readiness and memory for TissUUmaps."""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import time
import urllib.request
from pathlib import Path


IMAGE = "cavenel/tissuumaps@sha256:e2f7898c700cd550432cfe546c14546681a4235c27b40c49cbb019e97ab86c06"
VIEWER = "adrenal_hd_tissuumaps_tmap.h5ad"


def command(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, check=check, text=True, capture_output=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runs", type=int, default=6)
    return parser.parse_args()


def run_once(prepared: Path, repetition: int) -> dict[str, object]:
    name = f"tissuumaps_hd_startup_{repetition}"
    port = 18100 + repetition
    command("docker", "rm", "-f", name, check=False)
    started = time.perf_counter()
    command(
        "docker", "run", "-d", "--platform", "linux/amd64", "--name", name,
        "-p", f"127.0.0.1:{port}:80", "-e", "HDF5_DISABLE_VERSION_CHECK=1",
        "-v", f"{prepared.resolve()}:/mnt/data/shared", IMAGE,
    )
    deadline = time.perf_counter() + 180
    error = None
    try:
        while time.perf_counter() < deadline:
            try:
                with urllib.request.urlopen(
                    f"http://127.0.0.1:{port}/{VIEWER}", timeout=10
                ) as response:
                    body = response.read()
                    if response.status == 200 and b"Gene expression" in body:
                        elapsed = time.perf_counter() - started
                        break
            except Exception as exc:
                error = repr(exc)
                time.sleep(0.05)
        else:
            raise TimeoutError(error)
        time.sleep(2)
        memory = int(command("docker", "exec", name, "cat", "/sys/fs/cgroup/memory.current").stdout)
        return {
            "repetition": repetition,
            "startup_to_hd_viewer_http_200_s": elapsed,
            "steady_memory_bytes": memory,
            "viewer_document_bytes": len(body),
        }
    finally:
        command("docker", "rm", "-f", name, check=False)


def main() -> None:
    args = parse_args()
    rows = []
    for repetition in range(args.runs):
        row = run_once(args.prepared_dir, repetition)
        rows.append(row)
        print(json.dumps(row), flush=True)
    warm = rows[1:]
    starts = [float(row["startup_to_hd_viewer_http_200_s"]) for row in warm]
    memory = [int(row["steady_memory_bytes"]) for row in warm]
    payload = {
        "protocol": {
            "image": IMAGE,
            "version": "3.2.1.14",
            "host": "Apple M4 Pro",
            "platform": "linux/amd64 under Docker Desktop",
            "readiness": "wall time from docker run invocation until HTTP 200 for the prepared CSC HD viewer document containing the 15,109-gene selector",
            "memory": "container cgroup memory.current two seconds after readiness",
            "runs": args.runs,
            "warm_summary": "first run retained but excluded; median and range of the following five",
            "server_default": "8 Gunicorn gevent workers as supplied by the official image"
        },
        "runs": rows,
        "summary": {
            "warm_start_median_s": statistics.median(starts),
            "warm_start_min_s": min(starts),
            "warm_start_max_s": max(starts),
            "warm_memory_median_bytes": int(statistics.median(memory)),
            "warm_memory_min_bytes": min(memory),
            "warm_memory_max_bytes": max(memory),
            "viewer_document_bytes": rows[-1]["viewer_document_bytes"]
        }
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
