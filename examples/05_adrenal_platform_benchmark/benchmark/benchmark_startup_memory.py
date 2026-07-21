#!/usr/bin/env python3
"""Benchmark startup and steady memory for four adrenal-viewer workloads."""

from __future__ import annotations

import json
import statistics
import subprocess
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT.parents[1]
RESULTS = ROOT / "benchmark" / "startup_memory_results.json"
SPATIALATLAS_IMAGE = "spatialatlas-benchmark:amd64"
SHINY_IMAGE = "satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e"

SCENARIOS = {
    "seurat_shiny_visium": {
        "profiles": 992,
        "genes": 19_465,
        "kind": "shiny",
    },
    "spatialatlas_visium": {
        "profiles": 992,
        "genes": 19_465,
        "kind": "spatialatlas_classic",
    },
    "seurat_shiny_visium_hd": {
        "profiles": 143_112,
        "genes": 15_109,
        "kind": "shiny_hd",
    },
    "spatialatlas_visium_hd": {
        "profiles": 143_112,
        "genes": 15_109,
        "kind": "spatialatlas_hd",
    },
}


def command(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, check=check, text=True, capture_output=True)


def wait_http(port: int, timeout: float = 240.0) -> None:
    deadline = time.perf_counter() + timeout
    last_error = None
    while time.perf_counter() < deadline:
        try:
            with urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=1) as response:
                if response.status == 200:
                    return
        except Exception as exc:  # readiness polling deliberately tolerates connection errors
            last_error = repr(exc)
            time.sleep(0.05)
    raise TimeoutError(f"No HTTP 200 on port {port}: {last_error}")


def container_command(scenario: str, name: str, port: int) -> list[str]:
    kind = SCENARIOS[scenario]["kind"]
    base = [
        "docker", "run", "-d", "--platform", "linux/amd64",
        "--name", name, "-p", f"{port}:8050",
    ]
    if kind == "shiny":
        return base + [
            "-e", "SHINY_SEURAT_RDS=/work/source_data/adrenal_visium_seurat.rds",
            "-v", f"{ROOT}:/work:ro",
            "-w", "/work/benchmark/reference_app",
            SHINY_IMAGE,
            "R", "-q", "-e",
            "shiny::runApp('/work/benchmark/reference_app', host='0.0.0.0', port=8050, launch.browser=FALSE)",
        ]
    if kind == "shiny_hd":
        dataset = PACKAGE / "examples" / "02_visium_hd_8um" / "output"
        return base + [
            "-e", "SHINY_HD_RDS=/work/source_data/adrenal_visium_hd_shiny.rds",
            "-e", "SHINY_HD_IMAGE=/data/tissue_adrenal_female_hd_8um.png",
            "-v", f"{ROOT}:/work:ro",
            "-v", f"{dataset}:/data:ro",
            "-w", "/work/benchmark/reference_app_hd",
            SHINY_IMAGE,
            "R", "-q", "-e",
            "shiny::runApp('/work/benchmark/reference_app_hd', host='0.0.0.0', port=8050, launch.browser=FALSE)",
        ]
    if kind == "spatialatlas_classic":
        dataset = ROOT / "prepared"
        config = ROOT / "config.classic.docker.yaml"
    else:
        dataset = PACKAGE / "examples" / "02_visium_hd_8um" / "output"
        config = ROOT / "config.hd.docker.yaml"
    return base + [
        "-e", f"SPATIALATLAS_CONFIG=/work/{config.name}",
        "-v", f"{ROOT}:/work:ro",
        "-v", f"{dataset}:/data:ro",
        "-w", "/app",
        SPATIALATLAS_IMAGE,
        "gunicorn", "-b", "0.0.0.0:8050", "app:server",
        "--workers", "1", "--threads", "4", "--timeout", "180", "--preload",
    ]


def run_once(scenario: str, repetition: int, port: int) -> dict[str, object]:
    name = f"adrenal_bench_{scenario}_{repetition}"
    command("docker", "rm", "-f", name, check=False)
    start = time.perf_counter()
    command(*container_command(scenario, name, port))
    try:
        wait_http(port)
        startup = time.perf_counter() - start
        time.sleep(2)
        memory = int(command("docker", "exec", name, "cat", "/sys/fs/cgroup/memory.current").stdout)
        logs = command("docker", "logs", name, check=False)
        return {
            "scenario": scenario,
            "repetition": repetition,
            "startup_to_http_200_s": startup,
            "steady_memory_bytes": memory,
            "log_tail": (logs.stdout + logs.stderr)[-1500:],
        }
    finally:
        command("docker", "rm", "-f", name, check=False)


def summarize(rows: list[dict[str, object]], scenario: str) -> dict[str, object]:
    subset = [row for row in rows if row["scenario"] == scenario]
    warm = subset[1:]
    starts = [float(row["startup_to_http_200_s"]) for row in warm]
    memory = [int(row["steady_memory_bytes"]) for row in warm]
    return {
        "n_total": len(subset),
        "n_warm": len(warm),
        "first_start_s": subset[0]["startup_to_http_200_s"],
        "warm_start_median_s": statistics.median(starts),
        "warm_start_min_s": min(starts),
        "warm_start_max_s": max(starts),
        "warm_memory_median_bytes": int(statistics.median(memory)),
        "warm_memory_min_bytes": min(memory),
        "warm_memory_max_bytes": max(memory),
    }


def main() -> None:
    rows: list[dict[str, object]] = []
    names = list(SCENARIOS)
    for repetition in range(6):
        order = names[repetition % len(names):] + names[:repetition % len(names)]
        if repetition % 2:
            order.reverse()
        for offset, scenario in enumerate(order):
            row = run_once(scenario, repetition, 18400 + repetition * 10 + offset)
            rows.append(row)
            print(json.dumps({k: v for k, v in row.items() if k != "log_tail"}), flush=True)
    payload = {
        "protocol": {
            "host": "Apple M4 Pro",
            "platform": "linux/amd64 Docker for all four scenarios",
            "spatialatlas_image": SPATIALATLAS_IMAGE,
            "shiny_image": SHINY_IMAGE,
            "readiness": "wall time from docker run invocation to HTTP 200 at root",
            "memory": "container cgroup memory.current sampled 2 s after readiness",
            "runs_per_scenario": 6,
            "warm_summary": "first observation excluded; median and range of the next five",
            "spatialatlas_workers": 1,
        },
        "datasets": SCENARIOS,
        "runs": rows,
        "summary": {scenario: summarize(rows, scenario) for scenario in SCENARIOS},
    }
    RESULTS.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
