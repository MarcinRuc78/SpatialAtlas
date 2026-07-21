#!/usr/bin/env python3
"""Repeat container startup and steady-memory measurements for three examples."""

import json
import statistics
import subprocess
import time
import urllib.request
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[2]
BENCH = Path(__file__).resolve().parent
IMAGE = "spatialatlas-final:2026-07-20"
SCENARIOS = {
    "public_visium": (PACKAGE / "examples/01_public_visium", BENCH / "config_public_visium.yaml"),
    "visium_hd_8um": (PACKAGE / "examples/02_visium_hd_8um", BENCH / "config_visium_hd_8um.yaml"),
    "multisample_adrenal": (PACKAGE / "examples/03_multisample_adrenal", BENCH / "config_multisample_adrenal.yaml"),
}


def command(*args, check=True):
    return subprocess.run(args, check=check, text=True, capture_output=True)


def await_http(port, timeout=180):
    deadline = time.perf_counter() + timeout
    while time.perf_counter() < deadline:
        try:
            with urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=1) as response:
                if response.status == 200:
                    return
        except Exception:
            time.sleep(0.05)
    raise TimeoutError(f"No HTTP 200 on port {port}")


def run_once(scenario, repetition, port):
    dataset, config = SCENARIOS[scenario]
    name = f"atlas_bench_{scenario}_{repetition}"
    command("docker", "rm", "-f", name, check=False)
    args = [
        "docker", "run", "-d", "--name", name, "-p", f"{port}:8050",
        "-e", f"SPATIALATLAS_CONFIG=/bench/{config.name}",
        "-v", f"{dataset}:/dataset:ro", "-v", f"{BENCH}:/bench:ro",
        "-w", "/app", IMAGE, "gunicorn", "-b", "0.0.0.0:8050", "app:server",
        "--workers", "1", "--threads", "4", "--timeout", "120", "--preload",
    ]
    started = time.perf_counter()
    command(*args)
    try:
        await_http(port)
        startup = time.perf_counter() - started
        time.sleep(2)
        memory = int(command("docker", "exec", name, "cat", "/sys/fs/cgroup/memory.current").stdout)
        return {"scenario": scenario, "repetition": repetition,
                "startup_to_http_200_s": startup, "steady_memory_bytes": memory}
    finally:
        command("docker", "rm", "-f", name, check=False)


runs = []
for repetition in range(1, 6):
    order = list(SCENARIOS)
    if repetition % 2 == 0:
        order.reverse()
    for offset, scenario in enumerate(order):
        row = run_once(scenario, repetition, 18300 + repetition * 10 + offset)
        runs.append(row)
        print(json.dumps(row), flush=True)

summary = {}
for scenario in SCENARIOS:
    subset = [row for row in runs if row["scenario"] == scenario]
    summary[scenario] = {
        "n": len(subset),
        "startup_median_s": statistics.median(row["startup_to_http_200_s"] for row in subset),
        "startup_min_s": min(row["startup_to_http_200_s"] for row in subset),
        "startup_max_s": max(row["startup_to_http_200_s"] for row in subset),
        "steady_memory_median_bytes": int(statistics.median(row["steady_memory_bytes"] for row in subset)),
        "steady_memory_min_bytes": min(row["steady_memory_bytes"] for row in subset),
        "steady_memory_max_bytes": max(row["steady_memory_bytes"] for row in subset),
    }

(BENCH / "spatialatlas_three_examples_startup.json").write_text(
    json.dumps({"protocol": {"image": IMAGE, "workers": 1, "repetitions": 5,
                              "readiness": "docker run to HTTP 200",
                              "memory": "cgroup memory.current 2 s after readiness"},
                "runs": runs, "summary": summary}, indent=2) + "\n")
print(json.dumps(summary, indent=2))
