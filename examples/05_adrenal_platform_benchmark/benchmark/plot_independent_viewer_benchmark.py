#!/usr/bin/env python3
"""Generate the GraphPad-style independent-viewer and concurrency plot."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


BENCH = Path(__file__).resolve().parent
PLOTS = BENCH / "plots"
COLORS = ("#2E75B6", "#6A3D9A")


def graphpad_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11, "font.weight": "bold",
        "axes.labelweight": "bold", "axes.titleweight": "bold", "axes.linewidth": 2.2,
        "xtick.major.width": 2.0, "ytick.major.width": 2.0,
        "xtick.major.size": 6, "ytick.major.size": 6,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def style_axis(axis: plt.Axes, letter: str) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(2.2)
    axis.spines["bottom"].set_linewidth(2.2)
    axis.tick_params(direction="out", width=2.0, length=6)
    axis.grid(False)
    for label in axis.get_xticklabels() + axis.get_yticklabels():
        label.set_fontweight("bold")
    axis.text(-0.18, 1.08, letter, transform=axis.transAxes,
              fontsize=18, fontweight="bold", va="top")


def bars(axis: plt.Axes, groups: list[list[float]], ylabel: str, letter: str) -> None:
    medians = np.asarray([np.median(group) for group in groups])
    lower = medians - np.asarray([min(group) for group in groups])
    upper = np.asarray([max(group) for group in groups]) - medians
    x = np.arange(2)
    axis.bar(x, medians, width=0.62, color=COLORS, edgecolor="black", linewidth=1.9, zorder=2)
    axis.errorbar(x, medians, yerr=np.vstack([lower, upper]), fmt="none", color="black",
                  linewidth=1.8, capsize=5, capthick=1.8, zorder=4)
    rng = np.random.default_rng(20260721)
    for index, values in enumerate(groups):
        jitter = rng.uniform(-0.10, 0.10, len(values))
        axis.scatter(np.full(len(values), index) + jitter, values, s=30, facecolor="white",
                     edgecolor="black", linewidth=1.1, zorder=5)
        axis.text(index, max(values) * 1.04, f"{medians[index]:.2f}", ha="center", fontweight="bold")
    axis.set_xticks(x, ("SpatialAtlas", "TissUUmaps"))
    axis.set_ylabel(ylabel)
    axis.set_ylim(0, max(max(group) for group in groups) * 1.18)
    style_axis(axis, letter)


def main() -> None:
    PLOTS.mkdir(exist_ok=True)
    graphpad_style()
    spatial = json.loads((BENCH / "startup_memory_results.json").read_text())
    tissuu = json.loads((BENCH / "tissuumaps_startup_memory_results.json").read_text())
    browser = json.loads((BENCH / "browser_ready_results.json").read_text())
    users = json.loads((BENCH / "concurrent_users_results.json").read_text())

    spatial_rows = [row for row in spatial["runs"]
                    if row["scenario"] == "spatialatlas_visium_hd" and row["repetition"] > 0]
    startup = [
        [float(row["startup_to_http_200_s"]) for row in spatial_rows],
        [float(row["startup_to_hd_viewer_http_200_s"]) for row in tissuu["runs"][1:]],
    ]
    memory = [
        [float(row["steady_memory_bytes"]) / 2**20 for row in spatial_rows],
        [float(row["steady_memory_bytes"]) / 2**20 for row in tissuu["runs"][1:]],
    ]
    ready = [
        [value / 1000 for value in browser["spatialatlas"]["ready_ms"]],
        [value / 1000 for value in browser["tissuumaps_prepared_csc"]["ready_ms"]],
    ]

    figure, axes = plt.subplots(2, 2, figsize=(10.8, 8.3), constrained_layout=True)
    bars(axes[0, 0], startup, "Warm server startup (s)", "A")
    bars(axes[0, 1], memory, "Steady container memory (MiB)", "B")
    bars(axes[1, 0], ready, "Browser visual readiness (s)", "C")

    levels = np.asarray([1, 5, 10, 20])
    p95 = np.asarray([users["summary"][str(level)]["session_p95_s"] for level in levels])
    medians = np.asarray([users["summary"][str(level)]["session_median_s"] for level in levels])
    axis = axes[1, 1]
    axis.plot(levels, p95, "o-", color="#D95F02", markeredgecolor="black",
              markeredgewidth=1.3, markersize=8, linewidth=3.0, label="95th percentile")
    axis.plot(levels, medians, "s--", color="#2E75B6", markeredgecolor="black",
              markeredgewidth=1.2, markersize=7, linewidth=2.5, label="Median")
    axis.set_xlabel("Concurrent first-visit sessions")
    axis.set_ylabel("SpatialAtlas load bundle (s)")
    axis.set_xticks(levels)
    axis.set_ylim(bottom=0)
    axis.legend(frameon=False, prop={"weight": "bold", "size": 10})
    style_axis(axis, "D")

    for suffix in ("png", "svg"):
        kwargs = {"dpi": 300} if suffix == "png" else {}
        figure.savefig(PLOTS / f"independent_viewer_browser_benchmark.{suffix}",
                       bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(figure)


if __name__ == "__main__":
    main()
