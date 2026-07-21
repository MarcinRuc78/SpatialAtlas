#!/usr/bin/env python3
"""Generate GraphPad-style adrenal benchmark and resolution plots."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


BENCH = Path(__file__).resolve().parent
PACKAGE = BENCH.parents[2]
EXAMPLE = PACKAGE / "examples" / "05_adrenal_platform_benchmark"
HD = PACKAGE / "examples" / "02_visium_hd_8um" / "output"
PLOTS = BENCH / "plots"

SCENARIOS = (
    "seurat_shiny_visium",
    "spatialatlas_visium",
    "seurat_shiny_visium_hd",
    "spatialatlas_visium_hd",
)
LABELS = (
    "Seurat/Shiny\nVisium",
    "SpatialAtlas\nVisium",
    "R/Shiny\nVisium HD",
    "SpatialAtlas\nVisium HD",
)
COLORS = ("#777777", "#2E75B6", "#D98E5F", "#D95F02")


def graphpad_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.linewidth": 2.2,
        "xtick.major.width": 2.0,
        "ytick.major.width": 2.0,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def style_numeric_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_linewidth(2.2)
    axis.spines["left"].set_linewidth(2.2)
    axis.tick_params(axis="both", direction="out", width=2.0, length=6)
    axis.grid(False)
    for label in axis.get_xticklabels() + axis.get_yticklabels():
        label.set_fontweight("bold")


def panel_label(axis: plt.Axes, letter: str) -> None:
    axis.text(-0.17, 1.07, letter, transform=axis.transAxes,
              fontsize=18, fontweight="bold", va="top", ha="left")


def draw_bar_panel(axis: plt.Axes, values: list[list[float]], ylabel: str,
                   letter: str, log_scale: bool = False) -> None:
    medians = np.array([np.median(v) for v in values])
    low = medians - np.array([np.min(v) for v in values])
    high = np.array([np.max(v) for v in values]) - medians
    x = np.arange(len(values))
    global_max = max(max(observations) for observations in values)
    axis.bar(x, medians, width=0.66, color=COLORS, edgecolor="black", linewidth=1.8, zorder=2)
    axis.errorbar(x, medians, yerr=np.vstack([low, high]), fmt="none",
                  ecolor="black", elinewidth=1.8, capsize=5, capthick=1.8, zorder=4)
    rng = np.random.default_rng(20260720)
    for index, observations in enumerate(values):
        jitter = rng.uniform(-0.12, 0.12, len(observations))
        axis.scatter(np.full(len(observations), index) + jitter, observations,
                     s=22, facecolor="white", edgecolor="black", linewidth=1.0, zorder=5)
        if log_scale:
            label_y = medians[index] / 1.8
            label_color = "white"
            vertical = "top"
        else:
            label_y = max(observations) + 0.025 * global_max
            label_color = "black"
            vertical = "bottom"
        axis.text(index, label_y, f"{medians[index]:.1f}", ha="center", va=vertical,
                  color=label_color, fontweight="bold", zorder=6)
    axis.set_xticks(x, LABELS)
    axis.set_ylabel(ylabel)
    if log_scale:
        axis.set_yscale("log")
        axis.set_ylim(
            bottom=min(min(observations) for observations in values) * 0.6,
            top=global_max * 1.25,
        )
    else:
        axis.set_ylim(bottom=0)
    style_numeric_axis(axis)
    panel_label(axis, letter)


def read_core(path: Path, field: str, scale: float) -> list[float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [float(row[field]) * scale for row in csv.DictReader(handle)]


def generate_benchmark() -> None:
    PLOTS.mkdir(exist_ok=True)
    startup = json.loads((BENCH / "startup_memory_results.json").read_text(encoding="utf-8"))
    startup_values = []
    memory_values = []
    for scenario in SCENARIOS:
        rows = [row for row in startup["runs"] if row["scenario"] == scenario and row["repetition"] > 0]
        startup_values.append([row["startup_to_http_200_s"] for row in rows])
        memory_values.append([row["steady_memory_bytes"] / 2**20 for row in rows])

    core_files = (
        BENCH / "shiny_core_results.csv",
        BENCH / "spatialatlas_classic_core_results.csv",
        BENCH / "shiny_hd_core_results.csv",
        BENCH / "spatialatlas_hd_core_results.csv",
    )
    response_values = [read_core(path, "total_server_s", 1000) for path in core_files]
    payload_values = [read_core(path, "response_bytes", 1 / 1024) for path in core_files]

    figure, axes = plt.subplots(2, 2, figsize=(11.4, 8.6), constrained_layout=True)
    draw_bar_panel(axes[0, 0], startup_values, "Warm startup to HTTP 200 (s)", "A")
    draw_bar_panel(axes[0, 1], memory_values, "Steady container memory (MiB)", "B")
    draw_bar_panel(axes[1, 0], response_values, "Gene response generation (ms; log scale)", "C", log_scale=True)
    draw_bar_panel(axes[1, 1], payload_values, "Serialized response (KiB; log scale)", "D", log_scale=True)
    figure.savefig(PLOTS / "adrenal_operational_benchmark.png", dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(PLOTS / "adrenal_operational_benchmark.svg", bbox_inches="tight", facecolor="white")
    plt.close(figure)


def expression_vector(atlas, gene: str) -> np.ndarray:
    from scipy import sparse

    matrix = atlas[:, gene].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix).reshape(-1)


def style_image_axis(axis: plt.Axes, title: str, letter: str) -> None:
    axis.set_title(title, pad=7, fontsize=13, fontweight="bold")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2.2)
        spine.set_color("black")
    panel_label(axis, letter)


def show_histology(axis: plt.Axes, image: np.ndarray, title: str, letter: str) -> None:
    axis.imshow(image, origin="upper")
    style_image_axis(axis, title, letter)


def show_expression(axis: plt.Axes, image: np.ndarray, coordinates: np.ndarray,
                    values: np.ndarray, title: str, letter: str,
                    size: float, vmax: float):
    axis.imshow(image, origin="upper", alpha=0.55)
    positive = values > 0
    plot = axis.scatter(coordinates[positive, 0], coordinates[positive, 1],
                        c=values[positive], cmap="magma", vmin=0, vmax=vmax,
                        s=size, edgecolors="none", alpha=0.92, rasterized=True)
    style_image_axis(axis, title, letter)
    return plot


def generate_resolution_figure() -> None:
    PLOTS.mkdir(exist_ok=True)
    import anndata as ad
    from PIL import Image

    classic = ad.read_h5ad(EXAMPLE / "prepared" / "sample_adrenal_visium_classic.h5ad")
    hd = ad.read_h5ad(HD / "sample_adrenal_female_hd_8um.h5ad")
    classic_image = np.asarray(Image.open(EXAMPLE / "prepared" / "tissue_adrenal_visium_classic.png"))
    hd_image = np.asarray(Image.open(HD / "tissue_adrenal_female_hd_8um.png"))
    classic_xy = np.asarray(classic.obsm["spatial"])
    hd_xy = np.asarray(hd.obsm["spatial"]) * 0.205592
    classic_expr = expression_vector(classic, "Cyp11b2")
    hd_expr = expression_vector(hd, "Cyp11b2")
    positive = np.concatenate([classic_expr[classic_expr > 0], hd_expr[hd_expr > 0]])
    vmax = float(np.percentile(positive, 99.5))

    figure, axes = plt.subplots(2, 2, figsize=(10.2, 9.4), constrained_layout=True)
    show_histology(axes[0, 0], classic_image,
                   "Classic Visium adrenal (992 spots)", "A")
    show_histology(axes[0, 1], hd_image,
                   "Visium HD adrenal (143,112 bins)", "B")
    plot = show_expression(axes[1, 0], classic_image, classic_xy, classic_expr,
                           "Cyp11b2 — classic Visium", "C", 25, vmax)
    show_expression(axes[1, 1], hd_image, hd_xy, hd_expr,
                    "Cyp11b2 — Visium HD, 8 µm", "D", 1.7, vmax)
    colorbar = figure.colorbar(plot, ax=axes.ravel().tolist(), orientation="horizontal",
                              fraction=0.035, pad=0.04, aspect=40)
    colorbar.set_label("Cyp11b2 expression (log-normalized; shared scale)", fontweight="bold")
    colorbar.outline.set_linewidth(1.8)
    colorbar.ax.tick_params(width=1.8, length=5)
    for label in colorbar.ax.get_xticklabels():
        label.set_fontweight("bold")
    figure.savefig(PLOTS / "adrenal_resolution_comparison.png", dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(PLOTS / "adrenal_resolution_comparison.svg", bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    graphpad_style()
    generate_benchmark()
    generate_resolution_figure()
