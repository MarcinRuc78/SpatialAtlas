"""
SpatialAtlas — Interactive Web Viewer for Visium HD Spatial Transcriptomics Data
================================================================================
A configurable Dash/Plotly application for deploying spatial transcriptomic
atlases as self-hosted web services. Reads all parameters from config.yaml.

Supports:
  - 1 to 4 samples displayed simultaneously (auto-layout or manual override)
  - Standard Visium (SpaceRanger v1.x / v2.x)
  - Visium HD at 2µm, 8µm, 16µm bin resolution
  - Visium HD with cell segmentation (Cellpose, StarDist)
  - Pre-processed .h5ad files from GEO, 10x datasets, or custom pipelines

Usage:
    python app.py                        # local development
    gunicorn -b 0.0.0.0:8050 app:server  # production via Docker
"""

import gc
import os
import yaml
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, no_update
import dash_bootstrap_components as dbc
from flask import abort, send_file
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.preprocessing import StandardScaler
from plotly.subplots import make_subplots
from PIL import Image

# Plotly graph config — hide broken client-side camera icon for GL traces
GRAPH_CONFIG = {
    "modeBarButtonsToRemove": ["toImage"],
    "displaylogo": False,
}

# =============================================================================
# Load configuration
# =============================================================================
# Look for config.yaml: env var → local dir → Docker path
_default_config = "config.yaml" if os.path.exists("config.yaml") else "/app/config.yaml"
CONFIG_PATH = os.getenv("SPATIALATLAS_CONFIG", _default_config)
CONFIG_PATH = os.path.abspath(CONFIG_PATH)
CONFIG_DIR = os.path.dirname(CONFIG_PATH)

with open(CONFIG_PATH, "r") as fh:
    cfg = yaml.safe_load(fh)

atlas_cfg = cfg["atlas"]
samples_cfg = cfg["samples"]
display_cfg = cfg.get("display", {})
multi_cfg = cfg.get("multi_gene", {})
citation_cfg = cfg.get("citation", {})
server_cfg = cfg.get("server", {})
markers_csv_path = cfg.get("markers_csv")


def resolve_config_path(path):
    """Resolve relative paths from the directory containing config.yaml."""
    if not path or os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(CONFIG_DIR, path))


markers_csv_path = resolve_config_path(markers_csv_path)


# =============================================================================
# Helper functions
# =============================================================================
def image_dimensions(path):
    """Read image dimensions without embedding image bytes in callback JSON."""
    with Image.open(path) as image:
        return image.size


def image_mimetype(path):
    """Return an explicit media type for common histology image formats."""
    return {
        ".jpeg": "image/jpeg",
        ".jpg": "image/jpeg",
        ".png": "image/png",
        ".tif": "image/tiff",
        ".tiff": "image/tiff",
        ".webp": "image/webp",
    }.get(os.path.splitext(path)[1].lower())


def load_sample(sample_cfg):
    """Load a single sample: AnnData + tissue image + aligned coordinates."""
    h5ad_path = resolve_config_path(sample_cfg["h5ad"])
    image_path = resolve_config_path(sample_cfg["image"])
    display_image_path = resolve_config_path(
        sample_cfg.get("display_image", sample_cfg["image"])
    )
    if not os.path.isfile(h5ad_path):
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")
    if not os.path.isfile(image_path):
        raise FileNotFoundError(f"Tissue image not found: {image_path}")
    if not os.path.isfile(display_image_path):
        raise FileNotFoundError(f"Display image not found: {display_image_path}")
    adata = ad.read_h5ad(h5ad_path)
    missing = []
    if "cluster" not in adata.obs:
        missing.append(".obs['cluster']")
    if "spatial" not in adata.obsm:
        missing.append(".obsm['spatial']")
    if missing:
        raise ValueError(
            f"Sample {sample_cfg['name']!r} is missing required fields: "
            + ", ".join(missing))
    # Coordinate bounds follow the archival image. The optional display image
    # may use a browser-efficient encoding but must preserve its aspect ratio.
    img_size = image_dimensions(image_path)
    display_size = image_dimensions(display_image_path)
    source_ratio = img_size[0] / img_size[1]
    display_ratio = display_size[0] / display_size[1]
    if not np.isclose(source_ratio, display_ratio, rtol=1e-3):
        raise ValueError(
            f"Display image aspect ratio differs from the source image: "
            f"{display_size} versus {img_size}"
        )

    # Align spatial coordinates to image pixels
    spatial = np.asarray(adata.obsm["spatial"])
    sfx = sample_cfg.get("scale_factor_x", 1.0)
    sfy = sample_cfg.get("scale_factor_y", 1.0)
    ox = sample_cfg.get("offset_x", 0)
    oy = sample_cfg.get("offset_y", 0)
    coords = np.column_stack([
        spatial[:, 0] * sfx - ox,
        spatial[:, 1] * sfy - oy,
    ]).astype(np.float32)

    # CSC provides efficient per-gene slicing. Assigning it back to AnnData
    # releases the original CSR arrays instead of retaining both formats.
    if sp.issparse(adata.X):
        X = adata.X.tocsc(copy=False).astype(np.float32, copy=False)
    else:
        X = np.asarray(adata.X, dtype=np.float32)
    adata.X = X

    var_names = np.asarray(adata.var_names.astype(str))
    clusters = adata.obs["cluster"].astype(str).to_numpy(copy=True)
    umap = (np.asarray(adata.obsm["UMAP"], dtype=np.float32)
            if "UMAP" in adata.obsm else None)
    n_obs = adata.n_obs
    gene_index = {gene: index for index, gene in enumerate(var_names)}

    return {
        "name": sample_cfg["name"],
        "n_obs": n_obs,
        "clusters": clusters,
        "umap": umap,
        "image_path": image_path,
        "display_image_path": display_image_path,
        "coords": coords,
        "X": X,
        "gene_index": gene_index,
        "img_bounds": {
            "xmin": sample_cfg.get("img_xmin", 0),
            "xmax": sample_cfg.get("img_xmax", img_size[0]),
            "ymin": sample_cfg.get("img_ymin", 0),
            "ymax": sample_cfg.get("img_ymax", img_size[1]),
        },
    }


# =============================================================================
# Load all samples
# =============================================================================
print("Loading samples...", flush=True)
if not 1 <= len(samples_cfg) <= 4:
    raise ValueError("config.yaml must define between 1 and 4 samples")
samples = [load_sample(s) for s in samples_cfg]
gc.collect()
n_samples = len(samples)

# Adaptive default point size based on data density
_total_spots = sum(s["coords"].shape[0] for s in samples)
_total_sparse_bytes = sum(
    (s["X"].data.nbytes + s["X"].indices.nbytes + s["X"].indptr.nbytes)
    if sp.issparse(s["X"]) else s["X"].nbytes
    for s in samples)
print(f"  Expression cache memory: {_total_sparse_bytes / (1024 ** 3):.2f} GiB", flush=True)
_default_pt_size = max(0.5, min(4.0, 2.5 * (4000 / max(_total_spots, 1)) ** 0.5))
_default_pt_size = round(_default_pt_size * 2) / 2  # round to 0.5 step
print(f"  Total spots: {_total_spots}, default point size: {_default_pt_size}px")

# Build lightweight combined metadata without materializing a union expression
# matrix. A gene that is absent from a sample is reported as unavailable for
# that sample; it is never substituted by a measured expression value of zero.
all_genes = sorted(set().union(*[set(s["gene_index"]) for s in samples]))
n_obs_combined = sum(s["n_obs"] for s in samples)
clusters_combined = np.concatenate([s["clusters"] for s in samples])
sample_names_combined = np.concatenate([
    np.full(s["n_obs"], s["name"], dtype=object) for s in samples
])

# UMAP coordinates (if available)
has_umap = all(s["umap"] is not None for s in samples)
if has_umap:
    umap_coords = np.vstack([s["umap"] for s in samples])

# Marker genes table (optional)
markers_df = None
if markers_csv_path and os.path.exists(markers_csv_path):
    markers_df = pd.read_csv(markers_csv_path)

# =============================================================================
# Gene availability across samples
# =============================================================================
# Samples processed independently can retain different gene sets: a gene may be
# undetected in a section, or removed by the minimum-observation filter applied
# during preparation. Neither situation is a measurement of zero expression, so
# availability is tracked explicitly and propagated to every view.
gene_presence = {
    gene: tuple(gene in s["gene_index"] for s in samples)
    for gene in all_genes
}
partial_genes = {g for g, present in gene_presence.items() if not all(present)}


def gene_present_in(gene, sample_idx):
    """True when the sample matrix actually contains the gene."""
    presence = gene_presence.get(gene)
    return bool(presence) and presence[sample_idx]


def missing_samples_for(gene):
    """Names of samples whose matrix does not contain the gene."""
    presence = gene_presence.get(gene)
    if not presence:
        return []
    return [s["name"] for s, ok in zip(samples, presence) if not ok]


def natural_list(names):
    """Join names the way a sentence would: "A", "A and B", "A, B and C"."""
    if len(names) <= 1:
        return "".join(names)
    return f"{', '.join(names[:-1])} and {names[-1]}"


def gene_option_label(gene):
    """Dropdown label; partial genes carry an explicit availability count."""
    if gene not in partial_genes:
        return gene
    n_present = sum(gene_presence[gene])
    return f"{gene}  —  in {n_present} of {n_samples} samples"


# Wording shown to readers. The panels say only what is true of that panel;
# the full explanation is given once, in the notice under the gene selector.
UNAVAILABLE_NOTE = "no data for this gene in this sample"

# Gene dropdown options
gene_names = list(all_genes)
gene_set = set(gene_names)
gene_options = [{"label": gene_option_label(g), "value": g} for g in gene_names]
max_dropdown = display_cfg.get("max_gene_dropdown", 500)

print(f"Atlas ready: {n_obs_combined} cells, {len(all_genes)} genes, {n_samples} samples", flush=True)
if partial_genes:
    print(f"  {len(partial_genes)} of {len(all_genes)} genes are absent from at "
          f"least one sample matrix and are reported as unavailable, not zero",
          flush=True)

# =============================================================================
# Subplot layout (supports 1-4 samples)
# =============================================================================
layout_mode = display_cfg.get("layout", "auto")
if layout_mode == "auto":
    LAYOUT = {1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (2, 2)}[min(n_samples, 4)]
else:
    r, c = layout_mode.split("x")
    LAYOUT = (int(r), int(c))
N_ROWS, N_COLS = LAYOUT


def subplot_index(sample_idx):
    """Return (row, col) for make_subplots (1-based)."""
    if N_ROWS == 1:
        return 1, sample_idx + 1
    return (sample_idx // N_COLS) + 1, (sample_idx % N_COLS) + 1


def axis_ref(sample_idx):
    """Return xref/yref strings for layout images."""
    flat = sample_idx  # 0-based
    if flat == 0:
        return "x", "y"
    return f"x{flat + 1}", f"y{flat + 1}"


# =============================================================================
# Color scales
# =============================================================================
COLOR_SCALES = {
    "Viridis": "Viridis",
    "Plasma": "Plasma",
    "Cividis": "Cividis",
    "Turbo": "Turbo",
    "Transparent-Black":   [[0, "rgba(0,0,0,0)"],     [1, "rgba(0,0,0,1)"]],
    "Transparent-Green":   [[0, "rgba(0,255,0,0)"],   [1, "rgba(0,255,0,1)"]],
    "Transparent-Red":     [[0, "rgba(255,0,0,0)"],   [1, "rgba(255,0,0,1)"]],
    "Transparent-Blue":    [[0, "rgba(0,0,255,0)"],   [1, "rgba(0,0,255,1)"]],
    "Transparent-Magenta": [[0, "rgba(255,0,255,0)"], [1, "rgba(255,0,255,1)"]],
    "Transparent-Cyan":    [[0, "rgba(0,255,255,0)"], [1, "rgba(0,255,255,1)"]],
    "Transparent-Yellow":  [[0, "rgba(255,255,0,0)"], [1, "rgba(255,255,0,1)"]],
}

MULTI_GENE_COLORS = [
    "Transparent-Green", "Transparent-Red", "Transparent-Cyan", "Transparent-Yellow"
]

MULTI_GENE_LABELS = ["Gene A (green)", "Gene B (red)", "Gene C (cyan)", "Gene D (yellow)"]
MULTI_GENE_BG = ["#e8f5e9", "#ffebee", "#e0f7fa", "#fffde7"]


# =============================================================================
# Dash App Layout
# =============================================================================
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


@server.route("/_spatialatlas/tissue/<int:sample_index>")
def tissue_image(sample_index):
    """Serve configured tissue images with browser/proxy caching enabled."""
    if sample_index < 0 or sample_index >= len(samples):
        abort(404)
    return send_file(
        samples[sample_index]["display_image_path"],
        conditional=True,
        max_age=3600,
        mimetype=image_mimetype(samples[sample_index]["display_image_path"]),
    )


@server.route("/_spatialatlas/tissue-full/<int:sample_index>")
def full_tissue_image(sample_index):
    """Serve the archival histology image without changing display defaults."""
    if sample_index < 0 or sample_index >= len(samples):
        abort(404)
    return send_file(
        samples[sample_index]["image_path"],
        conditional=True,
        max_age=3600,
        mimetype=image_mimetype(samples[sample_index]["image_path"]),
    )

# --- Title ---
title_parts = [html.Span(atlas_cfg["title"])]
if atlas_cfg.get("subtitle"):
    title_parts += [html.Br(), html.Span(atlas_cfg["subtitle"])]

header = html.Div([
    html.H1(title_parts, style={
        "textAlign": "center", "color": "#2c3e50",
        "fontFamily": "Arial", "fontSize": "32px", "fontWeight": "bold"
    }),
])
_logo_path = resolve_config_path(atlas_cfg.get("logo"))
if _logo_path:
    # Resolve logo: check local path first, then /app/assets/
    import base64 as _b64
    _logo_local = _logo_path if os.path.exists(_logo_path) else None
    if not _logo_local and os.path.exists("assets/logo.png"):
        _logo_local = "assets/logo.png"
    if _logo_local:
        with open(_logo_local, "rb") as _lf:
            _logo_b64 = "data:image/png;base64," + _b64.b64encode(_lf.read()).decode()
        header.children.append(
            html.Img(src=_logo_b64, style={
                "position": "absolute", "top": "0px", "right": "30px", "height": "100px"
            })
        )
        header.style = {"position": "relative"}

# --- Gene selector ---
gene_selector = html.Div([
    html.Div([
        html.Label("Genes:", style={"fontWeight": "bold", "marginRight": "15px", "lineHeight": "40px"}),
        dcc.Dropdown(
            id="gene", options=gene_options[:max_dropdown], value=None,
            placeholder="Type gene name...", style={"width": "300px"},
            optionHeight=30, searchable=True
        )
    ], style={"display": "flex", "justifyContent": "center", "alignItems": "center"}),
    html.Div(id="gene-availability-note", style={
        "marginTop": "8px", "textAlign": "center", "fontSize": "13px",
    }),
], id="gene-selector-div", style={"marginBottom": "20px", "textAlign": "center"})

# --- Spatial tab ---
spatial_tab = dcc.Tab(label="Spatial", value="tab-spatial", children=[
    html.Div("Adjust point size, opacity and color scale:", style={"margin": "10px 0"}),
    html.Div([
        html.Div([
            html.Label("Point size (px)"),
            dcc.Slider(id="spatial-size", min=0.5, max=40, step=0.5,
                       value=display_cfg.get("default_point_size", _default_pt_size),
                       marks={i: str(i) for i in [1, 5, 10, 20, 30, 40]},
                       tooltip={"always_visible": True}),
        ], style={"width": "30%"}),
        html.Div([
            html.Label("Opacity"),
            dcc.Slider(id="spatial-opacity", min=0.1, max=1, step=0.1,
                       value=display_cfg.get("default_opacity", 0.8),
                       marks={i/10: str(i/10) for i in range(1, 11)},
                       tooltip={"always_visible": True}),
        ], style={"width": "30%"}),
        html.Div([
            html.Label("Color scale"),
            dcc.Dropdown(id="color-scale",
                         options=[{"label": k, "value": k} for k in COLOR_SCALES],
                         value=display_cfg.get("default_colorscale", "Transparent-Green"),
                         clearable=False),
        ], style={"width": "30%"}),
    ], style={"display": "flex", "justifyContent": "space-between", "gap": "5%"}),
    dcc.Graph(id="spatial-plot", config=GRAPH_CONFIG, style={"height": f"{display_cfg.get('spatial_plot_height', 850)}px"}),
])

# --- UMAP tab ---
umap_tab = dcc.Tab(label="UMAP", value="tab-umap", children=[
    dcc.Graph(id="umap-plot", config=GRAPH_CONFIG, style={"height": "500px"}),
]) if has_umap else None

# --- Violin tab ---
violin_tab = dcc.Tab(label="Violin", value="tab-violin", children=[
    dcc.Graph(id="violin-plot", config=GRAPH_CONFIG),
])

# --- Heatmap tab ---
heatmap_tab = dcc.Tab(label="Heatmap", value="tab-heatmap", children=[
    html.Div("Top N markers per cluster:", style={"margin": "10px 0"}),
    dcc.Slider(id="n-markers", min=2, max=20, step=1, value=5,
               marks={i: str(i) for i in range(2, 21)},
               tooltip={"always_visible": True}),
    dcc.Graph(id="heatmap-plot", config=GRAPH_CONFIG),
]) if markers_df is not None else None

# --- Multi-gene overlay tab ---
# Default genes follow the configuration, then the marker table, then gene order.
_cfg_genes = [multi_cfg.get("gene_a", ""), multi_cfg.get("gene_b", ""),
              multi_cfg.get("gene_c", ""), multi_cfg.get("gene_d", "")]
_cfg_genes = [g for g in _cfg_genes if g and g in gene_set]
if len(_cfg_genes) < 4:
    _extra = []
    if markers_df is not None and "gene" in markers_df.columns:
        # Pick top unique marker genes from different clusters
        _seen = set(_cfg_genes)
        for _, row in markers_df.iterrows():
            g = row["gene"]
            if g in gene_set and g not in _seen:
                _extra.append(g)
                _seen.add(g)
            if len(_extra) >= 4 - len(_cfg_genes):
                break
    if len(_cfg_genes) + len(_extra) < 4:
        # Keep the fallback deterministic so repeated deployments open identically.
        _candidates = [g for g in gene_names if g not in set(_cfg_genes + _extra)]
        _extra.extend(_candidates[:4 - len(_cfg_genes) - len(_extra)])
    _cfg_genes.extend(_extra[:4 - len(_cfg_genes)])
default_genes = (_cfg_genes + ["", "", "", ""])[:4]
print(f"  Multi-gene defaults: {default_genes}", flush=True)
multi_initial_options = list(gene_options[:max_dropdown])
for value in default_genes:
    if value and not any(option["value"] == value for option in multi_initial_options):
        multi_initial_options.append({"label": value, "value": value})

multi_gene_rows = []
for pair_start in range(0, 4, 2):
    row_children = []
    for idx in range(pair_start, min(pair_start + 2, 4)):
        row_children.append(html.Div([
            dcc.Checklist(id=f"mgene{idx+1}-enable",
                          options=[{"label": " ", "value": "on"}], value=["on"],
                          style={"display": "inline-block", "marginRight": "5px"}),
            html.Label(MULTI_GENE_LABELS[idx],
                       style={"fontWeight": "bold", "marginRight": "10px", "minWidth": "120px"}),
            dcc.Dropdown(id=f"mgene{idx+1}", options=multi_initial_options,
                         value=default_genes[idx] if default_genes[idx] in gene_set else None,
                         clearable=False,
                         style={"width": "180px", "display": "inline-block", "marginRight": "15px"},
                         optionHeight=50),
            html.Label("Opacity:", style={"marginRight": "10px", "whiteSpace": "nowrap"}),
            dcc.Input(id=f"mgene{idx+1}-opacity", type="number", min=0.1, max=1, step=0.1,
                      value=0.8, style={"width": "80px"}),
        ], style={"display": "flex", "alignItems": "center", "padding": "10px",
                  "backgroundColor": MULTI_GENE_BG[idx], "borderRadius": "5px", "width": "48%"}))
    multi_gene_rows.append(
        html.Div(row_children, style={
            "display": "flex", "justifyContent": "space-between",
            "gap": "2%", "marginBottom": "15px"
        })
    )

multi_tab = dcc.Tab(label="Multiple Genes", value="tab-multiple", children=[
    html.Div(" ", style={"margin": "10px 0"}),
    *multi_gene_rows,
    html.Div([
        html.Label("Point size (px)"),
        dcc.Slider(id="multi-size", min=0.5, max=40, step=0.5, value=3,
                   marks={i: str(i) for i in [1, 5, 10, 20, 30, 40]},
                   tooltip={"always_visible": True}),
    ], style={"marginTop": "20px", "marginBottom": "20px"}),
    dcc.Graph(id="multi-genes-plot", config=GRAPH_CONFIG, style={"height": f"{display_cfg.get('spatial_plot_height', 850)}px"}),
])

# --- Citation tab ---
citation_tab = dcc.Tab(label="Citation", value="tab-citation", children=[
    html.Div([
        html.H3("Citation", style={"marginTop": "20px"}),
        html.P(citation_cfg.get("text", "")),
        html.P([
            html.B(citation_cfg.get("paper_title", "")), html.Br(),
            citation_cfg.get("authors", ""), html.Br(),
            "DOI: ", html.A(citation_cfg.get("doi", ""), href=citation_cfg.get("doi", ""),
                            target="_blank"),
        ], style={"maxWidth": "800px"})
    ], style={"padding": "20px"})
])

# Assemble tabs
tabs_children = [spatial_tab]
if umap_tab is not None:
    tabs_children.append(umap_tab)
tabs_children.append(violin_tab)
if heatmap_tab is not None:
    tabs_children.append(heatmap_tab)
tabs_children.append(multi_tab)
tabs_children.append(citation_tab)

app.layout = html.Div([header, gene_selector,
                        dcc.Tabs(id="tabs-main", value="tab-spatial", children=tabs_children)])


# =============================================================================
# Helper: get expression vector for a gene
# =============================================================================
def get_expr(gene, missing=np.nan):
    """Combined expression vector for one gene.

    Observations belonging to a sample whose matrix does not contain the gene
    are filled with ``missing`` (NaN by default) rather than zero, so that
    "not retained in this sample" stays distinguishable from a measured zero
    in every downstream plot and summary.
    """
    if not gene or gene not in gene_set:
        return np.full(n_obs_combined, missing, dtype=np.float32)
    chunks = []
    for sample in samples:
        idx = sample["gene_index"].get(gene)
        if idx is None:
            chunks.append(np.full(sample["n_obs"], missing, dtype=np.float32))
            continue
        column = sample["X"][:, idx]
        if sp.issparse(column):
            column = column.toarray()
        chunks.append(np.asarray(column, dtype=np.float32).reshape(-1))
    return np.concatenate(chunks)


def expression_vmax(expr):
    """99th percentile of measured positive values, ignoring unavailable data."""
    measured = expr[np.isfinite(expr)]
    positive = measured[measured > 0]
    if positive.size == 0:
        return 1.0
    return max(float(np.percentile(positive, 99)), 0.1)


def add_extent_placeholder(fig, sample, row, col):
    """Keep the histology visible in a subplot that carries no data trace.

    Layout images are positioned in data coordinates, so a subplot without any
    trace would autoscale to an unrelated range and hide the tissue image.
    """
    bounds = sample["img_bounds"]
    fig.add_trace(go.Scattergl(
        x=[bounds["xmin"], bounds["xmax"]],
        y=[bounds["ymin"], bounds["ymax"]],
        mode="markers", marker=dict(size=0.1, opacity=0),
        hoverinfo="skip", showlegend=False,
    ), row=row, col=col)


def unavailable_annotation(fig, gene, row, col):
    """Mark a subplot whose sample has no data for the requested gene."""
    fig.add_annotation(
        text=f"<b>{gene}</b><br>{UNAVAILABLE_NOTE}",
        xref="x domain", yref="y domain",
        x=0.5, y=0.5, showarrow=False, align="center",
        font=dict(size=11, color="#5a3d00"),
        bgcolor="rgba(255,244,214,0.92)", bordercolor="#b8860b", borderwidth=1,
        borderpad=5, row=row, col=col,
    )


def subplot_title(sample, sample_idx, gene):
    """Subplot title, flagged when the sample does not carry the gene."""
    if gene and gene in partial_genes and not gene_present_in(gene, sample_idx):
        return (f"<b>{sample['name']}</b>"
                f"<br><span style='font-size:11px;color:#8a6d00'>"
                f"no data for {gene}</span>")
    return f"<b>{sample['name']}</b>"


def sample_slice(sample_idx):
    """Return (start, end) indices in combined vectors for a sample."""
    start = sum(s["n_obs"] for s in samples[:sample_idx])
    end = start + samples[sample_idx]["n_obs"]
    return start, end


# =============================================================================
# Callbacks
# =============================================================================

# --- Dynamic gene filtering ---
def gene_search_results(search, current_value):
    query = (search or "").strip().lower()
    if query:
        starts = [o for o in gene_options if o["value"].lower().startswith(query)]
        contains = [o for o in gene_options
                    if not o["value"].lower().startswith(query)
                    and query in o["value"].lower()]
        results = (starts + contains)[:200]
    else:
        results = list(gene_options[:max_dropdown])
    if current_value and not any(o["value"] == current_value for o in results):
        results = [{"label": gene_option_label(current_value),
                    "value": current_value}] + results
    return results


@app.callback(
    Output("gene", "options"),
    Output("gene", "value"),
    Input("gene", "search_value"),
    State("gene", "value"),
    prevent_initial_call=True
)
def filter_genes(search, current_value):
    s = (search or "").strip()
    s_lower = s.lower()
    if not s:
        return no_update, no_update

    results = gene_search_results(search, current_value)

    new_value = no_update
    if s in gene_set:
        new_value = s
    else:
        match = next((g for g in gene_set if g.lower() == s_lower), None)
        if match:
            new_value = match

    check_val = new_value if new_value is not no_update else current_value
    if check_val and not any(o["value"] == check_val for o in results):
        cur = next((o for o in gene_options if o["value"] == check_val), None)
        if cur:
            results = [cur] + results

    return results, new_value


for _index in range(1, 5):
    @app.callback(
        Output(f"mgene{_index}", "options"),
        Input(f"mgene{_index}", "search_value"),
        State(f"mgene{_index}", "value"),
        prevent_initial_call=True,
    )
    def filter_multi_gene(search, current_value):
        return gene_search_results(search, current_value)


# --- Toggle gene selector visibility ---
@app.callback(
    Output("gene-selector-div", "style"),
    Input("tabs-main", "value")
)
def toggle_gene_selector(active_tab):
    if active_tab in ["tab-spatial", "tab-umap", "tab-violin", "tab-heatmap"]:
        return {"marginBottom": "20px", "textAlign": "center"}
    return {"display": "none"}


# --- Gene availability banner ---
@app.callback(
    Output("gene-availability-note", "children"),
    Input("gene", "value"),
)
def gene_availability_note(gene):
    """Name the samples whose data does not contain the selected gene."""
    if not gene or gene not in partial_genes:
        return None
    absent = missing_samples_for(gene)
    return html.Div([
        html.B(f"{gene} is missing from {natural_list(absent)}"),
    ], style={
        "display": "inline-block", "maxWidth": "780px", "padding": "8px 12px",
        "backgroundColor": "#fff4d6", "border": "1px solid #b8860b",
        "borderRadius": "4px", "textAlign": "left", "color": "#5a3d00",
    })


# --- UMAP plot ---
if has_umap:
    @app.callback(
        Output("umap-plot", "figure"),
        Input("gene", "value"),
        Input("tabs-main", "value"),
    )
    def update_umap(gene, active_tab="tab-umap"):
        if active_tab != "tab-umap":
            return no_update
        expr = get_expr(gene)
        gene_label = gene or "expression"
        # Observations from samples without the gene are not colored on the
        # expression scale; they are drawn separately and labelled unavailable.
        measured = np.isfinite(expr)
        df = pd.DataFrame({
            "UMAP1": umap_coords[measured, 0], "UMAP2": umap_coords[measured, 1],
            "cluster": clusters_combined[measured], gene_label: expr[measured]
        })
        fig = px.scatter(df, x="UMAP1", y="UMAP2", color=gene_label,
                         hover_data=["cluster"], opacity=0.8,
                         color_continuous_scale="Viridis", render_mode="webgl")
        fig.update_traces(marker=dict(size=3))
        fig.update_coloraxes(colorbar_title="Expression<br>[arb. units]")
        fig.update_layout(template="plotly_white")
        if not measured.all():
            absent = ", ".join(missing_samples_for(gene))
            fig.add_trace(go.Scattergl(
                x=umap_coords[~measured, 0], y=umap_coords[~measured, 1],
                mode="markers", marker=dict(size=3, color="#c8c8c8", opacity=0.6),
                name=f"no data for {gene} ({absent})",
                hovertemplate=f"{gene}: {UNAVAILABLE_NOTE}<extra></extra>",
            ))
            fig.data = fig.data[-1:] + fig.data[:-1]
            fig.update_layout(showlegend=True, legend=dict(
                orientation="h", yanchor="bottom", y=1.02, x=0))
        # Centroids use every observation so that cluster labels stay in the
        # same place regardless of which samples carry the selected gene.
        cents = pd.DataFrame({
            "UMAP1": umap_coords[:, 0], "UMAP2": umap_coords[:, 1],
            "cluster": clusters_combined,
        }).groupby("cluster", observed=False)[["UMAP1", "UMAP2"]].median().reset_index()
        for _, row in cents.iterrows():
            fig.add_annotation(
                x=row["UMAP1"], y=row["UMAP2"], text=row["cluster"],
                showarrow=False, font=dict(size=14, color="black"),
                bgcolor="white", bordercolor="black"
            )
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)
        return fig


# --- Spatial plot ---
@app.callback(
    Output("spatial-plot", "figure"),
    [Input("gene", "value"), Input("spatial-size", "value"),
     Input("spatial-opacity", "value"), Input("color-scale", "value"),
     Input("tabs-main", "value")]
)
def update_spatial(gene, size, opacity, color_scale_name, active_tab="tab-spatial"):
    if active_tab != "tab-spatial":
        return no_update
    hspacing = 0.02 if N_COLS <= 2 else 0.01
    vspacing = 0.05 if N_ROWS > 1 else 0.02
    fig = make_subplots(
        rows=N_ROWS, cols=N_COLS,
        subplot_titles=[subplot_title(s, si, gene) for si, s in enumerate(samples)],
        horizontal_spacing=hspacing,
        vertical_spacing=vspacing,
    )

    if not gene:
        # Cluster mode
        all_clusters = sorted(set().union(*[set(s["clusters"]) for s in samples]))
        colors = px.colors.qualitative.Plotly + px.colors.qualitative.Dark24
        color_map = {c: colors[i % len(colors)] for i, c in enumerate(all_clusters)}
        shown = set()
        for si, s in enumerate(samples):
            row, col = subplot_index(si)
            clusters = s["clusters"]
            for cl in all_clusters:
                mask = clusters == cl
                if not mask.any():
                    continue
                fig.add_trace(go.Scattergl(
                    x=s["coords"][mask, 0], y=s["coords"][mask, 1], mode="markers",
                    marker=dict(size=size, color=color_map[cl], opacity=opacity),
                    name=cl, legendgroup=cl,
                    showlegend=(cl not in shown),
                    hovertemplate=f"<b>{cl}</b><extra></extra>"
                ), row=row, col=col)
                shown.add(cl)
    else:
        # Expression mode
        cscale = COLOR_SCALES[color_scale_name]
        expr = get_expr(gene)
        # Use 99th percentile as vmax to prevent outliers from washing out
        # signal. Unavailable samples contribute NaN and are excluded.
        vmax = expression_vmax(expr)
        available = [si for si in range(n_samples) if gene_present_in(gene, si)]
        scale_owner = available[-1] if available else None

        for si, s in enumerate(samples):
            row, col = subplot_index(si)
            if not gene_present_in(gene, si):
                # The matrix does not contain this gene. No expression is
                # drawn and the panel is labelled, so an unmeasured gene is
                # never shown as an all-zero expression profile.
                add_extent_placeholder(fig, s, row, col)
                unavailable_annotation(fig, gene, row, col)
                continue
            start, end = sample_slice(si)
            e = expr[start:end]
            mask = np.isfinite(e) & (e > 0)
            is_last = (si == scale_owner)
            fig.add_trace(go.Scattergl(
                x=s["coords"][mask, 0], y=s["coords"][mask, 1], mode="markers",
                marker=dict(size=size, color=e[mask], colorscale=cscale,
                            opacity=opacity, showscale=is_last,
                            cmin=0, cmax=vmax,
                            colorbar=dict(title="Expression<br>[arb. units]", x=1.01) if is_last else None),
                hoverinfo="skip", showlegend=False
            ), row=row, col=col)

    # Add tissue images and style
    for si, s in enumerate(samples):
        b = s["img_bounds"]
        xr, yr = axis_ref(si)
        fig.add_layout_image(dict(
            source=f"/_spatialatlas/tissue/{si}",
            x=b["xmin"], y=b["ymin"],
            sizex=b["xmax"] - b["xmin"], sizey=b["ymax"] - b["ymin"],
            xref=xr, yref=yr,
            sizing="stretch", layer="below"
        ))

    fig.update_traces(marker=dict(size=size))
    fig.update_xaxes(showgrid=False, showticklabels=False, title_text=None)
    fig.update_yaxes(showgrid=False, showticklabels=False, title_text=None, autorange="reversed")
    for si in range(n_samples):
        row, col = subplot_index(si)
        _, yr = axis_ref(si)
        fig.update_xaxes(scaleanchor=yr, scaleratio=1, row=row, col=col)

    base_height = display_cfg.get("spatial_plot_height", 850)
    plot_height = base_height if N_ROWS == 1 else int(base_height * 0.55 * N_ROWS)
    fig.update_layout(
        template="plotly_white", height=plot_height,
        margin=dict(l=10, r=10, t=40, b=10),
        legend=dict(x=1.01, y=1, xanchor="left")
    )
    return fig


# --- Violin plot ---
@app.callback(
    Output("violin-plot", "figure"),
    Input("gene", "value"),
    Input("tabs-main", "value"),
)
def update_violin(gene, active_tab="tab-violin"):
    if active_tab != "tab-violin":
        return no_update
    if not gene:
        return go.Figure()
    expr = get_expr(gene)
    clusters = clusters_combined
    sample_names = sample_names_combined
    cluster_order = sorted(set(clusters))
    unique_samples = [s["name"] for s in samples]
    palette = px.colors.qualitative.Plotly
    colors = {name: palette[i % len(palette)] for i, name in enumerate(unique_samples)}

    fig = go.Figure()
    sides = ["negative", "positive"] if len(unique_samples) == 2 else [None] * len(unique_samples)
    for si, sname in enumerate(unique_samples):
        # A sample whose matrix lacks the gene contributes no distribution.
        # Omitting it keeps an unmeasured gene out of the quantitative summary
        # instead of adding a spurious distribution concentrated at zero.
        if not gene_present_in(gene, si):
            continue
        mask_s = sample_names == sname
        for ci, cl in enumerate(cluster_order):
            mask = mask_s & (clusters == cl)
            vals = expr[mask]
            vals = vals[np.isfinite(vals)]
            fig.add_trace(go.Violin(
                x=[cl] * len(vals), y=vals,
                name=sname, legendgroup=sname,
                showlegend=(ci == 0),
                side=sides[si] if si < len(sides) else None,
                line_color=colors[sname], fillcolor=colors[sname],
                opacity=0.7, meanline_visible=True, points=False, box_visible=True,
            ))

    absent = missing_samples_for(gene)
    if absent:
        fig.add_annotation(
            text=f"<b>{gene}</b> is missing from {natural_list(absent)}",
            xref="paper", yref="paper", x=0.5, y=1.10, showarrow=False,
            font=dict(size=11, color="#5a3d00"),
            bgcolor="rgba(255,244,214,0.92)", bordercolor="#b8860b",
            borderwidth=1, borderpad=4,
        )

    fig.update_layout(
        violingap=0, violinmode="overlay",
        template="plotly_white", legend_title="Sample",
        yaxis_title="Expression [arb. units]", height=500,
        margin=dict(t=90 if absent else 40),
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    return fig


# --- Heatmap ---
if markers_df is not None:
    @app.callback(
        Output("heatmap-plot", "figure"),
        Input("n-markers", "value"),
        Input("tabs-main", "value"),
    )
    def update_heatmap(n_markers, active_tab="tab-heatmap"):
        if active_tab != "tab-heatmap":
            return no_update
        if "avg_log2FC" in markers_df.columns:
            top = markers_df.sort_values("avg_log2FC", ascending=False).groupby("cluster", observed=False).head(n_markers)
        else:
            top = markers_df.groupby("cluster", observed=False).head(n_markers)
        genes = [g for g in top["gene"].unique() if g in gene_set]
        if len(genes) < 2:
            return go.Figure().update_layout(
                template="plotly_white",
                title="Not enough marker genes overlap the expression matrix")
        # Unavailable samples enter as NaN, so each cluster mean is taken over
        # the samples that actually measured the gene rather than being diluted
        # by zeros substituted for absent genes.
        df = pd.DataFrame({gene: get_expr(gene) for gene in genes})
        df["cluster"] = clusters_combined
        mat = df.groupby("cluster", observed=False).mean()
        mat = mat.fillna(mat.mean())
        scaler = StandardScaler()
        mat_scaled = pd.DataFrame(scaler.fit_transform(mat), index=mat.index, columns=mat.columns)
        link_r = linkage(mat_scaled.values, method="ward")
        link_c = linkage(mat_scaled.T.values, method="ward")
        mat_scaled = mat_scaled.iloc[leaves_list(link_r), leaves_list(link_c)]
        # Genes measured in only part of the atlas are flagged on the axis so
        # that a partial cluster mean is not read as a complete one.
        column_labels = [f"{g} *" if g in partial_genes else g
                         for g in mat_scaled.columns]
        fig = go.Figure(data=go.Heatmap(
            z=mat_scaled.values, x=column_labels, y=mat_scaled.index,
            colorscale="Viridis", colorbar=dict(title="Expression<br>[arb. units]")
        ))
        title = f"Top {n_markers} markers heatmap (scaled by gene)"
        if any(g in partial_genes for g in mat_scaled.columns):
            title += ("<br><span style='font-size:11px;color:#8a6d00'>"
                      "* not present in every sample; the average uses only "
                      "the samples that have it</span>")
        fig.update_layout(template="plotly_white", title=title)
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)
        return fig


# --- Multi-gene overlay ---
@app.callback(
    Output("multi-genes-plot", "figure"),
    [Input(f"mgene{i+1}", "value") for i in range(4)] +
    [Input(f"mgene{i+1}-enable", "value") for i in range(4)] +
    [Input(f"mgene{i+1}-opacity", "value") for i in range(4)] +
    [Input("multi-size", "value"), Input("tabs-main", "value")]
)
def update_multi(*args):
    active_tab = args[13] if len(args) > 13 else "tab-multiple"
    if active_tab != "tab-multiple":
        return no_update
    genes = args[0:4]
    enables = args[4:8]
    opacities = args[8:12]
    size = args[12]

    hspacing = 0.02 if N_COLS <= 2 else 0.01
    vspacing = 0.05 if N_ROWS > 1 else 0.02
    fig = make_subplots(
        rows=N_ROWS, cols=N_COLS,
        subplot_titles=[f"<b>{s['name']}</b>" for s in samples],
        horizontal_spacing=hspacing,
        vertical_spacing=vspacing,
    )

    unavailable = {si: [] for si in range(n_samples)}
    for gi in range(4):
        gene = genes[gi]
        if not enables[gi] or "on" not in enables[gi]:
            continue
        expr = get_expr(gene or "")
        cscale = COLOR_SCALES[MULTI_GENE_COLORS[gi]]
        op = opacities[gi] or 0.8

        for si, s in enumerate(samples):
            row, col = subplot_index(si)
            if gene and not gene_present_in(gene, si):
                # No channel is drawn for a gene the sample does not contain.
                unavailable[si].append(gene)
                continue
            start, end = sample_slice(si)
            e = expr[start:end]
            mask = np.isfinite(e) & (e > 0)
            fig.add_trace(go.Scattergl(
                x=s["coords"][mask, 0], y=s["coords"][mask, 1], mode="markers",
                marker=dict(size=size, color=e[mask], colorscale=cscale,
                            opacity=op, showscale=False),
                name=gene, hoverinfo="skip"
            ), row=row, col=col)

    drawn = {trace.xaxis for trace in fig.data}
    for si, absent in unavailable.items():
        row, col = subplot_index(si)
        axis = "x" if si == 0 else f"x{si + 1}"
        if axis not in drawn:
            # Every requested channel is missing for this sample; keep the
            # histology framed instead of letting the subplot autoscale.
            add_extent_placeholder(fig, samples[si], row, col)
        if not absent:
            continue
        fig.add_annotation(
            text=f"no data here for: {', '.join(absent)}",
            xref="x domain", yref="y domain", x=0.5, y=0.02,
            showarrow=False, font=dict(size=10, color="#5a3d00"),
            bgcolor="rgba(255,244,214,0.92)", bordercolor="#b8860b",
            borderwidth=1, borderpad=3, row=row, col=col,
        )

    # Add tissue images
    for si, s in enumerate(samples):
        b = s["img_bounds"]
        xr, yr = axis_ref(si)
        fig.add_layout_image(dict(
            source=f"/_spatialatlas/tissue/{si}",
            x=b["xmin"], y=b["ymin"],
            sizex=b["xmax"] - b["xmin"], sizey=b["ymax"] - b["ymin"],
            xref=xr, yref=yr,
            sizing="stretch", layer="below"
        ))

    fig.update_yaxes(autorange="reversed", showticklabels=False, showgrid=False)
    fig.update_xaxes(showticklabels=False, showgrid=False)
    for si in range(n_samples):
        row, col = subplot_index(si)
        _, yr = axis_ref(si)
        fig.update_xaxes(scaleanchor=yr, scaleratio=1, row=row, col=col)

    base_height = display_cfg.get("spatial_plot_height", 850)
    plot_height = base_height if N_ROWS == 1 else int(base_height * 0.55 * N_ROWS)
    fig.update_layout(
        template="plotly_white",
        height=plot_height,
        margin=dict(l=10, r=10, t=40, b=10),
        showlegend=False
    )
    return fig


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    app.run(
        debug=server_cfg.get("debug", False),
        host=server_cfg.get("host", "0.0.0.0"),
        port=server_cfg.get("port", 8050),
    )
