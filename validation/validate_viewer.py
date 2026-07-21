#!/usr/bin/env python3
"""Import SpatialAtlas and exercise all visualization callbacks directly."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--app", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    os.environ["SPATIALATLAS_CONFIG"] = str(args.config.resolve())

    started = time.perf_counter()
    spec = importlib.util.spec_from_file_location("spatialatlas_app_under_test", args.app.resolve())
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    import_seconds = time.perf_counter() - started

    gene = next((g for g in module.default_genes if g in module.gene_set), module.gene_names[0])
    checks = {}

    def exercise(name, function, *call_args):
        before = time.perf_counter()
        figure = function(*call_args)
        elapsed = time.perf_counter() - before
        payload = figure.to_json()
        assert len(figure.data) > 0, f"{name} returned no traces"
        checks[name] = {
            "status": "PASS",
            "seconds": elapsed,
            "traces": len(figure.data),
            "json_bytes": len(payload.encode("utf-8")),
        }

    exercise("spatial_cluster", module.update_spatial, None, 2.5, 0.6, "Transparent-Green")
    exercise("spatial_gene", module.update_spatial, gene, 2.5, 0.6, "Transparent-Green")
    if module.has_umap:
        exercise("umap", module.update_umap, gene)
    exercise("violin", module.update_violin, gene)
    if module.markers_df is not None:
        exercise("heatmap", module.update_heatmap, 5)
    exercise(
        "multi_gene", module.update_multi,
        *module.default_genes,
        *(["on"] for _ in range(4)),
        *([0.8] * 4),
        2.5,
    )

    response = module.server.test_client().get("/")
    assert response.status_code == 200
    image_response = module.server.test_client().get("/_spatialatlas/tissue/0")
    assert image_response.status_code == 200
    assert image_response.headers.get("Cache-Control")
    checks["http"] = {
        "status": "PASS",
        "root_status": response.status_code,
        "image_status": image_response.status_code,
        "image_cache_control": image_response.headers.get("Cache-Control"),
    }

    report = {
        "status": "PASS",
        "config": str(args.config),
        "import_seconds": import_seconds,
        "n_obs": module.n_obs_combined,
        "n_genes": len(module.gene_names),
        "n_samples": module.n_samples,
        "expression_cache_gib": module._total_sparse_bytes / (1024 ** 3),
        "gene": gene,
        "checks": checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
