#!/usr/bin/env python3
"""Validate that the TissUUmaps adapter and CSC copy retain the HD data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


def digest(dataset: h5py.Dataset) -> str:
    hasher = hashlib.sha256()
    for start in range(0, dataset.shape[0], 1_000_000):
        values = dataset[start : start + 1_000_000]
        if values.dtype.kind in "OSU":
            for value in values.reshape(-1):
                encoded = value if isinstance(value, bytes) else str(value).encode("utf-8")
                hasher.update(len(encoded).to_bytes(8, "little"))
                hasher.update(encoded)
        else:
            hasher.update(np.ascontiguousarray(values).tobytes())
    return hasher.hexdigest()


def sparse_matrix(handle: h5py.File) -> sparse.spmatrix:
    group = handle["X"]
    arrays = (group["data"][:], group["indices"][:], group["indptr"][:])
    shape = tuple(int(value) for value in group.attrs["shape"])
    encoding = group.attrs["encoding-type"]
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    return sparse.csr_matrix(arrays, shape=shape) if encoding == "csr_matrix" else sparse.csc_matrix(arrays, shape=shape)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("original", type=Path)
    parser.add_argument("adapter", type=Path)
    parser.add_argument("converted", type=Path)
    parser.add_argument("output_json", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with h5py.File(args.original, "r") as original, h5py.File(args.adapter, "r") as adapter, h5py.File(args.converted, "r") as converted:
        exact_datasets = ["X/data", "X/indices", "X/indptr", "obs/_index", "var/_index", "obsm/spatial"]
        exact = {
            name: {
                "original_sha256": digest(original[name]),
                "adapter_sha256": digest(adapter[name]),
                "identical": digest(original[name]) == digest(adapter[name]),
            }
            for name in exact_datasets
        }
        source_matrix = sparse_matrix(original)
        converted_matrix = sparse_matrix(converted)
        difference = source_matrix - converted_matrix
        matrix_identical = difference.nnz == 0
        payload = {
            "original": str(args.original),
            "adapter": str(args.adapter),
            "tissuumaps_converted": str(args.converted),
            "shape": list(source_matrix.shape),
            "nonzero_values": int(source_matrix.nnz),
            "adapter_exact_dataset_checks": exact,
            "all_adapter_checks_passed": all(item["identical"] for item in exact.values()),
            "converted_encoding": str(converted["X"].attrs["encoding-type"]),
            "converted_matrix_numerically_identical": matrix_identical,
            "converted_difference_nonzero": int(difference.nnz),
            "converted_spatial_coordinates_identical": digest(original["obsm/spatial"]) == digest(converted["obsm/spatial"]),
            "converted_observation_ids_identical": digest(original["obs/_index"]) == digest(converted["obs/_index"]),
            "converted_gene_ids_identical": digest(original["var/_index"]) == digest(converted["var/_index"]),
        }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if not (payload["all_adapter_checks_passed"] and matrix_identical):
        raise SystemExit("Equivalence validation failed")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
