#!/usr/bin/env python3
"""Add only the Visium image metadata required by TissUUmaps.

The expression matrix, identifiers, annotations, embeddings and spatial
coordinates are copied without transformation.  Run this script in the pinned
TissUUmaps container, which supplies h5py, NumPy and pyvips.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import h5py
import numpy as np
import pyvips


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_h5ad", type=Path)
    parser.add_argument("histology_png", type=Path)
    parser.add_argument("output_h5ad", type=Path)
    parser.add_argument("--library-id", default="adrenal_hd_8um")
    parser.add_argument("--hires-scale", type=float, default=0.205592)
    return parser.parse_args()


def image_array(path: Path) -> np.ndarray:
    image = pyvips.Image.new_from_file(str(path), access="sequential")
    vips_types = {
        "uchar": np.uint8,
        "char": np.int8,
        "ushort": np.uint16,
        "short": np.int16,
        "uint": np.uint32,
        "int": np.int32,
        "float": np.float32,
        "double": np.float64,
    }
    array = np.ndarray(
        buffer=image.write_to_memory(),
        dtype=vips_types[image.format],
        shape=(image.height, image.width, image.bands),
    )
    return array.copy()


def main() -> None:
    args = parse_args()
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.input_h5ad, args.output_h5ad)
    image = image_array(args.histology_png)

    with h5py.File(args.output_h5ad, "r+") as handle:
        n_obs = int(handle["X"].attrs["shape"][0])
        if "library_id" in handle["obs"]:
            del handle["obs/library_id"]
        library = handle["obs"].create_group("library_id")
        library.attrs["encoding-type"] = "categorical"
        library.attrs["encoding-version"] = "0.2.0"
        library.attrs["ordered"] = False
        library.create_dataset("codes", data=np.zeros(n_obs, dtype=np.int8))
        string_type = h5py.string_dtype(encoding="utf-8")
        library.create_dataset(
            "categories", data=np.asarray([args.library_id], dtype=object), dtype=string_type
        )

        column_order = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in handle["obs"].attrs.get("column-order", [])
        ]
        if "library_id" not in column_order:
            column_order.append("library_id")
        if "column-order" in handle["obs"].attrs:
            del handle["obs"].attrs["column-order"]
        handle["obs"].attrs.create(
            "column-order", np.asarray(column_order, dtype=object), dtype=string_type
        )

        spatial = handle.require_group("uns").require_group("spatial")
        if args.library_id in spatial:
            del spatial[args.library_id]
        library_uns = spatial.create_group(args.library_id)
        images = library_uns.create_group("images")
        images.create_dataset(
            "hires", data=image, compression="gzip", compression_opts=4, shuffle=True
        )
        scalefactors = library_uns.create_group("scalefactors")
        scalefactors.create_dataset(
            "tissue_hires_scalef", data=np.float64(args.hires_scale)
        )

    print(
        f"Wrote {args.output_h5ad}: {n_obs:,} profiles; "
        f"histology {image.shape[1]} x {image.shape[0]}; scale {args.hires_scale}"
    )


if __name__ == "__main__":
    main()
