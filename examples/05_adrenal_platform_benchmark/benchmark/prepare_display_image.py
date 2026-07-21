#!/usr/bin/env python3
"""Create and validate the browser display image used by the HD benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def psnr(reference: Image.Image, candidate: Image.Image) -> float:
    squared_error = 0
    values = 0
    width, height = reference.size
    for y0 in range(0, height, 256):
        box = (0, y0, width, min(y0 + 256, height))
        left = np.asarray(reference.crop(box), dtype=np.int16)
        right = np.asarray(candidate.crop(box), dtype=np.int16)
        difference = left.astype(np.int32) - right.astype(np.int32)
        squared_error += int(np.square(difference, dtype=np.int64).sum())
        values += difference.size
    mean_squared_error = squared_error / values
    return 20 * math.log10(255.0) - 10 * math.log10(mean_squared_error)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("validation", type=Path)
    parser.add_argument("--quality", type=int, default=90)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(args.source) as image:
        reference = image.convert("RGB")
        reference.save(
            args.output,
            format="WEBP",
            quality=args.quality,
            method=6,
            optimize=True,
        )

    with Image.open(args.output) as image:
        candidate = image.convert("RGB")
        if candidate.size != reference.size:
            raise ValueError(
                f"Display image dimensions {candidate.size} differ from {reference.size}"
            )
        measured_psnr = psnr(reference, candidate)

    result = {
        "source": str(args.source),
        "display_image": str(args.output),
        "dimensions": list(reference.size),
        "source_bytes": args.source.stat().st_size,
        "display_bytes": args.output.stat().st_size,
        "size_reduction_percent": 100 * (
            1 - args.output.stat().st_size / args.source.stat().st_size
        ),
        "webp_quality": args.quality,
        "psnr_db": measured_psnr,
        "source_sha256": sha256(args.source),
        "display_sha256": sha256(args.output),
    }
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    args.validation.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
