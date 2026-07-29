#!/usr/bin/env python3
"""Compare two MicroSim VTK fragments field by field."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from analyze_pilot_morphology import read_fragment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    _, reference = read_fragment(args.reference)
    _, candidate = read_fragment(args.candidate)
    if reference.keys() != candidate.keys():
        raise ValueError(
            f"field mismatch: {sorted(reference)} != {sorted(candidate)}"
        )

    rows: list[dict[str, float | int | str]] = []
    for name in reference:
        if reference[name].shape != candidate[name].shape:
            raise ValueError(
                f"{name}: shape mismatch "
                f"{reference[name].shape} != {candidate[name].shape}"
            )
        difference = candidate[name] - reference[name]
        rows.append(
            {
                "field": name,
                "values": difference.size,
                "nonfinite_difference": int(
                    difference.size - np.count_nonzero(np.isfinite(difference))
                ),
                "maximum_absolute_difference": float(np.nanmax(np.abs(difference))),
                "rms_difference": float(np.sqrt(np.nanmean(difference**2))),
                "reference_l2": float(np.sqrt(np.mean(reference[name] ** 2))),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} field comparisons to {args.output}")


if __name__ == "__main__":
    main()
