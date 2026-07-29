#!/usr/bin/env python3
"""Report field ranges and non-finite counts in MicroSim rank-0 fragments."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np


FIELDS = ("FCC_A1", "LIQUID", "Mu_Si", "Composition_Si", "Temperature")


def timestep(path: Path) -> int:
    match = re.search(r"_(\d+)\.vtk$", path.name)
    if not match:
        raise ValueError(f"cannot extract timestep from {path}")
    return int(match.group(1))


def read_fragment(path: Path) -> dict[str, np.ndarray]:
    with path.open("rb") as handle:
        mesh_x = int(handle.readline())
        mesh_y = int(handle.readline())
        mesh_z = int(handle.readline())
        count = mesh_x * mesh_y * mesh_z
        arrays = {}
        for index, name in enumerate(FIELDS):
            arrays[name] = np.fromfile(handle, dtype=">f8", count=count).astype(
                np.float64
            )
            if arrays[name].size != count:
                raise ValueError(f"{path}: incomplete {name}")
            if index < len(FIELDS) - 1 and handle.read(1) != b"\n":
                raise ValueError(f"{path}: missing separator after {name}")
    return arrays


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    files = sorted(args.root.glob("**/Processor_0/*.vtk"), key=lambda p: (p.parent, timestep(p)))
    if not files:
        raise SystemExit(f"no rank fragments found below {args.root}")

    rows = []
    for path in files:
        row = {"file": str(path.relative_to(args.root)), "step": timestep(path)}
        arrays = read_fragment(path)
        for name, values in arrays.items():
            finite = values[np.isfinite(values)]
            row[f"{name}_nonfinite"] = int(values.size - finite.size)
            row[f"{name}_min"] = float(finite.min()) if finite.size else np.nan
            row[f"{name}_max"] = float(finite.max()) if finite.size else np.nan
        phase_sum = arrays["FCC_A1"] + arrays["LIQUID"]
        finite_sum = phase_sum[np.isfinite(phase_sum)]
        row["phase_sum_nonfinite"] = int(phase_sum.size - finite_sum.size)
        row["phase_sum_min"] = (
            float(finite_sum.min()) if finite_sum.size else np.nan
        )
        row["phase_sum_max"] = (
            float(finite_sum.max()) if finite_sum.size else np.nan
        )
        row["phase_sum_max_abs_error"] = (
            float(np.max(np.abs(finite_sum - 1.0)))
            if finite_sum.size
            else np.nan
        )
        rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} frame diagnostics to {args.output}")


if __name__ == "__main__":
    main()
