#!/usr/bin/env python3
"""Estimate seeded dendrite-tip radii from local quadratic interface fits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from analyze_pilot_morphology import read_fragment, subcell_interface_height, timestep


def fit_tip_radius(
    height: np.ndarray, center: int, half_window: int, fit_depth: float
) -> tuple[float, float]:
    nx = height.size
    offsets = np.arange(-half_window, half_window + 1)
    indices = (center + offsets) % nx
    local_height = height[indices]
    # Columns with no phi=0.5 crossing (deep intergranular liquid grooves) come
    # back as NaN from subcell_interface_height. np.argmax reports the NaN index
    # rather than the tip, so a single groove inside the window poisoned the whole
    # seed; nanargmax picks the true crest instead. The fit mask below already
    # drops NaNs, since any comparison against NaN is False.
    if not np.any(np.isfinite(local_height)):
        return float("nan"), float("nan")
    peak_index = int(np.nanargmax(local_height))
    peak_offset = float(offsets[peak_index])
    peak_height = float(local_height[peak_index])
    x = offsets.astype(float) - peak_offset
    mask = local_height >= peak_height - fit_depth
    if np.count_nonzero(mask) < 5:
        return float("nan"), peak_height
    a, _, _ = np.polyfit(x[mask], local_height[mask], 2)
    if not np.isfinite(a) or a >= 0:
        return float("nan"), peak_height
    return float(1.0 / abs(2.0 * a)), peak_height


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--centers", default="43,128,213")
    parser.add_argument("--half-window", type=int, default=30)
    parser.add_argument("--fit-depth-cells", type=float, default=6.0)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--dx", type=float, required=True)
    args = parser.parse_args()
    centers = [int(value) for value in args.centers.split(",")]

    paths = sorted(args.input_dir.glob("*.vtk"), key=timestep)
    if not paths:
        raise SystemExit(f"no VTK fragments in {args.input_dir}")

    rows: list[dict[str, float | int | str]] = []
    for path in paths:
        _, fields = read_fragment(path)
        height = subcell_interface_height(fields["FCC_A1"][:, :, 0])
        radii: list[float] = []
        tips: list[float] = []
        row: dict[str, float | int | str] = {
            "file": path.name,
            "step": timestep(path),
            "time_s": timestep(path) * args.dt,
        }
        for number, center in enumerate(centers, start=1):
            radius, tip = fit_tip_radius(
                height, center, args.half_window, args.fit_depth_cells
            )
            radii.append(radius)
            tips.append(tip)
            row[f"seed{number}_tip_cell"] = tip
            row[f"seed{number}_radius_cells"] = radius
            row[f"seed{number}_radius_um"] = radius * args.dx * 1e6
        finite_radii = np.asarray(radii)[np.isfinite(radii)]
        row["median_radius_cells"] = (
            float(np.median(finite_radii)) if finite_radii.size else float("nan")
        )
        row["median_radius_um"] = (
            float(row["median_radius_cells"]) * args.dx * 1e6
        )
        row["minimum_radius_over_dx"] = (
            float(np.min(finite_radii)) if finite_radii.size else float("nan")
        )
        row["mean_tip_cell"] = float(np.mean(tips))
        rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} tip-radius rows to {args.output}")


if __name__ == "__main__":
    main()
