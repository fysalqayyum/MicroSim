#!/usr/bin/env python3
"""Assert phase-fraction conservation: max |FCC_A1 + LIQUID - 1| over all frames.

This is the regression test MicroSim never had, and the one that would have
caught the addNoise defect in 2023. The upstream kernel perturbed each phase
independently, so every cell it touched drifted off sum(phi) = 1. With noise
off nothing violates it, which is exactly why the bug survived: the check only
bites when the feature under test is actually switched on.

Exits non-zero if any frame exceeds --tolerance.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from analyze_pilot_morphology import read_fragment, timestep


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("frame_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-12,
        help="max allowed |sum(phi) - 1| (default %(default)g)",
    )
    parser.add_argument("--phases", nargs="+", default=["FCC_A1", "LIQUID"])
    args = parser.parse_args()

    frames = sorted(args.frame_dir.glob("*.vtk"), key=timestep)
    if not frames:
        raise SystemExit(f"no .vtk frames in {args.frame_dir}")

    rows = []
    worst = 0.0
    worst_step = -1
    for frame in frames:
        _, fields = read_fragment(frame)
        missing = [p for p in args.phases if p not in fields]
        if missing:
            raise SystemExit(f"{frame.name}: missing phase field(s) {missing}")

        total = sum(fields[p] for p in args.phases)
        deviation = np.abs(total - 1.0)
        # Frames can carry NaN in unfilled regions; nanmax is the project rule.
        max_dev = float(np.nanmax(deviation))
        nonfinite = int(deviation.size - np.count_nonzero(np.isfinite(deviation)))
        step = timestep(frame)
        if max_dev > worst:
            worst, worst_step = max_dev, step
        rows.append(
            {
                "step": step,
                "max_abs_sum_minus_one": max_dev,
                "mean_abs_sum_minus_one": float(np.nanmean(deviation)),
                "nonfinite": nonfinite,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"  step {row['step']:>9}: max|sum(phi)-1| = "
            f"{row['max_abs_sum_minus_one']:.3e}  nonfinite={row['nonfinite']}"
        )
    print(
        f"WORST {worst:.3e} at step {worst_step} "
        f"against tolerance {args.tolerance:.1e}"
    )
    total_nonfinite = sum(r["nonfinite"] for r in rows)
    if total_nonfinite:
        print(f"FAIL: {total_nonfinite} non-finite values in the phase sum")
        raise SystemExit(2)
    if worst > args.tolerance:
        print(
            "FAIL: phase-fraction conservation violated. If this appears only "
            "with Noise_phasefield = 1, the noise kernel is perturbing phases "
            "independently instead of moving fraction between them."
        )
        raise SystemExit(1)
    print("PASS: phase fractions are conserved")


if __name__ == "__main__":
    main()
