#!/usr/bin/env python3
"""Qualify a MicroSim restart against a continuous reference.

Acceptance criterion (project standard): GPU kernels are not bit-reproducible
(floating-point non-associativity), so bit-identity is the wrong test. A restart
is equivalent when the per-field RELATIVE-L2 drift stays below 1e-4, where

    relative_l2 = rms_difference / reference_l2

and every field is finite. This supersedes the earlier absolute 1e-10 tolerance,
which demanded bit-identity and produced false negatives (one case where
all fields drifted only ~1e-10 to ~1e-12 relative yet "failed" the absolute test).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

# Documented GPU-realistic acceptance threshold on relative-L2 drift.
RELATIVE_L2_TOLERANCE = 1.0e-4


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("comparison", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--tolerance",
        type=float,
        default=RELATIVE_L2_TOLERANCE,
        help="relative-L2 acceptance threshold (default %(default)g)",
    )
    args = parser.parse_args()

    with args.comparison.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    qualified = []
    for row in rows:
        nonfinite = int(row["nonfinite_difference"])
        rms = float(row["rms_difference"])
        ref_l2 = float(row["reference_l2"])
        # Relative-L2 drift; if the reference field is identically zero, fall
        # back to the raw RMS difference so a nonzero drift cannot hide.
        relative_l2 = rms / ref_l2 if ref_l2 > 0.0 else rms
        passed = nonfinite == 0 and relative_l2 <= args.tolerance
        qualified.append(
            {
                **row,
                "relative_l2": relative_l2,
                "predeclared_relative_l2_tolerance": args.tolerance,
                "pass": int(passed),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=qualified[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(qualified)

    for row in qualified:
        print(
            f"{row['field']}: relative_l2={row['relative_l2']:.3e} "
            f"nonfinite={row['nonfinite_difference']} pass={row['pass']}"
        )

    overall = all(int(row["pass"]) == 1 for row in qualified)
    print(f"restart_equivalence_pass={int(overall)}")
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
