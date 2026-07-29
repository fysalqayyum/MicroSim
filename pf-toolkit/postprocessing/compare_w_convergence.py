#!/usr/bin/env python3
"""Interface-width convergence: is the answer insensitive to W?

Compares the base and refined runs at MATCHED PHYSICAL TIME. The two runs use
different dt and different SAVET chosen so that frame k of each sits at the same
physical time (base 50000*2.0e-9 = refined 112500*8.8889e-10 = 1.0e-4 s), so the
comparison is frame-index against frame-index after checking the times agree.

ACCEPT: < 5% change in tip radius AND tip undercooling.

Exits 0 if converged, 1 if not, 2 if the inputs are unusable. The exit code is
the point -- this is a publication gate, not a report.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def load(case_dir: Path) -> tuple[dict, list[dict], list[dict]]:
    manifest: dict[str, str] = {}
    man_path = case_dir / "case_manifest.txt"
    if man_path.exists():
        for line in man_path.read_text().splitlines():
            if "=" in line:
                key, _, value = line.partition("=")
                manifest[key.strip()] = value.strip()

    def read(name: str) -> list[dict]:
        path = case_dir / name
        if not path.exists():
            raise SystemExit(f"{path} missing - has {case_dir.name} finished?")
        with path.open() as handle:
            return list(csv.DictReader(handle))

    return manifest, read("tip_undercooling.csv"), read("tip_radius.csv")


def to_float(row: dict, key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=Path)
    parser.add_argument("refined_dir", type=Path)
    parser.add_argument("--tolerance", type=float, default=5.0, help="percent")
    args = parser.parse_args()

    b_man, b_dT, b_rad = load(args.base_dir)
    r_man, r_dT, r_rad = load(args.refined_dir)

    print(f"base    : {b_man.get('mesh','?')}  W/dx={b_man.get('W_over_dx','?')}  "
          f"W/d0={b_man.get('W_over_d0','?')}  dt={b_man.get('delta_t','?')}")
    print(f"refined : {r_man.get('mesh','?')}  W/dx={r_man.get('W_over_dx','?')}  "
          f"W/d0={r_man.get('W_over_d0','?')}  dt={r_man.get('delta_t','?')}")

    # The whole comparison rests on the frame times matching. Verify, don't assume.
    b_int = float(b_man.get("frame_interval_s", "nan"))
    r_int = float(r_man.get("frame_interval_s", "nan"))
    if not (abs(b_int - r_int) / b_int < 1e-3):
        raise SystemExit(
            f"frame intervals differ ({b_int:.6e} vs {r_int:.6e} s) - "
            "the runs are not comparable frame-by-frame")
    print(f"frame interval matched: {b_int:.6e} s\n")

    n = min(len(b_dT), len(r_dT), len(b_rad), len(r_rad))
    if n < 2:
        raise SystemExit("not enough frames to compare")

    print(f"{'frame':>6} {'t (s)':>11} {'dT_base':>9} {'dT_ref':>9} {'dT %':>8} "
          f"{'rad_base':>9} {'rad_ref':>9} {'rad %':>8}")
    print("-" * 78)

    worst_dT = worst_rad = 0.0
    n_dT = n_rad = 0          # valid (non-NaN) comparisons actually made
    for i in range(n):
        t_b = to_float(b_dT[i], "time_s")
        dT_b, dT_r = to_float(b_dT[i], "undercooling_K"), to_float(r_dT[i], "undercooling_K")
        rb = to_float(b_rad[i], "median_radius_um")
        rr = to_float(r_rad[i], "median_radius_um")
        d_dT = 100.0 * (dT_r - dT_b) / dT_b if dT_b else float("nan")
        d_rad = 100.0 * (rr - rb) / rb if rb else float("nan")
        # Skip the first frames: the seed has not yet grown a resolvable tip and
        # the percentage is dominated by the initial condition, not by W.
        if i >= 2:
            if d_dT == d_dT:
                worst_dT = max(worst_dT, abs(d_dT))
                n_dT += 1
            if d_rad == d_rad:
                worst_rad = max(worst_rad, abs(d_rad))
                n_rad += 1
        print(f"{i:>6} {t_b:>11.4e} {dT_b:>9.3f} {dT_r:>9.3f} {d_dT:>+8.2f} "
              f"{rb:>9.4f} {rr:>9.4f} {d_rad:>+8.2f}")

    print()
    print(f"worst |change| after frame 2:  undercooling {worst_dT:.2f}% "
          f"({n_dT} valid frames)   tip radius {worst_rad:.2f}% "
          f"({n_rad} valid frames)   (tolerance {args.tolerance:.1f}%)")

    # A metric with no valid comparisons must NOT pass. Before this guard,
    # worst_rad stayed at its 0.0 initialiser when every frame was NaN, and
    # 0.0 < tolerance reported CONVERGED on zero evidence. The tip-radius
    # extractor returns NaN once the doublon forms (a split tip defeats the
    # parabola fit), so this was the actual behaviour, not a hypothetical.
    MIN_VALID = 3
    insufficient = []
    if n_dT < MIN_VALID:
        insufficient.append(f"undercooling ({n_dT} valid, need {MIN_VALID})")
    if n_rad < MIN_VALID:
        insufficient.append(f"tip radius ({n_rad} valid, need {MIN_VALID})")
    if insufficient:
        print("VERDICT: INSUFFICIENT DATA - cannot certify convergence for: "
              + "; ".join(insufficient))
        print("  This is NOT a pass. Fix the extractor or supply more frames.")
        sys.exit(3)

    converged = worst_dT < args.tolerance and worst_rad < args.tolerance
    if converged:
        print("VERDICT: CONVERGED - results are insensitive to W at this resolution.")
        sys.exit(0)
    print("VERDICT: NOT CONVERGED - W is still affecting the answer.")
    print("  Production results at W = 6.4623e-8 are not publication-defensible")
    print("  until this passes. Next step is a third, finer level (W/2.25).")
    sys.exit(1)


if __name__ == "__main__":
    main()
