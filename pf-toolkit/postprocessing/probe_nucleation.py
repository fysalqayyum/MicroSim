#!/usr/bin/env python3
"""Undercooling-versus-height probe for directional-solidification frames.

Answers one question: is the liquid AHEAD of the front undercooled, and by how
much? Prints, for each requested height, the mean temperature, mean composition,
the liquidus of that local composition, and the resulting undercooling.

Typical uses
------------
* Deciding whether a morphology burst was spurious bulk nucleation or a genuine
  front instability -- probe the last clean frame and the first bad one.
* Sizing the liquid above the front. Note that in the KKS model implemented here
  bulk liquid at phi = 0 is INERT (h'(0) = h''(0) = 0, and the noise is
  interface-localised), so deep far-field undercooling is not by itself a fault.
  See pf-toolkit/LESSONS.md.

The liquidus is the linearised binary form
    T_liq(c) = T_liq(c0) + m_L * (wt%(c) - wt%(c0))
so --t-liq-c0, --m-l and --c0 must match the alloy. Defaults are Al-5.5wt%Si;
--solvent-mass / --solute-mass default to Al / Si.

Example
-------
    python probe_nucleation.py frame_1350000.vtk frame_1400000.vtk \\
        --heights 600 900 1200 1650 2000 --dx 1.0771e-8
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_pilot_morphology import read_fragment  # noqa: E402


def wt_solute(x, m_solute, m_solvent):
    """wt% solute of a binary liquid at mole fraction x."""
    return 100.0 * x * m_solute / (x * m_solute + (1.0 - x) * m_solvent)


def parse_args():
    p = argparse.ArgumentParser(
        description="Probe undercooling versus height in solidification frames.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("frames", nargs="+", type=Path, help="VTK frame(s) to probe")
    p.add_argument("--heights", type=int, nargs="+", default=None,
                   help="grid rows to report (default: 9 evenly spaced rows)")
    p.add_argument("--dx", type=float, default=1.0771e-8, help="cell size, m")
    p.add_argument("--c0", type=float, default=0.052953,
                   help="nominal composition, mole fraction solute")
    p.add_argument("--t-liq-c0", type=float, default=899.604509,
                   help="liquidus temperature at c0, K")
    p.add_argument("--m-l", type=float, default=-6.507,
                   help="liquidus slope, K per wt%% solute")
    p.add_argument("--solute-mass", type=float, default=28.0855,
                   help="solute molar mass, g/mol (default Si)")
    p.add_argument("--solvent-mass", type=float, default=26.9815,
                   help="solvent molar mass, g/mol (default Al)")
    p.add_argument("--solid-phase", default="FCC_A1",
                   help="VTK array name of the solid phase field")
    p.add_argument("--composition-array", default="Composition_Si",
                   help="VTK array name of the composition field")
    p.add_argument("--phi-threshold", type=float, default=0.5,
                   help="phi above which a cell counts as solid")
    return p.parse_args()


def main():
    args = parse_args()
    wt0 = wt_solute(args.c0, args.solute_mass, args.solvent_mass)

    for frame in args.frames:
        dims, arrays = read_fragment(frame)
        print(f"\n===== {frame.name}  dims={dims} =====")
        print("arrays:", sorted(arrays.keys()))

        if args.solid_phase not in arrays:
            print(f"  MISSING array '{args.solid_phase}' -- skipping frame")
            continue

        phi = arrays[args.solid_phase][:, :, 0]          # [x, y]
        xs = arrays[args.composition_array][:, :, 0]

        solid = phi >= args.phi_threshold
        # topmost solid row per column, -1 where the column is all liquid
        top = np.where(solid.any(axis=1),
                       solid.shape[1] - 1 - np.argmax(solid[:, ::-1], axis=1), -1)
        print(f"topmost solid cell: max={top.max()}  median={int(np.median(top))}")

        ny = phi.shape[1]
        heights = args.heights or list(np.linspace(0, ny - 1, 9, dtype=int))

        if "Temperature" not in arrays:
            print("  NO Temperature array in this VTK -- cannot evaluate undercooling.")
            print("  (Write the temperature field, or reconstruct it from the "
                  "gradient and pulling velocity in the input file.)")
            continue

        temperature = arrays["Temperature"][:, :, 0]
        print("  row     height       T (K)      x_solute    wt%     "
              "T_liq (K)   undercooling   phi_mean")
        for y in heights:
            if not 0 <= y < ny:
                continue
            t_mean = temperature[:, y].mean()
            x_mean = xs[:, y].mean()
            wt = wt_solute(x_mean, args.solute_mass, args.solvent_mass)
            t_liq = args.t_liq_c0 + args.m_l * (wt - wt0)
            print(f"  {y:5d} {y * args.dx * 1e6:8.2f} um {t_mean:10.3f} "
                  f"{x_mean:11.6f} {wt:7.3f} {t_liq:11.3f} "
                  f"{t_liq - t_mean:+13.3f} K {phi[:, y].mean():9.4f}")


if __name__ == "__main__":
    main()
