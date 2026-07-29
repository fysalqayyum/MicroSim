#!/usr/bin/env python3
"""Generate a broadband (multi-mode) FCC-Al substrate for Mullins-Sekerka selection.

WHY THIS EXISTS
---------------
The solver's own `Noise_phasefield` cannot seed a planar-front instability. In
`CL_Kim_Kernel.cl:6875` the phase-field noise increment is multiplied by the
phase index `ip`:

    phi[ip] -= (NoiseFac*noise) * phi[ip]*(1-phi[ip]) * ip / (npha*npha);

so for `ip = 0` -- which is FCC_A1, the solid, when NUMPHASES = 2 -- the
increment is identically zero. Only LIQUID is ever perturbed, which both leaves
the solid interface bit-exactly flat and breaks sum(phi) = 1. Separately, the
kernel's "random" value is a pure function of cell indices with no timestep
input and no carried state, so it is a static spatial bias rather than a noise
process. A dedicated gate run measured the consequence: interface amplitude
0.000000 -> -0.000000 um over 20000 steps.

`generate_sinusoidal_filling.py` is the wrong substitute: it seeds a SINGLE
mode, which re-imposes the primary spacing as an input. That is precisely the
defect that invalidated an early campaign (W/lambda1 = 1.00-1.03, spacing was an
input, not a result).

WHAT THIS DOES
--------------
Superposes `--n-modes` sinusoids with equal amplitude and pseudo-random phases,
normalised to a target RMS in cells, then rounds to integer cell heights (the
FILLCUBE grammar addresses cells, so sub-cell amplitudes cannot be expressed).
Every mode from 1 to n_modes is present at comparable amplitude, so the front
is free to select its own fastest-growing wavelength.

DETERMINISTIC BY CONSTRUCTION: phases come from a fixed `--seed`, so the run is
reproducible and the seed is a reportable input. This is a virtue over solver
RNG noise for a publication, not a compromise.

Choose `--n-modes` so the shortest seeded wavelength stays well resolved: at
nx = 3168, dx = 1.0771e-8 m, mode 64 has wavelength 5.33e-7 m = 8.2 W, which is
comfortably above the interface width. The script refuses a mode count whose
shortest wavelength falls under `--min-wavelength-over-w` interface widths when
`--dx` and `--epsilon` are supplied.
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--base-row", type=int, required=True)
    parser.add_argument("--amplitude-cells", type=float, default=1.0,
                        help="target RMS of the perturbation, in cells")
    parser.add_argument("--n-modes", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--dx", type=float, default=None)
    parser.add_argument("--epsilon", type=float, default=None)
    parser.add_argument("--min-wavelength-over-w", type=float, default=4.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.nx < 2:
        raise SystemExit("nx must be >= 2")
    if args.n_modes < 2:
        raise SystemExit("n_modes must be >= 2; a single mode imposes the spacing")
    if args.n_modes > args.nx // 2:
        raise SystemExit(f"n_modes {args.n_modes} exceeds Nyquist {args.nx // 2}")
    if args.base_row - math.ceil(3.0 * args.amplitude_cells) < 1:
        raise SystemExit("base_row too low: the perturbation would reach row 0")

    # Refuse a spectrum that reaches below the resolvable interface scale.
    if args.dx is not None and args.epsilon is not None:
        shortest = args.nx * args.dx / args.n_modes
        ratio = shortest / args.epsilon
        if ratio < args.min_wavelength_over_w:
            raise SystemExit(
                f"shortest seeded wavelength {shortest:.4e} m is only {ratio:.2f} W; "
                f"reduce --n-modes below {int(args.nx * args.dx / (args.min_wavelength_over_w * args.epsilon))}"
            )

    rng = random.Random(args.seed)
    phases = [rng.uniform(0.0, 2.0 * math.pi) for _ in range(args.n_modes)]

    # Equal-amplitude modes with random phase: the sum of n_modes unit cosines
    # has RMS sqrt(n_modes/2), so divide by that to hit the requested RMS.
    norm = args.amplitude_cells / math.sqrt(args.n_modes / 2.0)

    raw = []
    for x in range(args.nx):
        s = sum(
            math.cos(2.0 * math.pi * (m + 1) * x / args.nx + phases[m])
            for m in range(args.n_modes)
        )
        raw.append(norm * s)

    heights = [args.base_row + round(v) for v in raw]

    achieved_rms = math.sqrt(sum(v * v for v in raw) / len(raw))
    quantised_rms = math.sqrt(
        sum((h - args.base_row) ** 2 for h in heights) / len(heights)
    )

    header = [
        "## Generated broadband (multi-mode) FCC-Al substrate.",
        f"## nx={args.nx}; base_row={args.base_row}; n_modes={args.n_modes}; seed={args.seed}",
        f"## target_rms_cells={args.amplitude_cells}; continuous_rms_cells={achieved_rms:.4f};"
        f" quantised_rms_cells={quantised_rms:.4f}",
        f"## height_min={min(heights)}; height_max={max(heights)}",
    ]
    if args.dx is not None:
        header.append(
            f"## shortest_wavelength_m={args.nx * args.dx / args.n_modes:.4e};"
            f" longest_wavelength_m={args.nx * args.dx:.4e}"
        )

    lines = header + [
        f"FILLCUBE = {{0,{x},0,0,{x},{h},0}};" for x, h in enumerate(heights)
    ]
    args.output.write_text("\n".join(lines) + "\n")

    print(f"wrote {args.output}")
    print(f"  modes           = 1..{args.n_modes} (seed {args.seed})")
    print(f"  rms cells       = {quantised_rms:.4f} after integer quantisation")
    print(f"  height range    = {min(heights)} .. {max(heights)} (base {args.base_row})")
    if quantised_rms == 0.0:
        raise SystemExit(
            "quantised RMS is zero - the amplitude rounded away to a flat front; "
            "raise --amplitude-cells"
        )


if __name__ == "__main__":
    main()
