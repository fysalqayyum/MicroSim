#!/usr/bin/env python3
"""Print quantitative morphology metrics for only the latest stored frame."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import analyze_pilot_morphology as morphology


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--dt", type=float, default=4.0e-7)
    parser.add_argument("--dx", type=float, default=7.5e-7)
    parser.add_argument("--substrate-top-cell", type=int, default=23)
    args = parser.parse_args()

    morphology.DT_S = args.dt
    morphology.DX_M = args.dx
    morphology.SUBSTRATE_TOP_CELL = args.substrate_top_cell
    paths = sorted(args.input_dir.glob("*.vtk"), key=morphology.timestep)
    if not paths:
        raise SystemExit(f"no VTK fragments in {args.input_dir}")
    row, _ = morphology.morphology_row(paths[-1])
    print(json.dumps(row, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
