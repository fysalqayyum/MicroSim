#!/usr/bin/env python3
"""Measure dendrite/cell tip undercooling and time-to-steady from frame series.

WHY: it is tempting to size a run from the ASSUMED solute transient 20*D/V^2
and from a 2-D Ivantsov tip undercooling. Both propagate into every downstream
decision -- domain height and wall-clock cost, the latter scaling as 1/V^2 --
and both were wrong in practice. 20*D/V^2 in particular is falsified: fit
v_tip(t) = V(1 - B exp(-t/tau)) with the asymptote PINNED to V instead. This
script replaces the assumptions with measurements.

METHOD: Temperature is a stored field in every MicroSim frame, so the tip
temperature is read directly and needs NO moving-window shift bookkeeping --
the shift cancels because we sample T at the tip cell itself. Undercooling is

    dT_tip = T_liquidus(c0) - T(tip)

and steady state is reached when dT_tip stops drifting. That criterion is also
shift-free, which is why it is preferred here over a lab-frame tip velocity.

NaN handling: subcell_interface_height returns NaN for columns with no phi=0.5
crossing (deep grooves between cells). Bare argmax/max return the NaN index, so
every reduction here is the nan-safe variant.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from analyze_pilot_morphology import read_fragment, subcell_interface_height, timestep


def tip_state(
    phase: np.ndarray, temperature: np.ndarray, composition: np.ndarray
) -> dict[str, float]:
    """Locate the leading interface point and sample fields there."""
    height = subcell_interface_height(phase)
    if not np.any(np.isfinite(height)):
        return {key: float("nan") for key in (
            "tip_x", "tip_y_cell", "mean_y_cell", "amplitude_cells",
            "T_tip_K", "c_ahead_of_tip", "finite_columns",
        )}

    tip_x = int(np.nanargmax(height))
    tip_y = float(height[tip_x])

    # Linear interpolation of T between the two cells bracketing the phi=0.5
    # crossing. T varies linearly in y anyway, so this is exact up to roundoff.
    ny = temperature.shape[1]
    lower = int(math.floor(tip_y))
    frac = tip_y - lower
    if lower + 1 < ny:
        t_tip = float(
            temperature[tip_x, lower] * (1.0 - frac)
            + temperature[tip_x, lower + 1] * frac
        )
    else:
        t_tip = float(temperature[tip_x, lower])

    # Liquid composition just ahead of the tip - tests the Ivantsov c_l,tip
    # prediction. Sampled 2 cells ahead to clear the diffuse interface.
    probe = min(lower + 2, ny - 1)
    c_ahead = float(composition[tip_x, probe])

    return {
        "tip_x": float(tip_x),
        "tip_y_cell": tip_y,
        "mean_y_cell": float(np.nanmean(height)),
        "amplitude_cells": tip_y - float(np.nanmean(height)),
        "T_tip_K": t_tip,
        "c_ahead_of_tip": c_ahead,
        "finite_columns": float(np.count_nonzero(np.isfinite(height))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path, help="DATA/Processor_0 of one case")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--dx", type=float, required=True)
    parser.add_argument("--velocity", type=float, required=True, help="imposed V, m/s")
    parser.add_argument("--diffusivity", type=float, default=6.25e-9)
    parser.add_argument("--t-liquidus", type=float, default=899.604509)
    parser.add_argument("--label", default="")
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("*.vtk"), key=timestep)
    if not paths:
        raise SystemExit(f"no VTK fragments in {args.input_dir}")

    rows: list[dict[str, float | int | str]] = []
    for path in paths:
        row: dict[str, float | int | str] = {
            "file": path.name,
            "step": timestep(path),
            "time_s": timestep(path) * args.dt,
        }
        try:
            dimensions, fields = read_fragment(path)
        except ValueError as error:
            # read_fragment is also the finiteness gate; record and continue so
            # one bad frame cannot hide the rest of the series.
            row["error"] = str(error).replace(",", ";")
            rows.append(row)
            continue

        state = tip_state(
            fields["FCC_A1"][:, :, 0],
            fields["Temperature"][:, :, 0],
            fields["Composition_Si"][:, :, 0],
        )
        row.update(state)
        row["undercooling_K"] = args.t_liquidus - state["T_tip_K"]
        row["tip_y_um"] = state["tip_y_cell"] * args.dx * 1e6
        row["amplitude_um"] = state["amplitude_cells"] * args.dx * 1e6
        row["phi_max"] = float(np.max(fields["FCC_A1"][:, :, 0]))
        row["nx"] = int(dimensions[0])
        row["ny"] = int(dimensions[1])
        rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    # lineterminator="\n" is mandatory in this project: Python's default CRLF
    # once killed a job through a Slurm string assert.
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # ---- summary against the assumptions being tested -----------------------
    steps = np.array([r["step"] for r in rows if "undercooling_K" in r], dtype=float)
    times = np.array([r["time_s"] for r in rows if "undercooling_K" in r], dtype=float)
    dT = np.array([r["undercooling_K"] for r in rows if "undercooling_K" in r], dtype=float)
    good = np.isfinite(dT)
    steps, times, dT = steps[good], times[good], dT[good]

    tau_assumed = 20.0 * args.diffusivity / args.velocity**2

    print(f"===== {args.label or args.input_dir} =====")
    print(f"frames analysed        = {dT.size}")
    print(f"assumed transient 20D/V^2 = {tau_assumed:.4e} s "
          f"({tau_assumed/args.dt:.0f} steps)")
    if dT.size == 0:
        print("no finite tip undercooling recovered")
        return
    print(f"latest time reached    = {times[-1]:.4e} s "
          f"({times[-1]/tau_assumed:.2f} x the assumed transient)")
    print(f"undercooling first/last= {dT[0]:.2f} / {dT[-1]:.2f} K")
    print(f"undercooling min/max   = {np.min(dT):.2f} / {np.max(dT):.2f} K")

    # Drift over the last three intervals is the steady-state test.
    if dT.size >= 4:
        recent = dT[-4:]
        recent_t = times[-4:]
        drift = float(np.polyfit(recent_t, recent, 1)[0])
        span = float(np.max(recent) - np.min(recent))
        print(f"drift over last 4 frames = {drift:+.1f} K/s "
              f"(spread {span:.3f} K over {recent_t[-1]-recent_t[0]:.3e} s)")
        print("VERDICT: " + (
            "undercooling has PLATEAUED (spread < 0.5 K)" if span < 0.5
            else "still TRANSIENT - undercooling is still moving"
        ))
    print()


if __name__ == "__main__":
    main()
