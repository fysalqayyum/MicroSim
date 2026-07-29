#!/usr/bin/env python3
"""Fail-fast validation of the Al–Si MicroSim input package."""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MW_AL = 26.9815385
MW_SI = 28.0855


def scalar(text: str, key: str) -> float:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([^;]+);", text)
    if not match:
        raise AssertionError(f"missing {key}")
    return float(match.group(1).strip())


def tuple_values(text: str, key: str) -> list[list[float]]:
    matches = re.findall(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*\{{([^}}]+)\}}\s*;", text
    )
    return [[float(v.strip()) for v in row.split(",")] for row in matches]


def read_compositions() -> list[dict[str, float]]:
    path = ROOT / "tdbs_encrypted" / "Composition_FCC_A1.csv"
    with path.open(newline="") as handle:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(handle)]


def interp(rows: list[dict[str, float]], xkey: str, ykey: str, x: float) -> float:
    rows = sorted(rows, key=lambda row: row[xkey])
    for a, b in zip(rows, rows[1:]):
        if a[xkey] <= x <= b[xkey]:
            f = (x - a[xkey]) / (b[xkey] - a[xkey])
            return a[ykey] + f * (b[ykey] - a[ykey])
    raise AssertionError(f"{x} outside {xkey} table")


def validate(path: Path, moving: bool) -> None:
    text = path.read_text()
    nx = int(scalar(text, "MESH_X"))
    ny = int(scalar(text, "MESH_Y"))
    dx = scalar(text, "DELTA_X")
    dy = scalar(text, "DELTA_Y")
    dt = scalar(text, "DELTA_t")
    epsilon = scalar(text, "epsilon")
    nsteps = int(scalar(text, "NTIMESTEPS"))
    d_s, d_l = (row[2] for row in tuple_values(text, "DIFFUSIVITY"))
    tempgrad = tuple_values(text, "Tempgrady")[0]
    base_t, delta_t, distance, _, velocity = tempgrad
    gradient = delta_t / distance

    assert scalar(text, "DIMENSION") == 2
    assert scalar(text, "NUMPHASES") == 2
    assert scalar(text, "NUMCOMPONENTS") == 2
    assert scalar(text, "Function_F") == 4
    assert scalar(text, "ISOTHERMAL") == 0
    assert int(scalar(text, "Shift")) == int(moving)
    assert math.isclose(epsilon / dx, 6.0, rel_tol=1e-12)
    assert d_s < d_l
    assert dt <= dx * dx / (4.0 * d_l), "explicit diffusion stability bound exceeded"
    assert gradient == 1.0e4
    assert velocity == 9.0e-5
    assert base_t >= 850.1497982085
    assert base_t + gradient * ny * dy <= 930.0

    wt_si = 0.055
    x_si = (wt_si / MW_SI) / (wt_si / MW_SI + (1.0 - wt_si) / MW_AL)
    rows = read_compositions()
    liquidus = interp(
        sorted(rows, key=lambda r: r["X_SI_LIQUID"]),
        "X_SI_LIQUID",
        "T",
        x_si,
    )
    ceq = tuple_values(text, "ceq")
    assert math.isclose(ceq[1][2], x_si, rel_tol=2e-10)
    assert abs(scalar(text, "Equilibrium_temperature") - liquidus) < 0.01

    print(f"{path.name}: OK")
    print(f"  grid={nx}x{ny}, size={nx*dx*1e6:.1f}x{ny*dy*1e6:.1f} um")
    print(f"  interface={epsilon/dx:.1f} cells, dt={dt:.2e} s")
    print(f"  physical time={nsteps*dt:.3f} s, isotherm travel={nsteps*dt*velocity*1e6:.1f} um")
    print(f"  xSi={x_si:.12f}, COST507 liquidus={liquidus:.6f} K")
    print(f"  G={gradient:.0f} K/m, V={velocity:.2e} m/s, G*V={gradient*velocity:.3f} K/s")


if __name__ == "__main__":
    validate(ROOT / "Input_smoke.in", moving=False)
    validate(ROOT / "Input_production.in", moving=True)
