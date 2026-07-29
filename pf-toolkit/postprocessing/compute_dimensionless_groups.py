#!/usr/bin/env python3
"""Compute reproducible scale and stability groups from a MicroSim input."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def scalar(text: str, key: str) -> float:
    match = re.search(rf"(?m)^{re.escape(key)}\s*=\s*({FLOAT})\s*;", text)
    if match is None:
        raise ValueError(f"missing scalar {key}")
    return float(match.group(1))


def integer(text: str, key: str) -> int:
    return int(scalar(text, key))


def braced(text: str, key: str) -> list[float]:
    match = re.search(
        rf"(?m)^{re.escape(key)}\s*=\s*\{{([^}}]+)\}}\s*;", text
    )
    if match is None:
        raise ValueError(f"missing braced value {key}")
    return [float(value) for value in re.findall(FLOAT, match.group(1))]


def diffusivities(text: str) -> tuple[float, float]:
    rows = re.findall(r"(?m)^DIFFUSIVITY\s*=\s*\{([^}]+)\}\s*;", text)
    values: dict[int, float] = {}
    for row in rows:
        parsed = [float(value) for value in re.findall(FLOAT, row)]
        if len(parsed) != 3:
            raise ValueError(f"unexpected DIFFUSIVITY row: {row}")
        values[int(parsed[1])] = parsed[2]
    if 0 not in values or 1 not in values:
        raise ValueError("expected solid phase 0 and liquid phase 1 diffusivities")
    return values[0], values[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument(
        "--active-tau-estimate",
        type=float,
        default=None,
        help="optional diagnostic F4 TAU estimate in s",
    )
    args = parser.parse_args()

    text = args.input.read_text()
    dx = scalar(text, "DELTA_X")
    dy = scalar(text, "DELTA_Y")
    dt = scalar(text, "DELTA_t")
    width = scalar(text, "epsilon")
    nx = integer(text, "MESH_X")
    ny = integer(text, "MESH_Y")
    d_s, d_l = diffusivities(text)
    temp = braced(text, "Tempgrady")
    if len(temp) != 5:
        raise ValueError("Tempgrady must contain five values")
    gradient = temp[1] / temp[2]
    velocity = temp[4]
    diffusion_length = d_l / velocity

    groups: list[dict[str, object]] = [
        {
            "quantity": "interface_cells",
            "value": width / dx,
            "unit": "1",
            "meaning": "diffuse-interface sampling",
        },
        {
            "quantity": "liquid_diffusion_Fourier",
            "value": dt * d_l / dx**2,
            "unit": "1",
            "meaning": "explicit liquid diffusion number dt*D_L/dx^2",
        },
        {
            "quantity": "solid_diffusion_Fourier",
            "value": dt * d_s / dx**2,
            "unit": "1",
            "meaning": "explicit solid diffusion number dt*D_S/dx^2",
        },
        {
            "quantity": "liquid_diffusion_margin_vs_0p2",
            "value": 0.2 / (dt * d_l / dx**2),
            "unit": "1",
            "meaning": "margin relative to a published Al-Si 0.2 bound",
        },
        {
            "quantity": "cell_advective_CFL",
            "value": velocity * dt / dy,
            "unit": "1",
            "meaning": "thermal-field travel per step in cells",
        },
        {
            "quantity": "interface_Peclet",
            "value": velocity * width / d_l,
            "unit": "1",
            "meaning": "V*W/D_L",
        },
        {
            "quantity": "cell_Peclet",
            "value": velocity * dx / d_l,
            "unit": "1",
            "meaning": "V*dx/D_L",
        },
        {
            "quantity": "diffusion_length",
            "value": diffusion_length,
            "unit": "m",
            "meaning": "D_L/V",
        },
        {
            "quantity": "gradient_across_interface",
            "value": gradient * width,
            "unit": "K",
            "meaning": "G*W",
        },
        {
            "quantity": "cooling_per_step",
            "value": gradient * velocity * dt,
            "unit": "K",
            "meaning": "G*V*dt",
        },
        {
            "quantity": "domain_width_over_diffusion_length",
            "value": nx * dx / diffusion_length,
            "unit": "1",
            "meaning": "Lx/(D_L/V)",
        },
        {
            "quantity": "domain_height_over_diffusion_length",
            "value": ny * dy / diffusion_length,
            "unit": "1",
            "meaning": "Ly/(D_L/V)",
        },
        {
            "quantity": "steps_per_isotherm_cell",
            "value": dy / (velocity * dt),
            "unit": "steps",
            "meaning": "steps for imposed field to move one y cell",
        },
    ]
    if args.active_tau_estimate is not None:
        groups.append(
            {
                "quantity": "F4_phase_Euler_prefactor_estimate",
                "value": 2.0 * dt / args.active_tau_estimate,
                "unit": "1",
                "meaning": "2*dt/TAU; excludes dimensional driving-force scale",
            }
        )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("quantity", "value", "unit", "meaning"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(groups)
    args.json.write_text(
        json.dumps(
            {
                "input": str(args.input),
                "parameters": {
                    "dx_m": dx,
                    "dy_m": dy,
                    "dt_s": dt,
                    "interface_width_m": width,
                    "D_s_m2_per_s": d_s,
                    "D_l_m2_per_s": d_l,
                    "G_K_per_m": gradient,
                    "V_m_per_s": velocity,
                    "mesh": [nx, ny],
                },
                "groups": groups,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {args.csv} and {args.json}")


if __name__ == "__main__":
    main()
