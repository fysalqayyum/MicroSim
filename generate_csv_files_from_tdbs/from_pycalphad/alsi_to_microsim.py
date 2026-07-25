#!/usr/bin/env python3
"""Extract binary Al-Si equilibrium thermodynamics for MicroSim.

The source TDB is never copied to the output directory. MatCalc-specific TDBs are
projected in memory onto LIQUID, FCC_A1, and SI_DIAMOND_A4 before pycalphad reads
them. MicroSim CSVs use mole fraction Si and J/mol curvature units.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pycalphad
import scipy
from pycalphad import Database, Model, variables as v
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, least_squares
from symengine import Symbol, lambdify


M_AL = 26.98154
M_SI = 28.0855
SOLID = "FCC_A1"
LIQUID = "LIQUID"
SILICON = "SI_DIAMOND_A4"
TARGET_PHASES = (SOLID, LIQUID, SILICON)
ENERGY_SCALE = 1.0e5


@dataclass(frozen=True)
class PhaseFunctions:
    gm: callable
    chemical_potential_difference: callable
    hessian: callable


def wt_percent_to_mole_fraction_si(wt_percent_si: float) -> float:
    weight_si = wt_percent_si / 100.0
    return (weight_si / M_SI) / ((weight_si / M_SI) + ((1.0 - weight_si) / M_AL))


def mole_fraction_to_wt_percent_si(x_si: float) -> float:
    return 100.0 * x_si * M_SI / ((1.0 - x_si) * M_AL + x_si * M_SI)


def _decode_tdb(path: Path) -> str:
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "windows-1252", "iso-8859-1"):
        try:
            return raw.decode(encoding).replace("\r\n", "\n").replace("\r", "\n")
        except UnicodeDecodeError:
            pass
    raise UnicodeError(f"Could not decode {path}")


def _tdb_commands(text: str):
    """Yield TDB commands without comments or the free-form reference appendix."""
    command_lines: list[str] = []
    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("$"):
            continue
        if not command_lines and stripped.upper().startswith("LIST_OF_REFERENCES"):
            break
        command_lines.append(stripped)
        if "!" in stripped:
            command = " ".join(command_lines)
            yield command[: command.index("!") + 1]
            command_lines = []
    if command_lines:
        raise ValueError("Unterminated TDB command at end of thermodynamic section")


def project_tdb_for_alsi(text: str) -> str:
    """Keep commands required to model the three Al-Si phases.

    This removes MatCalc-only commands such as REFERENCE_ELEMENT,
    ATTACH_CONTRIBUTION, and vacancy-enthalpy HMVA parameters. The projection is
    held in memory because source-database redistribution terms may be restrictive.
    """
    selected: list[str] = []
    for command in _tdb_commands(text):
        upper = command.upper().lstrip()
        if upper.startswith(("ELEMENT ", "FUNCTION ")):
            selected.append(command)
        elif upper.startswith("TYPE_DEFINITION "):
            selected.append(command)
        elif upper.startswith("PHASE "):
            if any(re.match(rf"PHASE\s+{re.escape(phase)}\b", upper) for phase in TARGET_PHASES):
                selected.append(command)
        elif upper.startswith("CONSTITUENT "):
            if any(re.match(rf"CONSTITUENT\s+{re.escape(phase)}\b", upper) for phase in TARGET_PHASES):
                selected.append(command)
        elif upper.startswith("PARAMETER "):
            if upper.startswith("PARAMETER HMVA"):
                continue
            if any(re.match(rf"PARAMETER\s+[A-Z]+\({re.escape(phase)}(?:,|;)", upper) for phase in TARGET_PHASES):
                selected.append(command)
    projected = "\n".join(selected) + "\n"
    missing = [phase for phase in TARGET_PHASES if f"PHASE {phase} " not in projected.upper()]
    if missing:
        raise ValueError(f"Target phases missing from projected TDB: {missing}")
    return projected


def load_projected_database(tdb_path: Path) -> Database:
    projected = project_tdb_for_alsi(_decode_tdb(tdb_path))
    return Database.from_string(projected, fmt="tdb")


def _scalar_compiled(expression, x_symbol):
    compiled = lambdify([v.T, x_symbol], [expression], backend="llvm")
    return lambda temperature, composition: float(compiled(temperature, composition))


def build_phase_functions(db: Database, phase: str) -> PhaseFunctions:
    model = Model(db, ["AL", "SI", "VA"], phase)
    x_si = Symbol("X_SI")
    substitutions = {}
    for site_fraction in model.site_fractions:
        label = str(site_fraction)
        if label.endswith(",AL)"):
            substitutions[site_fraction] = 1.0 - x_si
        elif label.endswith(",SI)"):
            substitutions[site_fraction] = x_si
        elif label.endswith(",VA)"):
            substitutions[site_fraction] = 1.0
        else:
            raise ValueError(f"Unexpected active site fraction {site_fraction} in {phase}")
    gm = model.GM.subs(substitutions)
    chemical_potential_difference = gm.diff(x_si)
    hessian = chemical_potential_difference.diff(x_si)
    return PhaseFunctions(
        _scalar_compiled(gm, x_si),
        _scalar_compiled(chemical_potential_difference, x_si),
        _scalar_compiled(hessian, x_si),
    )


def _common_tangent_residual(
    phase_a: PhaseFunctions,
    phase_b: PhaseFunctions,
    temperature: float,
    x_a: float,
    x_b: float,
) -> np.ndarray:
    g_a = phase_a.gm(temperature, x_a)
    g_b = phase_b.gm(temperature, x_b)
    mu_a = phase_a.chemical_potential_difference(temperature, x_a)
    mu_b = phase_b.chemical_potential_difference(temperature, x_b)
    return np.asarray(
        [
            (mu_a - mu_b) / ENERGY_SCALE,
            ((g_a - x_a * mu_a) - (g_b - x_b * mu_b)) / ENERGY_SCALE,
        ]
    )


def find_eutectic(functions: dict[str, PhaseFunctions]):
    def residual(z):
        temperature, x_solid, x_liquid, x_silicon = z
        return np.concatenate(
            [
                _common_tangent_residual(
                    functions[SOLID], functions[LIQUID], temperature, x_solid, x_liquid
                ),
                _common_tangent_residual(
                    functions[LIQUID], functions[SILICON], temperature, x_liquid, x_silicon
                ),
            ]
        )

    result = least_squares(
        residual,
        x0=[850.0, 0.016, 0.122, 0.999],
        bounds=([800.0, 1e-10, 1e-10, 0.5], [900.0, 0.20, 0.40, 1.0 - 1e-10]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=5000,
    )
    max_residual_j_per_mol = float(np.max(np.abs(residual(result.x))) * ENERGY_SCALE)
    if not result.success or max_residual_j_per_mol > 1e-3:
        raise RuntimeError(
            f"Eutectic common-tangent solve failed: {result.message}; "
            f"max residual={max_residual_j_per_mol:.6g} J/mol"
        )
    return result.x, max_residual_j_per_mol


def solve_solid_liquid_tie_line(
    functions: dict[str, PhaseFunctions],
    temperature: float,
    initial_guess: tuple[float, float],
):
    def residual(z):
        return _common_tangent_residual(
            functions[SOLID], functions[LIQUID], temperature, z[0], z[1]
        )

    result = least_squares(
        residual,
        x0=initial_guess,
        bounds=([1e-10, 1e-10], [0.20, 0.40]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=2000,
    )
    max_residual_j_per_mol = float(np.max(np.abs(residual(result.x))) * ENERGY_SCALE)
    x_solid, x_liquid = map(float, result.x)
    if (
        not result.success
        or max_residual_j_per_mol > 1e-3
        or not 0.0 < x_solid < x_liquid < 1.0
    ):
        raise RuntimeError(
            f"FCC_A1/LIQUID tie-line solve failed at {temperature:g} K: "
            f"{result.message}; x={result.x}; residual={max_residual_j_per_mol:.6g} J/mol"
        )
    return x_solid, x_liquid, max_residual_j_per_mol


def temperature_grid(t_min: float, t_max: float, step: float, eutectic_temperature: float):
    if not t_min < t_max or step <= 0:
        raise ValueError("Require t_min < t_max and step > 0")
    stable_minimum = max(t_min, eutectic_temperature)
    regular_start = math.ceil(stable_minimum / step) * step
    regular = np.arange(regular_start, t_max + 0.5 * step, step)
    refined = np.asarray(
        [
            eutectic_temperature,
            eutectic_temperature + 0.1,
            eutectic_temperature + 0.2,
            eutectic_temperature + 0.5,
        ]
    )
    values = np.concatenate([regular, refined])
    values = values[(values >= stable_minimum - 1e-9) & (values <= t_max + 1e-9)]
    return np.unique(np.round(values, 10))


def validate_reference_eutectic(eutectic, expected_temperature, expected_liquid_wt, expected_solid_wt):
    temperature, x_solid, x_liquid, _ = eutectic
    errors = {
        "temperature_K": abs(float(temperature) - expected_temperature),
        "liquid_wt_percent_Si": abs(mole_fraction_to_wt_percent_si(float(x_liquid)) - expected_liquid_wt),
        "solid_wt_percent_Si": abs(mole_fraction_to_wt_percent_si(float(x_solid)) - expected_solid_wt),
    }
    tolerances = {
        "temperature_K": 2.0,
        "liquid_wt_percent_Si": 0.2,
        "solid_wt_percent_Si": 0.2,
    }
    failed = {key: value for key, value in errors.items() if value > tolerances[key]}
    if failed:
        raise RuntimeError(f"Eutectic validation failed: errors={failed}, tolerances={tolerances}")
    return errors, tolerances


def _write_csvs(output_dir: Path, diagnostics: pd.DataFrame):
    microsim_dir = output_dir / "tdbs_encrypted"
    microsim_dir.mkdir(parents=True, exist_ok=True)
    diagnostics[["T_K", "x_Si_FCC_A1", "x_Si_LIQUID"]].rename(
        columns={"T_K": "T", "x_Si_FCC_A1": "X_SI_FCC_A1", "x_Si_LIQUID": "X_SI_LIQUID"}
    ).to_csv(microsim_dir / "Composition_FCC_A1.csv", index=False)
    diagnostics[["T_K", "H_FCC_A1_J_per_mol"]].rename(
        columns={"T_K": "T", "H_FCC_A1_J_per_mol": "HSN(SI,SI)@FCC_A1"}
    ).to_csv(microsim_dir / "HSN_FCC_A1.csv", index=False)
    diagnostics[["T_K", "H_LIQUID_J_per_mol"]].rename(
        columns={"T_K": "T", "H_LIQUID_J_per_mol": "HSN(SI,SI)@LIQUID"}
    ).to_csv(microsim_dir / "HSN_LIQUID.csv", index=False)
    diagnostics.to_csv(output_dir / "AlSi_thermodynamic_diagnostics.csv", index=False)
    return microsim_dir


def validate_natural_cubic_interpolation(
    diagnostics: pd.DataFrame,
    functions: dict[str, PhaseFunctions],
    eutectic_guess: tuple[float, float],
):
    """Compare MicroSim-equivalent natural cubic splines with exact midpoint values."""
    first_regular_midpoint = math.ceil(float(diagnostics.T_K.min())) + 0.5
    midpoints = np.arange(first_regular_midpoint, float(diagnostics.T_K.max()), 1.0)
    columns = {
        "x_Si_FCC_A1": [],
        "x_Si_LIQUID": [],
        "H_FCC_A1_J_per_mol": [],
        "H_LIQUID_J_per_mol": [],
    }
    guess = eutectic_guess
    for temperature in midpoints:
        x_solid, x_liquid, _ = solve_solid_liquid_tie_line(functions, temperature, guess)
        guess = (x_solid, x_liquid)
        columns["x_Si_FCC_A1"].append(x_solid)
        columns["x_Si_LIQUID"].append(x_liquid)
        columns["H_FCC_A1_J_per_mol"].append(functions[SOLID].hessian(temperature, x_solid))
        columns["H_LIQUID_J_per_mol"].append(functions[LIQUID].hessian(temperature, x_liquid))
    errors = {}
    for name, exact in columns.items():
        spline = CubicSpline(diagnostics.T_K, diagnostics[name], bc_type="natural")
        exact_array = np.asarray(exact)
        interpolated = spline(midpoints)
        errors[name] = float(np.max(np.abs((interpolated - exact_array) / exact_array)))
    if max(errors.values()) > 0.01:
        raise RuntimeError(f"Natural-cubic interpolation validation exceeded 1%: {errors}")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdb", type=Path, required=True, help="User-supplied thermodynamic TDB")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--t-min", type=float, default=850.0)
    parser.add_argument("--t-max", type=float, default=930.0)
    parser.add_argument("--step", type=float, default=1.0)
    parser.add_argument("--kks-temperature", type=float, default=860.0)
    parser.add_argument(
        "--target-wt-percent-si",
        type=float,
        nargs="+",
        default=[3.5, 4.5, 5.5, 6.5, 7.5],
    )
    args = parser.parse_args()

    if not args.tdb.is_file():
        parser.error(f"TDB does not exist: {args.tdb}")
    args.output.mkdir(parents=True, exist_ok=True)
    db = load_projected_database(args.tdb)
    if not set(TARGET_PHASES).issubset(db.phases):
        raise RuntimeError(f"Parsed phases are {sorted(db.phases)}; expected {TARGET_PHASES}")
    functions = {phase: build_phase_functions(db, phase) for phase in TARGET_PHASES}

    eutectic, eutectic_residual = find_eutectic(functions)
    validation_errors, validation_tolerances = validate_reference_eutectic(
        eutectic, expected_temperature=850.0, expected_liquid_wt=12.6, expected_solid_wt=1.6
    )
    eutectic_temperature, eutectic_solid, eutectic_liquid, eutectic_silicon = map(float, eutectic)

    rows = []
    guess = (eutectic_solid, eutectic_liquid)
    for temperature in temperature_grid(args.t_min, args.t_max, args.step, eutectic_temperature):
        x_solid, x_liquid, tangent_residual = solve_solid_liquid_tie_line(
            functions, float(temperature), guess
        )
        guess = (x_solid, x_liquid)
        h_solid = functions[SOLID].hessian(float(temperature), x_solid)
        h_liquid = functions[LIQUID].hessian(float(temperature), x_liquid)
        if h_solid <= 0.0 or h_liquid <= 0.0:
            raise RuntimeError(f"Non-positive phase Hessian at {temperature:g} K")
        rows.append(
            {
                "T_K": float(temperature),
                "x_Si_FCC_A1": x_solid,
                "x_Si_LIQUID": x_liquid,
                "wt_percent_Si_FCC_A1": mole_fraction_to_wt_percent_si(x_solid),
                "wt_percent_Si_LIQUID": mole_fraction_to_wt_percent_si(x_liquid),
                "H_FCC_A1_J_per_mol": h_solid,
                "H_LIQUID_J_per_mol": h_liquid,
                "KKS_f0_FCC_A1_J_per_mol": 0.5 * h_solid,
                "KKS_f0_LIQUID_J_per_mol": 0.5 * h_liquid,
                "common_tangent_max_residual_J_per_mol": tangent_residual,
            }
        )
    diagnostics = pd.DataFrame(rows)
    interpolation_errors = validate_natural_cubic_interpolation(
        diagnostics, functions, (eutectic_solid, eutectic_liquid)
    )
    microsim_dir = _write_csvs(args.output, diagnostics)

    def exact_tie_line(temperature):
        x_solid, x_liquid, _ = solve_solid_liquid_tie_line(
            functions, temperature, (eutectic_solid, eutectic_liquid)
        )
        return x_solid, x_liquid

    target_rows = []
    upper_liquidus_limit = min(933.45, max(args.t_max, args.kks_temperature))
    for wt_percent_si in args.target_wt_percent_si:
        x_nominal = wt_percent_to_mole_fraction_si(wt_percent_si)

        def liquidus_objective(temperature):
            return exact_tie_line(temperature)[1] - x_nominal

        liquidus_temperature = brentq(
            liquidus_objective, eutectic_temperature + 1e-5, upper_liquidus_limit
        )
        x_solid_at_liquidus, _ = exact_tie_line(liquidus_temperature)
        target_rows.append(
            {
                "nominal_wt_percent_Si": wt_percent_si,
                "nominal_x_Si": x_nominal,
                "liquidus_T_K": liquidus_temperature,
                "x_Si_FCC_A1_at_liquidus": x_solid_at_liquidus,
                "wt_percent_Si_FCC_A1_at_liquidus": mole_fraction_to_wt_percent_si(
                    x_solid_at_liquidus
                ),
                "partition_coefficient_x_basis": x_solid_at_liquidus / x_nominal,
                "equilibrium_eutectic_fraction_mole_lever_rule": (
                    (x_nominal - eutectic_solid) / (eutectic_liquid - eutectic_solid)
                ),
                "equilibrium_eutectic_fraction_weight_lever_rule": (
                    (wt_percent_si - mole_fraction_to_wt_percent_si(eutectic_solid))
                    / (
                        mole_fraction_to_wt_percent_si(eutectic_liquid)
                        - mole_fraction_to_wt_percent_si(eutectic_solid)
                    )
                ),
            }
        )
    targets = pd.DataFrame(target_rows)
    targets.to_csv(args.output / "AlSi_target_alloys.csv", index=False)

    kks_x_solid, kks_x_liquid = exact_tie_line(args.kks_temperature)
    kks_h_solid = functions[SOLID].hessian(args.kks_temperature, kks_x_solid)
    kks_h_liquid = functions[LIQUID].hessian(args.kks_temperature, kks_x_liquid)
    kks_fragment = (
        "# Fixed-temperature KKS_CuFFT thermodynamic fragment.\n"
        "# Phase 0 = FCC_A1; phase 1 = LIQUID; independent component = Si.\n"
        "# f0 = H/2 in J/mol. KKS divides f0 by Vm; supply a validated Vm in m^3/mol.\n"
        f"# T = {args.kks_temperature:.10g} K. This is not a temperature-dependent reader.\n"
        f"f0 = {{0,{0.5 * kks_h_solid:.16g}}};\n"
        f"f0 = {{1,{0.5 * kks_h_liquid:.16g}}};\n"
        f"ceq = {{0,0,{kks_x_solid:.16g}}};\n"
        f"ceq = {{0,1,{kks_x_liquid:.16g}}};\n"
        f"ceq = {{1,0,{kks_x_solid:.16g}}};\n"
        f"ceq = {{1,1,{kks_x_liquid:.16g}}};\n"
    )
    (args.output / f"KKS_CuFFT_{args.kks_temperature:g}K_parameters.txt").write_text(
        kks_fragment, encoding="utf-8"
    )

    metadata = {
        "source": {
            "filename": args.tdb.name,
            "sha256": hashlib.sha256(args.tdb.read_bytes()).hexdigest(),
            "redistributed": False,
            "projection": "In-memory LIQUID/FCC_A1/SI_DIAMOND_A4 projection",
        },
        "software": {
            "python": platform.python_version(),
            "pycalphad": pycalphad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "basis_and_units": {
            "composition": "mole fraction Si; Al is the dependent solvent",
            "hessian": "d2(G_m)/d(x_Si)^2 = d(mu_Si-mu_Al)/d(x_Si), J/mol",
            "kks_f0": "H/2, J/mol; KKS_CuFFT divides by Vm",
            "temperature": "K",
            "molar_masses_g_per_mol": {"Al": M_AL, "Si": M_SI},
        },
        "requested_temperature_range_K": [args.t_min, args.t_max],
        "stable_two_phase_temperature_range_written_K": [
            float(diagnostics.T_K.min()),
            float(diagnostics.T_K.max()),
        ],
        "eutectic": {
            "T_K": eutectic_temperature,
            "x_Si_FCC_A1": eutectic_solid,
            "x_Si_LIQUID": eutectic_liquid,
            "x_Si_SI_DIAMOND_A4": eutectic_silicon,
            "wt_percent_Si_FCC_A1": mole_fraction_to_wt_percent_si(eutectic_solid),
            "wt_percent_Si_LIQUID": mole_fraction_to_wt_percent_si(eutectic_liquid),
            "wt_percent_Si_SI_DIAMOND_A4": mole_fraction_to_wt_percent_si(eutectic_silicon),
            "common_tangent_max_residual_J_per_mol": eutectic_residual,
        },
        "validation": {
            "reference_targets": {
                "T_K": 850.0,
                "wt_percent_Si_LIQUID": 12.6,
                "wt_percent_Si_FCC_A1": 1.6,
            },
            "absolute_errors": validation_errors,
            "tolerances": validation_tolerances,
            "natural_cubic_midpoint_max_relative_errors": interpolation_errors,
            "natural_cubic_midpoint_relative_tolerance": 0.01,
            "passed": True,
        },
        "microsim": {
            "grand_potential_mpi": {
                "Function_F": 4,
                "NUM_THERMO_PHASES": 2,
                "tdb_phases": [SOLID, LIQUID],
                "phase_order_requirement": "LIQUID must be the final simulation phase",
                "csv_directory": microsim_dir.name,
            },
            "kks_cufft": {
                "mapping": "isothermal f0=H/2 and diagonal ceq",
                "temperature_K": args.kks_temperature,
                "limitation": (
                    "Current KKS_CuFFT has no TDB/CSV thermodynamic reader and no "
                    "temperature-dependent directional-solidification free energy."
                ),
            },
        },
        "mobility": {
            "extracted": False,
            "reason": (
                "A thermodynamic TDB supplies Gibbs models, not atomic mobilities or "
                "phase diffusivities; a separately licensed mobility database or "
                "citable experimental diffusion model is required."
            ),
        },
    }
    (args.output / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    print(f"Wrote {len(diagnostics)} thermodynamic rows to {args.output}")
    print(
        "Eutectic: "
        f"{eutectic_temperature:.4f} K, "
        f"FCC_A1={mole_fraction_to_wt_percent_si(eutectic_solid):.4f} wt.% Si, "
        f"LIQUID={mole_fraction_to_wt_percent_si(eutectic_liquid):.4f} wt.% Si"
    )


if __name__ == "__main__":
    main()
