#!/usr/bin/env python3
"""Quantify and plot morphology in MicroSim rank-local binary fragments.

The MicroSim files analysed here are not standalone VTK files.  They begin with
three ASCII dimension lines followed by five big-endian float64 arrays.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


DT_S = 4.0e-7
DX_M = 7.5e-7
PHASE_THRESHOLD = 0.5
SUBSTRATE_TOP_CELL = 23
FIELDS = ("FCC_A1", "LIQUID", "Mu_Si", "Composition_Si", "Temperature")
OKABE_ITO = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
MM_TO_IN = 1.0 / 25.4


def timestep(path: Path) -> int:
    match = re.search(r"_(\d+)\.vtk$", path.name)
    if match is None:
        raise ValueError(f"cannot infer timestep from {path.name}")
    return int(match.group(1))


def read_fragment(path: Path) -> tuple[tuple[int, int, int], dict[str, np.ndarray]]:
    with path.open("rb") as handle:
        dimensions = tuple(int(handle.readline()) for _ in range(3))
        count = int(np.prod(dimensions))
        arrays: dict[str, np.ndarray] = {}
        for index, name in enumerate(FIELDS):
            values = np.fromfile(handle, dtype=">f8", count=count).astype(np.float64)
            if values.size != count:
                raise ValueError(f"{path}: incomplete {name} array")
            if not np.isfinite(values).all():
                raise ValueError(
                    f"{path}: {np.count_nonzero(~np.isfinite(values))} "
                    f"non-finite values in {name}"
                )
            # Writer ordering is x-major/y-minor, so y is the fastest index.
            arrays[name] = values.reshape(dimensions)
            if index < len(FIELDS) - 1 and handle.read(1) != b"\n":
                raise ValueError(f"{path}: missing separator after {name}")
    return dimensions, arrays


def interface_height(solid: np.ndarray) -> np.ndarray:
    """Highest solid cell at every x; -1 marks an empty column."""
    nx, ny = solid.shape
    result = np.full(nx, -1, dtype=int)
    for x in range(nx):
        occupied = np.flatnonzero(solid[x])
        if occupied.size:
            result[x] = int(occupied[-1])
    return result


def subcell_interface_height(phase: np.ndarray) -> np.ndarray:
    """Linearly interpolate the uppermost phi=0.5 crossing in every x column."""
    nx, ny = phase.shape
    result = np.full(nx, np.nan, dtype=float)
    for x in range(nx):
        crossings = np.flatnonzero(
            (phase[x, :-1] >= PHASE_THRESHOLD)
            & (phase[x, 1:] < PHASE_THRESHOLD)
        )
        if crossings.size:
            lower = int(crossings[-1])
            denominator = phase[x, lower] - phase[x, lower + 1]
            fraction = (
                (phase[x, lower] - PHASE_THRESHOLD) / denominator
                if denominator != 0
                else 0.0
            )
            result[x] = lower + fraction
    return result


def morphology_row(path: Path) -> tuple[dict[str, float | int | str], dict[str, np.ndarray]]:
    dimensions, arrays = read_fragment(path)
    if dimensions[2] != 1:
        raise ValueError("this analysis is defined for the 2-D case only")
    phase = arrays["FCC_A1"][:, :, 0]
    composition = arrays["Composition_Si"][:, :, 0]
    temperature = arrays["Temperature"][:, :, 0]
    solid = phase >= PHASE_THRESHOLD
    height = interface_height(solid)
    subcell_height = subcell_interface_height(phase)
    above_substrate = solid.copy()
    above_substrate[:, : SUBSTRATE_TOP_CELL + 1] = False

    labels, component_count = ndimage.label(above_substrate)
    component_sizes = np.bincount(labels.ravel())[1:]
    retained_components = int(np.count_nonzero(component_sizes >= 4))
    eroded = ndimage.binary_erosion(above_substrate)
    perimeter_cells = int(np.count_nonzero(above_substrate & ~eroded))
    active = height > SUBSTRATE_TOP_CELL
    active_height = height[active]
    distance_from_solid = ndimage.distance_transform_edt(~solid)
    liquid_core = phase < 0.1
    solid_core = phase > 0.9
    near_interface_liquid = (
        liquid_core & (distance_from_solid >= 1.0) & (distance_from_solid <= 12.0)
    )
    far_liquid = liquid_core & (distance_from_solid >= 40.0)

    def selected_mean(mask: np.ndarray) -> float:
        selected = composition[mask]
        return float(np.mean(selected)) if selected.size else float("nan")

    step = timestep(path)

    row: dict[str, float | int | str] = {
        "file": path.name,
        "step": step,
        "time_s": step * DT_S,
        "solid_area_fraction": float(np.mean(solid)),
        "solid_area_above_substrate_cells": int(np.count_nonzero(above_substrate)),
        "tip_cell": int(height.max()),
        "tip_um": float(height.max() * DX_M * 1e6),
        "tip_subcell": float(np.nanmax(subcell_height)),
        "tip_subcell_um": float(np.nanmax(subcell_height) * DX_M * 1e6),
        "phase_integral_cells": float(phase.sum()),
        "active_lateral_width_cells": int(np.count_nonzero(active)),
        "interface_mean_cell_active": float(np.mean(active_height)) if active_height.size else np.nan,
        "interface_std_cell_active": float(np.std(active_height)) if active_height.size else np.nan,
        "interface_total_variation_cells": int(np.abs(np.diff(height)).sum()),
        "retained_components_above_substrate": retained_components,
        "perimeter_cells_above_substrate": perimeter_cells,
        "phase_min": float(phase.min()),
        "phase_max": float(phase.max()),
        "xSi_min": float(composition.min()),
        "xSi_mean": float(composition.mean()),
        "xSi_max": float(composition.max()),
        "xSi_solid_core_mean": selected_mean(solid_core),
        "xSi_near_interface_liquid_mean": selected_mean(near_interface_liquid),
        "xSi_far_liquid_mean": selected_mean(far_liquid),
        "solute_rejection_delta_xSi": (
            selected_mean(near_interface_liquid) - selected_mean(far_liquid)
        ),
        "temperature_min_K": float(temperature.min()),
        "temperature_max_K": float(temperature.max()),
    }
    return row, arrays


def apply_style() -> None:
    try:
        import scienceplots  # noqa: F401

        plt.style.use(["science", "no-latex"])
    except ImportError:
        plt.style.use("default")
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def plot_phase_sequence(
    records: list[tuple[dict[str, float | int | str], dict[str, np.ndarray]]],
    output_dir: Path,
) -> None:
    fig, axes = plt.subplots(
        1,
        len(records),
        figsize=(180 * MM_TO_IN, 82 * MM_TO_IN),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    extent = [0, records[0][1]["FCC_A1"].shape[0] * DX_M * 1e6,
              0, records[0][1]["FCC_A1"].shape[1] * DX_M * 1e6]
    image = None
    for panel, (row, arrays) in enumerate(records):
        ax = np.atleast_1d(axes)[panel]
        phase = arrays["FCC_A1"][:, :, 0]
        image = ax.imshow(
            phase.T,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="viridis",
            vmin=0,
            vmax=1,
            interpolation="nearest",
            rasterized=True,
        )
        ax.contour(
            np.linspace(0, extent[1], phase.shape[0]),
            np.linspace(0, extent[3], phase.shape[1]),
            phase.T,
            levels=[PHASE_THRESHOLD],
            colors="white",
            linewidths=0.45,
        )
        ax.set_title(f"{float(row['time_s']):.3f} s")
        ax.set_xlabel(r"$x$ ($\mu$m)")
        ax.text(
            0.03,
            0.97,
            f"({chr(97 + panel)})",
            transform=ax.transAxes,
            va="top",
            color="white",
            fontweight="bold",
        )
    np.atleast_1d(axes)[0].set_ylabel(r"Growth direction $y$ ($\mu$m)")
    if image is not None:
        colorbar = fig.colorbar(image, ax=np.atleast_1d(axes), fraction=0.025, pad=0.02)
        colorbar.set_label(r"FCC-A1 phase field $\phi$")
    fig.savefig(output_dir / "pilot_phase_morphology.png", dpi=400)
    plt.close(fig)


def plot_final_fields(
    record: tuple[dict[str, float | int | str], dict[str, np.ndarray]],
    output_dir: Path,
) -> None:
    row, arrays = record
    fig, axes = plt.subplots(
        1, 2, figsize=(120 * MM_TO_IN, 98 * MM_TO_IN), constrained_layout=True
    )
    extent = [0, arrays["FCC_A1"].shape[0] * DX_M * 1e6,
              0, arrays["FCC_A1"].shape[1] * DX_M * 1e6]
    specifications = (
        ("FCC_A1", "viridis", 0, 1, r"FCC-A1 phase field $\phi$"),
        (
            "Composition_Si",
            "plasma",
            float(arrays["Composition_Si"].min()),
            float(arrays["Composition_Si"].max()),
            r"Si mole fraction $x_{\rm Si}$",
        ),
    )
    for panel, (ax, (name, cmap, vmin, vmax, label)) in enumerate(
        zip(axes, specifications)
    ):
        values = arrays[name][:, :, 0]
        image = ax.imshow(
            values.T,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            rasterized=True,
        )
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label(label)
        ax.set_xlabel(r"$x$ ($\mu$m)")
        ax.set_ylabel(r"Growth direction $y$ ($\mu$m)")
        ax.text(
            0.03,
            0.97,
            f"({chr(97 + panel)})",
            transform=ax.transAxes,
            va="top",
            color="white",
            fontweight="bold",
        )
    fig.suptitle(f"Production pilot at {float(row['time_s']):.3f} s")
    fig.savefig(output_dir / "pilot_final_phase_composition.png", dpi=400)
    plt.close(fig)


def plot_metrics(rows: list[dict[str, float | int | str]], output_dir: Path) -> None:
    time = np.asarray([float(row["time_s"]) for row in rows])
    tip = np.asarray([float(row["tip_um"]) for row in rows])
    area = np.asarray([float(row["solid_area_fraction"]) for row in rows])
    roughness = np.asarray([float(row["interface_std_cell_active"]) for row in rows])

    fig, axes = plt.subplots(
        1, 3, figsize=(180 * MM_TO_IN, 55 * MM_TO_IN), constrained_layout=True
    )
    series = (
        (tip, r"Tip position ($\mu$m)", "o-", OKABE_ITO[0]),
        (area, "Solid area fraction", "s--", OKABE_ITO[1]),
        (roughness, "Interface-height SD (cells)", "^-.", OKABE_ITO[2]),
    )
    for panel, (ax, (values, ylabel, style, color)) in enumerate(zip(axes, series)):
        ax.plot(time, values, style, color=color, markersize=4)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(ylabel)
        ax.minorticks_on()
        ax.text(
            -0.16,
            1.06,
            f"({chr(97 + panel)})",
            transform=ax.transAxes,
            va="top",
            fontweight="bold",
        )
    fig.savefig(output_dir / "pilot_morphology_metrics.pdf")
    fig.savefig(output_dir / "pilot_morphology_metrics.png", dpi=300)
    plt.close(fig)


def main() -> None:
    global DT_S, DX_M, SUBSTRATE_TOP_CELL

    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=DT_S, help="time step in s")
    parser.add_argument("--dx", type=float, default=DX_M, help="cell spacing in m")
    parser.add_argument(
        "--substrate-top-cell",
        type=int,
        default=SUBSTRATE_TOP_CELL,
        help="highest prescribed substrate row",
    )
    args = parser.parse_args()
    DT_S = args.dt
    DX_M = args.dx
    SUBSTRATE_TOP_CELL = args.substrate_top_cell

    paths = sorted(args.input_dir.glob("*.vtk"), key=timestep)
    if not paths:
        raise SystemExit(f"no rank fragments found in {args.input_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records = [morphology_row(path) for path in paths]
    rows = [record[0] for record in records]
    initial_phase = records[0][1]["FCC_A1"][:, :, 0]
    for row, (_, arrays) in zip(rows, records):
        difference = arrays["FCC_A1"][:, :, 0] - initial_phase
        row["phase_change_l1_cells"] = float(np.abs(difference).sum())
        row["phase_change_rms"] = float(np.sqrt(np.mean(difference**2)))
    with (args.output_dir / "pilot_morphology_metrics.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    apply_style()
    plot_phase_sequence(records, args.output_dir)
    plot_final_fields(records[-1], args.output_dir)
    plot_metrics(rows, args.output_dir)
    print(f"analysed {len(records)} frames in {args.output_dir}")


if __name__ == "__main__":
    main()
