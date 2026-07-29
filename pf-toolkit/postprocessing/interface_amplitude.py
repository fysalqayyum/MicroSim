#!/usr/bin/env python3
"""Is the noise-seeded flat front destabilising?

Reads MicroSim rank-local binary fragments (three ASCII dimension lines, then
five big-endian float64 arrays: FCC_A1, LIQUID, Mu_Si, Composition_Si,
Temperature) and tracks the amplitude of the solid-liquid interface h(x),
defined as the uppermost phi = 0.5 crossing in every x column, interpolated to
sub-cell precision.  Only the first field is read; the rest are skipped.

The x boundary is periodic (BOUNDARY = {phi,3,3,1,1,3,3}), so the discrete
Fourier transform of h(x) - <h> is the physically correct mode decomposition.

Reports, per frame: RMS amplitude, peak-to-valley, protrusion count, and the
dominant wavelength.  Fits sigma from ln(RMS) vs t for comparison with the
multimode-seed growth rate.

MOVING WINDOW.  h(x) is measured in window-local coordinates.  Once the window
starts shifting, lab-frame height = window-local height + cumulative shift, so
every absolute height (and the front velocity derived from it) needs shift.dat.
The amplitude quantities -- RMS, peak-to-valley, the spectrum -- are taken about
each frame's own mean and are shift-invariant, so they stay valid either way.
Pass --shift-file, or --no-shift to assert the window has not fired; asserting
that wrongly is a hard error rather than a silently wrong velocity.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
from scipy import signal

PHASE_THRESHOLD = 0.5
N_FIELDS = 5


def timestep(path: Path) -> int:
    match = re.search(r"_(\d+)\.vtk$", path.name)
    if match is None:
        raise ValueError(f"cannot infer timestep from {path.name}")
    return int(match.group(1))


def read_shift_map(path: Path) -> dict[int, int]:
    """step -> cumulative window shift in cells, as written by the solver.

    One line per SAVET step, "<step> <cumulative_cells>", starting at the first
    output step; step 0 is absent and means zero shift.
    """
    mapping: dict[int, int] = {}
    for line in path.read_text().split("\n"):
        parts = line.split()
        if len(parts) == 2:
            mapping[int(parts[0])] = int(parts[1])
    return mapping


def shift_at(shift_map: dict[int, int], step: int) -> int:
    """Cumulative shift at `step`: the last entry at or before it, else 0."""
    if step in shift_map:
        return shift_map[step]
    earlier = [s for s in shift_map if s <= step]
    return shift_map[max(earlier)] if earlier else 0


def read_phase(path: Path) -> np.ndarray:
    """Return FCC_A1 as (nx, ny); the remaining four fields are not read."""
    with path.open("rb") as handle:
        dimensions = tuple(int(handle.readline()) for _ in range(3))
        if dimensions[2] != 1:
            raise ValueError(f"{path}: this diagnostic is 2-D only")
        count = int(np.prod(dimensions))
        values = np.fromfile(handle, dtype=">f8", count=count)
    if values.size != count:
        raise ValueError(f"{path}: incomplete FCC_A1 array")
    nonfinite = int(np.count_nonzero(~np.isfinite(values)))
    if nonfinite:
        raise ValueError(f"{path}: {nonfinite} non-finite values in FCC_A1")
    return values.reshape(dimensions)[:, :, 0]


def interface_profile(phase: np.ndarray) -> np.ndarray:
    """Uppermost phi = 0.5 crossing per x column, linearly interpolated."""
    nx, _ = phase.shape
    result = np.full(nx, np.nan, dtype=float)
    below = phase[:, :-1] >= PHASE_THRESHOLD
    above = phase[:, 1:] < PHASE_THRESHOLD
    crossing = below & above
    for x in range(nx):
        indices = np.flatnonzero(crossing[x])
        if indices.size:
            lower = int(indices[-1])
            denominator = phase[x, lower] - phase[x, lower + 1]
            fraction = (
                (phase[x, lower] - PHASE_THRESHOLD) / denominator
                if denominator != 0.0
                else 0.0
            )
            result[x] = lower + fraction
    return result


def spectrum(deviation: np.ndarray, dx_m: float) -> tuple[np.ndarray, np.ndarray]:
    """One-sided power spectrum of a periodic profile; mode 0 dropped."""
    nx = deviation.size
    amplitude = np.abs(np.fft.rfft(deviation)) * 2.0 / nx
    modes = np.arange(amplitude.size)
    with np.errstate(divide="ignore"):
        wavelength_um = np.where(modes > 0, nx * dx_m * 1e6 / np.maximum(modes, 1), np.inf)
    return amplitude[1:], wavelength_um[1:]


def analyse(path: Path, dx_m: float, dt_s: float, substrate_top: int,
            shift_cells: int) -> dict:
    phase = read_phase(path)
    height = interface_profile(phase)
    valid = np.isfinite(height)
    if valid.sum() < phase.shape[0] // 2:
        raise ValueError(f"{path}: interface lost in {(~valid).sum()} of {valid.size} columns")

    # A column with no crossing cannot be left as NaN in a periodic FFT; the
    # front is continuous, so interpolate across the few gaps.
    if not valid.all():
        indices = np.arange(height.size)
        height = np.interp(indices, indices[valid], height[valid])

    deviation = height - height.mean()
    amplitude, wavelength_um = spectrum(deviation, dx_m)
    order = np.argsort(amplitude)[::-1]

    rms_cells = float(np.sqrt(np.mean(deviation**2)))
    prominence = max(0.5 * rms_cells, 0.25)
    peaks, _ = signal.find_peaks(deviation, prominence=prominence)

    # Below this the profile is flat to within float noise and every spectral
    # quantity is a code signature rather than a measurement.
    resolved = rms_cells > 1e-9

    def wavelength_of(rank: int) -> float:
        return float(wavelength_um[order[rank]]) if resolved else np.nan

    step = timestep(path)
    # Lab frame = window-local + cumulative shift. The substrate was seeded at
    # lab row `substrate_top` and never moves in the lab frame, so h - substrate
    # stays the total solidified distance once the window starts scrolling.
    h_mean_lab = float(height.mean()) + shift_cells
    h_max_lab = float(height.max()) + shift_cells
    h_min_lab = float(height.min()) + shift_cells
    return {
        "file": path.name,
        "step": step,
        "time_s": step * dt_s,
        "n_columns_no_crossing": int((~valid).sum()),
        "cumulative_shift_cells": int(shift_cells),
        "h_mean_cells": h_mean_lab,
        "h_mean_um": h_mean_lab * dx_m * 1e6,
        "h_tip_cells": h_max_lab,
        "h_tip_um": h_max_lab * dx_m * 1e6,
        "h_groove_cells": h_min_lab,
        "h_above_substrate_cells": h_mean_lab - substrate_top,
        "rms_amplitude_cells": rms_cells,
        "rms_amplitude_um": rms_cells * dx_m * 1e6,
        "peak_to_valley_cells": float(deviation.max() - deviation.min()),
        "peak_to_valley_um": float((deviation.max() - deviation.min()) * dx_m * 1e6),
        "protrusion_count": int(peaks.size),
        "mean_protrusion_spacing_um": (
            float(phase.shape[0] * dx_m * 1e6 / peaks.size) if peaks.size else np.nan
        ),
        "dominant_wavelength_um": wavelength_of(0),
        "dominant_mode_amplitude_cells": float(amplitude[order[0]]) if resolved else np.nan,
        "second_wavelength_um": wavelength_of(1),
        "third_wavelength_um": wavelength_of(2),
        "_deviation": deviation,
    }


def add_front_velocities(rows: list[dict], dx_m: float) -> None:
    """Backward-difference mean and tip front velocity, in mm/s, in place.

    The first frame has no predecessor, so its velocities are NaN rather than a
    zero that would read as a measurement.
    """
    previous = None
    for row in rows:
        if previous is None or row["time_s"] <= previous["time_s"]:
            row["v_mean_mm_s"] = np.nan
            row["v_tip_mm_s"] = np.nan
        else:
            span = row["time_s"] - previous["time_s"]
            scale = dx_m * 1e3 / span
            row["v_mean_mm_s"] = (row["h_mean_cells"] - previous["h_mean_cells"]) * scale
            row["v_tip_mm_s"] = (row["h_tip_cells"] - previous["h_tip_cells"]) * scale
        previous = row


def growth_rate(times: np.ndarray, rms: np.ndarray) -> tuple[float, float, int]:
    """Least-squares sigma from ln(RMS) = ln(A0) + sigma t over growing frames."""
    usable = np.isfinite(rms) & (rms > 0)
    if usable.sum() < 3:
        return np.nan, np.nan, int(usable.sum())
    t = times[usable]
    y = np.log(rms[usable])
    slope, intercept = np.polyfit(t, y, 1)
    residual = y - (slope * t + intercept)
    total = y - y.mean()
    r_squared = 1.0 - residual @ residual / (total @ total) if total @ total > 0 else np.nan
    return float(slope), float(r_squared), int(usable.sum())


# --- plotting -------------------------------------------------------------

MM_TO_IN = 1 / 25.4
OKABE_ITO = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")


def apply_style() -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    try:
        import scienceplots  # noqa: F401

        plt.style.use(["science", "no-latex"])
    except Exception:
        plt.style.use("default")
    mpl.rcParams.update({
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "legend.title_fontsize": 8,
        "lines.linewidth": 1.5, "lines.markersize": 5, "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "figure.dpi": 150, "text.usetex": False,
    })


def make_figure(rows: list[dict], dx_m: float, sigma: float, output_dir: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM_TO_IN, 70 * MM_TO_IN),
                             constrained_layout=True)

    time_ms = np.asarray([row["time_s"] for row in rows]) * 1e3
    rms_um = np.asarray([row["rms_amplitude_um"] for row in rows])
    ptv_um = np.asarray([row["peak_to_valley_um"] for row in rows])
    one_cell_um = dx_m * 1e6

    ax = axes[0]
    positive = rms_um > 0
    ax.plot(time_ms[positive], rms_um[positive], "o-", color=OKABE_ITO[0],
            label="RMS amplitude")
    ax.plot(time_ms[positive], ptv_um[positive], "s--", color=OKABE_ITO[1],
            label="peak-to-valley")
    ax.axhline(one_cell_um, color="0.4", ls=":", lw=1.0)
    ax.text(0.02, one_cell_um * 1.15, "one cell", transform=ax.get_yaxis_transform(),
            fontsize=7, color="0.4")
    if np.isfinite(sigma) and positive.sum() >= 3:
        t_fit = time_ms[positive]
        reference = rms_um[positive][0] * np.exp(sigma * (t_fit - t_fit[0]) * 1e-3)
        ax.plot(t_fit, reference, "-", color="0.5", lw=1.0, zorder=0,
                label=rf"fit $\sigma$ = {sigma:.3g} s$^{{-1}}$")
    ax.set_yscale("log")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(r"Interface amplitude ($\mu$m)")
    ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
    ax.legend(frameon=False, loc="best")

    ax = axes[1]
    nx = rows[0]["_deviation"].size
    x_um = np.arange(nx) * dx_m * 1e6
    for row, color, style in ((rows[0], OKABE_ITO[2], "-"), (rows[-1], OKABE_ITO[3], "-")):
        ax.plot(x_um, row["_deviation"] * dx_m * 1e6, style, color=color, lw=0.8,
                label=f"{row['time_s'] * 1e3:.2f} ms (step {row['step']:,})")
    ax.set_xlabel(r"$x$ ($\mu$m)")
    ax.set_ylabel(r"$h - \langle h \rangle$ ($\mu$m)")
    ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
    ax.legend(frameon=False, loc="best")

    for index, ax in enumerate(axes):
        ax.text(-0.16, 1.04, f"({chr(97 + index)})", transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")

    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png"):
        target = output_dir / f"interface_amplitude.{suffix}"
        fig.savefig(target, dpi=300)
        print(f"Saved: {target}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frames", nargs="+", type=Path)
    parser.add_argument("--dx", type=float, required=True, help="cell size (m)")
    parser.add_argument("--dt", type=float, required=True, help="timestep (s)")
    parser.add_argument("--substrate-top", type=int, required=True,
                        help="top y cell of the seeded flat front")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--no-plot", action="store_true")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--shift-file", type=Path,
                       help="DATA/shift.dat; required once the window has fired")
    group.add_argument("--no-shift", action="store_true",
                       help="assert the window has not shifted over these frames")
    args = parser.parse_args()

    paths = sorted(args.frames, key=timestep)

    # The window shift is not recoverable from the frames themselves, so a
    # missing shift.dat cannot be allowed to pass as "no shift" by default.
    if args.no_shift:
        shift_map: dict[int, int] = {}
    else:
        if not args.shift_file.is_file():
            raise SystemExit(
                f"--shift-file {args.shift_file} does not exist. The solver writes it "
                "only once the window fires; pass --no-shift to assert it has not.")
        shift_map = read_shift_map(args.shift_file)
        # Steps before the file's first entry predate the window and are
        # genuinely unshifted. A gap at or after it means the file is holed or
        # stops short of the frames, and the heights would be under-reported.
        if shift_map:
            first = min(shift_map)
            missing = [timestep(p) for p in paths
                       if timestep(p) > first and timestep(p) not in shift_map]
            if missing:
                raise SystemExit(
                    f"{args.shift_file} covers steps {first}..{max(shift_map)} but has "
                    f"no entry for {missing[:5]}{' ...' if len(missing) > 5 else ''}; "
                    "it does not cover these frames.")

    shifts = {timestep(p): shift_at(shift_map, timestep(p)) for p in paths}
    moved = max(shifts.values()) if shifts else 0
    if args.no_shift and args.shift_file is None:
        candidate = paths[0].parent.parent / "shift.dat"
        if candidate.is_file() and max(read_shift_map(candidate).values(), default=0) > 0:
            raise SystemExit(
                f"--no-shift was asserted but {candidate} records a nonzero shift. "
                "Re-run with --shift-file; the heights would be wrong.")
    print(f"window shift over these frames: 0 -> {moved} cells")

    rows = [analyse(path, args.dx, args.dt, args.substrate_top, shifts[timestep(path)])
            for path in paths]
    add_front_velocities(rows, args.dx)

    times = np.asarray([row["time_s"] for row in rows])
    rms = np.asarray([row["rms_amplitude_cells"] for row in rows])
    sigma, r_squared, n_fit = growth_rate(times, rms)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    columns = [key for key in rows[0] if not key.startswith("_")]
    csv_path = args.output_dir / "interface_amplitude.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {csv_path}")

    header = ("step", "t (ms)", "shift", "<h>-sub", "h_tip", "v_mean", "v_tip",
              "RMS (um)", "PtV (um)", "n_prot", "lam_dom")
    units = ("", "", "cell", "cell", "cell", "mm/s", "mm/s", "", "", "", "um")
    fmt = "{:>9} {:>8} {:>7} {:>9} {:>8} {:>8} {:>8} {:>10} {:>10} {:>7} {:>9}"
    print("\n" + fmt.format(*header))
    print(fmt.format(*units))
    for row in rows:
        print("{:>9} {:>8.3f} {:>7} {:>9.1f} {:>8.0f} {:>8.3f} {:>8.3f} "
              "{:>10.5f} {:>10.5f} {:>7} {:>9.3f}".format(
                  row["step"], row["time_s"] * 1e3, row["cumulative_shift_cells"],
                  row["h_above_substrate_cells"], row["h_tip_cells"],
                  row["v_mean_mm_s"], row["v_tip_mm_s"],
                  row["rms_amplitude_um"], row["peak_to_valley_um"],
                  row["protrusion_count"], row["dominant_wavelength_um"]))

    one_cell_um = args.dx * 1e6
    last = rows[-1]
    print(f"\none cell = {one_cell_um:.5f} um; domain width = "
          f"{last['_deviation'].size * one_cell_um:.3f} um")
    if np.isfinite(sigma):
        print(f"growth rate sigma = {sigma:.4g} 1/s over {n_fit} frames "
              f"(R^2 = {r_squared:.4f})")
    else:
        print(f"growth rate: INSUFFICIENT DATA ({n_fit} frames with RMS > 0)")

    verdict = ("DESTABILISING" if last["rms_amplitude_cells"] > 1.0
               else "still flat at sub-cell amplitude")
    print(f"VERDICT: final RMS = {last['rms_amplitude_cells']:.4f} cells "
          f"({last['rms_amplitude_um']:.5f} um) -> {verdict}")

    if not args.no_plot:
        make_figure(rows, args.dx, sigma, args.output_dir)


if __name__ == "__main__":
    main()
