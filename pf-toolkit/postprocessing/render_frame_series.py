"""Render a fixed-scale phase / composition frame series for one MicroSim case.

Designed to run ON THE HPCC so only the finished images are transferred.
Produces, in --output-dir:
  frames/frame_<idx>_step<step>.png   two-panel phase + composition per frame
  <prefix>_phase_evolution.gif        animation built from those panels
  <prefix>_phase_montage.png          single-sheet contact print of the phase field

Colour scales are fixed across every frame (phi in [0,1]; xSi over the global
range of the whole series), so panels are directly comparable frame to frame.

The vertical axis is WINDOW-LOCAL. The moving window shifts during the run, so
lab-frame y = window y + (cumulative shift). Each panel is annotated with its
own cumulative shift.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analyze_pilot_morphology import read_fragment, timestep  # noqa: E402

MM_TO_IN = 1 / 25.4
DOUBLE_COL_IN = 180 * MM_TO_IN


def apply_style() -> None:
    """SciencePlots 'science' if available, else an inline equivalent.

    The HPCC env has neither scienceplots nor LaTeX, so the fallback keeps the
    remote figures consistent with the locally styled ones apart from the font.
    """
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "no-latex"])
    except Exception:
        mpl.rcParams.update({
            "font.family": "serif",
            "mathtext.fontset": "dejavuserif",
            "axes.grid": False,
            "axes.prop_cycle": mpl.cycler(
                "color", ["#0C5DA5", "#00B945", "#FF9500", "#FF2C00", "#845B97"]
            ),
        })
    mpl.rcParams.update({
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.top": False, "ytick.right": False,
        "savefig.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "figure.dpi": 120, "text.usetex": False,
    })


def time_display(steps: list[int], dt: float) -> tuple[float, str, int]:
    """Pick a time unit and decimal count so frame-to-frame changes are visible.

    These runs are millisecond-scale (dt ~ 2e-9 s over ~1e6 steps), so the old
    fixed ".2f" in seconds printed every panel as "0.00 s". Scale to whichever
    unit puts the last frame at O(1) and keep enough decimals to resolve one
    save interval.
    """
    t_max = max(steps) * dt
    scale, label = 1e9, "ns"
    for candidate, name in ((1.0, "s"), (1e3, "ms"), (1e6, r"$\mu$s"), (1e9, "ns")):
        if t_max * candidate >= 1.0:
            scale, label = candidate, name
            break
    gaps = {b - a for a, b in zip(steps, steps[1:]) if b > a}
    step_gap = min(gaps) if gaps else (max(steps) or 1)
    gap = step_gap * dt * scale
    decimals = 2 if gap <= 0 else min(4, max(0, int(np.ceil(-np.log10(gap))) + 1))
    return scale, label, decimals


def read_shift_map(path: Path) -> dict[int, int]:
    """step -> cumulative shift in cells. Absent steps (e.g. 0) mean no shift."""
    mapping: dict[int, int] = {}
    if not path.is_file():
        return mapping
    for line in path.read_text().split("\n"):
        parts = line.split()
        if len(parts) == 2:
            mapping[int(parts[0])] = int(parts[1])
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path, help="DATA/Processor_0 of the case")
    parser.add_argument("--shift-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="case")
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--dx-um", type=float, required=True)
    parser.add_argument("--crop-y-um", type=float, default=105.0)
    parser.add_argument("--xsi-far-field", type=float, default=0.052952583584)
    parser.add_argument("--gif-ms", type=int, default=260)
    args = parser.parse_args()

    apply_style()
    frame_dir = args.output_dir / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted(args.input_dir.glob("*.vtk"), key=timestep)
    if not paths:
        raise SystemExit(f"no VTK fragments in {args.input_dir}")
    shift_map = read_shift_map(args.shift_file)

    # Pass 1: load once, keep cropped fields, and find the global xSi range so
    # every frame shares one colour scale.
    phis: list[np.ndarray] = []
    xsis: list[np.ndarray] = []
    steps: list[int] = []
    crop = int(round(args.crop_y_um / args.dx_um))
    for path in paths:
        dims, arrays = read_fragment(path)
        if dims[2] != 1:
            raise ValueError(f"{path}: expected a 2-D fragment")
        # stored [x, y]; imshow wants [row, col] = [y, x]
        phis.append(arrays["FCC_A1"][:, :, 0].T[:crop])
        xsis.append(arrays["Composition_Si"][:, :, 0].T[:crop])
        steps.append(timestep(path))
    print(f"loaded {len(paths)} frames, steps {steps[0]}..{steps[-1]}")

    xsi_lo = float(min(f.min() for f in xsis))
    xsi_hi = float(max(f.max() for f in xsis))
    print(f"global xSi range {xsi_lo:.6f} .. {xsi_hi:.6f}")

    tscale, tunit, tdec = time_display(steps, args.dt)
    print(f"time axis in {tunit} with {tdec} decimals "
          f"(last frame {steps[-1] * args.dt * tscale:.{tdec}f} {tunit})")

    ny, nx = phis[0].shape
    lx_um = nx * args.dx_um
    ly_um = ny * args.dx_um
    extent = [0.0, lx_um, 0.0, ly_um]

    panel_w_in = (DOUBLE_COL_IN - 0.85) / 2.0
    fig_h_in = panel_w_in / (lx_um / ly_um) + 0.62

    written: list[Path] = []
    for index, (phi, xsi, step) in enumerate(zip(phis, xsis, steps)):
        shift_cells = shift_map.get(step, 0)
        time_s = step * args.dt
        fig, axes = plt.subplots(
            1, 2, figsize=(DOUBLE_COL_IN, fig_h_in), constrained_layout=True
        )

        im0 = axes[0].imshow(
            phi, origin="lower", extent=extent, aspect="equal",
            cmap="viridis", vmin=0.0, vmax=1.0, interpolation="nearest",
        )
        cb0 = fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.02, aspect=14)
        cb0.set_label(r"FCC-Al phase fraction $\phi$")

        im1 = axes[1].imshow(
            xsi, origin="lower", extent=extent, aspect="equal",
            cmap="plasma", vmin=xsi_lo, vmax=xsi_hi, interpolation="nearest",
        )
        axes[1].contour(
            phi, levels=[0.5], colors="white", linewidths=0.5,
            extent=extent, origin="lower",
        )
        cb1 = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.02, aspect=14)
        cb1.set_label(r"Si mole fraction $x_\mathrm{Si}$")
        cb1.ax.axhline(args.xsi_far_field, color="white", linewidth=0.9)

        for label, ax in zip("ab", axes):
            ax.set_xlabel(r"$x$ ($\mu$m)")
            ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
            ax.text(
                0.022, 0.94, f"({label})", transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="top", color="white",
                bbox=dict(boxstyle="round,pad=0.15", fc="black", ec="none", alpha=0.55),
            )
        axes[0].set_ylabel(r"growth direction $y$ ($\mu$m)")

        fig.suptitle(
            f"step {step:,}   $t$ = {time_s * tscale:.{tdec}f} {tunit}   "
            f"window shifted {shift_cells * args.dx_um:.1f} " + r"$\mu$m",
            fontsize=9,
        )
        out = frame_dir / f"frame_{index:03d}_step{step:08d}.png"
        fig.savefig(out)
        plt.close(fig)
        written.append(out)
    print(f"wrote {len(written)} frame images to {frame_dir}")

    # Contact sheet of the phase field only.
    ncols = 6
    nrows = int(np.ceil(len(phis) / ncols))
    cell_w = DOUBLE_COL_IN / ncols
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(DOUBLE_COL_IN, nrows * cell_w / (lx_um / ly_um) + 0.4),
        constrained_layout=True,
    )
    for ax, phi, step in zip(axes.flat, phis, steps):
        ax.imshow(
            phi, origin="lower", extent=extent, aspect="equal",
            cmap="viridis", vmin=0.0, vmax=1.0, interpolation="nearest",
        )
        ax.set_title(f"{step * args.dt * tscale:.{tdec}f} {tunit}",
                     fontsize=6.5, pad=1.5)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in list(axes.flat)[len(phis):]:
        ax.axis("off")
    montage = args.output_dir / f"{args.prefix}_phase_montage.png"
    fig.savefig(montage, dpi=300)
    plt.close(fig)
    print(f"wrote {montage}")

    # Animation from the per-frame panels.
    try:
        from PIL import Image
        images = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in written]
        gif = args.output_dir / f"{args.prefix}_phase_evolution.gif"
        images[0].save(
            gif, save_all=True, append_images=images[1:],
            duration=args.gif_ms, loop=0,
        )
        print(f"wrote {gif}")
    except Exception as exc:  # pillow missing or palette failure
        print(f"GIF skipped: {exc}")


if __name__ == "__main__":
    main()
