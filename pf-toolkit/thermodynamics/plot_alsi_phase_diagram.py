#!/usr/bin/env python3
"""Al-Si binary phase diagram computed from the project's own COST507 TDB.

Runs REMOTELY in the `microsim-pp` env (pycalphad 0.11.2 -- the same version that
produced the frozen design constants); only the figure and CSVs are
pulled back. There is no LaTeX on the HPCC, so SciencePlots runs in `no-latex`.

Panel (a) the full binary, panel (b) the Al-rich region carrying this project's
operating point (5.5 wt.% Si) and the LINEARISED liquidus actually used to derive
m_L, k and dT0, so the linearisation error is visible rather than assumed.

All 25 phases that COST507 can populate from {AL, SI, VA} are offered to the
mapper; it returns only LIQUID, FCC_A1 and DIAMOND_A4 as stable, which is the
check that no spurious multicomponent phase leaks into the binary section.

Usage:
    python plot_alsi_phase_diagram.py --tdb <path> --outdir <dir>
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401

from pycalphad import Database, variables as v
from pycalphad.core.utils import filter_phases, unpack_species
from pycalphad.mapping import BinaryStrategy

# ----------------------------------------------------------------- constants
M_AL = 26.9815  # g/mol, TDB ELEMENT line
M_SI = 28.085   # g/mol, TDB ELEMENT line

# Frozen design values -- plotted here, never recomputed.
C0_WT = 5.5
T_LIQ_C0 = 899.604509   # K
M_L_LIN = -6.507        # K/wt%
K_PART = 0.10903        # -
T_EUT_USED = 850.0      # K, the isotherm below which this study claims nothing

MM_TO_IN = 1 / 25.4
DOUBLE_COL_IN = 180 * MM_TO_IN
OKABE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#0072B2", "#D55E00", "#CC79A7", "#000000"]


def x_to_wt(x_si):
    """Si mole fraction -> Si wt.%."""
    x_si = np.asarray(x_si, dtype=float)
    return 100.0 * x_si * M_SI / (x_si * M_SI + (1.0 - x_si) * M_AL)


def wt_to_x(w_si):
    """Si wt.% -> Si mole fraction."""
    w = np.asarray(w_si, dtype=float) / 100.0
    return (w / M_SI) / (w / M_SI + (1.0 - w) / M_AL)


def apply_style():
    plt.style.use(["science", "no-latex"])
    mpl.rcParams.update({
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 7, "legend.title_fontsize": 8,
        "lines.linewidth": 1.2, "lines.markersize": 4, "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "figure.dpi": 150, "text.usetex": False,
    })


# ------------------------------------------------------------------ mapping
def run_map(db, comps, phases, t_range, x_range):
    strat = BinaryStrategy(
        db, comps, phases,
        {v.T: t_range, v.X("SI"): x_range, v.P: 101325, v.N: 1},
    )
    strat.do_map()
    return strat


def extract_boundaries(strat):
    """-> list of segments, each {'phases': frozenset, 'data': {phase: (wt%, T)}}.

    Keying on the segment's own phase pair matters: inside the Al-rich window
    BOTH liquidus limbs (L+FCC_A1 and L+DIAMOND_A4) sit below 20 wt.% Si, so a
    composition filter cannot separate them and merging them corrupts any
    liquidus slope or partition coefficient taken from the result.
    """
    segments = []
    for seg in strat.get_tieline_data(v.X("SI"), v.T):
        data = {}
        for spd in seg.data:
            x = np.atleast_1d(np.asarray(spd.x, dtype=float))
            t = np.atleast_1d(np.asarray(spd.y, dtype=float))
            if x.size < 2:
                continue
            order = np.argsort(t)
            data[spd.phase] = (x_to_wt(x[order]), t[order])
        if data:
            segments.append({"phases": frozenset(data), "data": data})
    return segments


def eutectic_from(strat):
    """-> (T_eut, {phase: wt% Si}) from the mapper's invariant nodes."""
    inv = strat.get_invariant_data(v.X("SI"), v.T)
    if not inv:
        raise SystemExit("no invariant found - widen the mapped T range")
    # Several near-identical nodes are returned for one invariant; average them.
    temps, comps = [], {}
    for node in inv:
        for spd in node.data:
            temps.append(float(np.atleast_1d(spd.y)[0]))
            comps.setdefault(spd.phase, []).append(float(np.atleast_1d(spd.x)[0]))
    return float(np.mean(temps)), {p: x_to_wt(np.mean(c)) for p, c in comps.items()}


LIQ_FCC = frozenset({"LIQUID", "FCC_A1"})


def liquidus_solidus(segments, t_eut, t_hi):
    """Al-rich liquidus and solidus, taken ONLY from L+FCC_A1 segments."""
    curves = {}
    for phase in ("LIQUID", "FCC_A1"):
        xs, ts = [], []
        for seg in segments:
            if seg["phases"] != LIQ_FCC:
                continue
            w, t = seg["data"][phase]
            m = (t >= t_eut - 0.5) & (t <= t_hi)
            if m.sum() >= 2:
                xs.append(w[m])
                ts.append(t[m])
        if not xs:
            continue
        w = np.concatenate(xs)
        t = np.concatenate(ts)
        o = np.argsort(t)
        t, w = t[o], w[o]
        keep = np.concatenate(([True], np.diff(t) > 1e-9))
        curves[phase] = (w[keep], t[keep])
    return curves


# ------------------------------------------------------------------- figure
def make_figure(bounds_full, bounds_zoom, t_eut, comp_eut, stats, outdir):
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL_IN, 82 * MM_TO_IN), constrained_layout=True
    )

    # ---- (a) full binary -------------------------------------------------
    for seg in bounds_full:
        for w, t in seg["data"].values():
            ax_a.plot(w, t, "-", color="k", lw=1.0)
    w_fcc = comp_eut.get("FCC_A1", 1.55)
    ax_a.hlines(t_eut, w_fcc, 100.0, colors="k", lw=1.0)
    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(700, 1750)
    ax_a.set_xlabel("Si content (wt.\\%)" if mpl.rcParams["text.usetex"]
                    else "Si content (wt.%)")
    ax_a.set_ylabel("Temperature (K)")
    for x, y, s in [(32, 1480, "L"), (52, 1080, "L + Si(D)"),
                    (58, 770, r"$\alpha$ + Si(D)")]:
        ax_a.text(x, y, s, fontsize=8, ha="center", va="center")
    # The alpha and L+alpha fields are too narrow to label in place.
    for xy, xytext, s in [((0.8, 800), (14, 790), r"$\alpha$"),
                          ((4.0, 890), (20, 960), r"L + $\alpha$")]:
        ax_a.annotate(s, xy=xy, xytext=xytext, fontsize=8, ha="center",
                      va="center",
                      arrowprops=dict(arrowstyle="-", lw=0.5, color="k",
                                      shrinkA=2, shrinkB=1))
    ax_a.annotate(
        f"eutectic {t_eut:.1f} K, {comp_eut['LIQUID']:.2f} wt.% Si",
        xy=(comp_eut["LIQUID"], t_eut), xytext=(30, 900),
        fontsize=7, ha="left", va="center",
        arrowprops=dict(arrowstyle="-", lw=0.5, color="k",
                        shrinkA=2, shrinkB=1),
    )

    # ---- (b) Al-rich zoom ------------------------------------------------
    for seg in bounds_zoom:
        for w, t in seg["data"].values():
            ax_b.plot(w, t, "-", color="k", lw=1.0)
    ax_b.hlines(t_eut, w_fcc, 20.0, colors="k", lw=1.0)

    w_lin = np.linspace(0, 14, 200)
    t_lin = T_LIQ_C0 + M_L_LIN * (w_lin - C0_WT)
    ax_b.plot(w_lin, t_lin, "--", color=OKABE[5], lw=1.2,
              label=f"linearised liquidus\n$m_L$ = {M_L_LIN} K/wt.%")

    ax_b.axvline(C0_WT, color=OKABE[4], ls=":", lw=1.0)
    ax_b.plot([C0_WT], [T_LIQ_C0], "o", color=OKABE[4], ms=4.5, zorder=5,
              label=f"$c_0$ = {C0_WT} wt.% Si\n$T_{{liq}}$ = {T_LIQ_C0:.2f} K")
    ax_b.axhline(T_EUT_USED, color=OKABE[6], ls="-.", lw=1.0,
                 label=f"$T_{{eut}}$ = {T_EUT_USED:.0f} K (study floor)")

    ax_b.set_xlim(0, 20)
    ax_b.set_ylim(835, 940)
    ax_b.set_xlabel("Si content (wt.%)")
    ax_b.set_ylabel("Temperature (K)")
    ax_b.text(0.6, 858, r"$\alpha$", fontsize=8, ha="center", va="center")
    ax_b.text(11.5, 895, r"L", fontsize=8, ha="center", va="center")
    ax_b.text(7.6, 861, r"L + $\alpha$", fontsize=8, ha="center", va="center")
    ax_b.text(17.0, 843, r"$\alpha$ + Si(D)", fontsize=8, ha="center", va="center")
    if stats is not None:
        ax_b.text(0.4, 836.5,
                  f"linearisation error $\\leq$ {stats['max_dev']:.1f} K over "
                  f"$T_{{eut}}$–$T_{{liq}}(c_0)$\n"
                  f"$k$ from TDB: {stats['k_at_liq']:.4f} at $T_{{liq}}(c_0)$, "
                  f"{stats['k_hi']:.4f} at $T_{{eut}}$",
                  fontsize=6.5, ha="left", va="bottom")
    ax_b.legend(frameon=False, loc="upper right", handlelength=1.6,
                borderaxespad=0.4)

    for ax in (ax_a, ax_b):
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
        sec = ax.secondary_yaxis(
            "right", functions=(lambda t: t - 273.15, lambda c: c + 273.15))
        sec.set_ylabel("Temperature (°C)")
        sec.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
        ax.tick_params(right=False)

    for i, ax in enumerate((ax_a, ax_b)):
        ax.text(-0.16, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")

    for fmt in ("pdf", "png"):
        p = outdir / f"alsi_phase_diagram_COST507.{fmt}"
        fig.savefig(p, dpi=600)
        print(f"saved {p}")
    plt.close(fig)


# ------------------------------------------------- focused single-panel figure
def make_zoom_figure(bounds, t_eut, comp_eut, stats, outdir,
                     xmax=20.0, tmin=800.0, tmax=1000.0,
                     width_mm=110.0, height_mm=88.0):
    """Standalone Al-rich diagram over 0..xmax wt.% Si and tmin..tmax K."""
    fig, ax = plt.subplots(figsize=(width_mm * MM_TO_IN, height_mm * MM_TO_IN),
                           constrained_layout=True)

    for seg in bounds:
        for w, t in seg["data"].values():
            ax.plot(w, t, "-", color="k", lw=1.0)
    ax.hlines(t_eut, comp_eut.get("FCC_A1", 1.55), xmax, colors="k", lw=1.0)

    w_lin = np.linspace(0, 14, 200)
    ax.plot(w_lin, T_LIQ_C0 + M_L_LIN * (w_lin - C0_WT), "--", color=OKABE[5],
            lw=1.2, label=f"linearised liquidus, $m_L$ = {M_L_LIN} K/wt.%")
    # Clipped below the legend box so the guide line does not cross its text.
    ax.axvline(C0_WT, ymin=0.0, ymax=0.84, color=OKABE[4], ls=":", lw=1.0)
    ax.plot([C0_WT], [T_LIQ_C0], "o", color=OKABE[4], ms=4.5, zorder=5,
            label=f"$c_0$ = {C0_WT} wt.% Si, $T_{{liq}}$ = {T_LIQ_C0:.2f} K")
    ax.axhline(T_EUT_USED, color=OKABE[6], ls="-.", lw=1.0,
               label=f"$T_{{eut}}$ = {T_EUT_USED:.0f} K (study floor)")

    ax.set_xlim(0, xmax)
    ax.set_ylim(tmin, tmax)
    ax.set_xlabel("Si content (wt.%)")
    ax.set_ylabel("Temperature (K)")

    for x, y, s_ in [(14.0, 968, "L"), (6.0, 878, r"L + $\alpha$"),
                     (18.3, 878, "L + Si(D)"),
                     (9.6, 822, r"$\alpha$ + Si(D)")]:
        ax.text(x, y, s_, fontsize=8, ha="center", va="center")
    # The alpha field is ~1 wt.% wide; label it in place rather than with a
    # leader, which would have to terminate in a neighbouring two-phase field.
    ax.text(0.48, 822, r"$\alpha$", fontsize=7, ha="center", va="center")
    ax.annotate(f"eutectic {t_eut:.2f} K\n{comp_eut['LIQUID']:.2f} wt.% Si",
                xy=(comp_eut["LIQUID"], t_eut), xytext=(15.2, 838),
                fontsize=7, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="k",
                                shrinkA=2, shrinkB=1))
    if stats is not None:
        ax.text(0.98, 0.012,
                f"linearisation error $\\leq$ {stats['max_dev']:.1f} K over "
                f"$T_{{eut}}$–$T_{{liq}}(c_0)$\n"
                f"$k$ from TDB: {stats['k_at_liq']:.4f} at $T_{{liq}}(c_0)$ "
                f"$\\rightarrow$ {stats['k_hi']:.4f} at $T_{{eut}}$",
                transform=ax.transAxes, fontsize=6.5, ha="right", va="bottom",
                linespacing=1.4)

    ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
    sec = ax.secondary_yaxis("right",
                             functions=(lambda t: t - 273.15,
                                        lambda c: c + 273.15))
    sec.set_ylabel("Temperature (°C)")
    sec.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
    ax.tick_params(right=False)
    ax.legend(frameon=False, loc="upper left", handlelength=1.8,
              borderaxespad=0.5)

    for fmt in ("pdf", "png"):
        pth = outdir / f"alsi_phase_diagram_alrich.{fmt}"
        fig.savefig(pth, dpi=600)
        print(f"saved {pth}")
    plt.close(fig)


# ------------------------------------------------------- linearisation check
def validate_linearisation(curves, t_eut, outdir):
    """Quantify how far the frozen linear m_L and constant k drift from the TDB."""
    if "LIQUID" not in curves or "FCC_A1" not in curves:
        print("WARNING: could not isolate liquidus/solidus; skipping validation")
        return None
    w_l, t_l = curves["LIQUID"]
    w_s, t_s = curves["FCC_A1"]
    t_grid = np.linspace(max(t_eut, t_l.min(), t_s.min()) + 0.05,
                         min(t_l.max(), t_s.max()) - 0.05, 400)
    wl = np.interp(t_grid, t_l, w_l)
    ws = np.interp(t_grid, t_s, w_s)

    t_lin = T_LIQ_C0 + M_L_LIN * (wl - C0_WT)
    dev = t_lin - t_grid
    k_local = wt_to_x(ws) / wt_to_x(wl)
    m_local = np.gradient(t_grid, wl)

    hdr = ("T_K,w_liquidus_wtpct,w_solidus_wtpct,x_liquidus_molfrac,"
           "x_solidus_molfrac,k_local,m_L_local_K_per_wtpct,"
           "T_linearised_K,linearisation_error_K")
    arr = np.column_stack([t_grid, wl, ws, wt_to_x(wl), wt_to_x(ws),
                           k_local, m_local, t_lin, dev])
    p = outdir / "alsi_liquidus_solidus_COST507.csv"
    np.savetxt(p, arr, delimiter=",", header=hdr, comments="", fmt="%.8g")
    print(f"saved {p}")

    def at(t):
        i = int(np.argmin(np.abs(t_grid - t)))
        return k_local[i], m_local[i], dev[i], wl[i]

    print("\n--- linearisation check against the TDB "
          "(frozen: m_L = %.3f K/wt%%, k = %.5f) ---" % (M_L_LIN, K_PART))
    for label, t in [("at T_liq(c0)", T_LIQ_C0), ("at T_eut", t_eut + 0.1)]:
        k_, m_, d_, wl_ = at(t)
        print(f"  {label:14s} T = {t:8.3f} K  w_L = {wl_:6.3f} wt%  "
              f"k = {k_:.5f}  m_L = {m_:7.3f} K/wt%  "
              f"linearised liquidus off by {d_:+.2f} K")
    win = t_grid <= T_LIQ_C0
    max_dev = float(np.nanmax(np.abs(dev[win])))
    print(f"  max |linearisation error| over T_eut..T_liq(c0): {max_dev:.2f} K")
    print(f"  k range over the same window: "
          f"{np.nanmin(k_local[win]):.5f} .. {np.nanmax(k_local[win]):.5f}")
    return {"max_dev": max_dev,
            "k_lo": float(np.nanmin(k_local[win])),
            "k_hi": float(np.nanmax(k_local[win])),
            "k_at_liq": float(at(T_LIQ_C0)[0]),
            "m_at_liq": float(at(T_LIQ_C0)[1])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tdb", required=True)
    ap.add_argument("--outdir", default=".", type=Path)
    ap.add_argument("--xmax-wt", type=float, default=20.0,
                    help="right-hand limit of the focused panel, wt.% Si")
    ap.add_argument("--tmin", type=float, default=800.0)
    ap.add_argument("--tmax", type=float, default=1000.0)
    ap.add_argument("--skip-full", action="store_true",
                    help="skip the full-binary map and two-panel figure")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    db = Database(args.tdb)
    comps = ["AL", "SI", "VA"]
    phases = sorted(p for p in filter_phases(db, unpack_species(db, comps))
                    if p != "GAS")
    print(f"{len(phases)} candidate phases offered to the mapper")

    print("mapping full binary ...")
    s_full = run_map(db, comps, phases, (700, 1800, 10), (0, 1, 0.005))
    stable = s_full.get_all_phases()
    print("stable phases returned:", stable)
    if set(stable) != {"LIQUID", "FCC_A1", "DIAMOND_A4"}:
        print("WARNING: unexpected stable phase set for the Al-Si binary")

    t_eut, comp_eut = eutectic_from(s_full)
    print(f"\nEUTECTIC from this TDB: T = {t_eut:.3f} K ({t_eut - 273.15:.2f} °C)")
    for p in sorted(comp_eut):
        print(f"  {p:12s} {comp_eut[p]:8.4f} wt.% Si "
              f"({wt_to_x(comp_eut[p]):.6f} mol frac)")

    x_hi = min(0.999, float(wt_to_x(args.xmax_wt)) * 1.35)
    print(f"\nmapping Al-rich region: T {args.tmin}-{args.tmax} K, "
          f"x(Si) 0-{x_hi:.3f} ...")
    s_zoom = run_map(db, comps, phases,
                     (args.tmin - 5.0, args.tmax + 5.0, 2),
                     (0, x_hi, 0.002))

    bounds_full = extract_boundaries(s_full)
    bounds_zoom = extract_boundaries(s_zoom)
    curves = liquidus_solidus(bounds_zoom, t_eut, args.tmax + 5.0)

    stats = validate_linearisation(curves, t_eut, args.outdir)
    apply_style()
    make_zoom_figure(bounds_zoom, t_eut, comp_eut, stats, args.outdir,
                     xmax=args.xmax_wt, tmin=args.tmin, tmax=args.tmax)
    if not args.skip_full:
        make_figure(bounds_full, bounds_zoom, t_eut, comp_eut, stats,
                    args.outdir)

    with open(args.outdir / "alsi_phase_diagram_provenance.txt", "w") as fh:
        fh.write(f"tdb={args.tdb}\n")
        fh.write(f"candidate_phases={phases}\n")
        fh.write(f"stable_phases={stable}\n")
        fh.write(f"T_eutectic_K={t_eut:.6f}\n")
        for p in sorted(comp_eut):
            fh.write(f"x_eut_{p}_wtpct={comp_eut[p]:.6f}\n")


if __name__ == "__main__":
    main()
