#!/usr/bin/env python3
"""
Read one dT_tip shooting shot and say which way to move dT_init.

THE MEASUREMENT.  In the frozen-temperature approximation the tip undercooling
obeys

    d(dT_tip)/dt = G * (V - v_tip)

because the isotherms sweep past at V while the tip advances at v_tip.  So
`v_tip = V` and `dT_tip stationary` are the SAME condition, and the shot does
not need to reach steady state -- it only needs to show which side of it we are
on.  Integrating over the whole shot gives the primary verdict:

    delta_dT_tip = G * (V * t_shot - delta_front * dx)

This uses only the net front displacement, so the integer quantisation of
tip_row enters once over the whole run instead of once per sample.  It is far
less noisy than differencing v_tip, which is why it is the headline number and
the per-sample velocity is only a cross-check.

    delta_dT_tip > 0  ->  the tip got COLDER: v_tip < V, banking resumed,
                          dT_init was too WARM.  INCREASE dT_init.
    delta_dT_tip < 0  ->  the tip got WARMER: v_tip > V, undercooling unwound,
                          dT_init was too COLD.  DECREASE dT_init.
    |delta_dT_tip| small -> dT_init IS dT_tip.  That is deliverable (3).

Bias the bracket cold: overshooting cold is self-correcting (the front outruns
the isotherm and warms), overshooting warm re-enters banking, which is the
failure that cost two legs.

Reads DATA/farfield_guard.csv, which the solver writes every SAVET.
"""

import argparse
import csv
import math
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv", help="DATA/farfield_guard.csv from the shot")
    p.add_argument("--dt-init", type=float, required=True,
                   help="tip undercooling imposed at step 0, K")
    p.add_argument("--gradient", type=float, default=1.1e6, help="K/m")
    p.add_argument("--velocity", type=float, default=1.3333e-2, help="m/s")
    p.add_argument("--dx", type=float, default=1.0771e-8, help="m")
    p.add_argument("--dt", type=float, default=2.0e-9, help="s per step")
    p.add_argument("--skip-frac", type=float, default=0.34,
                   help="ignore this fraction of the shot at the start; the "
                        "seed's grooves remelt when the field is reset warmer "
                        "and that relaxation is not the signal")
    p.add_argument("--tol-v", type=float, default=0.02,
                   help="|1 - v_mean/V| below this counts as stationary. The "
                        "criterion is RELATIVE, not an absolute delta_dT_tip: "
                        "over a 190k-step window the entire measurable range is "
                        "only G*V*t = 3.5 K, so an absolute 1 K tolerance would "
                        "pass a front running at 0.8 V. 2% maps to ~1.1 K of "
                        "dT_tip via G*(V-v)*tau.")
    p.add_argument("--tau", type=float, default=3.855e-3,
                   help="interface response time, s, for the suggested dT_init "
                        "correction. Fitted on a PLANAR incubating front; a "
                        "developed cellular front should respond faster, so "
                        "treat the correction as an upper bound on the step.")
    args = p.parse_args()

    rows = []
    with open(args.csv) as fh:
        for r in csv.DictReader(fh):
            try:
                rows.append((int(r["step"]), int(r["front_abs_cells"]),
                             float(r["topband_dc_rel_c0"]),
                             float(r["v_tip_over_V"])))
            except (KeyError, ValueError):
                continue
    if len(rows) < 3:
        sys.exit(f"{args.csv}: need >= 3 usable rows, got {len(rows)}. "
                 f"Old-format CSV without the short-box columns?")

    rows.sort()
    i0 = int(len(rows) * args.skip_frac)
    seg = rows[i0:]
    step0, front0 = seg[0][0], seg[0][1]
    step1, front1 = seg[-1][0], seg[-1][1]

    t_shot = (step1 - step0) * args.dt
    d_front = (front1 - front0) * args.dx
    v_mean = d_front / t_shot if t_shot > 0 else float("nan")
    ddT = args.gradient * (args.velocity * t_shot - d_front)
    # One cell of tip_row quantisation, at each end of the segment.
    ddT_err = args.gradient * args.dx * 2.0

    dc = [r[2] for r in seg]
    vs = [r[3] for r in seg if not math.isnan(r[3])]

    print(f"shot from {args.csv}")
    print(f"  dT_init imposed          {args.dt_init:8.3f} K")
    print(f"  window analysed          steps {step0} -> {step1} "
          f"({len(seg)} samples, first {i0} skipped)")
    print(f"  front advance            {front1-front0} cells "
          f"= {d_front*1e6:.3f} um in {t_shot*1e3:.4f} ms")
    print(f"  mean v_tip               {v_mean*1e3:8.4f} mm/s "
          f"= {v_mean/args.velocity:.3f} x V")
    if vs:
        print(f"  per-sample v_tip/V       min {min(vs):.2f}  "
              f"max {max(vs):.2f}  (cross-check only)")
    print(f"  top-band dC/c0           {min(dc)*100:+.3f}% -> {max(dc)*100:+.3f}%"
          f"   (containment; must not climb above +2%)")
    print()
    print(f"  delta_dT_tip over window {ddT:+8.3f} +/- {ddT_err:.3f} K "
          f"(the whole window can only span {args.gradient*args.velocity*t_shot:.2f} K)")

    resid = v_mean / args.velocity - 1.0
    # Project the velocity deficit over the interface response time: that is
    # the undercooling the tip would still bank before settling.
    correction = args.gradient * (args.velocity - v_mean) * args.tau
    print(f"  velocity residual        {resid*100:+.2f}%  "
          f"(stationary if |.| <= {args.tol_v*100:.1f}%)")

    if abs(resid) <= args.tol_v:
        print(f"\n  VERDICT: STATIONARY. v_tip = V to within {args.tol_v*100:.1f}%.")
        print(f"  dT_tip = {args.dt_init:.3f} K "
              f"(+/- {abs(args.gradient*args.velocity*args.tau*args.tol_v):.1f} K "
              f"from the tolerance)  <-- deliverable (3).")
        print("  Confirm with the other seed, and with a longer shot, before "
              "believing it.")
    elif resid < 0:
        print("\n  VERDICT: v_tip < V -> the tip is banking undercooling -> "
              "dT_init is TOO WARM.")
        print(f"  Increase dT_init by roughly {correction:.1f} K "
              f"-> try {args.dt_init + correction:.1f} K.")
    else:
        print("\n  VERDICT: v_tip > V -> the tip is unwinding undercooling -> "
              "dT_init is TOO COLD.")
        print(f"  Decrease dT_init by roughly {abs(correction):.1f} K "
              f"-> try {args.dt_init + correction:.1f} K.")
        print("  (This is the safe side to be on: the front outruns the "
              "isotherm and warms.)")

    if abs(resid) > args.tol_v and not 0.0 < args.dt_init + correction < 31.468:
        print(f"\n  NOTE: the suggested {args.dt_init + correction:.1f} K is "
              "outside (0, dT0 = 31.468 K); solutal tip undercooling cannot "
              "exceed the freezing range.\n  Clamp the next shot inside the "
              "bracket -- the run script refuses dT_init >= dT0 (exit 29).")

    if max(dc) > 0.02:
        print("\n  WARNING: top-band composition exceeded +2% of c0. The solute "
              "layer is no longer contained and the alloy is enriching; the "
              "verdict above is not trustworthy. Raise ny.")


if __name__ == "__main__":
    main()
