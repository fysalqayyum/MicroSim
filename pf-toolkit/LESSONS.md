# Lessons — phase-field solidification with MicroSim KKS_OpenCl

Hard-won, each one paid for with GPU hours or a wrong conclusion. Read this
before designing a campaign. Nothing here is speculative; every item was
observed.

---

## 1. A run can finish "cleanly" and be meaningless

A directional-solidification run completed with `max|sum(phi) - 1| = 2.1e-13`,
zero non-finite values, and a moving-window mass audit exact to `1.7e-13` over
679 shift events. Its output was physically meaningless for the last 15% of its
duration: the front had banked undercooling well past the freezing range and
then released it at **15 × the pulling velocity** into a grid-scale structure.

**`sum(phi)` and `nonfinite` are not validity checks.** `sum(phi)` actually
*improved* across the failure (2.06e-13 → 6.27e-14).

> **A metric that cannot represent the failure mode is not a guard.**

This is why `functions/farfield_guard.h` exists. It watches quantities that can
represent the failure: composition drift, tip margin, front velocity.

## 2. Undercooled bulk liquid is INERT — do not size domains for thermal headroom

For `npha = 2` the interpolation is `h(φ) = 6φ⁵ − 15φ⁴ + 10φ³`, so
`h'(0) = h''(0) = 0`: the driving-force term vanishes to **second order** at
φ = 0. And `addNoise` returns immediately when `envelope = φ0·φ1 <= 0`, so bulk
liquid receives **exactly zero** stochastic forcing.

Uniform liquid at φ = 0 is therefore an exact stationary state at any
undercooling and **cannot nucleate**. Confirmed in practice: liquid sat **14.3 K
below its liquidus for 1.5e6 steps with zero nucleation**, peak excursion 26 K.

**Consequence:** sizing the liquid so the top row reaches the liquidus isotherm
wastes enormous compute. Size it for **solute containment** instead. The binding
constraint with a no-flux top is that the solute boundary layer must fit.

**Corollary — do not use `l_D = D/V` as the containment estimate.** Measured
22–40 `l_D` on fronts that were not yet at steady state, because `l_D` assumes
the front runs at `V`. Measure it.

## 3. Far-field undercooling is an algebraic identity — it tells you nothing

In a linear thermal gradient:

```
undercooling_farfield ≡ dT_tip − cells_above_tip·G·dx + (T_liq_band − T_liq_c0)
```

Verified over 45 steady-state samples: residual **0.0000 K, sd 0.0000 K**, exact
to machine precision. Any threshold on it is just a box-sizing rule restated in
kelvin, and in a short box it can never pass. An abort built on it would have
killed a valid run at 35% completion. It was removed.

**Banking of undercooling is a TIP quantity**, `d(dT_tip)/dt`. Measure that.

## 4. `shift.dat` is the only authoritative front position

Once the moving window fires, `tip_row` is pinned at `Shiftj` by construction.
Any velocity derived from it reads exactly zero, forever.

A guard that computed the "absolute" front position as `tip_row + shift_position`
was blinded for exactly this reason: **`shift_position` is set once when a
restart file is read and never updated during a run.** The live accumulator is
`shift_OFFSET`. The frozen value sat at 1112 for 1.4e6 steps while the true
front advanced 4687 cells — and the runaway guard could not fire.

```
velocity = Δcells × dx / (Δsteps × dt)      # from shift.dat
```

> **An exactly-constant diagnostic column is a bug until proven otherwise.**

## 5. Match your observation window to the transient timescale

Five runs of 190,000 steps each "showed" stalled fronts and prompted a whole
false hypothesis about contaminated initial conditions. The real transient was
**~1.9e6 steps** — the runs were stopped 10× too early to show anything. All
were still accelerating when they ended.

The incubation clock is the **morphological instability**, not solute diffusion:

| process | timescale | steps at dt = 2e-9 |
|---|---|---|
| solute layer `D/V²` | 3.52e-5 s | 17,579 |
| instability `1/σ` | 2.64e-4 s | 131,752 |
| observed lag to breakdown | 2.72e-3 s | ~1,360,000 |

The diffusion layer is built in 1.3% of the incubation. Pre-imposing an analytic
solute boundary layer addresses the wrong lever.

## 6. Fit a trend; never difference endpoints

When sample scatter is comparable to the drift, endpoint differencing gives the
**wrong sign**. Differencing `dT_tip` over 1e6 steps suggested `v/V = 0.97`;
windowed linear fits on the same data gave `v/V > 1`. The sd (0.37 K) exceeded
the trend.

## 7. `20 D/V²` is falsified as a transient timescale

Fit `v_tip(t) = V(1 − B e^{−t/τ})` with the asymptote **PINNED to V**. Never fit
the asymptote free — it produces plausible nonsense.

Banked undercooling for a cold start is `dT_bank = G·V·τ·B`. Sizing a domain to
absorb a cold-start transient is a losing game; seed a developed front instead.

## 8. Steady state needs BOTH halves of the criterion

`v_tip = V` **and** `dT_tip` stationary. Reading `v_tip = V` on a *rising*
`dT_tip` is an overshoot, not steady state. That misreading cost two full legs.

Since `d(dT_tip)/dt = G(V − v_tip)`, the two are the same condition — so a run
need only show which side of it you are on. Use a **relative** criterion
(`|1 − v/V| ≤ 2%`), not an absolute temperature tolerance: an absolute 1 K
tolerance passed a front running at **0.8 V**, because the whole measurable range
over the observation window was only 3.5 K.

Transient overshoot to ~1.6 × V is legitimate physics. A runaway reaches 15 × V.
Set velocity guards to distinguish them (3.0 works; 2.0 aborts on real physics).

## 9. Shell and build traps that cost real time

- **Never end a pipeline in `grep -q` under `set -o pipefail`.** `grep -q` exits
  at the first match, the upstream command dies of SIGPIPE, and `pipefail`
  returns **141**. A binary check written this way refused *every* binary,
  including ones that provably contained the symbol. Use `grep -c`, which drains
  stdin. Same family as "`| tail` swallows the exit code."
- **`clBuildProgram` failure does NOT abort — it only prints.** A broken `.cl`
  fails quietly and the run produces garbage. Validate kernels offline first:
  ```
  clang -x cl -cl-std=CL1.2 -fsyntax-only -Xclang -finclude-default-header \
        -target spir64 -I. solverloop/CL_Kim_Kernel.cl
  ```
  `-target spir64` is required (the default target lacks `cl_khr_fp64`).
  **The host `make` does not compile the `.cl`** — the OpenCL driver does, at run
  time.
- **A guard must fail LOUD, never fail closed.** Distinguish "condition absent"
  from "the detector could not run" with separate exit codes.
- **Dry-run every guard against pass AND fail inputs.** `bash -n` does not
  execute embedded interpreters, and `sed` is silent on a miss.

## 10. Solver behaviours worth knowing

- **`NTIMESTEPS` is relative to the leg.** Every leg length must be an exact
  multiple of `SAVET`, or the leg writes no final frame and cannot be resumed.
- **Restart files are not VTK** despite the extension. Layout:
  ```
  3 ASCII header lines (nx, ny, nz) -- READ BUT NEVER VALIDATED
  phi[a] : nx*ny*nz float64 BIG-ENDIAN per phase, each + b'\n'
  mu[k]  : same, NUMCOMPONENTS-1 blocks
  c[k]   : same, NUMCOMPONENTS-1 blocks
  T      : same, only when !ISOTHERMAL, no trailing separator
  ```
  Index order is x outer, z middle, y inner: `index = y + ny*(z + nz*x)`.
  Byte count `5*8*nx*ny*nz + len(header) + 4`. The reader accepts a frame it did
  not write — but the `b'\n'` separators are a real integrity check, so a wrong
  mesh fails loudly rather than silently.
- **Temperature is read from a restart and then immediately overwritten** by
  `apply_temperature_gradientY`. The thermal field is 100% controlled by the
  input file.
- **Clock reset** = `RESTART=1` + `STARTTIME=0` + `shift.dat` containing `0 0`.
  Setting `STARTTIME` to the real step count instead RESTORES all accumulated
  cooling — which is why you cannot simply "continue from frame N" after a
  thermal reset.
- **`fill_composition_cube` initialises liquid with `mu = Mu(c, Teq, liquid)` at
  `Teq`, not the local temperature.** Every liquid cell therefore starts
  thermodynamically inconsistent with its position in a gradient. Suspected
  cause of an unexplained far-field composition drift; unresolved.
- **There is NO Dirichlet BC** in this solver. The kernel set is `noflux` and
  `periodic` only, for φ, c, T and elasticity alike. With a no-flux top, a
  truncated boundary layer means rejected solute accumulates and the alloy
  silently enriches — and **the mass audit cannot detect this**, because its
  `expected_after` is *defined* to include refill at c0.
- **`.in` keys that no parser reads are silently ignored.** Two such dead keys
  were cited as real settings for weeks. Grep the source before trusting a key.

## 11. Resolution and domain-size traps

- **Throughput depends on `ny`, not `nx`.** Measured: `ny` 1888 → 1.64e8,
  2240 → 1.73e8, 2832 → 1.79e8 cells·steps/s; but 1984 → 5.7e7 and
  5600 → < 2.6e7. The dependence is **not monotonic** — benchmark any new `ny`
  before committing to it.
- **Periodic x quantises the primary spacing** to `L/n`. A 34 µm domain has a
  ±10% quantisation floor near 3 µm. Halving the domain roughly doubles that
  floor. Never report a spacing from a domain holding only 2–3 primary spacings.
- **Selected wavelength ≠ primary spacing.** Report `lambda1` only once the
  protrusion count has stabilised AND the FFT dominant mode agrees with
  `L/n_prot`. In one case they disagreed by 2.3×, and the real-space number was
  the honest one.
- **Raising anisotropy to force dendrites is a fitting trap.** At the physical
  value for Al (`dab = 0.01`) expect cells or doublon tips. Side branching should
  come from noise, not from an inflated `dab`. Observed: doublon at 0.01–0.02,
  side-branch onset 0.03, full dendrite 0.05.

## 12. Validate diagnostics against a known answer

`reembed_frame.py --self-test` round-trips a real frame **byte-identically**.
Do this after any edit to a file-format tool.

A convergence checker once passed silently on a `0.0` initialiser with zero
valid frames. **A metric with no valid data must FAIL, not pass; an
exactly-zero result is a code signature, not physics.**
