# pf-toolkit

Analysis, thermodynamics and launch tooling for MicroSim `KKS_OpenCl`
phase-field solidification runs, with an emphasis on **catching invalid runs
early**.

These tools were built while running two-phase directional solidification of a
hypoeutectic Al–Si alloy with `Function_F = 4`. They are alloy-agnostic wherever
that was practical; where a default is Al–Si it is documented as such and
overridable.

**Read [`LESSONS.md`](LESSONS.md) first.** It is the most valuable file here. It
records failures that each cost GPU-hours or a wrong conclusion, including the
one that matters most: *a run can satisfy every conventional validity check and
still be physically meaningless*.

## Layout

```
pf-toolkit/
├── LESSONS.md              <- start here
├── postprocessing/         <- frame analysis, diagnostics, restart tooling
├── thermodynamics/         <- CALPHAD -> MicroSim CSV generation (COST507 Al-Si)
└── slurm/                  <- launch template with pre-flight guards
```

## postprocessing/

All scripts take `--help`. Units are SI unless stated; `dx` in metres, `dt` in
seconds, compositions in mole fraction, temperatures in kelvin.

### Front and interface analysis
| script | what it does |
|---|---|
| `interface_amplitude.py` | interface amplitude, RMS, peak-to-valley, Fourier modes, instability growth rate σ, front and tip velocity. **Requires `--shift-file` or `--no-shift`** — absolute heights are wrong without it; RMS/PtV/spectrum are shift-invariant either way. Watch `n_columns_no_crossing`: once the window shift exceeds the groove depth, groove roots scroll out and are interpolated across. |
| `analyze_pilot_morphology.py` | shared VTK reader (`(nx, ny, nz)`, axis 1 is y) plus morphology metrics. Other scripts import `read_fragment` from here. |
| `tip_height.py` | topmost solid row per frame |
| `analyze_tip_undercooling.py` | tip undercooling and time-to-steady, measured rather than assumed |
| `analyze_seed_tip_radius.py` | parabola fit to the tip. Note: a doublon's split tip defeats the fit — check validity counts before quoting a radius. |
| `shot_verdict.py` | steady-state verdict from a guard CSV. Uses `d(dT_tip)/dt = G(V − v_tip)` and a **relative** criterion; integrates net displacement rather than differencing per-sample velocity. |
| `solute_profile.py` | x-averaged composition versus height; use it to **measure** the solute containment depth |
| `probe_nucleation.py` | undercooling versus height, for asking whether liquid ahead of the front could nucleate |

### Validity and convergence
| script | what it does |
|---|---|
| `check_phase_sum.py` | `sum(phi)` audit. **Necessary but nowhere near sufficient** — see LESSONS.md §1. |
| `check_xsi.py` | quick composition range check |
| `compare_vtk_fields.py` | field-by-field comparison of two frames |
| `compare_w_convergence.py` | interface-width convergence at matched physical time |
| `qualify_restart_equivalence.py` | restart equivalence by relative L2, not bit-identity (bit-identity produces false negatives) |
| `compute_dimensionless_groups.py` | Péclet, `d0`, `l_D`, `Γ`, stability and resolution ratios from raw inputs |
| `test_addnoise_hash.py` | reference vectors for the `addNoise` PRNG mapping |

### Files, restarts and rendering
| script | what it does |
|---|---|
| `reembed_frame.py` | crop/pad a restart frame into a new domain so the solver accepts it as step 0. Seam-matched x-cropping. **`--self-test` round-trips a real frame byte-identically — run it after any edit.** |
| `merge_rank_fragments.py` | merge per-rank VTK fragments |
| `diagnose_rank_fragments.py` | diagnose fragment coverage and gaps |
| `inspect_latest_frame.py` | quick look at the newest frame |
| `generate_multimode_filling.py` | multimode perturbation `Filling.in` generator |
| `render_frame_series.py` | frame series, GIF and montage. **Pass the whole frame set every time** — the colour scale is global over whatever is passed, so appending frames breaks comparability. |
| `make_paraview_series.py` | ParaView-friendly series |
| `watch_leg_complete.sh` | wait for a final frame to be *completely* written (exact byte count, twice), then run the interface diagnostic |

## thermodynamics/

`cost507_alsi/` contains everything needed to regenerate the MicroSim
`Function_F = 4` thermodynamic tables for Al–Si from an openly redistributable
CALPHAD database:

- `COST507_upstream_pycalphad.tdb` — untouched upstream database
- `COST507_MicroSim_compatible.tdb` — nomenclature-compatible variant
- `prepare_cost507_for_microsim.py` — the generator
- `Composition_FCC_A1.csv`, `HSN_FCC_A1.csv`, `HSN_LIQUID.csv` — the tables the
  solver reads (`HSN` = Hessian, `d²G_m/dx²` in J/mol)
- `metadata.json` — full provenance, software versions, validation gates
- `AlSi_thermodynamic_diagnostics.csv`, `SHA256SUMS.txt`, `README.md`

Validation against the accepted Al–Si eutectic: **850.150 K** (target 850, tol 2),
**12.502 wt% Si liquid** (target 12.6, tol 0.2), **1.549 wt% Si α** (target 1.6,
tol 0.2). Natural-cubic interpolation error < 0.01% on compositions and < 0.8%
on Hessians.

`plot_alsi_phase_diagram.py` renders the Al-rich phase diagram from the same TDB.

> **Licensing.** COST507 is redistributed here under its own terms: the TDB
> header permits free use at the user's own risk, and the European Commission
> report permits reproduction with source acknowledgement. This is *freely
> reusable with attribution*, not an OSI licence. Attribution details and
> upstream checksums are in `cost507_alsi/README.md`. **No commercially licensed
> database is included in this repository**, and no table derived from one.

The TDB supplies Gibbs energies only — **no mobilities, diffusivities, molar
volumes or interface properties.** Those must come from a separately licensed
mobility database or a citable experimental model.

## slurm/

`directional_solidification.slurm.example` is a template with no site-specific
paths. Its value is the **pre-flight guard structure**, not the parameters:
refuses unreadable or outdated binaries, leg lengths that are not a multiple of
`SAVET` (which silently produce unresumable runs), and domains too short to hold
the solute boundary layer. Each failure has a distinct exit code.

`Input_resolved_alsi55.in.template` is a token-substituted input file;
`validate_inputs.py` asserts no placeholder tokens survive substitution.

## Requirements

Python 3.11+, NumPy, SciPy, pandas, matplotlib. `pycalphad` 0.11.2 for the
thermodynamics tools only. No package installation step — the scripts are
standalone and import `read_fragment` from `analyze_pilot_morphology.py` by path.

## Contributing back

Bug reports and generalisations welcome, particularly making the alloy-specific
helpers (`farfield_wt_si`, the linearised liquidus) properly multicomponent.
