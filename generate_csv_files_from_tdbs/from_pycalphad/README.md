# Al-Si TDB to MicroSim

This directory provides an open-source, reproducible route from a user-supplied
Al-Si TDB to the positional CSV files read by MicroSim
`Grand_Potential_MPI` with `Function_F = 4`.

The route is:

```text
user-supplied TDB
  -> in-memory Al-Si phase projection
  -> pycalphad Gibbs models
  -> common-tangent equilibrium compositions
  -> constrained phase Hessians
  -> MicroSim CSVs
```

It does **not** use the `from_Pandat` scripts. Those scripts only reformat text
already calculated and exported by proprietary Pandat; they do not read a TDB
or calculate equilibria/Hessians. The Thermo-Calc route similarly requires the
proprietary `TC-Python` API. This implementation calculates and writes the CSVs
directly with open-source pycalphad.

## Environment

Tested with Python 3.11 and the pinned packages in `requirements.txt`.

```bash
python3.11 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Run

The database path must point to a TDB that you obtained and are permitted to
use. The source TDB is never copied to the output.

```bash
MPLCONFIGDIR=/tmp/mpl-alsi \
  .venv/bin/python alsi_to_microsim.py \
  --tdb /private/path/to/AlSi.tdb \
  --output /private/path/to/AlSi_850_930K \
  --t-min 850 --t-max 930 --step 1 \
  --kks-temperature 860
```

The two-phase table begins at the calculated eutectic temperature if the
requested lower bound is below it. It adds points at `Te`, `Te+0.1`, `Te+0.2`,
and `Te+0.5 K`, followed by the 1 K grid.

## Output and conventions

`tdbs_encrypted/` contains the three files consumed positionally by
`Grand_Potential_MPI/functions/functionF_04.h`:

- `Composition_FCC_A1.csv`: `T`, equilibrium `x_Si` in FCC_A1, equilibrium
  `x_Si` in LIQUID.
- `HSN_FCC_A1.csv`: `T`, `d²Gm_FCC_A1/dx_Si²`.
- `HSN_LIQUID.csv`: `T`, `d²Gm_LIQUID/dx_Si²`.

Al is the dependent solvent. Compositions are **mole fractions**, not wt.%.
The binary Hessian is

```text
d(mu_Si - mu_Al)/d(x_Si) = d²Gm/d(x_Si)²
```

and is written in `J/mol` per mole-fraction squared. For a binary system the
upper-triangular Hessian contains one value. The diagnostic and target-alloy
CSVs report both mole and weight bases explicitly.

The script also writes:

- `metadata.json`: source SHA-256, software versions, units, validation,
  phase/order requirements, and limitations.
- `AlSi_thermodynamic_diagnostics.csv`: both bases, Hessians, quadratic KKS
  coefficients, and common-tangent residuals.
- `AlSi_target_alloys.csv`: conversions, liquidus temperatures, partition
  coefficients, and lever-rule eutectic fractions for 3.5–7.5 wt.% Si.
- `KKS_CuFFT_860K_parameters.txt`: an isothermal parameter fragment.

## Grand-potential MicroSim mapping

Use:

```text
Function_F = 4;
num_thermo_phases = 2;
tdb_phases = {FCC_A1,LIQUID};
phase_map = {FCC_A1,LIQUID};
```

The liquid must be the final simulation phase because `functionF_04.h`
addresses `NUMPHASES-1` as liquid. Place the generated `tdbs_encrypted`
directory in the solver run directory. The current reader is positional, so do
not reorder columns.

The CSV free energy is a local quadratic approximation about each equilibrium
phase composition. It is not the full CALPHAD free-energy surface.

## KKS_CuFFT mapping and limitation

At a selected fixed temperature:

```text
ceq(FCC_A1) = equilibrium x_Si in FCC_A1
ceq(LIQUID) = equilibrium x_Si in LIQUID
f0(phase) = Hessian(phase) / 2
```

`f0` is in J/mol. KKS_CuFFT divides it by `Vm`, which therefore must be supplied
in `m³/mol` from a separately validated source.

The current KKS_CuFFT solver has no TDB/CSV reader and no
temperature-dependent directional-solidification free energy. The generated
fragment is valid for an **isothermal verification case only**. A production
thermal-gradient run requires either a temperature-dependent reader/model to
be added to KKS_CuFFT or use of the grand-potential solver that already reads
these tables.

## NIMS CPDDB access

CPDDB is part of the free NIMS Materials Database service, **MatNavi**. A DICE
account alone is not sufficient:

1. In the DICE User Portal, open `Applications` and then
   `Identity Proofing Status`. The Institution field must read `Confirmed`.
2. Under `Application Usage`, apply for
   `NIMS Materials Database (MatNavi)`, not AtomWork-Adv.
3. If the institution is not confirmed, the organization's domain
   administrator (or an authorized person with the administrator copied) must
   apply to register the institutional email domain.
4. After MatNavi approval, log in to MatNavi, open CPDDB, search `Al-Si`, and
   manually download one assessment TDB.
5. Record the assessment authors/year, database citation, access date, source
   URL, terms, and SHA-256. Do not scrape or redistribute the TDB.

Candidate assessments listed by CPDDB include COST 507 (1998), Mey and Hack
(1986), and Murray and McAlister (1984). Select one actual downloaded file only
after checking its supplied provenance and terms.

## Validation gates

The program stops unless:

- the calculated eutectic is within 2 K of 850 K;
- eutectic liquid and FCC_A1 compositions are each within 0.2 wt.% Si of
  12.6 and 1.6 wt.% Si, respectively;
- all phase Hessians are positive;
- common-tangent residuals are below `1e-3 J/mol`;
- natural-cubic interpolation errors at 0.5 K holdout points are below 1%.

These are acceptance checks for the initial Al-Si database, not a substitute
for citing the selected CALPHAD assessment.

## Mobility and diffusivity

A thermodynamic TDB contains Gibbs-energy models. It does not automatically
provide atomic mobilities or phase diffusivities. Those require a separate
mobility database (often licensed independently) or a citable experimental
diffusion model. This script intentionally does not invent or infer them.
