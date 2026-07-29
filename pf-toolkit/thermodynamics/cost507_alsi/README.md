# Open COST507 Al–Si thermodynamics for pycalphad and MicroSim

## Status

This directory contains an openly available, independently assessed Al–Si
thermodynamic cross-check derived from the COST 507 light-alloy database.
The 84-row extraction completed with pycalphad 0.11.2, passed the same
eutectic reference gates used for an independent commercial cross-check, and passed
the 1% natural-cubic interpolation gate for both compositions and Hessians.

This dataset is suitable for independent thermodynamic comparison and
controlled MicroSim `Function_F=4` verification. It does not contain mobility,
diffusivity, molar-volume, or interface-property data.

## Upstream sources and reuse terms

- OpenCalphad database page:
  <https://www.opencalphad.com/databases.php>
- Original COST507 archive:
  <https://www.opencalphad.com/databases/COST507.zip>
- European Commission COST 507 report:
  <https://www.opencalphad.com/databases/CGNA18499ENC_001.pdf>
- Pycalphad-maintained corrected TDB:
  <https://raw.githubusercontent.com/pycalphad/pycalphad/develop/pycalphad/tests/databases/COST507.tdb>

The TDB header states that the database “may be used freely” at the user's own
risk. The European Commission report states that reproduction is authorized
provided the source is acknowledged. These are explicit reuse permissions,
but they are not expressed as a modern SPDX or Open Data Commons license.
Therefore the cautious description is **freely reusable with attribution**,
not “OSI-licensed.”

The Al–Si assessment is attributed in the TDB to H. L. Lukas, updated
1994-12-14. COST 507 is the frozen Round II light-alloy database dated January
1999; the accompanying report was edited by I. Ansara, A. T. Dinsdale, and
M. H. Rand and published by the European Commission in 1998.

## File lineage and checksums

| Artifact | SHA-256 |
|---|---|
| Original OpenCalphad `COST507.zip` | `a5e5856a758c6881a692a0f554890d6a266d16bfb6a029d197f4b01df26141c0` |
| Original archived `COST507.tdb` | `80363e1de4aa1ee0750df28ef48e1524ef64b7988d6517649fff1c991e3ed641` |
| Local, untouched pycalphad-maintained `COST507_upstream_pycalphad.tdb` | `6565f4d67695fbbf779d68b9eeeb83e86e44a1a3dde9bafd612d99daffb9c34d` |
| Local `COST507_MicroSim_compatible.tdb` | `bfba3f8bae2b8f086792cb2ec758b70eff0ff5c2064bfe5394874f489d8a738d` |

The original OpenCalphad TDB does not parse directly in pycalphad 0.11.2
because it contains legacy syntax errors, beginning with a missing `E` in the
yttrium atomic-mass exponent. Pycalphad distributes a corrected copy in its
test database set. That corrected file parses successfully and is the upstream
source used here.

Use `COST507_upstream_pycalphad.tdb` directly for general pycalphad work. Use
`COST507_MicroSim_compatible.tdb` with the current MicroSim Al–Si extraction
script, whose phase-name and parser conventions require the three aliases
listed below. Both files are retained locally so future calculations do not
depend on network availability.

`prepare_cost507_for_microsim.py` makes only three compatibility edits:

1. expands the legacy `CONST` abbreviation to `CONSTITUENT`;
2. removes the legacy `:L` suffix from the LIQUID phase declaration; and
3. aliases `DIAMOND_A4` to the phase name `SI_DIAMOND_A4` expected by the
   existing MicroSim extractor.

No thermodynamic coefficient, function, phase model, constituent set, or
temperature interval is changed.

## Key thermodynamic results

| Quantity | COST507 result |
|---|---:|
| Eutectic temperature | 850.1497982085 K |
| FCC_A1 eutectic composition | 0.01488551278 mole fraction Si |
| FCC_A1 eutectic composition | 1.5485128054 wt.% Si |
| LIQUID eutectic composition | 0.12070131330 mole fraction Si |
| LIQUID eutectic composition | 12.5022425284 wt.% Si |
| Eutectic maximum tangent residual | 3.9314e-7 J/mol |
| Maximum two-phase temperature written | 930 K |
| Rows written | 84 |

At 860 K:

| Quantity | COST507 result |
|---|---:|
| `x_Si_FCC_A1` | 0.0129932392315 |
| `x_Si_LIQUID` | 0.1080164132301 |
| `H_FCC_A1` | 563178.606433 J/mol |
| `H_LIQUID` | 97947.408388 J/mol |

Natural-cubic midpoint maximum relative errors were
`1.2991e-5` for FCC_A1 composition, `3.3064e-5` for LIQUID composition,
`0.0079217` for the FCC_A1 Hessian, and `0.0077128` for the LIQUID Hessian.

## Comparison with an independent commercial assessment

As a sanity check the same extraction was run against a separately licensed
commercial Al–Mg–Si database. **No content from that database is redistributed
here** — neither the database itself nor any table derived from it. The five
scalar values below are reported solely as a validation comparison, in the same
sense that a paper quotes a literature value.

The two descriptions agree closely on eutectic temperature but are not
identical:

| Quantity | COST507 | Commercial assessment | Difference |
|---|---:|---:|---:|
| Eutectic T (K) | 850.149798 | 850.195902 | -0.046104 |
| FCC_A1 Si (wt.%) | 1.548513 | 1.729928 | -0.181415 |
| LIQUID Si (wt.%) | 12.502243 | 12.637072 | -0.134830 |
| 860 K FCC Hessian (J/mol) | 563178.606 | 506919.925 | +56258.682 |
| 860 K LIQUID Hessian (J/mol) | 97947.408 | 97244.929 | +702.480 |

The approximately 11.1% difference in the 860 K FCC_A1 Hessian is large enough
to affect the local quadratic free energy and should be included in a
thermodynamic sensitivity study.

## Reproduction

The upstream pycalphad-maintained TDB has already been archived locally. A
direct pycalphad parse can be checked from the project `DATA/` directory with:

```bash
../.venvs/pycalphad-0.11.2/bin/python -c \
  "from pycalphad import Database; db=Database('thermodynamics/open_cost507_alsi_850_930K/COST507_upstream_pycalphad.tdb'); print(len(db.phases), sorted({'LIQUID','FCC_A1','DIAMOND_A4'} & set(db.phases)))"
```

To reproduce the download and prepare the alias-compatible copy:

```bash
curl -L -o COST507_upstream_pycalphad.tdb \
  https://raw.githubusercontent.com/pycalphad/pycalphad/develop/pycalphad/tests/databases/COST507.tdb
python3 prepare_cost507_for_microsim.py \
  COST507_upstream_pycalphad.tdb COST507_MicroSim_compatible.tdb
```

From the project `DATA/` directory, regenerate the tables:

```bash
MPLCONFIGDIR=/tmp/mpl-cost507 python prepare_cost507_for_microsim.py \
  --tdb COST507_MicroSim_compatible.tdb \
  --output . \
  --t-min 850 --t-max 930 --step 1 --kks-temperature 860
```

Pycalphad emits legacy phase-type warnings for `X`, `R`, and `F` after the
extractor's in-memory projection. These warnings do not prevent model
construction or validation, but they should remain recorded rather than
silenced.
