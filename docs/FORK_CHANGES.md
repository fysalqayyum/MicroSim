# Maintained fork changes

This fork tracks
[ICME-India/MicroSim](https://github.com/ICME-India/MicroSim) and keeps its
history and GPL-3.0 attribution. Fork changes are organized as reviewable
commits on top of upstream rather than copied into a disconnected code dump.

## Upstream baseline

The first maintained release is based on upstream commit
`365aa4137e922db883bc3a91643964e2131fe0c2` from 25 July 2026.

The ten upstream commits incorporated while preparing this release were:

| Commit | Upstream change | Relevance to the maintained path |
|---|---|---|
| `c5bfd1a` | Added a new AMReX code revision | Reverted by the next commit |
| `4656634` | Reverted `c5bfd1a` | No net solver change retained from the add/revert pair |
| `62aaef2` | Reorganized and repaired duplicated AMReX content | Retained as upstream; not qualified here |
| `0a69004` | Renamed and reorganized OpenFOAM/docs content | Retained; paths in older instructions may have changed |
| `235ebc4` | Added Grand Potential SYCL as a submodule | Transitional commit |
| `214bf3a` | Removed the temporary SYCL submodule entry | Transitional commit |
| `9af1a0a` | Imported the Grand Potential SYCL beta source | Retained as upstream beta; not qualified here |
| `ce9a428` | Removed legacy examples | Retained |
| `a91ca87` | Added the coupled Bridgman/grain OpenFOAM module | Retained as upstream experimental code |
| `365aa41` | Added the OpenFOAM 8 electrochemistry module | Retained as upstream experimental code |

None of these commits modified the files carrying this fork's qualified
`KKS_OpenCl`, `KKS_CuFFT`, Grand Potential portability, or Al-Si extraction
patches. The local commits were therefore rebased without source conflicts.

## Fork-maintained areas

### `KKS_OpenCl`

- Corrected local temperature-array sizing and growth-direction indexing.
- Added a temperature-dependent `Function_F=4` lookup table sampled from the
  host GSL splines and evaluated on the device with four-point cubic
  interpolation.
- Corrected composition indexing during host thermodynamic-table generation.
- Restored the anisotropic phase-gradient calculation used by the active F4
  kernel.
- Made the anti-trapping term obey the parsed `atr` switch.
- Added moving-window detection, shifting, liquid refill, composition refill,
  and absolute-time temperature refill.
- Repaired long-number parsing in `Tempgrady`, `ceq`, `cfill`, `slopes`, and
  `c_guess`.
- Repaired rank-local binary composition output and restart input.
- Initialized restart ghost cells and displacement storage.
- Changed restart physics gates to use `t + STARTTIME`.
- Enabled OpenCL optimization by default after a controlled field-equivalence
  gate; `MICROSIM_DISABLE_CL_OPTIMIZER` remains available for diagnosis.
- Removed the unused duplicate `CL_Kim_Kernelx.cl` runtime variant. Its history
  remains available in Git.

### Build portability

- Added macOS endian compatibility to selected file writers.
- Made the Cahn-Hilliard FFT dependency roots configurable.
- Added current CUDA CUB compatibility and configurable `NVCC`/`ARCH` values to
  `KKS_CuFFT`.
- Made Grand Potential MPI endian and input handling more portable.
- Made the primary `KKS_OpenCl` dependency roots configurable.

### Thermodynamics

`generate_csv_files_from_tdbs/from_pycalphad` converts a user-supplied Al-Si
TDB into the equilibrium-composition and Hessian CSVs consumed by MicroSim.
The input database is neither copied nor redistributed.

## Deliberate non-claims

- The upstream SYCL, AMReX, Bridgman, electrochemistry, and CUDA-MPI modules
  have not been qualified by this fork.
- The current multi-rank/multi-GPU OpenCL path is not accepted for scientific
  runs.
- Passing numerical and restart gates does not establish experimental
  calibration or dendritic morphology.
- Commercial or access-controlled thermodynamic databases are not part of the
  repository.
