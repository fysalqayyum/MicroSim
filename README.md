# MicroSim

MicroSim is an open-source collection of phase-field solvers for CPU, MPI,
CUDA, OpenCL, SYCL, AMReX, and OpenFOAM workflows. The upstream project is a
collaboration among IISc Bangalore, IIT Hyderabad, IIT Bombay, IIT Madras,
Savitribai Phule Pune University, and C-DAC Pune.

This repository is a maintained fork of
[ICME-India/MicroSim](https://github.com/ICME-India/MicroSim). It follows the
upstream `main` branch while maintaining portability, GPU/HPC, restart, output,
and reproducibility repairs developed during real cluster use.

## What this fork adds

The maintained changes are deliberately narrow:

- portable macOS/Linux builds for selected CPU and MPI solvers;
- configurable CUDA architecture and current-toolkit CUB support for
  `KKS_CuFFT`;
- a repaired `KKS_OpenCl` `Function_F=4` path with temperature-dependent
  thermodynamics;
- correct moving-window refill and temperature evolution along the growth
  direction;
- restart-safe absolute timestep handling, ghost-cell initialization, and
  binary field restoration;
- correct rank-local binary composition output;
- safer input parsing for long floating-point values;
- default OpenCL compiler optimization with a diagnostic rollback switch;
- an open-source, user-supplied-TDB to MicroSim CSV workflow based on
  pycalphad;
- a reusable one-GPU Slurm example and explicit validation limitations.

See [Fork changes](docs/FORK_CHANGES.md), the
[OpenCL HPC guide](docs/KKS_OPENCL_HPC.md), and
[validation status](docs/VALIDATION.md).

## Module status

The repository contains modules at different maturity levels. A successful
build of one module does not validate the others.

| Module | Backend | Status in this fork |
|---|---|---|
| `Cahn_Hilliard_FFT_2D` | CPU, FFTW/GSL | Portable build fixes; local smoke-tested |
| `Grand_Potential_Serial` | CPU | Portable binary-output helpers; local smoke-tested |
| `Grand_Potential_MPI` | MPI/HDF5 | macOS/Linux portability and parser fixes |
| `KKS_CuFFT` | CUDA | Current CUB compatibility and selectable GPU architecture |
| `KKS_OpenCl` | MPI/OpenCL | Maintained fork path; one-rank/one-GPU F4 workflow qualified |
| `KKS_FD_CUDA_MPI` | CUDA MPI | Upstream module; not ported or qualified by this fork |
| `Grand_Potential_AMReX` | AMReX | Upstream module; not qualified by this fork |
| `Grand_Potential_OpenFOAM` | OpenFOAM | Upstream module; not qualified by this fork |
| `Grand_Potential_SYCL` | SYCL | Upstream beta; not qualified by this fork |
| `Bridgman_Grain` | OpenFOAM | Upstream experimental module; not qualified by this fork |
| `Electrochemistry_Module_OpenFoam8` | OpenFOAM 8 | Upstream experimental module; not qualified by this fork |

“Qualified” here means that a controlled software gate passed. It does not
mean that every material parameter or simulated morphology is experimentally
calibrated.

## Quick start

Clone this fork:

```bash
git clone git@github.com:fysalqayyum/MicroSim.git
cd MicroSim
```

### KKS OpenCL on a single GPU

The maintained HPC path uses exactly one MPI rank mapped to one OpenCL GPU.
Multi-rank OpenCL runs are not recommended until the halo/boundary path is
repaired and requalified.

Build with site-specific paths supplied on the command line:

```bash
cd KKS_OpenCl
make microsim_kks_opencl \
  CC=/path/to/mpicc \
  GSL_ROOT=/path/to/gsl \
  CUDA_ROOT=/path/to/cuda
```

Run from a case directory that contains `Input.in`, `Filling.in`,
`tdbs_encrypted/`, and the generated `solverloop/` files:

```bash
/path/to/MicroSim/KKS_OpenCl/microsim_kks_opencl \
  Input.in Filling.in output_prefix
```

A site-neutral Slurm template is available under
[`KKS_OpenCl/examples/single_gpu_slurm`](KKS_OpenCl/examples/single_gpu_slurm).
No research case or thermodynamic database is bundled. Generate tables from a
TDB that you are permitted to use with
[`generate_csv_files_from_tdbs/from_pycalphad`](generate_csv_files_from_tdbs/from_pycalphad).

### KKS CuFFT

```bash
cd KKS_CuFFT
make NVCC=/path/to/nvcc ARCH=-arch=sm_80
./microsim_kks_cufft Input.in Filling.in output_prefix
```

Choose `ARCH` for the target GPU; do not copy an architecture value from
another system without checking the device compute capability.

### Grand Potential serial

```bash
cd Grand_Potential_Serial
make
./microsim_gp Input.in Filling.in output_prefix
```

### Cahn-Hilliard FFT

```bash
cd Cahn_Hilliard_FFT_2D
make
./microsim_ch_fft Input.in Filling.in output_prefix
```

The exact dependency locations vary by operating system and HPC site. See
[`Installation_instructions.md`](Installation_instructions.md) and the README
inside each module before building.

## Thermodynamic data

This fork does not redistribute private or commercial thermodynamic
databases. The pycalphad converter accepts a user-supplied Al-Si TDB, performs
the equilibrium and Hessian calculations in memory, and writes the positional
CSV tables expected by MicroSim.

Thermodynamic data do not provide atomic mobility or phase diffusivity.
Mobility and diffusivity require independent, citable sources.

## Restart and output rules

For `KKS_OpenCl`:

- physics gated on the timestep must use absolute time,
  `restart-relative step + STARTTIME`;
- a restart must initialize ghost cells and device-side displacement storage;
- the rank-local binary reader restores phase, chemical potential,
  composition, and temperature into their matching arrays;
- rank fragments are MicroSim binary fragments and should be reconstructed or
  converted before opening them in ParaView;
- restart acceptance requires a field comparison against a continuous run at
  the same physical timestep.

Detailed procedures and tolerances are in [docs/VALIDATION.md](docs/VALIDATION.md).

## Keeping the fork current

The recommended remote layout is:

```text
origin    git@github.com:fysalqayyum/MicroSim.git
upstream  https://github.com/ICME-India/MicroSim.git
```

Update a clean branch with:

```bash
git fetch upstream
git rebase upstream/main
```

Resolve conflicts by preserving upstream changes first, then reapplying and
retesting the fork-specific patch. Never assume that a conflict-free rebase
means the solver remains scientifically equivalent.

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request. Keep
code, small inputs, scripts, and compact expected diagnostics in Git. Keep
binaries, raw VTK/HDF5 output, job logs, licensed databases, credentials, and
machine-specific configuration out of the repository.

Run the pre-publication check before pushing:

```bash
scripts/pre-publish-check.sh upstream/main
```

## Attribution and license

MicroSim is distributed under the
[GNU General Public License v3.0](LICENSE). The original module authors and
institutional contributors retain their upstream attribution. Fork-specific
changes are modifications of that GPL-licensed work; they do not replace the
original authorship record.

For upstream releases, workshops, and the full project history, visit
[ICME-India/MicroSim](https://github.com/ICME-India/MicroSim).
