# KKS OpenCL HPC guide

This guide covers the maintained `KKS_OpenCl` path for a single OpenCL GPU on
an HPC scheduler. Commands use placeholders so the repository does not encode
one institution's filesystem, account, or module names.

## Supported execution topology

Use one MPI rank and one GPU.

The present multi-rank OpenCL implementation has not passed field equivalence
at internal and periodic x boundaries. More ranks or GPUs may produce a result
that exits normally but differs from the one-rank reference. Do not use that
path for scientific production until halo exchange is repaired and the
equivalence gates are repeated.

## Dependencies

The maintained path requires:

- a C compiler wrapper from a compatible MPI installation;
- OpenCL headers and loader;
- a vendor OpenCL runtime for the selected GPU;
- GSL and GSL CBLAS;
- HDF5 compiler support only if building the reconstruction/HDF5 utilities;
- Python and pycalphad only when generating thermodynamic tables.

Record compiler, MPI, OpenCL runtime, GPU, and library versions with every
qualified executable.

## Build

From `KKS_OpenCl`:

```bash
make clean
make microsim_kks_opencl \
  CC=/opt/mpi/bin/mpicc \
  GSL_ROOT=/opt/gsl \
  CUDA_ROOT=/opt/cuda
```

`CUDA_ROOT` is the location providing OpenCL headers/libraries on NVIDIA
systems; other vendors may require site-specific include and library
overrides:

```bash
make microsim_kks_opencl \
  CC=/opt/mpi/bin/mpicc \
  CFLAGS="-I. -I/opt/gsl/include -I/opt/opencl/include" \
  LDIR="-L/opt/gsl/lib -L/opt/opencl/lib"
```

The OpenCL compiler optimizer is enabled by default. For a diagnostic
comparison only, rebuild with:

```bash
make clean
make microsim_kks_opencl \
  CC=/opt/mpi/bin/mpicc \
  CFLAGS="-I. -I/opt/gsl/include -I/opt/opencl/include -DMICROSIM_DISABLE_CL_OPTIMIZER" \
  LDIR="-L/opt/gsl/lib -L/opt/opencl/lib"
```

Treat the optimized and no-optimization binaries as distinct artifacts and
record a SHA-256 checksum for each.

## Case directory

Run from an isolated case directory containing:

```text
case/
├── Input.in
├── Filling.in
├── solverloop/
├── tdbs_encrypted/
│   ├── Composition_FCC_A1.csv
│   ├── HSN_FCC_A1.csv
│   └── HSN_LIQUID.csv
└── DATA/
    └── Processor_0/
```

The executable loads OpenCL sources and generated headers through relative
paths. Copy the solver's `solverloop/` directory into the case or run from a
controlled directory where those relative files resolve correctly.

Create a fresh directory for every scheduler job. Refuse overwrites. Preserve:

- exact input and filling files;
- source and executable SHA-256 checksums;
- build command and dependency versions;
- scheduler job state and exit code;
- selected GPU and live utilization;
- stdout/stderr;
- complete field diagnostics.

## Slurm

Use the template in
`KKS_OpenCl/examples/single_gpu_slurm/run_single_gpu.slurm`.
Site-specific partition names, accounts, modules, and paths are supplied
through the environment or edited in a local untracked copy.

The essential resource request is one task and one GPU:

```text
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
```

High host CPU usage does not show that the kernels ran on the CPU. Confirm the
reported OpenCL device and sample live GPU utilization.

## Temperature-dependent Function F4

The host reads equilibrium-composition and Hessian CSVs as GSL splines. It
samples a reference-anchored uniform temperature table, and the OpenCL kernels
use four-point cubic lookup for equilibrium, phase free energy, Hessian,
relaxation, diffusion, chemical potential, and thermal-composition terms.

For the established two-phase binary path:

- keep `LIQUID` as the final phase;
- use mole fractions in the composition tables;
- provide Hessians in `J/mol` per mole-fraction squared;
- ensure every required CSV has a common temperature interval;
- keep the simulated temperature range inside the available table range.

The working F4 phase relaxation is calculated inside the kernel from the
thermodynamic Hessian, inverse liquid diffusivity, equilibrium composition
gap, gradient coefficient, and well height. Historical input entries named
`tau` and `Tau` are not active controls in this path.

## Restart

A restart must behave like the continuous run at the same absolute physical
step.

The maintained implementation:

1. reads phase, chemical potential, composition, and temperature blocks into
   the correct arrays;
2. consumes the newline separator after each binary block;
3. fills ghost cells using the no-flux boundary state;
4. mirrors the restored new state into the old-state buffer;
5. zeros displacement storage not present in the restart file;
6. evaluates smoothing, anti-trapping, and anisotropy gates using
   `t + STARTTIME`.

First run a short finiteness smoke test. Then compare a restarted continuation
against a continuous reference at the same final timestep. Do not validate a
restart only by checking that the executable returned zero.

## Moving window

The moving-window trigger follows the leading phase interface. After a shift,
new rows must be:

- liquid;
- filled with the far-field composition;
- assigned temperatures using absolute time and cumulative displacement.

A forced shift validates data plumbing only. It does not validate the physical
trigger or morphology.

## Output

`DATA/Processor_*/...vtk` files are rank-local MicroSim binary fragments, not
ordinary standalone legacy VTK files. Use the repository reconstruction tools
or a parser that understands:

- the ASCII dimensions/header;
- big-endian float64 field blocks;
- block separators;
- the exact phase, chemical-potential, composition, and temperature order.

Audit every saved field for non-finite values before visualization.

## Promotion sequence

Qualify changes in increasing risk:

1. clean build and host/device layout;
2. thermodynamic-table audit;
3. temperature sensitivity;
4. short finite-field smoke;
5. timestep sweep;
6. restart equivalence;
7. forced moving-window plumbing;
8. production-size short pilot;
9. morphology and convergence gates;
10. full production.

See [VALIDATION.md](VALIDATION.md) for current evidence and limitations.
