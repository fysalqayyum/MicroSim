# Site-neutral single-GPU Slurm template

This directory contains scheduler plumbing only. It does not contain a
research input, thermodynamic table, cluster name, partition, account,
username, software path, or production parameter.

Prepare a separate case directory containing:

```text
Input.in
Filling.in
tdbs_encrypted/
```

Build `KKS_OpenCl/microsim_kks_opencl`, copy `run_single_gpu.slurm` into the
case directory, and submit it with site-specific options supplied outside the
tracked file:

```bash
export MICROSIM_ROOT=/path/to/MicroSim
sbatch --partition=your-partition run_single_gpu.slurm
```

If the site requires an allocation/account, provide it through the submission
command or a local untracked wrapper.

The script uses one MPI process and one GPU, creates an isolated run directory,
records source/input hashes, and refuses to overwrite an existing run.
