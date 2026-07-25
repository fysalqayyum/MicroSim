#!/bin/bash
#SBATCH --job-name=microsim-opencl
#SBATCH --output=microsim-opencl-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:1

set -euo pipefail

: "${MICROSIM_ROOT:?Set MICROSIM_ROOT to the MicroSim checkout}"

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
solver_dir="${MICROSIM_ROOT}/KKS_OpenCl"
solver="${solver_dir}/microsim_kks_opencl"
run_id="${SLURM_JOB_ID:-manual}"
run_dir="${submit_dir}/run_${run_id}"
output_prefix="${MICROSIM_OUTPUT_PREFIX:-microsim}"

if [[ ! -x "${solver}" ]]; then
  echo "Solver is not executable: ${solver}" >&2
  exit 2
fi
for required in Input.in Filling.in tdbs_encrypted; do
  if [[ ! -e "${submit_dir}/${required}" ]]; then
    echo "Missing case item: ${submit_dir}/${required}" >&2
    exit 2
  fi
done
if [[ -e "${run_dir}" ]]; then
  echo "Refusing to overwrite existing run: ${run_dir}" >&2
  exit 3
fi

mkdir -p "${run_dir}/DATA/Processor_0"
cp "${submit_dir}/Input.in" "${run_dir}/Input.in"
cp "${submit_dir}/Filling.in" "${run_dir}/Filling.in"
cp -R "${submit_dir}/tdbs_encrypted" "${run_dir}/tdbs_encrypted"
cp -R "${solver_dir}/solverloop" "${run_dir}/solverloop"

{
  echo "purpose=single_gpu_opencl_run"
  echo "tasks=${SLURM_NTASKS:-1}"
  sha256sum "${solver}"
  sha256sum "${run_dir}/Input.in"
  sha256sum "${run_dir}/Filling.in"
  sha256sum "${run_dir}/solverloop/CL_Kim_Kernel.cl"
} > "${run_dir}/manifest.txt"

cd "${run_dir}"
"${solver}" Input.in Filling.in "${output_prefix}"

test -n "$(find DATA -name '*.vtk' -print -quit)"
