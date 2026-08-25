#!/bin/bash
#
# Final normal_sigma_Id versus standardized multivariate-t production. The two
# scientific invocations run sequentially because nu depends on d:
#   d=2, nu=3; then d=5, nu=6.
# Run `mkdir -p logs` and the login-node preflight before sbatch.
#
#SBATCH --job-name=normal_sigma_Id_t
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/normal_sigma_Id_t_%j.out
#SBATCH --error=logs/normal_sigma_Id_t_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=07:15:00
#SBATCH --mem=6900M
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 0 ]]; then
  echo "Usage: sbatch scripts/run_c3_normal_sigma_Id_t_production.sbatch" >&2
  exit 2
fi

: "${SLURM_SUBMIT_DIR:?Submit this file with sbatch from the repository root}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"
if [[ "$SLURM_CPUS_PER_TASK" -ne 32 ]]; then
  echo "Production requires exactly 32 CPUs; Slurm assigned ${SLURM_CPUS_PER_TASK}." >&2
  exit 2
fi
cd "$SLURM_SUBMIT_DIR"

output_d2="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n50_100_200_400_M1000_B5000_Nderiv10000"
output_d5="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n50_100_200_400_M1000_B5000_Nderiv10000"

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

# The runner owns the 32-process outer queue; each bootstrap uses one core.
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "===== C3 normal_sigma_Id t production environment ====="
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "pwd=$(pwd)"
echo "commit=$(git rev-parse HEAD)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-not_set}"
echo "SLURM_JOB_NAME=${SLURM_JOB_NAME:-not_set}"
echo "SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "SLURM_NTASKS=${SLURM_NTASKS:-not_set}"
echo "SLURM_JOB_NUM_NODES=${SLURM_JOB_NUM_NODES:-not_set}"
echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-not_set}"
echo "SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU:-not_set}"
echo "R_LIBS_USER=${R_LIBS_USER}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
echo "MKL_NUM_THREADS=${MKL_NUM_THREADS}"
echo "output_d2=${output_d2}"
echo "output_d5=${output_d5}"
echo "===== Loaded modules ====="
module list 2>&1 || true
echo "===== R version ====="
R --version
echo "===== Slurm allocation ====="
scontrol show job "$SLURM_JOB_ID" || true

# Validate both destinations before starting either scientific invocation.
echo "===== In-job deployment and output preflight ====="
bash scripts/preflight_c3_normal_sigma_Id_t_production.sh

run_invocation() {
  local d="$1"
  local nu="$2"
  local output_dir="$3"
  local started_epoch
  local finished_epoch
  local runner_wall_seconds
  local runner_status
  local validation_status

  echo "===== normal_sigma_Id runner start: d=${d}, nu=${nu} ====="
  echo "date=$(date --iso-8601=seconds)"
  echo "output_dir=${output_dir}"
  started_epoch=$(date +%s)
  set +e
  Rscript --vanilla scripts/run_normal_sigma_Id_t_pilot.R \
    "--output_dir=${output_dir}" \
    --M=1000 --B=5000 "--dimensions=${d}" \
    --n_values=50,100,200,400 --beta_values=0,0.25,0.5,1 \
    "--nu=${nu}" --derivative_mc_size=10000 --cvm_block_size=50 \
    "--cores=${SLURM_CPUS_PER_TASK}" --checkpoint_results=256 \
    --seed=20260728 --show_progress=true
  runner_status=$?
  set -e
  finished_epoch=$(date +%s)
  runner_wall_seconds=$((finished_epoch - started_epoch))
  echo "runner_exit_code=${runner_status}"
  echo "runner_wall_seconds=${runner_wall_seconds}"
  if [[ $runner_status -ne 0 ]]; then
    echo "normal_sigma_Id runner failed for d=${d}, nu=${nu}." >&2
    return "$runner_status"
  fi

  echo "===== normal_sigma_Id final validation: d=${d}, nu=${nu} ====="
  set +e
  Rscript --vanilla scripts/check_c3_normal_sigma_Id_t_production.R \
    --mode=final "--output_dir=${output_dir}" "--d=${d}" "--nu=${nu}" \
    --M=1000 --B=5000 --seed=20260728 --derivative_mc_size=10000 \
    --cvm_block_size=50 "--runner_wall_seconds=${runner_wall_seconds}" \
    "--allocated_cores=${SLURM_CPUS_PER_TASK}" --allow_new=false
  validation_status=$?
  set -e
  echo "validation_exit_code=${validation_status}"
  if [[ $validation_status -ne 0 ]]; then
    echo "normal_sigma_Id validation failed for d=${d}, nu=${nu}." >&2
    return "$validation_status"
  fi
  echo "===== normal_sigma_Id completed: d=${d}, nu=${nu} ====="
}

total_started_epoch=$(date +%s)
run_invocation 2 3 "$output_d2"
run_invocation 5 6 "$output_d5"
total_finished_epoch=$(date +%s)

echo "===== C3 normal_sigma_Id t production end ====="
echo "date=$(date --iso-8601=seconds)"
echo "total_wall_seconds=$((total_finished_epoch - total_started_epoch))"
echo "runner_exit_code=0"
echo "validation_exit_code=0"
echo "completed_invocations=2/2"
