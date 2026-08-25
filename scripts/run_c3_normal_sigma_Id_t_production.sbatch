#!/bin/bash

#SBATCH --job-name=normal_sigma_Id_t
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/normal_sigma_Id_t_%j.out
#SBATCH --error=logs/normal_sigma_Id_t_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:40:00
#SBATCH --mem=5G
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/run_c3_normal_sigma_Id_t_production.sbatch BETA" >&2
  exit 2
fi

beta_value="$1"
case "$beta_value" in
  0|0.25|0.5|1) ;;
  *)
    echo "BETA must be one of: 0, 0.25, 0.5, 1" >&2
    exit 2
    ;;
esac

beta_label="${beta_value/./p}"

: "${SLURM_SUBMIT_DIR:?Submit from repository root}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"

if [[ "$SLURM_CPUS_PER_TASK" -ne 16 ]]; then
  echo "Production requires exactly 16 CPUs." >&2
  exit 2
fi

cd "$SLURM_SUBMIT_DIR"

output_d2="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n50_100_200_400_M1000_B5000_beta${beta_label}_Nderiv10000"
output_d5="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n50_100_200_400_M1000_B5000_beta${beta_label}_Nderiv10000"

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "===== C3 normal_sigma_Id t production environment ====="
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "commit=$(git rev-parse HEAD)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "beta=${beta_value}"
echo "output_d2=${output_d2}"
echo "output_d5=${output_d5}"

echo "===== In-job preflight ====="
bash scripts/preflight_c3_normal_sigma_Id_t_production.sh "$beta_value"

run_invocation() {
  local d="$1"
  local nu="$2"
  local output_dir="$3"

  echo "===== runner start: d=${d}, nu=${nu}, beta=${beta_value} ====="

  started_epoch=$(date +%s)

  Rscript --vanilla scripts/run_normal_sigma_Id_t_pilot.R \
    "--output_dir=${output_dir}" \
    --M=1000 --B=5000 "--dimensions=${d}" \
    --n_values=50,100,200,400 "--beta_values=${beta_value}" \
    "--nu=${nu}" \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    "--cores=${SLURM_CPUS_PER_TASK}" \
    --checkpoint_results=2000 \
    --seed=20260728 \
    --show_progress=true

  finished_epoch=$(date +%s)
  runner_wall_seconds=$((finished_epoch - started_epoch))

  echo "===== final validation: d=${d}, nu=${nu}, beta=${beta_value} ====="

  Rscript --vanilla scripts/check_c3_normal_sigma_Id_t_production.R \
    --mode=final \
    "--output_dir=${output_dir}" \
    "--d=${d}" "--nu=${nu}" "--beta=${beta_value}" \
    --M=1000 --B=5000 \
    --seed=20260728 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    "--runner_wall_seconds=${runner_wall_seconds}" \
    "--allocated_cores=${SLURM_CPUS_PER_TASK}" \
    --allow_new=false
}

run_invocation 2 3 "$output_d2"
run_invocation 5 6 "$output_d5"

echo "===== completed beta=${beta_value} ====="
