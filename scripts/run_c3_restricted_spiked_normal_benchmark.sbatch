#!/bin/bash
#
# Intermediate benchmark for one restricted-spiked mean configuration.
# Run mkdir -p logs before sbatch because Slurm opens log files first.
#
#SBATCH --job-name=rspike_bench
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/rspike_bench_%j.out
#SBATCH --error=logs/rspike_bench_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/run_c3_restricted_spiked_normal_benchmark.sbatch {axis_075|axis_100|diagonal_150|diagonal_100}" >&2
  exit 2
fi

mean_config="$1"
case "$mean_config" in
  axis_075|axis_100|diagonal_150|diagonal_100) ;;
  *)
    echo "Unsupported mean configuration '$mean_config'." >&2
    exit 2
    ;;
esac

: "${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the repository root}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"
if [[ "$SLURM_CPUS_PER_TASK" -ne 32 ]]; then
  echo "Benchmark requires exactly 32 allocated CPUs." >&2
  exit 2
fi
cd "$SLURM_SUBMIT_DIR"

output_dir="simulation_results/restricted_spiked_normal_c3_benchmark/${mean_config}_d2_d5_n50_100_200_400_M20_B5000_Nderiv10000"
if [[ -d "$output_dir" && ! -f "$output_dir/manifest.csv" ]]; then
  echo "Refusing nonempty or unrecognized benchmark directory: $output_dir" >&2
  exit 2
fi

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

echo "===== Restricted-spiked benchmark environment ====="
echo "mean_config=${mean_config}"
echo "output_dir=${output_dir}"
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "pwd=$(pwd)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-not_set}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-not_set}"
module list 2>&1 || true
R --version
scontrol show job "$SLURM_JOB_ID" || true

Rscript --vanilla scripts/check_c3_restricted_spiked_normal.R \
  --mode=preflight \
  "--mean_config=${mean_config}" \
  "--output_dir=${output_dir}" \
  --M=20 --B=5000 \
  --dimensions=2,5 --n_values=50,100,200,400 --beta_values=0,0.25,0.5,1 \
  --lambda=2 --derivative_mc_size=10000 --cvm_block_size=50 --allow_new=true

started_epoch=$(date +%s)
set +e
Rscript --vanilla scripts/run_restricted_spiked_normal_covariance_alternatives.R \
  "--mean_config=${mean_config}" \
  "--output_dir=${output_dir}" \
  --M=20 --B=5000 \
  --dimensions=2,5 --n_values=50,100,200,400 --beta_values=0,0.25,0.5,1 \
  --lambda=2 --derivative_mc_size=10000 --cvm_block_size=50 \
  "--cores=${SLURM_CPUS_PER_TASK}" --checkpoint_results=64 --show_progress=true
runner_status=$?
set -e
finished_epoch=$(date +%s)
runner_wall_seconds=$((finished_epoch - started_epoch))
if [[ $runner_status -ne 0 ]]; then
  echo "runner_exit_code=${runner_status}" >&2
  exit "$runner_status"
fi

set +e
Rscript --vanilla scripts/check_c3_restricted_spiked_normal.R \
  --mode=final \
  "--mean_config=${mean_config}" \
  "--output_dir=${output_dir}" \
  --M=20 --B=5000 \
  --dimensions=2,5 --n_values=50,100,200,400 --beta_values=0,0.25,0.5,1 \
  --lambda=2 --derivative_mc_size=10000 --cvm_block_size=50 \
  "--runner_wall_seconds=${runner_wall_seconds}" \
  "--allocated_cores=${SLURM_CPUS_PER_TASK}"
validation_status=$?
set -e

echo "===== Restricted-spiked benchmark end ====="
echo "mean_config=${mean_config}"
echo "runner_exit_code=${runner_status}"
echo "validation_exit_code=${validation_status}"
echo "runner_wall_seconds=${runner_wall_seconds}"
echo "output_dir=${output_dir}"
exit "$validation_status"
