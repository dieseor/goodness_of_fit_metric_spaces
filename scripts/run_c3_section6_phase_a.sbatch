#!/bin/bash
#
# Section 6 Phase A benchmark for one family. Submit this file twice, once with
# `normal` and once with `lg`. Slurm opens --output/--error before the script
# starts, so run `mkdir -p logs` from the repository root before calling sbatch.
#
# Production values retained:
#   B = 5000, derivative_mc_size = 1000, cvm_block_size = 50,
#   beta = c(0, 0.5, 1), both production scenarios, fused C++ fast bootstrap.
# Reduced benchmark values:
#   M = 2, dimensions = 10, n = 400.
# This gives 2 scenarios x 3 beta values x 2 replications = 12 jobs per family.
#
#SBATCH --job-name=section6_phase_a
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/section6_phase_a_%j.out
#SBATCH --error=logs/section6_phase_a_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/run_c3_section6_phase_a.sbatch {normal|lg}" >&2
  exit 2
fi

family="$1"
case "$family" in
  normal)
    seed=20260728
    output_dir="simulation_results/section6_new_scenarios/c3_phase_a/normal_d10_n400_M2_B5000"
    ;;
  lg)
    seed=20260727
    output_dir="simulation_results/section6_new_scenarios/c3_phase_a/lg_d10_n400_M2_B5000"
    ;;
  *)
    echo "Unsupported family '$family'; expected 'normal' or 'lg'." >&2
    exit 2
    ;;
esac

: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set; submit this file with sbatch}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"
cd "$SLURM_SUBMIT_DIR"

# This does not create the Slurm log paths retroactively. The submitter must
# create logs/ before sbatch, as noted above.
mkdir -p logs
mkdir -p "$(dirname "$output_dir")"

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

# Section 6 parallelism is the outer mcparallel queue. Every bootstrap uses
# n_cores=1, so threaded numerical libraries must remain serial.
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "===== C3 Section 6 Phase A environment ====="
echo "family=${family}"
echo "output_dir=${output_dir}"
echo "seed=${seed}"
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "pwd=$(pwd)"
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
echo "===== Loaded modules ====="
module list 2>&1 || true
echo "===== R version ====="
R --version
echo "===== Slurm allocation ====="
scontrol show job "$SLURM_JOB_ID" || true
echo "===== Phase A runner start ====="

phase_a_started_epoch=$(date +%s)
set +e
Rscript --vanilla scripts/run_section6_new_scenarios.R \
  "--family=${family}" \
  "--output_dir=${output_dir}" \
  --M=2 \
  --B=5000 \
  --dimensions=10 \
  --n_values=400 \
  --beta_values=0,0.5,1 \
  "--cores=${SLURM_CPUS_PER_TASK}" \
  --checkpoint_results=4 \
  "--seed=${seed}" \
  --derivative_mc_size=1000 \
  --cvm_block_size=50 \
  --show_progress=true
runner_status=$?
set -e

if [[ $runner_status -ne 0 ]]; then
  echo "Phase A runner failed with exit code ${runner_status}." >&2
  exit "$runner_status"
fi

echo "===== Phase A validation start ====="
set +e
Rscript --vanilla scripts/check_c3_section6_phase_a.R \
  "--family=${family}" \
  "--output_dir=${output_dir}" \
  "--seed=${seed}"
validation_status=$?
set -e

phase_a_finished_epoch=$(date +%s)
phase_a_wall_seconds=$((phase_a_finished_epoch - phase_a_started_epoch))

echo "===== C3 Section 6 Phase A end ====="
echo "family=${family}"
echo "date=$(date --iso-8601=seconds)"
echo "runner_exit_code=${runner_status}"
echo "validation_exit_code=${validation_status}"
echo "phase_a_wall_seconds=${phase_a_wall_seconds}"
echo "output_dir=${output_dir}"

exit "$validation_status"
