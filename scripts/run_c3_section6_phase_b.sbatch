#!/bin/bash
#
# Section 6 Phase B benchmark for one family. Submit this file twice, once with
# `normal` and once with `lg`. Slurm opens --output/--error before the script
# starts, so run `mkdir -p logs` from the repository root before calling sbatch.
#
# Production configuration retained:
#   B = 5000, dimensions = c(2, 10), n = c(50, 100, 200, 400),
#   beta = c(0, 0.5, 1), derivative_mc_size = 1000,
#   cvm_block_size = 50, both scenarios, fused C++ fast bootstrap.
# Only M is reduced from 1000 to 20:
#   2 scenarios x 2 dimensions x 4 n x 3 beta x 20 = 960 replications.
#
#SBATCH --job-name=section6_phase_b
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/section6_phase_b_%j.out
#SBATCH --error=logs/section6_phase_b_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=64G

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/run_c3_section6_phase_b.sbatch {normal|lg}" >&2
  exit 2
fi

family="$1"
case "$family" in
  normal)
    seed=20260728
    output_dir="simulation_results/section6_new_scenarios/c3_phase_b/normal_full_M20_B5000/job_${SLURM_JOB_ID}"
    completed_results=""
    ;;
  lg)
    seed=20260727
    output_dir="simulation_results/section6_new_scenarios/c3_phase_b/lg_full_M20_B5000/job_${SLURM_JOB_ID}"
    completed_results="simulation_results/section6_new_scenarios/final_lg_d2_d10_n50_100_200_400_M1000_B5000_auto/raw_results.csv"
    ;;
  *)
    echo "Unsupported family '$family'; expected 'normal' or 'lg'." >&2
    exit 2
    ;;
esac

: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set; submit this file with sbatch}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"
if [[ "$SLURM_CPUS_PER_TASK" -ne 32 ]]; then
  echo "Phase B requires exactly 32 CPUs; Slurm assigned ${SLURM_CPUS_PER_TASK}." >&2
  exit 2
fi
cd "$SLURM_SUBMIT_DIR"

# This does not create the Slurm log paths retroactively. The submitter must
# create logs/ before sbatch, as noted above.
mkdir -p logs
mkdir -p "$(dirname "$output_dir")"

if [[ -e "$output_dir" ]]; then
  echo "Refusing to reuse Phase B output directory: ${output_dir}" >&2
  exit 2
fi
if [[ -n "$completed_results" && ! -f "$completed_results" ]]; then
  echo "LG completed-results file is missing: ${completed_results}" >&2
  exit 2
fi

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

echo "===== C3 Section 6 Phase B environment ====="
echo "family=${family}"
echo "output_dir=${output_dir}"
echo "completed_results=${completed_results:-none}"
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
echo "===== Phase B runner start ====="

runner_started_epoch=$(date +%s)
set +e
Rscript --vanilla scripts/run_section6_new_scenarios.R \
  "--family=${family}" \
  "--output_dir=${output_dir}" \
  --M=20 \
  --B=5000 \
  --dimensions=2,10 \
  --n_values=50,100,200,400 \
  --beta_values=0,0.5,1 \
  "--cores=${SLURM_CPUS_PER_TASK}" \
  --checkpoint_results=64 \
  "--seed=${seed}" \
  --derivative_mc_size=1000 \
  --cvm_block_size=50 \
  --show_progress=true
runner_status=$?
set -e
runner_finished_epoch=$(date +%s)
runner_wall_seconds=$((runner_finished_epoch - runner_started_epoch))

if [[ $runner_status -ne 0 ]]; then
  echo "Phase B runner failed with exit code ${runner_status}." >&2
  exit "$runner_status"
fi
if [[ $runner_wall_seconds -le 0 ]]; then
  echo "Invalid Phase B runner wall time: ${runner_wall_seconds}." >&2
  exit 2
fi

checker_args=(
  "--family=${family}"
  "--output_dir=${output_dir}"
  "--seed=${seed}"
  "--cores=${SLURM_CPUS_PER_TASK}"
  "--runner_wall_seconds=${runner_wall_seconds}"
)
if [[ -n "$completed_results" ]]; then
  checker_args+=("--completed_results=${completed_results}")
fi

echo "===== Phase B validation and estimation start ====="
set +e
Rscript --vanilla scripts/check_c3_section6_phase_b.R "${checker_args[@]}"
validation_status=$?
set -e

echo "===== C3 Section 6 Phase B end ====="
echo "family=${family}"
echo "date=$(date --iso-8601=seconds)"
echo "runner_exit_code=${runner_status}"
echo "validation_exit_code=${validation_status}"
echo "runner_wall_seconds=${runner_wall_seconds}"
echo "output_dir=${output_dir}"

exit "$validation_status"
