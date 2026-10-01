#!/bin/bash
#SBATCH --job-name=paper_fast_matched
#SBATCH --account=uc3m-influentia-15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=2120M
#SBATCH --time=03:28:00
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: sbatch scripts/run_paper_fast_matched_batch.sbatch SCENARIO D REP_START REP_END" >&2
  exit 2
fi

: "${PAPER_FAST_MATCHED_OUTPUT_DIR:?Set PAPER_FAST_MATCHED_OUTPUT_DIR to a new timestamped directory}"
: "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
: "${SLURM_CPUS_PER_TASK:?Slurm did not report allocated CPUs}"
if [[ "$SLURM_CPUS_PER_TASK" -ne 16 ]]; then
  echo "Matched benchmark requires 16 CPUs per job." >&2
  exit 2
fi
cd "$SLURM_SUBMIT_DIR"

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_THREAD_LIMIT=1

printf 'job=%s host=%s cpus=%s output=%s args=%s\n' \
  "$SLURM_JOB_ID" "$(hostname)" "$SLURM_CPUS_PER_TASK" \
  "$PAPER_FAST_MATCHED_OUTPUT_DIR" "$*"
module list
env | grep -E '^(OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|BLIS_NUM_THREADS|VECLIB_MAXIMUM_THREADS|OMP_THREAD_LIMIT)='
Rscript --vanilla -e 'print(sessionInfo())'

Rscript --vanilla scripts/run_paper_fast_matched_batch.R "$@"
