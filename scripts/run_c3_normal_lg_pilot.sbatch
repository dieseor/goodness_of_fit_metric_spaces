#!/bin/bash
#
# IMPORTANT: Slurm opens --output/--error before this script starts.  Submit
# from the repository root and run `mkdir -p logs` before calling sbatch.
#
#SBATCH --job-name=normal_lg_pilot
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/normal_lg_pilot_%j.out
#SBATCH --error=logs/normal_lg_pilot_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mem=16G

set -euo pipefail

: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set; submit this file with sbatch}"
cd "$SLURM_SUBMIT_DIR"

# This is still useful for auxiliary files, but it is intentionally not relied
# upon for the Slurm output/error paths, which were opened before execution.
mkdir -p logs

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

# The production Section 6 runner uses up to SLURM_CPUS_PER_TASK outer forked
# workers and forces n_cores=1 inside every bootstrap.  Limit threaded numerical
# libraries so that four workers cannot each create additional BLAS/OpenMP
# threads.
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "===== C3 normal/LG pilot environment ====="
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "pwd=$(pwd)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-not_set}"
echo "SLURM_JOB_NAME=${SLURM_JOB_NAME:-not_set}"
echo "SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-not_set}"
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
echo "===== Pilot start ====="

set +e
Rscript --vanilla scripts/run_c3_normal_lg_pilot.R
r_status=$?
set -e

echo "===== Pilot end ====="
echo "date=$(date --iso-8601=seconds)"
echo "R_exit_code=${r_status}"
exit "$r_status"
