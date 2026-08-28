#!/bin/bash

#SBATCH --job-name=lg_ar1_n800
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/lg_t4_ar1_n800_%j.out
#SBATCH --error=logs/lg_t4_ar1_n800_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/run_c3_lg_t4_ar1_n800_production.sbatch BETA" >&2
  exit 2
fi

beta_value="$1"

case "$beta_value" in
  0|0.5|1) ;;
  *)
    echo "BETA must be one of: 0, 0.5, 1" >&2
    exit 2
    ;;
esac

beta_label="${beta_value/./p}"

: "${SLURM_SUBMIT_DIR:?Submit from repository root}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

output_dir="simulation_results/lg_t4_ar1_final_beta${beta_label}_d2_d5_n800_M1000_B1000_Nderiv10000"

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

echo "===== C3 Logistic-Gaussian t4 AR(1), n=800 ====="
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "commit=$(git rev-parse HEAD)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "beta=${beta_value}"
echo "output_dir=${output_dir}"

started_epoch=$(date +%s)

Rscript --vanilla \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
  "--output_dir=${output_dir}" \
  --scenario=t4 \
  --d=2,5 \
  --n=800 \
  "--beta=${beta_value}" \
  --M=1000 \
  --B=1000 \
  --nu=4 \
  --t_standardized=TRUE \
  --t_ar1_rho=0.5 \
  "--cores=${SLURM_CPUS_PER_TASK}" \
  --derivative_mc_size=10000 \
  --cvm_block_size=50 \
  --checkpoint_results=250 \
  --seed=20260903 \
  --show_progress=TRUE

finished_epoch=$(date +%s)
runner_wall_seconds=$((finished_epoch - started_epoch))

echo "===== final validation ====="

Rscript --vanilla - "$output_dir" "$beta_value" "$runner_wall_seconds" <<'RS'
args <- commandArgs(trailingOnly = TRUE)

output_dir <- args[[1L]]
beta_value <- as.numeric(args[[2L]])
runner_wall_seconds <- as.numeric(args[[3L]])

source("scripts/run_logistic_gaussian_sigma_shape_scenarios.R")

x <- read.csv(
  file.path(output_dir, "raw_results.csv"),
  stringsAsFactors = FALSE
)

expected_rows <- 2L * 1000L
ok <- conforming(x, 10000L)

if (nrow(x) != expected_rows) {
  stop(sprintf(
    "Expected %d rows, found %d.",
    expected_rows, nrow(x)
  ))
}

if (!all(ok)) {
  stop(sprintf(
    "Only %d/%d rows are conforming.",
    sum(ok), length(ok)
  ))
}

if (any(x$scenario != "t4") ||
    any(!x$d %in% c(2L, 5L)) ||
    any(x$n != 800L) ||
    any(x$beta != beta_value)) {
  stop("Final n=800 design validation failed.")
}

if ("rho_hat" %in% names(x)) {
  if (any(!is.finite(x$rho_hat))) {
    stop("Non-finite rho_hat found.")
  }
  if (any(abs(x$rho_hat) >= 0.995 + 1e-12)) {
    stop("rho_hat outside the fitted parameter interval.")
  }
}

cat(sprintf(
  "FINAL n=800 VALIDATION PASSED: %d/%d conforming rows.\n",
  sum(ok), length(ok)
))
cat(sprintf(
  "Runner wall time: %.1f minutes.\n",
  runner_wall_seconds / 60
))

print(
  read.csv(
    file.path(output_dir, "summary.csv"),
    stringsAsFactors = FALSE
  )
)
RS

echo "===== completed n=800 beta=${beta_value} ====="
