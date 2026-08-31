#!/bin/bash

#SBATCH --job-name=lg_dir15
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/lg_dir15_%j.out
#SBATCH --error=logs/lg_dir15_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=03:15:00
#SBATCH --mem=4G
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/run_c3_lg_dirichlet15_mu_only_production.sbatch BETA" >&2
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

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

output_dir="simulation_results/lg_dirichlet15_mu_only_final_beta${beta_label}_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000"

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

Rscript --vanilla \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
  "--output_dir=${output_dir}" \
  --scenario=dirichlet15 \
  --d=2,5 \
  --n=50,100,200,400 \
  "--beta=${beta_value}" \
  --M=1000 \
  --B=1000 \
  "--cores=${SLURM_CPUS_PER_TASK}" \
  --derivative_mc_size=10000 \
  --cvm_block_size=50 \
  --checkpoint_results=500 \
  --seed=20260904 \
  --show_progress=TRUE

Rscript --vanilla - "$output_dir" "$beta_value" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[[1L]]
beta_value <- as.numeric(args[[2L]])

source("scripts/run_logistic_gaussian_dirichlet15_mu_only.R")

x <- read.csv(
  file.path(output_dir, "raw_results.csv"),
  stringsAsFactors = FALSE
)

ok <- conforming(x, 10000L)

if (nrow(x) != 8000L)
  stop(sprintf("Expected 8000 rows, found %d.", nrow(x)))

if (!all(ok))
  stop(sprintf("Only %d/%d rows are conforming.", sum(ok), length(ok)))

if (any(x$scenario != "dirichlet15") ||
    any(!x$d %in% c(2L, 5L)) ||
    any(!x$n %in% c(50L, 100L, 200L, 400L)) ||
    any(abs(x$beta - beta_value) > 1e-12))
  stop("Final design validation failed.")

counts <- aggregate(replication ~ d + n, x, length)

if (nrow(counts) != 8L || any(counts$replication != 1000L))
  stop("Expected exactly 1000 replications per (d,n).")

cat(sprintf(
  "FINAL VALIDATION PASSED: %d/%d conforming rows.\n",
  sum(ok), length(ok)
))

print(read.csv(
  file.path(output_dir, "summary.csv"),
  stringsAsFactors = FALSE
))
RS
