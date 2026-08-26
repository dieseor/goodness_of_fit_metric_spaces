#!/bin/bash

#SBATCH --job-name=normal_Id_n800
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/normal_sigma_Id_t_n800_null_%j.out
#SBATCH --error=logs/normal_sigma_Id_t_n800_null_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:50:00
#SBATCH --mem=8G
#SBATCH --partition=batch

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

output_dir="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n800_M1000_B5000_beta0_Nderiv10000_cachetrue"

mkdir -p logs

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

echo "===== normal_sigma_Id n=800 null check ====="
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "commit=$(git rev-parse HEAD)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "output_dir=${output_dir}"

Rscript --vanilla scripts/run_normal_sigma_Id_t_pilot.R \
  "--output_dir=${output_dir}" \
  --M=1000 \
  --B=5000 \
  --dimensions=5 \
  --n_values=800 \
  --beta_values=0 \
  --nu=6 \
  --derivative_mc_size=10000 \
  --cvm_block_size=50 \
  --cache_corrections=true \
  "--cores=${SLURM_CPUS_PER_TASK}" \
  --checkpoint_results=500 \
  --seed=20260728 \
  --show_progress=true

Rscript --vanilla - "$output_dir" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[[1L]]

source("scripts/run_normal_sigma_Id_t_pilot.R")

x <- read.csv(file.path(output_dir, "raw_results.csv"),
              stringsAsFactors = FALSE)

ok <- conforming(x, 10000L)

if (nrow(x) != 1000L ||
    !all(ok) ||
    any(x$d != 5L) ||
    any(x$n != 800L) ||
    any(x$beta != 0)) {
  stop("Final n=800 validation failed.")
}

cat("Final validation passed: 1000/1000 conforming rows.\n")
RS

echo "===== END ====="
