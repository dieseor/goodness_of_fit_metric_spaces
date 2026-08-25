#!/bin/bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash scripts/preflight_c3_normal_sigma_Id_t_production.sh BETA" >&2
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

required_files=(
  bootstrap/model_specs.R
  bootstrap/normal_sigma_Id_model_spec.R
  bootstrap/normal_sigma_Id_bootstrap.R
  scripts/run_normal_sigma_Id_t_pilot.R
  scripts/check_c3_normal_sigma_Id_t_production.R
  scripts/preflight_c3_normal_sigma_Id_t_production.sh
  scripts/run_c3_normal_sigma_Id_t_production.sbatch
)

for path in "${required_files[@]}"; do
  if [[ ! -f "$path" ]]; then
    echo "Missing required deployed file: $path" >&2
    exit 2
  fi
  if ! git cat-file -e "HEAD:${path}" 2>/dev/null; then
    echo "Required file is not committed in deployed HEAD: $path" >&2
    exit 2
  fi
done

if ! git diff --quiet HEAD -- "${required_files[@]}"; then
  echo "A required normal_sigma_Id workflow file differs from deployed HEAD." >&2
  git status --short -- "${required_files[@]}" >&2
  exit 2
fi

head_commit=$(git rev-parse HEAD)
origin_commit=$(git rev-parse origin/main)

if [[ "$head_commit" != "$origin_commit" ]]; then
  echo "HEAD does not match origin/main." >&2
  echo "HEAD=$head_commit" >&2
  echo "origin/main=$origin_commit" >&2
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

output_d2="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n50_100_200_400_M1000_B5000_beta${beta_label}_Nderiv10000"
output_d5="simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n50_100_200_400_M1000_B5000_beta${beta_label}_Nderiv10000"

echo "===== C3 normal_sigma_Id t production preflight ====="
echo "commit=${head_commit}"
echo "beta=${beta_value}"
echo "output_d2=${output_d2}"
echo "output_d5=${output_d5}"

Rscript --vanilla scripts/check_c3_normal_sigma_Id_t_production.R \
  --mode=preflight "--output_dir=${output_d2}" --d=2 --nu=3 \
  "--beta=${beta_value}" \
  --M=1000 --B=5000 --seed=20260728 --derivative_mc_size=10000 \
  --cvm_block_size=50 --allow_new=true

Rscript --vanilla scripts/check_c3_normal_sigma_Id_t_production.R \
  --mode=preflight "--output_dir=${output_d5}" --d=5 --nu=6 \
  "--beta=${beta_value}" \
  --M=1000 --B=5000 --seed=20260728 --derivative_mc_size=10000 \
  --cvm_block_size=50 --allow_new=true

echo "normal_sigma_Id beta=${beta_value} preflights passed."
