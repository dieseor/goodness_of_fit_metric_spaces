#!/bin/bash

# Run on a C3 login node before submitting one restricted-spiked benchmark.
# This is read-only and does not create the output directory.

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash scripts/preflight_c3_restricted_spiked_normal_benchmark.sh {axis_075|axis_100|diagonal_150|diagonal_100}" >&2
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

if [[ ! -f scripts/run_restricted_spiked_normal_covariance_alternatives.R ||
      ! -f scripts/check_c3_restricted_spiked_normal.R ||
      ! -f bootstrap/restricted_spiked_normal_bootstrap.R ]]; then
  echo "Run this preflight from the repository root." >&2
  exit 2
fi

output_dir="simulation_results/restricted_spiked_normal_c3_benchmark/${mean_config}_d2_d5_n50_100_200_400_M20_B5000_Nderiv10000"

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

echo "===== Restricted-spiked C3 benchmark preflight ====="
echo "mean_config=${mean_config}"
echo "output_dir=${output_dir}"
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "pwd=$(pwd)"
echo "R_LIBS_USER=${R_LIBS_USER}"
module list 2>&1 || true
R --version

Rscript --vanilla scripts/check_c3_restricted_spiked_normal.R \
  --mode=preflight \
  "--mean_config=${mean_config}" \
  "--output_dir=${output_dir}" \
  --M=20 --B=5000 \
  --dimensions=2,5 --n_values=50,100,200,400 --beta_values=0,0.25,0.5,1 \
  --lambda=2 --derivative_mc_size=10000 --cvm_block_size=50 --allow_new=true
