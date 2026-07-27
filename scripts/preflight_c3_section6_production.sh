#!/bin/bash

# Run this on a C3 login node before sbatch. It catches missing unversioned
# result files immediately, without consuming a Slurm allocation.

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash scripts/preflight_c3_section6_production.sh {normal|lg}" >&2
  exit 2
fi

family="$1"
case "$family" in
  normal)
    seed=20260728
    output_dir="simulation_results/section6_new_scenarios/final_normal_d2_d10_n50_100_200_400_M1000_B5000_auto"
    ;;
  lg)
    seed=20260727
    output_dir="simulation_results/section6_new_scenarios/final_lg_d2_d10_n50_100_200_400_M1000_B5000_auto"
    ;;
  *)
    echo "Unsupported family '$family'; expected 'normal' or 'lg'." >&2
    exit 2
    ;;
esac

if [[ ! -f scripts/run_section6_new_scenarios.R ||
      ! -f scripts/check_c3_section6_production.R ]]; then
  echo "Run this preflight from the repository root." >&2
  exit 2
fi

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

echo "===== C3 Section 6 production preflight ====="
echo "family=${family}"
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "pwd=$(pwd)"
echo "output_dir=${output_dir}"
echo "R_LIBS_USER=${R_LIBS_USER}"
module list 2>&1 || true
R --version

Rscript --vanilla scripts/check_c3_section6_production.R \
  --mode=preflight \
  "--family=${family}" \
  "--output_dir=${output_dir}" \
  "--seed=${seed}"
