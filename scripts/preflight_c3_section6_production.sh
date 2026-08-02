#!/bin/bash

# Run this on a C3 login node before sbatch. It catches missing unversioned
# result files immediately, without consuming a Slurm allocation.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash scripts/preflight_c3_section6_production.sh {normal|lg} [2,10|5]" >&2
  exit 2
fi

family="$1"
dimensions="${2:-2,10}"
case "$dimensions" in
  2,10)
    dimension_tag="d2_d10"
    allow_new_output="true"
    ;;
  5)
    dimension_tag="d5"
    allow_new_output="true"
    ;;
  *)
    echo "Unsupported dimensions '$dimensions'; expected '2,10' or '5'." >&2
    exit 2
    ;;
esac

case "$family" in
  normal)
    seed=20260728
    output_dir="simulation_results/section6_new_scenarios/final_normal_${dimension_tag}_n50_100_200_400_M1000_B5000_quadrature"
    ;;
  lg)
    seed=20260727
    output_dir="simulation_results/section6_new_scenarios/final_lg_${dimension_tag}_n50_100_200_400_M1000_B5000_quadrature"
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
echo "dimensions=${dimensions}"
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
  "--seed=${seed}" \
  "--dimensions=${dimensions}" \
  --derivative_method=quadrature \
  "--allow_new=${allow_new_output}"
