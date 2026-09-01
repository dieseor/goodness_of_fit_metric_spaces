#!/bin/bash
# C3 pilot for one HvMF Section 6 scenario.
# Usage: sbatch scripts/run_c3_hvmf_section6_pilot.sbatch {hvmf_1_mixture|hvmf_2_angular|hvmf_1_dimension_scaled_location_mixture|hvmf_2_angular_sqrt_d_concentration} [2,5] [0,0.5]
# Restartable: the Section 6 runner reuses completed rows.
#SBATCH --job-name=section6_hvmf_pilot
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/section6_hvmf_pilot_%j.out
#SBATCH --error=logs/section6_hvmf_pilot_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:30:00
#SBATCH --mem=5G
#SBATCH --partition=batch

set -euo pipefail
if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: sbatch scripts/run_c3_hvmf_section6_pilot.sbatch {hvmf_1_mixture|hvmf_2_angular|hvmf_1_dimension_scaled_location_mixture|hvmf_2_angular_sqrt_d_concentration} [2,5] [0,0.5]" >&2
  exit 2
fi
scenario="$1"
dimensions="${2:-2,5}"
beta_values="${3:-0,0.5}"
case "$scenario" in
  hvmf_1_mixture|hvmf_2_angular|hvmf_1_dimension_scaled_location_mixture|hvmf_2_angular_sqrt_d_concentration) ;;
  *) echo "Unsupported scenario '$scenario'." >&2; exit 2 ;;
esac
case "$dimensions" in
  2|5|2,5) ;;
  *) echo "Unsupported dimensions '$dimensions'; expected 2, 5, or 2,5." >&2; exit 2 ;;
esac
case "$beta_values" in
  0|0.5|0,0.5) ;;
  *) echo "Unsupported beta grid '$beta_values'; expected 0, 0.5, or 0,0.5." >&2; exit 2 ;;
esac

seed=20260830
if [[ "$beta_values" == "0,0.5" ]]; then
  # Preserve the historical default location for the existing two pilots.
  output_dir="simulation_results/section6_new_scenarios/pilot_${scenario}_d${dimensions//,/_}_n50_100_200_400_M100_B1000_quadrature_F_integral"
else
  beta_tag="${beta_values//./p}"
  beta_tag="${beta_tag//,/__}"
  output_dir="simulation_results/section6_new_scenarios/pilot_${scenario}_d${dimensions//,/_}_n50_100_200_400_M100_B1000_beta${beta_tag}_quadrature_F_integral"
fi
: "${SLURM_SUBMIT_DIR:?Submit with sbatch from the repository root.}"
: "${SLURM_CPUS_PER_TASK:?Slurm did not provide a CPU allocation.}"
if [[ "$SLURM_CPUS_PER_TASK" -ne 16 ]]; then
  echo "This pilot requires 16 CPUs; Slurm assigned ${SLURM_CPUS_PER_TASK}." >&2
  exit 2
fi
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs "$(dirname "$output_dir")"

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1
export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

echo "scenario=${scenario} dimensions=${dimensions} output_dir=${output_dir}"
echo "M=100 B=1000 n=50,100,200,400 beta=${beta_values} derivative_method=quadrature"
echo "host=$(hostname) job=${SLURM_JOB_ID:-not_set} date=$(date --iso-8601=seconds)"
Rscript --vanilla scripts/run_section6_new_scenarios.R \
  "--family=hvmf" "--scenarios=${scenario}" "--output_dir=${output_dir}" \
  --M=100 --B=1000 "--dimensions=${dimensions}" \
  --n_values=50,100,200,400 "--beta_values=${beta_values}" \
  --cores="${SLURM_CPUS_PER_TASK}" --checkpoint_results=400 \
  "--seed=${seed}" --derivative_method=quadrature --derivative_mc_size=10000 \
  --cvm_block_size=50 --show_progress=true
