#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

# Respect the project renv so repo-declared packages are visible to each Rscript.

B="${B:-5000}"
N_CORES="${N_CORES:-8}"

run_cmd() {
  echo
  echo ">>> $*"
  "$@"
}

# Pending AoS reruns only.
# Explicitly excluded here:
# - giant H0 calibration table
# - Jones--Pewsey experiments
# - convergence / validation figures
# - Bahadur validation
# - rotational mixtures: beta_mixture2 and uniform_beta_mixture
# - items already completed in the previous rerun batch

# 1. Logistic Gaussian real-data screening
run_cmd Rscript scripts/run_logistic_gaussian_dataset_screening.R \
  "--output_dir=real_data/logistic_gaussian/screening/fast/paper_results_B5000_sampleks" \
  "--B=${B}" \
  "--n_cores=${N_CORES}" \
  "--bootstrap_mode=composite_multiplier" \
  "--omega_grid_type=sample_points" \
  "--t_grid_type=sample_distances"

# 2. Risoe 125m wind screening
run_cmd Rscript -e '
source("real_data/wind/run_risoe_125m_screening_ks_cvm.R")
run_risoe_125m_screening_ks_cvm(
  output_dir = "real_data/wind/month_diagnostics/screening_125m_ks_cvm_b5000/fast_sampleks",
  B = as.integer(Sys.getenv("B", "5000")),
  n_cores = as.integer(Sys.getenv("N_CORES", "8")),
  bootstrap_method = "fast_multiplier",
  ks_grid_mode = "sample_points_unique_distances"
)'

# 3. Sunspots rolling weighted-mixture GOF
run_cmd Rscript real_data/sunspots/run_sunspots_weighted_mixture_rolling_windows_gof.R \
  "--output_dir=real_data/sunspots/output/cycles21_23_weighted_mixture_rolling_windows_10cr_B5000/fast" \
  "--statistics=ks,cvm" \
  "--B=${B}" \
  "--n_cores=${N_CORES}" \
  "--bootstrap_method=fast_multiplier" \
  "--resume=TRUE"

echo
echo "Pending non-validation AoS reruns excluding rotational mixtures have been launched."
