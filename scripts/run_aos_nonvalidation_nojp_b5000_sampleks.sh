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

# 1. Comets: cardioid family (C1--C4 + uniform benchmark row already handled by the runner)
run_cmd Rscript -e '
source("scripts/run_comets_distance_profile_cardioid.R")
run_comets_distance_profile_cardioid(
  output_root = "real_data/comets/cardioid/paper_results_B5000_sampleks/fast",
  stages = c("oort_cvm", "oort_ks", "short_cvm", "short_ks"),
  cvm_B = as.integer(Sys.getenv("B", "5000")),
  ks_B = as.integer(Sys.getenv("B", "5000")),
  ks_grid_mode = "sample_points_unique_distances",
  n_cores = as.integer(Sys.getenv("N_CORES", "8")),
  bootstrap_method = "fast_multiplier",
  distance_type = "geodesic",
  control = list(cardioid_optim_control = list(maxit = 1000))
)'

# 2. Comets: small-circle model, short-period and long-period, KS and CvM
for dataset in short long; do
  for statistic in ks cvm; do
    if [[ "${dataset}" == "short" && "${statistic}" == "ks" ]]; then
      out_dir="real_data/comets/small_circle/short_comets_B5000_sampleks/fast"
    elif [[ "${dataset}" == "short" ]]; then
      out_dir="real_data/comets/small_circle/short_comets_B5000_samplecvm/fast"
    elif [[ "${statistic}" == "ks" ]]; then
      out_dir="real_data/comets/small_circle/long_comets_B5000_sampleks/fast"
    else
      out_dir="real_data/comets/small_circle/long_comets_B5000_samplecvm/fast"
    fi

    run_cmd Rscript scripts/run_comets_distance_profile_small_circle_benchmark.R \
      "--output_root=${out_dir}" \
      "--dataset=${dataset}" \
      "--B_values=${B}" \
      "--statistic=${statistic}" \
      "--n_cores=${N_CORES}" \
      "--bootstrap_method=fast_multiplier" \
      "--ks_grid_mode=sample_points_unique_distances" \
      "--distance_type=geodesic"
  done
done

# 3. Comets: beta-mixture model (exclude Jones--Pewsey and exclude uniform+beta here because it is not in the AoS paper table)
run_cmd Rscript scripts/run_comets_rotational_mixtures_short_long.R \
  "--output_root=real_data/comets/mixture/beta_mixture2_short_long_B5000/fast" \
  "--datasets=short,long" \
  "--models=beta_mixture2" \
  "--B=${B}" \
  "--statistics=ks,cvm" \
  "--n_cores=${N_CORES}" \
  "--bootstrap_method=fast_multiplier" \
  "--ks_grid_mode=sample_points_unique_distances"

# 4. Logistic Gaussian real-data screening
run_cmd Rscript scripts/run_logistic_gaussian_dataset_screening.R \
  "--output_dir=real_data/logistic_gaussian/screening/fast/paper_results_B5000_sampleks" \
  "--B=${B}" \
  "--n_cores=${N_CORES}" \
  "--bootstrap_mode=composite_multiplier" \
  "--omega_grid_type=sample_points" \
  "--t_grid_type=sample_distances"

# 5. Risoe 125m wind screening
run_cmd Rscript -e '
source("real_data/wind/run_risoe_125m_screening_ks_cvm.R")
run_risoe_125m_screening_ks_cvm(
  output_dir = "real_data/wind/month_diagnostics/screening_125m_ks_cvm_b5000/fast_sampleks",
  B = as.integer(Sys.getenv("B", "5000")),
  n_cores = as.integer(Sys.getenv("N_CORES", "8")),
  bootstrap_method = "fast_multiplier",
  ks_grid_mode = "sample_points_unique_distances"
)'

# 6. Sunspots rolling weighted-mixture GOF
run_cmd Rscript real_data/sunspots/run_sunspots_weighted_mixture_rolling_windows_gof.R \
  "--output_dir=real_data/sunspots/output/cycles21_23_weighted_mixture_rolling_windows_10cr_B5000/fast" \
  "--statistics=ks,cvm" \
  "--B=${B}" \
  "--n_cores=${N_CORES}" \
  "--bootstrap_method=fast_multiplier" \
  "--resume=TRUE"

echo
echo "All requested non-validation, non-JP AoS reruns have been launched."
