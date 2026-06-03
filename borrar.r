

#!/usr/bin/env Rscript

set.seed(123)

cat("\n============================================================\n")
cat("STEP 1/2: Real-data comets, beta_mixture2, short/long\n")
cat("============================================================\n")

source("scripts/run_comets_rotational_mixtures_short_long.R")

res_comets_beta <- run_comets_mixtures_short_long(
  output_root = file.path(
    "output", "real_data", "comets",
    "beta_mixture2_short_long_B1000_4cores"
  ),
  datasets = c("short", "long"),
  models = c("beta_mixture2"),
  B = 1000L,
  statistics = c("ks", "cvm"),
  n_cores = 4L,
  seed = 123L,
  M_value = 60L,
  ks_t_points = 200L,
  distance_type = "geodesic"
)

cat("\nFinished comets run. Summary saved in:\n")
cat(file.path(
  "output", "real_data", "comets",
  "beta_mixture2_short_long_B1000_4cores",
  "comets_rotational_mixtures_summary.csv"
), "\n")

gc(verbose = FALSE)

cat("\n============================================================\n")
cat("STEP 2/2: Beta-mixture composite calibration, M=B=250\n")
cat("============================================================\n")

source("bootstrap/calibration_study.R")

scenarios_beta_comp <- list(
  make_beta_mixture2_composite_calibration_scenario(
    weight1 = 0.4,
    alpha1 = 2,
    beta1 = 8,
    alpha2 = 8,
    beta2 = 2
  ),
  make_beta_mixture2_composite_calibration_scenario(
    weight1 = 0.55,
    alpha1 = 4,
    beta1 = 12,
    alpha2 = 10,
    beta2 = 3
  )
)

for (i in seq_along(scenarios_beta_comp)) {
  scenarios_beta_comp[[i]]$control <- modifyList(
    scenarios_beta_comp[[i]]$control %||% list(),
    list(
      beta_mixture2_profile_method = "legendre",
      beta_mixture2_quad_n = 100L
    )
  )
}

res_beta_comp_multi <- run_bootstrap_calibration_study(
  scenarios = scenarios_beta_comp,
  n_values = c(50, 100, 200),
  M_outer = 250,
  B = 250,
  alpha_nominal = 0.05,
  alphas = c(0.01, 0.05, 0.10),
  statistics = c("ks", "cvm"),
  n_cores_outer = 12,
  seed = 123,
  output_dir = file.path(
    "output", "calibration", "bootstrap", "mixtures",
    "beta_mixture2_composite_custom_multi_M250_B250_n50_100_200"
  ),
  show_progress = TRUE,
  verbose = TRUE
)

cat("\nFinished calibration run. Summary saved in:\n")
cat(file.path(
  "output", "calibration", "bootstrap", "mixtures",
  "beta_mixture2_composite_custom_multi_M250_B250_n50_100_200",
  "bootstrap_calibration_summary.csv"
), "\n")

cat("\nAll overnight beta-mixture jobs finished.\n")