source(file.path("bootstrap", "calibration_study.R"))

run_beta_composite_fast_20260602 <- function(
  output_root = file.path("output", "bootstrap_calibration", "rotational_beta_composite_fast_20260602"),
  n_values = c(50L, 100L),
  M_outer = 500L,
  B = 250L,
  statistics = c("ks", "cvm"),
  alphas = c(0.01, 0.05, 0.10),
  n_cores_outer = 12L,
  seed = 20260602L
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  scenario <- default_rotational_beta_mixture2_composite_calibration_scenarios()[[1L]]

  result <- run_bootstrap_calibration_study(
    scenarios = list(scenario),
    n_values = as.integer(n_values),
    M_outer = as.integer(M_outer),
    B = as.integer(B),
    alpha_nominal = 0.05,
    alphas = as.numeric(alphas),
    statistics = statistics,
    n_cores_outer = as.integer(n_cores_outer),
    seed = as.integer(seed),
    output_dir = output_root,
    show_progress = FALSE,
    verbose = TRUE
  )

  saveRDS(result, file = file.path(output_root, "run_result.rds"))
  invisible(result)
}

if (sys.nframe() == 0L) {
  run_beta_composite_fast_20260602()
}
