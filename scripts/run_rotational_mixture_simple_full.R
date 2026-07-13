source(file.path("bootstrap", "calibration_study.R"))

rotmix_simple_expected_raw_rows <- function(M_outer,
                                            n_values,
                                            statistics,
                                            n_scenarios) {
  as.integer(M_outer) * length(as.integer(n_values)) * length(statistics) * as.integer(n_scenarios)
}

rotmix_simple_expected_summary_rows <- function(n_values,
                                                statistics,
                                                alphas,
                                                n_scenarios) {
  length(as.integer(n_values)) * length(statistics) * length(as.numeric(alphas)) * as.integer(n_scenarios)
}

rotmix_simple_output_complete <- function(output_dir,
                                          expected_raw_rows,
                                          expected_summary_rows) {
  raw_csv <- file.path(output_dir, "bootstrap_calibration_raw.csv")
  summary_csv <- file.path(output_dir, "bootstrap_calibration_summary.csv")
  if (!file.exists(raw_csv) || !file.exists(summary_csv)) {
    return(FALSE)
  }

  raw_df <- try(utils::read.csv(raw_csv, stringsAsFactors = FALSE), silent = TRUE)
  summary_df <- try(utils::read.csv(summary_csv, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(raw_df, "try-error") || inherits(summary_df, "try-error")) {
    return(FALSE)
  }

  nrow(raw_df) >= expected_raw_rows && nrow(summary_df) >= expected_summary_rows
}

run_rotational_mixture_simple_full <- function(
  output_root = file.path("output", "calibration", "bootstrap", "rotational_mixtures_simple_full"),
  n_values = c(50L, 100L, 200L),
  M_outer = 500L,
  B = 5000L,
  statistics = c("ks", "cvm"),
  alphas = c(0.01, 0.05, 0.10),
  n_cores_outer = 12L,
  seed = 20260601L
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  scenarios <- list(
    default_beta_mixture2_simple_calibration_scenarios()[[1L]],
    default_logitnormal_mixture2_simple_calibration_scenarios()[[1L]]
  )

  expected_raw <- rotmix_simple_expected_raw_rows(
    M_outer = M_outer,
    n_values = n_values,
    statistics = statistics,
    n_scenarios = length(scenarios)
  )
  expected_summary <- rotmix_simple_expected_summary_rows(
    n_values = n_values,
    statistics = statistics,
    alphas = alphas,
    n_scenarios = length(scenarios)
  )

  if (rotmix_simple_output_complete(output_root, expected_raw, expected_summary)) {
    message("[SKIP] Full simple calibration already complete at ", output_root)
    return(invisible(list(output_root = output_root, skipped = TRUE)))
  }

  result <- run_bootstrap_calibration_study(
    scenarios = scenarios,
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
  run_rotational_mixture_simple_full()
}
