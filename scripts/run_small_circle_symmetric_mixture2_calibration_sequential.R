resolve_scsm2_calibration_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

calibration_study_script_scsm2 <- resolve_scsm2_calibration_path("bootstrap", "calibration_study.R")
source(calibration_study_script_scsm2)

scsm2_expected_raw_rows <- function(M_outer, n_values, statistics) {
  as.integer(M_outer) * length(as.integer(n_values)) * length(statistics)
}

scsm2_expected_summary_rows <- function(n_values, statistics, alphas) {
  length(as.integer(n_values)) * length(statistics) * length(as.numeric(alphas))
}

scsm2_is_scenario_output_complete <- function(output_dir, expected_raw_rows, expected_summary_rows) {
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

scsm2_run_single_scenario_with_checkpoint <- function(scenario,
                                                      output_dir,
                                                      n_values,
                                                      M_outer,
                                                      B,
                                                      statistics,
                                                      alphas,
                                                      n_cores_outer,
                                                      seed,
                                                      show_progress = TRUE,
                                                      verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  expected_raw <- scsm2_expected_raw_rows(M_outer = M_outer, n_values = n_values, statistics = statistics)
  expected_summary <- scsm2_expected_summary_rows(n_values = n_values, statistics = statistics, alphas = alphas)

  if (scsm2_is_scenario_output_complete(output_dir, expected_raw, expected_summary)) {
    message(sprintf("[SKIP] Scenario '%s' already complete at %s", scenario$id, output_dir))
    return(list(
      scenario_id = scenario$id,
      output_dir = output_dir,
      skipped = TRUE,
      raw_csv = file.path(output_dir, "bootstrap_calibration_raw.csv"),
      summary_csv = file.path(output_dir, "bootstrap_calibration_summary.csv")
    ))
  }

  message(sprintf("[RUN ] Scenario '%s' -> %s", scenario$id, output_dir))
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
    output_dir = output_dir,
    show_progress = isTRUE(show_progress),
    verbose = isTRUE(verbose)
  )

  list(
    scenario_id = scenario$id,
    output_dir = output_dir,
    skipped = FALSE,
    raw_csv = result$raw_csv,
    summary_csv = result$summary_csv
  )
}

scsm2_consolidate_sequential_outputs <- function(run_rows, consolidated_dir) {
  dir.create(consolidated_dir, recursive = TRUE, showWarnings = FALSE)

  summary_all <- do.call(rbind, lapply(run_rows$summary_csv, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$source_summary_csv <- path
    x
  }))
  raw_all <- do.call(rbind, lapply(run_rows$raw_csv, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$source_raw_csv <- path
    x
  }))

  summary_out <- file.path(consolidated_dir, "bootstrap_calibration_summary_all_scenarios.csv")
  raw_out <- file.path(consolidated_dir, "bootstrap_calibration_raw_all_scenarios.csv")
  manifest_out <- file.path(consolidated_dir, "sequential_manifest.csv")

  utils::write.csv(summary_all, summary_out, row.names = FALSE)
  utils::write.csv(raw_all, raw_out, row.names = FALSE)
  utils::write.csv(run_rows, manifest_out, row.names = FALSE)

  list(summary_csv = summary_out, raw_csv = raw_out, manifest_csv = manifest_out)
}

default_scsm2_selected_scenarios <- function(run_simple = TRUE,
                                             run_composite = TRUE,
                                             kappa = 12,
                                             nu = 0.6) {
  scenarios <- list()
  if (isTRUE(run_simple)) {
    scenarios[[length(scenarios) + 1L]] <- make_small_circle_symmetric_mixture2_simple_calibration_scenario(
      kappa = kappa,
      nu = nu
    )
  }
  if (isTRUE(run_composite)) {
    scenarios[[length(scenarios) + 1L]] <- make_small_circle_symmetric_mixture2_composite_calibration_scenario(
      kappa = kappa,
      nu = nu
    )
  }
  scenarios
}

run_small_circle_symmetric_mixture2_calibration_sequential <- function(
  output_root = file.path(
    "output",
    "bootstrap_calibration",
    "small_circle_symmetric_mixture2_simple_composite"
  ),
  n_values = c(50L, 100L, 200L),
  M_outer = 50L,
  B = 50L,
  n_cores_outer = 8L,
  statistics = c("ks", "cvm"),
  alphas = c(0.01, 0.05, 0.10),
  seed = 20260604L,
  show_progress = TRUE,
  verbose = TRUE,
  run_simple = TRUE,
  run_composite = TRUE,
  kappa = 12,
  nu = 0.6
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  scenarios <- default_scsm2_selected_scenarios(
    run_simple = isTRUE(run_simple),
    run_composite = isTRUE(run_composite),
    kappa = as.numeric(kappa),
    nu = as.numeric(nu)
  )
  run_rows <- list()

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    scenario_seed <- as.integer(seed) + i
    scenario_dir <- file.path(output_root, sprintf("%02d_%s", i, scenario$id))

    run_info <- scsm2_run_single_scenario_with_checkpoint(
      scenario = scenario,
      output_dir = scenario_dir,
      n_values = as.integer(n_values),
      M_outer = as.integer(M_outer),
      B = as.integer(B),
      statistics = statistics,
      alphas = as.numeric(alphas),
      n_cores_outer = as.integer(n_cores_outer),
      seed = scenario_seed,
      show_progress = isTRUE(show_progress),
      verbose = isTRUE(verbose)
    )

    run_rows[[length(run_rows) + 1L]] <- data.frame(
      order_id = i,
      scenario_id = run_info$scenario_id,
      output_dir = run_info$output_dir,
      skipped = run_info$skipped,
      raw_csv = run_info$raw_csv,
      summary_csv = run_info$summary_csv,
      stringsAsFactors = FALSE
    )
  }

  manifest_df <- do.call(rbind, run_rows)
  consolidated <- scsm2_consolidate_sequential_outputs(
    run_rows = manifest_df,
    consolidated_dir = file.path(output_root, "_consolidated")
  )

  result <- list(
    output_root = output_root,
    manifest = manifest_df,
    consolidated = consolidated,
    config = list(
      n_values = as.integer(n_values),
      M_outer = as.integer(M_outer),
      B = as.integer(B),
      n_cores_outer = as.integer(n_cores_outer),
      statistics = statistics,
      alphas = as.numeric(alphas),
      seed = as.integer(seed),
      run_simple = isTRUE(run_simple),
      run_composite = isTRUE(run_composite),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu)
    )
  )

  saveRDS(result, file = file.path(output_root, "run_result.rds"))
  result
}

if (sys.nframe() == 0L) {
  print(run_small_circle_symmetric_mixture2_calibration_sequential())
}
