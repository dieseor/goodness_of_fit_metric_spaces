resolve_jp_cardioid_calibration_path <- function(...) {
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

calibration_study_script_jp_cardioid <- resolve_jp_cardioid_calibration_path(
  "bootstrap",
  "calibration_study.R"
)
source(calibration_study_script_jp_cardioid)

expected_raw_rows_calibration <- function(M_outer,
                                          n_values,
                                          statistics) {
  as.integer(M_outer) * length(as.integer(n_values)) * length(statistics)
}

expected_summary_rows_calibration <- function(n_values,
                                              statistics,
                                              alphas) {
  length(as.integer(n_values)) * length(statistics) * length(alphas)
}

is_scenario_output_complete <- function(output_dir,
                                        expected_raw_rows,
                                        expected_summary_rows) {
  raw_csv <- file.path(output_dir, "bootstrap_calibration_raw.csv")
  summary_csv <- file.path(output_dir, "bootstrap_calibration_summary.csv")

  if (!file.exists(raw_csv) || !file.exists(summary_csv)) {
    return(FALSE)
  }

  raw_ok <- FALSE
  summary_ok <- FALSE

  raw_try <- try(utils::read.csv(raw_csv, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(raw_try, "try-error") && is.data.frame(raw_try)) {
    raw_ok <- nrow(raw_try) >= expected_raw_rows
  }

  summary_try <- try(utils::read.csv(summary_csv, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(summary_try, "try-error") && is.data.frame(summary_try)) {
    summary_ok <- nrow(summary_try) >= expected_summary_rows
  }

  raw_ok && summary_ok
}

build_jp_nm_boot_local_control <- function(jp_mle_maxit = 500L,
                                           jp_mle_reltol = 1e-10,
                                           jp_mle_max_abs_kappa_psi = 6) {
  list(
    jp_mle_method = "Nelder-Mead",
    jp_mle_maxit = as.integer(jp_mle_maxit),
    jp_mle_reltol = as.numeric(jp_mle_reltol),
    jp_mle_max_abs_kappa_psi = as.numeric(jp_mle_max_abs_kappa_psi),
    jp_profile_method = "tabulated",
    jp_profile_n_u = 1025L,
    jp_profile_n_delta = 257L,
    jp_vmf_switch_abs_kappa_psi = 1e-3
  )
}

build_cardioid_composite_control <- function(cardioid_maxit = 500L,
                                             cardioid_reltol = 1e-10) {
  # Cardioid MLE uses L-BFGS-B in fit_cardioid_theta; this controls that optimizer.
  list(
    cardioid_optim_control = list(
      maxit = as.integer(cardioid_maxit),
      reltol = as.numeric(cardioid_reltol)
    )
  )
}

run_single_scenario_with_checkpoint <- function(scenario,
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

  expected_raw <- expected_raw_rows_calibration(
    M_outer = M_outer,
    n_values = n_values,
    statistics = statistics
  )
  expected_summary <- expected_summary_rows_calibration(
    n_values = n_values,
    statistics = statistics,
    alphas = alphas
  )

  if (is_scenario_output_complete(
    output_dir = output_dir,
    expected_raw_rows = expected_raw,
    expected_summary_rows = expected_summary
  )) {
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

consolidate_sequential_outputs <- function(run_rows,
                                           consolidated_dir) {
  dir.create(consolidated_dir, recursive = TRUE, showWarnings = FALSE)

  summary_dfs <- lapply(run_rows$summary_csv, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$source_summary_csv <- path
    x
  })

  raw_dfs <- lapply(run_rows$raw_csv, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$source_raw_csv <- path
    x
  })

  summary_all <- do.call(rbind, summary_dfs)
  raw_all <- do.call(rbind, raw_dfs)

  summary_out <- file.path(consolidated_dir, "bootstrap_calibration_summary_all_scenarios.csv")
  raw_out <- file.path(consolidated_dir, "bootstrap_calibration_raw_all_scenarios.csv")
  manifest_out <- file.path(consolidated_dir, "sequential_manifest.csv")

  utils::write.csv(summary_all, summary_out, row.names = FALSE)
  utils::write.csv(raw_all, raw_out, row.names = FALSE)
  utils::write.csv(run_rows, manifest_out, row.names = FALSE)

  list(
    summary_csv = summary_out,
    raw_csv = raw_out,
    manifest_csv = manifest_out
  )
}

run_jp_cardioid_composite_calibration_sequential <- function(
  output_root = file.path("output", "calibration", "bootstrap", "jp_cardioid_composite_M500_B500_NM_bootlocal"),
  n_values = c(50L, 100L, 200L),
  M_outer = 500L,
  B = 500L,
  n_cores_outer = 12L,
  statistics = c("ks", "cvm"),
  alphas = c(0.01, 0.05, 0.10),
  seed = 20260530L,
  jp_mle_maxit = 500L,
  jp_mle_reltol = 1e-10,
  jp_mle_max_abs_kappa_psi = 6,
  cardioid_rho = 0.5,
  cardioid_maxit = 500L,
  cardioid_reltol = 1e-10,
  show_progress = TRUE,
  verbose = TRUE
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  jp_scenarios <- default_jp_composite_calibration_scenarios()
  jp_control <- build_jp_nm_boot_local_control(
    jp_mle_maxit = jp_mle_maxit,
    jp_mle_reltol = jp_mle_reltol,
    jp_mle_max_abs_kappa_psi = jp_mle_max_abs_kappa_psi
  )
  jp_scenarios <- lapply(jp_scenarios, function(scenario) {
    scenario$control <- modifyList(scenario$control %||% list(), jp_control)
    scenario
  })

  cardioid_scenarios <- default_cardioid_composite_calibration_scenarios(rho = cardioid_rho)
  cardioid_control <- build_cardioid_composite_control(
    cardioid_maxit = cardioid_maxit,
    cardioid_reltol = cardioid_reltol
  )
  cardioid_scenarios <- lapply(cardioid_scenarios, function(scenario) {
    scenario$control <- modifyList(scenario$control %||% list(), cardioid_control)
    scenario
  })

  # Execution order is explicit to guarantee at-least-one scenario is saved
  # even if the full sequence is interrupted.
  ordered_scenarios <- c(jp_scenarios, cardioid_scenarios)

  run_rows <- list()
  for (i in seq_along(ordered_scenarios)) {
    scenario <- ordered_scenarios[[i]]
    scenario_seed <- as.integer(seed) + i
    scenario_dir <- file.path(output_root, sprintf("%02d_%s", i, scenario$id))

    run_info <- run_single_scenario_with_checkpoint(
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
  consolidated <- consolidate_sequential_outputs(
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
      jp_control = jp_control,
      cardioid_control = cardioid_control
    )
  )

  saveRDS(result, file = file.path(output_root, "run_result.rds"))
  result
}
