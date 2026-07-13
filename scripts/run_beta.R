source(file.path("bootstrap", "calibration_study.R"))
source(file.path("scripts", "run_comets_rotational_mixtures_short_long.R"))

timestamp_tag_beta_ <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

expected_raw_rows_beta_ <- function(M_outer,
                                             n_values,
                                             statistics) {
  as.integer(M_outer) * length(as.integer(n_values)) * length(statistics)
}

expected_summary_rows_beta_ <- function(n_values,
                                                 statistics,
                                                 alphas) {
  length(as.integer(n_values)) * length(statistics) * length(as.numeric(alphas))
}

beta__output_complete <- function(output_dir,
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

estimate_beta_composite_hours <- function(n_values,
                                          B,
                                          M_outer = 500L,
                                          n_cores_outer = 12L) {
  # Benchmark anchor previously measured for the current beta composite code path:
  # n = 200, B = 500, KS + CvM jointly, about 6.59 minutes per outer replicate.
  per_outer_minutes_n200_B500 <- 6.59
  scale_B <- as.numeric(B) / 500
  per_outer_minutes <- per_outer_minutes_n200_B500 * (as.numeric(n_values) / 200) * scale_B
  batches <- ceiling(as.integer(M_outer) / as.integer(n_cores_outer))
  total_minutes <- batches * sum(per_outer_minutes)
  total_minutes / 60
}

choose_beta_composite_plan <- function(max_hours = 10,
                                       M_outer = 500L,
                                       n_cores_outer = 12L) {
  candidate_full <- list(
    n_values = c(50L, 100L, 200L),
    B = 5000L
  )
  candidate_short <- list(
    n_values = c(50L, 100L),
    B = 5000L
  )

  candidate_full$estimated_hours <- estimate_beta_composite_hours(
    n_values = candidate_full$n_values,
    B = candidate_full$B,
    M_outer = M_outer,
    n_cores_outer = n_cores_outer
  )
  if (candidate_full$estimated_hours <= max_hours) {
    candidate_full$decision <- "full_n_with_B500"
    return(candidate_full)
  }

  candidate_short$estimated_hours <- estimate_beta_composite_hours(
    n_values = candidate_short$n_values,
    B = candidate_short$B,
    M_outer = M_outer,
    n_cores_outer = n_cores_outer
  )
  if (candidate_short$estimated_hours <= max_hours) {
    candidate_short$decision <- "n50_n100_with_B500"
    return(candidate_short)
  }

  max_B <- floor(500 * (max_hours / candidate_short$estimated_hours))
  max_B <- max(50L, min(500L, as.integer(max_B)))

  list(
    n_values = c(50L, 100L),
    B = max_B,
    estimated_hours = estimate_beta_composite_hours(
      n_values = c(50L, 100L),
      B = max_B,
      M_outer = M_outer,
      n_cores_outer = n_cores_outer
    ),
    decision = "n50_n100_with_reduced_B"
  )
}

run_beta_composite_stage <- function(output_dir,
                                     n_values,
                                     B,
                                     M_outer = 500L,
                                     n_cores_outer = 12L,
                                     statistics = c("ks", "cvm"),
                                     alphas = c(0.01, 0.05, 0.10),
                                     seed = 20260602L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  expected_raw <- expected_raw_rows_beta(
    M_outer = M_outer,
    n_values = n_values,
    statistics = statistics
  )
  expected_summary <- expected_summary_rows_beta(
    n_values = n_values,
    statistics = statistics,
    alphas = alphas
  )

  if (beta__output_complete(
    output_dir = output_dir,
    expected_raw_rows = expected_raw,
    expected_summary_rows = expected_summary
  )) {
    message("[SKIP] Beta composite calibration already complete at ", output_dir)
    return(invisible(list(output_dir = output_dir, skipped = TRUE)))
  }

  scenario <- default_beta_mixture2_composite_calibration_scenarios()[[1L]]
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
    show_progress = FALSE,
    verbose = TRUE
  )

  saveRDS(result, file = file.path(output_dir, "run_result.rds"))
  invisible(result)
}

run_beta_comets_stage <- function(output_dir,
                                  B = 5000L,
                                  statistics = c("ks", "cvm"),
                                  n_cores = 12L,
                                  seed = 20260602L) {
  summary_path <- file.path(output_dir, "comets_rotational_mixtures_summary.csv")
  if (file.exists(summary_path)) {
    message("[SKIP] Beta comet analysis already complete at ", output_dir)
    return(invisible(list(output_root = output_dir, skipped = TRUE)))
  }

  run_comets_rotational_mixtures_short_long(
    output_root = output_dir,
    datasets = c("short", "long"),
    models = c("beta_mixture2"),
    B = as.integer(B),
    statistics = statistics,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    M_value = 60L,
    ks_t_points = 200L,
    distance_type = "geodesic"
  )
}

run_beta <- function(
  output_root = file.path("output", "beta"),
  max_beta_hours = 10,
  M_outer = 500L,
  n_cores_outer = 12L,
  comet_B = 1000L,
  seed = 20260602L
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  plan <- choose_beta_composite_plan(
    max_hours = max_beta_hours,
    M_outer = M_outer,
    n_cores_outer = n_cores_outer
  )
  plan_df <- data.frame(
    decision = plan$decision,
    n_values = paste(plan$n_values, collapse = ","),
    B = as.integer(plan$B),
    estimated_hours = as.numeric(plan$estimated_hours),
    M_outer = as.integer(M_outer),
    n_cores_outer = as.integer(n_cores_outer),
    comet_B = as.integer(comet_B),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(plan_df, file.path(output_root, "plan.csv"), row.names = FALSE)

  stage_rows <- list()

  beta_output_dir <- file.path(output_root, "01_beta_composite_calibration")
  beta_status <- tryCatch({
    t0 <- Sys.time()
    run_beta_composite_stage(
      output_dir = beta_output_dir,
      n_values = plan$n_values,
      B = plan$B,
      M_outer = M_outer,
      n_cores_outer = n_cores_outer,
      seed = seed + 1L
    )
    data.frame(
      stage = "beta_composite_calibration",
      status = "ok",
      started_at = format(t0, "%Y-%m-%d %H:%M:%S %Z"),
      ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      stage = "beta_composite_calibration",
      status = paste("error:", conditionMessage(e)),
      started_at = NA_character_,
      ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
  })
  stage_rows[[length(stage_rows) + 1L]] <- beta_status
  utils::write.csv(do.call(rbind, stage_rows), file.path(output_root, "stage_status.csv"), row.names = FALSE)

  comets_output_dir <- file.path(output_root, "02_comets_short_long_beta")
  comets_status <- tryCatch({
    t0 <- Sys.time()
    run_beta_comets_stage(
      output_dir = comets_output_dir,
      B = comet_B,
      statistics = c("ks", "cvm"),
      n_cores = n_cores_outer,
      seed = seed + 100L
    )
    data.frame(
      stage = "comets_short_long_beta",
      status = "ok",
      started_at = format(t0, "%Y-%m-%d %H:%M:%S %Z"),
      ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      stage = "comets_short_long_beta",
      status = paste("error:", conditionMessage(e)),
      started_at = NA_character_,
      ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
  })
  stage_rows[[length(stage_rows) + 1L]] <- comets_status
  utils::write.csv(do.call(rbind, stage_rows), file.path(output_root, "stage_status.csv"), row.names = FALSE)

  invisible(list(output_root = output_root, plan = plan_df, stages = do.call(rbind, stage_rows)))
}

if (sys.nframe() == 0L) {
  run_beta()
}
