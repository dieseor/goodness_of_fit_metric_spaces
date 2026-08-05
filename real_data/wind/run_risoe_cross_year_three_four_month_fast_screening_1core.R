#!/usr/bin/env Rscript

# Missing cyclic contiguous windows only: those crossing December--January.
# The outer loop is serial and every bootstrap uses exactly one core.
# No observations are imputed.  The data rule is inherited from the main
# runner: 77 m through 2004-12-31 and 125 m strictly before 2004-11-12.
#
# Example:
#   B=5000 Rscript real_data/wind/run_risoe_cross_year_three_four_month_fast_screening_1core.R

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

source(file.path(repo_root, "real_data", "wind", "run_risoe_contiguous_three_four_month_fast_screening_2004.R"))

make_cross_year_configs <- function() {
  missing_windows <- list(
    list(months = c(11L, 12L, 1L), window_length = 3L)
    # list(months = c(12L, 1L, 2L), window_length = 3L),
    # list(months = c(10L, 11L, 12L, 1L), window_length = 4L),
    # list(months = c(11L, 12L, 1L, 2L), window_length = 4L),
    # list(months = c(12L, 1L, 2L, 3L), window_length = 4L)
  )

  height_specs <- list(
    # list(height_m = 77L, speed_col = "ws77", direction_col = "wd77", cutoff = as.Date("2005-01-01")),
    list(
      height_m = 125L,
      speed_col = "ws125",
      direction_col = "wd125",
      cutoff = as.Date("2004-11-12")
    )
  )

  configs <- list()
  index <- 1L

  for (height_spec in height_specs) {
    for (window in missing_windows) {
      months_id <- paste(month_labels[window$months], collapse = "_")

      for (pattern in names(day_patterns)) {
        configs[[index]] <- c(
          height_spec,
          list(
            months = window$months,
            months_id = months_id,
            window_length = window$window_length,
            first_month = window$months[[1L]],
            last_month = window$months[[length(window$months)]],
            pattern = pattern,
            day_pattern = day_patterns[[pattern]],
            case_id = sprintf(
              "risoe_%dm_%s_%s_1996_2004",
              height_spec$height_m,
              months_id,
              pattern
            )
          )
        )

        index <- index + 1L
      }
    }
  }

  configs
}

run_cross_year_screening <- function(
    input_nc = file.path(
      repo_root,
      "real_data",
      "wind",
      "risoe_m_all.nc"
    ),
    B = read_positive_integer("B", 5000L),
    base_seed = read_positive_integer("BASE_SEED", 2026080400L),
    output_dir = Sys.getenv(
      "OUTPUT_DIR",
      unset = file.path(
        repo_root,
        "real_data",
        "wind",
        "month_diagnostics",
        sprintf(
          "screening_fast_cross_year_125m_nov_dec_jan_1996_2004_b%d_1core",
          B
        )
      )
    )) {

  n_cores <- 1L
  configs <- make_cross_year_configs()

  stopifnot(
    length(configs) == 4L,
    all(vapply(
      configs,
      function(config) identical(config$height_m, 125L),
      logical(1L)
    )),
    all(vapply(
      configs,
      function(config) identical(config$months, c(11L, 12L, 1L)),
      logical(1L)
    ))
  )

  if (!file.exists(input_nc)) {
    stop(
      sprintf("Input NetCDF file not found: %s", input_nc),
      call. = FALSE
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  checkpoint_path <- file.path(
    output_dir,
    sprintf(
      "risoe_fast_cross_year_125m_nov_dec_jan_1996_2004_b%d_1core_checkpoint.csv",
      B
    )
  )

  summary_path <- file.path(
    output_dir,
    sprintf(
      "risoe_fast_cross_year_125m_nov_dec_jan_1996_2004_b%d_1core_summary.csv",
      B
    )
  )
  
  results <- if (file.exists(checkpoint_path)) utils::read.csv(checkpoint_path, stringsAsFactors = FALSE, check.names = FALSE) else NULL
  selected_df <- select_noon_all_months(load_risoe_concurrent(input_nc, fixed_tz = "UTC"), fixed_tz = "UTC")

  cat(sprintf("Cross-year fast screening: %d datasets; B=%d; exactly one core per bootstrap.\n\n", length(configs), B))
  for (i in seq_along(configs)) {
    config <- configs[[i]]
    if (completed_case(results, config$case_id, B)) {
      cat(sprintf("[%d/%d] Resuming: %s already complete.\n", i, length(configs), config$case_id))
      next
    }
    seed <- as.integer(base_seed + i)
    cat(sprintf("[%d/%d] %dm | months=%s | %s | B=%d | cores=1 | seed=%d\n", i, length(configs), config$height_m, paste(config$months, collapse = ","), config$pattern, B, seed))
    flush.console()
    case_result <- tryCatch(
      run_case(selected_df, config, B, n_cores, seed),
      error = function(e) make_error_case_result(config, B, n_cores, seed, conditionMessage(e))
    )
    if (!is.null(results) && nrow(results) > 0L) {
      results <- results[!(results$case_id == config$case_id & results$B == B), , drop = FALSE]
      if (nrow(results) == 0L) results <- NULL
    }
    results <- if (is.null(results)) case_result else rbind(results, case_result)
    utils::write.csv(results, checkpoint_path, row.names = FALSE)
    utils::write.csv(make_summary(results), summary_path, row.names = FALSE)
    if (!all(case_result$status == "ok")) {
      stop(sprintf("Case %s failed and was checkpointed: %s", config$case_id, case_result$error_message[[1L]]), call. = FALSE)
    }
    cat(sprintf("          n=%d (2004: %d); p_KS=%.6f; p_CvM=%.6f\n", case_result$n_valid[[1L]], case_result$n_valid_2004[[1L]], case_result$p_value[case_result$statistic == "KS"], case_result$p_value[case_result$statistic == "CvM"]))
    flush.console()
  }
  cat(sprintf("Completed. Checkpoint: %s\nSummary: %s\n", normalizePath(checkpoint_path, winslash = "/", mustWork = TRUE), normalizePath(summary_path, winslash = "/", mustWork = TRUE)))
}

if (sys.nframe() == 0L) run_cross_year_screening()
