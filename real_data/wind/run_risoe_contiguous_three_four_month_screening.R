#!/usr/bin/env Rscript

# Sequential composite-HvMF screening for contiguous three- and four-month
# windows in the cleaned Risoe wind data.  Each bootstrap uses 10 cores, but
# the windows themselves are deliberately run one after another.
#
# The data construction matches the wind experiment in the manuscript:
#   * closest-to-noon observation on each selected date;
#   * years 1996--2003;
#   * strict four-day patterns start1,...,start4;
#   * the existing height-specific cleaning rule.
#
# By default, the script performs 2 x (10 + 9) x 4 = 152 composite tests:
# both heights, all non-wrapping contiguous windows of length 3 and 4, and
# all four day-pattern starts.  KS and CvM are calculated from the same
# composite multiplier-bootstrap run for each dataset.
#
# Run, for example, from any directory with
#   B=5000 Rscript real_data/wind/run_risoe_contiguous_three_four_month_screening.R
#
# `B` defaults to 5000.  `OUTPUT_DIR` and `BASE_SEED` may be set in the
# environment.  `N_CORES` is intentionally not configurable: every test is
# run with exactly 10 cores and the outer loop is sequential.

script_file_argument <- commandArgs(trailingOnly = FALSE)
script_file_argument <- script_file_argument[grepl("^--file=", script_file_argument)]
if (length(script_file_argument) != 1L) {
  stop("Run this file with Rscript so that its repository location is known.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_file_argument), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

# This runner supplies the data transformation, both composite tests, and the
# published height-specific cleaning dates used in the original wind analysis.
source(file.path(repo_root, "real_data", "wind", "run_risoe_125m_screening_ks_cvm.R"))

read_positive_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(sprintf("Environment variable %s must be a positive integer.", name), call. = FALSE)
  }
  parsed
}

month_abbreviations <- c(
  "jan", "feb", "mar", "apr", "may", "jun",
  "jul", "aug", "sep", "oct", "nov", "dec"
)

strict_step4_patterns <- lapply(seq_len(4L), function(start_day) {
  seq.int(start_day, 30L, by = 4L)
})
names(strict_step4_patterns) <- paste0("start", seq_len(4L))

make_contiguous_month_windows <- function(window_lengths = c(3L, 4L)) {
  window_lengths <- sort(unique(as.integer(window_lengths)))
  if (!identical(window_lengths, c(3L, 4L))) {
    stop("This screening is defined for exactly the three- and four-month windows.", call. = FALSE)
  }

  windows <- list()
  index <- 1L
  for (window_length in window_lengths) {
    for (first_month in seq_len(12L - window_length + 1L)) {
      months <- seq.int(first_month, length.out = window_length)
      windows[[index]] <- list(
        months = months,
        months_id = paste(month_abbreviations[months], collapse = "_"),
        window_length = window_length,
        first_month = first_month,
        last_month = max(months)
      )
      index <- index + 1L
    }
  }
  windows
}

make_screening_configs <- function() {
  height_specs <- list(
    list(height_m = 77L, speed_col = "ws77", direction_col = "wd77"),
    list(height_m = 125L, speed_col = "ws125", direction_col = "wd125")
  )
  month_windows <- make_contiguous_month_windows()

  configs <- list()
  index <- 1L
  for (height_spec in height_specs) {
    for (month_window in month_windows) {
      for (pattern_id in names(strict_step4_patterns)) {
        configs[[index]] <- c(
          height_spec,
          month_window,
          list(
            pattern = pattern_id,
            day_pattern = strict_step4_patterns[[pattern_id]],
            case_id = sprintf(
              "risoe_clean_%sm_%s_%s_1996_2003",
              height_spec$height_m,
              month_window$months_id,
              pattern_id
            )
          )
        )
        index <- index + 1L
      }
    }
  }
  configs
}

completed_case <- function(results, case_id, B) {
  if (is.null(results) || nrow(results) == 0L) {
    return(FALSE)
  }
  rows <- results[
    results$case_id == case_id &
      results$B == B &
      results$status == "ok" &
      results$statistic %in% c("KS", "CvM"),
    ,
    drop = FALSE
  ]
  identical(sort(unique(rows$statistic)), c("CvM", "KS"))
}

write_checkpoint <- function(results, path) {
  utils::write.csv(results, path, row.names = FALSE)
  invisible(path)
}

make_summary <- function(results, alpha = 0.05) {
  completed <- results[
    results$status == "ok" & results$statistic %in% c("KS", "CvM"),
    ,
    drop = FALSE
  ]
  if (nrow(completed) == 0L) {
    return(data.frame())
  }

  case_ids <- unique(completed$case_id)
  summary_rows <- lapply(case_ids, function(case_id) {
    rows <- completed[completed$case_id == case_id, , drop = FALSE]
    ks_row <- rows[rows$statistic == "KS", , drop = FALSE]
    cvm_row <- rows[rows$statistic == "CvM", , drop = FALSE]

    if (nrow(ks_row) != 1L || nrow(cvm_row) != 1L) {
      return(NULL)
    }

    data.frame(
      case_id = case_id,
      height_m = ks_row$height_m[[1L]],
      window_length = ks_row$window_length[[1L]],
      month_start = ks_row$month_start[[1L]],
      month_end = ks_row$month_end[[1L]],
      months = ks_row$months[[1L]],
      pattern = ks_row$pattern[[1L]],
      day_pattern = ks_row$day_pattern[[1L]],
      n = ks_row$n[[1L]],
      B = ks_row$B[[1L]],
      p_value_KS = ks_row$p_value[[1L]],
      p_value_CvM = cvm_row$p_value[[1L]],
      reject_KS_unadjusted_5pct = ks_row$p_value[[1L]] < alpha,
      reject_CvM_unadjusted_5pct = cvm_row$p_value[[1L]] < alpha,
      reject_either_unadjusted_5pct =
        ks_row$p_value[[1L]] < alpha || cvm_row$p_value[[1L]] < alpha,
      stringsAsFactors = FALSE
    )
  })

  summary_rows <- Filter(Negate(is.null), summary_rows)
  output <- do.call(rbind, summary_rows)
  output <- output[order(
    output$height_m,
    output$window_length,
    output$month_start,
    output$pattern
  ), , drop = FALSE]
  rownames(output) <- NULL
  output
}

run_risoe_contiguous_three_four_month_screening <- function(
    input_nc = file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"),
    B = read_positive_integer("B", 5000L),
    output_dir = {
      from_environment <- Sys.getenv("OUTPUT_DIR", unset = "")
      if (nzchar(from_environment)) {
        from_environment
      } else {
        file.path(
          repo_root,
          "real_data", "wind", "month_diagnostics",
          sprintf("screening_contiguous_three_four_month_b%d", B)
        )
      }
    },
    base_seed = read_positive_integer("BASE_SEED", 2026080300L)) {
  # This is an intentional invariant of the requested computational design.
  n_cores <- 10L
  configs <- make_screening_configs()

  if (!file.exists(input_nc)) {
    stop(sprintf("Input NetCDF file not found: %s", input_nc), call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  checkpoint_path <- file.path(
    output_dir,
    sprintf("risoe_contiguous_three_four_month_ks_cvm_b%d_checkpoint.csv", B)
  )
  summary_path <- file.path(
    output_dir,
    sprintf("risoe_contiguous_three_four_month_ks_cvm_b%d_summary.csv", B)
  )

  results <- if (file.exists(checkpoint_path)) {
    utils::read.csv(checkpoint_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }

  all_df <- load_risoe_concurrent(input_nc, fixed_tz = "UTC")
  selected_df <- select_noon_all_months(all_df, fixed_tz = "UTC")

  cat(sprintf(
    "Sequential screening: %d datasets; B=%d; exactly %d cores per bootstrap.\n",
    length(configs), B, n_cores
  ))
  cat("The outer loop is serial: the next month window starts only after the current bootstrap finishes.\n")
  cat("The reported rejection indicators are unadjusted screening decisions; windows overlap.\n\n")

  # Do not replace this loop with lapply/future/mclapply: it deliberately
  # reserves all ten cores for the one composite bootstrap currently running.
  for (index in seq_along(configs)) {
    config <- configs[[index]]
    if (completed_case(results, config$case_id, B)) {
      cat(sprintf("[%d/%d] Resuming: %s already complete.\n", index, length(configs), config$case_id))
      next
    }

    case_seed <- as.integer(base_seed + index)
    cat(sprintf(
      "[%d/%d] %s | months=%s | %s | B=%d | cores=%d | seed=%d\n",
      index,
      length(configs),
      config$height_m,
      paste(config$months, collapse = ","),
      config$pattern,
      B,
      n_cores,
      case_seed
    ))
    flush.console()

    case_result <- run_single_risoe_case(
      config = list(
        window = config$months_id,
        months = config$months,
        pattern = config$pattern,
        day_pattern = config$day_pattern
      ),
      selected_df = selected_df,
      speed_col = config$speed_col,
      direction_col = config$direction_col,
      height_m = config$height_m,
      B = B,
      n_cores = n_cores,
      bootstrap_method = "reestimated",
      fixed_tz = "UTC",
      ks_grid_mode = "sample_points_unique_distances",
      seed = case_seed
    )

    case_result$case_id <- config$case_id
    case_result$height_m <- config$height_m
    case_result$window_length <- config$window_length
    case_result$month_start <- config$first_month
    case_result$month_end <- config$last_month
    case_result$months <- paste(config$months, collapse = ",")
    case_result$day_pattern <- paste(config$day_pattern, collapse = ",")
    case_result$case_seed <- case_seed
    case_result$n_cores <- n_cores

    # A rerun after a partially written result replaces the previous rows for
    # this exact case and B, then updates the durable checkpoint immediately.
    if (!is.null(results) && nrow(results) > 0L) {
      keep <- !(results$case_id == config$case_id & results$B == B)
      results <- results[keep, , drop = FALSE]
    }
    results <- if (is.null(results)) case_result else rbind(results, case_result)
    rownames(results) <- NULL
    write_checkpoint(results, checkpoint_path)

    current_summary <- make_summary(results)
    utils::write.csv(current_summary, summary_path, row.names = FALSE)

    if (all(case_result$status == "ok")) {
      p_ks <- case_result$p_value[case_result$statistic == "KS"]
      p_cvm <- case_result$p_value[case_result$statistic == "CvM"]
      cat(sprintf("          n=%d; p_KS=%.6f; p_CvM=%.6f\n", case_result$n[[1L]], p_ks, p_cvm))
    } else {
      cat("          This case failed; its error message is recorded in the checkpoint CSV.\n")
    }
    flush.console()
  }

  final_summary <- make_summary(results)
  utils::write.csv(final_summary, summary_path, row.names = FALSE)
  cat("\nCheckpoint CSV: ", normalizePath(checkpoint_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("Summary CSV:    ", normalizePath(summary_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat(sprintf(
    "Completed datasets with both statistics: %d of %d.\n",
    nrow(final_summary), length(configs)
  ))

  invisible(list(
    results = results,
    summary = final_summary,
    checkpoint_path = checkpoint_path,
    summary_path = summary_path
  ))
}

if (sys.nframe() == 0L) {
  run_risoe_contiguous_three_four_month_screening()
}
