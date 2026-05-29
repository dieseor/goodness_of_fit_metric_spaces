source(file.path("wind", "preprocess_risoe_modern_hvmf.R"))
source(file.path("wind", "run_hvmf_real_data_cvm.R"))

jensen_like_day_patterns <- function() {
  list(
    set12 = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
    set3 = c(1L, 5L, 9L, 13L, 17L, 21L, 25L, 29L)
  )
}

rbind_fill_data_frames <- function(data_frames) {
  if (length(data_frames) == 0L) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(data_frames, names), use.names = FALSE))
  aligned <- lapply(data_frames, function(df) {
    missing_names <- setdiff(all_names, names(df))
    if (length(missing_names) > 0L) {
      for (name in missing_names) {
        df[[name]] <- NA
      }
    }
    df <- df[, all_names, drop = FALSE]
    rownames(df) <- NULL
    df
  })

  do.call(rbind, aligned)
}

count_pattern_days_in_month <- function(year, month, day_pattern, fixed_tz = "UTC") {
  month_start <- as.Date(sprintf("%04d-%02d-01", year, month), tz = fixed_tz)
  next_month <- seq(month_start, length.out = 2L, by = "1 month")[[2L]]
  month_end_day <- as.integer(format(next_month - 1, "%d"))
  sum(as.integer(day_pattern) <= month_end_day)
}

identify_complete_nov_dec_years <- function(df, fixed_tz = "UTC") {
  patterns <- jensen_like_day_patterns()
  candidate_years <- sort(unique(df$year[df$month %in% c(11L, 12L)]))
  complete_years <- integer()

  for (year in candidate_years) {
    is_complete <- TRUE

    for (month in c(11L, 12L)) {
      month_df <- df[df$year == year & df$month == month, , drop = FALSE]
      if (nrow(month_df) == 0L) {
        is_complete <- FALSE
        break
      }

      dates_present <- unique(as.Date(month_df$datetime, tz = fixed_tz))
      days_present <- as.integer(format(dates_present, "%d"))

      for (pattern in patterns) {
        if (!all(pattern %in% days_present)) {
          is_complete <- FALSE
          break
        }
      }

      if (!is_complete) {
        break
      }
    }

    if (is_complete) {
      complete_years <- c(complete_years, year)
    }
  }

  complete_years
}

make_risoe_all_window_configs <- function(complete_years) {
  complete_years <- sort(unique(as.integer(complete_years)))
  configs <- list()
  idx <- 1L

  if (length(complete_years) >= 10L) {
    first_10 <- complete_years[seq_len(10L)]
    configs[[idx]] <- list(
      window_id = sprintf(
        "10y_complete_%d_%d_excl_%s",
        min(first_10),
        max(first_10),
        paste(setdiff(seq(min(first_10), max(first_10)), first_10), collapse = "_")
      ),
      years = first_10,
      n_years_nominal = 10L
    )
    idx <- idx + 1L

    last_10 <- tail(complete_years, 10L)
    if (!identical(last_10, first_10)) {
      configs[[idx]] <- list(
        window_id = sprintf(
          "10y_complete_%d_%d_excl_%s",
          min(last_10),
          max(last_10),
          paste(setdiff(seq(min(last_10), max(last_10)), last_10), collapse = "_")
        ),
        years = last_10,
        n_years_nominal = 10L
      )
      idx <- idx + 1L
    }
  }

  if (length(complete_years) >= 20L) {
    first_20 <- complete_years[seq_len(20L)]
    configs[[idx]] <- list(
      window_id = sprintf(
        "20y_complete_%d_%d_excl_%s",
        min(first_20),
        max(first_20),
        paste(setdiff(seq(min(first_20), max(first_20)), first_20), collapse = "_")
      ),
      years = first_20,
      n_years_nominal = 20L
    )
    idx <- idx + 1L
  }

  if (length(complete_years) >= 30L) {
    first_30 <- complete_years[seq_len(30L)]
    configs[[idx]] <- list(
      window_id = sprintf(
        "30y_complete_%d_%d_excl_%s",
        min(first_30),
        max(first_30),
        paste(setdiff(seq(min(first_30), max(first_30)), first_30), collapse = "_")
      ),
      years = first_30,
      n_years_nominal = 30L
    )
    idx <- idx + 1L
  }

  if (length(complete_years) > 0L) {
    configs[[idx]] <- list(
      window_id = sprintf(
        "max_available_complete_%d_%d_excl_%s",
        min(complete_years),
        max(complete_years),
        paste(setdiff(seq(min(complete_years), max(complete_years)), complete_years), collapse = "_")
      ),
      years = complete_years,
      n_years_nominal = length(complete_years)
    )
  }

  configs
}

filter_window_pattern_rows <- function(selected_df, years, day_pattern) {
  keep <- selected_df$year %in% years &
    selected_df$month %in% c(11L, 12L) &
    selected_df$day %in% as.integer(day_pattern)
  output <- selected_df[keep, , drop = FALSE]
  rownames(output) <- NULL
  output
}

summarize_case_availability <- function(selected_df,
                                        years,
                                        day_pattern,
                                        pattern_name,
                                        speed_col,
                                        direction_col,
                                        height_m,
                                        window_id,
                                        dataset_id,
                                        source_file,
                                        fixed_tz = "UTC") {
  case_rows <- filter_window_pattern_rows(selected_df, years = years, day_pattern = day_pattern)
  summary_rows <- vector("list", length(years) * 2L)
  idx <- 1L

  for (year in years) {
    for (month in c(11L, 12L)) {
      month_rows <- case_rows[case_rows$year == year & case_rows$month == month, , drop = FALSE]
      expected_nominal <- count_pattern_days_in_month(year, month, day_pattern, fixed_tz = fixed_tz)
      dates_available <- length(unique(as.Date(month_rows$datetime, tz = fixed_tz)))
      valid <- is.finite(month_rows[[speed_col]]) &
        is.finite(month_rows[[direction_col]]) &
        month_rows[[speed_col]] > 0

      summary_rows[[idx]] <- data.frame(
        section = "window_year_month",
        source_file = source_file,
        dataset_id = dataset_id,
        window_id = window_id,
        pattern = pattern_name,
        height_m = as.numeric(height_m),
        year = as.integer(year),
        month = as.integer(month),
        expected_nominal = expected_nominal,
        dates_available = dates_available,
        selected_dates = nrow(month_rows),
        valid_after_filter = sum(valid),
        missing_or_removed = nrow(month_rows) - sum(valid),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_df
}

augment_hvmf_case_df <- function(df_case,
                                 dataset_id,
                                 pattern_name,
                                 window_id,
                                 source_file) {
  output <- data.frame(
    dataset_id = rep(dataset_id, nrow(df_case)),
    datetime = df_case$datetime,
    year = df_case$year,
    month = df_case$month,
    day = df_case$day,
    hour = df_case$hour,
    minute = df_case$minute,
    height_m = df_case$height_m,
    pattern = rep(pattern_name, nrow(df_case)),
    window_id = rep(window_id, nrow(df_case)),
    source_file = rep(source_file, nrow(df_case)),
    speed = df_case$speed,
    direction_deg = df_case$direction_deg,
    speed_mean = df_case$speed_mean_height,
    speed_scaled = df_case$speed_scaled,
    angle_rad = df_case$angle_rad,
    x0 = df_case$x0,
    x1 = df_case$x1,
    x2 = df_case$x2,
    minkowski_norm = df_case$minkowski_norm,
    stringsAsFactors = FALSE
  )

  attr(output, "dropped_days") <- attr(df_case, "dropped_days")
  output
}

write_augmented_hvmf_case_csv <- function(df, output_csv, fixed_tz = "UTC") {
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  output_df <- df
  output_df$datetime <- format_datetime_for_csv(output_df$datetime, fixed_tz = fixed_tz)
  write.csv(output_df, file = output_csv, row.names = FALSE)
  output_csv
}

make_risoe_all_case_configs <- function(window_configs,
                                        output_dir = file.path("wind", "jensen_like_all")) {
  patterns <- jensen_like_day_patterns()
  heights <- list(
    list(id = "77m", speed_col = "ws77", direction_col = "wd77", height_m = 77L),
    list(id = "125m", speed_col = "ws125", direction_col = "wd125", height_m = 125L)
  )

  configs <- list()
  idx <- 1L

  for (window in window_configs) {
    for (pattern_name in names(patterns)) {
      for (height in heights) {
        dataset_id <- sprintf(
          "risoe_all_%s_%s_nov_dec_%s",
          height$id,
          window$window_id,
          pattern_name
        )
        configs[[idx]] <- list(
          dataset_id = dataset_id,
          output_csv = file.path(output_dir, paste0(dataset_id, "_hvmf.csv")),
          pattern_name = pattern_name,
          day_pattern = patterns[[pattern_name]],
          speed_col = height$speed_col,
          direction_col = height$direction_col,
          height_m = height$height_m,
          window_id = window$window_id,
          years = window$years,
          n_years_nominal = window$n_years_nominal
        )
        idx <- idx + 1L
      }
    }
  }

  configs
}

compute_variable_missing_summary <- function(df, source_file) {
  vars <- c("ws77", "wd77", "ws125", "wd125")

  do.call(
    rbind,
    lapply(vars, function(var_name) {
      values <- as.numeric(df[[var_name]])
      data.frame(
        section = "variable_missing",
        source_file = source_file,
        variable = var_name,
        finite_count = sum(is.finite(values)),
        missing_count = sum(!is.finite(values)),
        zero_count = sum(values == 0, na.rm = TRUE),
        min_value = min(values, na.rm = TRUE),
        max_value = max(values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
}

compute_file_metadata_summary <- function(df, source_file, complete_years, metadata) {
  data.frame(
    section = "file_metadata",
    source_file = source_file,
    min_datetime = format_datetime_for_csv(min(df$datetime), fixed_tz = metadata$fixed_tz),
    max_datetime = format_datetime_for_csv(max(df$datetime), fixed_tz = metadata$fixed_tz),
    years_available = paste(sort(unique(df$year)), collapse = ","),
    complete_nov_dec_years = paste(complete_years, collapse = ","),
    time_units = metadata$time_units,
    time_calendar = ifelse(is.na(metadata$time_calendar), "NA", metadata$time_calendar),
    fixed_tz = metadata$fixed_tz,
    stringsAsFactors = FALSE
  )
}

compute_case_total_row <- function(case_summary_df,
                                   dataset_id,
                                   window_id,
                                   pattern_name,
                                   height_m,
                                   source_file) {
  data.frame(
    section = "window_total",
    source_file = source_file,
    dataset_id = dataset_id,
    window_id = window_id,
    pattern = pattern_name,
    height_m = as.numeric(height_m),
    year = NA_integer_,
    month = NA_integer_,
    expected_nominal = sum(case_summary_df$expected_nominal),
    dates_available = sum(case_summary_df$dates_available),
    selected_dates = sum(case_summary_df$selected_dates),
    valid_after_filter = sum(case_summary_df$valid_after_filter),
    missing_or_removed = sum(case_summary_df$missing_or_removed),
    stringsAsFactors = FALSE
  )
}

format_years_used <- function(years) {
  paste(as.integer(years), collapse = ",")
}

run_single_risoe_all_simple_case <- function(case_config,
                                             case_summary_df,
                                             processed_csv,
                                             results_dir,
                                             logs_dir,
                                             B = 5000L,
                                             n_cores = 10L,
                                             alpha = 0.05,
                                             seed = 63000L,
                                             tol = 1e-8) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  result_rds <- file.path(results_dir, paste0(case_config$dataset_id, "_cvm_simple_plugin_result.rds"))
  log_file <- file.path(logs_dir, paste0(case_config$dataset_id, "_log.txt"))

  start_time <- Sys.time()
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  warning_messages <- character()
  status <- "ok"
  error_message <- NA_character_
  output_row <- NULL

  log_line(log_con, "Dataset: ", case_config$dataset_id)
  log_line(log_con, "Start timestamp: ", format_timestamp(start_time))
  log_line(log_con, "Source file: wind/risoe_m_all.nc")
  log_line(log_con, "Processed CSV: ", processed_csv)
  log_line(log_con, "Output RDS: ", result_rds)
  log_line(log_con, "Window ID: ", case_config$window_id)
  log_line(log_con, "Years used: ", format_years_used(case_config$years))
  log_line(log_con, "Height: ", case_config$height_m, " m")
  log_line(log_con, "Pattern: ", case_config$pattern_name)
  log_line(log_con, "Months: 11,12")
  log_line(log_con, "Nominal expected n: ", sum(case_summary_df$expected_nominal))
  log_line(log_con, "Final n after filtering: ", sum(case_summary_df$valid_after_filter))
  log_line(log_con, "Counts by year-month:")
  for (line in capture.output(print(case_summary_df, row.names = FALSE))) {
    log_line(log_con, line, .console = FALSE)
  }

  tryCatch(
    {
      prepared <- prepare_existing_h2_csv_dataset(processed_csv, tol = tol)
      raw_df <- read.csv(processed_csv, stringsAsFactors = FALSE)
      speed_mean <- unique(raw_df$speed_mean)

      if (length(speed_mean) != 1L || !is.finite(speed_mean)) {
        stop("Processed CSV must contain a single finite `speed_mean` value.")
      }

      log_line(log_con, "Speed mean used for scaling: ", sprintf("%.10f", speed_mean))
      log_line(log_con, "Scaling computed after Jensen-like filtering: TRUE")
      log_line(log_con, "Max abs(minkowski_norm + 1): ", sprintf("%.10e", prepared$max_minkowski_error))

      if (prepared$n < 50L) {
        status <- "skipped_n_lt_50"
        error_message <- sprintf("n_final = %d < 50; GOF not run.", prepared$n)
        log_line(log_con, "Status: ", status)
        log_line(log_con, "Reason: ", error_message)
      } else {
        fit <- hvmf_mle_h2(prepared$data_matrix)
        theta0 <- list(mu = fit$mu, kappa = fit$kappa)

        log_line(log_con, "MLE initial parameters:")
        log_line(log_con, "  kappa_hat = ", sprintf("%.10f", fit$kappa))
        log_line(log_con, "  mu_hat = (", paste(sprintf("%.10f", fit$mu), collapse = ", "), ")")
        log_line(log_con, "  sinh_chi_hat = ", sprintf("%.10f", fit$sinh_chi))
        log_line(log_con, "  theta_deg_hat = ", sprintf("%.10f", fit$theta_deg))
        log_line(log_con, "Configuration: statistic=CvM, null=simple plug-in, B=", B, ", n_cores=", n_cores, ", alpha=", alpha, ", seed=", seed)

        bootstrap_result <- withCallingHandlers(
          multiplier_bootstrap_hvmf(
            data = prepared$data_matrix,
            null = list(type = "simple", theta = theta0),
            statistics = "cvm",
            B = B,
            alpha = alpha,
            n_cores = n_cores,
            seed = seed,
            keep = list(
              observed_process = FALSE,
              bootstrap_statistics = TRUE
            ),
            unknown_param = "both"
          ),
          warning = function(w) {
            warning_messages <<- c(warning_messages, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )

        saveRDS(
          list(
            case_config = case_config,
            case_summary = case_summary_df,
            fit = fit,
            bootstrap_result = bootstrap_result
          ),
          file = result_rds
        )

        elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        output_row <- data.frame(
          dataset_id = case_config$dataset_id,
          source_file = "wind/risoe_m_all.nc",
          height_m = case_config$height_m,
          window_id = case_config$window_id,
          year_start = min(case_config$years),
          year_end = max(case_config$years),
          n_years_nominal = case_config$n_years_nominal,
          pattern = case_config$pattern_name,
          months = "11,12",
          n_expected_nominal = sum(case_summary_df$expected_nominal),
          n_final = prepared$n,
          n_missing_or_removed = sum(case_summary_df$expected_nominal) - prepared$n,
          speed_mean = speed_mean,
          kappa_hat = fit$kappa,
          mu_hat_x0 = fit$mu[[1L]],
          mu_hat_x1 = fit$mu[[2L]],
          mu_hat_x2 = fit$mu[[3L]],
          sinh_chi_hat = fit$sinh_chi,
          theta_deg_hat = fit$theta_deg,
          statistic_CvM = bootstrap_result$observed$cvm$statistic,
          p_value = bootstrap_result$inference$cvm$p_value,
          B = B,
          n_cores = n_cores,
          elapsed_seconds = elapsed_seconds,
          processed_csv = processed_csv,
          result_rds = result_rds,
          log_file = log_file,
          status = status,
          error_message = NA_character_,
          stringsAsFactors = FALSE
        )

        log_line(log_con, "Observed CvM: ", sprintf("%.10f", output_row$statistic_CvM))
        log_line(log_con, "p-value: ", sprintf("%.10f", output_row$p_value))
        log_line(log_con, "Elapsed seconds: ", sprintf("%.3f", elapsed_seconds))
      }
    },
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      log_line(log_con, "ERROR: ", error_message)
    }
  )

  if (length(warning_messages) > 0L) {
    log_line(log_con, "Warnings:")
    for (warning_message in warning_messages) {
      log_line(log_con, "  - ", warning_message)
    }
  } else {
    log_line(log_con, "Warnings: none")
  }

  end_time <- Sys.time()
  elapsed_seconds_total <- as.numeric(difftime(end_time, start_time, units = "secs"))
  log_line(log_con, "End timestamp: ", format_timestamp(end_time))
  log_line(log_con, "Total elapsed seconds: ", sprintf("%.3f", elapsed_seconds_total))

  if (is.null(output_row)) {
    output_row <- data.frame(
      dataset_id = case_config$dataset_id,
      source_file = "wind/risoe_m_all.nc",
      height_m = case_config$height_m,
      window_id = case_config$window_id,
      year_start = min(case_config$years),
      year_end = max(case_config$years),
      n_years_nominal = case_config$n_years_nominal,
      pattern = case_config$pattern_name,
      months = "11,12",
      n_expected_nominal = sum(case_summary_df$expected_nominal),
      n_final = sum(case_summary_df$valid_after_filter),
      n_missing_or_removed = sum(case_summary_df$expected_nominal) - sum(case_summary_df$valid_after_filter),
      speed_mean = NA_real_,
      kappa_hat = NA_real_,
      mu_hat_x0 = NA_real_,
      mu_hat_x1 = NA_real_,
      mu_hat_x2 = NA_real_,
      sinh_chi_hat = NA_real_,
      theta_deg_hat = NA_real_,
      statistic_CvM = NA_real_,
      p_value = NA_real_,
      B = B,
      n_cores = n_cores,
      elapsed_seconds = elapsed_seconds_total,
      processed_csv = processed_csv,
      result_rds = result_rds,
      log_file = log_file,
      status = status,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }

  output_row
}

run_risoe_jensen_like_all <- function(input_nc = "wind/risoe_m_all.nc",
                                      output_dir = file.path("wind", "jensen_like_all"),
                                      B = 5000L,
                                      n_cores = 10L,
                                      fixed_tz = "UTC",
                                      tie_break = "earliest") {
  metadata <- read_risoe_nc_metadata(input_nc, fixed_tz = fixed_tz)
  all_df <- load_risoe_concurrent(input_nc, fixed_tz = fixed_tz)
  selected_df <- select_noon_nov_dec(all_df, tie_break = tie_break, fixed_tz = fixed_tz)
  complete_years <- identify_complete_nov_dec_years(all_df, fixed_tz = fixed_tz)
  window_configs <- make_risoe_all_window_configs(complete_years)

  if (length(window_configs) == 0L) {
    stop("No complete November-December year windows are available in risoe_m_all.nc.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results_dir <- file.path(output_dir, "results")
  logs_dir <- file.path(output_dir, "logs")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Using source file: ", input_nc, "\n", sep = "")
  cat("Fixed timezone convention: ", fixed_tz, " (timestamps used exactly as stored)\n", sep = "")
  cat("Only months 11 and 12 will be used.\n")
  cat("Day patterns are exactly Jensen set12 and set3.\n")
  cat("Nearest-to-noon selection uses tie_break = ", tie_break, ".\n", sep = "")
  cat("Speed scaling uses the mean speed of the final Jensen-like dataset after filtering.\n\n")

  file_metadata_summary <- compute_file_metadata_summary(all_df, input_nc, complete_years, metadata)
  variable_missing_summary <- compute_variable_missing_summary(all_df, input_nc)
  availability_rows <- list(file_metadata_summary, variable_missing_summary)

  configs <- make_risoe_all_case_configs(window_configs, output_dir = output_dir)
  summary_rows <- vector("list", length(configs))

  for (i in seq_along(configs)) {
    config <- configs[[i]]
    cat("\n=== risoe_m_all Jensen-like case ", i, "/", length(configs), ": ", config$dataset_id, " ===\n", sep = "")
    cat("Years used: ", format_years_used(config$years), "\n", sep = "")

    case_rows <- filter_window_pattern_rows(
      selected_df = selected_df,
      years = config$years,
      day_pattern = config$day_pattern
    )
    case_summary_df <- summarize_case_availability(
      selected_df = selected_df,
      years = config$years,
      day_pattern = config$day_pattern,
      pattern_name = config$pattern_name,
      speed_col = config$speed_col,
      direction_col = config$direction_col,
      height_m = config$height_m,
      window_id = config$window_id,
      dataset_id = config$dataset_id,
      source_file = input_nc,
      fixed_tz = fixed_tz
    )
    availability_rows[[length(availability_rows) + 1L]] <- case_summary_df
    availability_rows[[length(availability_rows) + 1L]] <- compute_case_total_row(
      case_summary_df = case_summary_df,
      dataset_id = config$dataset_id,
      window_id = config$window_id,
      pattern_name = config$pattern_name,
      height_m = config$height_m,
      source_file = input_nc
    )

    df_case <- build_hvmf_wind_set(
      df = case_rows,
      speed_col = config$speed_col,
      direction_col = config$direction_col,
      height_m = config$height_m,
      fixed_tz = fixed_tz
    )
    df_case_augmented <- augment_hvmf_case_df(
      df_case = df_case,
      dataset_id = config$dataset_id,
      pattern_name = config$pattern_name,
      window_id = config$window_id,
      source_file = "risoe_m_all.nc"
    )

    output_csv <- write_augmented_hvmf_case_csv(df_case_augmented, config$output_csv, fixed_tz = fixed_tz)

    cat("  n_expected_nominal = ", sum(case_summary_df$expected_nominal), "\n", sep = "")
    cat("  n_final = ", nrow(df_case_augmented), "\n", sep = "")
    cat("  speed_mean = ", sprintf("%.10f", unique(df_case_augmented$speed_mean)), "\n", sep = "")
    cat("  max |minkowski_norm + 1| = ", sprintf("%.10e", max(abs(df_case_augmented$minkowski_norm + 1))), "\n", sep = "")

    summary_rows[[i]] <- run_single_risoe_all_simple_case(
      case_config = config,
      case_summary_df = case_summary_df,
      processed_csv = output_csv,
      results_dir = results_dir,
      logs_dir = logs_dir,
      B = B,
      n_cores = n_cores,
      seed = 63000L + i
    )
  }

  availability_summary <- rbind_fill_data_frames(availability_rows)
  availability_csv <- file.path(output_dir, "risoe_m_all_availability_summary.csv")
  write.csv(availability_summary, availability_csv, row.names = FALSE)

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(summary_df$status != "ok", -summary_df$p_value, -summary_df$n_final), , drop = FALSE]
  rownames(summary_df) <- NULL
  summary_csv <- file.path(output_dir, "risoe_jensen_like_all_summary.csv")
  write.csv(summary_df, summary_csv, row.names = FALSE)

  stopifnot(file.exists(availability_csv))
  stopifnot(file.exists(summary_csv))
  stopifnot(all(file.exists(summary_df$processed_csv)))
  ok_rows <- summary_df$status == "ok"
  if (any(ok_rows)) {
    stopifnot(all(file.exists(summary_df$result_rds[ok_rows])))
    stopifnot(all(file.exists(summary_df$log_file[ok_rows])))
  }

  invisible(list(
    metadata = metadata,
    complete_years = complete_years,
    window_configs = window_configs,
    availability_summary = availability_summary,
    availability_csv = availability_csv,
    summary = summary_df,
    summary_csv = summary_csv
  ))
}

if (sys.nframe() == 0L) {
  run_risoe_jensen_like_all()
}
