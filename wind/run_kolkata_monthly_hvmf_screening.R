source(file.path("wind", "preprocess_kolkata_power_hvmf.R"))
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))

format_timestamp <- function(x, tz = "UTC") {
  format(x, "%Y-%m-%d %H:%M:%S", tz = tz)
}

log_line <- function(con, ..., .console = FALSE) {
  text <- paste0(..., collapse = "")
  writeLines(text, con = con)
  if (.console) {
    message(text)
  }
}

run_composite_benchmark_case <- function(df,
                                         dataset_id,
                                         output_dir,
                                         B = 500L,
                                         n_cores = 12L,
                                         profile_method = "tabulated") {
  results_dir <- file.path(output_dir, "results")
  logs_dir <- file.path(output_dir, "logs")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  csv_path <- file.path(output_dir, paste0(dataset_id, "_hvmf.csv"))
  write.csv(df, csv_path, row.names = FALSE)

  X <- as.matrix(df[, c("x0", "x1", "x2")])
  fit <- hvmf_mle_h2(X)
  start_time <- Sys.time()
  log_file <- file.path(logs_dir, sprintf("%s_composite_B%d_log.txt", dataset_id, as.integer(B)))
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  warning_messages <- character()
  log_line(log_con, "Dataset: ", dataset_id)
  log_line(log_con, "Start timestamp: ", format_timestamp(start_time))
  log_line(log_con, "Height: ", unique(df$height_m), " m")
  log_line(log_con, "Months id: ", unique(df$months_id))
  log_line(log_con, "Years used: ", paste(sort(unique(df$year)), collapse = ","))
  log_line(log_con, "n_final: ", nrow(df))
  log_line(log_con, "Speed mean used for scaling: ", sprintf("%.10f", unique(df$speed_mean)))
  log_line(log_con, "MLE kappa_hat = ", sprintf("%.10f", fit$kappa))
  log_line(log_con, "MLE mu_hat = (", paste(sprintf("%.10f", fit$mu), collapse = ", "), ")")
  log_line(log_con, "MLE theta_deg_hat = ", sprintf("%.10f", fit$theta_deg))
  log_line(log_con, "Configuration: composite CvM, B=", B, ", n_cores=", n_cores, ", hvmf_profile_method=", profile_method)

  result <- withCallingHandlers(
    multiplier_bootstrap_hvmf(
      data = X,
      null = list(type = "composite"),
      statistics = "cvm",
      B = B,
      alpha = 0.05,
      n_cores = n_cores,
      seed = 105000L + nrow(df) + unique(df$height_m),
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = list(hvmf_profile_method = profile_method),
      unknown_param = "both"
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  result_rds <- file.path(results_dir, sprintf("%s_cvm_composite_B%d_result.rds", dataset_id, as.integer(B)))
  saveRDS(
    list(
      dataset_id = dataset_id,
      fit = fit,
      result = result
    ),
    file = result_rds
  )

  log_line(log_con, "Observed CvM = ", sprintf("%.10f", result$observed$cvm$statistic))
  log_line(log_con, "p-value = ", sprintf("%.10f", result$inference$cvm$p_value))
  log_line(log_con, "Elapsed seconds = ", sprintf("%.3f", elapsed_seconds))
  if (length(warning_messages) == 0L) {
    log_line(log_con, "Warnings: none")
  } else {
    log_line(log_con, "Warnings:")
    for (msg in warning_messages) {
      log_line(log_con, "  - ", msg)
    }
  }
  log_line(log_con, "End timestamp: ", format_timestamp(Sys.time()))

  data.frame(
    dataset_id = dataset_id,
    height_m = unique(df$height_m),
    months_id = unique(df$months_id),
    n_final = nrow(df),
    year_range_used = paste(range(df$year), collapse = "-"),
    speed_mean = unique(df$speed_mean),
    kappa_hat = fit$kappa,
    theta_deg_hat = fit$theta_deg,
    statistic_CvM_composite = result$observed$cvm$statistic,
    p_value_composite = result$inference$cvm$p_value,
    B = B,
    n_cores = n_cores,
    output_plot = NA_character_,
    processed_csv = csv_path,
    result_rds = result_rds,
    log_file = log_file,
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

build_kolkata_month_case <- function(all_daily_df,
                                     month,
                                     start_day = 1L,
                                     step = 4L) {
  month <- as.integer(month)
  month_df <- filter_power_month(all_daily_df, month = month)
  month_diag <- diagnose_power_daily_data(month_df)

  clean_df <- month_df[
    is.finite(month_df$WD10M) &
      is.finite(month_df$WS10M) &
      month_df$WS10M > 0,
    ,
    drop = FALSE
  ]
  rownames(clean_df) <- NULL

  thinned_df <- thin_kolkata_by_day_pattern(
    df = clean_df,
    start_day = start_day,
    step = step,
    max_day = 31L
  )

  if (nrow(thinned_df) == 0L) {
    stop(sprintf("No rows remain for month %02d after day-pattern thinning.", month))
  }

  wind_df_thinned <- data.frame(
    datetime = as.POSIXct(thinned_df$date, tz = "UTC"),
    ws10m = thinned_df$WS10M,
    wd10m = thinned_df$WD10M,
    stringsAsFactors = FALSE
  )

  embedded_df <- build_hvmf_wind_set(
    df = wind_df_thinned,
    speed_col = "ws10m",
    direction_col = "wd10m",
    height_m = 10L,
    fixed_tz = "UTC"
  )

  embedded_df$months_id <- tolower(month.abb[month])
  embedded_df$speed_mean <- embedded_df$speed_mean_height
  embedded_df$dataset_source <- "NASA POWER CERES/MERRA2 daily point"
  embedded_df$start_day <- as.integer(start_day)
  embedded_df$step <- as.integer(step)

  stopifnot(all(embedded_df$day %in% make_day_pattern(start_day = start_day, step = step, max_day = 31L)))

  list(
    month_df = month_df,
    clean_df = clean_df,
    thinned_df = thinned_df,
    embedded_df = embedded_df,
    month_diag = month_diag
  )
}

failed_month_row <- function(month,
                             start_year,
                             end_year,
                             start_day,
                             step,
                             error_message,
                             n_raw_month = NA_integer_,
                             n_clean_month = NA_integer_,
                             n_thinned_month = NA_integer_) {
  data.frame(
    dataset_id = sprintf(
      "kolkata_10m_%s_%d_%d_start%d_step%d",
      tolower(month.abb[as.integer(month)]),
      as.integer(start_year),
      as.integer(end_year),
      as.integer(start_day),
      as.integer(step)
    ),
    month = as.integer(month),
    month_label = month.abb[as.integer(month)],
    n_raw_month = as.integer(n_raw_month),
    n_clean_month = as.integer(n_clean_month),
    n_thinned_month = as.integer(n_thinned_month),
    n_final = as.integer(n_thinned_month),
    kappa_hat = NA_real_,
    theta_deg_hat = NA_real_,
    statistic_CvM_composite = NA_real_,
    p_value_composite = NA_real_,
    elapsed_seconds = NA_real_,
    processed_csv = NA_character_,
    result_rds = NA_character_,
    log_file = NA_character_,
    status = "failed",
    error_message = as.character(error_message),
    stringsAsFactors = FALSE
  )
}

run_kolkata_monthly_hvmf_screening <- function(start_year = 1982L,
                                               end_year = 2022L,
                                               months = 1:12,
                                               start_day = 1L,
                                               step = 4L,
                                               B = 500L,
                                               n_cores = 10L,
                                               profile_method = "tabulated",
                                               output_dir = file.path("wind", "kolkata_monthly_screening")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  processed_dir <- file.path(output_dir, "processed")
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  raw_download <- load_kolkata_power_daily(
    start_year = start_year,
    end_year = end_year,
    latitude = 22.57,
    longitude = 88.36
  )
  all_daily_df <- raw_download$data

  all_daily_csv <- file.path(
    processed_dir,
    sprintf("kolkata_daily_%d_%d_all_months.csv", as.integer(start_year), as.integer(end_year))
  )
  utils::write.csv(all_daily_df, all_daily_csv, row.names = FALSE)

  qa_monthly <- monthly_power_quality_summary(
    df = all_daily_df,
    start_day = start_day,
    step = step,
    max_day = 31L
  )
  qa_monthly_csv <- file.path(
    output_dir,
    sprintf(
      "kolkata_monthly_quality_start%d_step%d_%d_%d.csv",
      as.integer(start_day),
      as.integer(step),
      as.integer(start_year),
      as.integer(end_year)
    )
  )
  utils::write.csv(qa_monthly, qa_monthly_csv, row.names = FALSE)

  summary_rows <- list()
  idx <- 1L

  for (month_value in as.integer(months)) {
    dataset_id <- sprintf(
      "kolkata_10m_%s_%d_%d_start%d_step%d",
      tolower(month.abb[month_value]),
      as.integer(start_year),
      as.integer(end_year),
      as.integer(start_day),
      as.integer(step)
    )
    cat(sprintf("Running Kolkata monthly case %s\n", dataset_id))

    month_result <- tryCatch(
      {
        case <- build_kolkata_month_case(
          all_daily_df = all_daily_df,
          month = month_value,
          start_day = start_day,
          step = step
        )

        raw_month_csv <- file.path(
          processed_dir,
          sprintf("%s_raw_daily.csv", dataset_id)
        )
        thinned_month_csv <- file.path(
          processed_dir,
          sprintf("%s_thinned.csv", dataset_id)
        )
        utils::write.csv(case$clean_df, raw_month_csv, row.names = FALSE)
        utils::write.csv(case$thinned_df, thinned_month_csv, row.names = FALSE)

        benchmark_row <- run_composite_benchmark_case(
          df = case$embedded_df,
          dataset_id = dataset_id,
          output_dir = output_dir,
          B = B,
          n_cores = n_cores,
          profile_method = profile_method
        )

        benchmark_row$month <- as.integer(month_value)
        benchmark_row$month_label <- month.abb[month_value]
        benchmark_row$n_raw_month <- nrow(case$month_df)
        benchmark_row$n_clean_month <- nrow(case$clean_df)
        benchmark_row$n_thinned_month <- nrow(case$thinned_df)
        benchmark_row$raw_month_csv <- raw_month_csv
        benchmark_row$thinned_month_csv <- thinned_month_csv
        benchmark_row$status <- "ok"
        benchmark_row$error_message <- NA_character_
        benchmark_row
      },
      error = function(e) {
        month_df <- filter_power_month(all_daily_df, month = month_value)
        clean_df <- month_df[
          is.finite(month_df$WD10M) &
            is.finite(month_df$WS10M) &
            month_df$WS10M > 0,
          ,
          drop = FALSE
        ]
        thinned_df <- tryCatch(
          thin_kolkata_by_day_pattern(
            df = clean_df,
            start_day = start_day,
            step = step,
            max_day = 31L
          ),
          error = function(...) data.frame()
        )

        failed_month_row(
          month = month_value,
          start_year = start_year,
          end_year = end_year,
          start_day = start_day,
          step = step,
          error_message = conditionMessage(e),
          n_raw_month = nrow(month_df),
          n_clean_month = nrow(clean_df),
          n_thinned_month = nrow(thinned_df)
        )
      }
    )

    summary_rows[[idx]] <- month_result
    idx <- idx + 1L
  }

  summary_df <- do.call(rbind, summary_rows)
  if ("p_value_composite" %in% names(summary_df)) {
    summary_df <- summary_df[order(summary_df$p_value_composite, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  rownames(summary_df) <- NULL

  summary_csv <- file.path(
    output_dir,
    sprintf(
      "kolkata_monthly_hvmf_screening_start%d_step%d_%d_%d.csv",
      as.integer(start_day),
      as.integer(step),
      as.integer(start_year),
      as.integer(end_year)
    )
  )
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  list(
    summary = summary_df,
    summary_csv = summary_csv,
    qa_monthly = qa_monthly,
    qa_monthly_csv = qa_monthly_csv,
    all_daily_csv = all_daily_csv,
    source_url = raw_download$url
  )
}

if (sys.nframe() == 0L) {
  result <- run_kolkata_monthly_hvmf_screening(
    start_year = 1982L,
    end_year = 2022L,
    months = 1:12,
    start_day = 1L,
    step = 4L,
    B = 200L,
    n_cores = 10L,
    profile_method = "tabulated"
  )

  cat("source_url=", result$source_url, "\n", sep = "")
  cat("all_daily_csv=", result$all_daily_csv, "\n", sep = "")
  cat("qa_monthly_csv=", result$qa_monthly_csv, "\n", sep = "")
  cat("summary_csv=", result$summary_csv, "\n", sep = "")
}
