resolve_hvmf_realdata_screening_runner_path <- function(...) {
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

  file.path(...)
}

utils_path_hvmf_realdata_screening <- resolve_hvmf_realdata_screening_runner_path(
  "hvmf_realdata_screening",
  "utils_hvmf_realdata_screening.R"
)
source(utils_path_hvmf_realdata_screening)

log_line_screening <- function(con, ...) {
  line <- paste0(...)
  writeLines(line, con = con)
  flush(con)
  cat(line, "\n", sep = "")
  flush.console()
}

failed_summary_row <- function(dataset_id,
                               source,
                               station = NA_character_,
                               window_start = NA_character_,
                               window_end = NA_character_,
                               phase = NA_character_,
                               n_raw = NA_integer_,
                               n_clean = NA_integer_,
                               n_final = NA_integer_,
                               scale = NA_real_,
                               B = NA_integer_,
                               n_cores = NA_integer_,
                               elapsed_seconds = NA_real_,
                               processed_csv = NA_character_,
                               result_rds = NA_character_,
                               log_file = NA_character_,
                               error_message = NA_character_) {
  summarize_hvmf_result(
    dataset_id = dataset_id,
    source = source,
    station = station,
    window_start = window_start,
    window_end = window_end,
    phase = phase,
    n_raw = n_raw,
    n_clean = n_clean,
    n_final = n_final,
    scale = scale,
    fit = NULL,
    bootstrap_result = NULL,
    B = B,
    n_cores = n_cores,
    elapsed_seconds = elapsed_seconds,
    processed_csv = processed_csv,
    result_rds = result_rds,
    log_file = log_file,
    status = "failed",
    error_message = error_message
  )
}

build_currents_processed_phase <- function(df_phase,
                                           dataset_id,
                                           source,
                                           station,
                                           phase_label,
                                           bin_value,
                                           projection,
                                           speed_threshold,
                                           preprocessing_notes) {
  if (nrow(df_phase) == 0L) {
    return(data.frame(
      dataset_id = character(0),
      source = character(0),
      datetime = as.POSIXct(character(0), tz = "UTC"),
      raw_speed = numeric(0),
      raw_direction_deg = numeric(0),
      v1 = numeric(0),
      v2 = numeric(0),
      scale = numeric(0),
      z1 = numeric(0),
      z2 = numeric(0),
      x0 = numeric(0),
      x1 = numeric(0),
      x2 = numeric(0),
      minkowski_norm = numeric(0),
      preprocessing_notes = character(0),
      station = character(0),
      phase = character(0),
      bin = integer(0),
      projection = numeric(0),
      speed_threshold = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  phase_speed <- sqrt(df_phase$v1^2 + df_phase$v2^2)
  scale_value <- stats::median(phase_speed, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) {
    stop("Phase-specific scale is not strictly positive and finite.")
  }

  embedded <- embed_components_to_h2(df_phase$v1, df_phase$v2, scale = scale_value)

  data.frame(
    dataset_id = dataset_id,
    source = source,
    datetime = df_phase$datetime,
    raw_speed = df_phase$raw_speed,
    raw_direction_deg = df_phase$raw_direction_deg,
    v1 = df_phase$v1,
    v2 = df_phase$v2,
    scale = embedded$scale,
    z1 = embedded$z1,
    z2 = embedded$z2,
    x0 = embedded$x0,
    x1 = embedded$x1,
    x2 = embedded$x2,
    minkowski_norm = embedded$minkowski_norm,
    preprocessing_notes = preprocessing_notes,
    station = station,
    phase = phase_label,
    bin = if (is.null(bin_value)) NA_integer_ else as.integer(bin_value),
    projection = as.numeric(projection),
    speed_threshold = rep_len(speed_threshold, nrow(df_phase)),
    stringsAsFactors = FALSE
  )
}

prepare_currents_phase_datasets <- function(raw_df,
                                            station,
                                            source,
                                            bin_value = NULL) {
  raw_df <- trim_character_columns(raw_df)
  if (!"date_time" %in% names(raw_df)) {
    stop("The currents dataset is missing the `Date Time` column.")
  }

  raw_df$datetime <- parse_coops_datetime(raw_df$date_time)
  if (any(!is.finite(as.numeric(raw_df$datetime)))) {
    stop("Could not parse the currents timestamps.")
  }

  n_raw <- nrow(raw_df)
  extracted <- extract_currents_components(raw_df)

  if (identical(extracted$method, "speed_direction")) {
    working_df <- extracted$df
    working_df$datetime <- parse_coops_datetime(working_df$date_time)
    working_df$v1 <- extracted$v1
    working_df$v2 <- extracted$v2
    phase_info <- split_phases_by_principal_axis(working_df$v1, working_df$v2)
    phase_method <- "principal_axis_projection"
    projection <- phase_info$projection
    positive_index <- phase_info$positive
    negative_index <- phase_info$negative
    components_note <- sprintf(
      "components from speed/direction using east=speed*sin(theta), north=speed*cos(theta); phase split by PCA principal axis"
    )
  } else {
    working_df <- raw_df
    working_df$v1 <- extracted$v1
    working_df$v2 <- extracted$v2
    working_df$raw_speed <- extracted$raw_speed
    working_df$raw_direction_deg <- extracted$raw_direction_deg
    keep <- is.finite(working_df$v1) &
      is.finite(working_df$v2) &
      is.finite(working_df$raw_speed) &
      (working_df$raw_speed > 0)
    working_df <- working_df[keep, , drop = FALSE]
    phase_method <- "along_channel_sign"
    projection <- working_df$v1
    positive_index <- projection > 0
    negative_index <- projection < 0
    components_note <- sprintf(
      "used along/cross-channel components directly; phase split by sign of along-channel component"
    )
  }

  if (nrow(working_df) == 0L) {
    stop("No finite current observations remained after the basic cleaning step.")
  }

  speed_threshold <- as.numeric(stats::quantile(working_df$raw_speed, probs = 0.25, names = FALSE, type = 8, na.rm = TRUE))
  if (!is.finite(speed_threshold)) {
    stop("Could not compute the 25th percentile speed threshold.")
  }

  keep_slack <- working_df$raw_speed > speed_threshold
  working_df <- working_df[keep_slack, , drop = FALSE]
  projection <- projection[keep_slack]
  positive_index <- positive_index[keep_slack]
  negative_index <- negative_index[keep_slack]
  n_clean <- nrow(working_df)

  if (n_clean == 0L) {
    stop("No current observations remained after removing slack water.")
  }

  positive_df <- working_df[positive_index, , drop = FALSE]
  negative_df <- working_df[negative_index, , drop = FALSE]

  positive_processed <- build_currents_processed_phase(
    df_phase = positive_df,
    dataset_id = if (is.null(bin_value)) {
      sprintf("%s_positive_phase", station)
    } else {
      sprintf("%s_bin%s_positive_phase", station, bin_value)
    },
    source = source,
    station = station,
    phase_label = "positive_phase",
    bin_value = bin_value,
    projection = projection[positive_index],
    speed_threshold = speed_threshold,
    preprocessing_notes = paste(
      components_note,
      sprintf("slack-water filter: raw_speed > %.10f", speed_threshold),
      sep = "; "
    )
  )

  negative_processed <- build_currents_processed_phase(
    df_phase = negative_df,
    dataset_id = if (is.null(bin_value)) {
      sprintf("%s_negative_phase", station)
    } else {
      sprintf("%s_bin%s_negative_phase", station, bin_value)
    },
    source = source,
    station = station,
    phase_label = "negative_phase",
    bin_value = bin_value,
    projection = projection[negative_index],
    speed_threshold = speed_threshold,
    preprocessing_notes = paste(
      components_note,
      sprintf("slack-water filter: raw_speed > %.10f", speed_threshold),
      sep = "; "
    )
  )

  list(
    n_raw = n_raw,
    n_clean = n_clean,
    speed_threshold = speed_threshold,
    phase_method = phase_method,
    subdatasets = list(
      positive = positive_processed,
      negative = negative_processed
    )
  )
}

prepare_ndbc_month_dataset <- function(raw_df,
                                       station,
                                       source,
                                       year,
                                       month) {
  month_df <- raw_df[as.integer(raw_df$month) == as.integer(month), , drop = FALSE]
  n_raw <- nrow(month_df)
  if (n_raw == 0L) {
    stop(sprintf("No raw NDBC rows were found for month %02d.", as.integer(month)))
  }

  month_df <- clean_speed_direction(
    df = month_df,
    speed_col = "wspd",
    direction_col = "wdir",
    sentinel_values = c(99, 999, 9999),
    direction_upper_exclusive = 360
  )
  n_clean <- nrow(month_df)
  if (n_clean == 0L) {
    stop(sprintf("No NDBC rows remained after cleaning for month %02d.", as.integer(month)))
  }

  components <- components_from_speed_direction(month_df$raw_speed, month_df$raw_direction_deg)
  month_df$v1 <- components$v1
  month_df$v2 <- components$v2

  scale_value <- stats::median(month_df$raw_speed, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) {
    stop(sprintf("The scale is not strictly positive for month %02d.", as.integer(month)))
  }

  embedded <- embed_components_to_h2(month_df$v1, month_df$v2, scale = scale_value)
  processed_df <- data.frame(
    dataset_id = sprintf("%s_%d_%02d", station, as.integer(year), as.integer(month)),
    source = source,
    datetime = month_df$datetime,
    raw_speed = month_df$raw_speed,
    raw_direction_deg = month_df$raw_direction_deg,
    v1 = month_df$v1,
    v2 = month_df$v2,
    scale = embedded$scale,
    z1 = embedded$z1,
    z2 = embedded$z2,
    x0 = embedded$x0,
    x1 = embedded$x1,
    x2 = embedded$x2,
    minkowski_norm = embedded$minkowski_norm,
    preprocessing_notes = "components from WSPD/WDIR using east=speed*sin(theta), north=speed*cos(theta); no phase split",
    station = station,
    year = rep_len(as.integer(year), n_clean),
    month = rep_len(as.integer(month), n_clean),
    stringsAsFactors = FALSE
  )

  list(
    n_raw = n_raw,
    n_clean = n_clean,
    processed_df = processed_df
  )
}

run_single_subdataset_analysis <- function(metadata,
                                           processed_df,
                                           directories,
                                           B,
                                           n_cores,
                                           seed,
                                           day_stride = 5L) {
  dataset_id <- metadata$dataset_id
  log_file <- file.path(directories$logs, sprintf("%s_log.txt", dataset_id))
  processed_csv <- file.path(directories$processed, sprintf("%s_h2.csv", dataset_id))
  result_rds <- file.path(directories$results, sprintf("%s_cvm_simple_plugin_result.rds", dataset_id))
  start_time <- Sys.time()
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  log_line_screening(log_con, "dataset_id: ", dataset_id)
  log_line_screening(log_con, "start_time: ", format_timestamp_screening(start_time))
  log_line_screening(log_con, "source: ", metadata$source)
  log_line_screening(log_con, "station: ", metadata$station %||% "")
  log_line_screening(log_con, "source_url: ", metadata$source_url %||% "")
  log_line_screening(log_con, "window_start: ", metadata$window_start %||% "")
  log_line_screening(log_con, "window_end: ", metadata$window_end %||% "")
  log_line_screening(log_con, "phase: ", metadata$phase %||% "")
  log_line_screening(log_con, "bin: ", if (is.null(metadata$bin)) "" else as.character(metadata$bin))
  log_line_screening(log_con, "download_parameters: ", metadata$download_parameters %||% "")
  log_line_screening(log_con, "n_raw: ", metadata$n_raw %||% NA_integer_)
  log_line_screening(log_con, "n_clean: ", metadata$n_clean %||% NA_integer_)
  log_line_screening(log_con, "n_after_phase_split: ", if (is.null(processed_df)) NA_integer_ else nrow(processed_df))
  log_line_screening(log_con, "configuration: statistic=CvM, null=simple plugin, B=", as.integer(B), ", n_cores=", as.integer(n_cores), ", seed=", as.integer(seed))

  summary_row <- NULL
  warnings_all <- character(0)

  tryCatch(
    {
      if (is.null(processed_df)) {
        stop(metadata$error_message %||% "Processed data is unavailable.")
      }

      processed_df <- processed_df[order(processed_df$datetime), , drop = FALSE]
      write_processed_dataset_csv(processed_df, processed_csv)

      if (nrow(processed_df) < 80L) {
        stop(sprintf("Only %d observations remained; fewer than the required 80.", nrow(processed_df)))
      }

      day_stride_result <- thin_by_day_stride_if_daily_coverage(
        df = processed_df,
        day_stride = day_stride,
        time_col = "datetime"
      )
      processed_df <- day_stride_result$data
      if (nrow(processed_df) < 80L) {
        stop(sprintf(
          "Only %d observations remained after the %d-day temporal filter; fewer than the required 80.",
          nrow(processed_df),
          as.integer(day_stride)
        ))
      }
      n_final <- nrow(processed_df)
      scale_value <- as.numeric(processed_df$scale[[1L]])

      stopifnot(nrow(processed_df) >= 80L)
      stopifnot(all(is.finite(processed_df$x0)))
      stopifnot(all(is.finite(processed_df$x1)))
      stopifnot(all(is.finite(processed_df$x2)))

      verification <- verify_h2_matrix_screening(processed_df[, c("x0", "x1", "x2"), drop = FALSE], tol = 1e-8)
      stopifnot(max(abs(verification$minkowski_norm + 1)) < 1e-8)
      stopifnot(all(processed_df$x0 > 0))

      analysis <- run_hvmf_simple_plugin_cvm(
        X = processed_df[, c("x0", "x1", "x2"), drop = FALSE],
        B = B,
        n_cores = n_cores,
        seed = seed
      )
      warnings_all <- unique(c(warnings_all, analysis$warnings))

      fit <- analysis$fit
      stopifnot(is.finite(fit$kappa), fit$kappa > 0)
      stopifnot(max(abs(-fit$mu[1]^2 + fit$mu[2]^2 + fit$mu[3]^2 + 1)) < 1e-8)

      saveRDS(
        list(
          metadata = metadata,
          fit = fit,
          theta0 = analysis$theta0,
          bootstrap_result = analysis$result,
          warnings = warnings_all
        ),
        file = result_rds
      )

      end_time <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
      aux <- compute_hvmf_auxiliary_parameters(fit$mu)

      log_line_screening(log_con, "processed_csv: ", processed_csv)
      log_line_screening(log_con, "result_rds: ", result_rds)
      log_line_screening(
        log_con,
        "day_stride_filter: ",
        if (isTRUE(day_stride_result$applied)) {
          sprintf(
            "applied_every_%d_days (unique_days %d -> %d)",
            as.integer(day_stride_result$day_stride),
            as.integer(day_stride_result$unique_days_before),
            as.integer(day_stride_result$unique_days_after)
          )
        } else {
          "not_applied"
        }
      )
      log_line_screening(log_con, "n_final: ", n_final)
      log_line_screening(log_con, "scale: ", sprintf("%.10f", scale_value))
      log_line_screening(log_con, "max_abs_minkowski_norm_plus_1: ", sprintf("%.12e", max(abs(verification$minkowski_norm + 1))))
      log_line_screening(
        log_con,
        "mle: mu=(",
        paste(sprintf("%.10f", fit$mu), collapse = ", "),
        "), kappa=",
        sprintf("%.10f", fit$kappa)
      )
      log_line_screening(log_con, "sinh_chi_hat: ", sprintf("%.10f", aux$sinh_chi_hat))
      log_line_screening(log_con, "theta_deg_hat: ", sprintf("%.10f", aux$theta_deg_hat))
      log_line_screening(log_con, "CvM_observed: ", sprintf("%.10f", analysis$result$observed$cvm$statistic))
      log_line_screening(log_con, "p_value: ", sprintf("%.10f", analysis$result$inference$cvm$p_value))
      log_line_screening(log_con, "B: ", as.integer(B))
      log_line_screening(log_con, "n_cores: ", as.integer(n_cores))
      log_line_screening(log_con, "elapsed_seconds: ", sprintf("%.3f", elapsed_seconds))
      if (length(warnings_all) > 0L) {
        log_line_screening(log_con, "warnings: ", paste(unique(warnings_all), collapse = " | "))
      } else {
        log_line_screening(log_con, "warnings: none")
      }
      log_line_screening(log_con, "status: ok")
      log_line_screening(log_con, "end_time: ", format_timestamp_screening(end_time))

      summary_row <- summarize_hvmf_result(
        dataset_id = dataset_id,
        source = metadata$source,
        station = metadata$station,
        window_start = metadata$window_start,
        window_end = metadata$window_end,
        phase = metadata$phase,
        n_raw = metadata$n_raw,
        n_clean = metadata$n_clean,
        n_final = n_final,
        scale = scale_value,
        fit = fit,
        bootstrap_result = analysis$result,
        B = B,
        n_cores = n_cores,
        elapsed_seconds = elapsed_seconds,
        processed_csv = processed_csv,
        result_rds = result_rds,
        log_file = log_file,
        status = "ok",
        error_message = NA_character_
      )
    },
    error = function(e) {
      end_time <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
      error_message <- conditionMessage(e)

      log_line_screening(log_con, "status: failed")
      log_line_screening(log_con, "error: ", error_message)
      if (length(warnings_all) > 0L) {
        log_line_screening(log_con, "warnings: ", paste(unique(warnings_all), collapse = " | "))
      } else {
        log_line_screening(log_con, "warnings: none")
      }
      log_line_screening(log_con, "elapsed_seconds: ", sprintf("%.3f", elapsed_seconds))
      log_line_screening(log_con, "end_time: ", format_timestamp_screening(end_time))

      summary_row <- failed_summary_row(
        dataset_id = dataset_id,
        source = metadata$source,
        station = metadata$station,
        window_start = metadata$window_start,
        window_end = metadata$window_end,
        phase = metadata$phase,
        n_raw = metadata$n_raw,
        n_clean = metadata$n_clean,
        n_final = if (is.null(processed_df)) NA_integer_ else min(500L, nrow(processed_df)),
        scale = if (is.null(processed_df) || !"scale" %in% names(processed_df) || nrow(processed_df) == 0L) NA_real_ else as.numeric(processed_df$scale[[1L]]),
        B = B,
        n_cores = n_cores,
        elapsed_seconds = elapsed_seconds,
        processed_csv = processed_csv,
        result_rds = result_rds,
        log_file = log_file,
        error_message = error_message
      )
    }
  )

  summary_row
}

attempt_coops_download <- function(station, windows, bins = list(NULL)) {
  attempts <- list()

  for (window in windows) {
    begin_date <- window[[1L]]
    end_date <- window[[2L]]
    for (bin_value in bins) {
      download_result <- download_coops_currents(
        station = station,
        begin_date = begin_date,
        end_date = end_date,
        bin = bin_value
      )
      attempts[[length(attempts) + 1L]] <- list(
        begin_date = begin_date,
        end_date = end_date,
        bin = bin_value,
        ok = download_result$ok,
        url = download_result$url,
        error_message = download_result$error_message
      )

      if (isTRUE(download_result$ok) && nrow(download_result$data) > 0L) {
        return(list(
          ok = TRUE,
          download = download_result,
          begin_date = begin_date,
          end_date = end_date,
          bin = bin_value,
          attempts = attempts
        ))
      }
    }
  }

  list(
    ok = FALSE,
    attempts = attempts,
    error_message = paste(
      vapply(attempts, function(attempt) {
        sprintf(
          "[%s,%s,bin=%s] %s",
          attempt$begin_date,
          attempt$end_date,
          if (is.null(attempt$bin)) "NA" else attempt$bin,
          attempt$error_message %||% "unknown error"
        )
      }, character(1)),
      collapse = " || "
    )
  )
}

run_hvmf_realdata_screening <- function(B = 1000L,
                                        n_cores = 3L,
                                        base_dir = "hvmf_realdata_screening",
                                        day_stride = 5L) {
  if (as.integer(n_cores) != 3L) {
    warning("The requested workflow was calibrated for `n_cores = 3`; continuing with the user-supplied value.")
  }

  directories <- ensure_screening_directories(base_dir = base_dir)
  summary_rows <- list()
  dataset_seed <- 2026052400L

  currents_candidates <- list(
    list(
      station = "n03020",
      source = "NOAA CO-OPS currents",
      windows = list(c("20250401", "20250414"), c("20250701", "20250714")),
      bins = list(NULL)
    ),
    list(
      station = "ca0101",
      source = "NOAA CO-OPS currents",
      windows = list(c("20250401", "20250414"), c("20250701", "20250714")),
      bins = list(NULL)
    ),
    list(
      station = "DEB2113",
      source = "NOAA CO-OPS historic currents",
      windows = list(c("20210701", "20210714")),
      bins = as.list(as.integer(1:10))
    )
  )

  for (candidate in currents_candidates) {
    cat(sprintf("Processing currents candidate %s\n", candidate$station))
    attempt <- attempt_coops_download(
      station = candidate$station,
      windows = candidate$windows,
      bins = candidate$bins
    )

    if (!isTRUE(attempt$ok)) {
      dataset_id <- candidate$station
      log_file <- file.path(directories$logs, sprintf("%s_log.txt", dataset_id))
      writeLines(
        c(
          sprintf("dataset_id: %s", dataset_id),
          sprintf("start_time: %s", format_timestamp_screening(Sys.time())),
          sprintf("source: %s", candidate$source),
          sprintf("status: failed"),
          sprintf("error: %s", attempt$error_message),
          sprintf("end_time: %s", format_timestamp_screening(Sys.time()))
        ),
        con = log_file
      )

      summary_rows[[length(summary_rows) + 1L]] <- failed_summary_row(
        dataset_id = dataset_id,
        source = candidate$source,
        station = candidate$station,
        B = B,
        n_cores = n_cores,
        log_file = log_file,
        error_message = attempt$error_message
      )
      next
    }

    prepared <- tryCatch(
      prepare_currents_phase_datasets(
        raw_df = attempt$download$data,
        station = candidate$station,
        source = candidate$source,
        bin_value = attempt$bin
      ),
      error = function(e) {
        e
      }
    )
    prepared_error <- NULL
    if (inherits(prepared, "error")) {
      prepared_error <- conditionMessage(prepared)
    }
    if (!is.null(prepared_error)) {
      dataset_id <- candidate$station
      log_file <- file.path(directories$logs, sprintf("%s_log.txt", dataset_id))
      writeLines(
        c(
          sprintf("dataset_id: %s", dataset_id),
          sprintf("start_time: %s", format_timestamp_screening(Sys.time())),
          sprintf("source: %s", candidate$source),
          sprintf("source_url: %s", attempt$download$url),
          sprintf("status: failed"),
          sprintf("error: %s", prepared_error),
          sprintf("end_time: %s", format_timestamp_screening(Sys.time()))
        ),
        con = log_file
      )

      summary_rows[[length(summary_rows) + 1L]] <- failed_summary_row(
        dataset_id = dataset_id,
        source = candidate$source,
        station = candidate$station,
        window_start = attempt$begin_date,
        window_end = attempt$end_date,
        B = B,
        n_cores = n_cores,
        log_file = log_file,
        error_message = prepared_error
      )
      next
    }

    for (phase_name in c("positive", "negative")) {
      dataset_seed <- dataset_seed + 1L
      processed_df <- prepared$subdatasets[[phase_name]]
      phase_label <- sprintf("%s_phase", phase_name)
      dataset_id <- if (is.null(attempt$bin)) {
        sprintf("%s_%s", candidate$station, phase_label)
      } else {
        sprintf("%s_bin%s_%s", candidate$station, attempt$bin, phase_label)
      }
      metadata <- list(
        dataset_id = dataset_id,
        source = candidate$source,
        station = candidate$station,
        source_url = attempt$download$url,
        window_start = attempt$begin_date,
        window_end = attempt$end_date,
        phase = phase_label,
        bin = attempt$bin,
        download_parameters = sprintf(
          "product=currents, station=%s, begin_date=%s, end_date=%s, units=metric, time_zone=gmt, format=csv%s",
          candidate$station,
          attempt$begin_date,
          attempt$end_date,
          if (is.null(attempt$bin)) "" else paste0(", bin=", attempt$bin)
        ),
        n_raw = prepared$n_raw,
        n_clean = prepared$n_clean
      )
      summary_rows[[length(summary_rows) + 1L]] <- run_single_subdataset_analysis(
        metadata = metadata,
        processed_df = processed_df,
        directories = directories,
        B = B,
        n_cores = n_cores,
        seed = dataset_seed,
        day_stride = day_stride
      )
    }
  }

  ndbc_candidates <- list(
    list(
      station = "51001",
      source = "NOAA NDBC stdmet",
      years = c(2024L, 2023L),
      months = c(7L, 8L)
    ),
    list(
      station = "46047",
      source = "NOAA NDBC stdmet",
      years = c(2024L, 2023L),
      months = c(4L, 5L)
    )
  )

  for (candidate in ndbc_candidates) {
    cat(sprintf("Processing NDBC candidate %s\n", candidate$station))
    year_download <- NULL
    year_used <- NULL

    for (year in candidate$years) {
      download_result <- download_ndbc_stdmet(station = candidate$station, year = year)
      if (isTRUE(download_result$ok)) {
        year_download <- download_result
        year_used <- year
        break
      }
    }

    if (is.null(year_download)) {
      dataset_id <- candidate$station
      log_file <- file.path(directories$logs, sprintf("%s_log.txt", dataset_id))
      writeLines(
        c(
          sprintf("dataset_id: %s", dataset_id),
          sprintf("start_time: %s", format_timestamp_screening(Sys.time())),
          sprintf("source: %s", candidate$source),
          sprintf("status: failed"),
          sprintf("error: could not download any requested NDBC year"),
          sprintf("end_time: %s", format_timestamp_screening(Sys.time()))
        ),
        con = log_file
      )

      summary_rows[[length(summary_rows) + 1L]] <- failed_summary_row(
        dataset_id = dataset_id,
        source = candidate$source,
        station = candidate$station,
        B = B,
        n_cores = n_cores,
        log_file = log_file,
        error_message = "Could not download any requested NDBC year."
      )
      next
    }

    raw_df <- tryCatch(
      parse_ndbc_stdmet(year_download$file),
      error = function(e) {
        e
      }
    )
    if (inherits(raw_df, "error")) {
      dataset_id <- candidate$station
      log_file <- file.path(directories$logs, sprintf("%s_log.txt", dataset_id))
      writeLines(
        c(
          sprintf("dataset_id: %s", dataset_id),
          sprintf("start_time: %s", format_timestamp_screening(Sys.time())),
          sprintf("source: %s", candidate$source),
          sprintf("source_url: %s", year_download$url),
          sprintf("status: failed"),
          sprintf("error: %s", conditionMessage(raw_df)),
          sprintf("end_time: %s", format_timestamp_screening(Sys.time()))
        ),
        con = log_file
      )

      summary_rows[[length(summary_rows) + 1L]] <- failed_summary_row(
        dataset_id = dataset_id,
        source = candidate$source,
        station = candidate$station,
        B = B,
        n_cores = n_cores,
        log_file = log_file,
        error_message = conditionMessage(raw_df)
      )
      next
    }
    for (month in candidate$months) {
      dataset_seed <- dataset_seed + 1L
      month_dataset_id <- sprintf("%s_%d_%02d", candidate$station, year_used, as.integer(month))
      prepared <- tryCatch(
        prepare_ndbc_month_dataset(
          raw_df = raw_df,
          station = candidate$station,
          source = candidate$source,
          year = year_used,
          month = month
        ),
        error = function(e) {
          list(error = conditionMessage(e))
        }
      )

      if (!is.null(prepared$error)) {
        log_file <- file.path(directories$logs, sprintf("%s_log.txt", month_dataset_id))
        writeLines(
          c(
            sprintf("dataset_id: %s", month_dataset_id),
            sprintf("start_time: %s", format_timestamp_screening(Sys.time())),
            sprintf("source: %s", candidate$source),
            sprintf("source_url: %s", year_download$url),
            sprintf("status: failed"),
            sprintf("error: %s", prepared$error),
            sprintf("end_time: %s", format_timestamp_screening(Sys.time()))
          ),
          con = log_file
        )

        summary_rows[[length(summary_rows) + 1L]] <- failed_summary_row(
          dataset_id = month_dataset_id,
          source = candidate$source,
          station = candidate$station,
          window_start = sprintf("%d-%02d-01", year_used, month),
          window_end = sprintf("%d-%02d-01", year_used, month),
          B = B,
          n_cores = n_cores,
          log_file = log_file,
          error_message = prepared$error
        )
        next
      }

      month_start <- format(min(prepared$processed_df$datetime), "%Y-%m-%d %H:%M:%S", tz = "UTC")
      month_end <- format(max(prepared$processed_df$datetime), "%Y-%m-%d %H:%M:%S", tz = "UTC")
      metadata <- list(
        dataset_id = prepared$processed_df$dataset_id[[1L]],
        source = candidate$source,
        station = candidate$station,
        source_url = year_download$url,
        window_start = month_start,
        window_end = month_end,
        phase = NA_character_,
        bin = NULL,
        download_parameters = sprintf("station=%s, year=%d, month=%02d", candidate$station, as.integer(year_used), as.integer(month)),
        n_raw = prepared$n_raw,
        n_clean = prepared$n_clean
      )

      summary_rows[[length(summary_rows) + 1L]] <- run_single_subdataset_analysis(
        metadata = metadata,
        processed_df = prepared$processed_df,
        directories = directories,
        B = B,
        n_cores = n_cores,
        seed = dataset_seed,
        day_stride = day_stride
      )
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- order_screening_summary(summary_df)
  summary_csv <- file.path(directories$results, "hvmf_realdata_screening_summary.csv")
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  list(
    summary = summary_df,
    summary_csv = summary_csv,
    directories = directories
  )
}

parse_simple_log_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  output <- list()

  for (line in lines) {
    if (!grepl(": ", line, fixed = TRUE)) {
      next
    }
    pieces <- strsplit(line, ": ", fixed = TRUE)[[1]]
    key <- pieces[[1L]]
    value <- paste(pieces[-1L], collapse = ": ")
    output[[key]] <- value
  }

  output
}

rebuild_hvmf_realdata_screening_summary <- function(base_dir = "hvmf_realdata_screening") {
  directories <- screening_directories(base_dir = base_dir)
  result_files <- list.files(directories$results, pattern = "_cvm_simple_plugin_result\\.rds$", full.names = TRUE)
  log_files <- list.files(directories$logs, pattern = "_log\\.txt$", full.names = TRUE)
  summary_rows <- list()

  for (result_file in result_files) {
    artifact <- readRDS(result_file)
    metadata <- artifact$metadata
    processed_csv <- file.path(directories$processed, sprintf("%s_h2.csv", metadata$dataset_id))
    log_file <- file.path(directories$logs, sprintf("%s_log.txt", metadata$dataset_id))
    processed_df <- utils::read.csv(processed_csv, stringsAsFactors = FALSE)

    summary_rows[[length(summary_rows) + 1L]] <- summarize_hvmf_result(
      dataset_id = metadata$dataset_id,
      source = metadata$source,
      station = metadata$station,
      window_start = metadata$window_start,
      window_end = metadata$window_end,
      phase = metadata$phase,
      n_raw = metadata$n_raw,
      n_clean = metadata$n_clean,
      n_final = nrow(processed_df),
      scale = as.numeric(processed_df$scale[[1L]]),
      fit = artifact$fit,
      bootstrap_result = artifact$bootstrap_result,
      B = artifact$bootstrap_result$bootstrap$B,
      n_cores = artifact$bootstrap_result$diagnostics$n_cores,
      elapsed_seconds = artifact$bootstrap_result$diagnostics$elapsed_seconds,
      processed_csv = processed_csv,
      result_rds = result_file,
      log_file = log_file,
      status = "ok",
      error_message = NA_character_
    )
  }

  success_dataset_ids <- vapply(summary_rows, function(row) row$dataset_id[[1L]], character(1))
  failed_logs <- log_files[!sub("_log\\.txt$", "", basename(log_files)) %in% success_dataset_ids]

  for (log_file in failed_logs) {
    parsed <- parse_simple_log_file(log_file)
    summary_rows[[length(summary_rows) + 1L]] <- failed_summary_row(
      dataset_id = parsed$dataset_id %||% sub("_log\\.txt$", "", basename(log_file)),
      source = parsed$source %||% NA_character_,
      station = parsed$station %||% NA_character_,
      window_start = parsed$window_start %||% NA_character_,
      window_end = parsed$window_end %||% NA_character_,
      phase = parsed$phase %||% NA_character_,
      n_raw = suppressWarnings(as.integer(parsed$n_raw %||% NA_character_)),
      n_clean = suppressWarnings(as.integer(parsed$n_clean %||% NA_character_)),
      n_final = suppressWarnings(as.integer(parsed$n_final %||% NA_character_)),
      scale = suppressWarnings(as.numeric(parsed$scale %||% NA_character_)),
      B = suppressWarnings(as.integer(parsed$B %||% 1000L)),
      n_cores = suppressWarnings(as.integer(parsed$n_cores %||% 3L)),
      elapsed_seconds = suppressWarnings(as.numeric(parsed$elapsed_seconds %||% NA_character_)),
      processed_csv = NA_character_,
      result_rds = NA_character_,
      log_file = log_file,
      error_message = parsed$error %||% "Unknown failure recorded in log."
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- order_screening_summary(summary_df)
  summary_csv <- file.path(directories$results, "hvmf_realdata_screening_summary.csv")
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  list(
    summary = summary_df,
    summary_csv = summary_csv
  )
}

if (sys.nframe() == 0L) {
  run_hvmf_realdata_screening()
}
