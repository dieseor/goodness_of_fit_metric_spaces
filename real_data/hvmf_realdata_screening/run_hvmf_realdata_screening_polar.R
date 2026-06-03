resolve_hvmf_realdata_screening_polar_path <- function(...) {
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

utils_path_hvmf_realdata_screening_polar <- resolve_hvmf_realdata_screening_polar_path(
  "hvmf_realdata_screening",
  "utils_hvmf_realdata_screening.R"
)
source(utils_path_hvmf_realdata_screening_polar)

polar_screening_directories <- function(base_dir = "hvmf_realdata_screening") {
  list(
    base = base_dir,
    processed = file.path(base_dir, "processed"),
    processed_polar = file.path(base_dir, "processed_polar"),
    results = file.path(base_dir, "results"),
    results_polar = file.path(base_dir, "results_polar"),
    logs_polar = file.path(base_dir, "logs_polar")
  )
}

ensure_polar_screening_directories <- function(base_dir = "hvmf_realdata_screening") {
  dirs <- polar_screening_directories(base_dir = base_dir)
  for (path in unname(unlist(dirs))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}

log_line_polar <- function(con, ...) {
  line <- paste0(...)
  writeLines(line, con = con)
  flush(con)
  cat(line, "\n", sep = "")
  flush.console()
}

polar_dataset_files <- function() {
  c(
    "n03020_negative_phase_h2.csv",
    "n03020_positive_phase_h2.csv",
    "ca0101_negative_phase_h2.csv",
    "ca0101_positive_phase_h2.csv",
    "51001_2024_07_h2.csv",
    "51001_2024_08_h2.csv",
    "46047_2024_04_h2.csv",
    "46047_2024_05_h2.csv"
  )
}

read_component_processed_dataset <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  df <- trim_character_columns(df)
  numeric_columns <- intersect(
    c(
      "raw_speed",
      "raw_direction_deg",
      "v1",
      "v2",
      "scale",
      "z1",
      "z2",
      "x0",
      "x1",
      "x2",
      "minkowski_norm",
      "projection",
      "speed_threshold",
      "year",
      "month"
    ),
    names(df)
  )
  df <- coerce_numeric_columns(df, numeric_columns)
  if ("datetime" %in% names(df)) {
    df$datetime <- as.POSIXct(df$datetime, tz = "UTC")
  }
  df
}

prepare_polar_embedding_dataset <- function(df_component) {
  required_columns <- c("dataset_id", "source", "datetime", "v1", "v2")
  missing_columns <- setdiff(required_columns, names(df_component))
  if (length(missing_columns) > 0L) {
    stop(sprintf(
      "The component CSV is missing required columns: %s",
      paste(missing_columns, collapse = ", ")
    ))
  }

  speed <- sqrt(df_component$v1^2 + df_component$v2^2)
  stopifnot(all(is.finite(speed)))
  stopifnot(all(speed > 0))

  unit1 <- df_component$v1 / speed
  unit2 <- df_component$v2 / speed

  speed_mean <- mean(speed, na.rm = TRUE)
  stopifnot(is.finite(speed_mean), speed_mean > 0)

  u <- speed / speed_mean
  x0_polar <- cosh(u)
  sinh_u <- sinh(u)
  x1_polar <- sinh_u * unit1
  x2_polar <- sinh_u * unit2
  minkowski_norm_polar <- -x0_polar^2 + x1_polar^2 + x2_polar^2

  stopifnot(max(abs(minkowski_norm_polar + 1)) < 1e-8)
  stopifnot(all(x0_polar > 0))

  output <- df_component

  rename_map <- intersect(
    c("x0", "x1", "x2", "minkowski_norm"),
    names(output)
  )
  for (column in rename_map) {
    names(output)[names(output) == column] <- paste0(column, "_component")
  }

  output$speed_polar <- speed
  output$speed_mean_polar <- rep_len(speed_mean, nrow(output))
  output$speed_scaled_polar <- u
  output$unit1 <- unit1
  output$unit2 <- unit2
  output$x0 <- x0_polar
  output$x1 <- x1_polar
  output$x2 <- x2_polar
  output$minkowski_norm <- minkowski_norm_polar
  output$embedding <- "polar_jensen"

  output
}

failed_polar_summary_row <- function(dataset_id,
                                     source = NA_character_,
                                     station = NA_character_,
                                     phase = NA_character_,
                                     n_final = NA_integer_,
                                     speed_mean_polar = NA_real_,
                                     B = NA_integer_,
                                     n_cores = NA_integer_,
                                     elapsed_seconds = NA_real_,
                                     processed_csv = NA_character_,
                                     result_rds = NA_character_,
                                     log_file = NA_character_,
                                     error_message = NA_character_) {
  data.frame(
    dataset_id = dataset_id,
    source = source,
    station = station,
    phase = phase,
    n_final = n_final,
    speed_mean_polar = speed_mean_polar,
    kappa_hat = NA_real_,
    mu_hat_x0 = NA_real_,
    mu_hat_x1 = NA_real_,
    mu_hat_x2 = NA_real_,
    sinh_chi_hat = NA_real_,
    theta_deg_hat = NA_real_,
    statistic_CvM = NA_real_,
    p_value = NA_real_,
    B = as.integer(B),
    n_cores = as.integer(n_cores),
    elapsed_seconds = elapsed_seconds,
    processed_csv = processed_csv,
    result_rds = result_rds,
    log_file = log_file,
    status = "failed",
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

ok_polar_summary_row <- function(dataset_id,
                                 source,
                                 station,
                                 phase,
                                 n_final,
                                 speed_mean_polar,
                                 fit,
                                 bootstrap_result,
                                 B,
                                 n_cores,
                                 elapsed_seconds,
                                 processed_csv,
                                 result_rds,
                                 log_file) {
  aux <- compute_hvmf_auxiliary_parameters(fit$mu)

  data.frame(
    dataset_id = dataset_id,
    source = source,
    station = station,
    phase = phase,
    n_final = n_final,
    speed_mean_polar = speed_mean_polar,
    kappa_hat = fit$kappa,
    mu_hat_x0 = fit$mu[1],
    mu_hat_x1 = fit$mu[2],
    mu_hat_x2 = fit$mu[3],
    sinh_chi_hat = aux$sinh_chi_hat,
    theta_deg_hat = aux$theta_deg_hat,
    statistic_CvM = bootstrap_result$observed$cvm$statistic,
    p_value = bootstrap_result$inference$cvm$p_value,
    B = as.integer(B),
    n_cores = as.integer(n_cores),
    elapsed_seconds = elapsed_seconds,
    processed_csv = processed_csv,
    result_rds = result_rds,
    log_file = log_file,
    status = "ok",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_single_polar_dataset <- function(component_csv,
                                     directories,
                                     B,
                                     n_cores,
                                     seed,
                                     day_stride = 5L) {
  start_time <- Sys.time()
  component_path <- file.path(directories$processed, component_csv)
  component_df <- read_component_processed_dataset(component_path)
  dataset_id <- as.character(component_df$dataset_id[[1L]])
  source_value <- if ("source" %in% names(component_df)) as.character(component_df$source[[1L]]) else NA_character_
  station_value <- if ("station" %in% names(component_df)) as.character(component_df$station[[1L]]) else NA_character_
  phase_value <- if ("phase" %in% names(component_df)) as.character(component_df$phase[[1L]]) else NA_character_

  processed_polar_csv <- file.path(directories$processed_polar, sprintf("%s_h2_polar.csv", dataset_id))
  result_rds <- file.path(directories$results_polar, sprintf("%s_cvm_simple_plugin_polar_result.rds", dataset_id))
  log_file <- file.path(directories$logs_polar, sprintf("%s_polar_log.txt", dataset_id))
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  log_line_polar(log_con, "dataset_id: ", dataset_id)
  log_line_polar(log_con, "start_time: ", format_timestamp_screening(start_time))
  log_line_polar(log_con, "embedding: polar_jensen")
  log_line_polar(log_con, "component_csv: ", component_path)
  log_line_polar(log_con, "n_final: ", nrow(component_df))
  log_line_polar(log_con, "configuration: statistic=CvM, null=simple plugin, B=", as.integer(B), ", n_cores=", as.integer(n_cores), ", seed=", as.integer(seed))

  warnings_all <- character(0)
  summary_row <- NULL

  tryCatch(
    {
      day_stride_result <- thin_by_day_stride_if_daily_coverage(
        df = component_df,
        day_stride = day_stride,
        time_col = "datetime"
      )
      component_df <- day_stride_result$data
      if (nrow(component_df) < 80L) {
        stop(sprintf(
          "Only %d observations remained after the %d-day temporal filter; fewer than the required 80.",
          nrow(component_df),
          as.integer(day_stride)
        ))
      }

      polar_df <- prepare_polar_embedding_dataset(component_df)
      speed_mean_polar <- as.numeric(polar_df$speed_mean_polar[[1L]])
      verification <- verify_h2_matrix_screening(polar_df[, c("x0", "x1", "x2"), drop = FALSE], tol = 1e-8)

      write_processed_dataset_csv(polar_df, processed_polar_csv)

      analysis <- run_hvmf_simple_plugin_cvm(
        X = polar_df[, c("x0", "x1", "x2"), drop = FALSE],
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
          dataset_id = dataset_id,
          embedding = "polar_jensen",
          component_csv = component_path,
          processed_polar_csv = processed_polar_csv,
          fit = fit,
          theta0 = analysis$theta0,
          bootstrap_result = analysis$result,
          warnings = warnings_all
        ),
        file = result_rds
      )

      aux <- compute_hvmf_auxiliary_parameters(fit$mu)
      end_time <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))

      log_line_polar(log_con, "processed_csv: ", processed_polar_csv)
      log_line_polar(log_con, "result_rds: ", result_rds)
      log_line_polar(
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
      log_line_polar(log_con, "n_final: ", nrow(polar_df))
      log_line_polar(log_con, "speed_mean_polar: ", sprintf("%.10f", speed_mean_polar))
      log_line_polar(log_con, "max_abs_minkowski_norm_plus_1: ", sprintf("%.12e", max(abs(verification$minkowski_norm + 1))))
      log_line_polar(
        log_con,
        "mu_hat: (",
        paste(sprintf("%.10f", fit$mu), collapse = ", "),
        ")"
      )
      log_line_polar(log_con, "kappa_hat: ", sprintf("%.10f", fit$kappa))
      log_line_polar(log_con, "sinh_chi_hat: ", sprintf("%.10f", aux$sinh_chi_hat))
      log_line_polar(log_con, "theta_deg_hat: ", sprintf("%.10f", aux$theta_deg_hat))
      log_line_polar(log_con, "CvM_observed: ", sprintf("%.10f", analysis$result$observed$cvm$statistic))
      log_line_polar(log_con, "p_value: ", sprintf("%.10f", analysis$result$inference$cvm$p_value))
      log_line_polar(log_con, "B: ", as.integer(B))
      log_line_polar(log_con, "n_cores: ", as.integer(n_cores))
      log_line_polar(log_con, "elapsed_seconds: ", sprintf("%.3f", elapsed_seconds))
      if (length(warnings_all) > 0L) {
        log_line_polar(log_con, "warnings: ", paste(unique(warnings_all), collapse = " | "))
      } else {
        log_line_polar(log_con, "warnings: none")
      }
      log_line_polar(log_con, "status: ok")
      log_line_polar(log_con, "end_time: ", format_timestamp_screening(end_time))

      summary_row <- ok_polar_summary_row(
        dataset_id = dataset_id,
        source = source_value,
        station = station_value,
        phase = phase_value,
        n_final = nrow(polar_df),
        speed_mean_polar = speed_mean_polar,
        fit = fit,
        bootstrap_result = analysis$result,
        B = B,
        n_cores = n_cores,
        elapsed_seconds = elapsed_seconds,
        processed_csv = processed_polar_csv,
        result_rds = result_rds,
        log_file = log_file
      )
    },
    error = function(e) {
      end_time <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
      speed_mean_value <- if (exists("speed_mean_polar", inherits = FALSE)) speed_mean_polar else NA_real_

      log_line_polar(log_con, "status: failed")
      log_line_polar(log_con, "error: ", conditionMessage(e))
      if (length(warnings_all) > 0L) {
        log_line_polar(log_con, "warnings: ", paste(unique(warnings_all), collapse = " | "))
      } else {
        log_line_polar(log_con, "warnings: none")
      }
      log_line_polar(log_con, "elapsed_seconds: ", sprintf("%.3f", elapsed_seconds))
      log_line_polar(log_con, "end_time: ", format_timestamp_screening(end_time))

      summary_row <- failed_polar_summary_row(
        dataset_id = dataset_id,
        source = source_value,
        station = station_value,
        phase = phase_value,
        n_final = nrow(component_df),
        speed_mean_polar = speed_mean_value,
        B = B,
        n_cores = n_cores,
        elapsed_seconds = elapsed_seconds,
        processed_csv = processed_polar_csv,
        result_rds = result_rds,
        log_file = log_file,
        error_message = conditionMessage(e)
      )
    }
  )

  summary_row
}

make_component_vs_polar_comparison <- function(component_summary_csv,
                                               polar_summary_df,
                                               output_csv) {
  component_summary <- utils::read.csv(component_summary_csv, stringsAsFactors = FALSE)
  component_summary <- component_summary[component_summary$status == "ok", , drop = FALSE]
  polar_ok <- polar_summary_df[polar_summary_df$status == "ok", , drop = FALSE]

  merged <- merge(
    component_summary[, c("dataset_id", "statistic_CvM", "p_value", "kappa_hat"), drop = FALSE],
    polar_ok[, c("dataset_id", "statistic_CvM", "p_value", "kappa_hat"), drop = FALSE],
    by = "dataset_id",
    suffixes = c("_component", "_polar"),
    all = FALSE
  )

  names(merged)[names(merged) == "statistic_CvM_component"] <- "CvM_component"
  names(merged)[names(merged) == "statistic_CvM_polar"] <- "CvM_polar"
  names(merged)[names(merged) == "kappa_hat_component"] <- "kappa_component"
  names(merged)[names(merged) == "kappa_hat_polar"] <- "kappa_polar"

  merged$delta_CvM <- merged$CvM_polar - merged$CvM_component
  merged$delta_p_value <- merged$p_value_polar - merged$p_value_component
  merged <- merged[order(-merged$p_value_polar, merged$dataset_id), , drop = FALSE]

  utils::write.csv(merged, output_csv, row.names = FALSE)
  merged
}

run_hvmf_realdata_screening_polar <- function(B = 1000L,
                                              n_cores = 3L,
                                              base_dir = "hvmf_realdata_screening",
                                              dataset_files = polar_dataset_files(),
                                              day_stride = 5L) {
  directories <- ensure_polar_screening_directories(base_dir = base_dir)
  summary_rows <- list()
  seed_base <- 2026052500L

  for (dataset_index in seq_along(dataset_files)) {
    component_csv <- dataset_files[[dataset_index]]
    cat(sprintf("Processing polar dataset %s\n", component_csv))
    summary_rows[[length(summary_rows) + 1L]] <- run_single_polar_dataset(
      component_csv = component_csv,
      directories = directories,
      B = B,
      n_cores = n_cores,
      seed = seed_base + dataset_index,
      day_stride = day_stride
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- order_screening_summary(summary_df)
  summary_csv <- file.path(directories$results_polar, "hvmf_realdata_screening_summary_polar.csv")
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  comparison_csv <- file.path(directories$results_polar, "hvmf_realdata_screening_component_vs_polar.csv")
  comparison_df <- make_component_vs_polar_comparison(
    component_summary_csv = file.path(directories$results, "hvmf_realdata_screening_summary.csv"),
    polar_summary_df = summary_df,
    output_csv = comparison_csv
  )

  list(
    summary = summary_df,
    summary_csv = summary_csv,
    comparison = comparison_df,
    comparison_csv = comparison_csv,
    directories = directories
  )
}

if (sys.nframe() == 0L) {
  run_hvmf_realdata_screening_polar()
}
