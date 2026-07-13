resolve_hvmf_real_data_path <- function(...) {
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

bootstrap_script_path <- resolve_hvmf_real_data_path("bootstrap", "multiplier_bootstrap.R")
if (!exists("multiplier_bootstrap_hvmf", mode = "function")) {
  source(bootstrap_script_path)
}

coerce_flag_values <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(logical(0))
  }
  if (is.logical(x)) {
    output <- x
    output[is.na(output)] <- default
    return(output)
  }
  if (is.numeric(x)) {
    output <- x != 0
    output[is.na(output)] <- default
    return(output)
  }

  values <- trimws(tolower(as.character(x)))
  output <- rep_len(default, length(values))
  true_values <- values %in% c("true", "t", "1", "yes", "y")
  false_values <- values %in% c("false", "f", "0", "no", "n", "", "na")
  unknown <- !(true_values | false_values)

  if (any(unknown)) {
    stop(sprintf("Could not coerce flag values: %s", paste(unique(values[unknown]), collapse = ", ")))
  }

  output[true_values] <- TRUE
  output[false_values] <- FALSE
  output
}

log_line <- function(con, ..., .console = TRUE) {
  line <- paste0(...)
  writeLines(line, con = con)
  flush(con)
  if (.console) {
    cat(line, "\n", sep = "")
    flush.console()
  }
}

format_timestamp <- function(x) {
  format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC", usetz = TRUE)
}

build_h2_from_speed_direction <- function(speed_scaled, direction_deg) {
  u <- as.numeric(speed_scaled)
  w <- as.numeric(direction_deg) * pi / 180

  x0 <- cosh(u)
  x1 <- sinh(u) * cos(w)
  x2 <- sinh(u) * sin(w)

  data.frame(
    x0 = x0,
    x1 = x1,
    x2 = x2,
    stringsAsFactors = FALSE
  )
}

verify_h2_matrix <- function(data_matrix, tol = 1e-8) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "double"

  minkowski_norm <- -data_matrix[, 1L]^2 + data_matrix[, 2L]^2 + data_matrix[, 3L]^2
  max_error <- max(abs(minkowski_norm + 1))

  if (!all(is.finite(data_matrix))) {
    stop("Data matrix contains non-finite entries.")
  }
  if (any(data_matrix[, 1L] <= 0)) {
    stop("Data matrix contains rows with x0 <= 0.")
  }
  if (max_error >= tol) {
    stop(sprintf("Data matrix is not on H^2 within tolerance %.1e; max error is %.10e.", tol, max_error))
  }

  normalize_hvmf_h2_data(data_matrix, tol = tol)

  list(
    data_matrix = data_matrix,
    minkowski_norm = minkowski_norm,
    max_minkowski_error = max_error
  )
}

prepare_jensen_hvmf_dataset <- function(path, tol = 1e-8) {
  if (!file.exists(path)) {
    stop(sprintf("Jensen dataset not found: %s", path))
  }

  raw_df <- read.csv(path, stringsAsFactors = FALSE)

  excluded_count <- 0L
  if ("excluded_manual" %in% names(raw_df)) {
    excluded <- coerce_flag_values(raw_df$excluded_manual, default = FALSE)
    excluded_count <- sum(excluded)
    raw_df <- raw_df[!excluded, , drop = FALSE]
  }

  if (all(c("x0", "x1", "x2") %in% names(raw_df))) {
    coord_df <- raw_df[, c("x0", "x1", "x2"), drop = FALSE]
  } else if (all(c("speed_scaled", "direction_deg") %in% names(raw_df))) {
    coord_df <- build_h2_from_speed_direction(
      speed_scaled = raw_df$speed_scaled,
      direction_deg = raw_df$direction_deg
    )
  } else {
    stop("Jensen dataset must contain either `x0`, `x1`, `x2` or `speed_scaled`, `direction_deg`.")
  }

  verification <- verify_h2_matrix(coord_df, tol = tol)
  user_marked_added <- if ("user_marked_added" %in% names(raw_df)) {
    sum(coerce_flag_values(raw_df$user_marked_added, default = FALSE))
  } else {
    0L
  }

  list(
    raw_df = raw_df,
    prepared_df = cbind(raw_df, as.data.frame(verification$data_matrix), minkowski_norm = verification$minkowski_norm),
    data_matrix = verification$data_matrix,
    n = nrow(verification$data_matrix),
    max_minkowski_error = verification$max_minkowski_error,
    notes = sprintf("included %d manually added points; excluded_manual removed %d rows", user_marked_added, excluded_count)
  )
}

prepare_existing_h2_csv_dataset <- function(path, tol = 1e-8) {
  if (!file.exists(path)) {
    stop(sprintf("Dataset not found: %s", path))
  }

  raw_df <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("x0", "x1", "x2")
  missing_cols <- setdiff(required_cols, names(raw_df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Dataset %s is missing columns: %s", path, paste(missing_cols, collapse = ", ")))
  }

  verification <- verify_h2_matrix(raw_df[, required_cols, drop = FALSE], tol = tol)

  list(
    raw_df = raw_df,
    prepared_df = cbind(raw_df, minkowski_norm = verification$minkowski_norm),
    data_matrix = verification$data_matrix,
    n = nrow(verification$data_matrix),
    max_minkowski_error = verification$max_minkowski_error,
    notes = "used existing x0, x1, x2 columns"
  )
}

prepare_real_hvmf_dataset <- function(dataset_type, path, tol = 1e-8) {
  switch(
    dataset_type,
    jensen = prepare_jensen_hvmf_dataset(path, tol = tol),
    modern = prepare_existing_h2_csv_dataset(path, tol = tol),
    stop(sprintf("Unsupported dataset_type: %s", dataset_type))
  )
}

fit_initial_hvmf_theta <- function(data_matrix) {
  fit_hvmf_theta(
    data = data_matrix,
    weights = NULL,
    null = list(type = "composite"),
    unknown_param = "both",
    control = list()
  )
}

make_real_data_analysis_configs <- function(results_dir = "wind/results") {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    list(
      dataset = "jensen_fig4",
      dataset_type = "jensen",
      input_csv = "wind/jensen_fig4_reconstructed_dataset1_v1_user_marked.csv",
      output_rds = file.path(results_dir, "cvm_hvmf_jensen_fig4_result.rds"),
      log_txt = file.path(results_dir, "log_cvm_hvmf_jensen_fig4.txt"),
      seed = 20260524L
    ),
    list(
      dataset = "risoe_modern_set_A_77m",
      dataset_type = "modern",
      input_csv = "wind/processed/risoe_modern_set_A_77m_hvmf.csv",
      output_rds = file.path(results_dir, "cvm_hvmf_risoe_modern_set_A_77m_result.rds"),
      log_txt = file.path(results_dir, "log_cvm_hvmf_risoe_modern_set_A_77m.txt"),
      seed = 20260525L
    ),
    list(
      dataset = "risoe_modern_set_B_125m",
      dataset_type = "modern",
      input_csv = "wind/processed/risoe_modern_set_B_125m_hvmf.csv",
      output_rds = file.path(results_dir, "cvm_hvmf_risoe_modern_set_B_125m_result.rds"),
      log_txt = file.path(results_dir, "log_cvm_hvmf_risoe_modern_set_B_125m.txt"),
      seed = 20260526L
    )
  )
}

run_single_hvmf_real_data_analysis <- function(config,
                                               B = 5000L,
                                               n_cores = 3L,
                                               alpha = 0.05,
                                               tol = 1e-8) {
  start_time <- Sys.time()
  log_con <- file(config$log_txt, open = "wt")
  on.exit(close(log_con), add = TRUE)

  warning_messages <- character()

  log_line(log_con, "Dataset: ", config$dataset)
  log_line(log_con, "Start timestamp: ", format_timestamp(start_time))
  log_line(log_con, "Input dataset: ", config$input_csv)
  log_line(log_con, "Output RDS: ", config$output_rds)
  log_line(log_con, "Configuration: statistic=CvM, null=composite, unknown_param=both, B=", B, ", n_cores=", n_cores, ", alpha=", alpha, ", seed=", config$seed)

  result_bundle <- tryCatch(
    withCallingHandlers({
      prepared <- prepare_real_hvmf_dataset(
        dataset_type = config$dataset_type,
        path = config$input_csv,
        tol = tol
      )
      data_matrix <- prepared$data_matrix
      initial_theta <- fit_initial_hvmf_theta(data_matrix)

      log_line(log_con, "Observations: ", prepared$n)
      log_line(log_con, "Max |<x,x>_M + 1|: ", sprintf("%.10e", prepared$max_minkowski_error))
      log_line(log_con, "Dataset notes: ", prepared$notes)
      log_line(
        log_con,
        "Initial MLE: kappa_hat=", sprintf("%.10f", initial_theta$kappa),
        ", mu_hat=(", sprintf("%.10f", initial_theta$mu[[1L]]), ", ",
        sprintf("%.10f", initial_theta$mu[[2L]]), ", ",
        sprintf("%.10f", initial_theta$mu[[3L]]), "), ",
        "sinh_chi_hat=", sprintf("%.10f", initial_theta$sinh_chi),
        ", theta_deg_hat=", sprintf("%.10f", initial_theta$theta_deg)
      )

      bootstrap_result <- multiplier_bootstrap_hvmf(
        data = data_matrix,
        null = list(type = "composite"),
        statistics = "cvm",
        B = B,
        alpha = alpha,
        n_cores = n_cores,
        seed = config$seed,
        keep = list(
          observed_process = FALSE,
          bootstrap_statistics = TRUE
        ),
        unknown_param = "both"
      )

      end_time <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
      theta_hat <- bootstrap_result$observed$theta_hat

      bundle <- list(
        dataset = config$dataset,
        dataset_type = config$dataset_type,
        input_csv = config$input_csv,
        output_rds = config$output_rds,
        log_txt = config$log_txt,
        start_time = start_time,
        end_time = end_time,
        elapsed_seconds = elapsed_seconds,
        seed = config$seed,
        B = B,
        n_cores = n_cores,
        alpha = alpha,
        max_minkowski_error = prepared$max_minkowski_error,
        n = prepared$n,
        initial_theta = initial_theta,
        warnings = warning_messages,
        prepared_notes = prepared$notes,
        result = bootstrap_result,
        summary = data.frame(
          dataset = config$dataset,
          n = prepared$n,
          kappa_hat = theta_hat$kappa,
          mu_hat_x0 = theta_hat$mu[[1L]],
          mu_hat_x1 = theta_hat$mu[[2L]],
          mu_hat_x2 = theta_hat$mu[[3L]],
          sinh_chi_hat = theta_hat$sinh_chi,
          theta_deg_hat = theta_hat$theta_deg,
          statistic_CvM = bootstrap_result$observed$cvm$statistic,
          p_value = bootstrap_result$inference$cvm$p_value,
          B = B,
          n_cores = n_cores,
          elapsed_seconds = elapsed_seconds,
          output_rds = config$output_rds,
          log_txt = config$log_txt,
          input_csv = config$input_csv,
          max_minkowski_error = prepared$max_minkowski_error,
          seed = config$seed,
          stringsAsFactors = FALSE
        )
      )

      saveRDS(bundle, file = config$output_rds)

      log_line(log_con, "Observed CvM statistic: ", sprintf("%.10f", bootstrap_result$observed$cvm$statistic))
      log_line(log_con, "Bootstrap p-value: ", sprintf("%.10f", bootstrap_result$inference$cvm$p_value))
      log_line(log_con, "Finish timestamp: ", format_timestamp(end_time))
      log_line(log_con, "Elapsed seconds: ", sprintf("%.3f", elapsed_seconds))
      if (length(warning_messages) == 0L) {
        log_line(log_con, "Warnings: none")
      } else {
        log_line(log_con, "Warnings:")
        for (warning_message in warning_messages) {
          log_line(log_con, "  - ", warning_message)
        }
      }

      bundle
    }, warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      end_time <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
      log_line(log_con, "Finish timestamp: ", format_timestamp(end_time))
      log_line(log_con, "Elapsed seconds: ", sprintf("%.3f", elapsed_seconds))
      log_line(log_con, "ERROR: ", conditionMessage(e))
      if (length(warning_messages) > 0L) {
        log_line(log_con, "Warnings seen before failure:")
        for (warning_message in warning_messages) {
          log_line(log_con, "  - ", warning_message)
        }
      }
      stop(e)
    }
  )

  result_bundle
}

write_real_data_summary_csv <- function(summary_rows,
                                        output_csv = "wind/results/cvm_hvmf_real_data_summary.csv") {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file = output_csv, row.names = FALSE)
  output_csv
}

run_hvmf_real_data_cvm_analyses <- function(B = 5000L,
                                            n_cores = 3L,
                                            results_dir = "wind/results",
                                            alpha = 0.05,
                                            tol = 1e-8) {
  configs <- make_real_data_analysis_configs(results_dir = results_dir)
  summary_rows <- vector("list", length(configs))

  for (i in seq_along(configs)) {
    config <- configs[[i]]
    cat("\n=== Running ", config$dataset, " (", i, "/", length(configs), ") ===\n", sep = "")
    bundle <- run_single_hvmf_real_data_analysis(
      config = config,
      B = B,
      n_cores = n_cores,
      alpha = alpha,
      tol = tol
    )
    summary_rows[[i]] <- bundle$summary
  }

  summary_csv <- write_real_data_summary_csv(
    summary_rows = summary_rows,
    output_csv = file.path(results_dir, "cvm_hvmf_real_data_summary.csv")
  )

  invisible(list(
    configs = configs,
    summary = do.call(rbind, summary_rows),
    summary_csv = summary_csv
  ))
}

if (sys.nframe() == 0L) {
  run_hvmf_real_data_cvm_analyses()
}
