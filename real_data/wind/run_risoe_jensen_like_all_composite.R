source(file.path("real_data", "wind", "run_hvmf_real_data_cvm.R"))

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

format_year_vector <- function(year_start, year_end, window_id = NULL) {
  if (!is.na(year_start) && !is.na(year_end) && grepl("excl_", window_id %||% "", fixed = TRUE)) {
    excluded <- sub("^.*excl_", "", window_id)
    return(sprintf("%d:%d excluding %s", year_start, year_end, gsub("_", ",", excluded, fixed = TRUE)))
  }
  sprintf("%d:%d", year_start, year_end)
}

safe_quantile <- function(x, probs) {
  if (length(x) == 0L || all(!is.finite(x))) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x[is.finite(x)], probs = probs, names = FALSE, type = 8))
}

extract_theta_star_diagnostics <- function(theta_star) {
  if (is.null(theta_star) || length(theta_star) == 0L) {
    return(list(
      kappa_values = numeric(),
      mu_norm_errors = numeric(),
      n_failed = NA_integer_
    ))
  }

  kappa_values <- vapply(theta_star, function(theta) {
    if (is.null(theta) || is.null(theta$kappa)) {
      return(NA_real_)
    }
    as.numeric(theta$kappa)
  }, numeric(1))

  mu_norm_errors <- vapply(theta_star, function(theta) {
    if (is.null(theta) || is.null(theta$mu)) {
      return(NA_real_)
    }
    mu <- as.numeric(theta$mu)
    abs(-mu[[1L]]^2 + mu[[2L]]^2 + mu[[3L]]^2 + 1)
  }, numeric(1))

  n_failed <- sum(!is.finite(kappa_values) | !is.finite(mu_norm_errors))

  list(
    kappa_values = kappa_values,
    mu_norm_errors = mu_norm_errors,
    n_failed = as.integer(n_failed)
  )
}

compute_bootstrap_cvm_summary <- function(bootstrap_values) {
  qs <- safe_quantile(bootstrap_values, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
  list(
    min = if (length(bootstrap_values) > 0L) min(bootstrap_values) else NA_real_,
    q01 = qs[[1L]],
    q05 = qs[[2L]],
    q25 = qs[[3L]],
    median = qs[[4L]],
    mean = if (length(bootstrap_values) > 0L) mean(bootstrap_values) else NA_real_,
    q75 = qs[[5L]],
    q95 = qs[[6L]],
    q99 = qs[[7L]],
    max = if (length(bootstrap_values) > 0L) max(bootstrap_values) else NA_real_
  )
}

stop_if_cvm_observed_mismatch <- function(dataset_id, cvm_simple, cvm_composite, tol = 1e-8) {
  if (!is.finite(cvm_simple) || !is.finite(cvm_composite)) {
    stop(sprintf("Observed CvM comparison failed for %s because one value is non-finite.", dataset_id))
  }

  if (abs(cvm_simple - cvm_composite) > tol) {
    stop(sprintf(
      "Observed CvM mismatch for %s: simple = %.10f, composite = %.10f, |diff| = %.3e > %.1e.",
      dataset_id,
      cvm_simple,
      cvm_composite,
      abs(cvm_simple - cvm_composite),
      tol
    ))
  }
}

composite_case_log_header <- function(log_con, config_row, B, n_cores, profile_method) {
  log_line(log_con, "Dataset: ", config_row$dataset_id)
  log_line(log_con, "Source file: ", config_row$source_file)
  log_line(log_con, "Processed CSV: ", config_row$processed_csv)
  log_line(log_con, "Window ID: ", config_row$window_id)
  log_line(log_con, "Years: ", format_year_vector(config_row$year_start, config_row$year_end, config_row$window_id))
  log_line(log_con, "Height: ", config_row$height_m, " m")
  log_line(log_con, "Pattern: ", config_row$pattern)
  log_line(log_con, "n_final: ", config_row$n_final)
  log_line(log_con, "speed_mean: ", sprintf("%.10f", config_row$speed_mean))
  log_line(log_con, "Simple plug-in reference: CvM = ", sprintf("%.10f", config_row$statistic_CvM), ", p = ", sprintf("%.10f", config_row$p_value))
  log_line(log_con, "Composite configuration: statistic=CvM, null=composite, unknown_param=both, B=", B, ", n_cores=", n_cores, ", hvmf_profile_method=", profile_method)
}

run_single_risoe_all_composite_case <- function(config_row,
                                                output_dir,
                                                B = 5000L,
                                                n_cores = 4L,
                                                alpha = 0.05,
                                                profile_method = "tabulated",
                                                tol = 1e-8,
                                                seed = 74000L) {
  results_dir <- file.path(output_dir, "results")
  logs_dir <- file.path(output_dir, "logs")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  result_rds <- file.path(results_dir, paste0(config_row$dataset_id, "_cvm_composite_result.rds"))
  log_file <- file.path(logs_dir, paste0(config_row$dataset_id, "_composite_log.txt"))

  start_time <- Sys.time()
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  warning_messages <- character()
  status <- "ok"
  error_message <- NA_character_
  summary_row <- NULL

  log_line(log_con, "Start timestamp: ", format_timestamp(start_time))
  composite_case_log_header(log_con, config_row, B = B, n_cores = n_cores, profile_method = profile_method)

  tryCatch(
    {
      prepared <- prepare_existing_h2_csv_dataset(config_row$processed_csv, tol = tol)
      X <- prepared$data_matrix
      minkowski_error <- prepared$max_minkowski_error
      stopifnot(all(X[, 1L] > 0))
      stopifnot(minkowski_error < tol)

      fit <- hvmf_mle_h2(X)
      log_line(log_con, "MLE on observed sample:")
      log_line(log_con, "  kappa_hat = ", sprintf("%.10f", fit$kappa))
      log_line(log_con, "  mu_hat = (", paste(sprintf("%.10f", fit$mu), collapse = ", "), ")")
      log_line(log_con, "  sinh_chi_hat = ", sprintf("%.10f", fit$sinh_chi))
      log_line(log_con, "  theta_deg_hat = ", sprintf("%.10f", fit$theta_deg))
      log_line(log_con, "  max abs(minkowski_norm + 1) = ", sprintf("%.10e", minkowski_error))

      bootstrap_result <- withCallingHandlers(
        multiplier_bootstrap_hvmf(
          data = X,
          null = list(type = "composite"),
          statistics = "cvm",
          B = B,
          alpha = alpha,
          n_cores = n_cores,
          seed = seed,
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

      cvm_composite <- bootstrap_result$observed$cvm$statistic
      stop_if_cvm_observed_mismatch(
        dataset_id = config_row$dataset_id,
        cvm_simple = config_row$statistic_CvM,
        cvm_composite = cvm_composite,
        tol = 1e-8
      )

      bootstrap_cvm <- bootstrap_result$bootstrap$statistics$cvm
      n_ge_observed <- sum(bootstrap_cvm >= cvm_composite)
      cvm_summary <- compute_bootstrap_cvm_summary(bootstrap_cvm)
      theta_diag <- extract_theta_star_diagnostics(bootstrap_result$bootstrap$theta_star)
      kappa_q <- safe_quantile(theta_diag$kappa_values, probs = c(0.05, 0.5, 0.95))

      elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      saveRDS(
        list(
          config_row = config_row,
          fit = fit,
          bootstrap_result = bootstrap_result,
          cvm_summary = cvm_summary,
          theta_star_diagnostics = theta_diag
        ),
        file = result_rds
      )

      log_line(log_con, "Observed CvM simple plug-in = ", sprintf("%.10f", config_row$statistic_CvM))
      log_line(log_con, "Observed CvM composite   = ", sprintf("%.10f", cvm_composite))
      log_line(log_con, "Observed CvM diff        = ", sprintf("%.3e", abs(config_row$statistic_CvM - cvm_composite)))
      log_line(log_con, "Bootstrap CvM summary:")
      log_line(log_con, "  min = ", sprintf("%.10f", cvm_summary$min))
      log_line(log_con, "  q01 = ", sprintf("%.10f", cvm_summary$q01))
      log_line(log_con, "  q05 = ", sprintf("%.10f", cvm_summary$q05))
      log_line(log_con, "  q25 = ", sprintf("%.10f", cvm_summary$q25))
      log_line(log_con, "  median = ", sprintf("%.10f", cvm_summary$median))
      log_line(log_con, "  mean = ", sprintf("%.10f", cvm_summary$mean))
      log_line(log_con, "  q75 = ", sprintf("%.10f", cvm_summary$q75))
      log_line(log_con, "  q95 = ", sprintf("%.10f", cvm_summary$q95))
      log_line(log_con, "  q99 = ", sprintf("%.10f", cvm_summary$q99))
      log_line(log_con, "  max = ", sprintf("%.10f", cvm_summary$max))
      log_line(log_con, "n_bootstrap_ge_observed = ", n_ge_observed)
      log_line(log_con, "p_value_composite = ", sprintf("%.10f", bootstrap_result$inference$cvm$p_value))
      log_line(log_con, "Bootstrap kappa_star summary:")
      log_line(log_con, "  min = ", sprintf("%.10f", min(theta_diag$kappa_values, na.rm = TRUE)))
      log_line(log_con, "  q05 = ", sprintf("%.10f", kappa_q[[1L]]))
      log_line(log_con, "  median = ", sprintf("%.10f", kappa_q[[2L]]))
      log_line(log_con, "  q95 = ", sprintf("%.10f", kappa_q[[3L]]))
      log_line(log_con, "  max = ", sprintf("%.10f", max(theta_diag$kappa_values, na.rm = TRUE)))
      log_line(log_con, "mu_star Minkowski norm error:")
      log_line(log_con, "  min = ", sprintf("%.10e", min(theta_diag$mu_norm_errors, na.rm = TRUE)))
      log_line(log_con, "  max = ", sprintf("%.10e", max(theta_diag$mu_norm_errors, na.rm = TRUE)))
      log_line(log_con, "n_failed_bootstrap_fits = ", theta_diag$n_failed)
      log_line(log_con, "elapsed_seconds = ", sprintf("%.3f", elapsed_seconds))

      summary_row <- data.frame(
        dataset_id = config_row$dataset_id,
        height_m = config_row$height_m,
        window_id = config_row$window_id,
        year_start = config_row$year_start,
        year_end = config_row$year_end,
        years = format_year_vector(config_row$year_start, config_row$year_end, config_row$window_id),
        pattern = config_row$pattern,
        n_final = config_row$n_final,
        speed_mean = config_row$speed_mean,
        kappa_hat = fit$kappa,
        mu_hat_x0 = fit$mu[[1L]],
        mu_hat_x1 = fit$mu[[2L]],
        mu_hat_x2 = fit$mu[[3L]],
        sinh_chi_hat = fit$sinh_chi,
        theta_deg_hat = fit$theta_deg,
        CvM_simple_plugin = config_row$statistic_CvM,
        p_value_simple_plugin = config_row$p_value,
        CvM_composite = cvm_composite,
        p_value_composite = bootstrap_result$inference$cvm$p_value,
        B = B,
        n_cores = n_cores,
        bootstrap_cvm_min = cvm_summary$min,
        bootstrap_cvm_q05 = cvm_summary$q05,
        bootstrap_cvm_median = cvm_summary$median,
        bootstrap_cvm_q95 = cvm_summary$q95,
        bootstrap_cvm_max = cvm_summary$max,
        n_bootstrap_ge_observed = n_ge_observed,
        kappa_star_min = min(theta_diag$kappa_values, na.rm = TRUE),
        kappa_star_q05 = kappa_q[[1L]],
        kappa_star_median = kappa_q[[2L]],
        kappa_star_q95 = kappa_q[[3L]],
        kappa_star_max = max(theta_diag$kappa_values, na.rm = TRUE),
        n_failed_bootstrap_fits = theta_diag$n_failed,
        elapsed_seconds = elapsed_seconds,
        processed_csv = config_row$processed_csv,
        result_rds = result_rds,
        log_file = log_file,
        status = status,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
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
  log_line(log_con, "End timestamp: ", format_timestamp(end_time))
  log_line(log_con, "Total elapsed seconds: ", sprintf("%.3f", as.numeric(difftime(end_time, start_time, units = "secs"))))

  if (is.null(summary_row)) {
    summary_row <- data.frame(
      dataset_id = config_row$dataset_id,
      height_m = config_row$height_m,
      window_id = config_row$window_id,
      year_start = config_row$year_start,
      year_end = config_row$year_end,
      years = format_year_vector(config_row$year_start, config_row$year_end, config_row$window_id),
      pattern = config_row$pattern,
      n_final = config_row$n_final,
      speed_mean = config_row$speed_mean,
      kappa_hat = NA_real_,
      mu_hat_x0 = NA_real_,
      mu_hat_x1 = NA_real_,
      mu_hat_x2 = NA_real_,
      sinh_chi_hat = NA_real_,
      theta_deg_hat = NA_real_,
      CvM_simple_plugin = config_row$statistic_CvM,
      p_value_simple_plugin = config_row$p_value,
      CvM_composite = NA_real_,
      p_value_composite = NA_real_,
      B = B,
      n_cores = n_cores,
      bootstrap_cvm_min = NA_real_,
      bootstrap_cvm_q05 = NA_real_,
      bootstrap_cvm_median = NA_real_,
      bootstrap_cvm_q95 = NA_real_,
      bootstrap_cvm_max = NA_real_,
      n_bootstrap_ge_observed = NA_integer_,
      kappa_star_min = NA_real_,
      kappa_star_q05 = NA_real_,
      kappa_star_median = NA_real_,
      kappa_star_q95 = NA_real_,
      kappa_star_max = NA_real_,
      n_failed_bootstrap_fits = NA_integer_,
      elapsed_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
      processed_csv = config_row$processed_csv,
      result_rds = result_rds,
      log_file = log_file,
      status = status,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }

  summary_row
}

compute_composite_cvm_observed <- function(X, theta, profile_method = "tabulated") {
  spec <- make_hvmf_spec(unknown_param = "both")
  prep <- prepare_cvm_observed_data(
    data = X,
    spec = spec,
    theta_hat = theta,
    control = list(hvmf_profile_method = profile_method)
  )
  prep$statistic
}

run_hvmf_parametric_composite_check <- function(config_row,
                                                output_dir,
                                                B_param = 1000L,
                                                n_cores = 4L,
                                                seed = 88000L,
                                                profile_method = "tabulated",
                                                tol = 1e-8) {
  parametric_dir <- file.path(output_dir, "parametric_check")
  dir.create(parametric_dir, recursive = TRUE, showWarnings = FALSE)

  result_rds <- file.path(parametric_dir, paste0(config_row$dataset_id, "_parametric_composite_check.rds"))
  log_file <- file.path(parametric_dir, paste0(config_row$dataset_id, "_parametric_composite_check.txt"))

  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)
  start_time <- Sys.time()
  log_line(log_con, "Start timestamp: ", format_timestamp(start_time))
  log_line(log_con, "Dataset: ", config_row$dataset_id)
  log_line(log_con, "Processed CSV: ", config_row$processed_csv)
  log_line(log_con, "Configuration: parametric composite check, B_param=", B_param, ", n_cores=", n_cores, ", profile_method=", profile_method)

  prepared <- prepare_existing_h2_csv_dataset(config_row$processed_csv, tol = tol)
  X <- prepared$data_matrix
  fit <- hvmf_mle_h2(X)
  theta_hat <- list(mu = fit$mu, kappa = fit$kappa)
  observed_cvm <- compute_composite_cvm_observed(X, theta = theta_hat, profile_method = profile_method)

  set.seed(seed)
  worker_fun <- function(b) {
    set.seed(seed + b)
    Xb <- rhvmf_h2_polar(n = nrow(X), mu = fit$mu, kappa = fit$kappa, check = TRUE)
    fit_b <- hvmf_mle_h2(Xb)
    theta_b <- list(mu = fit_b$mu, kappa = fit_b$kappa)
    cvm_b <- compute_composite_cvm_observed(Xb, theta = theta_b, profile_method = profile_method)
    list(
      cvm = cvm_b,
      kappa = fit_b$kappa,
      mu_norm_error = abs(-fit_b$mu[[1L]]^2 + fit_b$mu[[2L]]^2 + fit_b$mu[[3L]]^2 + 1)
    )
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parametric_results <- parallel::mclapply(seq_len(B_param), worker_fun, mc.cores = n_cores)
  } else {
    parametric_results <- lapply(seq_len(B_param), worker_fun)
  }

  cvm_star <- vapply(parametric_results, `[[`, numeric(1), "cvm")
  kappa_star <- vapply(parametric_results, `[[`, numeric(1), "kappa")
  mu_error_star <- vapply(parametric_results, `[[`, numeric(1), "mu_norm_error")
  p_value_parametric <- (1 + sum(cvm_star >= observed_cvm)) / (B_param + 1)

  output <- list(
    dataset_id = config_row$dataset_id,
    observed_cvm = observed_cvm,
    p_value_parametric = p_value_parametric,
    cvm_star = cvm_star,
    kappa_star = kappa_star,
    mu_error_star = mu_error_star,
    fit = fit,
    B_param = B_param,
    n_cores = n_cores,
    processed_csv = config_row$processed_csv
  )
  saveRDS(output, result_rds)

  log_line(log_con, "Observed CvM = ", sprintf("%.10f", observed_cvm))
  log_line(log_con, "Parametric p-value = ", sprintf("%.10f", p_value_parametric))
  log_line(log_con, "CvM star summary:")
  cvm_summary <- compute_bootstrap_cvm_summary(cvm_star)
  log_line(log_con, "  q05 = ", sprintf("%.10f", cvm_summary$q05))
  log_line(log_con, "  median = ", sprintf("%.10f", cvm_summary$median))
  log_line(log_con, "  q95 = ", sprintf("%.10f", cvm_summary$q95))
  log_line(log_con, "kappa star summary:")
  kappa_q <- safe_quantile(kappa_star, c(0.05, 0.5, 0.95))
  log_line(log_con, "  q05 = ", sprintf("%.10f", kappa_q[[1L]]))
  log_line(log_con, "  median = ", sprintf("%.10f", kappa_q[[2L]]))
  log_line(log_con, "  q95 = ", sprintf("%.10f", kappa_q[[3L]]))
  log_line(log_con, "mu norm error max = ", sprintf("%.10e", max(mu_error_star)))
  log_line(log_con, "End timestamp: ", format_timestamp(Sys.time()))

  output
}

run_risoe_jensen_like_all_composite <- function(simple_summary_csv = file.path("wind", "jensen_like_all", "risoe_jensen_like_all_summary.csv"),
                                                output_dir = file.path("wind", "jensen_like_all_composite"),
                                                B = 5000L,
                                                n_cores = 4L,
                                                profile_method = "tabulated",
                                                alpha = 0.05,
                                                run_parametric_if_needed = TRUE) {
  if (!file.exists(simple_summary_csv)) {
    stop(sprintf("Simple summary CSV not found: %s", simple_summary_csv))
  }

  simple_summary <- read.csv(simple_summary_csv, stringsAsFactors = FALSE)
  simple_summary <- simple_summary[simple_summary$status == "ok", , drop = FALSE]
  rownames(simple_summary) <- NULL

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_rows <- vector("list", nrow(simple_summary))

  for (i in seq_len(nrow(simple_summary))) {
    config_row <- simple_summary[i, , drop = FALSE]
    cat("\n=== Composite HvMF-CvM ", i, "/", nrow(simple_summary), ": ", config_row$dataset_id, " ===\n", sep = "")
    summary_rows[[i]] <- run_single_risoe_all_composite_case(
      config_row = config_row,
      output_dir = output_dir,
      B = B,
      n_cores = n_cores,
      alpha = alpha,
      profile_method = profile_method,
      seed = 74000L + i
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  safe_p <- ifelse(is.finite(summary_df$p_value_composite), summary_df$p_value_composite, -Inf)
  summary_df <- summary_df[order(summary_df$status != "ok", -safe_p, -summary_df$n_final), , drop = FALSE]
  rownames(summary_df) <- NULL

  summary_csv <- file.path(output_dir, "risoe_jensen_like_all_composite_summary.csv")
  write.csv(summary_df, summary_csv, row.names = FALSE)

  parametric_result <- NULL
  if (isTRUE(run_parametric_if_needed)) {
    subset_77 <- summary_df[summary_df$status == "ok" & summary_df$height_m == 77, , drop = FALSE]
    if (nrow(subset_77) > 0L && all(is.finite(subset_77$p_value_composite)) && all(subset_77$p_value_composite <= 0.01)) {
      target_id <- "risoe_all_77m_10y_complete_1997_2007_excl_2002_nov_dec_set12"
      target_row <- simple_summary[simple_summary$dataset_id == target_id, , drop = FALSE]
      if (nrow(target_row) == 1L) {
        parametric_result <- run_hvmf_parametric_composite_check(
          config_row = target_row,
          output_dir = output_dir,
          B_param = 1000L,
          n_cores = n_cores,
          profile_method = profile_method
        )
      }
    }
  }

  invisible(list(
    summary = summary_df,
    summary_csv = summary_csv,
    parametric_result = parametric_result
  ))
}

if (sys.nframe() == 0L) {
  run_risoe_jensen_like_all_composite()
}
