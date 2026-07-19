source(file.path("real_data", "wind", "preprocess_risoe_modern_hvmf.R"))
source(file.path("real_data", "wind", "run_risoe_jensen_like_all_composite.R"))

select_noon_all_months <- function(df, fixed_tz = "UTC") {
  df$date <- as.Date(df$datetime, tz = fixed_tz)
  noon_time <- as.POSIXct(sprintf("%s 12:00:00", df$date), tz = fixed_tz)
  df$abs_noon_diff <- abs(as.numeric(difftime(df$datetime, noon_time, units = "secs")))
  df$datetime_num <- as.numeric(df$datetime)

  ord <- order(df$date, df$abs_noon_diff, df$datetime_num)
  df_ord <- df[ord, , drop = FALSE]
  selected <- df_ord[!duplicated(df_ord$date), , drop = FALSE]
  selected <- selected[order(selected$datetime), , drop = FALSE]
  rownames(selected) <- NULL

  selected$abs_noon_diff <- NULL
  selected$datetime_num <- NULL
  selected
}

build_jensen_like_case <- function(selected_df,
                                   years,
                                   months,
                                   day_pattern,
                                   speed_col,
                                   direction_col,
                                   height_m) {
  keep <- selected_df$year %in% years &
    selected_df$month %in% months &
    selected_df$day %in% as.integer(day_pattern)

  rows <- selected_df[keep, , drop = FALSE]
  valid <- is.finite(rows[[speed_col]]) &
    is.finite(rows[[direction_col]]) &
    rows[[speed_col]] > 0
  rows <- rows[valid, , drop = FALSE]
  rownames(rows) <- NULL

  direction_deg <- ((as.numeric(rows[[direction_col]]) %% 360) + 360) %% 360
  speed <- as.numeric(rows[[speed_col]])
  speed_mean <- mean(speed)
  speed_scaled <- speed / speed_mean
  angle_rad <- direction_deg * pi / 180

  output <- data.frame(
    datetime = rows$datetime,
    year = rows$year,
    month = rows$month,
    day = rows$day,
    hour = rows$hour,
    minute = rows$minute,
    height_m = as.numeric(height_m),
    speed = speed,
    direction_deg = direction_deg,
    speed_mean = rep(speed_mean, length(speed)),
    speed_scaled = speed_scaled,
    angle_rad = angle_rad,
    x0 = cosh(speed_scaled),
    x1 = sinh(speed_scaled) * cos(angle_rad),
    x2 = sinh(speed_scaled) * sin(angle_rad),
    minkowski_norm = -cosh(speed_scaled)^2 + (sinh(speed_scaled) * cos(angle_rad))^2 + (sinh(speed_scaled) * sin(angle_rad))^2,
    stringsAsFactors = FALSE
  )

  stopifnot(all(abs(output$minkowski_norm + 1) < 1e-8))
  output
}

compute_direction_summary <- function(df) {
  split_key <- interaction(df$year, df$month, drop = TRUE, lex.order = TRUE)
  groups <- split(df, split_key)

  rows <- lapply(groups, function(g) {
    hvmf_fit <- hvmf_mle_h2(as.matrix(g[, c("x0", "x1", "x2")]))
    frechet_h2 <- fit_h2_frechet_mean(as.matrix(g[, c("x0", "x1", "x2")]))
    s1_intrinsic <- fit_s1_intrinsic_mean(g$angle_rad)

    data.frame(
      year = unique(g$year),
      month = unique(g$month),
      n = nrow(g),
      mean_speed_scaled = mean(g$speed_scaled),
      s1_intrinsic_mean_deg = s1_intrinsic$theta_deg,
      h2_frechet_theta_deg = frechet_h2$theta_deg,
      h2_frechet_chi = frechet_h2$chi,
      h2_frechet_sinh_chi = frechet_h2$sinh_chi,
      hvmf_kappa_hat = hvmf_fit$kappa,
      hvmf_theta_deg_hat = hvmf_fit$theta_deg,
      hvmf_sinh_chi_hat = hvmf_fit$sinh_chi,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, rows)
  summary_df[order(summary_df$year, summary_df$month), , drop = FALSE]
}

compute_month_aggregate_summary <- function(df) {
  groups <- split(df, df$month)
  rows <- lapply(groups, function(g) {
    hvmf_fit <- hvmf_mle_h2(as.matrix(g[, c("x0", "x1", "x2")]))
    frechet_h2 <- fit_h2_frechet_mean(as.matrix(g[, c("x0", "x1", "x2")]))
    s1_intrinsic <- fit_s1_intrinsic_mean(g$angle_rad)

    data.frame(
      month = unique(g$month),
      n = nrow(g),
      s1_intrinsic_mean_deg = s1_intrinsic$theta_deg,
      h2_frechet_theta_deg = frechet_h2$theta_deg,
      h2_frechet_chi = frechet_h2$chi,
      h2_frechet_sinh_chi = frechet_h2$sinh_chi,
      mean_speed_scaled = mean(g$speed_scaled),
      hvmf_kappa_hat = hvmf_fit$kappa,
      hvmf_theta_deg_hat = hvmf_fit$theta_deg,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, rows)
  summary_df[order(summary_df$month), , drop = FALSE]
}

plot_month_year_arrows <- function(summary_df, height_m, output_png) {
  months_to_plot <- c(10L, 11L, 12L)
  month_colors <- c("10" = "#c26a00", "11" = "#1f6f8b", "12" = "#8c1c13")
  month_labels <- c("10" = "Oct", "11" = "Nov", "12" = "Dec")
  years <- sort(unique(summary_df$year))

  grDevices::png(output_png, width = 1600, height = 1000, res = 160)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 5), mar = c(2.2, 2.2, 3.2, 1.2), oma = c(0, 0, 3, 0))

  theta <- seq(0, 2 * pi, length.out = 361L)
  max_sinh_chi <- max(summary_df$h2_frechet_sinh_chi, na.rm = TRUE)

  for (year in years) {
    panel <- summary_df[summary_df$year == year & summary_df$month %in% months_to_plot, , drop = FALSE]
    graphics::plot(
      cos(theta),
      sin(theta),
      type = "l",
      asp = 1,
      axes = FALSE,
      xlab = "",
      ylab = "",
      xlim = c(-1.15, 1.15),
      ylim = c(-1.15, 1.15),
      main = as.character(year),
      col = "#b9b9b9",
      lwd = 1.2
    )
    graphics::abline(h = 0, v = 0, col = "#efefef", lty = 3)

    if (nrow(panel) > 0L) {
      for (i in seq_len(nrow(panel))) {
        month_chr <- as.character(panel$month[[i]])
        angle <- panel$h2_frechet_theta_deg[[i]] * pi / 180
        length_r <- if (is.finite(max_sinh_chi) && max_sinh_chi > 0) {
          0.2 + 0.75 * panel$h2_frechet_sinh_chi[[i]] / max_sinh_chi
        } else {
          0.8
        }
        x1 <- length_r * cos(angle)
        y1 <- length_r * sin(angle)
        graphics::arrows(
          x0 = 0,
          y0 = 0,
          x1 = x1,
          y1 = y1,
          length = 0.08,
          lwd = 2.5,
          col = month_colors[[month_chr]]
        )
        graphics::points(x1, y1, pch = 16, cex = 0.9, col = month_colors[[month_chr]])
      }
    }
  }

  graphics::mtext(sprintf("Risoe %sm: monthly mean wind directions (arrow length = directional resultant)", height_m), side = 3, outer = TRUE, line = 1, cex = 1.1)
  graphics::mtext("Arrow direction from H^2 Fréchet mean; arrow length proportional to sinh(chi) of that mean", side = 3, outer = TRUE, line = -0.3, cex = 0.8)
  graphics::legend(
    "bottom",
    inset = -0.4,
    horiz = TRUE,
    legend = unname(month_labels),
    col = unname(month_colors),
    lwd = 2.5,
    bty = "n",
    xpd = TRUE
  )
}

h2_point_from_chi_theta <- function(chi, theta) {
  c(cosh(chi), sinh(chi) * cos(theta), sinh(chi) * sin(theta))
}

wrap_angle_to_pi <- function(theta) {
  atan2(sin(theta), cos(theta))
}

fit_s1_intrinsic_mean <- function(angle_rad) {
  angle_rad <- as.numeric(angle_rad)
  angle_rad <- angle_rad[is.finite(angle_rad)]
  if (length(angle_rad) == 0L) {
    stop("`angle_rad` must contain at least one finite value.")
  }

  objective <- function(theta) {
    delta <- wrap_angle_to_pi(angle_rad - theta[[1L]])
    mean(delta^2)
  }

  start <- atan2(mean(sin(angle_rad)), mean(cos(angle_rad)))
  fit <- stats::optim(par = start, fn = objective, method = "BFGS")
  theta_hat <- wrap_angle_to_pi(fit$par[[1L]])

  list(
    theta = theta_hat,
    theta_deg = (theta_hat * 180 / pi) %% 360,
    value = fit$value,
    convergence = fit$convergence
  )
}

fit_h2_frechet_mean <- function(X) {
  X <- as.matrix(X)
  start_fit <- hvmf_mle_h2(X)
  start_eta <- log(start_fit$chi + 1e-8)
  start_theta <- start_fit$theta

  objective <- function(par) {
    chi <- exp(par[[1L]])
    theta <- par[[2L]]
    mu <- h2_point_from_chi_theta(chi, theta)
    mean(apply(X, 1L, function(x) hyperbolic_geodesic_distance_h2(mu, x)^2))
  }

  fit <- stats::optim(
    par = c(start_eta, start_theta),
    fn = objective,
    method = "BFGS",
    control = list(reltol = 1e-12, maxit = 1000L)
  )

  chi_hat <- exp(fit$par[[1L]])
  theta_hat <- wrap_angle_to_pi(fit$par[[2L]])
  mu_hat <- h2_point_from_chi_theta(chi_hat, theta_hat)

  list(
    mu = mu_hat,
    chi = chi_hat,
    sinh_chi = sinh(chi_hat),
    theta = theta_hat,
    theta_deg = (theta_hat * 180 / pi) %% 360,
    value = fit$value,
    convergence = fit$convergence
  )
}

run_composite_benchmark_case <- function(df,
                                         dataset_id,
                                         output_dir,
                                         B = 5000L,
                                         n_cores = 4L,
                                         profile_method = "tabulated",
                                         seed = NULL,
                                         verbose = FALSE) {
  results_dir <- file.path(output_dir, "results")
  logs_dir <- file.path(output_dir, "logs")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  csv_path <- file.path(output_dir, paste0(dataset_id, "_hvmf.csv"))
  write.csv(df, csv_path, row.names = FALSE)

  X <- as.matrix(df[, c("x0", "x1", "x2")])
  fit <- hvmf_mle_h2(X)
  start_time <- Sys.time()
  bootstrap_seed <- if (is.null(seed)) {
    92000L + nrow(df)
  } else {
    as.integer(seed)
  }
  log_file <- file.path(logs_dir, sprintf("%s_composite_B%d_seed%d_log.txt", dataset_id, B, bootstrap_seed))
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)
  log_line <- function(...) {
    line <- paste0(...)
    if (isTRUE(verbose)) {
      cat(line, "\n")
    }
    cat(line, "\n", file = log_con)
  }

  warning_messages <- character()
  log_line("Dataset: ", dataset_id)
  log_line("Start timestamp: ", format_timestamp(start_time))
  log_line("n = ", nrow(df))
  log_line("Months = ", paste(sort(unique(df$month)), collapse = ","))
  log_line("Years = ", paste(sort(unique(df$year)), collapse = ","))
  log_line("Height = ", unique(df$height_m), " m")
  log_line("Speed mean = ", sprintf("%.10f", unique(df$speed_mean)))
  log_line("MLE kappa_hat = ", sprintf("%.10f", fit$kappa))
  log_line("MLE mu_hat = (", paste(sprintf("%.10f", fit$mu), collapse = ", "), ")")
  log_line("MLE theta_deg_hat = ", sprintf("%.10f", fit$theta_deg))
  log_line("Configuration: composite CvM, B=", B, ", n_cores=", n_cores, ", hvmf_profile_method=", profile_method)
  log_line("Bootstrap seed = ", bootstrap_seed)

  result <- withCallingHandlers(
    multiplier_bootstrap_hvmf(
      data = X,
      null = list(type = "composite"),
      statistics = "cvm",
      B = B,
      alpha = 0.05,
      n_cores = n_cores,
      seed = bootstrap_seed,
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
  saveRDS(
    list(
      dataset_id = dataset_id,
      fit = fit,
      result = result,
      seed = bootstrap_seed
    ),
    file = file.path(results_dir, sprintf("%s_cvm_composite_B%d_seed%d_result.rds", dataset_id, B, bootstrap_seed))
  )

  log_line("Observed CvM = ", sprintf("%.10f", result$observed$cvm$statistic))
  log_line("p-value = ", sprintf("%.10f", result$inference$cvm$p_value))
  log_line("Elapsed seconds = ", sprintf("%.3f", elapsed_seconds))
  if (length(warning_messages) == 0L) {
    log_line("Warnings: none")
  } else {
    log_line("Warnings:")
    for (msg in warning_messages) {
      log_line("  - ", msg)
    }
  }
  log_line("End timestamp: ", format_timestamp(Sys.time()))

  data.frame(
    dataset_id = dataset_id,
    height_m = unique(df$height_m),
    months = paste(sort(unique(df$month)), collapse = ","),
    n = nrow(df),
    speed_mean = unique(df$speed_mean),
    kappa_hat = fit$kappa,
    theta_deg_hat = fit$theta_deg,
    cvm_composite = result$observed$cvm$statistic,
    p_value_composite = result$inference$cvm$p_value,
    B = B,
    n_cores = n_cores,
    seed = bootstrap_seed,
    elapsed_seconds = elapsed_seconds,
    csv_path = csv_path,
    log_file = log_file,
    stringsAsFactors = FALSE
  )
}

run_risoe_oct_nov_dec_analysis <- function(input_nc = "real_data/wind/risoe_m_all.nc",
                                           output_dir = "real_data/wind/oct_nov_dec_analysis",
                                           years = c(1997:2001, 2003:2007),
                                           day_pattern = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
                                           B = 5000L,
                                           n_cores = 4L,
                                           profile_method = "tabulated",
                                           seed = 20260527L,
                                           verbose = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  all_df <- load_risoe_concurrent(input_nc)
  selected_df <- select_noon_all_months(all_df, fixed_tz = "UTC")

  height_specs <- list(
    list(height_m = 77L, speed_col = "ws77", direction_col = "wd77"),
    list(height_m = 125L, speed_col = "ws125", direction_col = "wd125")
  )

  benchmark_month_sets <- list(
    oct_nov = c(10L, 11L),
    nov_dec = c(11L, 12L),
    oct_dec = c(10L, 12L),
    oct_nov_dec = c(10L, 11L, 12L)
  )
  combined_benchmark_rows <- list()
  idx_combined <- 1L

  direction_rows <- list()
  month_rows <- list()
  benchmark_rows <- list()
  idx_dir <- 1L
  idx_month <- 1L
  idx_bench <- 1L

  for (spec in height_specs) {
    df_on_d <- build_jensen_like_case(
      selected_df = selected_df,
      years = years,
      months = c(10L, 11L, 12L),
      day_pattern = day_pattern,
      speed_col = spec$speed_col,
      direction_col = spec$direction_col,
      height_m = spec$height_m
    )

    direction_summary <- compute_direction_summary(df_on_d)
    direction_summary$height_m <- spec$height_m
    direction_rows[[idx_dir]] <- direction_summary
    idx_dir <- idx_dir + 1L

    month_summary <- compute_month_aggregate_summary(df_on_d)
    month_summary$height_m <- spec$height_m
    month_rows[[idx_month]] <- month_summary
    idx_month <- idx_month + 1L

    plot_month_year_arrows(
      summary_df = direction_summary,
      height_m = spec$height_m,
      output_png = file.path(output_dir, sprintf("risoe_%sm_oct_nov_dec_direction_arrows.png", spec$height_m))
    )

    for (benchmark_name in names(benchmark_month_sets)) {
      months <- benchmark_month_sets[[benchmark_name]]
      df_case <- build_jensen_like_case(
        selected_df = selected_df,
        years = years,
        months = months,
        day_pattern = day_pattern,
        speed_col = spec$speed_col,
        direction_col = spec$direction_col,
        height_m = spec$height_m
      )

      if (nrow(df_case) < 100L) {
        next
      }

      benchmark_rows[[idx_bench]] <- run_composite_benchmark_case(
        df = df_case,
        dataset_id = sprintf("risoe_%sm_%s_set12_1997_2007_excl_2002", spec$height_m, benchmark_name),
        output_dir = output_dir,
        B = B,
        n_cores = n_cores,
        profile_method = profile_method,
        seed = as.integer(seed + 100000L * spec$height_m + 1000L * min(months) + max(months)),
        verbose = verbose
      )
      idx_bench <- idx_bench + 1L
    }
  }

  combined_years <- c(1997:2001, 2003:2004)
  combined_77_df <- build_jensen_like_case(
    selected_df = selected_df,
    years = combined_years,
    months = c(10L, 11L, 12L),
    day_pattern = day_pattern,
    speed_col = "ws77",
    direction_col = "wd77",
    height_m = 77L
  )
  combined_125_df <- build_jensen_like_case(
    selected_df = selected_df,
    years = combined_years,
    months = c(10L, 11L, 12L),
    day_pattern = day_pattern,
    speed_col = "ws125",
    direction_col = "wd125",
    height_m = 125L
  )

  for (benchmark_name in names(benchmark_month_sets)) {
    months <- benchmark_month_sets[[benchmark_name]]
    combined_df <- rbind(
      combined_77_df[combined_77_df$month %in% months, , drop = FALSE],
      combined_125_df[combined_125_df$month %in% months, , drop = FALSE]
    )
    rownames(combined_df) <- NULL

    if (nrow(combined_df) < 100L) {
      next
    }

    combined_benchmark_rows[[idx_combined]] <- run_composite_benchmark_case(
      df = combined_df,
      dataset_id = sprintf("risoe_combined_77_125m_%s_set12_1997_2004_excl_2002", benchmark_name),
      output_dir = output_dir,
      B = B,
      n_cores = n_cores,
      profile_method = profile_method,
      seed = as.integer(seed + 9000000L + 1000L * min(months) + max(months)),
      verbose = verbose
    )
    idx_combined <- idx_combined + 1L
  }

  direction_summary_df <- do.call(rbind, direction_rows)
  direction_summary_csv <- file.path(output_dir, "direction_year_month_summary.csv")
  write.csv(direction_summary_df, direction_summary_csv, row.names = FALSE)

  month_summary_df <- do.call(rbind, month_rows)
  month_summary_csv <- file.path(output_dir, "direction_month_aggregate_summary.csv")
  write.csv(month_summary_df, month_summary_csv, row.names = FALSE)

  benchmark_summary_df <- do.call(rbind, benchmark_rows)
  benchmark_summary_df <- benchmark_summary_df[order(benchmark_summary_df$height_m, -benchmark_summary_df$p_value_composite), , drop = FALSE]
  benchmark_summary_csv <- file.path(output_dir, "composite_benchmark_B500_summary.csv")
  write.csv(benchmark_summary_df, benchmark_summary_csv, row.names = FALSE)

  combined_benchmark_summary_df <- do.call(rbind, combined_benchmark_rows)
  combined_benchmark_summary_df <- combined_benchmark_summary_df[order(-combined_benchmark_summary_df$p_value_composite), , drop = FALSE]
  combined_benchmark_summary_csv <- file.path(output_dir, "composite_benchmark_B500_combined_77_125_summary.csv")
  write.csv(combined_benchmark_summary_df, combined_benchmark_summary_csv, row.names = FALSE)

  invisible(list(
    direction_summary = direction_summary_df,
    direction_summary_csv = direction_summary_csv,
    month_summary = month_summary_df,
    month_summary_csv = month_summary_csv,
    benchmark_summary = benchmark_summary_df,
    benchmark_summary_csv = benchmark_summary_csv,
    combined_benchmark_summary = combined_benchmark_summary_df,
    combined_benchmark_summary_csv = combined_benchmark_summary_csv
  ))
}

if (sys.nframe() == 0L) {
  run_risoe_oct_nov_dec_analysis()
}
