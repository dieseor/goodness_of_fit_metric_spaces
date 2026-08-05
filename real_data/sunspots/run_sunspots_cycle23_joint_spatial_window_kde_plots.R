#!/usr/bin/env Rscript

# Spatial non-cumulative window diagnostics for cycle 23 under the joint model.
# Produces parametric conditional densities and directional KDEs side-by-side.

resolve_sunspots_joint_window_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_sunspots_joint_window_path(
  "real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"
))

sunspots_joint_require_dirstats <- function() {
  if (!requireNamespace("DirStats", quietly = TRUE)) {
    stop(
      paste(
        "The DirStats package is required for the window KDE plots.",
        "Install it with install.packages('DirStats')."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sunspots_joint_lon_lat_from_xyz <- function(x) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  data.frame(
    lon_deg = atan2(x[, 2L], x[, 1L]) * 180 / pi,
    lat_deg = asin(pmin(pmax(x[, 3L], -1), 1)) * 180 / pi,
    stringsAsFactors = FALSE
  )
}

sunspots_joint_xyz_from_lon_lat <- function(lon_deg, lat_deg) {
  lon_rad <- as.numeric(lon_deg) * pi / 180
  lat_rad <- as.numeric(lat_deg) * pi / 180
  x <- cos(lat_rad) * cos(lon_rad)
  y <- cos(lat_rad) * sin(lon_rad)
  z <- sin(lat_rad)
  cbind(x, y, z)
}

sunspots_joint_spatial_rank_windows <- function(s,
                                                lower_levels = c(0, 0.10, 0.20, 0.40, 0.60),
                                                upper_levels = c(0.10, 0.20, 0.40, 0.60, 0.80),
                                                center_levels = c(0.05, 0.15, 0.30, 0.50, 0.70)) {
  s <- as.numeric(s)
  n <- length(s)
  if (n < 1L || any(!is.finite(s))) stop("`s` must contain finite values.")
  if (!(length(lower_levels) == length(upper_levels) && length(lower_levels) == length(center_levels))) {
    stop("Window level vectors must have identical lengths.")
  }

  order_index <- order(s, seq_along(s), na.last = TRUE)
  windows <- vector("list", length(lower_levels))
  summary_rows <- vector("list", length(lower_levels))
  for (i in seq_along(lower_levels)) {
    start_rank <- floor(lower_levels[[i]] * n) + 1L
    end_rank <- floor(upper_levels[[i]] * n)
    center_rank <- min(max(as.integer(ceiling(center_levels[[i]] * n)), 1L), n)
    idx <- if (end_rank >= start_rank) order_index[start_rank:end_rank] else integer(0L)
    center_index <- order_index[[center_rank]]
    windows[[i]] <- list(
      index = idx,
      center_index = center_index,
      center_s = s[[center_index]],
      lower = lower_levels[[i]],
      upper = upper_levels[[i]],
      center_level = center_levels[[i]]
    )
    summary_rows[[i]] <- data.frame(
      window_id = i,
      lower_rank_level = lower_levels[[i]],
      upper_rank_level = upper_levels[[i]],
      center_rank_level = center_levels[[i]],
      center_rank = center_rank,
      n = length(idx),
      center_s = s[[center_index]],
      stringsAsFactors = FALSE
    )
  }

  list(summary = do.call(rbind, summary_rows), windows = windows)
}

sunspots_joint_hdr_thresholds <- function(density_values, area_weights,
                                          levels = c(0.50, 0.80, 0.95)) {
  density_values <- as.numeric(density_values)
  area_weights <- as.numeric(area_weights)
  levels <- as.numeric(levels)
  if (length(density_values) != length(area_weights)) {
    stop("`density_values` and `area_weights` must have the same length.")
  }
  if (any(!is.finite(density_values)) || any(!is.finite(area_weights)) || any(area_weights < 0)) {
    stop("Density values and area weights must be finite with nonnegative weights.")
  }
  if (sum(area_weights) <= 0) stop("Area weights must sum to a positive value.")
  if (length(levels) == 0L || any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) {
    stop("`levels` must lie strictly inside (0, 1).")
  }

  ord <- order(density_values, decreasing = TRUE)
  mass <- density_values[ord] * area_weights[ord]
  total_mass <- sum(mass)
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Weighted density mass must be finite and positive.")
  }
  cumulative <- cumsum(mass) / total_mass
  thresholds <- vapply(levels, function(level) {
    density_values[ord[[which(cumulative >= level)[[1L]]]]]
  }, numeric(1L))

  data.frame(level = levels, threshold = thresholds, stringsAsFactors = FALSE)
}

sunspots_joint_select_bandwidth_lcv_emi <- function(x_window,
                                                    bandwidth_seed = 20260805L,
                                                    h_grid = exp(seq(log(0.05), log(1.5), length.out = 100L))) {
  sunspots_joint_require_dirstats()
  x_window <- jp_normalize_unit_matrix(x_window, arg_name = "`x_window`", min_ncol = 3L)
  if (nrow(x_window) < 10L) {
    stop("Each spatial window must contain at least 10 observations for KDE bandwidth selection.")
  }

  warnings <- character(0L)
  warn_handler <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  lcv <- tryCatch(
    withCallingHandlers(
      DirStats::bw_dir_lcv(data = x_window, h_grid = h_grid, plot_it = FALSE, optim = TRUE),
      warning = warn_handler
    ),
    error = function(e) NULL
  )
  if (!is.null(lcv) && is.finite(lcv$h_opt) && lcv$h_opt > 0) {
    return(list(
      h = as.numeric(lcv$h_opt),
      method = "lcv",
      warnings = warnings
    ))
  }

  emi <- sunspots_joint_with_seed(as.integer(bandwidth_seed), tryCatch(
    withCallingHandlers(
      DirStats::bw_dir_emi(data = x_window, h_grid = h_grid, plot_it = FALSE, optim = TRUE),
      warning = warn_handler
    ),
    error = function(e) NULL
  ))
  if (!is.null(emi) && is.finite(emi$h_opt) && emi$h_opt > 0) {
    return(list(
      h = as.numeric(emi$h_opt),
      method = "emi_fallback",
      warnings = warnings
    ))
  }

  stop("Failed to obtain a finite positive KDE bandwidth with LCV and EMI fallback.")
}

sunspots_joint_parametric_conditional_density <- function(x, center_s, theta_hat) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  exp(sunspots_time_varying_log_density(
    x = x,
    u = rep(as.numeric(center_s), nrow(x)),
    theta = theta_hat
  ))
}

sunspots_joint_eval_window_densities <- function(x_grid,
                                                 x_window,
                                                 center_s,
                                                 theta_hat,
                                                 bandwidth_seed,
                                                 hdr_levels = c(0.50, 0.80, 0.95),
                                                 area_weights = NULL) {
  x_grid <- jp_normalize_unit_matrix(x_grid, arg_name = "`x_grid`", min_ncol = 3L)
  x_window <- jp_normalize_unit_matrix(x_window, arg_name = "`x_window`", min_ncol = 3L)

  parametric <- sunspots_joint_parametric_conditional_density(x_grid, center_s = center_s, theta_hat = theta_hat)

  bw <- sunspots_joint_select_bandwidth_lcv_emi(
    x_window = x_window,
    bandwidth_seed = as.integer(bandwidth_seed)
  )
  kde <- DirStats::kde_dir(x = x_grid, data = x_window, h = bw$h, L = NULL)

  if (is.null(area_weights)) area_weights <- rep(1 / length(parametric), length(parametric))
  parametric_hdr <- sunspots_joint_hdr_thresholds(parametric, area_weights = area_weights, levels = hdr_levels)
  kde_hdr <- sunspots_joint_hdr_thresholds(kde, area_weights = area_weights, levels = hdr_levels)

  list(
    parametric = as.numeric(parametric),
    kde = as.numeric(kde),
    bandwidth = bw,
    hdr_parametric = parametric_hdr,
    hdr_kde = kde_hdr
  )
}

sunspots_joint_lebedev_integral <- function(density_fun) {
  sunspots_joint_require_dirstats()
  xyz <- DirStats::lebedev$xyz
  w <- DirStats::lebedev$w
  values <- as.numeric(density_fun(xyz))
  if (length(values) != nrow(xyz) || any(!is.finite(values))) {
    return(NA_real_)
  }
  4 * pi * sum(w * values)
}

sunspots_joint_fit_from_mle_row <- function(mle_row) {
  if (nrow(mle_row) != 1L) stop("`mle_row` must have one row.")
  hemisphere_regression <- as.character(mle_row$hemisphere_regression[[1L]] %||% "asymmetric")
  theta <- list(
    a_N = as.numeric(mle_row$a_N[[1L]]),
    b_N = as.numeric(mle_row$b_N[[1L]]),
    a_S = as.numeric(mle_row$a_S[[1L]]),
    b_S = as.numeric(mle_row$b_S[[1L]]),
    c = as.numeric(mle_row$c[[1L]])
  )
  theta <- sunspots_time_varying_validate_theta(theta, hemisphere_regression = hemisphere_regression)
  eta <- sunspots_joint_time_canonicalize_eta(list(
    weight1 = as.numeric(mle_row$temporal_weight1[[1L]]),
    alpha1 = as.numeric(mle_row$temporal_alpha1[[1L]]),
    beta1 = as.numeric(mle_row$temporal_beta1[[1L]]),
    alpha2 = as.numeric(mle_row$temporal_alpha2[[1L]]),
    beta2 = as.numeric(mle_row$temporal_beta2[[1L]])
  ))
  list(theta_hat = theta, eta_hat = eta, hemisphere_regression = hemisphere_regression)
}

sunspots_joint_load_window_plot_inputs <- function(input_csv,
                                                   start_date,
                                                   end_date,
                                                   dequantization_seed,
                                                   hemisphere_regression,
                                                   control,
                                                   from_joint_output_dir = NULL) {
  if (!is.null(from_joint_output_dir)) {
    retained_path <- file.path(from_joint_output_dir, "cycle23_joint_time_space_retained_data.csv")
    mle_path <- file.path(from_joint_output_dir, "cycle23_joint_time_space_mle.csv")
    if (!file.exists(retained_path) || !file.exists(mle_path)) {
      stop("`from_joint_output_dir` must contain cycle23_joint_time_space_retained_data.csv and cycle23_joint_time_space_mle.csv")
    }
    retained <- utils::read.csv(retained_path, stringsAsFactors = FALSE)
    if (!all(c("x1", "x2", "x3", "s") %in% names(retained))) {
      stop("Retained data CSV is missing one of x1, x2, x3, s.")
    }
    retained$recorded_timestamp <- as.POSIXct(retained$recorded_timestamp, tz = "UTC")
    retained$calendar_day <- as.Date(retained$calendar_day)
    retained$dequantized_timestamp <- as.POSIXct(retained$dequantized_timestamp, tz = "UTC")
    mle <- utils::read.csv(mle_path, stringsAsFactors = FALSE)
    fit <- sunspots_joint_fit_from_mle_row(mle)
    return(list(retained = retained, fit = fit, source = "existing_joint_output"))
  }

  retained <- prepare_sunspots_cycle23_joint_time_space_data(
    input_csv = input_csv,
    start_date = start_date,
    end_date = end_date,
    dequantization_seed = as.integer(dequantization_seed)
  )
  data <- sunspots_joint_validate_data(as.matrix(retained[, c("x1", "x2", "x3")]), retained$s)
  eta_hat <- fit_sunspots_joint_time_beta_mixture2(data$s, control = control)
  theta_hat <- fit_sunspots_time_varying_asymmetric_mixture(
    x = data$x,
    u = data$s,
    hemisphere_regression = hemisphere_regression,
    control = control
  )
  list(retained = retained, fit = list(theta_hat = theta_hat, eta_hat = eta_hat), source = "refit")
}

run_sunspots_cycle23_joint_spatial_window_kde_plots <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle23_joint_spatial_window_kde"),
    start_date = "1997-06-01",
    end_date = "2006-01-01",
    dequantization_seed = 20260712L,
    hemisphere_regression = "asymmetric",
    bandwidth_seed = 20260806L,
    from_joint_output_dir = NULL,
    n_lon = 240L,
    n_lat = 121L,
    hdr_levels = c(0.50, 0.80, 0.95),
    control = list()) {
  sunspots_joint_require_dirstats()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  n_lon <- as.integer(n_lon)
  n_lat <- as.integer(n_lat)
  if (!is.finite(n_lon) || n_lon < 60L || !is.finite(n_lat) || n_lat < 31L) {
    stop("`n_lon` and `n_lat` must be reasonably large integers.")
  }

  loaded <- sunspots_joint_load_window_plot_inputs(
    input_csv = input_csv,
    start_date = start_date,
    end_date = end_date,
    dequantization_seed = dequantization_seed,
    hemisphere_regression = hemisphere_regression,
    control = control,
    from_joint_output_dir = from_joint_output_dir
  )
  retained <- loaded$retained
  theta_hat <- loaded$fit$theta_hat

  x <- jp_normalize_unit_matrix(as.matrix(retained[, c("x1", "x2", "x3")]), arg_name = "`retained`", min_ncol = 3L)
  windows <- sunspots_joint_spatial_rank_windows(retained$s)

  lon_seq <- seq(-180, 180, length.out = n_lon + 1L)[- (n_lon + 1L)]
  lat_seq <- seq(-90, 90, length.out = n_lat)
  grid_df <- expand.grid(lon_deg = lon_seq, lat_deg = lat_seq)
  x_grid <- sunspots_joint_xyz_from_lon_lat(grid_df$lon_deg, grid_df$lat_deg)
  area_weights <- pmax(cos(grid_df$lat_deg * pi / 180), 0)
  area_weights <- area_weights / sum(area_weights)

  per_window <- vector("list", nrow(windows$summary))
  bandwidth_rows <- vector("list", nrow(windows$summary))
  integral_rows <- vector("list", nrow(windows$summary))

  for (i in seq_len(nrow(windows$summary))) {
    window_meta <- windows$windows[[i]]
    idx <- window_meta$index
    if (length(idx) < 10L) {
      stop(sprintf("Window %d has only %d observations; expected at least 10.", i, length(idx)))
    }
    x_window <- x[idx, , drop = FALSE]
    eval_out <- sunspots_joint_eval_window_densities(
      x_grid = x_grid,
      x_window = x_window,
      center_s = window_meta$center_s,
      theta_hat = theta_hat,
      bandwidth_seed = as.integer(bandwidth_seed) + i - 1L,
      hdr_levels = hdr_levels,
      area_weights = area_weights
    )

    parametric_integral <- sunspots_joint_lebedev_integral(function(x_eval) {
      sunspots_joint_parametric_conditional_density(
        x = x_eval,
        center_s = window_meta$center_s,
        theta_hat = theta_hat
      )
    })
    kde_integral <- sunspots_joint_lebedev_integral(function(x_eval) {
      DirStats::kde_dir(x = x_eval, data = x_window, h = eval_out$bandwidth$h, L = NULL)
    })

    per_window[[i]] <- list(
      meta = window_meta,
      parametric = eval_out$parametric,
      kde = eval_out$kde,
      hdr_parametric = eval_out$hdr_parametric,
      hdr_kde = eval_out$hdr_kde,
      x_window = x_window
    )

    bandwidth_rows[[i]] <- data.frame(
      window_id = i,
      n = length(idx),
      center_s = window_meta$center_s,
      bandwidth = eval_out$bandwidth$h,
      bandwidth_method = eval_out$bandwidth$method,
      warnings = paste(unique(eval_out$bandwidth$warnings), collapse = " | "),
      stringsAsFactors = FALSE
    )

    integral_rows[[i]] <- data.frame(
      window_id = i,
      center_s = window_meta$center_s,
      n = length(idx),
      parametric_integral_lebedev = parametric_integral,
      kde_integral_lebedev = kde_integral,
      parametric_is_finite = is.finite(parametric_integral),
      kde_is_finite = is.finite(kde_integral),
      stringsAsFactors = FALSE
    )

    output_grid <- data.frame(
      window_id = i,
      lon_deg = grid_df$lon_deg,
      lat_deg = grid_df$lat_deg,
      parametric_density = eval_out$parametric,
      kde_density = eval_out$kde,
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      output_grid,
      file.path(output_dir, sprintf("cycle23_joint_spatial_window_%02d_density_grid.csv", i)),
      row.names = FALSE
    )
  }

  bandwidth_df <- do.call(rbind, bandwidth_rows)
  integrals_df <- do.call(rbind, integral_rows)

  zlim <- c(
    min(c(vapply(per_window, function(x) min(x$parametric), numeric(1L)),
          vapply(per_window, function(x) min(x$kde), numeric(1L)))),
    max(c(vapply(per_window, function(x) max(x$parametric), numeric(1L)),
          vapply(per_window, function(x) max(x$kde), numeric(1L))))
  )

  palette_colors <- hcl.colors(256L, palette = "YlOrRd", rev = FALSE)
  plot_path <- file.path(output_dir, "cycle23_joint_spatial_windows_parametric_vs_kde.png")
  grDevices::png(plot_path, width = 2400, height = 1200, res = 150)
  old_par <- par(mfrow = c(2, nrow(windows$summary)), mar = c(2.4, 2.6, 2.4, 1.2), oma = c(0, 0, 3, 0))
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  for (row_type in c("parametric", "kde")) {
    for (i in seq_len(nrow(windows$summary))) {
      density_values <- if (identical(row_type, "parametric")) per_window[[i]]$parametric else per_window[[i]]$kde
      zmat <- matrix(density_values, nrow = length(lon_seq), ncol = length(lat_seq), byrow = FALSE)
      image(
        x = lon_seq,
        y = lat_seq,
        z = zmat,
        col = palette_colors,
        zlim = zlim,
        xlim = c(-180, 180),
        ylim = c(-90, 90),
        xlab = if (identical(row_type, "kde")) "Longitude (degrees)" else "",
        ylab = if (i == 1L) {
          if (identical(row_type, "parametric")) "Parametric" else "KDE"
        } else "",
        main = sprintf("W%d (%.2f, %.2f] | n=%d", i,
                       windows$summary$lower_rank_level[[i]],
                       windows$summary$upper_rank_level[[i]],
                       windows$summary$n[[i]])
      )
      contour_levels <- if (identical(row_type, "parametric")) {
        per_window[[i]]$hdr_parametric$threshold
      } else {
        per_window[[i]]$hdr_kde$threshold
      }
      contour(
        x = lon_seq,
        y = lat_seq,
        z = zmat,
        levels = contour_levels,
        add = TRUE,
        drawlabels = FALSE,
        lwd = 1.1,
        col = c("#1f78b4", "#e31a1c", "#111111")
      )
      ll <- sunspots_joint_lon_lat_from_xyz(per_window[[i]]$x_window)
      points(ll$lon_deg, ll$lat_deg, pch = 16, cex = 0.22,
             col = grDevices::adjustcolor("#111111", alpha.f = 0.35))
      box()
    }
  }
  mtext("Cycle 23 spatial windows: conditional parametric density vs DirStats KDE", outer = TRUE, cex = 1.2)

  windows_out <- windows$summary
  windows_out$center_date <- if ("calendar_day" %in% names(retained)) {
    as.character(retained$calendar_day[vapply(windows$windows, `[[`, integer(1L), "center_index")])
  } else NA_character_
  windows_out$center_timestamp <- if ("dequantized_timestamp" %in% names(retained)) {
    as.character(retained$dequantized_timestamp[vapply(windows$windows, `[[`, integer(1L), "center_index")])
  } else NA_character_

  utils::write.csv(windows_out, file.path(output_dir, "cycle23_joint_spatial_window_definitions.csv"), row.names = FALSE)
  utils::write.csv(bandwidth_df, file.path(output_dir, "cycle23_joint_spatial_window_bandwidths.csv"), row.names = FALSE)
  utils::write.csv(integrals_df, file.path(output_dir, "cycle23_joint_spatial_window_integrals.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    source = loaded$source,
    input_csv = input_csv,
    start_date = start_date,
    end_date_exclusive = end_date,
    dequantization_seed = as.integer(dequantization_seed),
    bandwidth_seed = as.integer(bandwidth_seed),
    hemisphere_regression = hemisphere_regression,
    hdr_levels = paste(as.numeric(hdr_levels), collapse = ","),
    n_lon = n_lon,
    n_lat = n_lat,
    from_joint_output_dir = from_joint_output_dir %||% "",
    stringsAsFactors = FALSE
  ), file.path(output_dir, "cycle23_joint_spatial_window_metadata.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  invisible(list(
    windows = windows_out,
    bandwidths = bandwidth_df,
    integrals = integrals_df,
    plot_path = plot_path,
    output_dir = output_dir
  ))
}

parse_sunspots_joint_spatial_window_kde_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  out <- list()
  integer_keys <- c("dequantization_seed", "bandwidth_seed", "n_lon", "n_lat")
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --input_csv=PATH --output_dir=PATH --start_date=YYYY-MM-DD --end_date=YYYY-MM-DD ",
        "--dequantization_seed=INTEGER --bandwidth_seed=INTEGER --hemisphere_regression=asymmetric|shared ",
        "--from_joint_output_dir=PATH --n_lon=INTEGER --n_lat=INTEGER --hdr_levels=0.5,0.8,0.95\n"
      ))
      quit(save = "no", status = 0L)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop(sprintf("Invalid option: %s", arg))
    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key %in% c("input_csv", "output_dir", "start_date", "end_date", "hemisphere_regression", "from_joint_output_dir")) {
      out[[key]] <- value
    }
    if (key %in% integer_keys) out[[key]] <- as.integer(value)
    if (identical(key, "hdr_levels")) out[[key]] <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(run_sunspots_cycle23_joint_spatial_window_kde_plots, parse_sunspots_joint_spatial_window_kde_args())
}
