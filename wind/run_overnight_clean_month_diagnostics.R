`%||%` <- function(x, y) if (is.null(x)) y else x

this_file_path <- function() {
  frame_files <- Filter(
    Negate(is.null),
    lapply(sys.frames(), function(env) env$ofile %||% NULL)
  )
  if (length(frame_files) > 0L) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path("wind", "run_overnight_clean_month_diagnostics.R"), winslash = "/", mustWork = TRUE)
}

resolve_repo_path <- function(...) {
  script_dir <- dirname(this_file_path())
  repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
  file.path(repo_root, ...)
}

source(resolve_repo_path("wind", "analyze_risoe_oct_nov_dec_patterns.R"))
source(resolve_repo_path("wind", "plot_risoe_clean_temporal_heatmaps.R"))

build_temporal_scale <- function(months) {
  months <- sort(unique(as.integer(months)))
  start_date <- as.Date(sprintf("2001-%02d-01", min(months)))
  end_month <- max(months)
  end_date <- seq(start_date, by = "1 month", length.out = length(months) + 1L)[[length(months) + 1L]] - 1L

  anchor_dates <- as.Date(sprintf("2001-%02d-01", months))
  labels <- paste(as.integer(format(anchor_dates, "%d")), format(anchor_dates, "%b"))

  list(
    start_date = start_date,
    end_date = end_date,
    limits = c(1, as.integer(end_date - start_date) + 1L),
    breaks = as.integer(anchor_dates - start_date) + 1L,
    labels = labels,
    end_month = end_month
  )
}

month_name_short <- function(month) {
  c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")[month]
}

month_name_label <- function(month) {
  c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[month]
}

cardinal <- function(deg) {
  dirs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  idx <- floor(((deg %% 360) + 22.5) / 45) %% 8 + 1
  dirs[idx]
}

build_all_days_hvmf_case <- function(df,
                                     speed_col,
                                     direction_col,
                                     height_m,
                                     months,
                                     fixed_tz = "UTC") {
  df <- apply_height_cleaning(df[df$month %in% months, , drop = FALSE], height_m = height_m)
  valid <- is.finite(df[[speed_col]]) &
    is.finite(df[[direction_col]]) &
    df[[speed_col]] > 0
  rows <- df[valid, , drop = FALSE]
  rownames(rows) <- NULL

  direction_deg <- ((as.numeric(rows[[direction_col]]) %% 360) + 360) %% 360
  speed <- as.numeric(rows[[speed_col]])
  speed_mean <- mean(speed)
  speed_scaled <- speed / speed_mean
  angle_rad <- direction_deg * pi / 180

  out <- data.frame(
    datetime = rows$datetime,
    year = rows$year,
    month = rows$month,
    day = rows$day,
    hour = rows$hour,
    minute = rows$minute,
    height_m = as.integer(height_m),
    speed = speed,
    direction_deg = direction_deg,
    speed_mean = rep(speed_mean, length(speed)),
    speed_scaled = speed_scaled,
    angle_rad = angle_rad,
    x0 = cosh(speed_scaled),
    x1 = sinh(speed_scaled) * cos(angle_rad),
    x2 = sinh(speed_scaled) * sin(angle_rad),
    stringsAsFactors = FALSE
  )
  out$minkowski_norm <- -out$x0^2 + out$x1^2 + out$x2^2
  stopifnot(all(abs(out$minkowski_norm + 1) < 1e-8))
  out
}

compute_case_diagnostics <- function(df_case) {
  theta <- df_case$angle_rad
  C <- mean(cos(theta))
  S <- mean(sin(theta))
  R <- sqrt(C^2 + S^2)

  WC <- mean(df_case$speed_scaled * cos(theta))
  WS <- mean(df_case$speed_scaled * sin(theta))
  WR <- sqrt(WC^2 + WS^2)
  weighted_theta_deg <- (atan2(WS, WC) * 180 / pi) %% 360

  h2_frechet <- fit_h2_frechet_mean_safe(as.matrix(df_case[, c("x0", "x1", "x2")]))
  s1_intrinsic <- fit_s1_intrinsic_mean(theta)
  hvmf_fit <- hvmf_mle_h2(as.matrix(df_case[, c("x0", "x1", "x2")]))

  data.frame(
    n = nrow(df_case),
    speed_mean = unique(df_case$speed_mean)[1],
    circular_resultant = R,
    circular_sd = if (R > 0) sqrt(-2 * log(R)) else Inf,
    weighted_mean_deg = weighted_theta_deg,
    weighted_cardinal = cardinal(weighted_theta_deg),
    weighted_resultant = WR,
    s1_intrinsic_theta_deg = s1_intrinsic$theta_deg,
    h2_frechet_theta_deg = h2_frechet$theta_deg,
    h2_frechet_sinh_chi = h2_frechet$sinh_chi,
    hvmf_theta_deg_hat = hvmf_fit$theta_deg,
    hvmf_kappa_hat = hvmf_fit$kappa,
    stringsAsFactors = FALSE
  )
}

compute_case_yearly_diagnostics <- function(df_case) {
  groups <- split(df_case, df_case$year)
  rows <- lapply(groups, function(g) {
    diag <- compute_case_diagnostics(g)
    data.frame(year = unique(g$year), diag, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

fit_h2_frechet_mean_safe <- function(X) {
  X <- normalize_hvmf_h2_data(X)
  start_fit <- hvmf_mle_h2(X)
  start_eta <- log(pmax(start_fit$chi, 1e-8))
  start_theta <- start_fit$theta

  objective <- function(par) {
    eta <- par[[1L]]
    theta <- par[[2L]]

    if (!is.finite(eta) || !is.finite(theta)) {
      return(Inf)
    }

    chi <- exp(eta)
    if (!is.finite(chi) || chi > 25) {
      return(Inf)
    }

    mu <- h2_point_from_chi_theta(chi, theta)
    if (any(!is.finite(mu))) {
      return(Inf)
    }

    distances <- vapply(
      seq_len(nrow(X)),
      function(i) hyperbolic_geodesic_distance_h2(mu, X[i, ]),
      numeric(1)
    )

    if (any(!is.finite(distances))) {
      return(Inf)
    }

    mean(distances^2)
  }

  fit <- stats::optim(
    par = c(start_eta, start_theta),
    fn = objective,
    method = "L-BFGS-B",
    lower = c(log(1e-8), -pi),
    upper = c(log(25), pi),
    control = list(maxit = 1000L, factr = 1e7)
  )

  if (!is.finite(fit$value)) {
    stop("Failed to compute a finite H^2 Fréchet mean.")
  }

  chi_hat <- exp(fit$par[[1L]])
  theta_hat <- wrap_angle_to_pi(fit$par[[2L]])
  mu_hat <- h2_point_from_chi_theta(chi_hat, theta_hat)

  if (any(!is.finite(mu_hat))) {
    stop("Computed non-finite H^2 Fréchet mean.")
  }

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

plot_monthly_fréchet_circle <- function(summary_df, height_m, output_png) {
  colors <- grDevices::hcl.colors(12L, "Dark 3")
  names(colors) <- as.character(1:12)
  theta <- seq(0, 2 * pi, length.out = 361L)
  max_len <- max(summary_df$h2_frechet_sinh_chi, na.rm = TRUE)

  grDevices::png(output_png, width = 1200, height = 1200, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mar = c(2, 2, 4, 2))
  graphics::plot(
    cos(theta),
    sin(theta),
    type = "l",
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = c(-1.2, 1.2),
    ylim = c(-1.2, 1.2),
    col = "#c7c7c7",
    lwd = 1.2,
    main = sprintf("Risoe %sm | monthly H^2 Fréchet mean directions", height_m)
  )
  graphics::abline(h = 0, v = 0, col = "#eeeeee", lty = 3)

  for (i in seq_len(nrow(summary_df))) {
    month_chr <- as.character(summary_df$month[[i]])
    angle <- summary_df$h2_frechet_theta_deg[[i]] * pi / 180
    radius <- 0.18 + 0.76 * summary_df$h2_frechet_sinh_chi[[i]] / max_len
    x1 <- radius * cos(angle)
    y1 <- radius * sin(angle)
    graphics::arrows(
      x0 = 0, y0 = 0, x1 = x1, y1 = y1,
      length = 0.08,
      lwd = 2.6,
      col = colors[[month_chr]]
    )
    graphics::points(x1, y1, pch = 16, cex = 1, col = colors[[month_chr]])
    graphics::text(
      1.08 * x1,
      1.08 * y1,
      labels = month_name_label(summary_df$month[[i]]),
      col = colors[[month_chr]],
      cex = 0.85,
      font = 2
    )
  }

  graphics::mtext(
    "Arrow direction from H^2 Fréchet mean; arrow length proportional to sinh(chi)",
    side = 3,
    line = 0.6,
    cex = 0.85
  )
}

make_projection_plot <- function(df_case, title_prefix, output_png) {
  temporal_scale <- build_temporal_scale(sort(unique(df_case$month)))
  colors <- viridisLite::magma(256, begin = 0.1, end = 0.95, direction = -1)
  density_palette <- colorRampPalette(c("#f4f4f4", "#d6d6d6", "#9f9f9f", "#4f4f4f"))

  pseudo_date <- as.Date(sprintf("2001-%02d-%02d", df_case$month, df_case$day))
  temporal_position <- as.integer(pseudo_date - temporal_scale$start_date) + 1L

  subtitle_text <- sprintf(
    "n = %d | years = %s | speed mean = %.3f",
    nrow(df_case),
    paste(sort(unique(df_case$year)), collapse = ","),
    unique(df_case$speed_mean)
  )

  p <- ggplot2::ggplot(df_case, ggplot2::aes(x = x1, y = x2)) +
    ggplot2::stat_bin_2d(
      bins = 30,
      ggplot2::aes(fill = ggplot2::after_stat(count)),
      alpha = 0.65,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(colors = density_palette(256), trans = "log10") +
    ggplot2::geom_point(
      ggplot2::aes(color = temporal_position),
      size = 2.2,
      alpha = 0.85
    ) +
    ggplot2::scale_color_gradientn(
      colors = colors,
      limits = temporal_scale$limits,
      breaks = temporal_scale$breaks,
      labels = temporal_scale$labels,
      oob = scales::squish,
      name = "Calendar time"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title_prefix,
      subtitle = subtitle_text,
      x = expression(x[1]),
      y = expression(x[2])
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(output_png, plot = p, width = 8.5, height = 7, dpi = 180)
}

make_contiguous_two_month_windows <- function() {
  lapply(seq_len(11L), function(start_month) {
    months <- c(start_month, start_month + 1L)
    list(
      months = months,
      months_id = paste(month_name_short(months), collapse = "_"),
      window_len = 2L
    )
  })
}

make_day_patterns_step4 <- function() {
  lapply(seq_len(4L), function(start_day) {
    days <- seq.int(start_day, 30L, by = 4L)
    list(
      start_day = start_day,
      day_pattern = days,
      day_pattern_id = paste0("start", start_day)
    )
  })
}

run_set12_window_screening <- function(selected_df,
                                       years,
                                       output_dir,
                                       B = 500L,
                                       n_cores = 10L,
                                       profile_method = "tabulated",
                                       seed = 20260527L) {
  month_windows <- make_contiguous_two_month_windows()
  day_patterns <- make_day_patterns_step4()
  height_specs <- list(
    list(height_m = 77L, speed_col = "ws77", direction_col = "wd77"),
    list(height_m = 125L, speed_col = "ws125", direction_col = "wd125")
  )
  summary_rows <- list()
  idx <- 1L

  set.seed(seed)

  for (spec in height_specs) {
    for (window in month_windows) {
      for (day_spec in day_patterns) {
        case_seed <- as.integer(seed +
          100000L * spec$height_m +
          1000L * min(window$months) +
          day_spec$start_day)
        set.seed(case_seed)

        df_case <- build_clean_hvmf_case(
          selected_df = selected_df,
          years = years,
          months = window$months,
          months_id = window$months_id,
          day_pattern = day_spec$day_pattern,
          speed_col = spec$speed_col,
          direction_col = spec$direction_col,
          height_m = spec$height_m,
          fixed_tz = "UTC"
        )

        case_id <- sprintf(
          "risoe_clean_%sm_%s_%s_%s",
          spec$height_m,
          window$months_id,
          day_spec$day_pattern_id,
          paste(range(years), collapse = "_")
        )

        benchmark_row <- suppressMessages(
          run_composite_benchmark_case(
            df = df_case,
            dataset_id = case_id,
            output_dir = output_dir,
            B = B,
            n_cores = n_cores,
            profile_method = profile_method,
            seed = case_seed
          )
        )

        benchmark_row$window_len <- window$window_len
        benchmark_row$month_start <- min(window$months)
        benchmark_row$month_end <- max(window$months)
        benchmark_row$month_window <- window$months_id
        benchmark_row$day_start <- day_spec$start_day
        benchmark_row$day_pattern <- paste(day_spec$day_pattern, collapse = ",")
        benchmark_row$includes_apr <- any(window$months == 4L)
        benchmark_row$includes_nov <- any(window$months == 11L)
        benchmark_row$years_used <- paste(years, collapse = ",")
        benchmark_row$B <- B
        benchmark_row$n_cores <- n_cores
        benchmark_row$seed <- seed
        benchmark_row$case_seed <- case_seed
        summary_rows[[idx]] <- benchmark_row
        idx <- idx + 1L

        cat(
          sprintf(
            "%sm %s %s n=%d B=%d seed=%d p=%.4f\n",
            spec$height_m,
            window$months_id,
            day_spec$day_pattern_id,
            nrow(df_case),
            B,
            case_seed,
            benchmark_row$p_value_composite
          )
        )
        flush.console()
      }
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(
    summary_df$height_m,
    summary_df$month_start,
    summary_df$day_start
  ), , drop = FALSE]
  rownames(summary_df) <- NULL

  combination_summary <- aggregate(
    p_value_composite ~ height_m + month_start + month_end + month_window,
    data = summary_df,
    FUN = function(x) {
      min(x, na.rm = TRUE)
    }
  )
  names(combination_summary)[names(combination_summary) == "p_value_composite"] <- "min_p_value_composite"
  combination_summary$all_day_starts_above_0p05 <- combination_summary$min_p_value_composite > 0.05
  combination_summary <- combination_summary[order(
    combination_summary$height_m,
    combination_summary$month_start
  ), , drop = FALSE]
  rownames(combination_summary) <- NULL

  summary_csv <- file.path(
    output_dir,
    sprintf("contiguous_two_month_step4_screening_b%d_seed%d.csv", B, seed)
  )
  combination_csv <- file.path(
    output_dir,
    sprintf("contiguous_two_month_step4_screening_summary_b%d_seed%d.csv", B, seed)
  )
  write.csv(summary_df, summary_csv, row.names = FALSE)
  write.csv(combination_summary, combination_csv, row.names = FALSE)

  list(
    results = summary_df,
    combination_summary = combination_summary,
    result_csv = summary_csv,
    combination_csv = combination_csv
  )
}

run_overnight_clean_month_diagnostics <- function(input_nc = file.path("wind", "risoe_m_all.nc"),
                                                  output_dir = file.path("wind", "overnight_clean_month_diagnostics"),
                                                  years = 1996:2003,
                                                  B = 500L,
                                                  n_cores = 10L,
                                                  seed = 20260527L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata <- read_risoe_nc_metadata(input_nc, fixed_tz = "UTC")
  raw_df <- load_risoe_concurrent(input_nc, fixed_tz = metadata$fixed_tz)
  selected_df <- select_noon_all_months(raw_df, fixed_tz = metadata$fixed_tz)

  height_specs <- list(
    list(height_m = 77L, speed_col = "ws77", direction_col = "wd77"),
    list(height_m = 125L, speed_col = "ws125", direction_col = "wd125")
  )

  month_summary_rows <- list()
  idx <- 1L

  for (spec in height_specs) {
    month_rows <- list()
    for (month in 1:12) {
      df_case <- build_all_days_hvmf_case(
        df = raw_df[raw_df$year %in% years, , drop = FALSE],
        speed_col = spec$speed_col,
        direction_col = spec$direction_col,
        height_m = spec$height_m,
        months = month,
        fixed_tz = metadata$fixed_tz
      )
      diag <- compute_case_diagnostics(df_case)
      month_rows[[month]] <- data.frame(
        height_m = spec$height_m,
        month = month,
        month_label = month_name_label(month),
        dataset_type = "all_days",
        diag,
        stringsAsFactors = FALSE
      )
    }
    month_df <- do.call(rbind, month_rows)
    rownames(month_df) <- NULL
    month_summary_rows[[idx]] <- month_df
    idx <- idx + 1L

    plot_monthly_fréchet_circle(
      summary_df = month_df,
      height_m = spec$height_m,
      output_png = file.path(output_dir, sprintf("risoe_%sm_all_days_monthly_frechet_circle.png", spec$height_m))
    )
  }

  month_summary_df <- do.call(rbind, month_summary_rows)
  write.csv(month_summary_df, file.path(output_dir, "all_days_monthly_direction_summary.csv"), row.names = FALSE)

  special_rows <- list()
  special_yearly_rows <- list()
  idx <- 1L
  idx_year <- 1L
  day_pattern <- c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L)

  for (spec in height_specs) {
    for (month in c(4L, 12L)) {
      df_case <- build_clean_hvmf_case(
        selected_df = selected_df,
        years = years,
        months = month,
        months_id = month_name_short(month),
        day_pattern = day_pattern,
        speed_col = spec$speed_col,
        direction_col = spec$direction_col,
        height_m = spec$height_m,
        fixed_tz = metadata$fixed_tz
      )
      diag <- compute_case_diagnostics(df_case)
      yearly_diag <- compute_case_yearly_diagnostics(df_case)
      yearly_diag$height_m <- spec$height_m
      yearly_diag$month <- month
      yearly_diag$month_label <- month_name_label(month)

      special_rows[[idx]] <- data.frame(
        height_m = spec$height_m,
        month = month,
        month_label = month_name_label(month),
        dataset_type = "set12_noon",
        diag,
        stringsAsFactors = FALSE
      )
      special_yearly_rows[[idx_year]] <- yearly_diag
      idx <- idx + 1L
      idx_year <- idx_year + 1L

      make_projection_plot(
        df_case = df_case,
        title_prefix = sprintf("Risoe %sm | %s | set12 noon sample", spec$height_m, month_name_label(month)),
        output_png = file.path(output_dir, sprintf("risoe_%sm_%s_set12_projection.png", spec$height_m, month_name_short(month)))
      )
    }
  }

  special_summary_df <- do.call(rbind, special_rows)
  special_yearly_df <- do.call(rbind, special_yearly_rows)
  rownames(special_summary_df) <- NULL
  rownames(special_yearly_df) <- NULL
  write.csv(special_summary_df, file.path(output_dir, "april_december_set12_diagnostics.csv"), row.names = FALSE)
  write.csv(special_yearly_df, file.path(output_dir, "april_december_set12_yearly_diagnostics.csv"), row.names = FALSE)

  windows_dir <- file.path(output_dir, sprintf("windows_contiguous_two_month_step4_b%d_seed%d", B, seed))
  dir.create(windows_dir, recursive = TRUE, showWarnings = FALSE)
  window_summary_df <- run_set12_window_screening(
    selected_df = selected_df,
    years = years,
    output_dir = windows_dir,
    B = B,
    n_cores = n_cores,
    profile_method = "tabulated",
    seed = seed
  )

  invisible(list(
    monthly_summary = month_summary_df,
    special_summary = special_summary_df,
    special_yearly_summary = special_yearly_df,
    window_results = window_summary_df$results,
    window_combination_summary = window_summary_df$combination_summary
  ))
}

if (sys.nframe() == 0L) {
  run_overnight_clean_month_diagnostics()
}
