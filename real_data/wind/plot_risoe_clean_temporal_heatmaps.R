`%||%` <- function(x, y) if (is.null(x)) y else x

this_file_path <- function() {
  frame_files <- Filter(
    Negate(is.null),
    lapply(sys.frames(), function(env) env$ofile %||% NULL)
  )
  if (length(frame_files) > 0L) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path("real_data", "wind", "plot_risoe_clean_temporal_heatmaps.R"), winslash = "/", mustWork = TRUE)
}

resolve_repo_path <- function(...) {
  script_dir <- dirname(this_file_path())
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  file.path(repo_root, ...)
}

source(resolve_repo_path("real_data", "wind", "preprocess_risoe_modern_hvmf.R"))
source(resolve_repo_path("bootstrap", "model_specs.R"))
source(resolve_repo_path("bootstrap", "multiplier_bootstrap.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}
if (!requireNamespace("viridisLite", quietly = TRUE)) {
  stop("Package 'viridisLite' is required. Install it with install.packages('viridisLite').")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required. Install it with install.packages('scales').")
}

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

apply_height_cleaning <- function(df, height_m) {
  if (height_m == 77L) {
    keep <- as.Date(df$datetime, tz = "UTC") < as.Date("2007-08-01")
    return(df[keep, , drop = FALSE])
  }
  if (height_m == 125L) {
    keep <- as.Date(df$datetime, tz = "UTC") < as.Date("2004-12-01")
    return(df[keep, , drop = FALSE])
  }
  stop("Unsupported height for cleaning.")
}

month_set_definitions <- function() {
  list(
    oct = c(10L),
    nov = c(11L),
    dec = c(12L),
    oct_nov = c(10L, 11L),
    nov_dec = c(11L, 12L),
    oct_nov_dec = c(10L, 11L, 12L)
  )
}

month_set_label <- function(months_id) {
  labels <- c(
    oct = "October",
    nov = "November",
    dec = "December",
    oct_nov = "October + November",
    nov_dec = "November + December",
    oct_nov_dec = "October + November + December"
  )
  labels[[months_id]]
}

month_set_plot_tag <- function(months_id) {
  tags <- c(
    oct = "oct",
    nov = "nov",
    dec = "dec",
    oct_nov = "oct_nov",
    nov_dec = "nov_dec",
    oct_nov_dec = "oct_nov_dec"
  )
  tags[[months_id]]
}

build_temporal_scale <- function(months) {
  start_date <- as.Date(sprintf("2001-%02d-01", min(months)))
  end_day <- if (max(months) == 10L) {
    31L
  } else if (max(months) == 11L) {
    30L
  } else {
    31L
  }
  end_date <- as.Date(sprintf("2001-%02d-%02d", max(months), end_day))
  anchor_months <- c(10L, 11L, 12L)
  anchor_dates <- as.Date(sprintf("2001-%02d-01", anchor_months))
  keep_anchor <- anchor_dates >= start_date & anchor_dates <= end_date

  format_anchor <- function(x) {
    paste(as.integer(format(x, "%d")), format(x, "%b"))
  }

  list(
    start_date = start_date,
    end_date = end_date,
    limits = c(1, as.integer(end_date - start_date) + 1L),
    breaks = as.integer(anchor_dates[keep_anchor] - start_date) + 1L,
    labels = vapply(anchor_dates[keep_anchor], format_anchor, character(1))
  )
}

build_clean_hvmf_case <- function(selected_df,
                                  years,
                                  months,
                                  months_id,
                                  day_pattern,
                                  speed_col,
                                  direction_col,
                                  height_m,
                                  fixed_tz = "UTC") {
  keep <- selected_df$year %in% years &
    selected_df$month %in% months &
    selected_df$day %in% as.integer(day_pattern)
  rows <- selected_df[keep, , drop = FALSE]
  rows <- apply_height_cleaning(rows, height_m = height_m)

  valid <- is.finite(rows[[speed_col]]) &
    is.finite(rows[[direction_col]]) &
    rows[[speed_col]] > 0
  rows <- rows[valid, , drop = FALSE]
  rownames(rows) <- NULL

  stopifnot(!anyDuplicated(as.Date(rows$datetime, tz = fixed_tz)))

  direction_deg <- ((as.numeric(rows[[direction_col]]) %% 360) + 360) %% 360
  speed <- as.numeric(rows[[speed_col]])
  speed_mean <- mean(speed)
  speed_scaled <- speed / speed_mean
  angle_rad <- direction_deg * pi / 180

  pseudo_date <- as.Date(sprintf("2001-%02d-%02d", rows$month, rows$day))
  temporal_scale <- build_temporal_scale(months)
  temporal_position <- as.integer(pseudo_date - temporal_scale$start_date) + 1L

  output <- data.frame(
    datetime = rows$datetime,
    year = rows$year,
    month = rows$month,
    day = rows$day,
    hour = rows$hour,
    minute = rows$minute,
    height_m = as.integer(height_m),
    months_id = rep(months_id, nrow(rows)),
    speed = speed,
    direction_deg = direction_deg,
    speed_mean = rep(speed_mean, nrow(rows)),
    speed_scaled = speed_scaled,
    angle_rad = angle_rad,
    x0 = cosh(speed_scaled),
    x1 = sinh(speed_scaled) * cos(angle_rad),
    x2 = sinh(speed_scaled) * sin(angle_rad),
    temporal_position = temporal_position,
    temporal_label = format(pseudo_date, "%d %b"),
    stringsAsFactors = FALSE
  )
  output$minkowski_norm <- -output$x0^2 + output$x1^2 + output$x2^2

  stopifnot(all(abs(output$minkowski_norm + 1) < 1e-8))
  output
}

plot_clean_temporal_heatmap <- function(df, output_png, title_prefix) {
  temporal_scale <- build_temporal_scale(sort(unique(df$month)))
  colors <- viridisLite::magma(256, begin = 0.1, end = 0.95, direction = -1)
  density_palette <- colorRampPalette(c("#f4f4f4", "#d6d6d6", "#9f9f9f", "#4f4f4f"))

  subtitle_text <- sprintf(
    "n = %d | years = %s | speed mean = %.3f",
    nrow(df),
    paste(sort(unique(df$year)), collapse = ","),
    unique(df$speed_mean)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x1, y = x2)) +
    ggplot2::stat_bin_2d(
      bins = 30,
      ggplot2::aes(fill = ggplot2::after_stat(count)),
      alpha = 0.65,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(colors = density_palette(256), trans = "log10") +
    ggplot2::geom_point(
      ggplot2::aes(color = temporal_position),
      size = 2.4,
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

  ggplot2::ggsave(
    filename = output_png,
    plot = p,
    width = 8.5,
    height = 7,
    dpi = 180
  )
}

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
                                         B = 5000L,
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
  log_file <- file.path(logs_dir, paste0(dataset_id, "_composite_B500_log.txt"))
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
  result_rds <- file.path(results_dir, paste0(dataset_id, "_cvm_composite_B500_result.rds"))
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

run_risoe_clean_temporal_heatmaps <- function(input_nc = "real_data/wind/risoe_m_all.nc",
                                              output_dir = file.path("wind", "clean_temporal_heatmaps"),
                                              years = c(1997:2001, 2003:2007),
                                              day_pattern = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
                                              B = 5000L,
                                              n_cores = 12L,
                                              profile_method = "tabulated") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  metadata <- read_risoe_nc_metadata(input_nc)
  all_df <- load_risoe_concurrent(input_nc)
  selected_df <- select_noon_all_months(all_df, fixed_tz = metadata$fixed_tz)

  month_sets <- month_set_definitions()
  height_specs <- list(
    list(height_m = 77L, speed_col = "ws77", direction_col = "wd77"),
    list(height_m = 125L, speed_col = "ws125", direction_col = "wd125")
  )

  summary_rows <- list()
  idx <- 1L

  for (spec in height_specs) {
    for (months_id in names(month_sets)) {
      months <- month_sets[[months_id]]
      df_case <- build_clean_hvmf_case(
        selected_df = selected_df,
        years = years,
        months = months,
        months_id = months_id,
        day_pattern = day_pattern,
        speed_col = spec$speed_col,
        direction_col = spec$direction_col,
        height_m = spec$height_m,
        fixed_tz = metadata$fixed_tz
      )

      stopifnot(all(df_case$month %in% months))
      stopifnot(all(df_case$day %in% day_pattern))
      stopifnot(!anyDuplicated(as.Date(df_case$datetime, tz = metadata$fixed_tz)))
      if (spec$height_m == 77L) {
        stopifnot(all(as.Date(df_case$datetime, tz = metadata$fixed_tz) < as.Date("2007-08-01")))
      }
      if (spec$height_m == 125L) {
        stopifnot(all(as.Date(df_case$datetime, tz = metadata$fixed_tz) < as.Date("2004-12-01")))
      }

      output_plot <- file.path(
        output_dir,
        sprintf("risoe_%sm_%s_heatmap.png", spec$height_m, month_set_plot_tag(months_id))
      )
      plot_clean_temporal_heatmap(
        df = df_case,
        output_png = output_plot,
        title_prefix = sprintf("Risoe %sm | %s", spec$height_m, month_set_label(months_id))
      )

      benchmark_row <- run_composite_benchmark_case(
        df = df_case,
        dataset_id = sprintf("risoe_clean_%sm_%s_set12", spec$height_m, month_set_plot_tag(months_id)),
        output_dir = output_dir,
        B = B,
        n_cores = n_cores,
        profile_method = profile_method
      )
      benchmark_row$output_plot <- output_plot
      summary_rows[[idx]] <- benchmark_row
      idx <- idx + 1L
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(summary_df$height_m, summary_df$months_id), , drop = FALSE]
  rownames(summary_df) <- NULL
  summary_csv <- file.path(output_dir, "risoe_clean_temporal_heatmaps_summary.csv")
  write.csv(summary_df, summary_csv, row.names = FALSE)

  stopifnot(nrow(summary_df) == 12L)
  stopifnot(all(file.exists(summary_df$output_plot)))
  stopifnot(all(file.exists(summary_df$processed_csv)))
  stopifnot(all(file.exists(summary_df$result_rds)))
  stopifnot(all(file.exists(summary_df$log_file)))
  stopifnot(all(is.finite(summary_df$p_value_composite)))

  invisible(list(
    metadata = metadata,
    summary = summary_df,
    summary_csv = summary_csv
  ))
}

if (sys.nframe() == 0L) {
  run_risoe_clean_temporal_heatmaps()
}
