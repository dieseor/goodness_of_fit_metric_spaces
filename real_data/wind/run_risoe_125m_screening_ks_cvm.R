source(file.path("real_data", "wind", "preprocess_risoe_modern_hvmf.R"))
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("bootstrap", "calibration_study.R"))
source(file.path("scripts", "path_helpers.R"))

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
  selected$date <- NULL
  selected
}

apply_height_cleaning <- function(df, height_m, fixed_tz = "UTC") {
  date_vec <- as.Date(df$datetime, tz = fixed_tz)

  if (height_m == 125L) {
    return(df[date_vec < as.Date("2004-12-01"), , drop = FALSE])
  }

  stop("This runner only supports the cleaned 125m series.")
}

build_clean_hvmf_case <- function(selected_df,
                                  years,
                                  months,
                                  window_id,
                                  day_pattern,
                                  pattern_id,
                                  speed_col,
                                  direction_col,
                                  height_m,
                                  fixed_tz = "UTC") {
  keep <- selected_df$year %in% years &
    selected_df$month %in% as.integer(months) &
    selected_df$day %in% as.integer(day_pattern)
  rows <- selected_df[keep, , drop = FALSE]
  rows <- apply_height_cleaning(rows, height_m = height_m, fixed_tz = fixed_tz)

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

  out <- data.frame(
    dataset_id = sprintf("risoe_clean_125m_%s_%s_1996_2003", window_id, pattern_id),
    window = window_id,
    pattern = pattern_id,
    datetime = rows$datetime,
    year = rows$year,
    month = rows$month,
    day = rows$day,
    hour = rows$hour,
    minute = rows$minute,
    height_m = as.integer(height_m),
    speed = speed,
    direction_deg = direction_deg,
    speed_mean = speed_mean,
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

month_window_map <- function() {
  list(
    jan = 1L,
    feb = 2L,
    jan_feb = c(1L, 2L),
    feb_mar = c(2L, 3L),
    oct = 10L,
    nov = 11L,
    oct_nov = c(10L, 11L)
  )
}

day_pattern_map <- function() {
  list(
    start1 = c(1L, 5L, 9L, 13L, 17L, 21L, 25L, 29L),
    start2 = c(2L, 6L, 10L, 14L, 18L, 22L, 26L, 30L),
    start3 = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
    start4 = c(4L, 8L, 12L, 16L, 20L, 24L, 28L)
  )
}

make_risoe_125m_screening_configs <- function() {
  windows <- month_window_map()
  patterns <- day_pattern_map()
  configs <- list()
  idx <- 1L

  for (window_id in names(windows)) {
    for (pattern_id in names(patterns)) {
      configs[[idx]] <- list(
        dataset_id = sprintf("risoe_clean_125m_%s_%s_1996_2003", window_id, pattern_id),
        window = window_id,
        months = windows[[window_id]],
        pattern = pattern_id,
        day_pattern = patterns[[pattern_id]]
      )
      idx <- idx + 1L
    }
  }

  configs
}

safe_error_message <- function(e) {
  msg <- conditionMessage(e)
  if (!is.character(msg) || length(msg) != 1L || !nzchar(msg)) {
    return("Unknown error.")
  }
  msg
}

format_pvalue_tex <- function(p_value, digits = 4L) {
  cell <- sprintf(paste0("%.", digits, "f"), p_value)
  if (is.finite(p_value) && p_value >= 0.05) {
    return(sprintf("\\textbf{%s}", cell))
  }
  cell
}

make_ks_grid_for_case <- function(X,
                                  mu,
                                  ks_n_omega = 50L,
                                  ks_n_t = 100L,
                                  chi_margin = 0.25,
                                  t_margin = 0.25) {
  make_hvmf_ks_grid(
    data = X,
    mu = mu,
    n_omega = ks_n_omega,
    n_t = ks_n_t,
    chi_margin = chi_margin,
    t_margin = t_margin
  )
}

run_single_risoe_125m_case <- function(config,
                                       selected_df,
                                       B = 1000L,
                                       n_cores = 12L,
                                       bootstrap_method = "reestimated",
                                       fixed_tz = "UTC",
                                       ks_n_omega = 50L,
                                       ks_n_t = 100L,
                                       chi_margin = 0.25,
                                       t_margin = 0.25,
                                       seed = 91000L) {
  start_time <- Sys.time()

  tryCatch(
    {
      df_case <- build_clean_hvmf_case(
        selected_df = selected_df,
        years = 1996:2003,
        months = config$months,
        window_id = config$window,
        day_pattern = config$day_pattern,
        pattern_id = config$pattern,
        speed_col = "ws125",
        direction_col = "wd125",
        height_m = 125L,
        fixed_tz = fixed_tz
      )

      X <- as.matrix(df_case[, c("x0", "x1", "x2")])
      fit <- hvmf_mle_h2(X)
      ks_grid <- make_ks_grid_for_case(
        X = X,
        mu = fit$mu,
        ks_n_omega = ks_n_omega,
        ks_n_t = ks_n_t,
        chi_margin = chi_margin,
        t_margin = t_margin
      )

      result <- multiplier_bootstrap_hvmf(
        data = X,
        null = list(type = "composite"),
        statistics = c("ks", "cvm"),
        ks_grid = ks_grid,
        B = B,
        alpha = 0.05,
        n_cores = n_cores,
        seed = seed,
        bootstrap_method = bootstrap_method,
        keep = list(
          observed_process = FALSE,
          bootstrap_statistics = TRUE,
          bootstrap_thetas = FALSE
        ),
        control = list(hvmf_profile_method = "tabulated"),
        unknown_param = "both"
      )

      elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      common <- data.frame(
        window = config$window,
        pattern = config$pattern,
        B = B,
        n = nrow(df_case),
        kappa_hat = fit$kappa,
        theta_deg_hat = fit$theta_deg,
        ks_n_omega = ks_n_omega,
        ks_n_t = ks_n_t,
        elapsed_seconds = elapsed_seconds,
        status = "ok",
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )

      rbind(
        cbind(
          common,
          data.frame(
            statistic = "KS",
            p_value = result$inference$ks$p_value,
            observed_statistic = result$observed$ks$statistic,
            stringsAsFactors = FALSE
          )
        ),
        cbind(
          common,
          data.frame(
            statistic = "CvM",
            p_value = result$inference$cvm$p_value,
            observed_statistic = result$observed$cvm$statistic,
            stringsAsFactors = FALSE
          )
        )
      )
    },
    error = function(e) {
      elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      common <- data.frame(
        window = config$window,
        pattern = config$pattern,
        B = B,
        n = NA_integer_,
        kappa_hat = NA_real_,
        theta_deg_hat = NA_real_,
        ks_n_omega = ks_n_omega,
        ks_n_t = ks_n_t,
        elapsed_seconds = elapsed_seconds,
        status = "error",
        error_message = safe_error_message(e),
        stringsAsFactors = FALSE
      )

      rbind(
        cbind(
          common,
          data.frame(
            statistic = "KS",
            p_value = NA_real_,
            observed_statistic = NA_real_,
            stringsAsFactors = FALSE
          )
        ),
        cbind(
          common,
          data.frame(
            statistic = "CvM",
            p_value = NA_real_,
            observed_statistic = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      )
    }
  )
}

write_risoe_125m_latex_table <- function(results_df, output_txt) {
  windows_order <- c("jan", "feb", "jan_feb", "feb_mar", "oct", "nov", "oct_nov")
  patterns_order <- c("start1", "start2", "start3", "start4")

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Composite HvMF goodness-of-fit screening for the cleaned Ris{\\o} $125\\,\\mathrm{m}$ series (1996--2003). The table reports multiplier-bootstrap $p$-values with $B=1000$ for the KS and CvM statistics under four four-day subsampling patterns. Bold entries are not significant at the $5\\%$ level.}",
    "\\label{tab:risoe-wind-screening-125m}",
    "\\begin{tabular}{lcccccccc}",
    "\\toprule",
    " & \\multicolumn{2}{c}{\\texttt{start1}} & \\multicolumn{2}{c}{\\texttt{start2}} & \\multicolumn{2}{c}{\\texttt{start3}} & \\multicolumn{2}{c}{\\texttt{start4}} \\\\",
    "Month(s) & KS & CvM & KS & CvM & KS & CvM & KS & CvM \\\\",
    "\\midrule"
  )

  for (window_id in windows_order) {
    row_df <- results_df[results_df$window == window_id & results_df$status == "ok", , drop = FALSE]
    row_cells <- character()

    for (pattern_id in patterns_order) {
      ks_row <- row_df[row_df$pattern == pattern_id & row_df$statistic == "KS", , drop = FALSE]
      cvm_row <- row_df[row_df$pattern == pattern_id & row_df$statistic == "CvM", , drop = FALSE]

      ks_cell <- if (nrow(ks_row) == 1L) format_pvalue_tex(ks_row$p_value[[1L]]) else "--"
      cvm_cell <- if (nrow(cvm_row) == 1L) format_pvalue_tex(cvm_row$p_value[[1L]]) else "--"

      row_cells <- c(row_cells, ks_cell, cvm_cell)
    }

    lines <- c(
      lines,
      sprintf(
        "\\texttt{%s} & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        window_id,
        row_cells[[1L]], row_cells[[2L]],
        row_cells[[3L]], row_cells[[4L]],
        row_cells[[5L]], row_cells[[6L]],
        row_cells[[7L]], row_cells[[8L]]
      )
    )
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, con = output_txt)
  invisible(output_txt)
}

run_risoe_125m_screening_ks_cvm <- function(input_nc = file.path("wind", "risoe_m_all.nc"),
                                            output_dir = canonical_wind_screening_dir("slow"),
                                            B = 1000L,
                                            n_cores = 12L,
                                            bootstrap_method = "reestimated",
                                            fixed_tz = "UTC",
                                            ks_n_omega = 50L,
                                            ks_n_t = 100L,
                                            chi_margin = 0.25,
                                            t_margin = 0.25,
                                            seed = 91000L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  all_df <- load_risoe_concurrent(input_nc, fixed_tz = fixed_tz)
  selected_df <- select_noon_all_months(all_df, fixed_tz = fixed_tz)
  configs <- make_risoe_125m_screening_configs()

  rows <- vector("list", length(configs))

  for (i in seq_along(configs)) {
    config <- configs[[i]]
    cat(sprintf(
      "\n[%d/%d] %s | window=%s | pattern=%s\n",
      i, length(configs), config$dataset_id, config$window, config$pattern
    ))
    rows[[i]] <- run_single_risoe_125m_case(
      config = config,
      selected_df = selected_df,
      B = B,
      n_cores = n_cores,
      bootstrap_method = bootstrap_method,
      fixed_tz = fixed_tz,
      ks_n_omega = ks_n_omega,
      ks_n_t = ks_n_t,
      chi_margin = chi_margin,
      t_margin = t_margin,
      seed = seed + i
    )
  }

  results_df <- do.call(rbind, rows)
  rownames(results_df) <- NULL

  csv_path <- file.path(output_dir, "risoe_125m_ks_cvm_screening_b1000.csv")
  txt_path <- file.path(output_dir, "risoe_125m_ks_cvm_screening_table_for_tex.txt")
  write.csv(results_df, csv_path, row.names = FALSE)
  write_risoe_125m_latex_table(results_df, txt_path)

  cat("\nOutput CSV:\n", csv_path, "\n", sep = "")
  cat("Output TeX table:\n", txt_path, "\n", sep = "")

  invisible(list(results = results_df, csv = csv_path, txt = txt_path))
}

if (sys.nframe() == 0L) {
  run_risoe_125m_screening_ks_cvm()
}
