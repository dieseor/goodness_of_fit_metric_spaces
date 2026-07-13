#!/usr/bin/env Rscript

resolve_sunspots_axial_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

utils_path_sunspots_axial <- resolve_sunspots_axial_path("utils.R")
model_specs_path_sunspots_axial <- resolve_sunspots_axial_path("bootstrap", "model_specs.R")
bootstrap_path_sunspots_axial <- resolve_sunspots_axial_path("bootstrap", "multiplier_bootstrap.R")

source(utils_path_sunspots_axial)
source(model_specs_path_sunspots_axial)
source(bootstrap_path_sunspots_axial)

suppressPackageStartupMessages({
  library(stats)
  library(utils)
  library(grDevices)
  library(graphics)
})

plot_col_sc <- "#1f77b4"
plot_col_axial <- "#d95f02"
plot_col_ref <- "gray40"

parse_axial_bool <- function(value, key) {
  value_lower <- tolower(trimws(value))
  if (value_lower %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value_lower %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("Option '--%s' must be TRUE or FALSE.", key))
}

parse_axial_integer <- function(value, key, allow_inf = FALSE) {
  if (allow_inf && tolower(trimws(value)) %in% c("inf", "infinity")) return(Inf)
  if (!grepl("^-?[0-9]+$", value)) stop(sprintf("Option '--%s' must be an integer.", key))
  as.integer(value)
}

parse_axial_cycles <- function(value) {
  parts <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) stop("Option '--cycles' cannot be empty.")
  if (any(!grepl("^-?[0-9]+$", parts))) stop("Option '--cycles' must contain comma-separated integer cycle identifiers.")
  unique(as.integer(parts))
}

parse_axial_statistics <- function(value) {
  parts <- unique(tolower(trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))))
  parts <- parts[nzchar(parts)]
  if (!identical(parts, "ks")) stop("The axial rolling-window runner currently supports only KS.")
  parts
}

parse_axial_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    cycles = c(21L, 22L, 23L),
    B = 5000L,
    n_cores = 12L,
    seed = 123L,
    resume = TRUE,
    output_dir = file.path(
      "real_data", "sunspots", "output",
      "cycles21_23_axial_marginal_rolling_windows_10cr_B1000"
    ),
    sc_csv = file.path(
      "real_data", "sunspots", "output",
      "cycles21_23_weighted_mixture_rolling_windows_10cr_B1000_all_KS",
      "sunspots_cycles21_23_weighted_mixture_rolling_10cr_gof_results.csv"
    ),
    sunspots_csv = paste(
      file.path("real_data", "sunspots", "output", "sunspots_cycle21_s2_all.csv"),
      file.path("real_data", "sunspots", "output", "sunspots_cycle22_s2_all.csv"),
      file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
      sep = ","
    ),
    ks_omega_points = 60L,
    ks_t_points = 200L,
    start_window = 1L,
    end_window = Inf,
    statistics = "ks"
  )

  help_text <- paste(
    "Usage:",
    "  Rscript real_data/sunspots/run_sunspots_axial_marginal_rolling_windows_gof.R [options]",
    "",
    "Options:",
    "  --help",
    "  --cycles=21,22,23",
    "  --B=50000",
    "  --n_cores=12",
    "  --seed=123",
    "  --resume=TRUE",
    "  --output_dir=PATH",
    "  --sc_csv=PATH",
    "  --sunspots_csv=PATH1,PATH2,PATH3",
    "  --ks_omega_points=60",
    "  --ks_t_points=200",
    "  --start_window=1",
    "  --end_window=Inf",
    "  --statistics=KS",
    sep = "\n"
  )

  if (any(args %in% c("--help", "-h"))) {
    cat(help_text, sep = "\n")
    quit(save = "no", status = 0L)
  }

  if (length(args) == 0L) return(defaults)

  parsed <- list()
  for (arg in args) {
    key_value <- sub("^--", "", arg)
    parts <- strsplit(key_value, "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else NULL
    if (is.null(value)) stop(sprintf("Option '--%s' requires the form --%s=value.", key, key))
    parsed[[key]] <- value
  }

  if ("cycles" %in% names(parsed)) defaults$cycles <- parse_axial_cycles(parsed$cycles)
  if ("B" %in% names(parsed)) defaults$B <- parse_axial_integer(parsed$B, "B")
  if ("n_cores" %in% names(parsed)) defaults$n_cores <- parse_axial_integer(parsed$n_cores, "n_cores")
  if ("seed" %in% names(parsed)) defaults$seed <- parse_axial_integer(parsed$seed, "seed")
  if ("resume" %in% names(parsed)) defaults$resume <- parse_axial_bool(parsed$resume, "resume")
  if ("output_dir" %in% names(parsed)) defaults$output_dir <- parsed$output_dir
  if ("sc_csv" %in% names(parsed)) defaults$sc_csv <- parsed$sc_csv
  if ("sunspots_csv" %in% names(parsed)) defaults$sunspots_csv <- parsed$sunspots_csv
  if ("ks_omega_points" %in% names(parsed)) defaults$ks_omega_points <- parse_axial_integer(parsed$ks_omega_points, "ks_omega_points")
  if ("ks_t_points" %in% names(parsed)) defaults$ks_t_points <- parse_axial_integer(parsed$ks_t_points, "ks_t_points")
  if ("start_window" %in% names(parsed)) defaults$start_window <- parse_axial_integer(parsed$start_window, "start_window")
  if ("end_window" %in% names(parsed)) defaults$end_window <- parse_axial_integer(parsed$end_window, "end_window", allow_inf = TRUE)
  if ("statistics" %in% names(parsed)) defaults$statistics <- parse_axial_statistics(parsed$statistics)
  if ("statistic" %in% names(parsed)) defaults$statistics <- parse_axial_statistics(parsed$statistic)
  defaults
}

as_utc <- function(x) {
  as.POSIXct(x, tz = "UTC")
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[ok], y[ok])
}

read_axial_sc_results <- function(path, cycles) {
  if (!file.exists(path)) {
    stop("Small-circle rolling results CSV not found: ", path)
  }

  sc <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "cycle", "window_id", "start_date", "end_date", "n", "statistic_type",
    "p_value_raw", "pi_hat", "mu_1", "mu_2", "mu_3", "kappa1_hat", "nu1_hat", "kappa2_hat", "nu2_hat"
  )
  missing <- setdiff(required, names(sc))
  if (length(missing) > 0L) {
    stop(
      "The small-circle rolling CSV is missing required fitted-parameter columns: ",
      paste(missing, collapse = ", "),
      ". Save them in the S^2 rolling-window runner before running the axial diagnostic."
    )
  }

  if ("grid_group" %in% names(sc)) {
    sc <- sc[sc$grid_group == "all", , drop = FALSE]
  }

  sc$cycle <- as.integer(sc$cycle)
  sc$window_id <- as.integer(sc$window_id)
  sc$start_date <- as_utc(sc$start_date)
  sc$end_date <- as_utc(sc$end_date)
  if ("center_date" %in% names(sc)) {
    sc$center_date <- as_utc(sc$center_date)
  } else {
    sc$center_date <- as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")
  }

  sc <- sc[sc$cycle %in% cycles, , drop = FALSE]
  sc <- sc[order(sc$cycle, sc$window_id), , drop = FALSE]
  sc
}

read_axial_sunspots <- function(path_value, cycles) {
  candidate_paths <- trimws(unlist(strsplit(path_value, ",", fixed = TRUE), use.names = FALSE))
  candidate_paths <- candidate_paths[nzchar(candidate_paths)]
  if (length(candidate_paths) == 0L) stop("`sunspots_csv` cannot be empty.")
  missing_paths <- candidate_paths[!file.exists(candidate_paths)]
  if (length(missing_paths) > 0L) {
    stop("Sunspots CSV path(s) not found: ", paste(missing_paths, collapse = ", "))
  }

  data_list <- lapply(candidate_paths, function(one_path) {
    dat <- utils::read.csv(one_path, stringsAsFactors = FALSE)
    required <- c("cycle", "date", "x1", "x2", "x3")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
      stop("Sunspots CSV `", one_path, "` is missing columns: ", paste(missing, collapse = ", "))
    }
    dat$cycle <- as.integer(dat$cycle)
    dat$date <- as_utc(dat$date)
    dat
  })

  sunspots <- do.call(rbind, data_list)
  sunspots <- sunspots[sunspots$cycle %in% cycles, , drop = FALSE]
  sunspots[order(sunspots$cycle, sunspots$date), , drop = FALSE]
}

make_axial_ks_grid <- function(ks_omega_points = 60L, ks_t_points = 200L) {
  make_sample_unique_distance_ks_grid()
}

make_axial_start_theta <- function(sc_row) {
  normalize_axial_truncnorm_mixture2_theta(list(
    pi = as.numeric(sc_row$pi_hat),
    kappa1 = as.numeric(sc_row$kappa1_hat),
    nu1 = as.numeric(sc_row$nu1_hat),
    kappa2 = as.numeric(sc_row$kappa2_hat),
    nu2 = as.numeric(sc_row$nu2_hat)
  ))
}

axial_log_message <- function(log_path, message_text) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message_text), file = log_path, append = TRUE)
}

empty_axial_results <- function() {
  data.frame(
    cycle = integer(),
    window_id = integer(),
    start_date = character(),
    end_date = character(),
    center_date = character(),
    n = integer(),
    statistic_type = character(),
    test_statistic = double(),
    p_value_raw = double(),
    p_value_BY = double(),
    pi_axial = double(),
    kappa1_axial = double(),
    nu1_axial = double(),
    kappa2_axial = double(),
    nu2_axial = double(),
    mu1_sc = double(),
    mu2_sc = double(),
    mu3_sc = double(),
    pi_sc = double(),
    kappa1_sc = double(),
    nu1_sc = double(),
    kappa2_sc = double(),
    nu2_sc = double(),
    loglik_axial = double(),
    convergence_axial = integer(),
    elapsed_mle = double(),
    elapsed_bootstrap = double(),
    elapsed_total = double(),
    seed = integer(),
    mle_warnings = character(),
    bootstrap_warnings = character(),
    stringsAsFactors = FALSE
  )
}

read_existing_axial_results <- function(path) {
  if (!file.exists(path)) return(empty_axial_results())
  out <- utils::read.csv(path, stringsAsFactors = FALSE)
  out$cycle <- as.integer(out$cycle)
  out$window_id <- as.integer(out$window_id)
  out$p_value_raw <- as.numeric(out$p_value_raw)
  out$p_value_BY <- as.numeric(out$p_value_BY)
  out
}

recompute_axial_by <- function(results) {
  if (nrow(results) == 0L) return(results)
  results$p_value_BY <- NA_real_
  for (cycle_value in sort(unique(results$cycle))) {
    idx <- which(results$cycle == cycle_value)
    ok <- is.finite(results$p_value_raw[idx])
    if (any(ok)) {
      results$p_value_BY[idx[ok]] <- stats::p.adjust(results$p_value_raw[idx[ok]], method = "BY")
    }
  }
  results
}

write_axial_results <- function(results, path) {
  results <- recompute_axial_by(results)
  utils::write.csv(results, path, row.names = FALSE)
  results
}

compute_sc_by <- function(sc) {
  sc$p_SC_raw <- as.numeric(sc$p_value_raw)
  sc$p_SC_BY <- NA_real_
  for (cycle_value in sort(unique(sc$cycle))) {
    idx <- which(sc$cycle == cycle_value)
    ok <- is.finite(sc$p_SC_raw[idx])
    if (any(ok)) {
      sc$p_SC_BY[idx[ok]] <- stats::p.adjust(sc$p_SC_raw[idx[ok]], method = "BY")
    }
  }
  sc$neglog10_p_SC_raw <- -log10(pmax(sc$p_SC_raw, .Machine$double.xmin))
  sc$neglog10_p_SC_BY <- -log10(pmax(sc$p_SC_BY, .Machine$double.xmin))
  sc
}

make_axial_matched_table <- function(sc, axial_results) {
  sc_small <- compute_sc_by(sc[, c(
    "cycle", "window_id", "start_date", "end_date", "center_date", "p_value_raw"
  ), drop = FALSE])
  sc_small$start_date_chr <- format(sc_small$start_date, "%Y-%m-%d %H:%M:%S")
  sc_small$end_date_chr <- format(sc_small$end_date, "%Y-%m-%d %H:%M:%S")
  sc_small$center_date_chr <- format(sc_small$center_date, "%Y-%m-%d %H:%M:%S")

  axial_small <- axial_results
  axial_small$p_axial_raw <- as.numeric(axial_small$p_value_raw)
  axial_small$p_axial_BY <- as.numeric(axial_small$p_value_BY)
  axial_small$neglog10_p_axial_raw <- -log10(pmax(axial_small$p_axial_raw, .Machine$double.xmin))
  axial_small$neglog10_p_axial_BY <- -log10(pmax(axial_small$p_axial_BY, .Machine$double.xmin))

  matched <- merge(
    sc_small,
    axial_small[, c(
      "cycle", "window_id", "start_date", "end_date", "p_axial_raw", "p_axial_BY",
      "neglog10_p_axial_raw", "neglog10_p_axial_BY"
    ), drop = FALSE],
    by.x = c("cycle", "window_id", "start_date_chr", "end_date_chr"),
    by.y = c("cycle", "window_id", "start_date", "end_date"),
    all = FALSE
  )
  matched$center_date <- as_utc(matched$center_date_chr)
  matched[order(matched$cycle, matched$window_id), , drop = FALSE]
}

make_axial_correlations <- function(matched) {
  cycles <- c(sort(unique(matched$cycle)), NA_integer_)
  rows <- vector("list", 2L * length(cycles))
  k <- 0L
  for (cy in cycles) {
    df <- if (is.na(cy)) matched else matched[matched$cycle == cy, , drop = FALSE]
    cycle_label <- if (is.na(cy)) "all_cycles" else as.character(cy)

    k <- k + 1L
    rows[[k]] <- data.frame(
      cycle = cycle_label,
      comparison = "raw_vs_raw",
      n_matched = sum(is.finite(df$neglog10_p_SC_raw) & is.finite(df$neglog10_p_axial_raw)),
      cor_neglog10 = safe_cor(df$neglog10_p_SC_raw, df$neglog10_p_axial_raw),
      stringsAsFactors = FALSE
    )

    k <- k + 1L
    rows[[k]] <- data.frame(
      cycle = cycle_label,
      comparison = "BY_vs_BY",
      n_matched = sum(is.finite(df$neglog10_p_SC_BY) & is.finite(df$neglog10_p_axial_BY)),
      cor_neglog10 = safe_cor(df$neglog10_p_SC_BY, df$neglog10_p_axial_BY),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

plot_axial_timeplot <- function(matched, output_file, adjusted = FALSE) {
  png(output_file, width = 1800, height = 1200, res = 180)
  on.exit(dev.off(), add = TRUE)

  cycles <- sort(unique(matched$cycle))
  par(mfrow = c(length(cycles), 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  sc_col <- if (adjusted) "p_SC_BY" else "p_SC_raw"
  axial_col <- if (adjusted) "p_axial_BY" else "p_axial_raw"
  legend_labels <- if (adjusted) {
    c("Small circle (BY)", "Axial marginal (BY)", "0.05", "0.10")
  } else {
    c("Small circle (raw)", "Axial marginal (raw)", "0.05", "0.10")
  }

  for (cycle_value in cycles) {
    df <- matched[matched$cycle == cycle_value, , drop = FALSE]
    if (nrow(df) == 0L) next
    y_range <- range(c(df[[sc_col]], df[[axial_col]], 0.05, 0.10), finite = TRUE)
    plot(
      df$center_date,
      df[[sc_col]],
      type = "l",
      lwd = 2,
      col = plot_col_sc,
      ylim = y_range,
      xlab = "",
      ylab = "p-value",
      main = sprintf("Cycle %d", cycle_value)
    )
    lines(df$center_date, df[[axial_col]], lwd = 2, col = plot_col_axial)
    abline(h = 0.05, col = plot_col_ref, lty = 2)
    abline(h = 0.10, col = plot_col_ref, lty = 3)
    legend("topright", legend = legend_labels, col = c(plot_col_sc, plot_col_axial, plot_col_ref, plot_col_ref), lty = c(1, 1, 2, 3), lwd = c(2, 2, 1, 1), bty = "n")
  }

  mtext("Small-circle vs axial marginal rolling-window p-values", outer = TRUE, cex = 1.2)
}

run_axial_window <- function(sc_row, sunspots_cycle, args, log_path) {
  start_time_total <- Sys.time()
  in_window <- sunspots_cycle$date >= sc_row$start_date & sunspots_cycle$date < sc_row$end_date
  window_data <- sunspots_cycle[in_window, , drop = FALSE]
  if (nrow(window_data) != as.integer(sc_row$n)) {
    warning(sprintf(
      "Window (%d, %d) has n=%d in the rolling CSV but %d projected observations in the sunspots CSV subset.",
      sc_row$cycle, sc_row$window_id, sc_row$n, nrow(window_data)
    ))
  }

  mu_hat <- c(sc_row$mu_1, sc_row$mu_2, sc_row$mu_3)
  z <- as.numeric(as.matrix(window_data[, c("x1", "x2", "x3"), drop = FALSE]) %*% mu_hat)
  z <- pmin(1, pmax(-1, z))
  z <- normalize_axial_truncnorm_mixture2_data(z)

  start_theta <- make_axial_start_theta(sc_row)
  ks_grid <- make_axial_ks_grid(
    ks_omega_points = args$ks_omega_points,
    ks_t_points = args$ks_t_points
  )

  mle_warnings <- character()
  start_time_mle <- Sys.time()
  theta_hat <- withCallingHandlers(
    fit_axial_truncnorm_mixture2_theta(
      data = z,
      weights = NULL,
      null = list(type = "composite"),
      control = list(axial_truncnorm_mixture2_start_theta = start_theta)
    ),
    warning = function(w) {
      mle_warnings <<- c(mle_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed_mle <- as.numeric(difftime(Sys.time(), start_time_mle, units = "secs"))

  start_time_bootstrap <- Sys.time()
  bootstrap_result <- multiplier_bootstrap_axial_truncnorm_mixture2(
    data = z,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = ks_grid,
    B = args$B,
    n_cores = args$n_cores,
    seed = args$seed + as.integer(sc_row$cycle) * 100000L + as.integer(sc_row$window_id),
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      axial_truncnorm_mixture2_start_theta = theta_hat,
      axial_truncnorm_mixture2_bootstrap_optim_control = list(maxit = 80L, reltol = 1e-6)
    )
  )
  elapsed_bootstrap <- as.numeric(difftime(Sys.time(), start_time_bootstrap, units = "secs"))
  elapsed_total <- as.numeric(difftime(Sys.time(), start_time_total, units = "secs"))

  bootstrap_warnings <- unlist(lapply(bootstrap_result$bootstrap$theta_star %||% list(), function(x) x$warnings %||% NULL), use.names = FALSE)
  if (length(bootstrap_warnings) == 0L && !is.null(bootstrap_result$observed$theta_hat$opt$message)) {
    bootstrap_warnings <- character()
  }

  data.frame(
    cycle = as.integer(sc_row$cycle),
    window_id = as.integer(sc_row$window_id),
    start_date = format(sc_row$start_date, "%Y-%m-%d %H:%M:%S"),
    end_date = format(sc_row$end_date, "%Y-%m-%d %H:%M:%S"),
    center_date = format(sc_row$center_date, "%Y-%m-%d %H:%M:%S"),
    n = length(z),
    statistic_type = "KS_axial",
    test_statistic = as.numeric(bootstrap_result$inference$ks$observed),
    p_value_raw = as.numeric(bootstrap_result$inference$ks$p_value),
    p_value_BY = NA_real_,
    pi_axial = as.numeric(theta_hat$pi),
    kappa1_axial = as.numeric(theta_hat$kappa1),
    nu1_axial = as.numeric(theta_hat$nu1),
    kappa2_axial = as.numeric(theta_hat$kappa2),
    nu2_axial = as.numeric(theta_hat$nu2),
    mu1_sc = as.numeric(sc_row$mu_1),
    mu2_sc = as.numeric(sc_row$mu_2),
    mu3_sc = as.numeric(sc_row$mu_3),
    pi_sc = as.numeric(sc_row$pi_hat),
    kappa1_sc = as.numeric(sc_row$kappa1_hat),
    nu1_sc = as.numeric(sc_row$nu1_hat),
    kappa2_sc = as.numeric(sc_row$kappa2_hat),
    nu2_sc = as.numeric(sc_row$nu2_hat),
    loglik_axial = as.numeric(theta_hat$loglik),
    convergence_axial = as.integer(theta_hat$opt$convergence %||% NA_integer_),
    elapsed_mle = elapsed_mle,
    elapsed_bootstrap = elapsed_bootstrap,
    elapsed_total = elapsed_total,
    seed = as.integer(args$seed),
    mle_warnings = paste(unique(mle_warnings), collapse = " | "),
    bootstrap_warnings = paste(unique(bootstrap_warnings), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_axial_args()
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  log_path <- file.path(args$output_dir, "sunspots_axial_marginal_rolling_10cr.log")
  results_csv <- file.path(args$output_dir, "axial_marginal_rolling_10cr_gof_results.csv")
  matched_csv <- file.path(args$output_dir, "sc_vs_axial_marginal_window_matched_cycles21_23.csv")
  correlations_csv <- file.path(args$output_dir, "sc_vs_axial_marginal_correlations.csv")

  sc <- read_axial_sc_results(args$sc_csv, cycles = args$cycles)
  sunspots <- read_axial_sunspots(args$sunspots_csv, cycles = args$cycles)

  selected <- sc[
    sc$window_id >= args$start_window &
      (is.infinite(args$end_window) | sc$window_id <= args$end_window),
    ,
    drop = FALSE
  ]

  results <- if (isTRUE(args$resume)) read_existing_axial_results(results_csv) else empty_axial_results()
  completed_keys <- paste(results$cycle, results$window_id, results$statistic_type, sep = "::")

  axial_log_message(log_path, sprintf("Starting axial rolling-window GOF run for cycles: %s", paste(args$cycles, collapse = ",")))

  for (i in seq_len(nrow(selected))) {
    sc_row <- selected[i, , drop = FALSE]
    row_key <- paste(sc_row$cycle, sc_row$window_id, "KS_axial", sep = "::")
    if (row_key %in% completed_keys) {
      next
    }

    axial_log_message(log_path, sprintf("Processing cycle %d window %d", sc_row$cycle, sc_row$window_id))
    sunspots_cycle <- sunspots[sunspots$cycle == sc_row$cycle, , drop = FALSE]
    result_row <- run_axial_window(sc_row = sc_row, sunspots_cycle = sunspots_cycle, args = args, log_path = log_path)
    results <- rbind(results, result_row)
    results <- results[order(results$cycle, results$window_id), , drop = FALSE]
    results <- write_axial_results(results, results_csv)
    completed_keys <- paste(results$cycle, results$window_id, results$statistic_type, sep = "::")
  }

  results <- write_axial_results(results, results_csv)

  matched <- make_axial_matched_table(sc = sc, axial_results = results)
  utils::write.csv(matched, matched_csv, row.names = FALSE)

  correlations <- make_axial_correlations(matched)
  utils::write.csv(correlations, correlations_csv, row.names = FALSE)

  plot_axial_timeplot(
    matched = matched,
    output_file = file.path(args$output_dir, "sc_vs_axial_marginal_pvalue_timeplot_raw.png"),
    adjusted = FALSE
  )
  plot_axial_timeplot(
    matched = matched,
    output_file = file.path(args$output_dir, "sc_vs_axial_marginal_pvalue_timeplot_BY.png"),
    adjusted = TRUE
  )

  writeLines(capture.output(utils::sessionInfo()), con = file.path(args$output_dir, "sessionInfo_axial_marginal.txt"))
  axial_log_message(log_path, "Finished axial rolling-window GOF run.")
}

main()
