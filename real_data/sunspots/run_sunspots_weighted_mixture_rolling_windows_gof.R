#!/usr/bin/env Rscript

resolve_sunspots_weighted_rolling_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

utils_path_sunspots_weighted_rolling <- resolve_sunspots_weighted_rolling_path("utils.R")
model_specs_path_sunspots_weighted_rolling <- resolve_sunspots_weighted_rolling_path("bootstrap", "model_specs.R")
bootstrap_path_sunspots_weighted_rolling <- resolve_sunspots_weighted_rolling_path("bootstrap", "multiplier_bootstrap.R")
spec_path_sunspots_weighted_rolling <- resolve_sunspots_weighted_rolling_path("bootstrap", "small_circle_weighted_mixture2_model_spec.R")
prep_path_sunspots_weighted_rolling <- resolve_sunspots_weighted_rolling_path("real_data", "sunspots", "sunspots.R")
fmgp_rdata_dir_sunspots_weighted_rolling <- resolve_sunspots_weighted_rolling_path("real_data", "sunspots", "fmgp_rdata")

source(utils_path_sunspots_weighted_rolling)
source(model_specs_path_sunspots_weighted_rolling)
source(spec_path_sunspots_weighted_rolling)
source(bootstrap_path_sunspots_weighted_rolling)

default_weighted_cycle24_theta_start <- function() {
  small_circle_weighted_mixture2_normalize_theta(list(
    pi = 0.529731,
    mu = c(-0.001131, -0.004108, 0.999991),
    kappa1 = 26.806931,
    nu1 = 0.237492,
    kappa2 = 24.109693,
    nu2 = 0.268405
  ))
}

normalize_weighted_rolling_statistics <- function(statistics) {
  statistics <- unique(tolower(trimws(unlist(strsplit(paste(statistics, collapse = ","), ",", fixed = TRUE), use.names = FALSE))))
  statistics <- statistics[nzchar(statistics)]
  if (length(statistics) == 0L) stop("At least one statistic must be supplied.")
  if (!all(statistics %in% c("ks", "cvm"))) stop("The weighted mixture rolling-window runner currently supports only KS and CvM.")
  statistics
}

parse_weighted_rolling_bool <- function(value, key) {
  value_lower <- tolower(trimws(value))
  if (value_lower %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value_lower %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("Option '--%s' must be TRUE or FALSE.", key))
}

parse_weighted_rolling_integer <- function(value, key, allow_inf = FALSE) {
  if (allow_inf && tolower(trimws(value)) %in% c("inf", "infinity")) return(Inf)
  if (!grepl("^-?[0-9]+$", value)) stop(sprintf("Option '--%s' must be an integer.", key))
  as.integer(value)
}

parse_weighted_rolling_cycles <- function(value) {
  parts <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) stop("Option '--cycles' cannot be empty.")
  if (any(!grepl("^-?[0-9]+$", parts))) stop("Option '--cycles' must contain comma-separated integer cycle identifiers.")
  unique(as.integer(parts))
}

parse_weighted_rolling_grid_groups <- function(value) {
  parts <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts)]
  supported <- c("all", "N", "S")
  if (length(parts) == 0L) stop("Option '--grid_groups' cannot be empty.")
  if (any(!parts %in% supported)) stop("Option '--grid_groups' must use values in {all,N,S}.")
  unique(parts)
}

parse_sunspots_weighted_rolling_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  help_text <- paste(
    "Usage:",
    "  Rscript real_data/sunspots/run_sunspots_weighted_mixture_rolling_windows_gof.R [options]",
    "",
    "Options:",
    "  --help",
    "  --cycles=21,22,23",
    "  --grid_groups=all",
    "  --input_csv=PATH",
    "  --output_dir=PATH",
    "  --B=1000",
    "  --n_cores=12",
    "  --seed=123",
    "  --distance_type=geodesic",
    "  --window_rotations=10",
    "  --step_rotations=1",
    "  --min_n=50",
    "  --start_window=1",
    "  --end_window=Inf",
    "  --resume=TRUE",
    "  --windows_only=FALSE",
    "  --statistics=KS|CvM|KS,CvM",
    sep = "\n"
  )
  if (any(args %in% c("--help", "-h"))) {
    cat(help_text, sep = "\n")
    quit(save = "no", status = 0L)
  }

  parsed <- list()
  for (arg in args) {
    key_value <- sub("^--", "", arg)
    parts <- strsplit(key_value, "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else NULL
    if (is.null(value)) stop(sprintf("Option '--%s' requires the form --%s=value.", key, key))
    parsed[[key]] <- value
  }

  out <- list()
  if ("cycles" %in% names(parsed)) out$cycles <- parse_weighted_rolling_cycles(parsed$cycles)
  if ("grid_groups" %in% names(parsed)) out$grid_groups <- parse_weighted_rolling_grid_groups(parsed$grid_groups)
  if ("input_csv" %in% names(parsed)) out$input_csv <- parsed$input_csv
  if ("output_dir" %in% names(parsed)) out$output_dir <- parsed$output_dir
  if ("B" %in% names(parsed)) out$B <- parse_weighted_rolling_integer(parsed$B, "B")
  if ("n_cores" %in% names(parsed)) out$n_cores <- parse_weighted_rolling_integer(parsed$n_cores, "n_cores")
  if ("seed" %in% names(parsed)) out$seed <- parse_weighted_rolling_integer(parsed$seed, "seed")
  if ("distance_type" %in% names(parsed)) out$distance_type <- parsed$distance_type
  if ("window_rotations" %in% names(parsed)) out$window_rotations <- as.numeric(parsed$window_rotations)
  if ("step_rotations" %in% names(parsed)) out$step_rotations <- as.numeric(parsed$step_rotations)
  if ("min_n" %in% names(parsed)) out$min_n <- parse_weighted_rolling_integer(parsed$min_n, "min_n")
  if ("start_window" %in% names(parsed)) out$start_window <- parse_weighted_rolling_integer(parsed$start_window, "start_window")
  if ("end_window" %in% names(parsed)) out$end_window <- parse_weighted_rolling_integer(parsed$end_window, "end_window", allow_inf = TRUE)
  if ("resume" %in% names(parsed)) out$resume <- parse_weighted_rolling_bool(parsed$resume, "resume")
  if ("windows_only" %in% names(parsed)) out$windows_only <- parse_weighted_rolling_bool(parsed$windows_only, "windows_only")
  if ("statistics" %in% names(parsed)) out$statistics <- normalize_weighted_rolling_statistics(parsed$statistics)
  if ("statistic" %in% names(parsed)) out$statistics <- normalize_weighted_rolling_statistics(parsed$statistic)
  out
}

sunspots_weighted_rolling_log_message <- function(log_path, message_text) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message_text), file = log_path, append = TRUE)
}

build_cycle_specific_rolling_windows <- function(df,
                                                 cycle,
                                                 grid_group = "all",
                                                 window_rotations = 10,
                                                 step_rotations = 1,
                                                 carrington_days = 27.2753,
                                                 min_n = 50L,
                                                 global_window_start = 1L) {
  cycle_df_all <- df[df$cycle == cycle, , drop = FALSE]
  if (grid_group == "N") {
    cycle_df_grid <- cycle_df_all[cycle_df_all$hemisphere == "N", , drop = FALSE]
  } else if (grid_group == "S") {
    cycle_df_grid <- cycle_df_all[cycle_df_all$hemisphere == "S", , drop = FALSE]
  } else {
    cycle_df_grid <- cycle_df_all
  }

  if (nrow(cycle_df_all) == 0L || nrow(cycle_df_grid) == 0L) {
    empty_windows <- data.frame(
      cycle = integer(),
      grid_group = character(),
      window_id = integer(),
      global_window_id = integer(),
      start_date = character(),
      end_date = character(),
      center_date = character(),
      n = integer(),
      n_north = integer(),
      n_south = integer(),
      first_grid_date = character(),
      last_grid_date = character(),
      first_cycle_date_all = character(),
      last_cycle_date_all = character(),
      window_rotations = double(),
      step_rotations = double(),
      carrington_days = double(),
      window_complete_by_fmgp_rule = logical(),
      keep_window = logical(),
      stringsAsFactors = FALSE
    )
    return(list(cycle_data_all = cycle_df_all, windows = empty_windows, next_global_window_id = global_window_start))
  }

  dates_all <- as.POSIXct(cycle_df_all$date, tz = "UTC")
  dates_grid <- as.POSIXct(cycle_df_grid$date, tz = "UTC")
  cycle_df_all <- cycle_df_all[order(dates_all), , drop = FALSE]
  dates_all <- as.POSIXct(cycle_df_all$date, tz = "UTC")
  cycle_df_grid <- cycle_df_grid[order(dates_grid), , drop = FALSE]
  dates_grid <- as.POSIXct(cycle_df_grid$date, tz = "UTC")

  window_seconds <- as.numeric(window_rotations) * carrington_days * 24 * 60 * 60
  step_seconds <- as.numeric(step_rotations) * carrington_days * 24 * 60 * 60
  first_time <- as.numeric(min(dates_grid))
  last_time <- as.numeric(max(dates_grid))
  total_windows <- floor((last_time - (first_time + window_seconds)) / step_seconds)

  first_grid_date <- format(as.POSIXct(first_time, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H:%M:%S")
  last_grid_date <- format(as.POSIXct(last_time, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H:%M:%S")
  first_cycle_date_all <- format(min(dates_all), "%Y-%m-%d %H:%M:%S")
  last_cycle_date_all <- format(max(dates_all), "%Y-%m-%d %H:%M:%S")

  if (!is.finite(total_windows) || total_windows < 1L) {
    empty_windows <- data.frame(
      cycle = integer(),
      grid_group = character(),
      window_id = integer(),
      global_window_id = integer(),
      start_date = character(),
      end_date = character(),
      center_date = character(),
      n = integer(),
      n_north = integer(),
      n_south = integer(),
      first_grid_date = character(),
      last_grid_date = character(),
      first_cycle_date_all = character(),
      last_cycle_date_all = character(),
      window_rotations = double(),
      step_rotations = double(),
      carrington_days = double(),
      window_complete_by_fmgp_rule = logical(),
      keep_window = logical(),
      stringsAsFactors = FALSE
    )
    return(list(cycle_data_all = cycle_df_all, windows = empty_windows, next_global_window_id = global_window_start))
  }

  start_times <- first_time + seq.int(0L, total_windows - 1L) * step_seconds
  window_list <- vector("list", length(start_times))
  for (i in seq_along(start_times)) {
    start_time <- as.POSIXct(start_times[[i]], origin = "1970-01-01", tz = "UTC")
    end_time <- as.POSIXct(start_times[[i]] + window_seconds, origin = "1970-01-01", tz = "UTC")
    in_window_all <- dates_all >= start_time & dates_all < end_time
    window_df_all <- cycle_df_all[in_window_all, , drop = FALSE]
    n_window <- nrow(window_df_all)
    window_list[[i]] <- data.frame(
      cycle = cycle,
      grid_group = grid_group,
      window_id = i,
      global_window_id = global_window_start + i - 1L,
      start_date = format(start_time, "%Y-%m-%d %H:%M:%S"),
      end_date = format(end_time, "%Y-%m-%d %H:%M:%S"),
      center_date = format(as.POSIXct(start_times[[i]] + window_seconds / 2, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H:%M:%S"),
      n = n_window,
      n_north = sum(window_df_all$hemisphere == "N"),
      n_south = sum(window_df_all$hemisphere == "S"),
      first_grid_date = first_grid_date,
      last_grid_date = last_grid_date,
      first_cycle_date_all = first_cycle_date_all,
      last_cycle_date_all = last_cycle_date_all,
      window_rotations = as.numeric(window_rotations),
      step_rotations = as.numeric(step_rotations),
      carrington_days = as.numeric(carrington_days),
      window_complete_by_fmgp_rule = TRUE,
      keep_window = n_window >= as.integer(min_n),
      stringsAsFactors = FALSE
    )
  }

  list(
    cycle_data_all = cycle_df_all,
    windows = do.call(rbind, window_list),
    next_global_window_id = global_window_start + length(start_times)
  )
}

build_weighted_rolling_windows_all_cycles <- function(df,
                                                      cycles,
                                                      grid_groups = "all",
                                                      window_rotations = 10,
                                                      step_rotations = 1,
                                                      carrington_days = 27.2753,
                                                      min_n = 50L) {
  next_global_window_id <- 1L
  cycle_data_list <- list()
  windows_list <- list()
  counter <- 1L
  for (cycle in cycles) {
    for (grid_group in grid_groups) {
      out <- build_cycle_specific_rolling_windows(
        df = df,
        cycle = cycle,
        grid_group = grid_group,
        window_rotations = window_rotations,
        step_rotations = step_rotations,
        carrington_days = carrington_days,
        min_n = min_n,
        global_window_start = next_global_window_id
      )
      cycle_data_list[[paste(cycle, grid_group, sep = "::")]] <- out$cycle_data_all
      windows_list[[counter]] <- out$windows
      next_global_window_id <- out$next_global_window_id
      counter <- counter + 1L
    }
  }
  windows_df <- if (length(windows_list) == 0L) data.frame() else do.call(rbind, windows_list)
  list(cycle_data = cycle_data_list, windows = windows_df)
}

theta_from_weighted_results_row <- function(row_df) {
  small_circle_weighted_mixture2_normalize_theta(list(
    pi = row_df$pi_hat[[1L]],
    mu = c(row_df$mu_1[[1L]], row_df$mu_2[[1L]], row_df$mu_3[[1L]]),
    kappa1 = row_df$kappa1_hat[[1L]],
    nu1 = row_df$nu1_hat[[1L]],
    kappa2 = row_df$kappa2_hat[[1L]],
    nu2 = row_df$nu2_hat[[1L]]
  ))
}

lookup_cycle_full_theta_start <- function(cycle, output_root = file.path("real_data", "sunspots", "output")) {
  if (cycle == 24L) return(default_weighted_cycle24_theta_start())
  theta_path <- file.path(output_root, sprintf("cycle%d_small_circle_weighted_mixture", cycle), "sunspots_cycle24_theta_hat.csv")
  if (!file.exists(theta_path)) return(NULL)
  theta_df <- utils::read.csv(theta_path, stringsAsFactors = FALSE)
  theta_from_weighted_results_row(theta_df[1L, , drop = FALSE])
}

load_fmgp_alignment_table <- function(cycles, grid_groups, alignment_output_path) {
  group_to_tag <- c(all = "all", N = "N", S = "S")
  rows <- list()
  counter <- 1L
  for (cycle in cycles) {
    for (grid_group in grid_groups) {
      fpath <- file.path(
        fmgp_rdata_dir_sunspots_weighted_rolling,
        sprintf("sunspots_cir_test_%s_10_results_cycle%d.RData", group_to_tag[[grid_group]], cycle)
      )
      if (!file.exists(fpath)) next
      e <- new.env(parent = emptyenv())
      load(fpath, envir = e)
      if (!exists("test.results", envir = e, inherits = FALSE)) next
      test_results <- get("test.results", envir = e, inherits = FALSE)
      if (!is.data.frame(test_results) || !"date" %in% names(test_results)) next
      rows[[counter]] <- data.frame(
        cycle = cycle,
        grid_group = grid_group,
        window_id = seq_len(nrow(test_results)),
        fmgp_date = as.character(test_results$date),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

compute_fmgp_alignment <- function(windows_df, cycles, grid_groups, output_dir) {
  fmgp_df <- load_fmgp_alignment_table(cycles = cycles, grid_groups = grid_groups, alignment_output_path = output_dir)
  if (is.null(fmgp_df) || nrow(windows_df) == 0L) return(NULL)

  merged <- merge(
    windows_df[, c("cycle", "grid_group", "window_id", "end_date")],
    fmgp_df,
    by = c("cycle", "grid_group", "window_id"),
    all = TRUE,
    sort = TRUE
  )
  merged$our_end_date <- merged$end_date
  merged$diff_days <- as.numeric(as.POSIXct(merged$our_end_date, tz = "UTC") - as.POSIXct(merged$fmgp_date, tz = "UTC")) / (24 * 60 * 60)
  merged <- merged[, c("cycle", "grid_group", "window_id", "our_end_date", "fmgp_date", "diff_days")]
  utils::write.csv(
    merged,
    file.path(output_dir, "sunspots_cycles21_23_rolling_10cr_windows_fmgp_alignment.csv"),
    row.names = FALSE
  )

  summary_rows <- do.call(rbind, lapply(split(merged, list(merged$cycle, merged$grid_group), drop = TRUE), function(d) {
    data.frame(
      cycle = d$cycle[[1L]],
      grid_group = d$grid_group[[1L]],
      n_our = sum(!is.na(d$our_end_date)),
      n_fmgp = sum(!is.na(d$fmgp_date)),
      max_abs_diff_days = if (all(is.na(d$diff_days))) NA_real_ else max(abs(d$diff_days), na.rm = TRUE),
      median_abs_diff_days = if (all(is.na(d$diff_days))) NA_real_ else stats::median(abs(d$diff_days), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(
    summary_rows,
    file.path(output_dir, "sunspots_cycles21_23_rolling_10cr_windows_fmgp_alignment_summary.csv"),
    row.names = FALSE
  )
  list(alignment = merged, summary = summary_rows)
}

run_weighted_window_gof_once <- function(window_df,
                                         window_meta,
                                         statistics = "ks",
                                         B = 1000L,
                                         n_cores = 12L,
                                         bootstrap_method = "reestimated",
                                         seed = 123L,
                                         distance_type = "geodesic",
                                         theta_start = NULL) {
  x <- as.matrix(window_df[, c("x1", "x2", "x3")])
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  control <- list(
    small_circle_weighted_mixture2_profile_method = "legendre",
    small_circle_weighted_mixture2_L_max = 200L,
    small_circle_weighted_mixture2_quad_n = 400L,
    small_circle_weighted_mixture2_tol = 1e-10,
    small_circle_weighted_mixture2_optim_control = list(maxit = 400L, reltol = 1e-9),
    small_circle_weighted_mixture2_n_starts = 8L
  )
  if (!is.null(theta_start)) {
    control$small_circle_weighted_mixture2_start_theta <- theta_start
    control$small_circle_weighted_mixture2_warm_start_only <- TRUE
  }

  mle_warning_log <- character()
  mle_start <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    small_circle_weighted_mixture2_mle_s2_weighted(x = x, control = control),
    warning = function(w) {
      mle_warning_log <<- c(mle_warning_log, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed_mle <- proc.time()[["elapsed"]] - mle_start
  theta_hat <- small_circle_weighted_mixture2_normalize_theta(list(
    pi = fit$pi, mu = fit$mu, kappa1 = fit$kappa1, nu1 = fit$nu1, kappa2 = fit$kappa2, nu2 = fit$nu2
  ))
  loglik <- fit$loglik
  n <- nrow(x)
  aic <- 2 * 7L - 2 * loglik
  bic <- log(n) * 7L - 2 * loglik

  statistics <- normalize_weighted_rolling_statistics(statistics)
  ks_grid <- list(omega_grid = generate_canonical_lattice(60L, dim = 3), t_grid = seq(1e-8, pi - 1e-8, length.out = 200L))
  spec <- make_small_circle_weighted_mixture2_spec(distance_type = distance_type)
  bootstrap_warning_log <- character()
  bootstrap_start <- proc.time()[["elapsed"]]
  gof_result <- withCallingHandlers(
    multiplier_bootstrap_gof(
      data = x,
      spec = spec,
      null = list(type = "composite"),
      statistics = statistics,
      ks_grid = if ("ks" %in% statistics) ks_grid else NULL,
      B = as.integer(B),
      alpha = 0.05,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      bootstrap_method = bootstrap_method,
      observed_theta_hat = theta_hat,
      keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE),
      control = control
    ),
    warning = function(w) {
      bootstrap_warning_log <<- c(bootstrap_warning_log, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed_bootstrap <- proc.time()[["elapsed"]] - bootstrap_start
  elapsed_total <- elapsed_mle + elapsed_bootstrap

  list(
    result_row = do.call(rbind, lapply(statistics, function(stat_name) {
      data.frame(
        cycle = window_meta$cycle[[1L]],
        grid_group = window_meta$grid_group[[1L]],
        window_id = window_meta$window_id[[1L]],
        global_window_id = window_meta$global_window_id[[1L]],
        start_date = window_meta$start_date[[1L]],
        end_date = window_meta$end_date[[1L]],
        center_date = window_meta$center_date[[1L]],
        n = n,
        n_north = sum(window_df$hemisphere == "N"),
        n_south = sum(window_df$hemisphere == "S"),
        model = "small_circle_weighted_mixture2",
        statistic_type = stat_name,
        test_statistic = gof_result$inference[[stat_name]]$observed,
        p_value_raw = gof_result$inference[[stat_name]]$p_value,
        loglik = loglik,
        aic = aic,
        bic = bic,
        pi_hat = theta_hat$pi,
        mu_1 = theta_hat$mu[[1L]],
        mu_2 = theta_hat$mu[[2L]],
        mu_3 = theta_hat$mu[[3L]],
        kappa1_hat = theta_hat$kappa1,
        nu1_hat = theta_hat$nu1,
        kappa2_hat = theta_hat$kappa2,
        nu2_hat = theta_hat$nu2,
        convergence = fit$opt$convergence %||% NA_integer_,
        elapsed_mle = elapsed_mle,
        elapsed_bootstrap = elapsed_bootstrap,
        elapsed_total = elapsed_total,
        seed = as.integer(seed),
        mle_warnings = paste(unique(mle_warning_log), collapse = " | "),
        bootstrap_warnings = paste(unique(bootstrap_warning_log), collapse = " | "),
        stringsAsFactors = FALSE
      )
    })),
    theta_hat = theta_hat
  )
}

run_sunspots_weighted_mixture_rolling_windows_gof <- function(
    cycles = c(21L, 22L, 23L),
    grid_groups = "all",
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycles21_23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycles21_23_weighted_mixture_rolling_windows_10cr_B1000"),
    statistics = "ks",
    B = 1000L,
    n_cores = 12L,
    bootstrap_method = "reestimated",
    seed = 123L,
    distance_type = "geodesic",
    window_rotations = 10,
    step_rotations = 1,
    min_n = 50L,
    start_window = 1L,
    end_window = Inf,
    resume = TRUE,
    windows_only = FALSE,
    carrington_days = 27.2753) {
  cycles <- unique(as.integer(cycles))
  grid_groups <- unique(as.character(grid_groups))
  statistics <- normalize_weighted_rolling_statistics(statistics)
  if (!file.exists(input_csv)) source(prep_path_sunspots_weighted_rolling)
  if (!file.exists(input_csv)) stop(sprintf("Input CSV not found: %s", input_csv))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sessionInfo.txt"))
  log_path <- file.path(output_dir, "sunspots_cycles21_23_weighted_mixture_rolling_10cr.log")
  start_elapsed <- proc.time()[["elapsed"]]
  sunspots_weighted_rolling_log_message(log_path, sprintf("Starting rolling-window GOF run for cycles=%s and grid_groups=%s", paste(cycles, collapse = ","), paste(grid_groups, collapse = ",")))

  df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  df <- df[df$cycle %in% cycles, , drop = FALSE]
  windows_out <- build_weighted_rolling_windows_all_cycles(
    df = df,
    cycles = cycles,
    grid_groups = grid_groups,
    window_rotations = window_rotations,
    step_rotations = step_rotations,
    carrington_days = carrington_days,
    min_n = min_n
  )
  windows_df <- windows_out$windows
  if (nrow(windows_df) > 0L) {
    windows_df$selected_window <- windows_df$window_id >= as.integer(start_window) &
      windows_df$window_id <= if (is.infinite(end_window)) Inf else as.integer(end_window)
  }
  windows_path <- file.path(output_dir, "sunspots_cycles21_23_rolling_10cr_windows.csv")
  utils::write.csv(windows_df, windows_path, row.names = FALSE)
  alignment_out <- compute_fmgp_alignment(windows_df = windows_df, cycles = cycles, grid_groups = grid_groups, output_dir = output_dir)

  if (isTRUE(windows_only)) {
    sunspots_weighted_rolling_log_message(log_path, "windows_only=TRUE, stopping after window construction and FMGP alignment.")
    return(list(windows = windows_df, alignment = alignment_out, results = NULL, elapsed_total = proc.time()[["elapsed"]] - start_elapsed))
  }

  results_path <- file.path(output_dir, "sunspots_cycles21_23_weighted_mixture_rolling_10cr_gof_results.csv")
  if (file.exists(results_path)) {
    results_df <- utils::read.csv(results_path, stringsAsFactors = FALSE)
  } else {
    results_df <- data.frame(
      cycle = integer(), grid_group = character(), window_id = integer(), global_window_id = integer(),
      start_date = character(), end_date = character(), center_date = character(),
      n = integer(), n_north = integer(), n_south = integer(),
      model = character(), statistic_type = character(), test_statistic = double(), p_value_raw = double(),
      loglik = double(), aic = double(), bic = double(), pi_hat = double(),
      mu_1 = double(), mu_2 = double(), mu_3 = double(),
      kappa1_hat = double(), nu1_hat = double(), kappa2_hat = double(), nu2_hat = double(),
      convergence = integer(), elapsed_mle = double(), elapsed_bootstrap = double(), elapsed_total = double(),
      seed = integer(), mle_warnings = character(), bootstrap_warnings = character(), stringsAsFactors = FALSE
    )
  }

  selected_windows <- windows_df[windows_df$selected_window, , drop = FALSE]
  cycle_grid_theta_start <- list()
  for (cycle in cycles) {
    for (grid_group in grid_groups) {
      key <- paste(cycle, grid_group, sep = "::")
      completed <- results_df[
        results_df$cycle == cycle &
          results_df$grid_group == grid_group &
          results_df$statistic_type == statistics[[1L]] &
          is.finite(results_df$p_value_raw),
        ,
        drop = FALSE
      ]
      if (nrow(completed) > 0L) {
        completed <- completed[order(completed$window_id), , drop = FALSE]
        cycle_grid_theta_start[[key]] <- theta_from_weighted_results_row(completed[nrow(completed), , drop = FALSE])
      } else {
        cycle_grid_theta_start[[key]] <- lookup_cycle_full_theta_start(cycle)
      }
    }
  }

  for (i in seq_len(nrow(selected_windows))) {
    window_meta <- selected_windows[i, , drop = FALSE]
    if (!isTRUE(window_meta$keep_window[[1L]])) {
      sunspots_weighted_rolling_log_message(log_path, sprintf("Skipping cycle %d group %s window %d because n=%d < min_n=%d.", window_meta$cycle[[1L]], window_meta$grid_group[[1L]], window_meta$window_id[[1L]], window_meta$n[[1L]], as.integer(min_n)))
      next
    }

    key <- paste(window_meta$cycle[[1L]], window_meta$grid_group[[1L]], sep = "::")
    existing <- results_df[
      results_df$cycle == window_meta$cycle[[1L]] &
        results_df$grid_group == window_meta$grid_group[[1L]] &
        results_df$window_id == window_meta$window_id[[1L]] &
        is.finite(results_df$p_value_raw),
      ,
      drop = FALSE
    ]
    missing_statistics <- setdiff(statistics, existing$statistic_type)
    if (isTRUE(resume) && length(missing_statistics) == 0L) {
      sunspots_weighted_rolling_log_message(log_path, sprintf("Skipping cycle %d group %s window %d because requested statistics are already completed.", window_meta$cycle[[1L]], window_meta$grid_group[[1L]], window_meta$window_id[[1L]]))
      existing <- existing[order(existing$window_id), , drop = FALSE]
      cycle_grid_theta_start[[key]] <- theta_from_weighted_results_row(existing[1L, , drop = FALSE])
      next
    }

    cycle_df_all <- windows_out$cycle_data[[key]]
    cycle_dates_all <- as.POSIXct(cycle_df_all$date, tz = "UTC")
    start_time <- as.POSIXct(window_meta$start_date[[1L]], tz = "UTC")
    end_time <- as.POSIXct(window_meta$end_date[[1L]], tz = "UTC")
    in_window <- cycle_dates_all >= start_time & cycle_dates_all < end_time
    window_df <- cycle_df_all[in_window, , drop = FALSE]
    window_seed <- as.integer(seed) + as.integer(window_meta$global_window_id[[1L]])

    sunspots_weighted_rolling_log_message(log_path, sprintf("Cycle %d group %s window %d (global %d): [%s, %s) | n=%d | seed=%d", window_meta$cycle[[1L]], window_meta$grid_group[[1L]], window_meta$window_id[[1L]], window_meta$global_window_id[[1L]], window_meta$start_date[[1L]], window_meta$end_date[[1L]], nrow(window_df), window_seed))
    result <- try(
      run_weighted_window_gof_once(
        window_df = window_df,
        window_meta = window_meta,
        statistics = if (isTRUE(resume)) missing_statistics else statistics,
        B = as.integer(B),
        n_cores = as.integer(n_cores),
        bootstrap_method = bootstrap_method,
        seed = window_seed,
        distance_type = distance_type,
        theta_start = cycle_grid_theta_start[[key]]
      ),
      silent = TRUE
    )
    if (inherits(result, "try-error")) {
      sunspots_weighted_rolling_log_message(log_path, sprintf("Cycle %d group %s window %d failed: %s", window_meta$cycle[[1L]], window_meta$grid_group[[1L]], window_meta$window_id[[1L]], as.character(result)))
      next
    }

    duplicate_rows <- results_df$cycle == window_meta$cycle[[1L]] &
      results_df$grid_group == window_meta$grid_group[[1L]] &
      results_df$window_id == window_meta$window_id[[1L]] &
      results_df$statistic_type %in% result$result_row$statistic_type
    if (any(duplicate_rows)) results_df <- results_df[!duplicate_rows, , drop = FALSE]
    results_df <- rbind(results_df, result$result_row)
    results_df <- results_df[order(results_df$cycle, results_df$grid_group, results_df$window_id, results_df$statistic_type), , drop = FALSE]
    utils::write.csv(results_df, results_path, row.names = FALSE)
    cycle_grid_theta_start[[key]] <- result$theta_hat
    sunspots_weighted_rolling_log_message(log_path, sprintf("Cycle %d group %s window %d completed: KS=%.6f | p-value=%.6f | convergence=%s", window_meta$cycle[[1L]], window_meta$grid_group[[1L]], window_meta$window_id[[1L]], result$result_row$test_statistic[[1L]], result$result_row$p_value_raw[[1L]], as.character(result$result_row$convergence[[1L]])))
  }

  summary_rows <- do.call(rbind, lapply(statistics, function(stat_name) {
    do.call(rbind, lapply(split(windows_df, list(windows_df$cycle, windows_df$grid_group), drop = TRUE), function(wg) {
      cycle_value <- wg$cycle[[1L]]
      grid_group_value <- wg$grid_group[[1L]]
      completed <- results_df[
        results_df$cycle == cycle_value &
          results_df$grid_group == grid_group_value &
          results_df$statistic_type == stat_name &
          results_df$window_id %in% wg$window_id[wg$keep_window] &
          is.finite(results_df$p_value_raw),
        ,
        drop = FALSE
      ]
      data.frame(
        cycle = cycle_value,
        grid_group = grid_group_value,
        statistic_type = stat_name,
        total_windows = nrow(wg),
        windows_below_min_n = sum(!wg$keep_window),
        completed_windows = nrow(completed),
        min_p_value = if (nrow(completed) > 0L) min(completed$p_value_raw) else NA_real_,
        median_p_value = if (nrow(completed) > 0L) stats::median(completed$p_value_raw) else NA_real_,
        prop_p_lt_0_05 = if (nrow(completed) > 0L) mean(completed$p_value_raw < 0.05) else NA_real_,
        prop_p_lt_0_01 = if (nrow(completed) > 0L) mean(completed$p_value_raw < 0.01) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }))
  utils::write.csv(summary_rows, file.path(output_dir, "sunspots_cycles21_23_weighted_mixture_rolling_10cr_summary_by_cycle.csv"), row.names = FALSE)
  elapsed_total <- proc.time()[["elapsed"]] - start_elapsed
  sunspots_weighted_rolling_log_message(log_path, sprintf("Run finished in %.3f seconds.", elapsed_total))

  list(windows = windows_df, alignment = alignment_out, results = results_df, summary_by_cycle = summary_rows, elapsed_total = elapsed_total)
}

if (sys.nframe() == 0L) {
  cli_args <- parse_sunspots_weighted_rolling_args()
  output <- do.call(run_sunspots_weighted_mixture_rolling_windows_gof, cli_args)
  if (!is.null(output$alignment$summary)) print(output$alignment$summary)
  if (!is.null(output$summary_by_cycle)) print(output$summary_by_cycle)
  cat(sprintf("Total elapsed seconds: %.3f\n", output$elapsed_total))
}
