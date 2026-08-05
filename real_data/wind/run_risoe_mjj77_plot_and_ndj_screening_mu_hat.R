#!/usr/bin/env Rscript

# Unified fast screening for two 77 m Risoe windows using the AoS wind-data rule:
# - Windows: May-June-July and November-December-January.
# - Years: 1996-2004.
# - Daily sampling: closest-to-noon timestamp on days 4, 8, 12, 16, 20, 24, 28.
# - Quality filter: finite speed/direction and speed > 0.
# - No imputation; duplicated daily selections are rejected.
#
# Run from any directory, for example:
#   B=5000 Rscript real_data/wind/run_risoe_mjj77_plot_and_ndj_screening_mu_hat.R

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
} else if (length(script_argument) == 0L) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("More than one --file argument was supplied.", call. = FALSE)
}
setwd(repo_root)

source(file.path(repo_root, "real_data", "wind", "run_risoe_125m_screening_ks_cvm.R"))

read_positive_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("%s must be a positive integer.", name), call. = FALSE)
  }
  value
}

mu_hat_string <- function(mu) {
  sprintf("(%.10f, %.10f, %.10f)", unname(mu[[1L]]), unname(mu[[2L]]), unname(mu[[3L]]))
}

scalar_diag <- function(diagnostics, name) {
  value <- diagnostics[[name]]
  if (is.null(value) || length(value) != 1L) NA else value
}

assert_requested_fast_implementation <- function(result) {
  d <- result$diagnostics
  violations <- c(
    effective_bootstrap_method = !identical(as.character(scalar_diag(d, "effective_bootstrap_method")), "fast_multiplier"),
    fallback_to_reestimated = isTRUE(scalar_diag(d, "fallback_to_reestimated")),
    derivative_method = !identical(as.character(scalar_diag(d, "derivative_method_effective")), "quadrature"),
    distance_profile_backend = !identical(as.character(scalar_diag(d, "distance_profile_backend_effective")), "r"),
    fast_multiplier_backend = !identical(as.character(scalar_diag(d, "fast_multiplier_backend_effective")), "cpp"),
    fast_multiplier_cpp_kernel = !identical(as.character(scalar_diag(d, "fast_multiplier_cpp_kernel_effective")), "contiguous_double"),
    fused_KS_CvM = !isTRUE(scalar_diag(d, "fast_multiplier_fuse_ks_cvm_effective"))
  )

  if (any(violations)) {
    stop(
      sprintf(
        "Requested fast implementation was not obtained: %s.",
        paste(names(violations)[violations], collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

make_configs <- function() {
  list(
    list(
      case_id = "risoe_77m_may_jun_jul_start4_1996_2004",
      window_id = "may_jun_jul",
      months = c(5L, 6L, 7L),
      pattern = "start4",
      day_pattern = c(4L, 8L, 12L, 16L, 20L, 24L, 28L),
      year_start = 1996L,
      year_end = 2004L,
      cutoff_exclusive = as.Date("2005-01-01"),
      speed_col = "ws77",
      direction_col = "wd77",
      height_m = 77L
    ),
    list(
      case_id = "risoe_77m_nov_dec_jan_start4_1996_2004",
      window_id = "nov_dec_jan",
      months = c(11L, 12L, 1L),
      pattern = "start4",
      day_pattern = c(4L, 8L, 12L, 16L, 20L, 24L, 28L),
      year_start = 1996L,
      year_end = 2004L,
      cutoff_exclusive = as.Date("2005-01-01"),
      speed_col = "ws77",
      direction_col = "wd77",
      height_m = 77L
    )
  )
}

build_case_data <- function(selected_df, config, fixed_tz = "UTC") {
  dates <- as.Date(selected_df$datetime, tz = fixed_tz)

  selected <- selected_df[
    selected_df$year >= config$year_start &
      selected_df$year <= config$year_end &
      selected_df$month %in% config$months &
      selected_df$day %in% config$day_pattern &
      dates < config$cutoff_exclusive,
    , drop = FALSE
  ]

  if (nrow(selected) == 0L) {
    stop("No rows remain after year/month/day filtering.", call. = FALSE)
  }

  duplicate_dates_selected <- sum(duplicated(as.Date(selected$datetime, tz = fixed_tz)))
  if (duplicate_dates_selected > 0L) {
    stop("Daily selection produced duplicate dates before value filtering.", call. = FALSE)
  }

  valid <- is.finite(selected[[config$speed_col]]) &
    is.finite(selected[[config$direction_col]]) &
    selected[[config$speed_col]] > 0

  rows <- selected[valid, , drop = FALSE]
  if (nrow(rows) < 2L) {
    stop("Fewer than two valid observations remain after quality filtering.", call. = FALSE)
  }

  duplicate_dates_valid <- sum(duplicated(as.Date(rows$datetime, tz = fixed_tz)))
  if (duplicate_dates_valid > 0L) {
    stop("Daily selection produced duplicate dates after value filtering.", call. = FALSE)
  }

  speed <- as.numeric(rows[[config$speed_col]])
  direction_deg <- ((as.numeric(rows[[config$direction_col]]) %% 360) + 360) %% 360
  speed_scaled <- speed / mean(speed)
  angle <- direction_deg * pi / 180

  X <- cbind(
    cosh(speed_scaled),
    sinh(speed_scaled) * cos(angle),
    sinh(speed_scaled) * sin(angle)
  )

  if (any(abs(-X[, 1]^2 + X[, 2]^2 + X[, 3]^2 + 1) > 1e-8)) {
    stop("The transformed observations are not on the unit hyperboloid.", call. = FALSE)
  }

  list(
    X = X,
    n_selected = nrow(selected),
    n_valid = nrow(rows),
    n_removed_missing_or_invalid = nrow(selected) - nrow(rows),
    n_selected_2004 = sum(selected$year == 2004L),
    n_valid_2004 = sum(rows$year == 2004L),
    years_valid = paste(sort(unique(rows$year)), collapse = ","),
    first_selected_datetime = as.character(min(selected$datetime)),
    last_selected_datetime = as.character(max(selected$datetime)),
    first_valid_datetime = as.character(min(rows$datetime)),
    last_valid_datetime = as.character(max(rows$datetime))
  )
}

run_case <- function(selected_df, config, B, n_cores, seed) {
  started <- Sys.time()
  case_data <- build_case_data(selected_df, config)
  fit <- hvmf_mle_h2(case_data$X)

  result <- multiplier_bootstrap_hvmf(
    data = case_data$X,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE),
    control = list(
      derivative_method = "quadrature",
      hvmf_profile_method = "tabulated",
      hvmf_profile_n_y = 4097L,
      fast_multiplier_cvm_block_size = 50L,
      fast_bootstrap_chunk_size = 100L,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto"
    ),
    unknown_param = "both",
    distance_profile_backend = "r",
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fuse_ks_cvm = TRUE,
    cache_block_corrections = "auto"
  )

  assert_requested_fast_implementation(result)
  d <- result$diagnostics

  common <- data.frame(
    case_id = config$case_id,
    window_id = config$window_id,
    height_m = config$height_m,
    years = sprintf("%d-%d", config$year_start, config$year_end),
    months = paste(config$months, collapse = ","),
    pattern = config$pattern,
    day_pattern = paste(config$day_pattern, collapse = ","),
    cutoff_exclusive = as.character(config$cutoff_exclusive),
    n_selected_before_quality_filter = case_data$n_selected,
    n_valid = case_data$n_valid,
    n_removed_missing_or_invalid = case_data$n_removed_missing_or_invalid,
    n_selected_2004_before_quality_filter = case_data$n_selected_2004,
    n_valid_2004 = case_data$n_valid_2004,
    years_valid = case_data$years_valid,
    first_selected_datetime = case_data$first_selected_datetime,
    last_selected_datetime = case_data$last_selected_datetime,
    first_valid_datetime = case_data$first_valid_datetime,
    last_valid_datetime = case_data$last_valid_datetime,
    B = B,
    n_cores = n_cores,
    seed = seed,
    kappa_hat = fit$kappa,
    theta_deg_hat = fit$theta_deg,
    mu1_hat = unname(fit$mu[[1L]]),
    mu2_hat = unname(fit$mu[[2L]]),
    mu3_hat = unname(fit$mu[[3L]]),
    mu_hat = mu_hat_string(fit$mu),
    ks_grid_mode = "sample_points_unique_distances",
    bootstrap_method_requested = "fast_multiplier",
    bootstrap_method_effective = as.character(scalar_diag(d, "effective_bootstrap_method")),
    derivative_method_requested = "quadrature",
    derivative_method_effective = as.character(scalar_diag(d, "derivative_method_effective")),
    distance_profile_backend_requested = "r",
    distance_profile_backend_effective = as.character(scalar_diag(d, "distance_profile_backend_effective")),
    fast_multiplier_backend_effective = as.character(scalar_diag(d, "fast_multiplier_backend_effective")),
    fast_multiplier_cpp_kernel_effective = as.character(scalar_diag(d, "fast_multiplier_cpp_kernel_effective")),
    fast_multiplier_fuse_ks_cvm_effective = isTRUE(scalar_diag(d, "fast_multiplier_fuse_ks_cvm_effective")),
    fallback_to_reestimated = isTRUE(scalar_diag(d, "fallback_to_reestimated")),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    status = "ok",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  rbind(
    cbind(common, data.frame(statistic = "KS", observed_statistic = result$observed$ks$statistic, p_value = result$inference$ks$p_value)),
    cbind(common, data.frame(statistic = "CvM", observed_statistic = result$observed$cvm$statistic, p_value = result$inference$cvm$p_value))
  )
}

make_error_case_result <- function(config, B, n_cores, seed, message) {
  data.frame(
    case_id = config$case_id,
    window_id = config$window_id,
    height_m = config$height_m,
    years = sprintf("%d-%d", config$year_start, config$year_end),
    months = paste(config$months, collapse = ","),
    pattern = config$pattern,
    day_pattern = paste(config$day_pattern, collapse = ","),
    cutoff_exclusive = as.character(config$cutoff_exclusive),
    n_selected_before_quality_filter = NA_integer_,
    n_valid = NA_integer_,
    n_removed_missing_or_invalid = NA_integer_,
    n_selected_2004_before_quality_filter = NA_integer_,
    n_valid_2004 = NA_integer_,
    years_valid = NA_character_,
    first_selected_datetime = NA_character_,
    last_selected_datetime = NA_character_,
    first_valid_datetime = NA_character_,
    last_valid_datetime = NA_character_,
    B = B,
    n_cores = n_cores,
    seed = seed,
    kappa_hat = NA_real_,
    theta_deg_hat = NA_real_,
    mu1_hat = NA_real_,
    mu2_hat = NA_real_,
    mu3_hat = NA_real_,
    mu_hat = NA_character_,
    ks_grid_mode = "sample_points_unique_distances",
    bootstrap_method_requested = "fast_multiplier",
    bootstrap_method_effective = NA_character_,
    derivative_method_requested = "quadrature",
    derivative_method_effective = NA_character_,
    distance_profile_backend_requested = "r",
    distance_profile_backend_effective = NA_character_,
    fast_multiplier_backend_effective = NA_character_,
    fast_multiplier_cpp_kernel_effective = NA_character_,
    fast_multiplier_fuse_ks_cvm_effective = NA,
    fallback_to_reestimated = NA,
    elapsed_seconds = NA_real_,
    status = "error",
    error_message = message,
    statistic = "ERROR",
    observed_statistic = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )
}

completed_case <- function(results, case_id, B) {
  if (is.null(results) || nrow(results) == 0L) return(FALSE)

  rows <- results[
    results$case_id == case_id &
      results$B == B &
      results$status == "ok",
    , drop = FALSE
  ]

  if (nrow(rows) != 2L) return(FALSE)
  identical(sort(as.character(rows$statistic)), c("CvM", "KS"))
}

make_summary <- function(results) {
  good <- results[results$status == "ok", , drop = FALSE]
  if (nrow(good) == 0L) return(data.frame())

  ids <- unique(good$case_id)
  rows <- lapply(ids, function(id) {
    z <- good[good$case_id == id, , drop = FALSE]
    ks <- z[z$statistic == "KS", , drop = FALSE]
    cvm <- z[z$statistic == "CvM", , drop = FALSE]
    if (nrow(ks) != 1L || nrow(cvm) != 1L) return(NULL)

    data.frame(
      case_id = id,
      window_id = ks$window_id,
      height_m = ks$height_m,
      years = ks$years,
      months = ks$months,
      pattern = ks$pattern,
      day_pattern = ks$day_pattern,
      n_valid = ks$n_valid,
      n_removed_missing_or_invalid = ks$n_removed_missing_or_invalid,
      first_valid_datetime = ks$first_valid_datetime,
      last_valid_datetime = ks$last_valid_datetime,
      kappa_hat = ks$kappa_hat,
      theta_deg_hat = ks$theta_deg_hat,
      mu1_hat = ks$mu1_hat,
      mu2_hat = ks$mu2_hat,
      mu3_hat = ks$mu3_hat,
      mu_hat = ks$mu_hat,
      p_value_KS = ks$p_value,
      p_value_CvM = cvm$p_value,
      reject_KS_unadjusted_5pct = ks$p_value < .05,
      reject_CvM_unadjusted_5pct = cvm$p_value < .05,
      reject_either_unadjusted_5pct = ks$p_value < .05 || cvm$p_value < .05,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, Filter(Negate(is.null), rows))
  out[order(out$window_id), , drop = FALSE]
}

run_screening_two_windows <- function(
    input_nc = file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"),
    B = read_positive_integer("B", 5000L),
    base_seed = read_positive_integer("BASE_SEED", 2026080400L),
    output_dir = Sys.getenv(
      "OUTPUT_DIR",
      unset = file.path(
        repo_root,
        "real_data",
        "wind",
        "month_diagnostics",
        sprintf("screening_fast_mjj_ndj_77m_start4_1996_2004_b%d_mu_hat", B)
      )
    )) {
  n_cores <- read_positive_integer("N_CORES", 6L)
  if (!file.exists(input_nc)) {
    stop(sprintf("Input NetCDF file not found: %s", input_nc), call. = FALSE)
  }

  configs <- make_configs()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  checkpoint_path <- file.path(output_dir, sprintf("risoe_fast_mjj_ndj_77m_start4_1996_2004_b%d_checkpoint.csv", B))
  summary_path <- file.path(output_dir, sprintf("risoe_fast_mjj_ndj_77m_start4_1996_2004_b%d_summary.csv", B))

  results <- if (file.exists(checkpoint_path)) {
    utils::read.csv(checkpoint_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }

  selected_df <- select_noon_all_months(load_risoe_concurrent(input_nc, fixed_tz = "UTC"), fixed_tz = "UTC")

  cat("AoS wind-data rule: years 1996-2004, noon-nearest on days 4,8,12,16,20,24,28, no imputation.\n")
  cat(sprintf("Running %d cases at 77 m with B=%d and cores=%d.\n\n", length(configs), B, n_cores))

  for (i in seq_along(configs)) {
    config <- configs[[i]]

    if (completed_case(results, config$case_id, B)) {
      cat(sprintf("[%d/%d] Resuming: %s already complete.\n", i, length(configs), config$case_id))
      next
    }

    seed <- as.integer(base_seed + i)
    cat(sprintf(
      "[%d/%d] %s | months=%s | pattern=%s | B=%d | cores=%d | seed=%d\n",
      i,
      length(configs),
      config$case_id,
      paste(config$months, collapse = ","),
      config$pattern,
      B,
      n_cores,
      seed
    ))
    flush.console()

    case_result <- tryCatch(
      run_case(selected_df, config, B, n_cores, seed),
      error = function(e) make_error_case_result(config, B, n_cores, seed, conditionMessage(e))
    )

    if (!is.null(results) && nrow(results) > 0L) {
      results <- results[!(results$case_id == config$case_id & results$B == B), , drop = FALSE]
      if (nrow(results) == 0L) results <- NULL
    }

    results <- if (is.null(results)) case_result else rbind(results, case_result)
    utils::write.csv(results, checkpoint_path, row.names = FALSE)
    utils::write.csv(make_summary(results), summary_path, row.names = FALSE)

    if (all(case_result$status == "ok")) {
      cat(sprintf(
        "          n=%d; removed=%d; first=%s; last=%s; p_KS=%.6f; p_CvM=%.6f; mu_hat=%s\n",
        case_result$n_valid[[1L]],
        case_result$n_removed_missing_or_invalid[[1L]],
        case_result$first_valid_datetime[[1L]],
        case_result$last_valid_datetime[[1L]],
        case_result$p_value[case_result$statistic == "KS"],
        case_result$p_value[case_result$statistic == "CvM"],
        case_result$mu_hat[[1L]]
      ))
    } else {
      stop(sprintf("Case %s failed and was checkpointed: %s", config$case_id, case_result$error_message[[1L]]), call. = FALSE)
    }

    flush.console()
  }

  cat(sprintf(
    "\nCompleted. Checkpoint: %s\nSummary: %s\n",
    normalizePath(checkpoint_path, winslash = "/", mustWork = TRUE),
    normalizePath(summary_path, winslash = "/", mustWork = TRUE)
  ))

  summary_df <- make_summary(results)
  if (nrow(summary_df) > 0L) {
    print(summary_df, row.names = FALSE, digits = 6)
  }

  invisible(list(checkpoint_path = checkpoint_path, summary_path = summary_path, summary = summary_df))
}

if (sys.nframe() == 0L) run_screening_two_windows()
