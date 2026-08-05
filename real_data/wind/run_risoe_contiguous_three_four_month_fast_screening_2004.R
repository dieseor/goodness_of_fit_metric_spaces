#!/usr/bin/env Rscript

# Sequential screening of the Risoe wind series with the fast composite-HvMF
# multiplier bootstrap.  The outer loop is serial; each individual test uses
# exactly ten cores.  The multiplier kernel is C++; the HvMF distance-profile
# backend remains R because its C++ implementation is correctly gated out by
# the project when its exactness/performance checks do not pass.
#
# Data rule: use every finite speed/direction pair in 1996--2004 after the
# closest-to-noon daily selection.  The 77 m series is used through 2004-12-31.
# The 125 m series is used only before 2004-11-12, the first day at which its
# wind direction becomes frozen; this prevents the defective tail entering any
# month window.  Missing observations are discarded, never imputed.
#
# Run from any directory, for example:
#   B=5000 Rscript real_data/wind/run_risoe_contiguous_three_four_month_fast_screening_2004.R

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
} else if (length(script_argument) == 0L) {
  # This branch is for programmatic inspection with `source()` from the root.
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

month_labels <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
day_patterns <- setNames(lapply(1:4, function(start) seq.int(start, 30L, by = 4L)), paste0("start", 1:4))
make_configs <- function() {
  configs <- list()
  index <- 1L

  for (height_spec in list(
    # list(height_m = 77L, speed_col = "ws77", direction_col = "wd77", cutoff = as.Date("2005-01-01")),
    list(
      height_m = 125L,
      speed_col = "ws125",
      direction_col = "wd125",
      cutoff = as.Date("2004-11-12")
    )
  )) {
    for (window_length in c(3L)) {
      for (first_month in 5L) {
        months <- seq.int(first_month, length.out = window_length)

        for (pattern in names(day_patterns)) {
          configs[[index]] <- c(
            height_spec,
            list(
              months = months,
              months_id = paste(month_labels[months], collapse = "_"),
              window_length = window_length,
              first_month = first_month,
              last_month = max(months),
              pattern = pattern,
              day_pattern = day_patterns[[pattern]],
              case_id = sprintf(
                "risoe_%dm_%s_%s_1996_2004",
                height_spec$height_m,
                paste(month_labels[months], collapse = "_"),
                pattern
              )
            )
          )

          index <- index + 1L
        }
      }
    }
  }

  configs
}

build_case_data <- function(selected_df, config, fixed_tz = "UTC") {
  dates <- as.Date(selected_df$datetime, tz = fixed_tz)
  selected <- selected_df[
    selected_df$year >= 1996L & selected_df$year <= 2004L &
      selected_df$month %in% config$months & selected_df$day %in% config$day_pattern &
      dates < config$cutoff,
    , drop = FALSE
  ]
  n_selected <- nrow(selected)
  n_selected_2004 <- sum(selected$year == 2004L)
  valid <- is.finite(selected[[config$speed_col]]) &
    is.finite(selected[[config$direction_col]]) & selected[[config$speed_col]] > 0
  rows <- selected[valid, , drop = FALSE]
  if (nrow(rows) < 2L) {
    stop("Fewer than two valid observations remain after the prespecified data rule.", call. = FALSE)
  }
  if (anyDuplicated(as.Date(rows$datetime, tz = fixed_tz))) {
    stop("Daily selection produced duplicate dates.", call. = FALSE)
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
    n_selected = n_selected,
    n_valid = nrow(rows),
    n_selected_2004 = n_selected_2004,
    n_valid_2004 = sum(rows$year == 2004L),
    years_valid = paste(sort(unique(rows$year)), collapse = ",")
  )
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
    stop(sprintf("Requested fast implementation was not obtained: %s.", paste(names(violations)[violations], collapse = ", ")), call. = FALSE)
  }
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
    case_id = config$case_id, height_m = config$height_m,
    window_length = config$window_length, month_start = config$first_month,
    month_end = config$last_month, months = paste(config$months, collapse = ","),
    pattern = config$pattern, day_pattern = paste(config$day_pattern, collapse = ","),
    data_cutoff_exclusive = as.character(config$cutoff),
    n_selected_before_value_filter = case_data$n_selected,
    n_valid = case_data$n_valid,
    n_selected_2004_before_value_filter = case_data$n_selected_2004,
    n_valid_2004 = case_data$n_valid_2004,
    years_valid = case_data$years_valid,
    B = B, n_cores = n_cores, seed = seed,
    kappa_hat = fit$kappa, theta_deg_hat = fit$theta_deg,
    mu1_hat = unname(fit$mu[[1L]]),
    mu2_hat = unname(fit$mu[[2L]]),
    mu3_hat = unname(fit$mu[[3L]]),
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
    status = "ok", error_message = NA_character_, stringsAsFactors = FALSE
  )
  rbind(
    cbind(common, data.frame(statistic = "KS", observed_statistic = result$observed$ks$statistic, p_value = result$inference$ks$p_value)),
    cbind(common, data.frame(statistic = "CvM", observed_statistic = result$observed$cvm$statistic, p_value = result$inference$cvm$p_value))
  )
}

completed_case <- function(results, case_id, B) {
  if (is.null(results) || nrow(results) == 0L) {
    return(FALSE)
  }

  rows <- results[
    results$case_id == case_id &
      results$B == B &
      results$status == "ok",
    ,
    drop = FALSE
  ]

  if (nrow(rows) != 2L) {
    return(FALSE)
  }

  identical(
    sort(as.character(rows$statistic)),
    c("CvM", "KS")
  )
}

make_summary <- function(results) {
  good <- results[results$status == "ok", , drop = FALSE]
  if (nrow(good) == 0L) return(data.frame())
  ids <- unique(good$case_id)
  out <- lapply(ids, function(id) {
    z <- good[good$case_id == id, , drop = FALSE]
    ks <- z[z$statistic == "KS", , drop = FALSE]
    cvm <- z[z$statistic == "CvM", , drop = FALSE]
    if (nrow(ks) != 1L || nrow(cvm) != 1L) return(NULL)
    data.frame(case_id = id, height_m = ks$height_m, window_length = ks$window_length,
      months = ks$months, pattern = ks$pattern, n_valid = ks$n_valid, n_valid_2004 = ks$n_valid_2004,
      p_value_KS = ks$p_value, p_value_CvM = cvm$p_value,
      reject_KS_unadjusted_5pct = ks$p_value < .05, reject_CvM_unadjusted_5pct = cvm$p_value < .05,
      reject_either_unadjusted_5pct = ks$p_value < .05 || cvm$p_value < .05, mu1_hat = ks$mu1_hat,
  mu2_hat = ks$mu2_hat,
  mu3_hat = ks$mu3_hat, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  out[order(out$height_m, out$window_length, out$months, out$pattern), , drop = FALSE]
}

make_error_case_result <- function(config, B, n_cores, seed, message) {
  data.frame(
    case_id = config$case_id, height_m = config$height_m,
    window_length = config$window_length, month_start = config$first_month,
    month_end = config$last_month, months = paste(config$months, collapse = ","),
    pattern = config$pattern, day_pattern = paste(config$day_pattern, collapse = ","),
    data_cutoff_exclusive = as.character(config$cutoff),
    n_selected_before_value_filter = NA_integer_, n_valid = NA_integer_,
    n_selected_2004_before_value_filter = NA_integer_, n_valid_2004 = NA_integer_,
    years_valid = NA_character_, B = B, n_cores = n_cores, seed = seed,
    kappa_hat = NA_real_, theta_deg_hat = NA_real_,
    mu1_hat = NA_real_,
mu2_hat = NA_real_,
mu3_hat = NA_real_,
    ks_grid_mode = "sample_points_unique_distances",
    bootstrap_method_requested = "fast_multiplier", bootstrap_method_effective = NA_character_,
    derivative_method_requested = "quadrature", derivative_method_effective = NA_character_,
    distance_profile_backend_requested = "r", distance_profile_backend_effective = NA_character_, fast_multiplier_backend_effective = NA_character_,
    fast_multiplier_cpp_kernel_effective = NA_character_, fast_multiplier_fuse_ks_cvm_effective = NA,
    fallback_to_reestimated = NA, elapsed_seconds = NA_real_, status = "error",
    error_message = message, statistic = "ERROR", observed_statistic = NA_real_, p_value = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_screening <- function(
    input_nc = file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"),
    B = read_positive_integer("B", 5000L),
    base_seed = read_positive_integer("BASE_SEED", 2026080300L),
    output_dir = Sys.getenv("OUTPUT_DIR", unset = file.path(repo_root, "real_data", "wind", "month_diagnostics", sprintf("screening_fast_contiguous_three_four_month_1996_2004_b%d", B)))) {
  n_cores <- 6L
  configs <- make_configs()
  if (!file.exists(input_nc)) stop(sprintf("Input NetCDF file not found: %s", input_nc), call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_path <- file.path(output_dir, sprintf("risoe_fast_contiguous_three_four_month_1996_2004_b%d_checkpoint.csv", B))
  summary_path <- file.path(output_dir, sprintf("risoe_fast_contiguous_three_four_month_1996_2004_b%d_summary.csv", B))
  results <- if (file.exists(checkpoint_path)) utils::read.csv(checkpoint_path, stringsAsFactors = FALSE, check.names = FALSE) else NULL
  selected_df <- select_noon_all_months(load_risoe_concurrent(input_nc, fixed_tz = "UTC"), fixed_tz = "UTC")
    # Upgrade old checkpoints by recomputing only the HvMF MLE.
  mu_columns <- c("mu1_hat", "mu2_hat", "mu3_hat")

  if (!is.null(results) && nrow(results) > 0L) {

    # Old checkpoints do not contain these columns.
    for (column in setdiff(mu_columns, names(results))) {
      results[[column]] <- NA_real_
    }

    for (config in configs) {
      rows <- results$case_id == config$case_id &
        results$B == B &
        results$status == "ok"

      if (!any(rows)) {
        next
      }

      mu_already_available <-
        all(is.finite(results$mu1_hat[rows])) &&
        all(is.finite(results$mu2_hat[rows])) &&
        all(is.finite(results$mu3_hat[rows]))

      if (mu_already_available) {
        next
      }

      case_data <- build_case_data(selected_df, config)
      fit <- hvmf_mle_h2(case_data$X)

      results$kappa_hat[rows] <- unname(fit$kappa)
      results$theta_deg_hat[rows] <- unname(fit$theta_deg)
      results$mu1_hat[rows] <- unname(fit$mu[[1L]])
      results$mu2_hat[rows] <- unname(fit$mu[[2L]])
      results$mu3_hat[rows] <- unname(fit$mu[[3L]])

      cat(sprintf(
        "Updated MLE: %s; mu_hat=(%.8f, %.8f, %.8f).\n",
        config$case_id,
        fit$mu[[1L]],
        fit$mu[[2L]],
        fit$mu[[3L]]
      ))
    }

    utils::write.csv(
      results,
      checkpoint_path,
      row.names = FALSE
    )

    utils::write.csv(
      make_summary(results),
      summary_path,
      row.names = FALSE
    )
  }
  cat(sprintf("Fast sequential screening: %d datasets; B=%d; ten cores in each bootstrap.\\n", length(configs), B))
  cat("77 m: all valid records through 2004-12-31. 125 m: all valid records before 2004-11-12.\\n\\n")
  for (i in seq_along(configs)) {
    config <- configs[[i]]
    if (completed_case(results, config$case_id, B)) {
      cat(sprintf("[%d/%d] Resuming: %s already complete.\\n", i, length(configs), config$case_id))
      next
    }
    seed <- as.integer(base_seed + i)
    cat(sprintf("[%d/%d] %dm | months=%s | %s | B=%d | cores=%d | seed=%d\\n", i, length(configs), config$height_m, paste(config$months, collapse = ","), config$pattern, B, n_cores, seed))
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
      cat(sprintf("          n=%d (2004: %d); p_KS=%.6f; p_CvM=%.6f\\n", case_result$n_valid[[1L]], case_result$n_valid_2004[[1L]], case_result$p_value[case_result$statistic == "KS"], case_result$p_value[case_result$statistic == "CvM"]))
    } else {
      stop(sprintf("Case %s failed and was checkpointed: %s", config$case_id, case_result$error_message[[1L]]), call. = FALSE)
    }
    flush.console()
  }
  cat(sprintf("Completed. Checkpoint: %s\\nSummary: %s\\n", normalizePath(checkpoint_path, winslash = "/", mustWork = TRUE), normalizePath(summary_path, winslash = "/", mustWork = TRUE)))
}

if (sys.nframe() == 0L) run_screening()
