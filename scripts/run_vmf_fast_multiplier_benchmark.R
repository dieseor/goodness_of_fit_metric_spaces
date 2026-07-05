resolve_vmf_fast_multiplier_benchmark_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

source(resolve_vmf_fast_multiplier_benchmark_path("bootstrap", "calibration_study.R"))
source(resolve_vmf_fast_multiplier_benchmark_path("scripts", "path_helpers.R"))

build_vmf_fast_multiplier_summary <- function(raw_df) {
  groups <- split(raw_df, list(raw_df$model, raw_df$scenario, raw_df$n, raw_df$statistic), drop = TRUE)
  rows <- lapply(groups, function(df) {
    data.frame(
      model = df$model[[1L]],
      scenario = df$scenario[[1L]],
      n = df$n[[1L]],
      statistic = df$statistic[[1L]],
      M_outer = nrow(df),
      old_p_value_mean = mean(df$old_p_value),
      fast_p_value_mean = mean(df$fast_p_value),
      p_value_correlation = stats::cor(df$old_p_value, df$fast_p_value),
      abs_diff_mean = mean(df$abs_p_value_diff),
      abs_diff_median = stats::median(df$abs_p_value_diff),
      old_rejection_rate_0_05 = mean(df$old_reject_0_05),
      fast_rejection_rate_0_05 = mean(df$fast_reject_0_05),
      old_runtime_mean = mean(df$old_total_seconds),
      fast_runtime_mean = mean(df$fast_total_seconds),
      median_speedup = stats::median(df$speedup_factor),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_vmf_fast_multiplier_scaling_summary <- function(raw_df) {
  groups <- split(raw_df, list(raw_df$B, raw_df$statistic), drop = TRUE)
  rows <- lapply(groups, function(df) {
    data.frame(
      B = df$B[[1L]],
      statistic = df$statistic[[1L]],
      median_old_prep_seconds = stats::median(df$old_prep_seconds),
      median_old_loop_seconds = stats::median(df$old_loop_seconds),
      median_old_total_seconds = stats::median(df$old_total_seconds),
      median_fast_prep_seconds = stats::median(df$fast_prep_seconds),
      median_fast_loop_seconds = stats::median(df$fast_loop_seconds),
      median_fast_total_seconds = stats::median(df$fast_total_seconds),
      median_speedup = stats::median(df$speedup_factor),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_one_vmf_fast_multiplier_case <- function(scenario,
                                             n,
                                             replicate_id,
                                             B,
                                             sample_seed,
                                             bootstrap_seed,
                                             derivative_mc_seed,
                                             statistics,
                                             derivative_mc_size) {
  set.seed(sample_seed)
  data <- simulate_h0_sample(
    scenario = scenario,
    n = n,
    replicate_id = replicate_id
  )

  old_result <- multiplier_bootstrap_vmf(
    data = data,
    null = scenario$null,
    statistics = statistics,
    ks_grid = scenario$ks_grid,
    B = B,
    alpha = 0.05,
    seed = bootstrap_seed,
    n_cores = 1,
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    distance_type = scenario$distance_type,
    unknown_param = scenario$unknown_param,
    bootstrap_method = "reestimated"
  )
  fast_result <- multiplier_bootstrap_vmf(
    data = data,
    null = scenario$null,
    statistics = statistics,
    ks_grid = scenario$ks_grid,
    B = B,
    alpha = 0.05,
    seed = bootstrap_seed,
    n_cores = 1,
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    distance_type = scenario$distance_type,
    unknown_param = scenario$unknown_param,
    bootstrap_method = "fast_multiplier",
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = derivative_mc_size,
      derivative_mc_seed = derivative_mc_seed,
      fast_multiplier_cvm_block_size = n
    )
  )

  do.call(rbind, lapply(statistics, function(stat_name) {
    old_inf <- old_result$inference[[stat_name]]
    fast_inf <- fast_result$inference[[stat_name]]
    data.frame(
      model = scenario$model,
      scenario = scenario$id,
      n = n,
      B = B,
      statistic = stat_name,
      replicate_id = replicate_id,
      observed_statistic = old_inf$observed,
      derivative_method = fast_result$diagnostics$derivative_method,
      derivative_mc_size = fast_result$diagnostics$derivative_mc_size,
      derivative_mc_seed = fast_result$diagnostics$derivative_mc_seed,
      common_observed_seconds = fast_result$diagnostics$common_observed_seconds,
      old_p_value = old_inf$p_value,
      fast_p_value = fast_inf$p_value,
      abs_p_value_diff = abs(old_inf$p_value - fast_inf$p_value),
      old_prep_seconds = old_result$diagnostics$old_prep_seconds,
      old_loop_seconds = old_result$diagnostics$old_loop_seconds,
      old_total_seconds = old_result$diagnostics$old_total_seconds,
      fast_prep_seconds = fast_result$diagnostics$fast_prep_seconds,
      fast_loop_seconds = fast_result$diagnostics$fast_loop_seconds,
      fast_total_seconds = fast_result$diagnostics$fast_total_seconds,
      speedup_factor = old_result$diagnostics$old_total_seconds / fast_result$diagnostics$fast_total_seconds,
      old_reject_0_01 = as.integer(old_inf$p_value <= 0.01),
      old_reject_0_05 = as.integer(old_inf$p_value <= 0.05),
      old_reject_0_10 = as.integer(old_inf$p_value <= 0.10),
      fast_reject_0_01 = as.integer(fast_inf$p_value <= 0.01),
      fast_reject_0_05 = as.integer(fast_inf$p_value <= 0.05),
      fast_reject_0_10 = as.integer(fast_inf$p_value <= 0.10),
      stringsAsFactors = FALSE
    )
  }))
}

run_vmf_fast_multiplier_benchmark <- function(mode = c("pilot", "scaling"),
                                              output_root = NULL,
                                             derivative_mc_size = 1000L,
                                              base_seed = 20260612L) {
  mode <- match.arg(mode)
  scenario <- make_vmf_composite_calibration_scenario(2.0)
  n_value <- 50L
  statistics <- c("ks", "cvm")
  if (is.null(output_root)) {
    output_root <- if (identical(mode, "pilot")) {
      canonical_calibration_bootstrap_dir("vmf", "fast", "pilot_kappa2_n50_B1000")
    } else {
      canonical_calibration_bootstrap_dir("vmf", "fast", "scaling_kappa2_n50_B100_300_1000")
    }
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  if (identical(mode, "pilot")) {
    B_values <- 1000L
    M_outer <- 20L
  } else {
    B_values <- c(100L, 300L, 1000L)
    M_outer <- 3L
  }

  rows <- list()
  idx <- 1L
  for (B_value in B_values) {
    for (replicate_id in seq_len(M_outer)) {
      sample_seed <- as.integer(base_seed + 10000L * B_value + replicate_id)
      bootstrap_seed <- as.integer(base_seed + 20000L * B_value + replicate_id)
      derivative_mc_seed <- as.integer(base_seed + 30000L * B_value + replicate_id)
      message(sprintf("[vMF fast multiplier %s] B=%d replicate=%d/%d", mode, B_value, replicate_id, M_outer))
      rows[[idx]] <- run_one_vmf_fast_multiplier_case(
        scenario = scenario,
        n = n_value,
        replicate_id = replicate_id,
        B = B_value,
        sample_seed = sample_seed,
        bootstrap_seed = bootstrap_seed,
        derivative_mc_seed = derivative_mc_seed,
        statistics = statistics,
        derivative_mc_size = derivative_mc_size
      )
      idx <- idx + 1L
    }
  }

  raw_df <- do.call(rbind, rows)
  summary_df <- if (identical(mode, "pilot")) {
    build_vmf_fast_multiplier_summary(raw_df)
  } else {
    build_vmf_fast_multiplier_scaling_summary(raw_df)
  }

  utils::write.csv(raw_df, file.path(output_root, "benchmark_raw.csv"), row.names = FALSE)
  utils::write.csv(summary_df, file.path(output_root, "benchmark_summary.csv"), row.names = FALSE)

  list(output_root = output_root, raw = raw_df, summary = summary_df)
}

parse_named_args_vmf_fast_multiplier <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_vmf_fast_multiplier(commandArgs(trailingOnly = TRUE))
  run_vmf_fast_multiplier_benchmark(
    mode = as.character(args$mode %||% "pilot"),
    output_root = args$output_root %||% NULL,
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L),
    base_seed = as.integer(args$base_seed %||% 20260612L)
  )
}
