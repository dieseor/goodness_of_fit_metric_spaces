resolve_c2_short_joint_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_c2_short_joint_path("scripts", "run_comets_distance_profile_cardioid.R"))

parse_c2_short_joint_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- if (length(parts) > 1L) {
      paste(parts[-1L], collapse = "=")
    } else {
      TRUE
    }
  }
  out
}

run_comets_cardioid_c2_short_joint_fast <- function(
    output_dir = "real_data/comets/cardioid/recheck_C2_short_B5000_4cores_joint_fast",
    B = 5000L,
    n_cores = 4L,
    seed = 20260713L) {
  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  seed <- as.integer(seed)
  if (!is.finite(B) || B <= 0L) stop("`B` must be a positive integer.")
  if (!is.finite(n_cores) || n_cores <= 0L) stop("`n_cores` must be a positive integer.")
  if (!is.finite(seed)) stop("`seed` must be a finite integer.")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  comets_data <- load_comets_distance_profile_data()
  data_matrix <- comets_data$short$normal
  model <- make_comet_cardioid_models("C2")[[1L]]
  spec <- make_cardioid_spec(k = model$k, distance_type = "geodesic", unknown_param = "both")

  start_time <- Sys.time()
  result <- multiplier_bootstrap_gof(
    data = data_matrix,
    spec = spec,
    null = model$null,
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      cardioid_optim_control = list(maxit = 1000L),
      derivative_method = "score_mc",
      derivative_mc_size = 1000L,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto"
    )
  )
  finished_at <- Sys.time()
  result$diagnostics$requested_bootstrap_method <- "fast_multiplier"

  required_checks <- c(
    effective_fast = identical(result$diagnostics$effective_bootstrap_method, "fast_multiplier"),
    no_reestimated_fallback = !isTRUE(result$diagnostics$fallback_to_reestimated),
    lightweight_ks = isTRUE(result$diagnostics$lightweight_ks_prep),
    lightweight_cvm = isTRUE(result$diagnostics$lightweight_cvm_prep),
    shared_sample_prep = isTRUE(result$diagnostics$shared_sample_ks_cvm_cache),
    cpp_backend = identical(result$diagnostics$fast_multiplier_backend_effective, "cpp"),
    fused_ks_cvm = isTRUE(result$diagnostics$fast_multiplier_fuse_ks_cvm_effective),
    streamed_ks = identical(
      result$diagnostics$fast_ks_mode,
      "sample_points_unique_distances_streamed"
    ),
    streamed_cvm = identical(
      result$diagnostics$fast_cvm_mode,
      "sample_points_unique_distances_sorted_rows"
    )
  )
  if (!all(required_checks)) {
    stop(sprintf(
      "The joint fast KS/CvM contract failed: %s",
      paste(names(required_checks)[!required_checks], collapse = ", ")
    ))
  }

  summary <- do.call(rbind, lapply(c("ks", "cvm"), function(statistic) {
    summarize_cardioid_model_result(
      result = result,
      model = model,
      dataset_label = "Short-period cardioid KS+CvM",
      statistic = statistic
    )
  }))
  summary$joint_elapsed_seconds <- as.numeric(difftime(finished_at, start_time, units = "secs"))

  saveRDS(result, file.path(output_dir, "c2_short_joint_ks_cvm_result.rds"))
  utils::write.csv(summary, file.path(output_dir, "c2_short_joint_ks_cvm_summary.csv"), row.names = FALSE)
  saveRDS(
    list(
      output_dir = output_dir,
      B = B,
      n_cores = n_cores,
      seed = seed,
      statistics = c("ks", "cvm"),
      ks_grid_mode = "sample_points_unique_distances",
      bootstrap_method = "fast_multiplier",
      keep_observed_process = FALSE,
      required_checks = required_checks,
      started_at = start_time,
      finished_at = finished_at
    ),
    file.path(output_dir, "run_config.rds")
  )
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(list(result = result, summary = summary, checks = required_checks))
}

if (sys.nframe() == 0L) {
  args <- parse_c2_short_joint_args(commandArgs(trailingOnly = TRUE))
  run_comets_cardioid_c2_short_joint_fast(
    output_dir = args$output_dir %||%
      "real_data/comets/cardioid/recheck_C2_short_B5000_4cores_joint_fast",
    B = as.integer(args$B %||% 5000L),
    n_cores = as.integer(args$n_cores %||% 4L),
    seed = as.integer(args$seed %||% 20260713L)
  )
}
