#!/usr/bin/env Rscript

# Joint time--space GOF analysis for cycle-23 sunspot-group first records.
# This is intentionally separate from the legacy rank-based conditional runner.

resolve_sunspots_joint_runner_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

sunspots_joint_module_path <- resolve_sunspots_joint_runner_path(
  "real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"
)
sunspots_joint_model_spec_path <- resolve_sunspots_joint_runner_path(
  "bootstrap", "sunspots_joint_time_space_model_spec.R"
)
multiplier_bootstrap_engine_path <- resolve_sunspots_joint_runner_path(
  "bootstrap", "multiplier_bootstrap.R"
)
source(sunspots_joint_module_path)
source(sunspots_joint_model_spec_path)
source(multiplier_bootstrap_engine_path)

sunspots_joint_parse_statistics <- function(statistics) {
  normalize_requested_statistics(statistics)
}

sunspots_joint_summary_row <- function(statistic_name, bootstrap_values,
                                       observed_value, fit, fast_prep, data,
                                       center_indices, settings, timing) {
  theta <- fit$theta_hat
  eta <- fit$eta_hat
  data.frame(
    model = sprintf("joint_time_space_%s_small_circle_beta_mixture", settings$hemisphere_regression),
    hemisphere_regression = settings$hemisphere_regression,
    n_parameters = fit$n_parameters,
    statistic_type = if (identical(statistic_name, "ks")) {
      "ks_sample_centers_product_metric"
    } else {
      "cvm_sample_centers_product_metric"
    },
    bootstrap_method = "fast_multiplier_joint_time_space",
    profile_scope = "joint_time_space_product_law",
    product_metric = "0.5 * (geodesic/pi + absolute_time_difference)",
    n = nrow(data$x), n_sample_centers = length(center_indices), B = settings$B,
    derivative_mc_size = fast_prep$derivative_mc_size,
    n_cores = settings$n_cores,
    observed_profile_n_cores = settings$observed_profile_n_cores,
    start_date = settings$start_date, end_date_exclusive = settings$end_date,
    center_seed = settings$center_seed,
    dequantization_seed = settings$dequantization_seed,
    derivative_mc_seed = settings$derivative_mc_seed,
    bootstrap_seed = settings$bootstrap_seed,
    time_quad_n_per_component = settings$time_quad_n,
    profile_l_max = settings$profile_l_max,
    spatial_quad_n = settings$spatial_quad_n,
    distance_profile_backend = settings$distance_profile_backend,
    boundary_fast_override = settings$allow_boundary_fast,
    statistic = observed_value,
    critical_value_0.95 = as.numeric(stats::quantile(bootstrap_values, probs = 0.95, names = FALSE, type = 8)),
    p_value = (1 + sum(bootstrap_values >= observed_value)) / (length(bootstrap_values) + 1),
    temporal_weight1 = eta$weight1,
    temporal_alpha1 = eta$alpha1,
    temporal_beta1 = eta$beta1,
    temporal_alpha2 = eta$alpha2,
    temporal_beta2 = eta$beta2,
    temporal_mean1 = eta$mean1,
    temporal_mean2 = eta$mean2,
    temporal_loglik = fit$temporal_loglik,
    conditional_loglik = fit$conditional_loglik,
    joint_loglik = fit$loglik,
    aic = fit$aic,
    bic = fit$bic,
    a_N = theta$a_N, b_N = theta$b_N, a_S = theta$a_S, b_S = theta$b_S, c = theta$c,
    temporal_convergence = eta$opt$convergence,
    spatial_convergence = theta$opt$convergence,
    temporal_boundary_weight = eta$boundary_flags$weight,
    temporal_boundary_shape_lower = eta$boundary_flags$shape_lower,
    temporal_boundary_shape_upper = eta$boundary_flags$shape_upper,
    vhat_rcond = fast_prep$diagnostics$Vhat_rcond,
    vhat_condition_number = fast_prep$diagnostics$Vhat_condition_number %||% NA_real_,
    vhat_min_positive_eigenvalue = fast_prep$Vhat_min_positive_eigenvalue %||% NA_real_,
    vhat_max_positive_eigenvalue = fast_prep$Vhat_max_positive_eigenvalue %||% NA_real_,
    vhat_inversion_method = fast_prep$inversion_method %||% "",
    vhat_regularization_added = fast_prep$regularization_added %||% NA_real_,
    correction_all_finite = fast_prep$correction_all_finite %||% NA,
    temporal_mle_seconds = timing$temporal_mle_seconds,
    spatial_mle_seconds = timing$spatial_mle_seconds,
    observed_profile_seconds = timing$observed_profile_seconds,
    fast_preparation_seconds = timing$fast_preparation_seconds,
    bootstrap_seconds = timing$bootstrap_seconds,
    total_seconds = timing$total_seconds,
    stringsAsFactors = FALSE
  )
}

run_sunspots_cycle23_joint_time_space_gof <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle23_joint_time_space_gof_sampleks_fast"),
    start_date = "1997-06-01",
    end_date = "2006-01-01",
    statistics = "ks",
    hemisphere_regression = "asymmetric",
    B = 1000L,
    n_cores = 12L,
    observed_profile_n_cores = n_cores,
    center_seed = 20260711L,
    dequantization_seed = 20260712L,
    derivative_mc_seed = 20260713L,
    bootstrap_seed = 20260714L,
    n_sample_centers = 100L,
    derivative_mc_size = 5000L,
    bootstrap_block_size = 25L,
    time_quad_n = 64L,
    allow_boundary_fast = FALSE,
    control = list()) {
  statistics <- sunspots_joint_parse_statistics(statistics)
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  observed_profile_n_cores <- as.integer(observed_profile_n_cores)
  if (!is.finite(B) || B < 1L) stop("`B` must be a positive integer.")
  if (!is.finite(n_cores) || n_cores < 1L) stop("`n_cores` must be a positive integer.")
  if (!is.finite(observed_profile_n_cores) || observed_profile_n_cores < 1L) {
    stop("`observed_profile_n_cores` must be a positive integer.")
  }
  if (!is.finite(derivative_mc_size) || derivative_mc_size < 10L) {
    stop("`derivative_mc_size` must be an integer of at least 10.")
  }
  allow_boundary_fast <- isTRUE(allow_boundary_fast) || isTRUE(control$allow_boundary_fast)
  control <- utils::modifyList(control, list(
    hemisphere_regression = hemisphere_regression,
    allow_boundary_fast = allow_boundary_fast
  ))
  profile_l_max <- as.integer(control$profile_l_max %||% 100L)
  spatial_quad_n <- as.integer(control$profile_quad_n %||% 400L)
  center_block_size <- as.integer(control$center_block_size %||% 8L)
  distance_profile_backend <- sunspots_joint_normalize_backend(
    control$distance_profile_backend %||% "cpp"
  )
  effective_backend <- sunspots_joint_effective_backend(distance_profile_backend)
  settings <- list(
    hemisphere_regression = hemisphere_regression,
    B = B,
    n_cores = n_cores,
    observed_profile_n_cores = observed_profile_n_cores,
    start_date = start_date, end_date = end_date,
    center_seed = as.integer(center_seed),
    dequantization_seed = as.integer(dequantization_seed),
    derivative_mc_seed = as.integer(derivative_mc_seed),
    bootstrap_seed = as.integer(bootstrap_seed),
    time_quad_n = as.integer(time_quad_n), profile_l_max = profile_l_max,
    spatial_quad_n = spatial_quad_n, distance_profile_backend = effective_backend,
    allow_boundary_fast = allow_boundary_fast
  )

  total_start <- proc.time()[["elapsed"]]
  retained <- prepare_sunspots_cycle23_joint_time_space_data(
    input_csv = input_csv, start_date = start_date, end_date = end_date,
    dequantization_seed = dequantization_seed
  )
  data <- sunspots_joint_validate_data(
    as.matrix(retained[, c("x1", "x2", "x3")]), retained$s
  )

  temporal_start <- proc.time()[["elapsed"]]
  eta_hat <- fit_sunspots_joint_time_beta_mixture2(data$s, control)
  temporal_mle_seconds <- proc.time()[["elapsed"]] - temporal_start
  spatial_start <- proc.time()[["elapsed"]]
  theta_hat <- fit_sunspots_time_varying_asymmetric_mixture(
    data$x, data$s, hemisphere_regression = hemisphere_regression, control = control
  )
  spatial_mle_seconds <- proc.time()[["elapsed"]] - spatial_start
  fit <- fit_sunspots_cycle23_joint_time_space(
    data$x, data$s, hemisphere_regression = hemisphere_regression, control = control,
    eta_hat = eta_hat, theta_hat = theta_hat
  )

  center_indices <- seq_len(nrow(data$x))
  observed_profile_seconds <- NA_real_
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- plot_sunspots_cycle23_joint_time_space_diagnostics(data, fit, output_dir)
  retained$temporal_fitted_density <- sunspots_joint_time_density(data$s, eta_hat, control)
  retained$temporal_pit <- diagnostics$time_pit
  retained$conditional_axial_pit <- diagnostics$axial_pit
  utils::write.csv(retained, file.path(output_dir, "cycle23_joint_time_space_retained_data.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    hemisphere_regression = hemisphere_regression, n_parameters = fit$n_parameters,
    temporal_weight1 = eta_hat$weight1, temporal_alpha1 = eta_hat$alpha1,
    temporal_beta1 = eta_hat$beta1, temporal_alpha2 = eta_hat$alpha2,
    temporal_beta2 = eta_hat$beta2, temporal_mean1 = eta_hat$mean1,
    temporal_mean2 = eta_hat$mean2, temporal_loglik = fit$temporal_loglik,
    conditional_loglik = fit$conditional_loglik, joint_loglik = fit$loglik,
    aic = fit$aic, bic = fit$bic, a_N = theta_hat$a_N, b_N = theta_hat$b_N,
    a_S = theta_hat$a_S, b_S = theta_hat$b_S, c = theta_hat$c,
    temporal_convergence = eta_hat$opt$convergence, spatial_convergence = theta_hat$opt$convergence,
    temporal_boundary_weight = eta_hat$boundary_flags$weight,
    temporal_boundary_shape_lower = eta_hat$boundary_flags$shape_lower,
    temporal_boundary_shape_upper = eta_hat$boundary_flags$shape_upper,
    stringsAsFactors = FALSE
  ), file.path(output_dir, "cycle23_joint_time_space_mle.csv"), row.names = FALSE)
  utils::write.csv(
    eta_hat$boundary_diagnostics,
    file.path(output_dir, "cycle23_joint_time_space_temporal_boundary_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(data.frame(
    center_rank = seq_along(center_indices), center_index = center_indices,
    recorded_timestamp = retained$recorded_timestamp[center_indices],
    calendar_day = retained$calendar_day[center_indices],
    dequantization_jitter_day = retained$dequantization_jitter_day[center_indices],
    s = data$s[center_indices], x1 = data$x[center_indices, 1L],
    x2 = data$x[center_indices, 2L], x3 = data$x[center_indices, 3L]
  ), file.path(output_dir, "cycle23_joint_time_space_centers.csv"), row.names = FALSE)
  if (isTRUE(eta_hat$boundary_flags$weight) ||
      isTRUE(eta_hat$boundary_flags$shape_lower) ||
      isTRUE(eta_hat$boundary_flags$shape_upper)) {
    stop(
      paste(
        "Temporal beta-mixture MLE is on an admissible boundary;",
        "fast multiplier GOF for the joint model is not valid in this regime."
      ),
      call. = FALSE
    )
  }

  fast_preparation_seconds <- NA_real_
  bootstrap_start <- proc.time()[["elapsed"]]
  control <- utils::modifyList(control, list(
    hemisphere_regression = hemisphere_regression,
    derivative_mc_size = as.integer(derivative_mc_size),
    derivative_mc_seed = as.integer(derivative_mc_seed),
    fast_multiplier_vhat_rcond_tol = as.numeric(control$fast_multiplier_vhat_rcond_tol %||% 1e-12),
    distance_profile_backend = effective_backend,
    center_block_size = center_block_size,
    observed_profile_n_cores = observed_profile_n_cores,
    time_quad_n = as.integer(time_quad_n),
    profile_l_max = profile_l_max,
    profile_quad_n = spatial_quad_n,
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fast_multiplier_fuse_ks_cvm = TRUE,
    fast_multiplier_cache_corrections = "auto"
  ))
  joint_data_matrix <- cbind(data$x, data$s)
  joint_spec <- make_sunspots_joint_time_space_spec(
    hemisphere_regression = hemisphere_regression
  )
  generic_result <- multiplier_bootstrap_gof(
    data = joint_data_matrix,
    spec = joint_spec,
    null = list(type = "composite"),
    statistics = statistics,
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = n_cores,
    seed = as.integer(bootstrap_seed),
    observed_theta_hat = fit,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = control,
    distance_profile_backend = effective_backend
  )
  if (!identical(generic_result$diagnostics$effective_bootstrap_method, "fast_multiplier")) {
    stop(
      sprintf(
        paste(
          "Fast multiplier preparation was not valid for the joint model.",
          "effective_bootstrap_method = %s, fallback_to_reestimated = %s"
        ),
        as.character(generic_result$diagnostics$effective_bootstrap_method %||% NA_character_),
        as.character(generic_result$diagnostics$fallback_to_reestimated %||% NA)
      ),
      call. = FALSE
    )
  }
  bootstrap_statistics <- generic_result$bootstrap$statistics
  bootstrap_seconds <- proc.time()[["elapsed"]] - bootstrap_start
  fast_preparation_seconds <- as.numeric(generic_result$diagnostics$fast_prep_seconds %||% NA_real_)
  observed_profile_seconds <- as.numeric(generic_result$diagnostics$common_observed_seconds %||% NA_real_)
  numerical_diagnostics <- data.frame(
    n = nrow(data$x),
    n_sample_centers = length(center_indices),
    observed_profile_n_cores = observed_profile_n_cores,
    derivative_mc_size = generic_result$diagnostics$derivative_mc_size,
    Vhat_rcond = generic_result$diagnostics$Vhat_rcond,
    Vhat_condition_number = generic_result$diagnostics$Vhat_condition_number,
    Vhat_min_positive_eigenvalue = {
      eig <- generic_result$diagnostics$Vhat_eigenvalues
      eig <- eig[is.finite(eig) & eig > 0]
      if (length(eig) == 0L) NA_real_ else min(eig)
    },
    Vhat_max_positive_eigenvalue = {
      eig <- generic_result$diagnostics$Vhat_eigenvalues
      eig <- eig[is.finite(eig) & eig > 0]
      if (length(eig) == 0L) NA_real_ else max(eig)
    },
    vhat_inversion_method = "solve",
    vhat_regularization_added = 0,
    correction_all_finite = TRUE,
    correction_any_nonfinite = FALSE,
    correction_nonfinite_centers = 0,
    temporal_near_any_bound = any(eta_hat$boundary_diagnostics$near_any_bound),
    temporal_min_distance_to_lower_bound = min(eta_hat$boundary_diagnostics$distance_to_lower),
    temporal_min_distance_to_upper_bound = min(eta_hat$boundary_diagnostics$distance_to_upper),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    numerical_diagnostics,
    file.path(output_dir, "cycle23_joint_time_space_numerical_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      center_rank = integer(0),
      center_index = integer(0),
      correction_all_finite = logical(0),
      correction_max_abs = numeric(0),
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "cycle23_joint_time_space_center_correction_diagnostics.csv"),
    row.names = FALSE
  )
  total_seconds <- proc.time()[["elapsed"]] - total_start
  observed_statistics <- c(
    ks = as.numeric(generic_result$observed$ks$statistic %||% NA_real_),
    cvm = as.numeric(generic_result$observed$cvm$statistic %||% NA_real_)
  )[statistics]
  timing <- list(
    temporal_mle_seconds = temporal_mle_seconds, spatial_mle_seconds = spatial_mle_seconds,
    observed_profile_seconds = observed_profile_seconds, fast_preparation_seconds = fast_preparation_seconds,
    bootstrap_seconds = bootstrap_seconds, total_seconds = total_seconds
  )

  summary <- do.call(rbind, lapply(statistics, function(statistic_name) {
    sunspots_joint_summary_row(
      statistic_name, bootstrap_statistics[[statistic_name]], observed_statistics[[statistic_name]],
      fit,
      list(
        derivative_mc_size = generic_result$diagnostics$derivative_mc_size,
        diagnostics = list(
          Vhat_rcond = generic_result$diagnostics$Vhat_rcond,
          Vhat_condition_number = generic_result$diagnostics$Vhat_condition_number
        ),
        Vhat_min_positive_eigenvalue = numerical_diagnostics$Vhat_min_positive_eigenvalue,
        Vhat_max_positive_eigenvalue = numerical_diagnostics$Vhat_max_positive_eigenvalue,
        inversion_method = "solve",
        regularization_added = 0,
        correction_all_finite = TRUE
      ),
      data, center_indices, settings, timing
    )
  }))
  utils::write.csv(summary, file.path(output_dir, "cycle23_joint_time_space_samplegof_fast_results.csv"), row.names = FALSE)
  for (statistic_name in statistics) {
    utils::write.csv(data.frame(
      replicate = seq_along(bootstrap_statistics[[statistic_name]]),
      statistic = bootstrap_statistics[[statistic_name]]
    ), file.path(output_dir, sprintf("cycle23_joint_time_space_samplegof_fast_bootstrap_%s.csv", statistic_name)), row.names = FALSE)
  }
  utils::write.csv(as.data.frame(timing), file.path(output_dir, "cycle23_joint_time_space_timing.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(list(
    summary = summary, fit = fit, observed_statistics = observed_statistics,
    bootstrap_statistics = bootstrap_statistics, centers = center_indices,
    timing = timing, output_dir = output_dir
  ))
}

parse_sunspots_joint_time_space_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  out <- list()
  integer_keys <- c(
    "B", "n_cores", "observed_profile_n_cores", "center_seed", "dequantization_seed", "derivative_mc_seed", "bootstrap_seed", "derivative_mc_size",
    "bootstrap_block_size", "time_quad_n"
  )
  character_keys <- c(
    "input_csv", "output_dir", "start_date", "end_date", "hemisphere_regression"
  )
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --statistics=ks,cvm --hemisphere_regression=asymmetric|shared ",
        "--B=INTEGER --n_cores=INTEGER --observed_profile_n_cores=INTEGER --n_sample_centers=INTEGER|all ",
        "--derivative_mc_size=INTEGER --bootstrap_block_size=INTEGER --time_quad_n=INTEGER ",
        "--center_seed=INTEGER --dequantization_seed=INTEGER --derivative_mc_seed=INTEGER --bootstrap_seed=INTEGER ",
        "--start_date=YYYY-MM-DD ",
        "--end_date=YYYY-MM-DD --output_dir=PATH ",
        "--allow_boundary_fast=true|false (exploratory only)\\n"
      ))
      quit(save = "no", status = 0L)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop(sprintf("Invalid option: %s", arg))
    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key %in% integer_keys) out[[key]] <- as.integer(value)
    if (identical(key, "n_sample_centers")) {
      out[[key]] <- if (tolower(value) %in% c("all", "inf", "infinity")) Inf else as.integer(value)
    }
    if (identical(key, "statistics")) {
      out[[key]] <- strsplit(tolower(value), ",", fixed = TRUE)[[1L]]
    }
    if (identical(key, "allow_boundary_fast")) {
      value <- tolower(value)
      if (!value %in% c("true", "false")) stop("`allow_boundary_fast` must be true or false.")
      out[[key]] <- identical(value, "true")
    }
    if (key %in% character_keys) out[[key]] <- value
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(run_sunspots_cycle23_joint_time_space_gof, parse_sunspots_joint_time_space_args())
}
