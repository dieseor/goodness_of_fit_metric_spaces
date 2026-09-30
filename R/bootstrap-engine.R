multiplier_bootstrap_gof <- function(data,
                                     spec,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 5000,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     observed_theta_hat = NULL,
                                     bootstrap_method = c("reestimated", "fast_multiplier"),
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list()) {
  validate_model_spec(spec)
  null <- validate_null_object(null)
  statistics <- normalize_requested_statistics(statistics)
  bootstrap_method <- match.arg(bootstrap_method)
  keep <- normalize_keep_options(keep)

  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  if (!is.finite(B) || B <= 0) {
    stop("`B` must be a strictly positive integer.")
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must belong to (0, 1).")
  }
  if (!is.finite(n_cores) || n_cores <= 0) {
    stop("`n_cores` must be a strictly positive integer.")
  }

  want_ks <- "ks" %in% statistics
  want_cvm <- "cvm" %in% statistics

  data_normalized <- spec_normalize_data(spec, data, control)
  n <- spec_n_obs(spec, data_normalized, control)
  use_lightweight_sample_ks_prep <- identical(bootstrap_method, "fast_multiplier") &&
    want_ks &&
    is_sample_unique_distance_ks_grid(ks_grid %||% list()) &&
    !keep$observed_process
  use_lightweight_cvm_prep <- identical(bootstrap_method, "fast_multiplier") &&
    want_cvm &&
    !keep$observed_process

  if (want_ks && is.null(ks_grid)) {
    ks_grid <- make_sample_unique_distance_ks_grid()
    use_lightweight_sample_ks_prep <- identical(bootstrap_method, "fast_multiplier") &&
      is_sample_unique_distance_ks_grid(ks_grid) &&
      !keep$observed_process
  }

  multiplier_spec <- resolve_multiplier_spec(multipliers)
  scale_factor <- multiplier_spec$mean / multiplier_spec$sd

  start_time <- Sys.time()
  common_observed_start <- proc.time()[["elapsed"]]

  theta_hat <- if (is.null(observed_theta_hat)) {
    spec$fit_theta(
      data = data_normalized,
      weights = NULL,
      null = null,
      control = control
    )
  } else {
    observed_theta_hat
  }

  ks_prep <- if (want_ks) {
    prepare_ks_observed_data(
      data = data_normalized,
      spec = spec,
      theta_hat = theta_hat,
      ks_grid = ks_grid,
      control = control,
      light = use_lightweight_sample_ks_prep,
      share_cvm_statistic = use_lightweight_sample_ks_prep &&
        use_lightweight_cvm_prep && want_cvm
    )
  } else {
    NULL
  }

  cvm_prep <- if (want_cvm) {
    if (isTRUE(use_lightweight_cvm_prep) &&
        isTRUE(ks_prep$light) &&
        identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
      prepare_cvm_observed_data_from_sample_ks(
        data = data_normalized,
        spec = spec,
        theta_hat = theta_hat,
        ks_prep = ks_prep,
        control = control
      )
    } else if (identical(bootstrap_method, "reestimated") &&
               isTRUE(want_ks) && !isTRUE(ks_prep$light) &&
               identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances") &&
               !is.function(spec$cvm_prepare) &&
               !is.function(spec$cvm_bootstrap_stat) &&
               isTRUE(control$reestimated_fuse_ks_cvm %||% TRUE)) {
      prepare_cvm_observed_data_from_full_sample_ks(ks_prep)
    } else {
      prepare_cvm_observed_data(
        data = data_normalized,
        spec = spec,
        theta_hat = theta_hat,
        control = control,
        light = use_lightweight_cvm_prep
      )
    }
  } else {
    NULL
  }

  reestimated_fuse_sample_ks_cvm <- if (identical(bootstrap_method, "reestimated")) {
    reestimated_sample_ks_cvm_fusion_eligible(
      spec = spec,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      control = control
    )
  } else FALSE

  raw_multiplier_matrix <- generate_multiplier_matrix(
    B = B,
    n = n,
    multiplier_spec = multiplier_spec,
    seed = seed
  )
  normalized_multiplier_matrix <- raw_multiplier_matrix / rowMeans(raw_multiplier_matrix)
  common_observed_seconds <- proc.time()[["elapsed"]] - common_observed_start

  n_cores_effective <- min(n_cores, B)

  if (identical(bootstrap_method, "fast_multiplier")) {
    chunk_results <- list(run_fast_multiplier_bootstrap(
      weight_matrix = normalized_multiplier_matrix,
      spec = spec,
      data = data_normalized,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      theta_hat = theta_hat,
      keep_bootstrap_thetas = keep$bootstrap_thetas,
      n_cores = n_cores_effective
    ))
  } else {
    chunk_results <- run_reestimated_bootstrap_chunks(
      weight_matrix = normalized_multiplier_matrix,
      spec = spec,
      data = data_normalized,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      fuse_sample_ks_cvm = reestimated_fuse_sample_ks_cvm,
      keep_bootstrap_thetas = keep$bootstrap_thetas,
      theta_hat = theta_hat,
      n_cores = n_cores_effective
    )
  }

  bootstrap_statistics_internal <- list()
  if (want_ks) {
    bootstrap_statistics_internal$ks <- unlist(lapply(chunk_results, `[[`, "ks"), use.names = FALSE)
  }
  if (want_cvm) {
    bootstrap_statistics_internal$cvm <- unlist(lapply(chunk_results, `[[`, "cvm"), use.names = FALSE)
  }
  bootstrap_theta_internal <- if (keep$bootstrap_thetas &&
    identical(null$type, "composite") &&
    identical(bootstrap_method, "reestimated")) {
    unlist(lapply(chunk_results, `[[`, "theta"), recursive = FALSE, use.names = FALSE)
  } else {
    NULL
  }

  observed_statistics <- list()
  if (want_ks) {
    observed_statistics$ks <- ks_prep$statistic
  }
  if (want_cvm) {
    observed_statistics$cvm <- cvm_prep$statistic
  }

  inference <- compute_inference_summary(
    observed_statistics = observed_statistics,
    bootstrap_statistics = bootstrap_statistics_internal,
    alpha = alpha
  )

  end_time <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  branch_prep_seconds <- sum(vapply(chunk_results, function(x) as.numeric(x$prep_seconds %||% 0), numeric(1)))
  branch_loop_seconds <- sum(vapply(chunk_results, function(x) as.numeric(x$loop_seconds %||% 0), numeric(1)))
  fallback_to_reestimated <- any(vapply(chunk_results, function(x) isTRUE(x$fallback_to_reestimated), logical(1)))
  effective_bootstrap_method <- chunk_results[[1L]]$effective_bootstrap_method %||% bootstrap_method
  fallback_reason <- chunk_results[[1L]]$fallback_reason %||% NA_character_

  result <- list(
    observed = build_observed_output(
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      keep_options = keep
    ),
    bootstrap = list(
      statistics = if (keep$bootstrap_statistics) bootstrap_statistics_internal else NULL,
      theta_star = bootstrap_theta_internal,
      multiplier = list(
        name = multiplier_spec$name,
        mean = multiplier_spec$mean,
        sd = multiplier_spec$sd
      ),
      B = B
    ),
    inference = inference,
    grid = if (want_ks) ks_grid else NULL,
    diagnostics = list(
      n = n,
      B = B,
      alpha = alpha,
      seed = seed,
      n_cores = n_cores_effective,
      null_type = null$type,
      spec_name = spec$name,
      bootstrap_method = bootstrap_method,
      effective_bootstrap_method = effective_bootstrap_method,
      fallback_to_reestimated = fallback_to_reestimated,
      fallback_reason = fallback_reason,
      engine = "multiplier_bootstrap_gof",
      method = "distance_profiles",
      weighted_mle = isTRUE(spec$weighted_mle),
      lightweight_ks_prep = isTRUE(ks_prep$light),
      lightweight_cvm_prep = isTRUE(cvm_prep$light),
      shared_sample_ks_cvm_cache = isTRUE(cvm_prep$shared_with_ks),
      reestimated_fuse_ks_cvm_requested =
        if (identical(bootstrap_method, "reestimated")) {
          isTRUE(control$reestimated_fuse_ks_cvm %||% TRUE)
        } else NA,
      reestimated_fuse_ks_cvm_effective =
        if (identical(bootstrap_method, "reestimated")) {
          reestimated_fuse_sample_ks_cvm
        } else NA,
      ks_prep_bytes = if (!is.null(ks_prep)) as.numeric(object.size(ks_prep)) else NA_real_,
      cvm_prep_bytes = if (!is.null(cvm_prep)) as.numeric(object.size(cvm_prep)) else NA_real_,
      derivative_method = chunk_results[[1L]]$derivative_method %||% NA_character_,
      derivative_method_requested =
        chunk_results[[1L]]$derivative_method_requested %||% NA_character_,
      derivative_method_effective =
        chunk_results[[1L]]$derivative_method_effective %||%
          chunk_results[[1L]]$derivative_method %||% NA_character_,
      derivative_method_selection_source =
        chunk_results[[1L]]$derivative_method_selection_source %||% NA_character_,
      derivative_mc_size = chunk_results[[1L]]$derivative_mc_size %||% NA_integer_,
      derivative_mc_seed = chunk_results[[1L]]$derivative_mc_seed %||% NA_integer_,
      quadrature_algorithm =
        chunk_results[[1L]]$quadrature_algorithm %||% NA_character_,
      quadrature_abs_tol =
        chunk_results[[1L]]$quadrature_abs_tol %||% NA_real_,
      quadrature_max_terms =
        chunk_results[[1L]]$quadrature_max_terms %||% NA_integer_,
      quadrature_initial_upper =
        chunk_results[[1L]]$quadrature_initial_upper %||% NA_real_,
      quadrature_max_upper =
        chunk_results[[1L]]$quadrature_max_upper %||% NA_real_,
      quadrature_tail_consecutive =
        chunk_results[[1L]]$quadrature_tail_consecutive %||% NA_integer_,
      quadrature_eigen_rel_tol =
        chunk_results[[1L]]$quadrature_eigen_rel_tol %||% NA_real_,
      quadrature_clip_tol =
        chunk_results[[1L]]$quadrature_clip_tol %||% NA_real_,
      quadrature_center_evaluations =
        chunk_results[[1L]]$quadrature_center_evaluations %||% NA_integer_,
      quadrature_max_terms_used =
        chunk_results[[1L]]$quadrature_max_terms_used %||% NA_integer_,
      quadrature_max_residual_error_estimate =
        chunk_results[[1L]]$quadrature_max_residual_error_estimate %||% NA_real_,
      quadrature_max_condition_number =
        chunk_results[[1L]]$quadrature_max_condition_number %||% NA_real_,
      quadrature_max_propagated_error_estimate =
        chunk_results[[1L]]$quadrature_max_propagated_error_estimate %||% NA_real_,
      quadrature_max_upper_used =
        chunk_results[[1L]]$quadrature_max_upper_used %||% NA_real_,
      quadrature_max_evaluations =
        chunk_results[[1L]]$quadrature_max_evaluations %||% NA_integer_,
      vhat_method = chunk_results[[1L]]$vhat_method %||% NA_character_,
      S_obs_dim = chunk_results[[1L]]$S_obs_dim %||% NA_integer_,
      Psi_aux_dim = chunk_results[[1L]]$Psi_aux_dim %||% NA_integer_,
      D_ks_dim = chunk_results[[1L]]$D_ks_dim %||% NA_integer_,
      D_cvm_dim = chunk_results[[1L]]$D_cvm_dim %||% NA_integer_,
      fast_ks_mode = chunk_results[[1L]]$fast_ks_mode %||% NA_character_,
      fast_cvm_mode = chunk_results[[1L]]$fast_cvm_mode %||% NA_character_,
      sample_correction_cache_bytes = chunk_results[[1L]]$sample_correction_cache_bytes %||% NA_real_,
      shared_sample_correction_cache = chunk_results[[1L]]$shared_sample_correction_cache %||% NA,
      fast_multiplier_backend_requested =
        chunk_results[[1L]]$fast_multiplier_backend_requested %||% NA_character_,
      fast_multiplier_backend_effective =
        chunk_results[[1L]]$fast_multiplier_backend_effective %||% NA_character_,
      fast_multiplier_cpp_kernel_requested =
        chunk_results[[1L]]$fast_multiplier_cpp_kernel_requested %||% NA_character_,
      fast_multiplier_cpp_kernel_effective =
        chunk_results[[1L]]$fast_multiplier_cpp_kernel_effective %||% NA_character_,
      fast_multiplier_fuse_ks_cvm_requested =
        chunk_results[[1L]]$fast_multiplier_fuse_ks_cvm_requested %||% NA,
      fast_multiplier_fuse_ks_cvm_effective =
        chunk_results[[1L]]$fast_multiplier_fuse_ks_cvm_effective %||% NA,
      fast_multiplier_cache_corrections_requested =
        chunk_results[[1L]]$fast_multiplier_cache_corrections_requested %||% NA_character_,
      fast_multiplier_cache_corrections_effective =
        chunk_results[[1L]]$fast_multiplier_cache_corrections_effective %||% NA,
      Vhat_dim = chunk_results[[1L]]$Vhat_dim %||% NA_integer_,
      score_mean_aux = chunk_results[[1L]]$score_mean_aux %||% NA_real_,
      score_mean_aux_norm = chunk_results[[1L]]$score_mean_aux_norm %||% NA_real_,
      Vhat_eigenvalues = chunk_results[[1L]]$Vhat_eigenvalues %||% NA_real_,
      Vhat_rcond = chunk_results[[1L]]$Vhat_rcond %||% NA_real_,
      Vhat_condition_number = chunk_results[[1L]]$Vhat_condition_number %||% NA_real_,
      fast_parameter_summary = chunk_results[[1L]]$fast_parameter_summary %||% NA_real_,
      common_observed_seconds = common_observed_seconds,
      old_prep_seconds = if (identical(bootstrap_method, "reestimated")) branch_prep_seconds else NA_real_,
      old_loop_seconds = if (identical(bootstrap_method, "reestimated")) branch_loop_seconds else NA_real_,
      old_total_seconds = if (identical(bootstrap_method, "reestimated")) common_observed_seconds + branch_prep_seconds + branch_loop_seconds else NA_real_,
      fast_prep_seconds = if (identical(bootstrap_method, "fast_multiplier")) branch_prep_seconds else NA_real_,
      fast_loop_seconds = if (identical(bootstrap_method, "fast_multiplier")) branch_loop_seconds else NA_real_,
      fast_total_seconds = if (identical(bootstrap_method, "fast_multiplier")) common_observed_seconds + branch_prep_seconds + branch_loop_seconds else NA_real_,
      elapsed_seconds = elapsed_seconds
    )
  )

  class(result) <- c("multiplier_bootstrap_gof_result", "list")
  result
}
.multiplier_bootstrap_gof_backend_implementation <- multiplier_bootstrap_gof
multiplier_bootstrap_gof <- function(data,
                                     spec,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 5000,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     observed_theta_hat = NULL,
                                     bootstrap_method = c("reestimated", "fast_multiplier"),
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list(),
                                     distance_profile_backend = c("r", "cpp")) {
  backend <- normalize_distance_profile_backend(distance_profile_backend)
  spec_name <- as.character(spec$name)
  is_jones_pewsey <- length(spec_name) == 1L && !is.na(spec_name) && grepl("^jp_", spec_name)
  if (is_jones_pewsey && identical(backend, "r")) {
    return(.multiplier_bootstrap_gof_backend_implementation(
      data = data,
      spec = spec,
      null = null,
      statistics = statistics,
      ks_grid = ks_grid,
      B = B,
      alpha = alpha,
      multipliers = multipliers,
      n_cores = n_cores,
      seed = seed,
      observed_theta_hat = observed_theta_hat,
      bootstrap_method = bootstrap_method,
      keep = keep,
      control = control
    ))
  }
  if (identical(backend, "cpp")) {
    assert_distance_profile_cpp_spec_available(spec$name)
  }
  control$distance_profile_backend <- backend
  result <- with_distance_profile_backend(
    backend,
    .multiplier_bootstrap_gof_backend_implementation(
      data = data,
      spec = spec,
      null = null,
      statistics = statistics,
      ks_grid = ks_grid,
      B = B,
      alpha = alpha,
      multipliers = multipliers,
      n_cores = n_cores,
      seed = seed,
      observed_theta_hat = observed_theta_hat,
      bootstrap_method = bootstrap_method,
      keep = keep,
      control = control
    )
  )
  result$diagnostics$distance_profile_backend_requested <- backend
  result$diagnostics$distance_profile_backend_effective <- backend
  result
}
