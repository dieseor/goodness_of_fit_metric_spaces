run_bootstrap_chunk <- function(weight_chunk,
                                spec,
                                data,
                                null,
                                control,
                                scale_factor,
                                ks_prep = NULL,
                                cvm_prep = NULL,
                                want_ks = FALSE,
                                want_cvm = FALSE,
                                fuse_sample_ks_cvm = FALSE,
                                keep_bootstrap_thetas = FALSE,
                                theta_start = NULL,
                                replicate_indices = NULL) {
  scoped_backend <- distance_profile_backend_from_control(control)
  previous_backend <- distance_profile_backend_current()
  if (!identical(scoped_backend, previous_backend)) {
    if (identical(scoped_backend, "cpp")) ensure_distance_profile_cpp_loaded()
    .distance_profile_cpp_state$active_backend <- scoped_backend
    on.exit({
      .distance_profile_cpp_state$active_backend <- previous_backend
    }, add = TRUE)
  }
  n_reps <- nrow(weight_chunk)
  ks_values <- if (want_ks) numeric(n_reps) else NULL
  cvm_values <- if (want_cvm) numeric(n_reps) else NULL
  prep_seconds_total <- 0
  loop_seconds_total <- 0
  theta_values <- if (keep_bootstrap_thetas && identical(null$type, "composite")) {
    vector("list", n_reps)
  } else {
    NULL
  }
  n <- spec_n_obs(spec, data, control)

  for (b in seq_len(n_reps)) {
    replicate_start <- proc.time()[["elapsed"]]
    normalized_weights <- weight_chunk[b, ]
    replicate_index <- if (is.null(replicate_indices)) b else as.integer(replicate_indices[[b]])
    debug_memory_log(
      control,
      sprintf("run_bootstrap_chunk: start replicate %d/%d", replicate_index, n_reps),
      list(
        weight_chunk = weight_chunk,
        normalized_weights = normalized_weights
      )
    )

    theta_star <- NULL
    bootstrap_fit_warnings <- character()
    theta_star_loglik <- NA_real_
    theta_star_convergence <- NA_integer_
    if (identical(null$type, "composite")) {
      bootstrap_control <- control
      if (grepl("^restricted_spiked_normal_", spec$name)) {
        # A finite-sample weighted MLE may attain its supremum at lambda = 0.
        # Allow this boundary value only for restricted-spiked bootstrap refits.
        bootstrap_control$restricted_spiked_allow_boundary_lambda_zero <- TRUE
      }
      if (!is.null(theta_start) && grepl("^jp_", spec$name)) {
        # JP composite bootstrap refits use a warm-started local re-optimization.
        # Together with the logic in jp_mle_s2_weighted(), this keeps the refit
        # on the observed sign branch of psi unless the caller explicitly
        # overrides it. This is a stabilization device for the JP optimizer, not
        # the fully unconstrained composite re-fit.
        bootstrap_control$jp_mle_start_theta <- theta_start
        bootstrap_control$jp_mle_warm_start_only <- TRUE
        bootstrap_control$jp_mle_bootstrap_refit <- TRUE
      } else if (!is.null(theta_start) && grepl("^beta_mixture2_", spec$name)) {
        bootstrap_control$beta_mixture2_start_theta <- theta_start
        bootstrap_control$beta_mixture2_warm_start_only <- TRUE
        bootstrap_control$beta_mixture2_n_starts <- bootstrap_control$beta_mixture2_bootstrap_n_starts %||% 1L
        bootstrap_control$beta_mixture2_optim_control <- bootstrap_control$beta_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^uniform_beta_mixture_", spec$name)) {
        bootstrap_control$uniform_beta_mixture_start_theta <- theta_start
        bootstrap_control$uniform_beta_mixture_warm_start_only <- TRUE
        bootstrap_control$uniform_beta_mixture_n_starts <- bootstrap_control$uniform_beta_mixture_bootstrap_n_starts %||% 1L
        bootstrap_control$uniform_beta_mixture_optim_control <- bootstrap_control$uniform_beta_mixture_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^small_circle_symmetric_mixture2_", spec$name)) {
        bootstrap_control$small_circle_symmetric_mixture2_start_theta <- theta_start
        bootstrap_control$small_circle_symmetric_mixture2_warm_start_only <- TRUE
        bootstrap_control$small_circle_symmetric_mixture2_n_starts <-
          bootstrap_control$small_circle_symmetric_mixture2_bootstrap_n_starts %||% 1L
        bootstrap_control$small_circle_symmetric_mixture2_optim_control <-
          bootstrap_control$small_circle_symmetric_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^small_circle_weighted_mixture2_", spec$name)) {
        bootstrap_control$small_circle_weighted_mixture2_start_theta <- theta_start
        bootstrap_control$small_circle_weighted_mixture2_warm_start_only <- TRUE
        bootstrap_control$small_circle_weighted_mixture2_n_starts <-
          bootstrap_control$small_circle_weighted_mixture2_bootstrap_n_starts %||% 1L
        bootstrap_control$small_circle_weighted_mixture2_optim_control <-
          bootstrap_control$small_circle_weighted_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^axial_truncnorm_mixture2_", spec$name)) {
        bootstrap_control$axial_truncnorm_mixture2_start_theta <- theta_start
        bootstrap_control$axial_truncnorm_mixture2_optim_control <-
          bootstrap_control$axial_truncnorm_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^logitnormal_mixture2_", spec$name)) {
        bootstrap_control$logitnormal_mixture2_start_theta <- theta_start
        bootstrap_control$logitnormal_mixture2_warm_start_only <- TRUE
      } else if (!is.null(theta_start) && is.null(bootstrap_control$jp_mle_start_theta)) {
        bootstrap_control$jp_mle_start_theta <- theta_start
      }
      theta_star <- tryCatch(
        withCallingHandlers(
          spec$fit_theta(
            data = data,
            weights = normalized_weights,
            null = null,
            control = bootstrap_control
          ),
          warning = function(w) {
            bootstrap_fit_warnings <<- c(bootstrap_fit_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          bootstrap_fit_warnings <<- c(
            bootstrap_fit_warnings,
            paste0("Bootstrap MLE error: ", conditionMessage(e))
          )
          NULL
        }
      )
      theta_star_loglik <- as.numeric(theta_star$loglik %||% NA_real_)
      theta_star_convergence <- as.integer(theta_star$opt$convergence %||% NA_integer_)
      if (grepl("^small_circle_symmetric_mixture2_", spec$name) &&
          (is.null(theta_star) ||
             any(!is.finite(as.numeric(c(theta_star$mu, theta_star$kappa, theta_star$nu)))))) {
        bootstrap_fit_warnings <- c(
          bootstrap_fit_warnings,
          "Bootstrap MLE returned non-finite theta_star; falling back to observed theta_hat."
        )
        theta_star <- theta_start
        theta_star_loglik <- as.numeric(theta_start$loglik %||% NA_real_)
        theta_star_convergence <- as.integer(theta_start$opt$convergence %||% NA_integer_)
      }
      if (grepl("^axial_truncnorm_mixture2_", spec$name) &&
          (is.null(theta_star) ||
             any(!is.finite(as.numeric(c(
               theta_star$pi,
               theta_star$kappa1,
               theta_star$nu1,
               theta_star$kappa2,
               theta_star$nu2
             )))))) {
        bootstrap_fit_warnings <- c(
          bootstrap_fit_warnings,
          "Bootstrap axial MLE returned an invalid theta_star; falling back to observed theta_hat."
        )
        theta_star <- theta_start
        theta_star_loglik <- as.numeric(theta_start$loglik %||% NA_real_)
        theta_star_convergence <- as.integer(theta_start$opt$convergence %||% NA_integer_)
      }
      if (!is.null(theta_values)) {
        theta_values[[b]] <- theta_star
      }
    }

    if (isTRUE(fuse_sample_ks_cvm)) {
      f_star_sample <- compute_weighted_sample_profile_matrix(
        order_matrix = ks_prep$order_matrix,
        rank_linear_index = ks_prep$rank_linear_index,
        normalized_weights = normalized_weights
      )
      if (identical(null$type, "simple")) {
        process_star_sample <- scale_factor * sqrt(n) * (
          f_star_sample - ks_prep$empirical_profile
        )
      } else {
        f_theta_star_sample <- compute_theoretical_sample_profile_matrix(
          spec = spec,
          data = data,
          distance_matrix = ks_prep$distance_matrix,
          theta = theta_star,
          control = control
        )
        process_star_sample <- scale_factor * sqrt(n) * (
          (f_star_sample - f_theta_star_sample) -
            (ks_prep$empirical_profile - ks_prep$theoretical_profile)
        )
      }
      ks_values[b] <- max(abs(process_star_sample))
      cvm_values[b] <- mean(process_star_sample^2)
    } else if (want_ks) {
      if (identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
        ks_values[b] <- compute_ks_sample_stat_blocked(
          spec = spec,
          normalized_data = data,
          ks_prep = ks_prep,
          normalized_weights = normalized_weights,
          scale_factor = scale_factor,
          theta_star = if (identical(null$type, "composite")) theta_star else NULL,
          null_type = null$type,
          control = control
        )
      } else {
        f_star_grid <- compute_grid_weighted_profile(
          distance_matrix = ks_prep$distance_matrix,
          t_grid = ks_prep$t_grid,
          normalized_weights = normalized_weights,
          sorted_distance_matrix = ks_prep$sorted_distance_matrix,
          order_matrix = ks_prep$order_matrix,
          threshold_index_matrix = ks_prep$threshold_index_matrix
        )

        if (identical(null$type, "simple")) {
          process_star_grid <- scale_factor * sqrt(n) * (f_star_grid - ks_prep$empirical_profile)
        } else {
          f_theta_star <- compute_theoretical_profile_matrix(
            spec = spec,
            omega_grid = ks_prep$omega_grid,
            t_grid = ks_prep$t_grid,
            theta = theta_star,
            control = control
          )
          process_star_grid <- scale_factor * sqrt(n) * (
            (f_star_grid - f_theta_star) -
              (ks_prep$empirical_profile - ks_prep$theoretical_profile)
          )
        }

        ks_values[b] <- max(abs(process_star_grid))
      }
    }

    if (want_cvm && !isTRUE(fuse_sample_ks_cvm)) {
      cvm_control <- utils::modifyList(
        control,
        list(
          small_circle_symmetric_mixture2_bootstrap_replicate_index = replicate_index,
          small_circle_symmetric_mixture2_bootstrap_warnings = bootstrap_fit_warnings,
          small_circle_symmetric_mixture2_bootstrap_loglik = theta_star_loglik,
          small_circle_symmetric_mixture2_bootstrap_convergence = theta_star_convergence
        )
      )
      cvm_stat_fast <- spec_cvm_bootstrap_stat(
        spec = spec,
        data = data,
        normalized_weights = normalized_weights,
        theta_star = theta_star,
        cvm_prep = cvm_prep,
        null = null,
        control = cvm_control,
        scale_factor = scale_factor
      )
      if (!is.null(cvm_stat_fast)) {
        cvm_values[b] <- cvm_stat_fast
      } else {
        f_star_sample <- compute_weighted_sample_profile_matrix(
          order_matrix = cvm_prep$order_matrix,
          rank_linear_index = cvm_prep$rank_linear_index,
          normalized_weights = normalized_weights
        )
        debug_memory_log(control, sprintf("run_bootstrap_chunk: replicate %d after f_star_sample", b), list(f_star_sample = f_star_sample))

        if (identical(null$type, "simple")) {
          process_star_sample <- scale_factor * sqrt(n) * (f_star_sample - cvm_prep$empirical_profile)
        } else {
          f_theta_star_sample <- compute_theoretical_sample_profile_matrix(
            spec = spec,
            data = data,
            distance_matrix = cvm_prep$distance_matrix,
            theta = theta_star,
            control = cvm_control
          )
          debug_memory_log(
            control,
            sprintf("run_bootstrap_chunk: replicate %d after f_theta_star_sample", b),
            list(f_theta_star_sample = f_theta_star_sample)
          )
          process_star_sample <- scale_factor * sqrt(n) * (
            (f_star_sample - f_theta_star_sample) -
              (cvm_prep$empirical_profile - cvm_prep$theoretical_profile)
          )
        }
        debug_memory_log(
          control,
          sprintf("run_bootstrap_chunk: replicate %d after process_star_sample", b),
          list(process_star_sample = process_star_sample)
        )

        cvm_values[b] <- mean(process_star_sample^2)
      }
    }
    loop_seconds_total <- loop_seconds_total + (proc.time()[["elapsed"]] - replicate_start)
  }

  list(
    ks = ks_values,
    cvm = cvm_values,
    theta = theta_values,
    prep_seconds = prep_seconds_total,
    loop_seconds = loop_seconds_total
  )
}
compute_inference_summary <- function(observed_statistics, bootstrap_statistics, alpha) {
  output <- list()

  for (stat_name in names(observed_statistics)) {
    observed_value <- observed_statistics[[stat_name]]
    bootstrap_values <- bootstrap_statistics[[stat_name]]

    critical_value <- as.numeric(stats::quantile(
      bootstrap_values,
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    ))
    p_value <- (1 + sum(bootstrap_values >= observed_value)) / (length(bootstrap_values) + 1)

    output[[stat_name]] <- list(
      observed = observed_value,
      critical_value = critical_value,
      p_value = p_value,
      reject = isTRUE(p_value <= alpha)
    )
  }

  output
}
build_observed_output <- function(theta_hat, ks_prep, cvm_prep, keep_options) {
  output <- list(theta_hat = theta_hat)

  if (!is.null(ks_prep)) {
    output$ks <- list(statistic = ks_prep$statistic)
    if (keep_options$observed_process) {
      output$ks$process_matrix <- ks_prep$process_matrix
      output$ks$empirical_profile <- ks_prep$empirical_profile
      output$ks$theoretical_profile <- ks_prep$theoretical_profile
    }
  }

  if (!is.null(cvm_prep)) {
    output$cvm <- list(statistic = cvm_prep$statistic)
    if (keep_options$observed_process && !is.null(cvm_prep$process_matrix)) {
      output$cvm$process_matrix <- cvm_prep$process_matrix
      output$cvm$distance_matrix <- cvm_prep$distance_matrix
      output$cvm$empirical_profile <- cvm_prep$empirical_profile
      output$cvm$theoretical_profile <- cvm_prep$theoretical_profile
    }
  }

  output
}
