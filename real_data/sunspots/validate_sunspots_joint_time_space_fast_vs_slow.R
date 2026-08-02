#!/usr/bin/env Rscript

# Small-sample validation: the fast multiplier law is compared with a
# parametric bootstrap that simulates from the fitted joint law and refits both
# temporal and conditional components in every replicate.

resolve_sunspots_joint_validation_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_sunspots_joint_validation_path(
  "real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"
))

validate_sunspots_joint_time_space_fast_vs_slow <- function(
    output_dir = file.path("real_data", "sunspots", "output", "joint_time_space_fast_vs_slow_validation"),
    hemisphere_regression = "asymmetric", n = 100L, n_sample_centers = 12L,
    B = 100L, derivative_mc_size = 3000L, seed = 20260711L,
    n_cores = 1L, time_quad_n = 32L, control = list()) {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  control <- utils::modifyList(control, list(hemisphere_regression = hemisphere_regression))
  theta <- if (identical(hemisphere_regression, "shared")) {
    list(a_N = 0.58, b_N = -0.22, a_S = 0.58, b_S = -0.22, c = 18)
  } else {
    list(a_N = 0.60, b_N = -0.23, a_S = 0.52, b_S = -0.16, c = 18)
  }
  eta <- sunspots_joint_time_canonicalize_eta(list(
    weight1 = 0.45, alpha1 = 4.5, beta1 = 10, alpha2 = 10, beta2 = 4.5
  ), control)
  target_fit <- list(eta_hat = eta, theta_hat = theta)
  target_par <- sunspots_joint_pack_par(target_fit, control)
  observed_data <- sunspots_joint_with_seed(seed, sample_sunspots_joint_time_space(n, target_par, control))
  observed_fit <- fit_sunspots_cycle23_joint_time_space(
    observed_data$x, observed_data$s, hemisphere_regression = hemisphere_regression, control = control
  )
  if (any(unlist(observed_fit$eta_hat$boundary_flags, use.names = FALSE))) {
    stop("The validation sample has a temporal boundary MLE; rerun with a different `seed`.")
  }
  center_indices <- sunspots_joint_select_centers(n, n_sample_centers, seed + 1L)
  prepared <- sunspots_joint_prepare_centers(
    observed_data, observed_fit, center_indices, time_quad_n = time_quad_n,
    l_max = as.integer(control$profile_l_max %||% 60L),
    spatial_quad_n = as.integer(control$profile_quad_n %||% 200L),
    distance_profile_backend = control$distance_profile_backend %||% "r"
  )
  fast <- sunspots_joint_prepare_fast_corrections(
    observed_data, observed_fit, prepared$centers, derivative_mc_size = derivative_mc_size,
    seed = seed + 2L, control = control
  )
  fast_statistics <- sunspots_time_gof_fast_statistics(
    fast$score_observed, fast$centers, statistics = c("ks", "cvm"), B = B,
    seed = seed + 3L, n_cores = n_cores, bootstrap_block_size = 10L
  )
  slow_statistics <- sunspots_joint_slow_reestimated_statistics(
    n = n, par = sunspots_joint_pack_par(observed_fit, control), center_indices = center_indices,
    B = B, seed = seed + 4L, hemisphere_regression = hemisphere_regression,
    statistics = c("ks", "cvm"), time_quad_n = time_quad_n,
    l_max = as.integer(control$profile_l_max %||% 60L),
    spatial_quad_n = as.integer(control$profile_quad_n %||% 200L),
    distance_profile_backend = control$distance_profile_backend %||% "r", control = control
  )
  slow_boundary_flags <- attr(slow_statistics, "temporal_boundary_flags")
  slow_boundary_replicates <- sum(rowSums(slow_boundary_flags) > 0L)
  observed_statistics <- c(ks = prepared$ks_statistic, cvm = prepared$cvm_statistic)
  summary <- do.call(rbind, lapply(c("ks", "cvm"), function(statistic_name) {
    fast_values <- fast_statistics[[statistic_name]]
    slow_values <- slow_statistics[[statistic_name]]
    data.frame(
      hemisphere_regression = hemisphere_regression, statistic_type = statistic_name,
      observed_statistic = observed_statistics[[statistic_name]], B = B,
      slow_temporal_boundary_replicates = slow_boundary_replicates,
      fast_p_value = (1 + sum(fast_values >= observed_statistics[[statistic_name]])) / (B + 1),
      slow_p_value = (1 + sum(slow_values >= observed_statistics[[statistic_name]])) / (B + 1),
      fast_q05 = stats::quantile(fast_values, 0.05, names = FALSE),
      slow_q05 = stats::quantile(slow_values, 0.05, names = FALSE),
      fast_q50 = stats::quantile(fast_values, 0.50, names = FALSE),
      slow_q50 = stats::quantile(slow_values, 0.50, names = FALSE),
      fast_q95 = stats::quantile(fast_values, 0.95, names = FALSE),
      slow_q95 = stats::quantile(slow_values, 0.95, names = FALSE),
      stringsAsFactors = FALSE
    )
  }))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, file.path(output_dir, "joint_time_space_fast_vs_slow_summary.csv"), row.names = FALSE)
  utils::write.csv(data.frame(replicate = seq_len(B), fast_ks = fast_statistics$ks,
                              fast_cvm = fast_statistics$cvm,
                              slow_ks = slow_statistics$ks,
                              slow_cvm = slow_statistics$cvm),
                    file.path(output_dir, "joint_time_space_fast_vs_slow_replicates.csv"), row.names = FALSE)
  invisible(list(summary = summary, fast = fast_statistics, slow = slow_statistics,
                 observed = observed_statistics, output_dir = output_dir))
}

if (sys.nframe() == 0L) {
  validate_sunspots_joint_time_space_fast_vs_slow()
}
