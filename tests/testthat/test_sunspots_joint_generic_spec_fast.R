library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"))
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("real_data", "sunspots", "run_sunspots_cycle23_time_varying_asymmetric_mixture_gof.R"))

joint_generic_fixture <- function(hemisphere_regression = "asymmetric",
                                  n = 36L,
                                  n_sample_centers = Inf,
                                  seed = 20260805L,
                                  derivative_mc_size = 500L,
                                  control = list()) {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  control <- utils::modifyList(control, list(
    hemisphere_regression = hemisphere_regression,
    time_beta_n_starts = 2L,
    time_beta_nelder_mead_control = list(maxit = 500L, reltol = 1e-7),
    optim_control = list(maxit = 100L, reltol = 1e-7),
    profile_l_max = 32L,
    profile_quad_n = 120L,
    time_quad_n = 20L,
    center_block_size = 4L,
    derivative_mc_size = as.integer(derivative_mc_size),
    derivative_mc_seed = as.integer(seed) + 1L,
    distance_profile_backend = "r",
    fast_multiplier_backend = "r",
    fast_multiplier_fuse_ks_cvm = TRUE,
    fast_multiplier_cache_corrections = "auto"
  ))

  theta <- if (identical(hemisphere_regression, "shared")) {
    list(a_N = 0.56, b_N = -0.18, a_S = 0.56, b_S = -0.18, c = 16)
  } else {
    list(a_N = 0.59, b_N = -0.21, a_S = 0.51, b_S = -0.15, c = 16)
  }
  eta <- sunspots_joint_time_canonicalize_eta(list(
    weight1 = 0.42,
    alpha1 = 3.5,
    beta1 = 8.5,
    alpha2 = 8,
    beta2 = 3.2
  ), control = control)
  target_fit <- list(eta_hat = eta, theta_hat = theta)
  par <- sunspots_joint_pack_par(target_fit, control = control)

  simulated <- sunspots_joint_with_seed(seed, sample_sunspots_joint_time_space(n, par, control = control))
  fit <- suppressWarnings(fit_sunspots_cycle23_joint_time_space(
    simulated$x,
    simulated$s,
    hemisphere_regression = hemisphere_regression,
    control = control
  ))
  centers <- sunspots_joint_select_centers(n, n_sample_centers = n_sample_centers, seed = seed + 2L)
  prepared <- sunspots_joint_prepare_centers(
    data = simulated,
    fit = fit,
    center_indices = centers,
    time_quad_n = as.integer(control$time_quad_n),
    l_max = as.integer(control$profile_l_max),
    spatial_quad_n = as.integer(control$profile_quad_n),
    center_block_size = as.integer(control$center_block_size),
    distance_profile_backend = control$distance_profile_backend
  )
  fast <- sunspots_joint_prepare_fast_corrections(
    data = simulated,
    fit = fit,
    centers = prepared$centers,
    derivative_mc_size = as.integer(control$derivative_mc_size),
    seed = as.integer(control$derivative_mc_seed),
    control = control
  )

  list(
    data = simulated,
    fit = fit,
    centers = centers,
    prepared = prepared,
    fast = fast,
    control = control,
    hemisphere_regression = hemisphere_regression
  )
}

test_that("joint generic spec matches custom MLE, distances, observed stats and bootstrap stats (asymmetric)", {
  fixture <- joint_generic_fixture(hemisphere_regression = "asymmetric")
  spec <- make_sunspots_joint_time_space_spec(hemisphere_regression = "asymmetric")

  fit_spec <- spec$fit_theta(
    data = fixture$data,
    weights = NULL,
    null = list(type = "composite"),
    control = fixture$control
  )
  expect_equal(fit_spec$eta_hat$weight1, fixture$fit$eta_hat$weight1, tolerance = 1e-8)
  expect_equal(fit_spec$eta_hat$alpha1, fixture$fit$eta_hat$alpha1, tolerance = 1e-8)
  expect_equal(fit_spec$eta_hat$beta1, fixture$fit$eta_hat$beta1, tolerance = 1e-8)
  expect_equal(fit_spec$eta_hat$alpha2, fixture$fit$eta_hat$alpha2, tolerance = 1e-8)
  expect_equal(fit_spec$eta_hat$beta2, fixture$fit$eta_hat$beta2, tolerance = 1e-8)
  expect_equal(fit_spec$theta_hat$a_N, fixture$fit$theta_hat$a_N, tolerance = 1e-8)
  expect_equal(fit_spec$theta_hat$b_N, fixture$fit$theta_hat$b_N, tolerance = 1e-8)
  expect_equal(fit_spec$theta_hat$a_S, fixture$fit$theta_hat$a_S, tolerance = 1e-8)
  expect_equal(fit_spec$theta_hat$b_S, fixture$fit$theta_hat$b_S, tolerance = 1e-8)
  expect_equal(fit_spec$theta_hat$c, fixture$fit$theta_hat$c, tolerance = 1e-8)

  z <- cbind(fixture$data$x, fixture$data$s)
  d_spec <- spec$distance_matrix(z, z, control = fixture$control)
  d_ref <- matrix(0, nrow = nrow(fixture$data$x), ncol = nrow(fixture$data$x))
  for (j in seq_len(nrow(fixture$data$x))) {
    d_ref[, j] <- sunspots_joint_distance(
      x = fixture$data$x,
      s = fixture$data$s,
      omega = fixture$data$x[j, ],
      center_s = fixture$data$s[j]
    )
  }
  expect_equal(d_spec, d_ref, tolerance = 1e-12)

  result <- multiplier_bootstrap_gof(
    data = z,
    spec = spec,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 10L,
    alpha = 0.05,
    n_cores = 1L,
    seed = 20260831L,
    observed_theta_hat = fixture$fit,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = TRUE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE),
    control = fixture$control,
    distance_profile_backend = "r"
  )

  expect_equal(result$observed$ks$statistic, fixture$prepared$ks_statistic, tolerance = 1e-8)
  expect_equal(result$observed$cvm$statistic, fixture$prepared$cvm_statistic, tolerance = 1e-8)

  custom_boot <- sunspots_time_gof_fast_statistics(
    score_observed = fixture$fast$score_observed,
    centers = fixture$fast$centers,
    statistics = c("ks", "cvm"),
    B = 10L,
    seed = 20260831L,
    n_cores = 1L,
    bootstrap_block_size = 5L
  )
  expect_equal(result$bootstrap$statistics$ks, custom_boot$ks, tolerance = 1e-8)
  expect_equal(result$bootstrap$statistics$cvm, custom_boot$cvm, tolerance = 1e-8)
})

test_that("joint generic spec works for shared and errors explicitly on boundary fast-invalid fit", {
  fixture_shared <- joint_generic_fixture(hemisphere_regression = "shared")
  spec_shared <- make_sunspots_joint_time_space_spec(hemisphere_regression = "shared")
  z_shared <- cbind(fixture_shared$data$x, fixture_shared$data$s)

  result_shared <- multiplier_bootstrap_gof(
    data = z_shared,
    spec = spec_shared,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 6L,
    alpha = 0.05,
    n_cores = 1L,
    seed = 20260901L,
    observed_theta_hat = fixture_shared$fit,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE),
    control = fixture_shared$control,
    distance_profile_backend = "r"
  )
  expect_true(is.numeric(result_shared$bootstrap$statistics$ks))
  expect_true(is.numeric(result_shared$bootstrap$statistics$cvm))

  boundary_fit <- fixture_shared$fit
  boundary_fit$eta_hat$boundary_flags <- list(weight = TRUE, shape_lower = FALSE, shape_upper = FALSE)

  expect_error(
    multiplier_bootstrap_gof(
      data = z_shared,
      spec = spec_shared,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = 3L,
      alpha = 0.05,
      n_cores = 1L,
      seed = 20260902L,
      observed_theta_hat = boundary_fit,
      bootstrap_method = "fast_multiplier",
      keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE),
      control = fixture_shared$control,
      distance_profile_backend = "r"
    ),
    regexp = "boundary"
  )
})
