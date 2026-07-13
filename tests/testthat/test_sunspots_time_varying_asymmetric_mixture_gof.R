library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "sunspots", "run_sunspots_cycle23_time_varying_asymmetric_mixture_gof.R"))

test_that("conditional averaged profile is a valid rotational probability law", {
  theta <- list(a_N = 0.55, b_N = -0.20, a_S = 0.50, b_S = -0.15, c = 18)
  coefficients <- sunspots_time_gof_average_legendre_coefficients(
    theta, u_reference = seq(0.05, 0.95, length.out = 30), l_max = 60L, quad_n = 250L
  )
  expect_equal(coefficients[[1L]], 1, tolerance = 1e-10)
  profile <- sunspots_time_gof_profile_row(c(0, 0, 1), seq(0, pi, length.out = 101), coefficients)
  expect_true(all(profile >= 0 & profile <= 1))
  expect_true(all(diff(profile) >= -1e-10))
  expect_equal(profile[[1L]], 0, tolerance = 1e-12)
  expect_equal(profile[[101L]], 1, tolerance = 1e-12)
})

test_that("a common time point recovers the static symmetric-mixture profile", {
  theta <- list(a_N = 0.47, b_N = -0.10, a_S = 0.47, b_S = -0.10, c = 18)
  coefficients <- sunspots_time_gof_average_legendre_coefficients(
    theta, u_reference = rep(0.5, 30), l_max = 100L, quad_n = 400L
  )
  omega <- c(0.3, -0.4, sqrt(0.75))
  thresholds <- seq(0.05, pi - 0.05, length.out = 31)
  observed <- sunspots_time_gof_profile_row(omega, thresholds, coefficients)
  expected <- distance_profile_small_circle_symmetric_mixture2(
    omega = omega,
    t_values = thresholds,
    mu = c(0, 0, 1),
    kappa = theta$c,
    nu = theta$a_N + 0.5 * theta$b_N,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 100L,
    quad_n = 400L
  )
  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("fast conditional sample KS returns finite bootstrap inference", {
  set.seed(123)
  theta <- list(a_N = 0.50, b_N = -0.15, a_S = 0.48, b_S = -0.12, c = 15)
  u <- (seq_len(36) - 0.5) / 36
  par <- sunspots_time_varying_pack_theta(theta)
  simulated <- sample_sunspots_time_gof_conditional(36L, u, par)
  data <- normalize_sunspots_time_gof_data(simulated$x, simulated$u)
  centers <- sunspots_time_gof_prepare_centers(data, theta, center_indices = 1:6, l_max = 40L, quad_n = 200L)
  fast <- sunspots_time_gof_prepare_fast_corrections(
    data, theta, centers$centers, derivative_mc_size = 120L, seed = 321L
  )
  values <- sunspots_time_gof_fast_ks_statistics(
    fast$score_observed, fast$centers, B = 6L, seed = 456L, n_cores = 1L, bootstrap_block_size = 2L
  )
  expect_equal(length(values), 6L)
  expect_true(all(is.finite(values)))
  expect_true(is.finite(centers$statistic))
})

test_that("fast conditional sample CvM shares the KS bootstrap preparation", {
  set.seed(124)
  theta <- list(a_N = 0.50, b_N = -0.15, a_S = 0.48, b_S = -0.12, c = 15)
  u <- (seq_len(30) - 0.5) / 30
  simulated <- sample_sunspots_time_gof_conditional(
    30L, u, sunspots_time_varying_pack_theta(theta)
  )
  data <- normalize_sunspots_time_gof_data(simulated$x, simulated$u)
  centers <- sunspots_time_gof_prepare_centers(data, theta, center_indices = 1:8, l_max = 40L, quad_n = 200L)
  fast <- sunspots_time_gof_prepare_fast_corrections(
    data, theta, centers$centers, derivative_mc_size = 120L, seed = 322L
  )
  both <- sunspots_time_gof_fast_statistics(
    fast$score_observed, fast$centers, statistics = c("ks", "cvm"),
    B = 6L, seed = 457L, n_cores = 1L, bootstrap_block_size = 2L
  )
  ks_only <- sunspots_time_gof_fast_ks_statistics(
    fast$score_observed, fast$centers,
    B = 6L, seed = 457L, n_cores = 1L, bootstrap_block_size = 2L
  )

  expect_equal(both$ks, ks_only, tolerance = 1e-12)
  expect_equal(length(both$cvm), 6L)
  expect_true(all(is.finite(both$cvm)))
  expect_true(is.finite(centers$cvm_statistic))
})

test_that("shared hemisphere regression uses a three-dimensional fast correction", {
  set.seed(125)
  theta <- list(a_N = 0.50, b_N = -0.15, a_S = 0.50, b_S = -0.15, c = 15)
  control <- list(hemisphere_regression = "shared")
  u <- (seq_len(32) - 0.5) / 32
  par <- sunspots_time_varying_pack_theta(theta, hemisphere_regression = "shared")
  simulated <- sample_sunspots_time_gof_conditional(32L, u, par, control)
  data <- normalize_sunspots_time_gof_data(simulated$x, simulated$u)
  centers <- sunspots_time_gof_prepare_centers(data, theta, center_indices = 1:8, l_max = 40L, quad_n = 200L)
  fast <- sunspots_time_gof_prepare_fast_corrections(
    data, theta, centers$centers, derivative_mc_size = 120L, seed = 323L, control = control
  )
  values <- sunspots_time_gof_fast_statistics(
    fast$score_observed, fast$centers, statistics = c("ks", "cvm"),
    B = 6L, seed = 458L, n_cores = 1L, bootstrap_block_size = 2L
  )

  expect_identical(dim(fast$score_observed), c(32L, 3L))
  expect_identical(dim(fast$vhat), c(3L, 3L))
  expect_true(all(is.finite(values$ks)))
  expect_true(all(is.finite(values$cvm)))
})

test_that("shared fast score matches the three-coordinate numerical derivative", {
  theta <- list(a_N = 0.52, b_N = -0.16, a_S = 0.52, b_S = -0.16, c = 15)
  control <- list(hemisphere_regression = "shared")
  par <- sunspots_time_varying_pack_theta(theta, hemisphere_regression = "shared")
  set.seed(126)
  u <- (seq_len(16) - 0.5) / 16
  simulated <- sample_sunspots_time_gof_conditional(16L, u, par, control)
  data <- normalize_sunspots_time_gof_data(simulated$x, simulated$u)
  analytic <- sunspots_time_gof_score_matrix(data, par, control)
  step <- 1e-6
  numeric <- vapply(seq_along(par), function(index) {
    plus <- par
    minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (
      sunspots_time_varying_log_density(
        data$x, data$u,
        sunspots_time_varying_unpack_par(plus, hemisphere_regression = "shared")
      ) - sunspots_time_varying_log_density(
        data$x, data$u,
        sunspots_time_varying_unpack_par(minus, hemisphere_regression = "shared")
      )
    ) / (2 * step)
  }, numeric(nrow(data$x)))

  expect_equal(as.numeric(analytic), as.numeric(numeric), tolerance = 1e-5)
})
