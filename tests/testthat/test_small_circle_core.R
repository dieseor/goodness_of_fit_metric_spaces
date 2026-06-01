library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)

test_that("Small Circle Legendre coefficients recover the uniform case", {
  coeffs <- small_circle_legendre_coefficients(kappa = 0, nu = 0.3, l_max = 8L)
  expect_equal(coeffs$coefficients[[1L]], 1, tolerance = 1e-14)
  expect_equal(coeffs$coefficients[-1L], rep(0, 8), tolerance = 1e-14)
})

test_that("Small Circle profile recovers the uniform S2 law", {
  mu <- c(0, 0, 1)
  omega <- c(1, 0, 0)
  t_grid <- seq(0, pi, length.out = 21)

  observed <- distance_profile_small_circle(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    kappa = 0,
    nu = 0.4,
    distance_type = "geodesic"
  )

  expect_equal(observed, (1 - cos(t_grid)) / 2, tolerance = 1e-12)
})

test_that("Small Circle special cases omega = +/- mu match the axis CDF", {
  mu <- c(0, 0, 1)
  t_grid <- seq(0.15, pi - 0.15, length.out = 11)
  z_threshold <- cos(t_grid)
  g_values <- small_circle_axis_cdf(z_threshold, kappa = 5, nu = 0.3)

  profile_mu <- distance_profile_small_circle(
    omega = mu,
    t_values = t_grid,
    mu = mu,
    kappa = 5,
    nu = 0.3,
    distance_type = "geodesic"
  )
  profile_minus_mu <- distance_profile_small_circle(
    omega = -mu,
    t_values = t_grid,
    mu = mu,
    kappa = 5,
    nu = 0.3,
    distance_type = "geodesic"
  )

  expect_equal(profile_mu, 1 - g_values, tolerance = 5e-8)
  expect_equal(profile_minus_mu, small_circle_axis_cdf(-z_threshold, kappa = 5, nu = 0.3), tolerance = 5e-8)
})

test_that("Small Circle profile is invariant to the order of t_values", {
  mu <- c(0, 0, 1)
  omega <- jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  t_sorted <- seq(0.1, pi - 0.1, length.out = 17)
  permutation <- c(5, 1, 17, 3, 8, 2, 12, 6, 10, 14, 4, 11, 7, 15, 9, 13, 16)
  t_unsorted <- t_sorted[permutation]

  profile_sorted <- distance_profile_small_circle(
    omega = omega,
    t_values = t_sorted,
    mu = mu,
    kappa = 5,
    nu = 0.3,
    distance_type = "geodesic"
  )
  profile_unsorted <- distance_profile_small_circle(
    omega = omega,
    t_values = t_unsorted,
    mu = mu,
    kappa = 5,
    nu = 0.3,
    distance_type = "geodesic"
  )

  expect_equal(profile_unsorted, profile_sorted[permutation], tolerance = 1e-10)
})

test_that("Small Circle Legendre and integral profiles agree across representative cases", {
  mu <- c(0, 0, 1)
  complement <- jp_orthonormal_complement(mu)
  omega_list <- list(
    mu,
    -mu,
    complement[, 1L],
    jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  )
  t_grid <- seq(0.05, pi - 0.05, length.out = 31)
  parameter_grid <- expand.grid(
    kappa = c(0, 1, 5, 20, 50),
    nu = c(0, 0.3, 0.7, 0.9),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(parameter_grid))) {
    comparison <- small_circle_compare_profile_methods(
      mu = mu,
      kappa = parameter_grid$kappa[[i]],
      nu = parameter_grid$nu[[i]],
      omega_list = omega_list,
      t_grid = t_grid,
      distance_type = "geodesic",
      l_max = 200L,
      quad_n = 500L,
      tol = 1e-10
    )

    expect_lte(max(comparison$max_abs_diff), 5e-4)
  }
})

test_that("Small Circle sampler returns unit vectors and MLE respects weighted replication", {
  set.seed(20260531)
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)
  sample <- r_sph_small_circle(n = 160, mu = mu, kappa = 20, nu = 0.7)
  expect_equal(dim(sample), c(160, 3))
  expect_lt(max(abs(sqrt(rowSums(sample^2)) - 1)), 1e-8)

  weights <- rep(c(1, 2, 3, 1), length.out = nrow(sample))
  fit_weighted <- small_circle_mle_s2_weighted(
    sample,
    weights = weights,
    control = list(small_circle_optim_control = list(maxit = 200L, reltol = 1e-9))
  )
  sample_replicated <- sample[rep(seq_len(nrow(sample)), times = weights), , drop = FALSE]
  fit_replicated <- small_circle_mle_s2_weighted(
    sample_replicated,
    control = list(
      small_circle_mle_start_theta = fit_weighted,
      small_circle_optim_control = list(maxit = 200L, reltol = 1e-9)
    )
  )

  expect_true(sum(fit_weighted$mu * fit_replicated$mu) > 0.98)
  expect_equal(fit_weighted$kappa, fit_replicated$kappa, tolerance = 0.5)
  expect_equal(fit_weighted$nu, fit_replicated$nu, tolerance = 0.03)
})