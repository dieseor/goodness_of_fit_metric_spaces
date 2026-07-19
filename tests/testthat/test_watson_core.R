library(testthat)

utils_path <- if (file.exists("utils.R")) "utils.R" else file.path("..", "..", "utils.R")
source(utils_path)

test_that("Watson agrees with Small Circle at nu = 0", {
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)
  omega <- jp_normalize_unit_vector(c(-0.5, 0.2, 0.84), arg_name = "omega", min_length = 3L)
  x <- rbind(mu, -mu, c(1, 0, 0), c(0, 1, 0))
  t_grid <- seq(0.08, pi - 0.08, length.out = 13L)

  for (kappa in c(0, 0.2, 3, 18)) {
    z <- as.numeric(x %*% mu)
    expect_equal(
      watson_axis_density(z, kappa),
      small_circle_axis_density(z, kappa, nu = 0),
      tolerance = 1e-12
    )
    expect_equal(watson_axis_cdf(z, kappa), small_circle_axis_cdf(z, kappa, nu = 0), tolerance = 1e-12)
    expect_equal(
      d_sph_watson_s2(x, mu, kappa),
      2 * small_circle_axis_density(z, kappa, nu = 0),
      tolerance = 1e-12
    )
    expect_equal(
      distance_profile_watson(omega, t_grid, mu, kappa, l_max = 100L, quad_n = 300L),
      distance_profile_small_circle(omega, t_grid, mu, kappa, nu = 0, l_max = 100L, quad_n = 300L),
      tolerance = 2e-9
    )
  }
})

test_that("Watson profiles, CvM grid and sampler match their nu = 0 references", {
  mu <- c(0, 0, 1)
  omega_grid <- rbind(mu, -mu, c(1, 0, 0))
  t_grid <- seq(0.1, pi - 0.1, length.out = 9L)
  set.seed(418)
  x_watson <- r_sph_watson(24, mu, 7)
  set.seed(418)
  x_small_circle <- r_sph_small_circle(24, mu, 7, nu = 0)
  expect_identical(x_watson, x_small_circle)
  expect_equal(
    distance_profile_watson_grid(omega_grid, mu, 7, t_grid, l_max = 100L, quad_n = 300L),
    distance_profile_small_circle_grid(omega_grid, mu, 7, nu = 0, t_grid = t_grid, l_max = 100L, quad_n = 300L),
    tolerance = 2e-9
  )
  expect_equal(
    distance_profile_watson_cvm_grid(x_watson, mu, 7, l_max = 100L, quad_n = 300L),
    distance_profile_small_circle_cvm_grid(x_watson, mu, 7, nu = 0, l_max = 100L, quad_n = 300L),
    tolerance = 2e-9
  )
  coeffs <- watson_legendre_coefficients(7, l_max = 20L, quad_n = 200L)$coefficients
  expect_equal(coeffs[seq.int(2L, length(coeffs), by = 2L)], rep(0, 10L), tolerance = 0)
})

test_that("Watson weighted MLE maximizes the constrained Small Circle likelihood", {
  set.seed(419)
  x <- r_sph_watson(240, c(0, 0, 1), 9)
  weights <- rep(c(1, 2, 3), length.out = nrow(x))
  fit <- watson_mle_s2_weighted(x, weights)
  expect_identical(fit$boundary, "interior")
  expect_equal(
    fit$loglik,
    small_circle_weighted_loglik_s2(fit$mu, fit$kappa, nu = 0, x = x, prob_weights = weights),
    tolerance = 1e-12
  )
  expect_lt(abs(watson_axis_second_moment(fit$kappa) - fit$q_hat), 1e-8)

  isotropic <- rbind(diag(3), -diag(3))
  uniform_fit <- watson_mle_s2_weighted(isotropic)
  expect_identical(uniform_fit$boundary, "uniform")
  expect_equal(uniform_fit$kappa, 0, tolerance = 0)
})
