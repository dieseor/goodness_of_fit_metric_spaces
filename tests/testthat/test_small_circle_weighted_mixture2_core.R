library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)

test_that("weighted small-circle-mixture canonicalization orients mu_3 >= 0 and swaps components", {
  theta <- small_circle_weighted_mixture2_canonicalize_theta(list(
    mu = c(0.1, 0.2, -0.97),
    pi = 0.3,
    kappa1 = 7,
    nu1 = 0.2,
    kappa2 = 11,
    nu2 = 0.6
  ))

  expect_gte(theta$mu[[3L]], 0)
  expect_equal(theta$pi, 0.7, tolerance = 1e-12)
  expect_equal(theta$kappa1, 11, tolerance = 1e-12)
  expect_equal(theta$nu1, 0.6, tolerance = 1e-12)
  expect_equal(theta$kappa2, 7, tolerance = 1e-12)
  expect_equal(theta$nu2, 0.2, tolerance = 1e-12)
  expect_equal(sqrt(sum(theta$mu^2)), 1, tolerance = 1e-12)
})

test_that("weighted small-circle-mixture reduces to the symmetric model", {
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)
  x <- r_sph_small_circle_symmetric_mixture2(n = 20, mu = mu, kappa = 12, nu = 0.35)
  omega <- jp_normalize_unit_vector(c(-0.4, 0.2, 0.894), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 17)

  density_weighted <- d_sph_small_circle_weighted_mixture2_s2(
    x = x,
    mu = mu,
    pi = 0.5,
    kappa1 = 12,
    nu1 = 0.35,
    kappa2 = 12,
    nu2 = 0.35
  )
  density_symmetric <- d_sph_small_circle_symmetric_mixture2_s2(
    x = x,
    mu = mu,
    kappa = 12,
    nu = 0.35
  )
  profile_weighted <- distance_profile_small_circle_weighted_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    pi = 0.5,
    kappa1 = 12,
    nu1 = 0.35,
    kappa2 = 12,
    nu2 = 0.35,
    method = "legendre",
    l_max = 120L,
    quad_n = 400L
  )
  profile_symmetric <- distance_profile_small_circle_symmetric_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    kappa = 12,
    nu = 0.35,
    method = "legendre",
    l_max = 120L,
    quad_n = 400L
  )

  expect_equal(density_weighted, density_symmetric, tolerance = 1e-12)
  expect_equal(profile_weighted, profile_symmetric, tolerance = 1e-10)
})

test_that("weighted small-circle-mixture sampler returns unit vectors and respects weights roughly", {
  set.seed(20260604)
  x <- r_sph_small_circle_weighted_mixture2(
    n = 400,
    mu = c(0, 0, 1),
    pi = 0.7,
    kappa1 = 15,
    nu1 = 0.5,
    kappa2 = 10,
    nu2 = 0.2
  )

  expect_equal(dim(x), c(400, 3))
  expect_lt(max(abs(sqrt(rowSums(x^2)) - 1)), 1e-8)
  expect_gt(mean(x[, 3] > 0), 0.55)
})

test_that("weighted small-circle-mixture weighted MLE agrees approximately with replicated sample fit", {
  set.seed(20260604)
  mu <- jp_normalize_unit_vector(c(0.15, -0.2, 0.968), arg_name = "mu", min_length = 3L)
  x <- r_sph_small_circle_weighted_mixture2(
    n = 120,
    mu = mu,
    pi = 0.65,
    kappa1 = 14,
    nu1 = 0.45,
    kappa2 = 9,
    nu2 = 0.25
  )
  weights <- rep(c(1, 2, 3), length.out = nrow(x))

  fit_weighted <- small_circle_weighted_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = list(
      small_circle_weighted_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9),
      small_circle_weighted_mixture2_n_starts = 4L
    )
  )
  x_rep <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_rep <- small_circle_weighted_mixture2_mle_s2_weighted(
    x = x_rep,
    control = list(
      small_circle_weighted_mixture2_start_theta = fit_weighted,
      small_circle_weighted_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9),
      small_circle_weighted_mixture2_n_starts = 1L
    )
  )

  expect_gt(sum(fit_weighted$mu * fit_rep$mu), 0.95)
  expect_equal(fit_weighted$pi, fit_rep$pi, tolerance = 0.12)
  expect_equal(fit_weighted$nu1, fit_rep$nu1, tolerance = 0.12)
  expect_equal(fit_weighted$nu2, fit_rep$nu2, tolerance = 0.12)
})
