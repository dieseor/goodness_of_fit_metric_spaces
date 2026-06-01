library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)

test_that("generic rotational Legendre coefficients recover the uniform law", {
  coeffs <- rotational_legendre_coefficients(
    density_h = function(z) rep(1, length(z)),
    Lmax = 8L,
    quad_n = 200L
  )

  expect_equal(coeffs$coefficients[[1L]], 1, tolerance = 1e-14)
  expect_equal(coeffs$coefficients[-1L], rep(0, 8L), tolerance = 1e-12)
})

test_that("rotational beta-mixture coefficients and special profiles are correct", {
  theta <- rotational_beta_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.4,
    alpha1 = 2,
    beta1 = 7,
    alpha2 = 9,
    beta2 = 3
  ))
  coeffs <- rotational_beta_mixture2_legendre_coefficients(theta, l_max = 40L, quad_n = 400L)
  expect_lt(coeffs$a0_error, 1e-10)

  t_grid <- seq(0, pi, length.out = 31)
  profile_mu <- distance_profile_rotational_beta_mixture2(
    omega = theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 120L,
    quad_n = 500L
  )
  profile_minus_mu <- distance_profile_rotational_beta_mixture2(
    omega = -theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 120L,
    quad_n = 500L
  )
  thresholds <- cos(t_grid)
  expect_equal(
    profile_mu,
    1 - rotational_beta_mixture2_cdf_y(
      y = (thresholds + 1) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    ),
    tolerance = 1e-9
  )
  expect_equal(
    profile_minus_mu,
    rotational_beta_mixture2_cdf_y(
      y = (1 - thresholds) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    ),
    tolerance = 1e-9
  )
  expect_equal(profile_mu[[1L]], 0, tolerance = 1e-12)
  expect_equal(profile_mu[[length(profile_mu)]], 1, tolerance = 1e-12)
  expect_true(all(diff(profile_mu) >= -1e-10))
  expect_true(all(profile_mu >= -1e-12 & profile_mu <= 1 + 1e-12))
})

test_that("rotational beta-mixture Legendre and integral profiles agree", {
  theta <- rotational_beta_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.45,
    alpha1 = 3,
    beta1 = 10,
    alpha2 = 12,
    beta2 = 4
  ))
  omega <- jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 25)

  legendre <- distance_profile_rotational_beta_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    method = "legendre",
    l_max = 150L,
    quad_n = 600L
  )
  integral <- distance_profile_rotational_beta_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    method = "integral",
    quad_n = 600L
  )

  expect_lt(max(abs(legendre - integral)), 5e-4)
})

test_that("rotational beta-mixture sampler and weighted MLE behave coherently", {
  set.seed(20260601)
  theta <- rotational_beta_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
    weight1 = 0.4,
    alpha1 = 2.5,
    beta1 = 8,
    alpha2 = 9,
    beta2 = 2.5
  ))
  x <- r_sph_rotational_beta_mixture2(
    n = 180,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )
  expect_lt(max(abs(sqrt(rowSums(x^2)) - 1)), 1e-8)

  weights <- rep(c(1, 2, 3), length.out = nrow(x))
  fit_weighted <- rotational_beta_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = list(rotational_beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
  )
  x_rep <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_rep <- rotational_beta_mixture2_mle_s2_weighted(
    x = x_rep,
    control = list(
      rotational_beta_mixture2_start_theta = fit_weighted,
      rotational_beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(sum(fit_weighted$mu * fit_rep$mu) > 0.95)
  expect_equal(fit_weighted$weight1, fit_rep$weight1, tolerance = 0.08)
  expect_equal(fit_weighted$alpha1, fit_rep$alpha1, tolerance = 1.5)
  expect_equal(fit_weighted$beta1, fit_rep$beta1, tolerance = 2.0)

  y <- (as.numeric(x %*% theta$mu) + 1) / 2
  ecdf_y <- stats::ecdf(y)
  grid <- seq(0.05, 0.95, length.out = 41)
  fitted_cdf <- rotational_beta_mixture2_cdf_y(
    y = grid,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )
  expect_lt(max(abs(ecdf_y(grid) - fitted_cdf)), 0.12)
})

test_that("rotational logit-normal-mixture coefficients and special profiles are correct", {
  theta <- rotational_logitnormal_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.45,
    mean1 = -1.1,
    sd1 = 0.45,
    mean2 = 1.0,
    sd2 = 0.55
  ))
  coeffs <- rotational_logitnormal_mixture2_legendre_coefficients(theta, l_max = 40L, quad_n = 400L)
  expect_lt(coeffs$a0_error, 1e-10)

  t_grid <- seq(0, pi, length.out = 31)
  profile_mu <- distance_profile_rotational_logitnormal_mixture2(
    omega = theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 120L,
    quad_n = 500L
  )
  profile_minus_mu <- distance_profile_rotational_logitnormal_mixture2(
    omega = -theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 120L,
    quad_n = 500L
  )
  thresholds <- cos(t_grid)
  expect_equal(
    profile_mu,
    1 - rotational_logitnormal_mixture2_cdf_y(
      y = (thresholds + 1) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2
    ),
    tolerance = 1e-9
  )
  expect_equal(
    profile_minus_mu,
    rotational_logitnormal_mixture2_cdf_y(
      y = (1 - thresholds) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2
    ),
    tolerance = 1e-9
  )
  expect_equal(profile_mu[[1L]], 0, tolerance = 1e-12)
  expect_equal(profile_mu[[length(profile_mu)]], 1, tolerance = 1e-12)
  expect_true(all(diff(profile_mu) >= -1e-10))
  expect_true(all(profile_mu >= -1e-12 & profile_mu <= 1 + 1e-12))
})

test_that("rotational logit-normal-mixture Legendre and integral profiles agree", {
  theta <- rotational_logitnormal_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.35,
    mean1 = -0.6,
    sd1 = 0.4,
    mean2 = 1.3,
    sd2 = 0.5
  ))
  omega <- jp_normalize_unit_vector(c(-0.45, 0.2, 0.87), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 25)

  legendre <- distance_profile_rotational_logitnormal_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2,
    method = "legendre",
    l_max = 150L,
    quad_n = 600L
  )
  integral <- distance_profile_rotational_logitnormal_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2,
    method = "integral",
    quad_n = 600L
  )

  expect_lt(max(abs(legendre - integral)), 5e-4)
})

test_that("rotational logit-normal-mixture sampler and weighted MLE behave coherently", {
  set.seed(20260602)
  theta <- rotational_logitnormal_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(-0.25, 0.3, 0.92), arg_name = "mu", min_length = 3L),
    weight1 = 0.45,
    mean1 = -1.0,
    sd1 = 0.45,
    mean2 = 1.15,
    sd2 = 0.55
  ))
  x <- r_sph_rotational_logitnormal_mixture2(
    n = 180,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2
  )
  expect_lt(max(abs(sqrt(rowSums(x^2)) - 1)), 1e-8)

  weights <- rep(c(1, 3, 2), length.out = nrow(x))
  fit_weighted <- rotational_logitnormal_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = list(rotational_logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
  )
  x_rep <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_rep <- rotational_logitnormal_mixture2_mle_s2_weighted(
    x = x_rep,
    control = list(
      rotational_logitnormal_mixture2_start_theta = fit_weighted,
      rotational_logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(sum(fit_weighted$mu * fit_rep$mu) > 0.95)
  expect_equal(fit_weighted$weight1, fit_rep$weight1, tolerance = 0.08)
  expect_equal(fit_weighted$mean1, fit_rep$mean1, tolerance = 0.25)
  expect_equal(fit_weighted$mean2, fit_rep$mean2, tolerance = 0.25)
  cdf_grid <- seq(0.1, 0.9, length.out = 21)
  cdf_weighted <- rotational_logitnormal_mixture2_cdf_y(
    y = cdf_grid,
    weight1 = fit_weighted$weight1,
    mean1 = fit_weighted$mean1,
    sd1 = fit_weighted$sd1,
    mean2 = fit_weighted$mean2,
    sd2 = fit_weighted$sd2
  )
  cdf_rep <- rotational_logitnormal_mixture2_cdf_y(
    y = cdf_grid,
    weight1 = fit_rep$weight1,
    mean1 = fit_rep$mean1,
    sd1 = fit_rep$sd1,
    mean2 = fit_rep$mean2,
    sd2 = fit_rep$sd2
  )
  expect_lt(max(abs(cdf_weighted - cdf_rep)), 0.15)

  y <- (as.numeric(x %*% theta$mu) + 1) / 2
  ecdf_y <- stats::ecdf(y)
  grid <- seq(0.05, 0.95, length.out = 41)
  fitted_cdf <- rotational_logitnormal_mixture2_cdf_y(
    y = grid,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2
  )
  expect_lt(max(abs(ecdf_y(grid) - fitted_cdf)), 0.12)
})

test_that("mixture canonicalization swaps labels correctly", {
  beta_theta <- rotational_beta_mixture2_canonicalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.7,
    alpha1 = 8,
    beta1 = 2,
    alpha2 = 2,
    beta2 = 8
  ))
  expect_lt(beta_theta$alpha1 / (beta_theta$alpha1 + beta_theta$beta1),
            beta_theta$alpha2 / (beta_theta$alpha2 + beta_theta$beta2))
  expect_equal(beta_theta$weight1, 0.3, tolerance = 1e-12)

  logit_theta <- rotational_logitnormal_mixture2_canonicalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.7,
    mean1 = 2,
    sd1 = 0.4,
    mean2 = -1,
    sd2 = 0.5
  ))
  expect_lt(logit_theta$mean1, logit_theta$mean2)
  expect_equal(logit_theta$weight1, 0.3, tolerance = 1e-12)
})
