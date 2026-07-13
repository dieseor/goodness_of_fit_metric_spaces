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
  theta <- beta_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.4,
    alpha1 = 2,
    beta1 = 7,
    alpha2 = 9,
    beta2 = 3
  ))
  coeffs <- beta_mixture2_legendre_coefficients(theta, l_max = 40L, quad_n = 1000L)
  expect_lt(coeffs$a0_error, 1e-10)

  t_grid <- seq(0, pi, length.out = 31)
  profile_mu <- distance_profile_beta_mixture2(
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
  profile_minus_mu <- distance_profile_beta_mixture2(
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
    1 - beta_mixture2_cdf_y(
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
    beta_mixture2_cdf_y(
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
  theta <- beta_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.45,
    alpha1 = 3,
    beta1 = 10,
    alpha2 = 12,
    beta2 = 4
  ))
  omega <- jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 25)

  legendre <- distance_profile_beta_mixture2(
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
  integral <- distance_profile_beta_mixture2(
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

test_that("rotational beta-mixture grid and CvM helpers match naive row-wise evaluation", {
  theta <- beta_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
    weight1 = 0.4,
    alpha1 = 2.5,
    beta1 = 8,
    alpha2 = 9,
    beta2 = 2.5
  ))
  omega_grid <- rbind(
    theta$mu,
    -theta$mu,
    jp_normalize_unit_vector(c(0.4, 0.2, 0.89), arg_name = "omega", min_length = 3L)
  )
  t_grid <- seq(0.1, pi - 0.1, length.out = 11)
  x <- r_sph_beta_mixture2(
    n = 6,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )
  dot_products <- pmin(pmax(x %*% t(x), -1), 1)

  for (method in c("legendre", "integral")) {
    tol_match <- 1e-8
    grid_naive <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
      distance_profile_beta_mixture2(
        omega = omega_grid[i, ],
        t_values = t_grid,
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
    }, numeric(length(t_grid))))
    grid_fast <- distance_profile_beta_mixture2_grid(
      omega_grid = omega_grid,
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      distance_type = "geodesic",
      t_grid = t_grid,
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10
    )
    expect_equal(grid_fast, grid_naive, tolerance = tol_match)

    cvm_naive <- t(vapply(seq_len(nrow(x)), function(i) {
      distance_profile_beta_mixture2(
        omega = x[i, ],
        t_values = acos(dot_products[i, ]),
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
    }, numeric(nrow(x))))
    cvm_fast <- distance_profile_beta_mixture2_cvm_grid(
      X = x,
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10
    )
    expect_equal(cvm_fast, cvm_naive, tolerance = tol_match)
  }
})

test_that("rotational beta-mixture sampler and weighted MLE behave coherently", {
  set.seed(20260601)
  theta <- beta_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
    weight1 = 0.4,
    alpha1 = 2.5,
    beta1 = 8,
    alpha2 = 9,
    beta2 = 2.5
  ))
  x <- r_sph_beta_mixture2(
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
  fit_weighted <- beta_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = list(beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
  )
  x_rep <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_rep <- beta_mixture2_mle_s2_weighted(
    x = x_rep,
    control = list(
      beta_mixture2_start_theta = fit_weighted,
      beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(sum(fit_weighted$mu * fit_rep$mu) > 0.95)
  expect_equal(fit_weighted$weight1, fit_rep$weight1, tolerance = 0.08)
  expect_equal(fit_weighted$alpha1, fit_rep$alpha1, tolerance = 1.5)
  expect_equal(fit_weighted$beta1, fit_rep$beta1, tolerance = 2.0)

  y <- (as.numeric(x %*% theta$mu) + 1) / 2
  ecdf_y <- stats::ecdf(y)
  grid <- seq(0.05, 0.95, length.out = 41)
  fitted_cdf <- beta_mixture2_cdf_y(
    y = grid,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )
  expect_lt(max(abs(ecdf_y(grid) - fitted_cdf)), 0.12)
})

test_that("rotational beta-mixture density is numerically stable at y endpoints", {
  log_density <- beta_mixture2_density_y(
    y = c(0, 1),
    weight1 = 0.4,
    alpha1 = 0.4,
    beta1 = 2.5,
    alpha2 = 3.0,
    beta2 = 0.3,
    log = TRUE
  )

  expect_true(all(is.finite(log_density)))

  x <- rbind(c(0, 0, 1), c(0, 0, -1))
  log_s2 <- d_sph_beta_mixture2_s2(
    x = x,
    mu = c(0, 0, 1),
    weight1 = 0.4,
    alpha1 = 0.4,
    beta1 = 2.5,
    alpha2 = 3.0,
    beta2 = 0.3,
    log = TRUE
  )

  expect_true(all(is.finite(log_s2)))
})

test_that("rotational beta-mixture axial and S2 densities are numerically normalized", {
  axial_integral <- integrate(
    f = function(z) beta_mixture2_density_gz(
      z = z,
      weight1 = 0.4,
      alpha1 = 2.5,
      beta1 = 8,
      alpha2 = 9,
      beta2 = 2.5
    ),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  s2_integral <- integrate(
    f = function(z) {
      x <- cbind(sqrt(pmax(0, 1 - z^2)), 0, z)
      2 * pi * d_sph_beta_mixture2_s2(
        x = x,
        mu = c(0, 0, 1),
        weight1 = 0.4,
        alpha1 = 2.5,
        beta1 = 8,
        alpha2 = 9,
        beta2 = 2.5
      )
    },
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  expect_equal(axial_integral, 1, tolerance = 1e-8)
  expect_equal(s2_integral, 1, tolerance = 1e-8)
})

test_that("rotational uniform-beta-mixture coefficients and special profiles are correct", {
  theta <- uniform_beta_mixture_normalize_theta(list(
    mu = c(0, 0, 1),
    weight_uniform = 0.2,
    alpha = 10,
    beta = 3
  ))
  coeffs <- uniform_beta_mixture_legendre_coefficients(theta, l_max = 40L, quad_n = 1000L)
  expect_lt(coeffs$a0_error, 1e-10)

  t_grid <- seq(0, pi, length.out = 31)
  profile_mu <- distance_profile_uniform_beta_mixture(
    omega = theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 120L,
    quad_n = 500L
  )
  profile_minus_mu <- distance_profile_uniform_beta_mixture(
    omega = -theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta,
    distance_type = "geodesic",
    method = "legendre",
    l_max = 120L,
    quad_n = 500L
  )
  thresholds <- cos(t_grid)
  expect_equal(
    profile_mu,
    1 - uniform_beta_mixture_cdf_y(
      y = (thresholds + 1) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    ),
    tolerance = 1e-9
  )
  expect_equal(
    profile_minus_mu,
    uniform_beta_mixture_cdf_y(
      y = (1 - thresholds) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    ),
    tolerance = 1e-9
  )
  expect_equal(profile_mu[[1L]], 0, tolerance = 1e-12)
  expect_equal(profile_mu[[length(profile_mu)]], 1, tolerance = 1e-12)
  expect_true(all(diff(profile_mu) >= -1e-10))
  expect_true(all(profile_mu >= -1e-12 & profile_mu <= 1 + 1e-12))
})

test_that("rotational uniform-beta-mixture Legendre and integral profiles agree", {
  theta <- uniform_beta_mixture_normalize_theta(list(
    mu = c(0, 0, 1),
    weight_uniform = 0.25,
    alpha = 12,
    beta = 4
  ))
  omega <- jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 25)

  legendre <- distance_profile_uniform_beta_mixture(
    omega = omega,
    t_values = t_grid,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta,
    method = "legendre",
    l_max = 150L,
    quad_n = 600L
  )
  integral <- distance_profile_uniform_beta_mixture(
    omega = omega,
    t_values = t_grid,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta,
    method = "integral",
    quad_n = 600L
  )

  expect_lt(max(abs(legendre - integral)), 5e-4)
})

test_that("rotational uniform-beta-mixture grid and CvM helpers match naive row-wise evaluation", {
  theta <- uniform_beta_mixture_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
    weight_uniform = 0.3,
    alpha = 9,
    beta = 2.5
  ))
  omega_grid <- rbind(
    theta$mu,
    -theta$mu,
    jp_normalize_unit_vector(c(0.4, 0.2, 0.89), arg_name = "omega", min_length = 3L)
  )
  t_grid <- seq(0.1, pi - 0.1, length.out = 11)
  x <- r_sph_uniform_beta_mixture(
    n = 6,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
  dot_products <- pmin(pmax(x %*% t(x), -1), 1)

  for (method in c("legendre", "integral")) {
    grid_naive <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
      distance_profile_uniform_beta_mixture(
        omega = omega_grid[i, ],
        t_values = t_grid,
        mu = theta$mu,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
    }, numeric(length(t_grid))))
    grid_fast <- distance_profile_uniform_beta_mixture_grid(
      omega_grid = omega_grid,
      mu = theta$mu,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      t_grid = t_grid,
      distance_type = "geodesic",
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10
    )
    expect_equal(grid_fast, grid_naive, tolerance = 1e-8)

    cvm_naive <- t(vapply(seq_len(nrow(x)), function(i) {
      distance_profile_uniform_beta_mixture(
        omega = x[i, ],
        t_values = acos(dot_products[i, ]),
        mu = theta$mu,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
    }, numeric(nrow(x))))
    cvm_fast <- distance_profile_uniform_beta_mixture_cvm_grid(
      X = x,
      mu = theta$mu,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10
    )
    expect_equal(cvm_fast, cvm_naive, tolerance = 1e-8)
  }
})

test_that("rotational uniform-beta-mixture sampler and weighted MLE behave coherently", {
  set.seed(20260602)
  theta <- uniform_beta_mixture_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
    weight_uniform = 0.25,
    alpha = 10,
    beta = 3
  ))
  x <- r_sph_uniform_beta_mixture(
    n = 180,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
  expect_lt(max(abs(sqrt(rowSums(x^2)) - 1)), 1e-8)

  weights <- rep(c(1, 2, 3), length.out = nrow(x))
  fit_weighted <- uniform_beta_mixture_mle_s2_weighted(
    x = x,
    weights = weights,
    control = list(uniform_beta_mixture_optim_control = list(maxit = 250L, reltol = 1e-9))
  )
  x_rep <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_rep <- uniform_beta_mixture_mle_s2_weighted(
    x = x_rep,
    control = list(
      uniform_beta_mixture_start_theta = fit_weighted,
      uniform_beta_mixture_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(sum(fit_weighted$mu * fit_rep$mu) > 0.95)
  expect_equal(fit_weighted$weight_uniform, fit_rep$weight_uniform, tolerance = 0.08)
  expect_equal(fit_weighted$alpha, fit_rep$alpha, tolerance = 2.0)
  expect_equal(fit_weighted$beta, fit_rep$beta, tolerance = 1.0)

  y <- (as.numeric(x %*% theta$mu) + 1) / 2
  ecdf_y <- stats::ecdf(y)
  grid <- seq(0.05, 0.95, length.out = 41)
  fitted_cdf <- uniform_beta_mixture_cdf_y(
    y = grid,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
  expect_lt(max(abs(ecdf_y(grid) - fitted_cdf)), 0.12)
})

test_that("rotational uniform-beta-mixture density is numerically stable at y endpoints", {
  log_density <- uniform_beta_mixture_density_y(
    y = c(0, 1),
    weight_uniform = 0.2,
    alpha = 0.4,
    beta = 2.5,
    log = TRUE
  )

  expect_true(all(is.finite(log_density)))

  x <- rbind(c(0, 0, 1), c(0, 0, -1))
  log_s2 <- d_sph_uniform_beta_mixture_s2(
    x = x,
    mu = c(0, 0, 1),
    weight_uniform = 0.2,
    alpha = 0.4,
    beta = 2.5,
    log = TRUE
  )

  expect_true(all(is.finite(log_s2)))
})

test_that("rotational uniform-beta-mixture axial and S2 densities are numerically normalized", {
  axial_integral <- integrate(
    f = function(z) uniform_beta_mixture_density_gz(
      z = z,
      weight_uniform = 0.2,
      alpha = 10,
      beta = 3
    ),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  s2_integral <- integrate(
    f = function(z) {
      x <- cbind(sqrt(pmax(0, 1 - z^2)), 0, z)
      2 * pi * d_sph_uniform_beta_mixture_s2(
        x = x,
        mu = c(0, 0, 1),
        weight_uniform = 0.2,
        alpha = 10,
        beta = 3
      )
    },
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  expect_equal(axial_integral, 1, tolerance = 1e-8)
  expect_equal(s2_integral, 1, tolerance = 1e-8)
})

test_that("rotational logit-normal-mixture coefficients and special profiles are correct", {
  theta <- logitnormal_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.45,
    mean1 = -1.1,
    sd1 = 0.45,
    mean2 = 1.0,
    sd2 = 0.55
  ))
  coeffs <- logitnormal_mixture2_legendre_coefficients(theta, l_max = 40L, quad_n = 400L)
  expect_lt(coeffs$a0_error, 1e-10)

  t_grid <- seq(0, pi, length.out = 31)
  profile_mu <- distance_profile_logitnormal_mixture2(
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
  profile_minus_mu <- distance_profile_logitnormal_mixture2(
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
    1 - logitnormal_mixture2_cdf_y(
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
    logitnormal_mixture2_cdf_y(
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
  theta <- logitnormal_mixture2_normalize_theta(list(
    mu = c(0, 0, 1),
    weight1 = 0.35,
    mean1 = -0.6,
    sd1 = 0.4,
    mean2 = 1.3,
    sd2 = 0.5
  ))
  omega <- jp_normalize_unit_vector(c(-0.45, 0.2, 0.87), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 25)

  legendre <- distance_profile_logitnormal_mixture2(
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
  integral <- distance_profile_logitnormal_mixture2(
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

test_that("rotational logit-normal-mixture grid and CvM helpers match naive row-wise evaluation", {
  theta <- logitnormal_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(-0.25, 0.3, 0.92), arg_name = "mu", min_length = 3L),
    weight1 = 0.45,
    mean1 = -1.0,
    sd1 = 0.45,
    mean2 = 1.15,
    sd2 = 0.55
  ))
  omega_grid <- rbind(
    theta$mu,
    -theta$mu,
    jp_normalize_unit_vector(c(-0.35, 0.15, 0.925), arg_name = "omega", min_length = 3L)
  )
  t_grid <- seq(0.1, pi - 0.1, length.out = 11)
  x <- r_sph_logitnormal_mixture2(
    n = 6,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2
  )
  dot_products <- pmin(pmax(x %*% t(x), -1), 1)

  for (method in c("legendre", "integral")) {
    tol_match <- if (identical(method, "legendre")) 1e-8 else 1e-12
    grid_naive <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
      distance_profile_logitnormal_mixture2(
        omega = omega_grid[i, ],
        t_values = t_grid,
        mu = theta$mu,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10,
        eps = 1e-12
      )
    }, numeric(length(t_grid))))
    grid_fast <- distance_profile_logitnormal_mixture2_grid(
      omega_grid = omega_grid,
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      t_grid = t_grid,
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10,
      eps = 1e-12
    )
    expect_equal(grid_fast, grid_naive, tolerance = tol_match)

    cvm_naive <- t(vapply(seq_len(nrow(x)), function(i) {
      distance_profile_logitnormal_mixture2(
        omega = x[i, ],
        t_values = acos(dot_products[i, ]),
        mu = theta$mu,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10,
        eps = 1e-12
      )
    }, numeric(nrow(x))))
    cvm_fast <- distance_profile_logitnormal_mixture2_cvm_grid(
      X = x,
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10,
      eps = 1e-12
    )
    expect_equal(cvm_fast, cvm_naive, tolerance = tol_match)
  }
})

test_that("rotational logit-normal-mixture sampler and weighted MLE behave coherently", {
  set.seed(20260602)
  theta <- logitnormal_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(-0.25, 0.3, 0.92), arg_name = "mu", min_length = 3L),
    weight1 = 0.45,
    mean1 = -1.0,
    sd1 = 0.45,
    mean2 = 1.15,
    sd2 = 0.55
  ))
  x <- r_sph_logitnormal_mixture2(
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
  fit_weighted <- logitnormal_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = list(logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
  )
  x_rep <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_rep <- logitnormal_mixture2_mle_s2_weighted(
    x = x_rep,
    control = list(
      logitnormal_mixture2_start_theta = fit_weighted,
      logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(sum(fit_weighted$mu * fit_rep$mu) > 0.95)
  expect_equal(fit_weighted$weight1, fit_rep$weight1, tolerance = 0.08)
  expect_equal(fit_weighted$mean1, fit_rep$mean1, tolerance = 0.25)
  expect_equal(fit_weighted$mean2, fit_rep$mean2, tolerance = 0.25)
  cdf_grid <- seq(0.1, 0.9, length.out = 21)
  cdf_weighted <- logitnormal_mixture2_cdf_y(
    y = cdf_grid,
    weight1 = fit_weighted$weight1,
    mean1 = fit_weighted$mean1,
    sd1 = fit_weighted$sd1,
    mean2 = fit_weighted$mean2,
    sd2 = fit_weighted$sd2
  )
  cdf_rep <- logitnormal_mixture2_cdf_y(
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
  fitted_cdf <- logitnormal_mixture2_cdf_y(
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
  beta_theta <- beta_mixture2_canonicalize_theta(list(
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

  logit_theta <- logitnormal_mixture2_canonicalize_theta(list(
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

test_that("logit-normal unpack leaves means unrestricted by default", {
  theta <- logitnormal_mixture2_unpack_par(
    par = c(stats::qlogis(0.2), 12, log(0.2), -11, log(0.4), 0, 0, 1),
    control = list(
      logitnormal_mixture2_mean_lower = -8,
      logitnormal_mixture2_mean_upper = 8,
      logitnormal_mixture2_sd_lower = 0.05,
      logitnormal_mixture2_sd_upper = 5,
      logitnormal_mixture2_weight_eps = 0.01
    )
  )

  expect_equal(theta$mean1, -11, tolerance = 1e-12)
  expect_equal(theta$mean2, 12, tolerance = 1e-12)

  clipped <- logitnormal_mixture2_unpack_par(
    par = c(stats::qlogis(0.2), 12, log(0.2), -11, log(0.4), 0, 0, 1),
    control = list(
      logitnormal_mixture2_clip_means = TRUE,
      logitnormal_mixture2_mean_lower = -8,
      logitnormal_mixture2_mean_upper = 8,
      logitnormal_mixture2_sd_lower = 0.05,
      logitnormal_mixture2_sd_upper = 5,
      logitnormal_mixture2_weight_eps = 0.01
    )
  )

  expect_equal(clipped$mean1, -8, tolerance = 1e-12)
  expect_equal(clipped$mean2, 8, tolerance = 1e-12)
})
