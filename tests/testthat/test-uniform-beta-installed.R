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


test_that("rotational uniform-beta-mixture bootstrap supports simple and composite nulls", {
  set.seed(202606031)
  theta <- list(
    mu = c(0, 0, 1),
    weight_uniform = 0.2,
    alpha = 8,
    beta = 2
  )
  x <- r_sph_uniform_beta_mixture(
    n = 12,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(3, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 3)
  )

  result_1 <- multiplier_bootstrap_uniform_beta_mixture(
    data = x,
    null = list(type = "simple", theta = theta),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 2,
    seed = 142,
    n_cores = 1
  )
  result_2 <- multiplier_bootstrap_uniform_beta_mixture(
    data = x,
    null = list(type = "simple", theta = theta),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 2,
    seed = 142,
    n_cores = 1
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)
  expect_true(isTRUE(result_1$diagnostics$weighted_mle))

  composite_result <- multiplier_bootstrap_uniform_beta_mixture(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 1,
    seed = 184,
    n_cores = 1,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      uniform_beta_mixture_n_starts = 1L,
      uniform_beta_mixture_optim_control = list(maxit = 50L, reltol = 1e-7)
    )
  )

  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_length(composite_result$bootstrap$theta_star, 1)
})
