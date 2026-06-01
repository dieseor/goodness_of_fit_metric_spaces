library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)
model_specs_path <- if (file.exists(file.path("bootstrap", "model_specs.R"))) {
  file.path("bootstrap", "model_specs.R")
} else {
  file.path("..", "..", "bootstrap", "model_specs.R")
}
source(model_specs_path)

finite_difference_gradient <- function(fn, par, eps = 1e-6) {
  grad <- numeric(length(par))
  for (j in seq_along(par)) {
    step <- rep(0, length(par))
    step[[j]] <- eps
    grad[[j]] <- (fn(par + step) - fn(par - step)) / (2 * eps)
  }
  grad
}

test_that("spherical Cauchy rho zero recovers the uniform S2 profiles", {
  mu <- c(0, 0, 1)
  omega <- c(1, 0, 0)
  t_geo <- c(0.2, 0.7, 1.2)
  t_ch <- c(0.3, 0.9, 1.4)

  observed_geo <- distance_profile_spherical_cauchy(
    omega = omega,
    t_values = t_geo,
    mu = mu,
    rho = 0,
    distance_type = "geodesic",
    warn = FALSE
  )
  observed_ch <- distance_profile_spherical_cauchy(
    omega = omega,
    t_values = t_ch,
    mu = mu,
    rho = 0,
    distance_type = "chordal",
    warn = FALSE
  )

  expect_equal(observed_geo, (1 - cos(t_geo)) / 2, tolerance = 1e-12)
  expect_equal(observed_ch, (t_ch^2) / 4, tolerance = 1e-12)
})

test_that("spherical Cauchy projected CDF matches the axial closed form", {
  mu <- c(0, 0, 1)
  x_grid <- seq(-0.8, 0.8, by = 0.2)

  observed <- spherical_cauchy_projected_cdf(
    x = x_grid,
    omega = mu,
    mu = mu,
    rho = 0.4,
    warn = FALSE
  )
  expected <- spherical_cauchy_axis_projected_cdf(x_grid, rho = 0.4)

  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("spherical Cauchy geodesic profiles are monotone in t", {
  mu <- c(0, 0, 1)
  omega <- c(1, 0, 0)
  t_grid <- seq(0, pi, length.out = 81)

  profile <- distance_profile_spherical_cauchy(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    rho = 0.8,
    distance_type = "geodesic",
    warn = FALSE
  )

  expect_true(all(diff(profile) >= -1e-10))
  expect_equal(profile[[1]], 0, tolerance = 1e-12)
  expect_equal(profile[[length(profile)]], 1, tolerance = 1e-12)
})

test_that("spherical Cauchy sampler returns unit vectors and matches the axial CDF", {
  set.seed(123)
  mu <- c(0, 0, 1)
  sample <- r_sph_spherical_cauchy(n = 2500, mu = mu, rho = 0.35)
  projections <- as.numeric(sample %*% mu)
  x_grid <- seq(-1, 1, length.out = 81)
  empirical_cdf <- vapply(x_grid, function(x0) mean(projections <= x0), numeric(1))
  theoretical_cdf <- spherical_cauchy_axis_projected_cdf(x_grid, rho = 0.35)

  expect_equal(dim(sample), c(2500, 3))
  expect_true(max(abs(sqrt(rowSums(sample^2)) - 1)) < 1e-8)
  expect_lt(max(abs(empirical_cdf - theoretical_cdf)), 0.04)
})

test_that("spherical Cauchy analytic phi-gradient matches finite differences", {
  x <- jp_normalize_unit_matrix(rbind(
    c(0.0, 0.1, 1.0),
    c(0.1, 0.0, 1.0),
    c(-0.1, 0.05, 1.0),
    c(0.05, -0.1, 1.0)
  ))
  prob_weights <- rep(1 / nrow(x), nrow(x))
  phi <- c(0.1, -0.05, 0.2)

  observed <- spherical_cauchy_weighted_loglik_phi_grad(phi, x, prob_weights)
  expected <- finite_difference_gradient(
    fn = function(phi_arg) spherical_cauchy_weighted_loglik_phi(phi_arg, x, prob_weights),
    par = phi,
    eps = 1e-6
  )

  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("spherical Cauchy weighted MLE matches equal-weight and replicated fits", {
  set.seed(2026)
  x <- r_sph_spherical_cauchy(24, mu = c(0, 0, 1), rho = 0.3)

  fit_unweighted <- spherical_cauchy_mle_s2_weighted(
    x,
    weights = NULL,
    control = list(spherical_cauchy_maxit = 300L)
  )
  fit_equal <- spherical_cauchy_mle_s2_weighted(
    x,
    weights = rep(1, nrow(x)),
    control = list(
      spherical_cauchy_maxit = 300L,
      theta_start = fit_unweighted
    )
  )

  weights <- c(2, 1, 3, 2, rep(1, nrow(x) - 4L))
  expanded_x <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_weighted <- spherical_cauchy_mle_s2_weighted(
    x,
    weights = weights,
    control = list(spherical_cauchy_maxit = 300L)
  )
  fit_replicated <- spherical_cauchy_mle_s2_weighted(
    expanded_x,
    weights = NULL,
    control = list(
      spherical_cauchy_maxit = 300L,
      theta_start = fit_weighted
    )
  )

  expect_equal(fit_equal$rho, fit_unweighted$rho, tolerance = 1e-6)
  expect_equal(fit_equal$mu, fit_unweighted$mu, tolerance = 1e-5)
  expect_equal(fit_weighted$rho, fit_replicated$rho, tolerance = 1e-5)
  expect_equal(fit_weighted$mu, fit_replicated$mu, tolerance = 1e-5)
})

test_that("spherical Cauchy grid profile path matches scalar evaluation", {
  mu <- c(0, 0, 1)
  rho <- 0.72
  omega_grid <- generate_canonical_lattice(6, dim = 3)
  t_grid <- seq(0, pi, length.out = 17)

  observed <- distance_profile_spherical_cauchy_grid(
    omega_grid = omega_grid,
    mu = mu,
    rho = rho,
    t_grid = t_grid,
    distance_type = "geodesic",
    warn = FALSE
  )
  expected <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
    distance_profile_spherical_cauchy(
      omega = omega_grid[i, ],
      t_values = t_grid,
      mu = mu,
      rho = rho,
      distance_type = "geodesic",
      warn = FALSE
    )
  }, numeric(length(t_grid))))

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("spherical Cauchy sample profile matrix path matches scalar evaluation", {
  set.seed(2027)
  x <- r_sph_spherical_cauchy(12, mu = c(0, 0, 1), rho = 0.45)
  distance_matrix <- acos(pmin(pmax(x %*% t(x), -1), 1))

  observed <- distance_profile_spherical_cauchy_cvm_grid(
    data = x,
    mu = c(0, 0, 1),
    rho = 0.45,
    distance_matrix = distance_matrix,
    distance_type = "geodesic",
    warn = FALSE
  )
  expected <- t(vapply(seq_len(nrow(x)), function(i) {
    distance_profile_spherical_cauchy(
      omega = x[i, ],
      t_values = distance_matrix[i, ],
      mu = c(0, 0, 1),
      rho = 0.45,
      distance_type = "geodesic",
      warn = FALSE
    )
  }, numeric(nrow(x))))

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("spherical Cauchy model spec fast paths match generic profile matrices", {
  set.seed(2028)
  spec <- make_spherical_cauchy_spec(distance_type = "geodesic")
  theta <- list(mu = c(0, 0, 1), rho = 0.55)
  omega_grid <- generate_canonical_lattice(5, dim = 3)
  t_grid <- seq(0.05, pi - 0.05, length.out = 9)
  data <- r_sph_spherical_cauchy(10, mu = theta$mu, rho = theta$rho)
  distance_matrix <- spec$distance_matrix(data, data)

  fast_grid <- spec$profile_matrix_eval(omega_grid, t_grid, theta, control = list(spherical_cauchy_profile_warn = FALSE))
  generic_grid <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
    spec$profile_eval(omega_grid[i, ], t_grid, theta, control = list(spherical_cauchy_profile_warn = FALSE))
  }, numeric(length(t_grid))))

  fast_sample <- spec$sample_profile_matrix_eval(data, distance_matrix, theta, control = list(spherical_cauchy_profile_warn = FALSE))
  generic_sample <- t(vapply(seq_len(nrow(data)), function(i) {
    spec$profile_eval(data[i, ], distance_matrix[i, ], theta, control = list(spherical_cauchy_profile_warn = FALSE))
  }, numeric(nrow(data))))

  expect_equal(fast_grid, generic_grid, tolerance = 1e-12)
  expect_equal(fast_sample, generic_sample, tolerance = 1e-12)
})
