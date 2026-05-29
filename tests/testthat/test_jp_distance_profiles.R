library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)

test_that("c_jp_sphere at alpha zero recovers the uniform constant", {
  expect_equal(
    c_jp_sphere(q = 2, alpha = 0, beta = 2),
    1 / sphere_surface_area(2),
    tolerance = 1e-14
  )
  expect_equal(
    c_jp_sphere(q = 3, alpha = 0, beta = -1.5),
    1 / sphere_surface_area(3),
    tolerance = 1e-14
  )
})

test_that("d_proj_jp integrates to one", {
  mu2 <- c(0, 0, 1)
  mu3 <- c(0, 0, 0, 1)

  integral_q2 <- integrate(
    f = function(t) d_proj_jp(t, mu = mu2, kappa = 2, psi = 0.5),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value
  integral_q3 <- integrate(
    f = function(t) d_proj_jp(t, mu = mu3, kappa = 1.5, psi = -0.5),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  expect_equal(integral_q2, 1, tolerance = 1e-6)
  expect_equal(integral_q3, 1, tolerance = 1e-6)
})

test_that("p_proj_jp integrates to one", {
  mu2 <- c(0, 0, 1)
  omega2 <- c(1, 0, 0)
  mu3 <- c(0, 0, 0, 1)
  omega3 <- c(1, 0, 0, 0)

  integral_q2 <- integrate(
    f = function(s) p_proj_jp(s, omega = omega2, mu = mu2, kappa = 2, psi = 0.5),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value
  integral_q3 <- integrate(
    f = function(s) p_proj_jp(s, omega = omega3, mu = mu3, kappa = 1.5, psi = -0.5),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  expect_equal(integral_q2, 1, tolerance = 2e-3)
  expect_equal(integral_q3, 1, tolerance = 2e-3)
})

test_that("p_proj_jp with omega equal to mu matches d_proj_jp", {
  mu <- c(0, 0, 1)
  s_grid <- seq(-0.9, 0.9, by = 0.2)

  expect_equal(
    p_proj_jp(s_grid, omega = mu, mu = mu, kappa = 2, psi = 0.5),
    d_proj_jp(s_grid, mu = mu, kappa = 2, psi = 0.5),
    tolerance = 1e-8
  )
})

test_that("uniform S2 geodesic distance profile is recovered exactly", {
  mu <- c(0, 0, 1)
  omega <- c(1, 0, 0)
  t_values <- c(0.1, 0.7, 1.2, 2.5)

  observed <- distance_profile_jp(
    omega = omega,
    t_values = t_values,
    mu = mu,
    kappa = 0,
    psi = 1,
    distance_type = "geodesic"
  )
  expected <- (1 - cos(t_values)) / 2

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("psi zero delegates exactly to the existing vMF profile", {
  if (!requireNamespace("rotasym", quietly = TRUE)) {
    skip("rotasym not installed")
  }

  mu <- c(0, 0, 1)
  omega <- c(1, 0, 0)
  t_values <- c(0.2, 0.8, 1.4)

  observed <- distance_profile_jp(
    omega = omega,
    t_values = t_values,
    mu = mu,
    kappa = 2,
    psi = 0,
    distance_type = "geodesic"
  )
  expected <- theoretical_distance_profile_vmf_s2_fast(
    omega = omega,
    mu = mu,
    kappa = 2,
    t_values = t_values,
    distance_type = "geodesic"
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("r_sph_jp returns unit vectors and matches the inversion CDF", {
  set.seed(123)
  mu <- c(0, 0, 1)
  sample <- r_sph_jp(n = 2000, mu = mu, kappa = 2, psi = 0.5)

  expect_equal(dim(sample), c(2000, 3))
  expect_true(max(abs(sqrt(rowSums(sample^2)) - 1)) < 1e-8)

  cdf_table <- build_jp_axis_cdf_table(mu = mu, kappa = 2, psi = 0.5, grid_size = 8193L)
  projections <- as.numeric(sample %*% mu)
  t_grid <- seq(-1, 1, length.out = 81)
  empirical_cdf <- vapply(t_grid, function(t0) mean(projections <= t0), numeric(1))
  theoretical_cdf <- jp_interpolate_cdf(t_grid, x_grid = cdf_table$u, cdf_grid = cdf_table$cdf)

  expect_lt(max(abs(empirical_cdf - theoretical_cdf)), 0.04)
})

test_that("alpha zero sampler recovers uniformity on S2", {
  set.seed(456)
  mu <- c(0, 0, 1)
  sample <- r_sph_jp(n = 2000, mu = mu, kappa = 0, psi = 1)
  projections <- as.numeric(sample %*% mu)
  t_grid <- seq(-1, 1, length.out = 81)
  empirical_cdf <- vapply(t_grid, function(t0) mean(projections <= t0), numeric(1))
  theoretical_cdf <- (t_grid + 1) / 2

  expect_true(max(abs(sqrt(rowSums(sample^2)) - 1)) < 1e-8)
  expect_lt(max(abs(empirical_cdf - theoretical_cdf)), 0.04)
})
