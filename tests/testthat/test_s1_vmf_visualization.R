library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("utils.R"))
source(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))
source(file.path("convergence_empirical_process", "gaussian_process_s1_visualization.R"))

test_that("generate_circle_grid returns deterministic unit vectors in angular order", {
  grid <- generate_circle_grid(6)
  expect_true(is.data.frame(grid))
  expect_equal(nrow(grid), 6)
  expect_equal(colnames(grid), c("theta", "x", "y", "label"))
  expect_true(all(diff(grid$theta) > 0))
  expect_equal(grid$theta[1], 0, tolerance = 1e-12)
  norms <- sqrt(grid$x^2 + grid$y^2)
  expect_equal(norms, rep(1, 6), tolerance = 1e-10)
})

test_that("S1 chordal distance profile is monotone with correct endpoints", {
  mu <- c(1, 0)
  omega <- c(1, 0)
  t_values <- seq(0, 2, length.out = 9)
  profile <- theoretical_distance_profile_vmf_s1_chordal(
    omega = omega,
    mu = mu,
    kappa = 2,
    t_values = t_values,
    cdf_grid_size = 4097
  )
  expect_equal(profile[1], 0, tolerance = 1e-12)
  expect_equal(profile[length(profile)], 1, tolerance = 1e-12)
  expect_true(all(diff(profile) >= -1e-8))
})

test_that("S1 chordal inversion is consistent with the deterministic CDF", {
  mu <- c(1, 0)
  omega <- c(cos(pi / 4), sin(pi / 4))
  cdf_object <- build_vmf_s1_cdf(mu, kappa = 2, n_grid = 4097)
  u_values <- seq(0.1, 0.9, by = 0.2)
  t_values <- invert_distance_profile_vmf_s1_chordal(
    omega = omega,
    mu = mu,
    kappa = 2,
    u_values = u_values,
    cdf_object = cdf_object,
    tol = 1e-8
  )
  recovered <- theoretical_distance_profile_vmf_s1_chordal(
    omega = omega,
    mu = mu,
    kappa = 2,
    t_values = t_values,
    cdf_object = cdf_object
  )
  expect_equal(recovered, u_values, tolerance = 1e-4)
})

test_that("Exact S1 joint probabilities are symmetric", {
  mu <- c(1, 0)
  omega1 <- c(1, 0)
  omega2 <- c(cos(2 * pi / 3), sin(2 * pi / 3))
  cdf_object <- build_vmf_s1_cdf(mu, kappa = 2, n_grid = 4097)
  p12 <- joint_probability_vmf_s1_chordal_exact(
    omega1 = omega1,
    t1 = 0.75,
    omega2 = omega2,
    t2 = 1.10,
    mu = mu,
    kappa = 2,
    cdf_object = cdf_object
  )
  p21 <- joint_probability_vmf_s1_chordal_exact(
    omega1 = omega2,
    t1 = 1.10,
    omega2 = omega1,
    t2 = 0.75,
    mu = mu,
    kappa = 2,
    cdf_object = cdf_object
  )
  expect_equal(p12, p21, tolerance = 1e-10)
})

test_that("Exact S1 covariance route returns a finite valid matrix on a small grid", {
  mu <- c(1, 0)
  omega_grid <- as.matrix(generate_circle_grid(4)[, c("x", "y")])
  t_grid <- c(0.3, 0.9, 1.5)
  sigma_exact <- cov_vmf(
    omega_grid = omega_grid,
    t_grid = t_grid,
    mu = mu,
    kappa = 2,
    distance_type = "chordal",
    h0 = "simple",
    cov_method = "exact_s1_simple",
    n_cores = 1,
    cdf_grid_size = 4097
  )
  expect_true(all(is.finite(sigma_exact)))
  diag_exact <- validate_covariance_matrix(sigma_exact, stop_on_failure = FALSE)
  expect_true(isTRUE(diag_exact$valid))
})

test_that("S1 visualization works in the exact covariance mode", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("mvtnorm")

  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 4,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 30,
    cdf_grid_size = 4097
  )
  expect_true(!is.null(res$plot))
  expect_true(is.data.frame(res$curve_data))
  expect_true(nrow(res$curve_data) > 0)
})

test_that("S1 visualization works in the Monte Carlo covariance mode on a stable small grid", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("mvtnorm")

  omega_grid <- as.matrix(generate_circle_grid(2)[, c("x", "y")])
  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    omega_grid = omega_grid,
    t_grid = c(0.5, 1.2),
    cov_method = "mc",
    n_mc_samples = 20000,
    n_cores = 1,
    seed = 123,
    curve_points = 25,
    cdf_grid_size = 4097
  )
  expect_true(!is.null(res$plot))
  expect_true(is.data.frame(res$curve_data))
  expect_true(nrow(res$curve_data) > 0)
})
