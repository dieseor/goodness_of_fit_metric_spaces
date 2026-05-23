library(testthat)
library(rotasym)
# When executing from tests/testthat, ensure working directory is repository root so sourced
# files that use 'source("utils.R")' resolve correctly.
oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)
source(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))
source(file.path("utils.R"))

test_that("Edge case handling for distance profiles - chordal/geodesic", {
  mu <- c(0, 0, 1)
  kappa <- 2.0
  omega <- c(0, 0, 1)
  t_values_chordal <- c(1e-5, 0.01, 0.5, 1.0, 1.5, 1.99, 2.0 - 1e-6)
  F_chordal <- theoretical_distance_profile_vmf(omega, mu, kappa, t_values_chordal, "chordal")
  expect_lt(F_chordal[1], 0.01)
  expect_gt(F_chordal[length(F_chordal)], 0.99)

  t_values_geodesic <- c(1e-5, 0.01, 0.5, 1.0, 2.0, 3.0, pi - 1e-6)
  F_geodesic <- theoretical_distance_profile_vmf(omega, mu, kappa, t_values_geodesic, "geodesic")
  expect_lt(F_geodesic[1], 0.01)
  expect_gt(F_geodesic[length(F_geodesic)], 0.99)
})

test_that("Vectorization correctness: row_cov_vmf properties", {
  mu <- c(0, 0, 1)
  kappa <- 2.0
  omega_grid <- generate_canonical_lattice(5)
  t_grid <- seq(0.1, 1.0, length.out = 5)
  mc_samples <- rotasym::r_vMF(2000, mu = mu, kappa = kappa)
  q_sphere <- length(mu) - 1
  A_q_kappa <- besselI(kappa, nu = (q_sphere + 1) / 2, expon.scaled = TRUE) / 
               besselI(kappa, nu = (q_sphere - 1) / 2, expon.scaled = TRUE)
  q_ambient <- length(mu)
  scalar_coef <- 1 - A_q_kappa^2 - ((q_ambient) * A_q_kappa / kappa)
  var_X <- (A_q_kappa / kappa) * diag(q_ambient) + scalar_coef * outer(mu, mu)
  row_vectorized <- row_cov_vmf(h0 = 'simple', unknown_param = NULL, idx = 1, omega_grid = omega_grid, t_grid = t_grid, mu = mu, kappa = kappa, distance_type = "chordal", mc_samples = mc_samples, A_q_kappa = A_q_kappa, var_X = var_X)
  expect_length(row_vectorized, 25)
  expect_true(all(is.finite(row_vectorized)))
})

test_that("Symmetry and diagonal positivity checks", {
  mu <- c(0, 0, 1)
  kappa <- 2.0
  omega_grid <- generate_canonical_lattice(5)
  t_grid <- seq(0.1, 1.0, length.out = 5)
  mc_samples <- rotasym::r_vMF(2000, mu = mu, kappa = kappa)
  q_sphere <- length(mu) - 1
  A_q_kappa <- besselI(kappa, nu = (q_sphere + 1) / 2, expon.scaled = TRUE) / 
               besselI(kappa, nu = (q_sphere - 1) / 2, expon.scaled = TRUE)
  q_ambient <- length(mu)
  scalar_coef <- 1 - A_q_kappa^2 - ((q_ambient) * A_q_kappa / kappa)
  var_X <- (A_q_kappa / kappa) * diag(q_ambient) + scalar_coef * outer(mu, mu)
  rv1 <- row_cov_vmf(h0 = 'simple', unknown_param = NULL, idx = 1, omega_grid = omega_grid, t_grid = t_grid, mu = mu, kappa = kappa, distance_type = "chordal", mc_samples = mc_samples, A_q_kappa = A_q_kappa, var_X = var_X)
  rv2 <- row_cov_vmf(h0 = 'simple', unknown_param = NULL, idx = 2, omega_grid = omega_grid, t_grid = t_grid, mu = mu, kappa = kappa, distance_type = "chordal", mc_samples = mc_samples, A_q_kappa = A_q_kappa, var_X = var_X)
  expect_true(rv1[1] > 0)
  expect_lt(abs(rv1[2] - rv2[1]), 1e-6)
})
