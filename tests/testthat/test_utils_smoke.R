library(testthat)
source(file.path("..", "..", "utils.R"))

test_that("generate_canonical_lattice returns unit vectors", {
  mat <- generate_canonical_lattice(20, dim = 3)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 20)
  expect_equal(ncol(mat), 3)
  norms <- sqrt(rowSums(mat^2))
  expect_true(all(abs(norms - 1) < 1e-8))
})

test_that("check_dot_products validates valid/invalid inputs", {
  v <- c(-1, -0.5, 0, 0.5, 1)
  r <- check_dot_products(v)
  expect_equal(v, r)
  expect_error(check_dot_products(c(1.1)))
})

test_that("compute_precomp_vmf returns expected structure when rotasym present", {
  if (!requireNamespace("rotasym", quietly = TRUE)) skip("rotasym not installed")
  set.seed(123)
  mu <- c(0, 0, 1)
  kappa <- 1
  n_mc <- 30
  mc_samples <- matrix(rnorm(n_mc * 3), nrow = n_mc, ncol = 3)
  mc_samples <- t(apply(mc_samples, 1, function(v) v / sqrt(sum(v^2))))
  omega_grid <- generate_canonical_lattice(5, dim = 3)
  t_grid <- c(0.1, 0.5)
  q <- length(mu) - 1
  A_q_kappa <- besselI(kappa, nu = (q + 1) / 2, expon.scaled = TRUE) / besselI(kappa, nu = (q - 1) / 2, expon.scaled = TRUE)
  res <- compute_precomp_vmf(mc_samples, omega_grid, t_grid, mu, kappa, "chordal", A_q_kappa)
  expect_true(is.list(res))
  expect_true(all(c("dists_all", "F2_matrix", "E2_array", "E2_mat", "F2_vec", "m2_mat", "in_ball_list") %in% names(res)))
})
