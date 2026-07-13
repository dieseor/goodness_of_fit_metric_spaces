project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
old_wd <- getwd()
setwd(project_root)
on.exit(setwd(old_wd), add = TRUE)
source("convergence_empirical_process/gaussian_process_normal.R")
source("convergence_empirical_process/gaussian_process_vmf.R")

test_that("normal sample KS helper matches brute-force simple null", {
  sample_x <- c(-1.1, -0.3, 0.2, 1.4)
  mu <- 0.1
  sigma <- 1.3

  brute_force <- function() {
    n <- length(sample_x)
    max_diff <- 0
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        radius_ij <- abs(sample_x[[j]] - sample_x[[i]])
        empirical_ij <- mean(abs(sample_x - sample_x[[i]]) <= radius_ij)
        theoretical_ij <- theoretical_distance_profile_normal(sample_x[[i]], mu, sigma, radius_ij)
        max_diff <- max(max_diff, abs(empirical_ij - theoretical_ij))
      }
    }
    sqrt(n) * max_diff
  }

  expect_equal(
    compute_sample_ks_sup_normal(sample_x, mu = mu, sigma = sigma, h0 = "simple"),
    brute_force(),
    tolerance = 1e-12
  )
})

test_that("normal sample KS helper matches brute-force composite null", {
  sample_x <- c(-1.1, -0.3, 0.2, 1.4)

  brute_force <- function() {
    n <- length(sample_x)
    mu_hat <- mean(sample_x)
    sigma_hat <- sqrt(mean((sample_x - mu_hat)^2))
    max_diff <- 0
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        radius_ij <- abs(sample_x[[j]] - sample_x[[i]])
        empirical_ij <- mean(abs(sample_x - sample_x[[i]]) <= radius_ij)
        theoretical_ij <- theoretical_distance_profile_normal(sample_x[[i]], mu_hat, sigma_hat, radius_ij)
        max_diff <- max(max_diff, abs(empirical_ij - theoretical_ij))
      }
    }
    sqrt(n) * max_diff
  }

  expect_equal(
    compute_sample_ks_sup_normal(sample_x, mu = 0, sigma = 1, h0 = "composite", unknown_param = "both"),
    brute_force(),
    tolerance = 1e-12
  )
})

test_that("vmf sample KS helper matches brute-force simple null", {
  sample_data <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 1, 1) / sqrt(3)
  )
  mu <- c(0, 0, 1)
  kappa <- 1.2

  brute_force <- function() {
    n <- nrow(sample_data)
    dot_products <- sample_data %*% t(sample_data)
    dot_products <- pmax(pmin(dot_products, 1), -1)
    distance_matrix <- acos(dot_products)
    max_diff <- 0
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        radius_ij <- distance_matrix[j, i]
        empirical_ij <- mean(distance_matrix[, i] <= radius_ij)
        theoretical_ij <- theoretical_distance_profile_vmf(sample_data[i, ], mu, kappa, radius_ij, "geodesic")
        max_diff <- max(max_diff, abs(empirical_ij - theoretical_ij))
      }
    }
    sqrt(n) * max_diff
  }

  expect_equal(
    compute_sample_ks_sup_vmf(sample_data, mu = mu, kappa = kappa, distance_type = "geodesic", h0 = "simple"),
    brute_force(),
    tolerance = 1e-10
  )
})

test_that("empirical process simulators run in sample mode", {
  normal_sup <- simulate_empirical_process_normal(
    omega_grid = seq(-1, 1, length.out = 5),
    t_grid = seq(0, 2, length.out = 5),
    n = 8,
    mu = 0,
    sigma = 1,
    M = 2,
    n_cores = 1,
    h0 = "simple",
    empirical_ks_mode = "sample"
  )
  expect_length(normal_sup, 2)
  expect_true(all(is.finite(normal_sup)))

  vmf_sup <- simulate_empirical_process_vmf(
    omega_grid = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)),
    t_grid = seq(0.2, 2, length.out = 4),
    n = 6,
    mu = c(0, 0, 1),
    kappa = 1,
    distance_type = "geodesic",
    M = 2,
    n_cores = 1,
    seed = 123,
    h0 = "simple",
    empirical_ks_mode = "sample"
  )
  expect_length(vmf_sup, 2)
  expect_true(all(is.finite(vmf_sup)))
})
