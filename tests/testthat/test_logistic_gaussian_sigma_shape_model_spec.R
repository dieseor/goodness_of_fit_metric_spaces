library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "logistic_gaussian_sigma_shape_bootstrap.R"))
source(file.path("scripts", "run_logistic_gaussian_sigma_shape_scenarios.R"))

test_that("Maphosa M2 covariance shape is correct in repository ilr coordinates", {
  for (d in c(2L, 5L)) {
    A <- logistic_gaussian_maphosa_shape(d)
    expect_equal(diag(A), c(rep(1, d - 1L), 1 / (d + 1)), tolerance = 0)
    expect_equal(A %*% logistic_gaussian_maphosa_shape_inverse(d), diag(d), tolerance = 1e-14)
  }
})

test_that("restricted Logistic-Gaussian weighted MLE has the closed form", {
  set.seed(101)
  d <- 2L
  x <- rlogistic_gaussian_sigma_shape(
    n = 20,
    mu_ilr = c(0.3, -0.4),
    sigma = 0.8
  )
  weights <- seq_len(nrow(x))
  probability_weights <- weights / sum(weights)
  z <- logistic_gaussian_ilr_matrix(x)
  mu_expected <- colSums(z * probability_weights)
  centered <- sweep(z, 2L, mu_expected, FUN = "-")
  Ainv <- logistic_gaussian_maphosa_shape_inverse(d)
  sigma2_expected <- sum(
    probability_weights * rowSums((centered %*% Ainv) * centered)
  ) / d

  fit <- fit_logistic_gaussian_sigma_shape_theta(
    data = x,
    weights = weights,
    null = list(type = "composite")
  )

  expect_equal(fit$mu_ilr, mu_expected, tolerance = 1e-12)
  expect_equal(fit$sigma^2, sigma2_expected, tolerance = 1e-12)
  expect_equal(fit$Sigma_ilr, sigma2_expected * logistic_gaussian_maphosa_shape(d), tolerance = 1e-12)
})

test_that("restricted score vanishes at the unweighted MLE", {
  set.seed(202)
  x <- rlogistic_gaussian_sigma_shape(
    n = 30,
    mu_ilr = c(-0.2, 0.4),
    sigma = 1.1
  )
  fit <- fit_logistic_gaussian_sigma_shape_theta(
    x,
    null = list(type = "composite")
  )
  score <- logistic_gaussian_sigma_shape_score_matrix(x, fit)
  expect_lt(max(abs(colMeans(score))), 1e-10)
})

test_that("Dirichlet M5 alpha has the intended pattern and concentration", {
  for (d in c(2L, 5L)) {
    alpha <- maphosa_dirichlet_alpha(d)
    D <- d + 1L
    expect_equal(alpha, 10 * seq_len(D) / (D * (D + 1L)), tolerance = 0)
    expect_equal(sum(alpha), 5, tolerance = 1e-14)

    alpha_twice <- maphosa_dirichlet_alpha(d, concentration_multiplier = 2)
    expect_equal(alpha_twice / sum(alpha_twice), alpha / sum(alpha), tolerance = 1e-14)
    expect_equal(sum(alpha_twice), 10, tolerance = 1e-14)
  }
})

test_that("t-logistic generator produces valid simplex observations", {
  skip_if_not_installed("mvtnorm")
  set.seed(303)
  x <- r_maphosa_t_logistic(n = 10, d = 5, nu = 4, standardized = FALSE)
  expect_equal(dim(x), c(10L, 6L))
  expect_true(all(x > 0))
  expect_equal(rowSums(x), rep(1, 10), tolerance = 1e-12)
})

test_that("Dirichlet generator uses gtools and produces valid simplex observations", {
  skip_if_not_installed("gtools")
  set.seed(404)
  x <- r_maphosa_dirichlet(n = 10, d = 2)
  expect_equal(dim(x), c(10L, 3L))
  expect_true(all(x > 0))
  expect_equal(rowSums(x), rep(1, 10), tolerance = 1e-12)
})

test_that("restricted Logistic-Gaussian bootstrap runs in smoke size", {
  set.seed(505)
  x <- rlogistic_gaussian_sigma_shape(
    n = 20,
    mu_ilr = c(0, 0),
    sigma = 1
  )
  result <- multiplier_bootstrap_logistic_gaussian_sigma_shape(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 5,
    seed = 123,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 100,
      derivative_mc_seed = 321,
      fast_multiplier_backend = "r"
    ),
    fast_multiplier_backend = "r",
    fuse_ks_cvm = TRUE
  )

  expect_true(result$inference$ks$p_value >= 0 && result$inference$ks$p_value <= 1)
  expect_true(result$inference$cvm$p_value >= 0 && result$inference$cvm$p_value <= 1)
  expect_equal(length(result$observed$theta_hat$mu_ilr), 2L)
  expect_true(is.finite(result$observed$theta_hat$sigma))
})
