library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "restricted_spiked_normal_openmx.R"))
source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))

test_that("restricted-spiked covariance has the prescribed spectrum", {
  theta <- c(1, -2, 0.5, 0.75, -1.25)
  lambda <- 3.5
  fitted <- normalize_restricted_spiked_normal_theta(list(theta = theta, lambda = lambda))
  eig <- eigen(fitted$Sigma, symmetric = TRUE)
  expect_equal(sort(eig$values), sort(c(1 + lambda, rep(1, 4L))), tolerance = 1e-12)
  expect_equal(drop(fitted$Sigma %*% fitted$u), (1 + lambda) * fitted$u, tolerance = 1e-12)
  expect_equal(fitted$Sigma, t(fitted$Sigma), tolerance = 1e-14)
  expect_true(all(eig$values > 0))
})

test_that("restricted-spiked simulation has the specified moments", {
  set.seed(4101)
  theta <- c(0.8, -0.4)
  lambda <- 1.75
  x <- rrestricted_spiked_normal(50000L, theta, lambda)
  target <- restricted_spiked_normal_covariance(theta, lambda)
  expect_equal(colMeans(x), theta, tolerance = 0.025)
  expect_equal(stats::cov(x) * (nrow(x) - 1) / nrow(x), target, tolerance = 0.04)
})

test_that("restricted-spiked custom MLE agrees with weighted replicated data", {
  set.seed(4102)
  x <- rrestricted_spiked_normal(20L, c(0.75, -0.3), 1.2)
  weights <- c(1, 2, 3, 1, 2, 1, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3, 1, 2, 3, 1)
  weighted <- fit_restricted_spiked_normal_theta(
    x, weights = weights, null = list(type = "composite")
  )
  replicated <- fit_restricted_spiked_normal_theta(
    x[rep(seq_len(nrow(x)), weights), , drop = FALSE],
    null = list(type = "composite")
  )
  expect_equal(weighted$theta, replicated$theta, tolerance = 1e-8)
  expect_equal(weighted$lambda, replicated$lambda, tolerance = 1e-7)
  expect_equal(weighted$Sigma, replicated$Sigma, tolerance = 1e-7)
})

test_that("restricted-spiked score agrees with numerical log-likelihood derivatives", {
  set.seed(4103)
  x <- rrestricted_spiked_normal(25L, c(0.7, -0.5), 1.3)
  parameter <- c(theta = c(0.6, -0.4), lambda = 0.9)
  score_mean <- colMeans(restricted_spiked_normal_score_matrix(
    x, list(theta = parameter[1:2], lambda = parameter[[3L]])
  ))
  step <- 1e-6
  numeric_gradient <- vapply(seq_along(parameter), function(j) {
    left <- parameter; right <- parameter
    left[[j]] <- left[[j]] - step
    right[[j]] <- right[[j]] + step
    (restricted_spiked_normal_loglik(x, list(theta = right[1:2], lambda = right[[3L]])) -
       restricted_spiked_normal_loglik(x, list(theta = left[1:2], lambda = left[[3L]]))) / (2 * step)
  }, numeric(1L))
  expect_equal(unname(score_mean), unname(numeric_gradient), tolerance = 3e-6)
})

test_that("restricted-spiked custom MLE attains the grid-checked profile maximum", {
  set.seed(4104)
  x <- rrestricted_spiked_normal(80L, c(0.8, -0.25, 0.5), 2)
  fitted <- fit_restricted_spiked_normal_theta(
    x, null = list(type = "composite"),
    control = list(restricted_spiked_profile_grid_size = 801L)
  )
  xbar <- colMeans(x)
  S <- crossprod(sweep(x, 2L, xbar, FUN = "-")) / nrow(x)
  diagnostic <- restricted_spiked_normal_profile_maximize(xbar, S, grid_size = 10001L)
  expect_equal(fitted$fit_diagnostics$profile_value, diagnostic$best$value, tolerance = 2e-7)
  expect_equal(fitted$lambda, diagnostic$best$lambda, tolerance = 2e-5)
})

test_that("restricted-spiked MLE stops at a numerically null profiled radius", {
  x <- rbind(c(-1, 0), c(1, 0), c(0, -2), c(0, 2))
  expect_error(
    fit_restricted_spiked_normal_theta(x, null = list(type = "composite")),
    "profiled radius"
  )
})

test_that("restricted-spiked model excludes the lambda boundary", {
  expect_error(
    normalize_restricted_spiked_normal_theta(list(theta = c(0.8, -0.2), lambda = 0)),
    "strictly positive"
  )
  x <- rbind(c(1.05, 0), c(0.55, 0), c(0.8, 0.25), c(0.8, -0.25))
  expect_error(
    fit_restricted_spiked_normal_theta(x, null = list(type = "composite")),
    "excluded boundary"
  )
})

test_that("restricted-spiked likelihood equals mvtnorm density", {
  set.seed(4105)
  theta <- c(0.5, -0.6, 0.3)
  lambda <- 1.1
  x <- matrix(rnorm(30L), ncol = 3L)
  fitted <- list(theta = theta, lambda = lambda)
  expected <- mean(mvtnorm::dmvnorm(
    x, mean = theta, sigma = restricted_spiked_normal_covariance(theta, lambda), log = TRUE
  ))
  expect_equal(restricted_spiked_normal_loglik(x, fitted), expected, tolerance = 1e-12)
})

test_that("restricted-spiked model is orthogonally equivariant", {
  set.seed(4106)
  theta <- c(0.9, -0.5, 0.4)
  lambda <- 1.4
  x <- rrestricted_spiked_normal(120L, theta, lambda)
  Q <- qr.Q(qr(matrix(rnorm(9L), 3L, 3L)))
  fit <- fit_restricted_spiked_normal_theta(x, null = list(type = "composite"))
  rotated_fit <- fit_restricted_spiked_normal_theta(x %*% Q, null = list(type = "composite"))
  expect_equal(rotated_fit$theta, drop(t(Q) %*% fit$theta), tolerance = 1e-6)
  expect_equal(rotated_fit$lambda, fit$lambda, tolerance = 1e-6)
  expect_equal(rotated_fit$Sigma, t(Q) %*% fit$Sigma %*% Q, tolerance = 1e-6)
})

test_that("restricted-spiked custom MLE agrees with the independent OpenMx reference", {
  skip_if_not_installed("OpenMx")
  set.seed(4107)
  x <- rrestricted_spiked_normal(100L, c(0.85, -0.35), 1.25)
  custom <- fit_restricted_spiked_normal_theta(x, null = list(type = "composite"))
  reference <- fit_restricted_spiked_normal_openmx(x)
  expect_equal(custom$loglik * nrow(x), reference$loglik, tolerance = 5e-5)
  expect_equal(custom$loglik * nrow(x), reference$loglik_openmx, tolerance = 5e-5)
  expect_equal(custom$theta, reference$theta, tolerance = 2e-4)
  expect_equal(custom$lambda, reference$lambda, tolerance = 2e-4)
  expect_equal(custom$Sigma, reference$Sigma, tolerance = 2e-4)
})

test_that("restricted-spiked analytical score matches numerical gradients across dimensions", {
  configurations <- list(
    list(theta = c(0.7, -0.4), lambda = 1e-4),
    list(theta = c(0.8, 0.2), lambda = 0.35),
    list(theta = c(0.6, -0.3, 0.2, 0.1, -0.4), lambda = 2.5)
  )
  step <- 1e-6
  for (configuration in configurations) {
    d <- length(configuration$theta)
    points <- rrestricted_spiked_normal(
      4L, configuration$theta, configuration$lambda
    )
    for (i in seq_len(nrow(points))) {
      point <- points[i, , drop = FALSE]
      parameter <- c(configuration$theta, configuration$lambda)
      numeric_gradient <- vapply(seq_along(parameter), function(j) {
        left <- parameter; right <- parameter
        left[[j]] <- left[[j]] - step
        right[[j]] <- right[[j]] + step
        (
          restricted_spiked_normal_loglik(point, list(theta = right[seq_len(d)], lambda = right[[d + 1L]])) -
            restricted_spiked_normal_loglik(point, list(theta = left[seq_len(d)], lambda = left[[d + 1L]]))
        ) / (2 * step)
      }, numeric(1L))
      score <- restricted_spiked_normal_score_matrix(
        point, list(theta = configuration$theta, lambda = configuration$lambda)
      )
      expect_equal(as.numeric(score), numeric_gradient, tolerance = 5e-6)
    }
  }
})

test_that("restricted-spiked V is minus the analytical Fisher information", {
  theta <- list(theta = c(0.7, -0.2, 0.4), lambda = 1.3)
  fisher <- restricted_spiked_normal_fisher_information(theta)
  V <- restricted_spiked_normal_score_jacobian_V(theta)
  expect_equal(V, -fisher, tolerance = 1e-14)
  expect_equal(fisher, t(fisher), tolerance = 1e-14)
  expect_true(all(eigen(fisher, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(fisher[1:3, 4L], rep(0, 3L), tolerance = 1e-14)
  expect_equal(fisher[4L, 1:3], rep(0, 3L), tolerance = 1e-14)
})

test_that("restricted-spiked fast preparation uses its own score and Fisher information", {
  set.seed(4108)
  # This test requires an interior fitted spike because the model now rejects
  # the non-identified lambda = 0 boundary explicitly.
  x <- rrestricted_spiked_normal(80L, c(0.8, -0.3), 2)
  spec <- make_restricted_spiked_normal_spec()
  theta_hat <- spec$fit_theta(x, null = list(type = "composite"))
  ks_prep <- prepare_ks_observed_data(
    x, spec, theta_hat, make_sample_unique_distance_ks_grid()
  )
  cvm_prep <- prepare_cvm_observed_data(x, spec, theta_hat, control = list())
  prep <- spec_fast_multiplier_prepare(
    spec, x, theta_hat, ks_prep, cvm_prep,
    control = list(derivative_method = "score_mc", derivative_mc_size = 400L,
                   derivative_mc_seed = 4109L, fast_multiplier_cvm_block_size = 10L)
  )
  expected_score <- restricted_spiked_normal_score_matrix(x, theta_hat)
  expected_fisher <- restricted_spiked_normal_fisher_information(theta_hat)
  expect_equal(prep$S_obs, expected_score, tolerance = 1e-12)
  expect_equal(prep$Vhat, expected_fisher, tolerance = 1e-12)
  expect_equal(prep$paper_Vhat, -expected_fisher, tolerance = 1e-12)
  expect_identical(prep$vhat_method, "restricted_spiked_analytic_fisher")
  expect_identical(prep$auxiliary_sampler, "restricted_spiked_normal")
  expect_identical(prep$score_parameterization, "beta=(theta,lambda)")
  expect_equal(dim(prep$S_obs), c(nrow(x), ncol(x) + 1L))
  expect_equal(dim(prep$Psi_aux), c(400L, ncol(x) + 1L))
})
