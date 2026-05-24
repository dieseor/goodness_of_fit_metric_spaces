library(testthat)
source(file.path("..", "..", "utils.R"))

hyperboloid_point <- function(chi, theta) {
  c(cosh(chi), sinh(chi) * cos(theta), sinh(chi) * sin(theta))
}

test_that("HvMF H^2 MLE matches the closed-form unweighted formula", {
  data <- rbind(
    hyperboloid_point(0.7, 0.1),
    hyperboloid_point(0.9, 0.5),
    hyperboloid_point(1.1, -0.3),
    hyperboloid_point(0.8, 0.9)
  )

  fit <- hvmf_mle_h2(data)

  S <- colSums(data)
  R <- sqrt(-minkowski_inner_product(S, S))
  xi_hat <- S / R
  kappa_hat <- nrow(data) / (R - nrow(data))

  expect_equal(fit$xi, xi_hat, tolerance = 1e-12)
  expect_equal(fit$mu, xi_hat, tolerance = 1e-12)
  expect_equal(fit$kappa, kappa_hat, tolerance = 1e-12)
  expect_equal(fit$R, R, tolerance = 1e-12)
  expect_equal(fit$W, nrow(data), tolerance = 1e-12)
  expect_equal(fit$xi_inner, -1, tolerance = 1e-10)
  expect_true(fit$xi[1] > 0)
})

test_that("HvMF H^2 weighted MLE matches the closed-form weighted formula", {
  data <- rbind(
    hyperboloid_point(0.6, -0.2),
    hyperboloid_point(0.8, 0.4),
    hyperboloid_point(1.0, 0.7),
    hyperboloid_point(0.9, -0.5)
  )
  weights <- c(0.5, 2, 1.5, 3)

  fit <- hvmf_mle_h2(data, weights = weights)

  W <- sum(weights)
  S_w <- colSums(data * weights)
  R_w <- sqrt(-minkowski_inner_product(S_w, S_w))
  xi_hat_w <- S_w / R_w
  kappa_hat_w <- W / (R_w - W)

  expect_equal(fit$xi, xi_hat_w, tolerance = 1e-12)
  expect_equal(fit$kappa, kappa_hat_w, tolerance = 1e-12)
  expect_equal(fit$R, R_w, tolerance = 1e-12)
  expect_equal(fit$W, W, tolerance = 1e-12)
  expect_equal(fit$xi_inner, -1, tolerance = 1e-10)

  fit_rescaled <- hvmf_mle_h2(data, weights = 10 * weights)
  expect_equal(fit_rescaled$xi, fit$xi, tolerance = 1e-12)
  expect_equal(fit_rescaled$kappa, fit$kappa, tolerance = 1e-12)
})

test_that("HvMF H^2 MLE rejects degenerate and invalid inputs", {
  x <- hyperboloid_point(0.75, 0.3)
  degenerate_data <- rbind(x, x, x)

  expect_error(
    hvmf_mle_h2(degenerate_data),
    "Degenerate or near-degenerate HvMF MLE"
  )

  invalid_data <- rbind(
    hyperboloid_point(0.75, 0.3),
    c(2, 0, 0)
  )
  expect_error(
    hvmf_mle_h2(invalid_data),
    "<x, x>_M = -1"
  )

  valid_data <- rbind(
    hyperboloid_point(0.7, 0.1),
    hyperboloid_point(0.9, 0.5)
  )
  expect_error(
    hvmf_mle_h2(valid_data, weights = c(1, -1)),
    "finite and nonnegative"
  )
})
