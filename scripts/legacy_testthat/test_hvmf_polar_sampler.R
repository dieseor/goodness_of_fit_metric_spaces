library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source("utils.R")

hvmf_polar_test_mu <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))

test_that("polar HvMF sampler returns the requested H2 sample", {
  set.seed(1001)
  x <- rhvmf_h2_polar(25, mu = hvmf_polar_test_mu, kappa = 5)

  expect_equal(dim(x), c(25, 3))
  expect_true(all(x[, 1] > 0))
  expect_lt(max(abs(-x[, 1]^2 + rowSums(x[, -1, drop = FALSE]^2) + 1)), 1e-8)
})

test_that("delta zero reproduces the polar null sampler", {
  set.seed(1002)
  x0 <- rhvmf_h2_polar(12, mu = hvmf_polar_test_mu, kappa = 25)
  set.seed(1002)
  x_delta0 <- rhvmf_h2_angular_mixture(12, mu = hvmf_polar_test_mu, kappa = 25, delta = 0)

  expect_equal(x_delta0, x0, tolerance = 0)
})

test_that("angular mixture preserves the polar radial regularisation", {
  set.seed(1003)
  x0 <- rhvmf_h2_polar(300, mu = hvmf_polar_test_mu, kappa = 25)
  set.seed(1004)
  x1 <- rhvmf_h2_angular_mixture(300, mu = hvmf_polar_test_mu, kappa = 25, delta = 0.1)

  u0 <- acosh(x0[, 1])
  u1 <- acosh(x1[, 1])
  expect_gt(stats::ks.test(u0, u1)$p.value, 1e-4)
})

test_that("angular mixture has the prescribed first conditional moment", {
  set.seed(1005)
  kappa <- 25
  delta <- 0.1
  x <- rhvmf_h2_angular_mixture(1000, mu = hvmf_polar_test_mu, kappa = kappa, delta = delta)
  u <- acosh(x[, 1])
  phi <- atan2(x[, 3], x[, 2])
  theta <- pi / 4
  lambda <- kappa * sinh(acosh(sqrt(2))) * sinh(u)
  a1 <- besselI(lambda, nu = 1, expon.scaled = TRUE) / besselI(lambda, nu = 0, expon.scaled = TRUE)

  observed_cos <- mean(cos(phi - theta))
  expected_cos <- mean(a1 * cos(delta))
  expect_lt(abs(observed_cos - expected_cos), 0.08)
  expect_lt(abs(mean(sin(phi - theta))), 0.08)
})

test_that("polar HvMF samples yield finite MLEs", {
  set.seed(1006)
  x <- rhvmf_h2_polar(100, mu = hvmf_polar_test_mu, kappa = 5)
  fit <- hvmf_mle_h2(x)

  expect_true(all(is.finite(fit$mu)))
  expect_true(is.finite(fit$kappa) && fit$kappa > 0)
  expect_lt(abs(-fit$mu[[1]]^2 + sum(fit$mu[-1]^2) + 1), 1e-8)
})
