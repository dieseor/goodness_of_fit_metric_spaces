library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source("utils.R")

test_that("GIG HvMF sampler returns points on H^2", {
  skip_if_not_installed("GIGrvg")

  mu <- c(
    cosh(0.6),
    sinh(0.6) / sqrt(2),
    sinh(0.6) / sqrt(2)
  )

  set.seed(123)
  samples <- rhvmf_h2_gig(25, mu = mu, kappa = 50)

  expect_type(samples, "double")
  expect_equal(dim(samples), c(25, 3))
  expect_true(all(samples[, 1L] > 0))
  expect_true(max(abs(-samples[, 1L]^2 + rowSums(samples[, -1L, drop = FALSE]^2) + 1)) < 1e-8)
})

test_that("GIG HvMF sampler validates its inputs", {
  skip_if_not_installed("GIGrvg")

  expect_error(
    rhvmf_h2_gig(5, mu = c(1, 0, 0), kappa = -1),
    "strictly positive finite scalar"
  )
  expect_error(
    rhvmf_h2_gig(5, mu = c(1, 0, 0, 0), kappa = 1),
    "length 3"
  )
})