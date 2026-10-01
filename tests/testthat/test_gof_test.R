library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("utils.R"))

set.seed(41)
x <- rnorm(30)

test_that("gof_test gives the same result as the model wrapper", {
  fast <- gof_test(x, h0 = "normal", B = 49, seed = 7)
  direct <- multiplier_bootstrap_normal(
    x, null = list(type = "composite"), unknown_param = "both", B = 49, seed = 7,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE)
  )
  expect_identical(fast$inference, direct$inference)

  slow <- gof_test(x, h0 = "normal", bootstrap = "reestimated", B = 19, seed = 7)
  direct <- multiplier_bootstrap_normal(
    x, null = list(type = "composite"), unknown_param = "both", B = 19, seed = 7,
    bootstrap_method = "reestimated"
  )
  expect_identical(slow$inference, direct$inference)
})

test_that("gof_test uses theta for a simple null and passes model arguments on", {
  theta <- list(mu = 0, sigma = 1)
  simple <- gof_test(x, h0 = "normal", theta = theta, bootstrap = "reestimated", B = 19, seed = 3)
  direct <- multiplier_bootstrap_normal(
    x, null = list(type = "simple", theta = theta), B = 19, seed = 3, bootstrap_method = "reestimated"
  )
  expect_identical(simple$inference, direct$inference)
  mu_only <- gof_test(x, h0 = "normal", unknown_param = "mu",
                      null = list(type = "composite", fixed = list(sigma = 1)), B = 19, seed = 3)
  expect_true(is.finite(mu_only$inference$ks$p_value))
})

test_that("gof_test rejects unknown null families", {
  expect_error(gof_test(x, h0 = "not_a_model"), "must be one of")
  expect_error(gof_test(x), "must be one of")
})
