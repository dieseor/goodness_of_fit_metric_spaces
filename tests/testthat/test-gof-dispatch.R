library(gofmetric)

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

test_that("gof_test reports which bootstrap was used", {
  fast <- gof_test(x, h0 = "normal", B = 19, seed = 5)
  expect_identical(fast$bootstrap$used, "fast")
  expect_false(fast$bootstrap$fallback)

  slow <- gof_test(x, h0 = "normal", bootstrap = "reestimated", B = 19, seed = 5)
  expect_identical(slow$bootstrap$used, "reestimated")
  expect_false(slow$bootstrap$fallback)

  set.seed(6)
  y <- gofmetric:::r_sph_car(30, c(0, 0, 1), 0.4, 2)
  expect_warning(
    forced <- gof_test(y, h0 = "cardioid", k = 2, B = 9, seed = 7,
                       control = list(cardioid_fast_boundary_eps = 1)),
    "re-estimated bootstrap was used"
  )
  expect_identical(forced$bootstrap$used, "reestimated")
  expect_true(forced$bootstrap$fallback)
  expect_identical(forced$bootstrap$fallback_reason, "cardioid_rho_zero_nonidentification")
})
