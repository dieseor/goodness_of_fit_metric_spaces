library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "model_specs.R"))

test_that("Axial truncated-normal component is numerically normalized in representative regimes", {
  parameter_grid <- list(
    list(kappa = 0.5, mean_value = 0.2),
    list(kappa = 10, mean_value = 0.85),
    list(kappa = 500, mean_value = -0.9),
    list(kappa = 1e6, mean_value = 0.999),
    list(kappa = 1e6, mean_value = -0.999)
  )

  for (params in parameter_grid) {
    integral_value <- integrate(
      f = function(z) {
        exp(axial_truncnorm_component_log_density(
          z = z,
          kappa = params$kappa,
          mean_value = params$mean_value
        ))
      },
      lower = -1,
      upper = 1,
      rel.tol = 1e-8,
      abs.tol = 1e-10,
      subdivisions = 500L
    )$value

    expect_equal(integral_value, 1, tolerance = 1e-6)
  }
})

test_that("Axial truncated-normal CDF remains valid for large concentration", {
  z_grid <- seq(-1, 1, length.out = 41)
  cdf_values <- axial_truncnorm_component_cdf(
    z = z_grid,
    kappa = 1e6,
    mean_value = 0.999
  )

  expect_equal(cdf_values[[1L]], 0, tolerance = 1e-12)
  expect_equal(cdf_values[[length(cdf_values)]], 1, tolerance = 1e-12)
  expect_true(all(diff(cdf_values) >= -1e-12))
  expect_true(all(cdf_values >= 0 & cdf_values <= 1))
})

test_that("Axial truncated-normal decode clips extreme log-kappa values to finite bounds", {
  theta <- axial_truncnorm_decode_theta(c(0, 1000, 0, -1000, 0))

  expect_true(is.finite(theta$kappa1))
  expect_true(is.finite(theta$kappa2))
  expect_equal(theta$kappa1, 1e8, tolerance = 1e-12)
  expect_equal(theta$kappa2, 1e-8, tolerance = 1e-20)
})

test_that("Axial truncated-normal mixture fit respects the hemispheric nu parameterization", {
  set.seed(20260608)
  z <- c(
    pmin(1, pmax(-1, rnorm(60, mean = 0.65, sd = 0.08))),
    pmin(1, pmax(-1, rnorm(40, mean = -0.55, sd = 0.07)))
  )

  theta_start <- list(
    pi = 0.6,
    kappa1 = 20,
    nu1 = 0.6,
    kappa2 = 18,
    nu2 = 0.55
  )

  fit <- fit_axial_truncnorm_mixture2_theta(
    data = z,
    null = list(type = "composite"),
    control = list(
      axial_truncnorm_mixture2_start_theta = theta_start,
      axial_truncnorm_mixture2_optim_control = list(maxit = 100L, reltol = 1e-7)
    )
  )

  expect_true(is.finite(fit$loglik))
  expect_true(fit$pi > 0 && fit$pi < 1)
  expect_true(fit$kappa1 > 0)
  expect_true(fit$kappa2 > 0)
  expect_true(fit$nu1 >= 0 && fit$nu1 < 1)
  expect_true(fit$nu2 >= 0 && fit$nu2 < 1)
})
