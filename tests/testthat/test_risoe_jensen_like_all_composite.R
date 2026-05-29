library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("wind", "run_risoe_jensen_like_all_composite.R"))

test_that("theta-star diagnostics summarize kappa and Minkowski errors", {
  theta_star <- list(
    list(mu = c(1, 0, 0), kappa = 1.2),
    list(mu = c(cosh(0.3), sinh(0.3), 0), kappa = 1.5)
  )

  diag <- extract_theta_star_diagnostics(theta_star)

  expect_equal(diag$kappa_values, c(1.2, 1.5))
  expect_true(all(diag$mu_norm_errors < 1e-12))
  expect_identical(diag$n_failed, 0L)
})

test_that("observed CvM mismatch guard tolerates machine noise and rejects large gaps", {
  expect_no_error(
    stop_if_cvm_observed_mismatch(
      dataset_id = "demo",
      cvm_simple = 0.123456789,
      cvm_composite = 0.123456789 + 1e-12,
      tol = 1e-8
    )
  )

  expect_error(
    stop_if_cvm_observed_mismatch(
      dataset_id = "demo",
      cvm_simple = 0.12,
      cvm_composite = 0.13,
      tol = 1e-8
    ),
    "Observed CvM mismatch"
  )
})
