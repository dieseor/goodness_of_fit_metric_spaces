library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path(
  "real_data",
  "sunspots",
  "run_sunspots_cycle23_temporal_uniform_beta_gof.R"
))

temporal_uniform_beta_test_eta <- function() {
  sunspots_joint_parsimonious_time_validate_eta(
    list(
      time_model = "uniform_beta",
      beta_weight = 0.68,
      alpha = 2.8,
      beta = 4.1
    ),
    time_model = "uniform_beta"
  )
}

test_that("temporal KS and CvM use the fitted-CDF formulas", {
  s <- c(0.08, 0.21, 0.37, 0.56, 0.83)
  eta <- temporal_uniform_beta_test_eta()
  result <- sunspots_temporal_uniform_beta_statistics(s, eta)

  u <- sunspots_joint_parsimonious_time_cdf(
    sort(s),
    eta,
    time_model = "uniform_beta"
  )
  n <- length(s)
  index <- seq_len(n)

  expected_ks <- sqrt(n) * max(
    max(index / n - u),
    max(u - (index - 1) / n)
  )
  expected_cvm <- 1 / (12 * n) + sum(
    (u - (2 * index - 1) / (2 * n))^2
  )

  expect_equal(result$ks, expected_ks, tolerance = 1e-14)
  expect_equal(result$cvm, expected_cvm, tolerance = 1e-14)
  expect_equal(result$fitted_cdf, u, tolerance = 1e-14)
})

test_that("temporal parametric bootstrap re-estimates every replicate", {
  set.seed(20260812)
  eta <- temporal_uniform_beta_test_eta()
  s <- sample_sunspots_joint_parsimonious_time(
    300L,
    eta,
    time_model = "uniform_beta"
  )
  fit <- sunspots_temporal_uniform_beta_refit(s)

  bootstrap <- sunspots_temporal_uniform_beta_bootstrap(
    s = s,
    eta_hat = fit,
    B = 8L,
    n_cores = 1L,
    seed = 20260813L
  )

  expect_identical(nrow(bootstrap), 8L)
  expect_identical(bootstrap$replicate, seq_len(8L))
  expect_true(all(is.finite(bootstrap$ks)))
  expect_true(all(is.finite(bootstrap$cvm)))
  expect_true(all(is.finite(bootstrap$beta_weight)))
  expect_true(all(is.finite(bootstrap$alpha)))
  expect_true(all(is.finite(bootstrap$beta)))
  expect_true(all(bootstrap$selected_converged))
})

test_that("temporal bootstrap is deterministic for a fixed seed", {
  set.seed(20260814)
  eta <- temporal_uniform_beta_test_eta()
  s <- sample_sunspots_joint_parsimonious_time(
    250L,
    eta,
    time_model = "uniform_beta"
  )
  fit <- sunspots_temporal_uniform_beta_refit(s)

  first <- sunspots_temporal_uniform_beta_bootstrap(
    s = s,
    eta_hat = fit,
    B = 5L,
    n_cores = 1L,
    seed = 20260815L
  )
  repeated <- sunspots_temporal_uniform_beta_bootstrap(
    s = s,
    eta_hat = fit,
    B = 5L,
    n_cores = 1L,
    seed = 20260815L
  )

  expect_equal(first, repeated, tolerance = 0)
})

test_that("a converged exploratory fit survives failed refinements", {
  set.seed(20260816)
  eta <- temporal_uniform_beta_test_eta()
  s <- sample_sunspots_joint_parsimonious_time(
    350L,
    eta,
    time_model = "uniform_beta"
  )

  fit <- suppressWarnings(
    fit_sunspots_joint_parsimonious_time(
      s,
      time_model = "uniform_beta",
      control = list(
        parsimonious_time_n_starts = 4L,
        parsimonious_time_nelder_mead_control =
          list(maxit = 2500L, reltol = 1e-10),
        parsimonious_time_optim_control =
          list(maxit = 0L)
      )
    )
  )

  expect_true(is.finite(fit$loglik))
  expect_true(fit$selected_converged)
  expect_identical(fit$opt$convergence, 0L)
})

