library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

jp_test_control <- list(
  jp_mle_maxit = 250L,
  jp_mle_reltol = 1e-9,
  jp_mle_psi_min = 1e-3,
  jp_mle_psi_abs_starts = c(0.25, 0.5, 1, 2)
)

test_that("closed-form JP normalizing constant on S^2 matches numerical integration", {
  parameter_grid <- list(
    c(alpha = 0, beta = 2),
    c(alpha = 0.2, beta = 1.5),
    c(alpha = -0.4, beta = -0.75),
    c(alpha = 0.35, beta = -1),
    c(alpha = -0.6, beta = -0.999999)
  )

  for (params in parameter_grid) {
    expected <- c_jp_sphere(q = 2, alpha = params[["alpha"]], beta = params[["beta"]])
    observed <- exp(jp_log_norm_constant_s2(params[["alpha"]], params[["beta"]]))
    expect_equal(observed, expected, tolerance = 1e-8)
  }
})

test_that("JP weighted log-likelihood on S^2 is finite for valid parameters", {
  set.seed(101)
  x <- r_sph_jp(40, mu = c(0, 0, 1), kappa = 1.25, psi = 0.5)
  prob_weights <- rep.int(1 / nrow(x), nrow(x))

  expect_true(is.finite(jp_weighted_loglik_s2(c(0, 0, 1), 1.25, 0.5, x, prob_weights)))
  expect_true(is.finite(jp_weighted_loglik_s2(c(0, 0, 1), 1.5, -0.5, x, prob_weights)))
  expect_true(is.finite(jp_weighted_loglik_s2(c(0, 0, 1), 1.25, 0, x, prob_weights)))
})

test_that("JP composite weighted fit agrees with replicated-data fit", {
  set.seed(202)
  x <- r_sph_jp(24, mu = c(0, 0, 1), kappa = 1.4, psi = 0.5)
  weights <- c(1, 2, 1, 3, 2, 1, 2, 3, 1, 2, 1, 2, 3, 1, 2, 1, 3, 2, 1, 2, 1, 3, 2, 1)

  fit_weighted <- jp_mle_s2_weighted(x, weights = weights, control = jp_test_control)
  x_replicated <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_replicated <- jp_mle_s2_weighted(x_replicated, control = jp_test_control)

  expect_true(sum(fit_weighted$mu * fit_replicated$mu) > 0.995)
  expect_equal(fit_weighted$kappa, fit_replicated$kappa, tolerance = 0.12)
  expect_equal(fit_weighted$psi, fit_replicated$psi, tolerance = 0.12)
})

test_that("JP composite weighted fit accepts observed-fit warm starts", {
  set.seed(252)
  x <- r_sph_jp(36, mu = c(0, 0, 1), kappa = 1.4, psi = 0.5)

  fast_control <- modifyList(jp_test_control, list(
    jp_mle_maxit = 80L,
    jp_mle_reltol = 1e-6,
    jp_mle_psi_abs_starts = c(0.5),
    jp_mle_sign_branches = c(1L)
  ))

  theta_start <- jp_mle_s2_weighted(x, control = fast_control)
  fit_warm <- jp_mle_s2_weighted(
    x,
    weights = runif(nrow(x), min = 0.75, max = 1.25),
    control = modifyList(fast_control, list(jp_mle_start_theta = theta_start))
  )

  expect_true(is.finite(fit_warm$kappa))
  expect_true(is.finite(fit_warm$psi))
  if (!is.null(fit_warm$loglik)) {
    expect_true(is.finite(fit_warm$loglik))
  }

  expect_equal(sum(fit_warm$mu^2), 1, tolerance = 1e-8)
})

test_that("JP composite MLE recovers the sign of psi", {
  set.seed(303)
  mu <- c(0, 0, 1)

  x_pos <- r_sph_jp(160, mu = mu, kappa = 1.5, psi = 0.5)
  fit_pos <- jp_mle_s2_weighted(x_pos, control = jp_test_control)
  expect_true(sum(fit_pos$mu * mu) > 0.9)
  expect_true(fit_pos$kappa > 0)
  expect_true(fit_pos$psi > 0)

  x_neg <- r_sph_jp(160, mu = mu, kappa = 1.8, psi = -0.5)
  fit_neg <- jp_mle_s2_weighted(x_neg, control = jp_test_control)
  expect_true(sum(fit_neg$mu * mu) > 0.9)
  expect_true(fit_neg$kappa > 0)
  expect_true(fit_neg$psi < 0)
})

test_that("JP composite fit stays close to the vMF submodel on vMF data", {
  set.seed(404)
  mu <- c(0, 0, 1)
  x <- rotasym::r_vMF(n = 240, mu = mu, kappa = 2)

  fit <- jp_mle_s2_weighted(
    x,
    control = modifyList(jp_test_control, list(jp_mle_psi_min = 0.05))
  )

  expect_true(abs(fit$psi) < 0.2)
  expect_true(sum(fit$mu * mu) > 0.9)
  expect_true(fit$kappa > 0)
})
