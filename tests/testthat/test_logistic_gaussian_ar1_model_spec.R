library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "logistic_gaussian_ar1_bootstrap.R"))
source(file.path("scripts", "run_logistic_gaussian_sigma_shape_scenarios.R"))

test_that("AR(1) covariance has the intended structure", {
  R <- logistic_gaussian_ar1_covariance(5, 0.5)

  expect_equal(diag(R), rep(1, 5), tolerance = 0)
  expect_equal(R[1, ], c(1, 0.5, 0.25, 0.125, 0.0625),
               tolerance = 1e-14)
  expect_gt(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("analytic AR(1) covariance derivative matches finite differences", {
  d <- 5L
  rho <- 0.35
  h <- 1e-6

  analytic <- logistic_gaussian_ar1_covariance_derivative(d, rho)
  numeric <- (
    logistic_gaussian_ar1_covariance(d, rho + h) -
      logistic_gaussian_ar1_covariance(d, rho - h)
  ) / (2 * h)

  expect_equal(analytic, numeric, tolerance = 1e-7)
})

test_that("AR(1) MLE estimates weighted mean and gives nearly zero score", {
  set.seed(610)
  d <- 5L

  x <- rlogistic_gaussian_ar1(
    n = 200,
    mu_ilr = c(0.2, -0.1, 0.3, 0, 0.15),
    rho = 0.5
  )

  weights <- seq_len(nrow(x))
  probability_weights <- weights / sum(weights)
  z <- logistic_gaussian_ilr_matrix(x)

  fit <- fit_logistic_gaussian_ar1_theta(
    data = x,
    weights = weights,
    null = list(type = "composite")
  )

  expect_equal(
    fit$mu_ilr,
    colSums(z * probability_weights),
    tolerance = 1e-12
  )
  expect_true(is.finite(fit$rho))
  expect_lt(abs(fit$rho), 1)

  score <- logistic_gaussian_ar1_score_matrix(x, fit)
  weighted_score_mean <- colSums(
    sweep(score, 1L, probability_weights, FUN = "*")
  )
  expect_lt(max(abs(weighted_score_mean)), 1e-5)
})

test_that("AR(1) Fisher information is positive definite", {
  theta <- list(
    mu_ilr = rep(0, 5),
    rho = 0.5
  )

  information <- logistic_gaussian_ar1_fisher_information(theta)
  eig <- eigen(information, symmetric = TRUE, only.values = TRUE)$values

  expect_gt(min(eig), 0)
})

test_that("t scenario null generator uses AR(1) covariance", {
  set.seed(620)
  d <- 5L
  rho0 <- 0.5

  x <- generate_lg_sigma_shape_sample(
    scenario = "t4",
    n = 4000,
    d = d,
    beta = 0,
    nu = 4,
    t_standardized = TRUE,
    t_ar1_rho = rho0
  )

  z <- logistic_gaussian_ilr_matrix(x)
  expect_equal(
    stats::cov(z),
    logistic_gaussian_ar1_covariance(d, rho0),
    tolerance = 0.08
  )
})

test_that("AR(1) t contaminant is standardized to the null covariance", {
  skip_if_not_installed("mvtnorm")
  set.seed(630)
  d <- 5L
  rho0 <- 0.5

  x <- r_t_logistic_ar1(
    n = 8000,
    d = d,
    rho = rho0,
    nu = 4,
    standardized = TRUE
  )

  z <- logistic_gaussian_ilr_matrix(x)
  expect_equal(
    stats::cov(z),
    logistic_gaussian_ar1_covariance(d, rho0),
    tolerance = 0.12
  )
})

test_that("AR(1) fast bootstrap runs in smoke size", {
  set.seed(640)

  x <- rlogistic_gaussian_ar1(
    n = 25,
    mu_ilr = c(0, 0),
    rho = 0.4
  )

  result <- multiplier_bootstrap_logistic_gaussian_ar1(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 5,
    seed = 123,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 100,
      derivative_mc_seed = 321,
      fast_multiplier_backend = "r"
    ),
    fast_multiplier_backend = "r",
    fuse_ks_cvm = TRUE
  )

  expect_true(result$inference$ks$p_value >= 0 &&
              result$inference$ks$p_value <= 1)
  expect_true(result$inference$cvm$p_value >= 0 &&
              result$inference$cvm$p_value <= 1)
  expect_equal(length(result$observed$theta_hat$mu_ilr), 2L)
  expect_true(is.finite(result$observed$theta_hat$rho))
  expect_lt(abs(result$observed$theta_hat$rho), 1)
})

test_that("AR(1) scalar profile likelihood equals the matrix likelihood", {
  set.seed(701)

  for (d in c(2L, 5L)) {
    x <- rlogistic_gaussian_ar1(
      n = 80,
      mu_ilr = rep(0, d),
      rho = 0.5
    )

    z <- logistic_gaussian_ilr_matrix(x)
    mu_hat <- colMeans(z)
    centered <- sweep(z, 2L, mu_hat, FUN = "-")
    weights <- rep(1 / nrow(z), nrow(z))

    moments <- logistic_gaussian_ar1_profile_moments(
      centered,
      weights
    )

    for (rho in c(-0.7, -0.2, 0, 0.3, 0.5, 0.8)) {
      Sigma <- logistic_gaussian_ar1_covariance(d, rho)
      Sigma_inv <- solve(Sigma)

      matrix_value <- -0.5 * (
        as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) +
          sum(
            weights *
              rowSums((centered %*% Sigma_inv) * centered)
          )
      )

      scalar_value <-
        logistic_gaussian_ar1_profile_loglik_from_moments(
          rho,
          moments
        )

      expect_equal(
        scalar_value,
        matrix_value,
        tolerance = 1e-11
      )
    }
  }
})

test_that("AR(1) cubic roots contain the profile-likelihood maximizer", {
  set.seed(702)

  for (d in c(2L, 5L)) {
    for (rho0 in c(-0.6, 0, 0.3, 0.7)) {
      x <- rlogistic_gaussian_ar1(
        n = 120,
        mu_ilr = rep(0, d),
        rho = rho0
      )

      z <- logistic_gaussian_ilr_matrix(x)
      centered <- sweep(z, 2L, colMeans(z), FUN = "-")
      weights <- rep(1 / nrow(z), nrow(z))

      global <- fit_logistic_gaussian_ar1_rho_global(
        centered,
        weights
      )

      legacy <- fit_logistic_gaussian_ar1_rho_optimize(
        centered,
        weights
      )

      # Global candidate evaluation cannot be worse than the legacy
      # one-dimensional optimizer on the same interval.
      moments <- logistic_gaussian_ar1_profile_moments(
        centered,
        weights
      )

      legacy_exact <-
        logistic_gaussian_ar1_profile_loglik_from_moments(
          legacy$rho,
          moments
        )

      expect_gte(
        global$objective + 1e-11,
        legacy_exact
      )

      # In ordinary interior cases both approaches should agree closely.
      if (!global$at_boundary) {
        expect_equal(
          global$rho,
          legacy$rho,
          tolerance = 1e-6
        )
      }
    }
  }
})

test_that("AR(1) global MLE satisfies the rho score equation when interior", {
  set.seed(703)

  for (d in c(2L, 5L)) {
    x <- rlogistic_gaussian_ar1(
      n = 200,
      mu_ilr = rep(0, d),
      rho = 0.5
    )

    fit <- fit_logistic_gaussian_ar1_theta(
      x,
      null = list(type = "composite")
    )

    expect_false(fit$fit_diagnostics$rho_at_boundary)

    score <- logistic_gaussian_ar1_score_matrix(x, fit)

    expect_lt(
      abs(mean(score[, "rho"])),
      1e-9
    )
  }
})

test_that("AR(1) weighted global MLE satisfies weighted score equations", {
  set.seed(704)

  for (d in c(2L, 5L)) {
    x <- rlogistic_gaussian_ar1(
      n = 150,
      mu_ilr = rep(0, d),
      rho = 0.4
    )

    weights <- rexp(nrow(x))

    fit <- fit_logistic_gaussian_ar1_theta(
      x,
      weights = weights,
      null = list(type = "composite")
    )

    probability_weights <- weights / sum(weights)
    score <- logistic_gaussian_ar1_score_matrix(x, fit)

    weighted_score <- colSums(
      sweep(score, 1L, probability_weights, FUN = "*")
    )

    expect_lt(
      max(abs(weighted_score)),
      1e-8
    )
  }
})
