library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

clr_distance_manual <- function(x, y) {
  clr_x <- log(x) - mean(log(x))
  clr_y <- log(y) - mean(log(y))
  sqrt(sum((clr_x - clr_y)^2))
}

test_that("logistic Gaussian ilr transform is invertible and matches Aitchison distance", {
  x3 <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.5, 0.2, 0.3),
    c(0.25, 0.35, 0.40)
  )
  z3 <- logistic_gaussian_ilr_matrix(x3)
  x3_back <- logistic_gaussian_ilr_to_simplex(z3, ambient_dim = 3)
  expect_equal(x3_back, x3, tolerance = 1e-12)

  x4 <- rbind(
    c(0.10, 0.20, 0.30, 0.40),
    c(0.25, 0.25, 0.25, 0.25)
  )
  z4 <- logistic_gaussian_ilr_matrix(x4)
  x4_back <- logistic_gaussian_ilr_to_simplex(z4, ambient_dim = 4)
  expect_equal(x4_back, x4, tolerance = 1e-12)

  spec <- make_logistic_gaussian_spec()
  observed_distances <- spec$distance_matrix(x3, x3, control = list())
  expected_distances <- outer(
    seq_len(nrow(x3)),
    seq_len(nrow(x3)),
    Vectorize(function(i, j) clr_distance_manual(x3[i, ], x3[j, ]))
  )

  expect_equal(observed_distances, expected_distances, tolerance = 1e-12)
})

test_that("logistic Gaussian weighted fit agrees with replicated-data fit", {
  set.seed(101)
  x <- rlogistic_gaussian_simplex(
    n = 12,
    mu_ilr = c(0.3, -0.4),
    Sigma_ilr = matrix(c(0.40, 0.10, 0.10, 0.25), nrow = 2L)
  )
  weights <- c(1, 2, 1, 3, 2, 1, 2, 1, 3, 2, 1, 2)

  fit_weighted <- fit_logistic_gaussian_theta(
    data = x,
    weights = weights,
    null = list(type = "composite"),
    unknown_param = "both"
  )
  x_replicated <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_replicated <- fit_logistic_gaussian_theta(
    data = x_replicated,
    null = list(type = "composite"),
    unknown_param = "both"
  )

  expect_equal(fit_weighted$mu_ilr, fit_replicated$mu_ilr, tolerance = 1e-12)
  expect_equal(fit_weighted$Sigma_ilr, fit_replicated$Sigma_ilr, tolerance = 1e-12)
})

test_that("logistic Gaussian bootstrap supports simple and composite nulls", {
  set.seed(202)
  mu_ilr <- c(0.2, -0.3)
  Sigma_ilr <- matrix(c(0.35, 0.08, 0.08, 0.30), nrow = 2L)
  x <- rlogistic_gaussian_simplex(20, mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr)
  ks_grid <- make_logistic_gaussian_ks_grid(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr)

  simple_null <- list(type = "simple", theta = list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr))
  result_1 <- multiplier_bootstrap_logistic_gaussian(
    data = x,
    null = simple_null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1
  )
  result_2 <- multiplier_bootstrap_logistic_gaussian(
    data = x,
    null = simple_null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)
  expect_true(isTRUE(result_1$diagnostics$weighted_mle))

  composite_result <- multiplier_bootstrap_logistic_gaussian(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 84,
    n_cores = 1,
    unknown_param = "both"
  )

  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_equal(length(composite_result$observed$theta_hat$mu_ilr), 2L)
  expect_equal(dim(composite_result$observed$theta_hat$Sigma_ilr), c(2L, 2L))
})

test_that("logistic Gaussian quadratic-form dispatcher returns finite probabilities", {
  normal_prob <- logistic_gaussian_quadform_tail_probability(
    q = 1,
    lambda = c(0.35, 0.25),
    h = c(1, 1),
    delta = c(0.1, 0.2),
    control = list()
  )
  expect_true(is.finite(normal_prob))
  expect_true(normal_prob >= 0 && normal_prob <= 1)

  ill_prob <- logistic_gaussian_quadform_tail_probability(
    q = 247,
    lambda = c(
      8.88904185274295,
      2.81096184712308,
      0.888904185274295,
      0.281096184712308,
      0.0888904185274296,
      0.0281096184712308,
      0.00888904185274296,
      0.00281096184712308,
      0.000888904185274296,
      0.000281096184712308,
      0.0000888904185274297,
      0.0000281096184712308,
      0.00000888904185274296
    ),
    h = rep(1, 13),
    delta = rep(75, 13),
    control = list()
  )
  expect_true(is.finite(ill_prob))
  expect_true(ill_prob >= 0 && ill_prob <= 1)
})

test_that("logistic Gaussian calibration scenarios run in smoke size", {
  simple_output_dir <- file.path(tempdir(), "lg_bootstrap_calibration_smoke")
  simple_result <- run_bootstrap_calibration_study(
    scenarios = default_logistic_gaussian_simple_calibration_scenarios()[1],
    n_values = 24,
    M_outer = 2,
    B = 5,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 321,
    output_dir = simple_output_dir,
    show_progress = FALSE,
    verbose = FALSE
  )

  expect_equal(unique(simple_result$raw_results$model), "logistic_gaussian")
  expect_equal(unique(simple_result$raw_results$scenario), "lg_simple_delta2_balanced")
  expect_equal(sort(unique(simple_result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(simple_result$raw_results), 4)
  expect_true(file.exists(file.path(simple_output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(simple_output_dir, "bootstrap_calibration_summary.csv")))

  composite_output_dir <- file.path(tempdir(), "lg_bootstrap_calibration_composite_smoke")
  composite_result <- run_bootstrap_calibration_study(
    scenarios = default_logistic_gaussian_composite_calibration_scenarios()[1],
    n_values = 24,
    M_outer = 2,
    B = 5,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 654,
    output_dir = composite_output_dir,
    show_progress = FALSE,
    verbose = FALSE
  )

  expect_equal(unique(composite_result$raw_results$model), "logistic_gaussian")
  expect_equal(unique(composite_result$raw_results$scenario), "lg_composite_delta2_balanced")
  expect_equal(sort(unique(composite_result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(composite_result$raw_results), 4)
  expect_true(file.exists(file.path(composite_output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(composite_output_dir, "bootstrap_calibration_summary.csv")))
})
