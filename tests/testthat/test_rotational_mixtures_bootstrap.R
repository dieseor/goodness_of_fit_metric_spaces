library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

test_that("rotational beta-mixture bootstrap supports simple and composite nulls", {
  set.seed(20260603)
  theta <- list(
    mu = c(0, 0, 1),
    weight1 = 0.4,
    alpha1 = 2,
    beta1 = 8,
    alpha2 = 8,
    beta2 = 2
  )
  x <- r_sph_rotational_beta_mixture2(
    n = 24,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )

  result_1 <- multiplier_bootstrap_rotational_beta_mixture2(
    data = x,
    null = list(type = "simple", theta = theta),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 42,
    n_cores = 1
  )
  result_2 <- multiplier_bootstrap_rotational_beta_mixture2(
    data = x,
    null = list(type = "simple", theta = theta),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 42,
    n_cores = 1
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)
  expect_true(isTRUE(result_1$diagnostics$weighted_mle))

  composite_result <- multiplier_bootstrap_rotational_beta_mixture2(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 5,
    seed = 84,
    n_cores = 1,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      rotational_beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_length(composite_result$bootstrap$theta_star, 5)
})

test_that("rotational logit-normal-mixture bootstrap supports simple and composite nulls", {
  set.seed(20260604)
  theta <- list(
    mu = c(0, 0, 1),
    weight1 = 0.45,
    mean1 = -1.1,
    sd1 = 0.45,
    mean2 = 1.0,
    sd2 = 0.55
  )
  x <- r_sph_rotational_logitnormal_mixture2(
    n = 24,
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2
  )
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )

  result_1 <- multiplier_bootstrap_rotational_logitnormal_mixture2(
    data = x,
    null = list(type = "simple", theta = theta),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 24,
    n_cores = 1
  )
  result_2 <- multiplier_bootstrap_rotational_logitnormal_mixture2(
    data = x,
    null = list(type = "simple", theta = theta),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 24,
    n_cores = 1
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)
  expect_true(isTRUE(result_1$diagnostics$weighted_mle))

  composite_result <- multiplier_bootstrap_rotational_logitnormal_mixture2(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 5,
    seed = 48,
    n_cores = 1,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      rotational_logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_length(composite_result$bootstrap$theta_star, 5)
})

test_that("rotational beta-mixture calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "rotational_beta_mixture2_bootstrap_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_rotational_beta_mixture2_simple_calibration_scenarios()[1],
    n_values = 24,
    M_outer = 2,
    B = 4,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 321,
    output_dir = output_dir,
    show_progress = FALSE,
    verbose = FALSE
  )

  expect_equal(unique(result$raw_results$model), "rotational_beta_mixture2")
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})

test_that("rotational logit-normal-mixture calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "rotational_logitnormal_mixture2_bootstrap_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_rotational_logitnormal_mixture2_composite_calibration_scenarios()[1],
    n_values = 24,
    M_outer = 2,
    B = 4,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 654,
    output_dir = output_dir,
    show_progress = FALSE,
    verbose = FALSE
  )

  expect_equal(unique(result$raw_results$model), "rotational_logitnormal_mixture2")
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})
