library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

test_that("HvMF calibration simulates polar H^2 samples in memory", {
  scenario <- make_hvmf_composite_calibration_scenario(5)
  sample_data <- simulate_h0_sample(scenario = scenario, n = 50)

  expect_true(is.matrix(sample_data))
  expect_equal(dim(sample_data), c(50, 3))
  minkowski_norms <- -sample_data[, 1]^2 + rowSums(sample_data[, -1, drop = FALSE]^2)
  expect_true(all(abs(minkowski_norms + 1) < 1e-10))
  expect_true(all(sample_data[, 1] > 0))
})

test_that("HvMF composite CvM calibration runs with polar samples", {
  scenario <- make_hvmf_composite_calibration_scenario(5)
  results <- run_calibration_scenario(
    scenario = scenario,
    n_values = 50,
    M_outer = 2,
    B = 2,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = "cvm",
    n_cores_outer = 1,
    seed = 123
  )

  expect_equal(nrow(results), 2)
  expect_true(all(results$model == "hvmf"))
  expect_true(all(results$statistic == "cvm"))
  expect_true(all(results$null_type == "composite"))
  expect_true(all(results$unknown_param == "both"))
  expect_true(all(results$p_value >= 0 & results$p_value <= 1))
})

test_that("HvMF simple KS and CvM calibration runs with polar samples", {
  scenario <- make_hvmf_simple_calibration_scenario(5)
  sample_data <- simulate_h0_sample(scenario = scenario, n = 50)
  ks_grid <- make_hvmf_ks_grid(sample_data, mu = scenario$sample_params$mu)

  expect_equal(dim(ks_grid$omega_grid), c(10, 3))
  expect_equal(length(ks_grid$t_grid), 10)
  expect_true(max(acosh(sample_data[, 1])) <= max(acosh(ks_grid$omega_grid[, 1])) + 1e-12)

  results <- run_calibration_scenario(
    scenario = scenario,
    n_values = 50,
    M_outer = 1,
    B = 2,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 123
  )

  expect_equal(sort(unique(results$statistic)), c("cvm", "ks"))
  expect_true(all(results$model == "hvmf"))
  expect_true(all(results$null_type == "simple"))
  expect_true(all(results$p_value >= 0 & results$p_value <= 1))
})
