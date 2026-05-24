library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

test_that("HvMF calibration loader reads prepared H^2 samples from disk", {
  scenario <- make_hvmf_composite_calibration_scenario(50)
  sample_data <- load_hvmf_calibration_sample(
    scenario = scenario,
    n = 50,
    replicate_id = 1
  )

  expect_true(is.matrix(sample_data))
  expect_equal(dim(sample_data), c(50, 3))
  minkowski_norms <- -sample_data[, 1]^2 + rowSums(sample_data[, -1, drop = FALSE]^2)
  expect_true(all(abs(minkowski_norms + 1) < 1e-10))
  expect_true(all(sample_data[, 1] > 0))
})

test_that("HvMF composite CvM calibration runs on prepared disk samples", {
  scenario <- make_hvmf_composite_calibration_scenario(50)
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
