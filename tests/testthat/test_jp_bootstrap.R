library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

jp_bootstrap_test_control <- list(
  jp_mle_maxit = 250L,
  jp_mle_reltol = 1e-9,
  jp_mle_psi_min = 1e-3,
  jp_mle_psi_abs_starts = c(0.25, 0.5, 1, 2)
)

test_that("JP bootstrap supports both simple and composite nulls", {
  mu <- c(0, 0, 1)
  x <- r_sph_jp(24, mu = mu, kappa = 1, psi = 0.5)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )
  null <- list(type = "simple", theta = list(mu = mu, kappa = 1, psi = 0.5))

  result_1 <- multiplier_bootstrap_jp(
    data = x,
    null = null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1,
    distance_type = "geodesic"
  )

  result_2 <- multiplier_bootstrap_jp(
    data = x,
    null = null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1,
    distance_type = "geodesic"
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)
  expect_equal(result_1$diagnostics$n, nrow(x))

  composite_result <- multiplier_bootstrap_jp(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 24,
    n_cores = 1,
    distance_type = "geodesic",
    control = jp_bootstrap_test_control
  )

  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_equal(composite_result$diagnostics$n, nrow(x))
})

test_that("JP calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "jp_bootstrap_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_jp_simple_calibration_scenarios()[1],
    n_values = 30,
    M_outer = 2,
    B = 5,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 321,
    output_dir = output_dir,
    show_progress = FALSE,
    verbose = FALSE
  )

  expect_equal(unique(result$raw_results$model), "jp")
  expect_equal(unique(result$raw_results$scenario), "jp_simple_s2_geodesic_kappa_1_psi_0p5")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})

test_that("JP composite calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "jp_bootstrap_calibration_composite_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_jp_composite_calibration_scenarios()[1],
    n_values = 30,
    M_outer = 2,
    B = 5,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1,
    seed = 654,
    output_dir = output_dir,
    show_progress = FALSE,
    verbose = FALSE
  )

  expect_equal(unique(result$raw_results$model), "jp")
  expect_equal(unique(result$raw_results$scenario), "jp_composite_s2_geodesic_kappa_1_psi_0p5")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})
