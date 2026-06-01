library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

test_that("Small Circle bootstrap supports simple and composite nulls", {
  set.seed(20260601)
  mu <- c(0, 0, 1)
  x <- r_sph_small_circle(30, mu = mu, kappa = 8, nu = 0.4)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )

  result_simple <- multiplier_bootstrap_small_circle(
    data = x,
    null = list(type = "simple", theta = list(mu = mu, kappa = 8, nu = 0.4)),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 42,
    n_cores = 1,
    distance_type = "geodesic",
    control = list(small_circle_L_max = 200L, small_circle_quad_n = 500L)
  )

  result_composite <- multiplier_bootstrap_small_circle(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 24,
    n_cores = 1,
    distance_type = "geodesic",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      small_circle_L_max = 200L,
      small_circle_quad_n = 500L,
      small_circle_optim_control = list(maxit = 200L, reltol = 1e-9)
    )
  )

  expect_equal(length(result_simple$bootstrap$statistics$ks), 6)
  expect_equal(length(result_simple$bootstrap$statistics$cvm), 6)
  expect_true(result_simple$inference$ks$p_value >= 0 && result_simple$inference$ks$p_value <= 1)
  expect_true(result_composite$inference$cvm$p_value >= 0 && result_composite$inference$cvm$p_value <= 1)
  expect_equal(length(result_composite$bootstrap$theta_star), 6)
  expect_true(all(vapply(result_composite$bootstrap$theta_star, function(theta) theta$nu >= 0, logical(1))))
})

test_that("Small Circle calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "small_circle_bootstrap_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_small_circle_simple_calibration_scenarios(kappa_values = 5, nu_values = 0.3)[1],
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

  expect_equal(unique(result$raw_results$model), "small_circle")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
})

test_that("Small Circle composite calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "small_circle_bootstrap_calibration_composite_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_small_circle_composite_calibration_scenarios(kappa_values = 5, nu_values = 0.3)[1],
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

  expect_equal(unique(result$raw_results$model), "small_circle")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})