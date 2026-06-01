library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

test_that("spherical Cauchy bootstrap supports simple and composite nulls", {
  mu <- c(0, 0, 1)
  x <- r_sph_spherical_cauchy(24, mu = mu, rho = 0.35)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )

  result_1 <- multiplier_bootstrap_spherical_cauchy(
    data = x,
    null = list(type = "simple", theta = list(mu = mu, rho = 0.35)),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1,
    distance_type = "geodesic"
  )
  result_2 <- multiplier_bootstrap_spherical_cauchy(
    data = x,
    null = list(type = "simple", theta = list(mu = mu, rho = 0.35)),
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

  composite_result <- multiplier_bootstrap_spherical_cauchy(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 24,
    n_cores = 1,
    distance_type = "geodesic",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      spherical_cauchy_maxit = 300L,
      spherical_cauchy_profile_warn = FALSE
    )
  )

  rho_star <- vapply(composite_result$bootstrap$theta_star, `[[`, numeric(1), "rho")
  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_equal(composite_result$diagnostics$n, nrow(x))
  expect_length(composite_result$bootstrap$theta_star, 8)
  expect_true(diff(range(rho_star)) > 0)
})

test_that("spherical Cauchy calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "spherical_cauchy_bootstrap_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_spherical_cauchy_simple_calibration_scenarios(rho_values = 0.4)[1],
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

  expect_equal(unique(result$raw_results$model), "spherical_cauchy")
  expect_equal(unique(result$raw_results$scenario), "spherical_cauchy_simple_s2_geodesic_rho_0p4")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})

test_that("spherical Cauchy composite calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "spherical_cauchy_bootstrap_calibration_composite_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_spherical_cauchy_composite_calibration_scenarios(rho_values = 0.4)[1],
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

  expect_equal(unique(result$raw_results$model), "spherical_cauchy")
  expect_equal(unique(result$raw_results$scenario), "spherical_cauchy_composite_s2_geodesic_rho_0p4")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})