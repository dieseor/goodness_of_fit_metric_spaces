library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("bootstrap", "calibration_study.R"))

test_that("symmetric small-circle-mixture bootstrap supports simple and composite nulls", {
  set.seed(20260603)
  mu <- c(0, 0, 1)
  x <- r_sph_small_circle_symmetric_mixture2(30, mu = mu, kappa = 8, nu = 0.35)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )

  result_simple <- multiplier_bootstrap_small_circle_symmetric_mixture2(
    data = x,
    null = list(type = "simple", theta = list(mu = mu, kappa = 8, nu = 0.35)),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 42,
    n_cores = 1,
    distance_type = "geodesic",
    control = list(
      small_circle_symmetric_mixture2_L_max = 200L,
      small_circle_symmetric_mixture2_quad_n = 500L
    )
  )

  result_composite <- multiplier_bootstrap_small_circle_symmetric_mixture2(
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
      small_circle_symmetric_mixture2_L_max = 200L,
      small_circle_symmetric_mixture2_quad_n = 500L,
      small_circle_symmetric_mixture2_optim_control = list(maxit = 200L, reltol = 1e-9)
    )
  )

  expect_equal(length(result_simple$bootstrap$statistics$ks), 6)
  expect_equal(length(result_simple$bootstrap$statistics$cvm), 6)
  expect_true(result_simple$inference$ks$p_value >= 0 && result_simple$inference$ks$p_value <= 1)
  expect_true(result_composite$inference$cvm$p_value >= 0 && result_composite$inference$cvm$p_value <= 1)
  expect_equal(length(result_composite$bootstrap$theta_star), 6)
  expect_true(all(vapply(result_composite$bootstrap$theta_star, function(theta) theta$nu >= 0, logical(1))))
  expect_true(all(vapply(result_composite$bootstrap$theta_star, function(theta) {
    theta$mu[which(abs(theta$mu) > 1e-12)[1]] >= 0
  }, logical(1))))
})

test_that("fast symmetric small-circle CvM light prep matches the dense route", {
  set.seed(20260712)
  theta_hat <- list(mu = c(0, 0, 1), kappa = 10, nu = 0.35)
  x <- r_sph_small_circle_symmetric_mixture2(
    36L, mu = theta_hat$mu, kappa = theta_hat$kappa, nu = theta_hat$nu
  )
  spec <- make_small_circle_symmetric_mixture2_spec(distance_type = "geodesic")
  control <- list(
    derivative_mc_size = 300L,
    derivative_mc_seed = 20260713L,
    fast_bootstrap_chunk_size = 3L,
    fast_multiplier_cache_corrections = FALSE,
    small_circle_symmetric_mixture2_profile_method = "legendre",
    small_circle_symmetric_mixture2_L_max = 80L,
    small_circle_symmetric_mixture2_quad_n = 200L
  )
  run_fast <- function(keep_process) {
    multiplier_bootstrap_gof(
      data = x,
      spec = spec,
      null = list(type = "composite"),
      statistics = "cvm",
      B = 6L,
      seed = 20260714L,
      n_cores = 1L,
      observed_theta_hat = theta_hat,
      bootstrap_method = "fast_multiplier",
      keep = list(observed_process = keep_process, bootstrap_statistics = TRUE),
      control = control
    )
  }

  dense <- run_fast(TRUE)
  light <- run_fast(FALSE)

  expect_identical(dense$diagnostics$fast_cvm_mode, "dense_matrix")
  expect_identical(light$diagnostics$fast_cvm_mode, "sample_points_unique_distances_sorted_rows")
  expect_true(light$diagnostics$lightweight_cvm_prep)
  expect_equal(dense$inference$cvm$observed, light$inference$cvm$observed, tolerance = 1e-10)
  expect_equal(dense$bootstrap$statistics$cvm, light$bootstrap$statistics$cvm, tolerance = 1e-10)
})

test_that("symmetric small-circle-mixture simple calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "small_circle_symmetric_mixture2_simple_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_small_circle_symmetric_mixture2_simple_calibration_scenarios(
      kappa_values = 12,
      nu_values = 0.6
    )[1],
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

  expect_equal(unique(result$raw_results$model), "small_circle_symmetric_mixture2")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_raw.csv")))
})

test_that("symmetric small-circle-mixture composite calibration scenarios run in smoke size", {
  output_dir <- file.path(tempdir(), "small_circle_symmetric_mixture2_composite_calibration_smoke")

  result <- run_bootstrap_calibration_study(
    scenarios = default_small_circle_symmetric_mixture2_composite_calibration_scenarios(
      kappa_values = 12,
      nu_values = 0.6
    )[1],
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

  expect_equal(unique(result$raw_results$model), "small_circle_symmetric_mixture2")
  expect_equal(sort(unique(result$raw_results$statistic)), c("cvm", "ks"))
  expect_equal(nrow(result$raw_results), 4)
  expect_true(file.exists(file.path(output_dir, "bootstrap_calibration_summary.csv")))
})
