library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)
source(file.path("bootstrap", "multiplier_bootstrap.R"))

test_that("Watson simple bootstrap agrees with Small Circle nu = 0", {
  set.seed(420)
  mu <- c(0, 0, 1)
  x <- r_sph_watson(28, mu, 6)
  ks_grid <- list(omega_grid = generate_canonical_lattice(4L, dim = 3L), t_grid = c(0.4, 0.9))
  common <- list(
    data = x, null = list(type = "simple", theta = list(mu = mu, kappa = 6)),
    statistics = "ks", ks_grid = ks_grid, B = 5L, n_cores = 1L, seed = 421L,
    bootstrap_method = "reestimated", distance_type = "geodesic"
  )
  watson <- do.call(multiplier_bootstrap_watson, common)
  small_circle_args <- common
  small_circle_args$null <- list(type = "simple", theta = list(mu = mu, kappa = 6, nu = 0))
  small_circle <- do.call(multiplier_bootstrap_small_circle, small_circle_args)
  expect_equal(watson$observed$ks$statistic, small_circle$observed$ks$statistic, tolerance = 2e-9)
  expect_equal(watson$bootstrap$statistics$ks, small_circle$bootstrap$statistics$ks, tolerance = 2e-9)
})

test_that("Watson composite bootstrap uses the three-parameter fast branch", {
  set.seed(422)
  x <- r_sph_watson(36, c(0, 0, 1), 8)
  result <- multiplier_bootstrap_watson(
    data = x, null = list(type = "composite"), statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(4L, dim = 3L), t_grid = c(0.4, 0.9)),
    B = 3L, n_cores = 1L, seed = 423L, bootstrap_method = "fast_multiplier",
    distance_type = "geodesic", control = list(
      derivative_mc_size = 500L, derivative_mc_seed = 424L,
      fast_multiplier_vhat_method = "numeric_jacobian"
    )
  )
  expect_identical(result$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_equal(result$diagnostics$S_obs_dim[[2L]], 3L)
  expect_true(all(is.finite(result$bootstrap$statistics$ks)))
})

test_that("Watson fast bootstrap rejects the nonregular uniform fit", {
  isotropic <- rbind(diag(3), -diag(3))
  expect_error(
    multiplier_bootstrap_watson(
      data = isotropic, null = list(type = "composite"), statistics = "ks",
      ks_grid = list(omega_grid = generate_canonical_lattice(3L, dim = 3L), t_grid = 0.8),
      B = 2L, n_cores = 1L, bootstrap_method = "fast_multiplier", distance_type = "geodesic"
    ),
    "not regular"
  )
})
