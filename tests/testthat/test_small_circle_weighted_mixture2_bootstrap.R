library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

test_that("weighted small-circle-mixture KS bootstrap supports simple and composite nulls", {
  set.seed(20260604)
  mu <- c(0, 0, 1)
  x <- r_sph_small_circle_weighted_mixture2(
    n = 32,
    mu = mu,
    pi = 0.6,
    kappa1 = 10,
    nu1 = 0.45,
    kappa2 = 8,
    nu2 = 0.25
  )
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(5, dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 5)
  )

  result_simple <- multiplier_bootstrap_small_circle_weighted_mixture2(
    data = x,
    null = list(type = "simple", theta = list(
      mu = mu, pi = 0.6, kappa1 = 10, nu1 = 0.45, kappa2 = 8, nu2 = 0.25
    )),
    statistics = "ks",
    ks_grid = ks_grid,
    B = 4,
    seed = 42,
    n_cores = 1,
    distance_type = "geodesic",
    control = list(
      small_circle_weighted_mixture2_L_max = 120L,
      small_circle_weighted_mixture2_quad_n = 300L
    )
  )

  result_composite <- multiplier_bootstrap_small_circle_weighted_mixture2(
    data = x,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = ks_grid,
    B = 4,
    seed = 24,
    n_cores = 1,
    distance_type = "geodesic",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      small_circle_weighted_mixture2_L_max = 120L,
      small_circle_weighted_mixture2_quad_n = 300L,
      small_circle_weighted_mixture2_optim_control = list(maxit = 200L, reltol = 1e-9),
      small_circle_weighted_mixture2_n_starts = 4L
    )
  )

  expect_equal(length(result_simple$bootstrap$statistics$ks), 4)
  expect_true(result_simple$inference$ks$p_value >= 0 && result_simple$inference$ks$p_value <= 1)
  expect_true(result_composite$inference$ks$p_value >= 0 && result_composite$inference$ks$p_value <= 1)
  expect_equal(length(result_composite$bootstrap$theta_star), 4)
  expect_true(all(vapply(result_composite$bootstrap$theta_star, function(theta) {
    all(is.finite(c(theta$pi, theta$kappa1, theta$nu1, theta$kappa2, theta$nu2)))
  }, logical(1))))
})
