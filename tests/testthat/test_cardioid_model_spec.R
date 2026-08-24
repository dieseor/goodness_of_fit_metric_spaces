library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("bootstrap", "cardioid_model_spec.R"))

test_that("cardioid geodesic profile matches projected CDF complement", {
  spec <- make_cardioid_spec(k = 2, distance_type = "geodesic")
  theta <- list(mu = c(0, 0, 1), rho = 0.35, k = 2)
  omega <- c(1, 0, 0)
  t_values <- c(0.2, 0.8, 1.4)

  expected <- 1 - p_proj_car_gamma(
    x = cos(t_values),
    rho = theta$rho,
    k = theta$k,
    p = 3,
    mu = theta$mu,
    gamma = omega
  )

  observed <- spec$profile_eval(omega, t_values, theta, control = list())
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("cardioid rho zero recovers the S2 uniform geodesic profile", {
  spec <- make_cardioid_spec(k = 1, distance_type = "geodesic")
  theta <- list(mu = c(0, 0, 1), rho = 0, k = 1)
  omega <- c(1, 0, 0)
  t_values <- c(0.1, 0.7, 1.2, 2.5)

  observed <- spec$profile_eval(omega, t_values, theta, control = list())
  expected <- (1 - cos(t_values)) / 2

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("cardioid rho domain depends on the parity of k", {
  even_theta <- normalize_cardioid_theta(
    list(mu = c(0, 0, 1), rho = -0.4, k = 2)
  )
  expect_equal(even_theta$rho, -0.4)

  expect_error(
    normalize_cardioid_theta(list(mu = c(0, 0, 1), rho = -0.4, k = 3)),
    "rho.*\\[0, 1\\]"
  )
})

test_that("weighted C2 MLE can select the negative rho branch", {
  set.seed(923)
  x <- r_sph_car(n = 500L, mu = c(0, 0, 1), rho = -0.6, k = 2)
  fit <- mle_sph_car_weighted(x, k = 2)

  expect_lt(fit$rho, 0)
})

test_that("odd-order cardioid sampler has the stated Gegenbauer moment", {
  set.seed(847)
  rho <- 0.4

  for (k in c(1L, 3L)) {
    x <- r_sph_car(n = 40000L, mu = c(0, 0, 1), rho = rho, k = k)
    z <- x[, 3L]
    qk <- drop(sphunif::Gegen_polyn(theta = acos(z), k = k, p = 3)) /
      drop(sphunif::Gegen_polyn(theta = 0, k = k, p = 3))

    expect_equal(mean(qk), rho / sphunif::d_p_k(p = 3, k = k), tolerance = 0.012)
  }
})

test_that("cardioid profiles respect geodesic and chordal boundaries", {
  theta <- list(mu = c(0, 0, 1), rho = 0.25, k = 2)
  omega <- c(0, 1, 0)

  spec_geo <- make_cardioid_spec(k = 2, distance_type = "geodesic")
  geo_values <- spec_geo$profile_eval(omega, c(-0.1, 0, pi / 2, pi, pi + 0.1), theta, control = list())
  expect_equal(geo_values[[1]], 0, tolerance = 1e-12)
  expect_equal(geo_values[[2]], 0, tolerance = 1e-12)
  expect_equal(geo_values[[4]], 1, tolerance = 1e-12)
  expect_equal(geo_values[[5]], 1, tolerance = 1e-12)

  spec_chordal <- make_cardioid_spec(k = 2, distance_type = "chordal")
  chordal_values <- spec_chordal$profile_eval(omega, c(-0.1, 0, sqrt(2), 2, 2.1), theta, control = list())
  expect_equal(chordal_values[[1]], 0, tolerance = 1e-12)
  expect_equal(chordal_values[[2]], 0, tolerance = 1e-12)
  expect_equal(chordal_values[[4]], 1, tolerance = 1e-12)
  expect_equal(chordal_values[[5]], 1, tolerance = 1e-12)
})

test_that("weighted cardioid MLE matches the unweighted MLE for equal weights", {
  X <- normalize_cardioid_data(rbind(
    c(0.0, 0.1, 1.0),
    c(0.1, 0.0, 1.0),
    c(-0.1, 0.05, 1.0),
    c(0.05, -0.1, 1.0)
  ))

  common_start <- list(mu = c(0, 0, 1), rho = 0.2, k = 1)
  fit_weighted <- mle_sph_car_weighted(
    X = X,
    k = 1,
    weights = rep(1, nrow(X)),
    theta_start = common_start,
    control = list(cardioid_optim_control = list(maxit = 1000))
  )
  fit_unweighted <- mle_sph_car(
    X = X,
    k = 1,
    mu0 = common_start$mu,
    rho0 = common_start$rho,
    control = list(maxit = 1000)
  )

  expect_equal(fit_weighted$rho, fit_unweighted$rho, tolerance = 1e-6)
  expect_equal(fit_weighted$mu, fit_unweighted$mu, tolerance = 1e-6)
  expect_equal(fit_weighted$ll, fit_unweighted$ll, tolerance = 1e-6)
})

test_that("integer-weight cardioid MLE matches the replicated-sample MLE", {
  X <- normalize_cardioid_data(rbind(
    c(0.0, 0.0, 1.0),
    c(0.1, 0.0, 1.0),
    c(0.0, 0.1, 1.0),
    c(-0.1, 0.0, 1.0)
  ))
  weights <- c(2, 1, 3, 2)
  expanded_X <- X[rep(seq_len(nrow(X)), times = weights), , drop = FALSE]
  common_start <- list(mu = c(0, 0, 1), rho = 0.25, k = 1)

  fit_weighted <- mle_sph_car_weighted(
    X = X,
    k = 1,
    weights = weights,
    theta_start = common_start,
    control = list(cardioid_optim_control = list(maxit = 1000))
  )
  fit_replicated <- mle_sph_car(
    X = expanded_X,
    k = 1,
    mu0 = common_start$mu,
    rho0 = common_start$rho,
    control = list(maxit = 1000)
  )

  expect_equal(fit_weighted$rho, fit_replicated$rho, tolerance = 1e-6)
  expect_equal(fit_weighted$mu, fit_replicated$mu, tolerance = 1e-6)
  expect_equal(fit_weighted$ll / nrow(X), fit_replicated$ll / nrow(expanded_X), tolerance = 1e-6)
})

test_that("cardioid composite bootstrap stores varying theta stars", {
  set.seed(321)
  X <- r_sph_car(
    n = 18,
    mu = c(0, 0, 1),
    rho = 0.45,
    k = 2
  )

  result <- multiplier_bootstrap_gof(
    data = X,
    spec = make_cardioid_spec(k = 2, distance_type = "geodesic"),
    null = list(type = "composite"),
    statistics = "cvm",
    B = 8,
    seed = 99,
    n_cores = 1,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      cardioid_optim_control = list(maxit = 1000)
    )
  )

  theta_star <- result$bootstrap$theta_star
  rho_star <- vapply(theta_star, `[[`, numeric(1), "rho")
  mu_star_1 <- vapply(theta_star, function(theta) theta$mu[[1]], numeric(1))

  expect_length(theta_star, 8)
  expect_true(diff(range(rho_star)) > 0)
  expect_true(diff(range(mu_star_1)) > 0)
  expect_identical(result$diagnostics$engine, "multiplier_bootstrap_gof")
  expect_identical(result$diagnostics$method, "distance_profiles")
  expect_true(isTRUE(result$diagnostics$weighted_mle))
})

test_that("joint fast KS and CvM refuse an unshared dense evaluation", {
  set.seed(812)
  X <- r_sph_car(n = 12L, mu = c(0, 0, 1), rho = 0.35, k = 2L)

  expect_error(
    multiplier_bootstrap_gof(
      data = X,
      spec = make_cardioid_spec(k = 2L, distance_type = "geodesic"),
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = 2L,
      seed = 91L,
      n_cores = 1L,
      bootstrap_method = "fast_multiplier",
      keep = list(
        observed_process = TRUE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = list(derivative_mc_size = 50L)
    ),
    "shared sample preparation and fused evaluation are not active",
    fixed = TRUE
  )
})

test_that("joint fast KS and CvM use the shared fused C++ path", {
  set.seed(813)
  X <- r_sph_car(n = 12L, mu = c(0, 0, 1), rho = 0.35, k = 2L)

  result <- multiplier_bootstrap_gof(
    data = X,
    spec = make_cardioid_spec(k = 2L, distance_type = "geodesic"),
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 3L,
    seed = 92L,
    n_cores = 1L,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_mc_size = 50L,
      fast_multiplier_backend = "cpp",
      fast_multiplier_fuse_ks_cvm = TRUE
    )
  )

  expect_true(isTRUE(result$diagnostics$shared_sample_ks_cvm_cache))
  expect_true(isTRUE(result$diagnostics$fast_multiplier_fuse_ks_cvm_effective))
  expect_identical(result$diagnostics$fast_multiplier_backend_effective, "cpp")
  expect_identical(
    result$diagnostics$fast_cvm_mode,
    "sample_points_unique_distances_sorted_rows"
  )
})
