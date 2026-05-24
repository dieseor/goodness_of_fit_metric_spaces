library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

hyperboloid_point <- function(chi, theta) {
  c(cosh(chi), sinh(chi) * cos(theta), sinh(chi) * sin(theta))
}

test_that("HvMF normalization and closed-form weighted fit are consistent", {
  data <- rbind(
    hyperboloid_point(0.7, 0.1),
    hyperboloid_point(0.9, 0.5),
    hyperboloid_point(1.1, -0.3),
    hyperboloid_point(0.8, 0.9)
  )
  weights <- c(0.5, 2, 1.5, 3)

  normalized_data <- normalize_hvmf_data(data)
  expect_equal(normalized_data, data, tolerance = 1e-12)

  theta <- normalize_hvmf_theta(list(mu = data[1, ], kappa = 7))
  expect_equal(theta$mu, data[1, ], tolerance = 1e-12)
  expect_equal(theta$kappa, 7, tolerance = 1e-12)
  expect_equal(minkowski_inner_product(theta$mu, theta$mu), -1, tolerance = 1e-10)

  fit <- fit_hvmf_theta(
    data = data,
    weights = weights,
    null = list(type = "composite")
  )
  manual <- hvmf_mle_h2(data, weights = weights)

  expect_equal(fit$mu, manual$mu, tolerance = 1e-12)
  expect_equal(fit$kappa, manual$kappa, tolerance = 1e-12)
  expect_equal(minkowski_inner_product(fit$mu, fit$mu), -1, tolerance = 1e-10)
})

test_that("HvMF distance helper and model spec use the hyperbolic geodesic distance", {
  omega_grid <- rbind(
    hyperboloid_point(0.6, 0.0),
    hyperboloid_point(0.8, 0.6)
  )
  data <- rbind(
    hyperboloid_point(0.7, 0.1),
    hyperboloid_point(0.9, 0.5),
    hyperboloid_point(1.0, -0.2)
  )

  helper_distances <- hvmf_distance_matrix(omega_grid, data)
  expect_equal(dim(helper_distances), c(2, 3))
  expect_equal(
    helper_distances[1, 1],
    hyperbolic_geodesic_distance_h2(omega_grid[1, ], data[1, ]),
    tolerance = 1e-12
  )

  spec <- make_hvmf_spec()
  spec_distances <- spec$distance_matrix(data, omega_grid, control = list())
  expect_equal(dim(spec_distances), c(3, 2))
  expect_equal(spec_distances, t(helper_distances), tolerance = 1e-12)
})

test_that("HvMF composite bootstrap runs through the generic infrastructure", {
  data <- rbind(
    hyperboloid_point(0.7, 0.1),
    hyperboloid_point(0.9, 0.5),
    hyperboloid_point(1.1, -0.3),
    hyperboloid_point(0.8, 0.9)
  )

  result <- multiplier_bootstrap_gof(
    data = data,
    spec = make_hvmf_spec(),
    null = list(type = "composite"),
    statistics = "cvm",
    B = 2,
    seed = 123,
    n_cores = 1,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE
    )
  )

  expect_true(is.list(result$observed$theta_hat))
  expect_true(all(is.finite(result$observed$theta_hat$mu)))
  expect_true(is.finite(result$observed$theta_hat$kappa))
  expect_true(result$observed$theta_hat$kappa > 0)
  expect_length(result$bootstrap$statistics$cvm, 2)
  expect_true(all(is.finite(result$bootstrap$statistics$cvm)))
  expect_true(result$inference$cvm$p_value >= 0 && result$inference$cvm$p_value <= 1)
})

test_that("HvMF tabulated CvM matrix matches the exact matrix closely for kappa=50", {
  data <- rbind(
    hyperboloid_point(0.55, -0.1),
    hyperboloid_point(0.70, 0.25),
    hyperboloid_point(0.85, -0.35),
    hyperboloid_point(0.95, 0.60),
    hyperboloid_point(0.65, 1.05),
    hyperboloid_point(1.05, -0.75)
  )
  theta <- hvmf_mle_h2(data)
  theta$kappa <- 50

  spec <- make_hvmf_spec()
  distance_matrix <- spec$distance_matrix(data, data, control = list())

  exact_matrix <- compute_theoretical_sample_profile_matrix(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(hvmf_profile_method = "exact")
  )
  fast_matrix <- compute_theoretical_sample_profile_matrix(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(hvmf_profile_method = "tabulated", hvmf_profile_n_y = 4097L)
  )

  expect_true(max(abs(exact_matrix - fast_matrix)) < 1e-3)
})

test_that("HvMF tabulated CvM matrix matches the exact matrix closely for kappa=200", {
  data <- rbind(
    hyperboloid_point(0.55, -0.1),
    hyperboloid_point(0.70, 0.25),
    hyperboloid_point(0.85, -0.35),
    hyperboloid_point(0.95, 0.60),
    hyperboloid_point(0.65, 1.05),
    hyperboloid_point(1.05, -0.75)
  )
  theta <- hvmf_mle_h2(data)
  theta$kappa <- 200

  spec <- make_hvmf_spec()
  distance_matrix <- spec$distance_matrix(data, data, control = list())

  exact_matrix <- compute_theoretical_sample_profile_matrix(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(hvmf_profile_method = "exact")
  )
  fast_matrix <- compute_theoretical_sample_profile_matrix(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(hvmf_profile_method = "tabulated", hvmf_profile_n_y = 4097L)
  )

  expect_true(max(abs(exact_matrix - fast_matrix)) < 1e-3)
})
