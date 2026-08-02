library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("bootstrap", "profile_lookup_interpolation.R"))

test_that("C++ cache lock releases an obsolete compilation lock", {
  cache_dir <- tempfile("distance-profile-cache-lock-")
  dir.create(cache_dir)
  lock_dir <- file.path(cache_dir, ".distance_profile_sourcecpp.lock")
  dir.create(lock_dir)
  Sys.setFileTime(lock_dir, Sys.time() - 2)

  expect_true(with_distance_profile_cpp_cache_lock(
    cache_dir,
    TRUE,
    timeout_seconds = 1,
    stale_seconds = 0.1
  ))
  expect_false(dir.exists(lock_dir))
})

test_that("C++ tensor lookup is exact on its local polynomial space", {
  t_grid <- seq(0, 1, length.out = 7L)
  geometry_grid <- seq(-1, 1, length.out = 8L)
  kappa_grid <- seq(1, 3, length.out = 9L)
  grid <- expand.grid(
    t = t_grid,
    geometry = geometry_grid,
    kappa = kappa_grid
  )
  values <- cbind(
    grid$t^5 + 2 * grid$geometry^4 - 0.5 * grid$kappa^3,
    grid$t^2 * grid$geometry^3 * grid$kappa^4
  )
  query <- data.frame(
    t = c(0, 0.13, 0.79, 1),
    geometry = c(-1, -0.37, 0.62, 1),
    kappa = c(1, 1.42, 2.71, 3)
  )

  observed <- distance_profile_cpp_call(
    "cpp_profile_lookup_tensor_local_polynomial",
    values,
    t_grid,
    geometry_grid,
    kappa_grid,
    query$t,
    query$geometry,
    query$kappa,
    6L
  )
  expected <- cbind(
    query$t^5 + 2 * query$geometry^4 - 0.5 * query$kappa^3,
    query$t^2 * query$geometry^3 * query$kappa^4
  )

  expect_equal(observed, expected, tolerance = 2e-12)
})

test_that("C++ tensor lookup refuses extrapolation instead of clipping", {
  grid <- seq(0, 1, length.out = 6L)
  values <- matrix(0, nrow = length(grid)^3, ncol = 1L)

  expect_error(
    distance_profile_cpp_call(
      "cpp_profile_lookup_tensor_local_polynomial",
      values,
      grid,
      grid,
      grid,
      1 + 1e-8,
      0.5,
      0.5,
      6L
    ),
    "outside its certified grid"
  )
})

test_that("invariant lookup reconstructs vMF and HvMF derivatives at table nodes", {
  configurations <- list(
    vmf = list(
      kappa_grid = seq(1, 2, length.out = 6L),
      geometry_grid = seq(-1, 1, length.out = 6L),
      t_grid = seq(0, pi, length.out = 6L)
    ),
    hvmf = list(
      kappa_grid = seq(1, 2, length.out = 6L),
      geometry_grid = seq(0, 1, length.out = 6L),
      t_grid = seq(0, 2, length.out = 6L)
    )
  )

  for (model in names(configurations)) {
    configuration <- configurations[[model]]
    table <- profile_lookup_build(
      model = model,
      q = 2L,
      kappa_grid = configuration$kappa_grid,
      geometry_grid = configuration$geometry_grid,
      t_grid = configuration$t_grid,
      integration_grid_size = 257L,
      cores = 1L
    )
    kappa <- configuration$kappa_grid[[4L]]
    geometry <- configuration$geometry_grid[[3L]]
    thresholds <- configuration$t_grid[c(2L, 5L)]
    mu <- c(1, 0, 0)
    omega <- if (identical(model, "vmf")) {
      c(geometry, sqrt(1 - geometry^2), 0)
    } else {
      c(cosh(geometry), sinh(geometry), 0)
    }
    observed <- profile_lookup_evaluate(
      table = table,
      xi = kappa * mu,
      centers = matrix(omega, nrow = 1L),
      thresholds = thresholds
    )
    expected_full <- if (identical(model, "vmf")) {
      vmf_profile_and_derivative_xi(
        omega = omega,
        xi = kappa * mu,
        t_values = configuration$t_grid,
        distance_type = "geodesic",
        grid_size = 257L
      )
    } else {
      hvmf_profile_and_derivative_xi(
        omega = omega,
        xi = kappa * mu,
        t_values = configuration$t_grid,
        grid_size = 257L
      )
    }
    expected_indices <- c(2L, 5L)

    expect_equal(
      drop(observed$F),
      expected_full$F[expected_indices],
      tolerance = 2e-12
    )
    expect_equal(
      observed$derivative_sorted,
      expected_full$derivative[expected_indices, , drop = FALSE],
      tolerance = 2e-12
    )
  }
})
