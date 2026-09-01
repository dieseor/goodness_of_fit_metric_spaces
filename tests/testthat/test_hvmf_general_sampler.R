library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source("utils.R")

hvmf_section6_mu <- function(q) {
  c(sqrt(2), 1, rep.int(0, q - 1L))
}

hvmf_spatial_direction <- function(x) {
  spatial <- x[, -1L, drop = FALSE]
  spatial / sqrt(rowSums(spatial^2))
}

vmf_mean_resultant_ratio <- function(lambda, q) {
  lambda <- as.numeric(lambda)
  output <- numeric(length(lambda))
  positive <- lambda > sqrt(.Machine$double.eps)
  output[positive] <- besselI(lambda[positive], nu = q / 2, expon.scaled = TRUE) /
    besselI(lambda[positive], nu = q / 2 - 1, expon.scaled = TRUE)
  output
}

hvmf_regularised_population_distance <- function(q,
                                                  kappa,
                                                  chi = asinh(1),
                                                  p_max = 0.999) {
  table <- hvmf_build_radial_quantile_table(
    q = q,
    kappa = kappa,
    chi = chi,
    p_max = p_max,
    probability_step = 0.01,
    upper = 3
  )
  probability <- seq(0, p_max, length.out = 2001L)
  interpolated_quantile <- stats::approx(
    x = table$base,
    y = table$quantile_values,
    xout = probability
  )$y
  true_cdf <- hvmf_radial_cdf(
    u = interpolated_quantile,
    q = q,
    kappa = kappa,
    chi = chi
  )

  max(abs(probability / p_max - true_cdf), 1 - p_max)
}

test_that("the default HvMF polar cutoff is 0.999", {
  expect_equal(formals(hvmf_build_radial_quantile_table)$p_max, 0.999)
  expect_equal(formals(rhvmf_polar)$p_max, 0.999)
  expect_equal(formals(rhvmf_angular_mixture)$p_max, 0.999)
})

test_that("general polar sampler returns samples on H^2 and H^10", {
  for (q in c(2L, 10L)) {
    set.seed(2100 + q)
    x <- rhvmf_polar(
      n = 40,
      mu = hvmf_section6_mu(q),
      kappa = q
    )

    expect_equal(dim(x), c(40, q + 1L))
    expect_true(all(x[, 1L] > 0))
    minkowski_norm <- -x[, 1L]^2 + rowSums(x[, -1L, drop = FALSE]^2)
    expect_lt(max(abs(minkowski_norm + 1)), 1e-8)
  }
})

test_that("H2 wrappers are exactly compatible with the general sampler", {
  mu <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))

  set.seed(2201)
  general_null <- rhvmf_polar(30, mu = mu, kappa = 5)
  set.seed(2201)
  legacy_null <- rhvmf_h2_polar(30, mu = mu, kappa = 5)
  expect_equal(legacy_null, general_null, tolerance = 0)

  set.seed(2202)
  general_alternative <- rhvmf_angular_mixture(
    30,
    mu = mu,
    kappa = 5,
    delta = 0.2
  )
  set.seed(2202)
  legacy_alternative <- rhvmf_h2_angular_mixture(
    30,
    mu = mu,
    kappa = 5,
    delta = 0.2
  )
  expect_equal(legacy_alternative, general_alternative, tolerance = 0)
})

test_that("zero angular displacement reproduces the general null sampler", {
  mu <- hvmf_section6_mu(10L)

  set.seed(2301)
  null_sample <- rhvmf_polar(25, mu = mu, kappa = 10)
  set.seed(2301)
  zero_displacement <- rhvmf_angular_mixture(
    25,
    mu = mu,
    kappa = 10,
    delta = 0
  )

  expect_equal(zero_displacement, null_sample, tolerance = 0)
})

test_that("general sampler validates dimension, concentration, and tangent", {
  mu10 <- hvmf_section6_mu(10L)

  expect_error(rhvmf_polar(5, mu = c(1, 0), kappa = 1), "q \\+ 1")
  expect_error(rhvmf_polar(5, mu = mu10, kappa = 0), "strictly positive")
  expect_error(
    rhvmf_polar(5, mu = c(2, rep(0, 10)), kappa = 1),
    "<x, x>_M = -1"
  )
  expect_error(
    rhvmf_angular_mixture(5, mu = mu10, kappa = 10, delta = 0.2),
    "tangent.*required"
  )
  expect_error(
    rhvmf_angular_mixture(
      5,
      mu = mu10,
      kappa = 10,
      delta = 0.2,
      tangent = c(1, rep(0, 9))
    ),
    "nonzero component orthogonal"
  )
})

test_that("general HvMF MLE and tabulated profiles work in H^10", {
  q <- 10L
  mu <- hvmf_section6_mu(q)
  set.seed(2401)
  x <- rhvmf_polar(180L, mu = mu, kappa = 10)
  fit <- hvmf_mle_hq(x)

  expect_identical(fit$q, q)
  expect_true(is.finite(fit$kappa) && fit$kappa > 0)
  expect_lt(abs(hvmf_minkowski_inner_product(fit$mu, fit$mu) + 1), 1e-8)

  radii <- seq(0, 2, length.out = 17L)
  omega <- c(1, rep.int(0, q))
  chi <- acosh(-hvmf_minkowski_inner_product(mu, omega))
  tabulated <- hvmf_distance_profile_tabulated(
    radii, q = q, kappa = 10, chi = chi, grid_size = 4097L
  )
  direct <- hvmf_radial_cdf(radii, q = q, kappa = 10, chi = chi)
  expect_lt(max(abs(tabulated - direct)), 5e-4)
})

test_that("general q=2 radial density equals the legacy closed form", {
  kappa <- 3
  chi <- asinh(1)
  u <- seq(0.05, 2.8, length.out = 100L)
  lambda <- kappa * sinh(chi) * sinh(u)
  legacy_density <-
    kappa * exp(kappa) * sinh(u) *
    exp(-kappa * cosh(chi) * cosh(u)) *
    besselI(lambda, nu = 0)

  expect_equal(
    hvmf_radial_density(u, q = 2, kappa = kappa, chi = chi),
    legacy_density,
    tolerance = 1e-11,
    scale = 1
  )
})

test_that("radial and mean-projection densities integrate to one", {
  for (configuration in list(
    list(q = 2L, kappa = 2),
    list(q = 2L, kappa = 3),
    list(q = 10L, kappa = 10),
    list(q = 10L, kappa = 15)
  )) {
    radial_mass <- stats::integrate(
      hvmf_radial_density,
      lower = 0,
      upper = Inf,
      q = configuration$q,
      kappa = configuration$kappa,
      chi = asinh(1),
      rel.tol = 1e-8
    )$value
    projection_mass <- stats::integrate(
      hvmf_mean_projection_density,
      lower = 1,
      upper = Inf,
      q = configuration$q,
      kappa = configuration$kappa,
      rel.tol = 1e-8
    )$value

    expect_equal(radial_mass, 1, tolerance = 1e-7)
    expect_equal(projection_mass, 1, tolerance = 1e-7)
  }
})

test_that("H2 mean projection density is the shifted exponential density", {
  y <- seq(1, 5, length.out = 100L)
  kappa <- 2.5

  expect_equal(
    hvmf_mean_projection_density(y, q = 2, kappa = kappa),
    kappa * exp(-kappa * (y - 1)),
    tolerance = 1e-12
  )
})

test_that("regularisation and interpolation stay within 0.015 of the true radial CDF", {
  configurations <- list(
    list(q = 2L, kappa = 2),
    list(q = 2L, kappa = 3),
    list(q = 10L, kappa = 10),
    list(q = 10L, kappa = 15)
  )

  distances <- vapply(configurations, function(configuration) {
    hvmf_regularised_population_distance(
      q = configuration$q,
      kappa = configuration$kappa
    )
  }, numeric(1))

  expect_true(all(distances <= 0.015), info = paste(distances, collapse = ", "))
})

test_that("radial quantile builder expands an insufficient initial bracket", {
  table <- hvmf_build_radial_quantile_table(
    q = 2,
    kappa = 2,
    chi = asinh(1),
    upper = 0.5
  )

  expect_gt(attr(table, "upper"), 0.5)
  expect_equal(tail(table$base, 1L), 0.999)
  expect_equal(nrow(table), 101L)
})

test_that("conditional vMF first moment is correct in H^10", {
  q <- 10L
  kappa <- 15
  mu <- hvmf_section6_mu(q)
  mean_direction <- c(1, rep.int(0, q - 1L))

  set.seed(2401)
  x <- rhvmf_polar(1800, mu = mu, kappa = kappa)
  u <- acosh(x[, 1L])
  angular <- hvmf_spatial_direction(x)
  lambda <- kappa * sinh(asinh(1)) * sinh(u)
  expected_resultant <- vmf_mean_resultant_ratio(lambda, q)

  expect_lt(
    abs(mean(angular %*% mean_direction) - mean(expected_resultant)),
    0.035
  )
  expect_lt(max(abs(colMeans(angular[, -1L, drop = FALSE]))), 0.035)
})

test_that("H10 angular mixture has the prescribed symmetric first moment", {
  q <- 10L
  kappa <- 15
  delta <- pi / 5
  mu <- hvmf_section6_mu(q)
  mean_direction <- c(1, rep.int(0, q - 1L))
  tangent <- c(0, 1, rep.int(0, q - 2L))

  set.seed(2501)
  x <- rhvmf_angular_mixture(
    2400,
    mu = mu,
    kappa = kappa,
    delta = delta,
    tangent = tangent
  )
  u <- acosh(x[, 1L])
  angular <- hvmf_spatial_direction(x)
  lambda <- kappa * sinh(asinh(1)) * sinh(u)
  expected_resultant <- vmf_mean_resultant_ratio(lambda, q)

  expect_lt(
    abs(mean(angular %*% mean_direction) -
          mean(expected_resultant * cos(delta))),
    0.035
  )
  expect_lt(abs(mean(angular %*% tangent)), 0.035)
  expect_lt(max(abs(colMeans(angular[, -(1:2), drop = FALSE]))), 0.035)
})

test_that("quantile cache is isolated across intrinsic dimensions", {
  set.seed(2601)
  x2 <- rhvmf_polar(8, mu = hvmf_section6_mu(2L), kappa = 10)
  set.seed(2602)
  x10 <- rhvmf_polar(8, mu = hvmf_section6_mu(10L), kappa = 10)

  expect_equal(ncol(x2), 3L)
  expect_equal(ncol(x10), 11L)
  expect_true(all(is.finite(x2)))
  expect_true(all(is.finite(x10)))
})
