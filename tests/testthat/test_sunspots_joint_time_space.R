library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"))

joint_test_eta <- function() {
  sunspots_joint_time_canonicalize_eta(list(
    weight1 = 0.42, alpha1 = 3.5, beta1 = 8.5, alpha2 = 8, beta2 = 3.2
  ))
}

joint_test_theta <- function(shared = FALSE) {
  if (isTRUE(shared)) {
    return(list(a_N = 0.56, b_N = -0.18, a_S = 0.56, b_S = -0.18, c = 16))
  }
  list(a_N = 0.59, b_N = -0.21, a_S = 0.51, b_S = -0.15, c = 16)
}

joint_test_fit <- function(shared = FALSE) {
  list(eta_hat = joint_test_eta(), theta_hat = joint_test_theta(shared))
}

test_that("dequantized first-record days are reproducible and remain inside the fixed window", {
  temporary_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temporary_csv), add = TRUE)
  dates <- as.POSIXct(c(
    "1997-05-31", "1997-06-01", "1997-06-01", "1998-01-15", "2005-12-31", "2006-01-01"
  ), tz = "UTC")
  z <- c(0.1, 0.2, 0.3, -0.2, -0.3, 0.4)
  example_data <- data.frame(
    cycle = 23L, date = format(dates, tz = "UTC", usetz = TRUE), NOAA = seq_along(dates),
    x1 = 0, x2 = sqrt(1 - z^2), x3 = z
  )
  utils::write.csv(example_data, temporary_csv, row.names = FALSE)

  first <- prepare_sunspots_cycle23_joint_time_space_data(
    temporary_csv, "1997-06-01", "2006-01-01", dequantization_seed = 11L
  )
  repeated <- prepare_sunspots_cycle23_joint_time_space_data(
    temporary_csv, "1997-06-01", "2006-01-01", dequantization_seed = 11L
  )
  changed <- prepare_sunspots_cycle23_joint_time_space_data(
    temporary_csv, "1997-06-01", "2006-01-01", dequantization_seed = 12L
  )

  expect_equal(nrow(first), 4L)
  expect_equal(first$recorded_timestamp, repeated$recorded_timestamp)
  expect_equal(first$dequantization_jitter_day, repeated$dequantization_jitter_day)
  expect_false(isTRUE(all.equal(first$dequantization_jitter_day, changed$dequantization_jitter_day)))
  expect_true(all(first$s > 0 & first$s < 1))
  expect_true(all(first$calendar_day >= as.Date("1997-06-01")))
  expect_true(all(first$calendar_day < as.Date("2006-01-01")))
  expect_equal(first$date, example_data$date[match(first$NOAA, example_data$NOAA)])
  same_day <- which(first$calendar_day == as.Date("1997-06-01"))
  expect_gt(abs(as.numeric(diff(first$dequantized_timestamp[same_day]))), 0)
})

test_that("two-beta time law is normalized, canonically ordered, and has analytic scores", {
  eta <- joint_test_eta()
  quadrature <- sunspots_joint_time_quadrature(eta, n_nodes = 64L)
  expect_equal(sum(quadrature$weights), 1, tolerance = 1e-12)
  expect_lt(quadrature$mass_error, 1e-12)
  expect_lte(eta$mean1, eta$mean2)

  high_precision <- integrate(
    function(s) exp(s) * sunspots_joint_time_density(s, eta), lower = 0, upper = 1,
    subdivisions = 2000L, rel.tol = 1e-11
  )$value
  expect_equal(sum(quadrature$weights * exp(quadrature$nodes)), high_precision, tolerance = 1e-10)

  set.seed(330)
  s <- sample_sunspots_joint_time_beta_mixture2(60L, eta)
  par <- sunspots_joint_time_pack_eta(eta)
  analytic <- sunspots_joint_time_score_matrix(s, par)
  step <- 1e-6
  numeric <- vapply(seq_along(par), function(index) {
    plus <- par
    minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (sunspots_joint_time_log_density(s, sunspots_joint_time_unpack_eta(plus)) -
       sunspots_joint_time_log_density(s, sunspots_joint_time_unpack_eta(minus))) / (2 * step)
  }, numeric(length(s)))
  expect_equal(as.numeric(analytic), as.numeric(numeric), tolerance = 1e-5)

  fit <- suppressWarnings(fit_sunspots_joint_time_beta_mixture2(
    sample_sunspots_joint_time_beta_mixture2(120L, eta),
    control = list(time_beta_n_starts = 3L, time_beta_nelder_mead_control = list(maxit = 1000L, reltol = 1e-8))
  ))
  expect_true(is.finite(fit$loglik))
  expect_gte(fit$n_successful_starts, 1L)
})

test_that("the product metric and joint profile have the required endpoint and reduction properties", {
  omega <- c(0.3, -0.4, sqrt(0.75))
  expect_equal(sunspots_joint_distance(rbind(omega, omega), c(0.4, 0.4), omega, 0.4), c(0, 0), tolerance = 1e-14)
  expect_equal(sunspots_joint_distance(rbind(c(0, 0, -1)), 0, c(0, 0, 1), 1), 1, tolerance = 1e-14)

  fit <- joint_test_fit(shared = TRUE)
  quadrature <- sunspots_joint_time_quadrature(fit$eta_hat, n_nodes = 24L)
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat, quadrature$nodes, l_max = 50L, quad_n = 200L
  )
  radii <- matrix(seq(0, 1, length.out = 101L), nrow = 1L)
  profile <- sunspots_joint_profile_block_r(
    radii, rho = omega[[3L]], center_s = 0.5, time_nodes = quadrature$nodes,
    time_weights = quadrature$weights, coefficients = coefficients
  )
  expect_true(all(profile >= 0 & profile <= 1))
  expect_true(all(diff(profile[1, ]) >= -1e-10))
  expect_equal(profile[1, 1], 0, tolerance = 1e-12)
  expect_equal(profile[1, 101], 1, tolerance = 1e-12)

  fixed_time <- 0.5
  fixed_coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat, fixed_time, l_max = 80L, quad_n = 300L
  )
  product_radii <- matrix(seq(0, 0.5, length.out = 31L), nrow = 1L)
  product_profile <- sunspots_joint_profile_block_r(
    product_radii, rho = omega[[3L]], center_s = fixed_time, time_nodes = fixed_time,
    time_weights = 1, coefficients = fixed_coefficients
  )
  expected <- distance_profile_small_circle_symmetric_mixture2(
    omega = omega, t_values = 2 * pi * product_radii[1, ], mu = c(0, 0, 1),
    kappa = fit$theta_hat$c,
    nu = fit$theta_hat$a_N + fixed_time * fit$theta_hat$b_N,
    distance_type = "geodesic", method = "legendre", l_max = 80L, quad_n = 300L
  )
  expect_equal(as.numeric(product_profile), expected, tolerance = 1e-10)
})

test_that("the compiled joint profile agrees with the R reference by blocks", {
  skip_if_not_installed("Rcpp")
  eta <- joint_test_eta()
  theta <- joint_test_theta()
  quadrature <- sunspots_joint_time_quadrature(eta, n_nodes = 18L)
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    theta, quadrature$nodes, l_max = 40L, quad_n = 180L
  )
  set.seed(331)
  radii <- matrix(runif(35L), nrow = 5L)
  rho <- runif(5L, -1, 1)
  center_s <- runif(5L)
  reference <- sunspots_joint_profile_block(
    radii, rho, center_s, quadrature$nodes, quadrature$weights, coefficients, backend = "r"
  )
  compiled <- sunspots_joint_profile_block(
    radii, rho, center_s, quadrature$nodes, quadrature$weights, coefficients, backend = "cpp"
  )
  expect_equal(compiled, reference, tolerance = 1e-12)
})

test_that("joint scores and fast KS/CvM work for both parameter dimensions", {
  for (shared in c(FALSE, TRUE)) {
    control <- if (shared) list(hemisphere_regression = "shared") else list()
    fit <- joint_test_fit(shared)
    par <- sunspots_joint_pack_par(fit, control)
    set.seed(if (shared) 333 else 332)
    simulated <- sample_sunspots_joint_time_space(40L, par, control)
    analytic <- sunspots_joint_score_matrix(simulated, par, control)
    step <- 1e-6
    numeric <- vapply(seq_along(par), function(index) {
      plus <- par
      minus <- par
      plus[[index]] <- plus[[index]] + step
      minus[[index]] <- minus[[index]] - step
      plus_state <- sunspots_joint_state_from_par(plus, control)
      minus_state <- sunspots_joint_state_from_par(minus, control)
      plus_density <- sunspots_joint_time_log_density(simulated$s, plus_state$eta) +
        sunspots_time_varying_log_density(simulated$x, simulated$s, plus_state$theta)
      minus_density <- sunspots_joint_time_log_density(simulated$s, minus_state$eta) +
        sunspots_time_varying_log_density(simulated$x, simulated$s, minus_state$theta)
      (plus_density - minus_density) / (2 * step)
    }, numeric(nrow(simulated$x)))
    expect_equal(as.numeric(analytic), as.numeric(numeric), tolerance = 1e-5)

    prepared <- sunspots_joint_prepare_centers(
      simulated, fit, center_indices = 1:6, time_quad_n = 12L,
      l_max = 30L, spatial_quad_n = 120L, distance_profile_backend = "r"
    )
    fast <- sunspots_joint_prepare_fast_corrections(
      simulated, fit, prepared$centers, derivative_mc_size = 300L,
      seed = if (shared) 335 else 334, control = control
    )
    statistics <- sunspots_time_gof_fast_statistics(
      fast$score_observed, fast$centers, statistics = c("ks", "cvm"),
      B = 5L, seed = if (shared) 337 else 336, n_cores = 1L, bootstrap_block_size = 2L
    )
    dimension <- if (shared) 8L else 10L
    expect_identical(dim(fast$score_observed), c(40L, dimension))
    expect_identical(dim(fast$vhat), c(dimension, dimension))
    expect_true(all(is.finite(statistics$ks)))
    expect_true(all(is.finite(statistics$cvm)))
  }
})

test_that("the slow parametric reference refits the joint model on small samples", {
  control <- list(time_beta_n_starts = 2L,
                  time_beta_nelder_mead_control = list(maxit = 500L, reltol = 1e-7),
                  optim_control = list(maxit = 100L, reltol = 1e-7))
  fit <- joint_test_fit()
  values <- sunspots_joint_slow_reestimated_statistics(
    n = 32L, par = sunspots_joint_pack_par(fit), center_indices = 1:4,
    B = 2L, seed = 338L, statistics = c("ks", "cvm"), time_quad_n = 10L,
    l_max = 24L, spatial_quad_n = 80L, control = control
  )
  expect_identical(names(values), c("ks", "cvm"))
  expect_true(all(is.finite(as.matrix(values))))
  expect_identical(dim(attr(values, "temporal_boundary_flags")), c(2L, 3L))
})
