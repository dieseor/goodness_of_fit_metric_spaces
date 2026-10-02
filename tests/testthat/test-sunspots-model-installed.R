sunspots_joint_time_density <- function(s, eta) exp(sunspots_joint_time_log_density(s, eta))

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

test_that("temporal beta-mixture selection only accepts converged fits", {
  fits <- list(
    list(value = 1, convergence = 1L, par = rep(0, 5L), message = "iteration limit"),
    list(value = 2, convergence = 0L, par = rep(0, 5L), message = "converged")
  )

  selected <- sunspots_joint_time_select_fit(fits)

  expect_identical(selected$fit$convergence, 0L)
  expect_equal(selected$fit$value, 2)
  expect_identical(selected$n_finite_fits, 2L)
  expect_identical(selected$n_converged_fits, 1L)
  expect_true(selected$selected_converged)
})

test_that("temporal beta-mixture selection fails if no fit converges", {
  fits <- list(
    list(value = 1, convergence = 1L, par = rep(0, 5L), message = "iteration limit"),
    list(value = 2, convergence = 52L, par = rep(0, 5L), message = "line-search failure")
  )

  expect_error(
    sunspots_joint_time_select_fit(fits),
    "No temporal two-beta-mixture optimization converged"
  )
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

test_that("two-beta temporal formulas support positive shapes below one", {
  eta <- sunspots_joint_time_canonicalize_eta(list(
    weight1 = 0.4,
    alpha1 = 0.7,
    beta1 = 2.5,
    alpha2 = 3.0,
    beta2 = 0.8
  ))

  normalization <- stats::integrate(
    function(s) sunspots_joint_time_density(s, eta),
    lower = 0,
    upper = 1,
    subdivisions = 4000L,
    rel.tol = 1e-10,
    abs.tol = 1e-12
  )$value
  expect_equal(normalization, 1, tolerance = 1e-8)

  quadrature <- sunspots_joint_time_quadrature(eta, n_nodes = 96L)
  expect_equal(sum(quadrature$weights), 1, tolerance = 1e-12)
  expect_lt(quadrature$mass_error, 1e-12)

  set.seed(20260806)
  s <- sample_sunspots_joint_time_beta_mixture2(80L, eta)
  par <- sunspots_joint_time_pack_eta(eta)
  analytic <- sunspots_joint_time_score_matrix(s, par)

  step <- 1e-6
  numeric <- vapply(seq_along(par), function(index) {
    plus <- par
    minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (
      sunspots_joint_time_log_density(
        s,
        sunspots_joint_time_unpack_eta(plus)
      ) -
      sunspots_joint_time_log_density(
        s,
        sunspots_joint_time_unpack_eta(minus)
      )
    ) / (2 * step)
  }, numeric(length(s)))

  expect_equal(
    as.numeric(analytic),
    as.numeric(numeric),
    tolerance = 1e-5
  )
})
