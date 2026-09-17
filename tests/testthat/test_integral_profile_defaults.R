library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

integral_profile_unit <- function(x) {
  x <- as.numeric(x)
  x / sqrt(sum(x^2))
}

integral_profile_hpoint <- function(chi, angle = 0) {
  c(
    cosh(chi),
    sinh(chi) * cos(angle),
    sinh(chi) * sin(angle)
  )
}

legacy_vmf_direct_profile <- function(
    omega,
    mu,
    kappa,
    t_values,
    distance_type = "geodesic") {
  omega <- integral_profile_unit(omega)
  mu <- integral_profile_unit(mu)
  rho <- sum(mu * omega)
  q <- length(mu)

  vapply(t_values, function(t) {
    threshold <- if (identical(distance_type, "chordal")) {
      1 - t^2 / 2
    } else {
      cos(t)
    }

    density_T <- function(s) {
      log_c_q <- rotasym::c_vMF(
        p = q,
        kappa = kappa,
        log = TRUE
      )
      log_exp_term <- kappa * rho * s
      log_power_term <- ((q - 3) / 2) * log(1 - s^2)
      log_numerator <- log_c_q + log_exp_term + log_power_term
      kappa_term <- kappa * sqrt(
        (1 - s^2) * (1 - rho^2)
      )
      log_denominator <- rotasym::c_vMF(
        p = q - 1,
        kappa = kappa_term,
        log = TRUE
      )

      exp(log_numerator - log_denominator)
    }

    1 - stats::integrate(
      density_T,
      lower = -1 + 1e-8,
      upper = threshold,
      rel.tol = 1e-8,
      abs.tol = 1e-10
    )$value
  }, numeric(1))
}

test_that("shared vMF log normalizer agrees with rotasym", {
  for (q in c(2L, 5L)) {
    for (kappa in c(0.5, 2, 10)) {
      expect_equal(
        vmf_log_normalizing_constant_intrinsic(q, kappa),
        rotasym::c_vMF(p = q + 1L, kappa = kappa, log = TRUE),
        tolerance = 2e-12
      )
    }
  }
})

test_that("general vMF direct integral preserves the old S5 values", {
  mu <- integral_profile_unit(c(1, -0.4, 0.7, 0.2, -0.5, 0.9))
  omega <- integral_profile_unit(c(-0.3, 0.8, 0.1, -0.6, 0.4, 1.1))
  radii <- c(0.35, 0.9, 1.4, 2.1)

  for (kappa in c(2, 10)) {
    legacy <- legacy_vmf_direct_profile(
      omega = omega,
      mu = mu,
      kappa = kappa,
      t_values = radii,
      distance_type = "geodesic"
    )
    current <- theoretical_distance_profile_vmf(
      omega = omega,
      mu = mu,
      kappa = kappa,
      t_values = radii,
      distance_type = "geodesic"
    )

    expect_equal(current, legacy, tolerance = 2e-9)
  }
})

test_that("vMF S2 defaults to integral while tabulation remains explicit", {
  x <- rbind(
    integral_profile_unit(c(1, 0.2, -0.1)),
    integral_profile_unit(c(0.5, 0.8, 0.1)),
    integral_profile_unit(c(-0.2, 0.4, 1))
  )
  theta <- list(
    mu = integral_profile_unit(c(1, 0.3, -0.2)),
    kappa = 2
  )
  spec <- make_vmf_spec(distance_type = "geodesic")
  distances <- spec$distance_matrix(x, x, control = list())

  expect_null(spec$sample_profile_matrix_eval(
    data = x,
    distance_matrix = distances,
    theta = theta,
    control = list()
  ))

  expect_true(is.matrix(spec$sample_profile_matrix_eval(
    data = x,
    distance_matrix = distances,
    theta = theta,
    control = list(
      vmf_profile_method = "tabulated",
      vmf_profile_n_u = 257L
    )
  )))

  radii <- c(0.3, 0.8, 1.5)
  expect_equal(
    spec$profile_eval(x[1, ], radii, theta, control = list()),
    spec$profile_eval(
      x[1, ],
      radii,
      theta,
      control = list(vmf_profile_method = "integral")
    ),
    tolerance = 0
  )
})

test_that("HvMF H2 defaults to integral while tabulation remains explicit", {
  x <- rbind(
    integral_profile_hpoint(0.5, -0.2),
    integral_profile_hpoint(0.8, 0.4),
    integral_profile_hpoint(1.0, -0.5)
  )
  theta <- list(
    mu = integral_profile_hpoint(0.7, 0.1),
    kappa = 4
  )
  spec <- make_hvmf_spec(unknown_param = "both")
  distances <- spec$distance_matrix(x, x, control = list())

  expect_null(spec$sample_profile_matrix_eval(
    data = x,
    distance_matrix = distances,
    theta = theta,
    control = list()
  ))

  expect_true(is.matrix(spec$sample_profile_matrix_eval(
    data = x,
    distance_matrix = distances,
    theta = theta,
    control = list(
      hvmf_profile_method = "tabulated",
      hvmf_profile_n_y = 257L
    )
  )))

  radii <- c(0.2, 0.7, 1.2)
  default_profile <- spec$profile_eval(
    x[1, ],
    radii,
    theta,
    control = list()
  )
  explicit_integral <- spec$profile_eval(
    x[1, ],
    radii,
    theta,
    control = list(hvmf_profile_method = "integral")
  )

  expect_equal(default_profile, explicit_integral, tolerance = 0)
})
