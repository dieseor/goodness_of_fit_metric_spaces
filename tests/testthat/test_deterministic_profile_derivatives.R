library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

det_unit <- function(x) x / sqrt(sum(x^2))

det_hpoint <- function(rho, angle = 0, q = 2L) {
  spatial <- rep.int(0, q)
  spatial[[1L]] <- sinh(rho) * cos(angle)
  if (q >= 2L) spatial[[2L]] <- sinh(rho) * sin(angle)
  c(cosh(rho), spatial)
}

det_exact_hvmf_h2_sample <- function(n, mu, kappa) {
  z <- stats::rexp(n, rate = kappa)
  angle <- stats::runif(n, min = -pi, max = pi)
  base <- cbind(
    1 + z,
    sqrt(z * (z + 2)) * cos(angle),
    sqrt(z * (z + 2)) * sin(angle)
  )
  spatial <- mu[-1L]
  spatial_norm <- sqrt(sum(spatial^2))
  if (spatial_norm <= sqrt(.Machine$double.eps)) return(base)
  direction <- spatial / spatial_norm
  boost <- matrix(0, 3L, 3L)
  boost[1L, 1L] <- mu[[1L]]
  boost[1L, -1L] <- spatial
  boost[-1L, 1L] <- spatial
  boost[-1L, -1L] <- diag(2L) +
    (mu[[1L]] - 1) * tcrossprod(direction)
  base %*% t(boost)
}

det_vmf_strict_F <- function(xi, omega, t, distance_type = "geodesic") {
  threshold <- if (identical(distance_type, "geodesic")) cos(t) else 1 - t^2 / 2
  if (threshold <= -1) return(1)
  if (threshold >= 1) return(0)
  stats::integrate(
    function(s) vmf_projected_density_canonical(s, xi, omega)$density,
    lower = threshold,
    upper = 1,
    rel.tol = 1e-11,
    abs.tol = 1e-12,
    subdivisions = 2000L
  )$value
}

det_hvmf_strict_F <- function(xi, omega, t) {
  kappa <- sqrt(-hvmf_minkowski_inner_product(xi, xi))
  a <- hvmf_minkowski_inner_product(xi, omega)
  b_sq <- profile_derivative_nonnegative_square(
    a^2 - kappa^2,
    max(a^2, kappa^2),
    "test HvMF b^2"
  )
  chi <- asinh(sqrt(b_sq) / kappa)
  stats::integrate(
    function(r) hvmf_radial_density(
      r,
      q = length(xi) - 1L,
      kappa = kappa,
      chi = chi
    ),
    lower = 0,
    upper = t,
    rel.tol = 2e-10,
    abs.tol = 1e-12,
    subdivisions = 2000L
  )$value
}

det_componentwise_difference <- function(fun, x, relative_step = 1e-5) {
  vapply(seq_along(x), function(j) {
    step <- relative_step * max(1, abs(x[[j]]))
    plus <- minus <- x
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step
    (fun(plus) - fun(minus)) / (2 * step)
  }, numeric(1))
}

test_that("vMF deterministic canonical derivative matches strict finite differences", {
  xi <- c(2.4, -0.7, 1.1)
  omega <- det_unit(c(0.2, 0.8, -0.3))
  for (configuration in list(
    list(t = 1.2, distance = "geodesic"),
    list(t = 1.0, distance = "chordal")
  )) {
    analytic <- vmf_profile_and_derivative_xi(
      omega,
      xi,
      configuration$t,
      configuration$distance,
      grid_size = 16385L
    )$derivative[1L, ]
    numerical <- det_componentwise_difference(
      function(value) det_vmf_strict_F(
        value,
        omega,
        configuration$t,
        configuration$distance
      ),
      xi
    )
    expect_equal(analytic, numerical, tolerance = 3e-7)
  }
})

test_that("HvMF deterministic canonical derivative matches strict finite differences", {
  mu <- det_hpoint(0.6, 0.4)
  xi <- 8 * mu
  omega <- det_hpoint(0.9, -0.3)
  analytic <- hvmf_profile_and_derivative_xi(
    omega,
    xi,
    1.1,
    grid_size = 16385L
  )$derivative[1L, ]
  numerical <- det_componentwise_difference(
    function(value) det_hvmf_strict_F(value, omega, 1.1),
    xi,
    relative_step = 2e-6
  )
  expect_equal(analytic, numerical, tolerance = 4e-7)
})

test_that("HvMF code score is the Euclidean-coordinate score for canonical xi", {
  q <- 2L
  xi <- 6 * det_hpoint(0.55, 0.3)
  x <- det_hpoint(0.8, -0.2)
  log_density <- function(value) {
    kappa <- sqrt(-hvmf_minkowski_inner_product(value, value))
    hvmf_log_normalizing_constant(q, kappa) +
      hvmf_minkowski_inner_product(x, value)
  }
  numerical_score <- det_componentwise_difference(
    log_density,
    xi,
    relative_step = 2e-6
  )
  analytic_score <- drop(hvmf_canonical_score_matrix(
    matrix(x, nrow = 1L),
    xi
  ))
  expect_equal(analytic_score, numerical_score, tolerance = 2e-8)

  numerical_score_jacobian <- vapply(seq_along(xi), function(j) {
    step <- 2e-6 * max(1, abs(xi[[j]]))
    plus <- minus <- xi
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step
    (
      drop(hvmf_canonical_score_matrix(matrix(x, nrow = 1L), plus)) -
        drop(hvmf_canonical_score_matrix(matrix(x, nrow = 1L), minus))
    ) / (2 * step)
  }, numeric(length(xi)))
  expect_equal(
    hvmf_canonical_information(xi),
    -numerical_score_jacobian,
    tolerance = 2e-8
  )
})

test_that("HvMF canonical score satisfies both information identities", {
  set.seed(7261)
  mu <- det_hpoint(0.45, -0.25)
  kappa <- 20
  xi <- kappa * mu
  sample <- det_exact_hvmf_h2_sample(30000L, mu = mu, kappa = kappa)
  score <- hvmf_canonical_score_matrix(sample, xi)
  information <- hvmf_canonical_information(xi)

  score_mean <- colMeans(score)
  score_mean_se <- sqrt(diag(information) / nrow(score))
  expect_true(all(abs(score_mean) <= 4.5 * score_mean_se))

  score_outer <- crossprod(score) / nrow(score)
  expect_lt(
    sqrt(sum((score_outer - information)^2)) /
      sqrt(sum(information^2)),
    0.08
  )

  score_jacobian <- fast_multiplier_numeric_jacobian(
    function(value) colMeans(hvmf_canonical_score_matrix(sample, value)),
    xi,
    step = pmax(1e-5, abs(xi) * 1e-6)
  )
  expect_equal(-score_jacobian, information, tolerance = 2e-8)

  mean_log_likelihood <- function(value) {
    value_kappa <- sqrt(-hvmf_minkowski_inner_product(value, value))
    hvmf_log_normalizing_constant(2L, value_kappa) +
      mean(sample %*% diag(c(-1, 1, 1)) %*% value)
  }
  p <- length(xi)
  step <- pmax(2e-4, abs(xi) * 2e-5)
  f0 <- mean_log_likelihood(xi)
  hessian <- matrix(0, p, p)
  for (j in seq_len(p)) {
    delta <- rep.int(0, p)
    delta[[j]] <- step[[j]]
    hessian[j, j] <- (
      mean_log_likelihood(xi + delta) - 2 * f0 +
        mean_log_likelihood(xi - delta)
    ) / step[[j]]^2
    if (j < p) {
      for (k in (j + 1L):p) {
        delta_k <- rep.int(0, p)
        delta_k[[k]] <- step[[k]]
        hessian[j, k] <- (
          mean_log_likelihood(xi + delta + delta_k) -
            mean_log_likelihood(xi + delta - delta_k) -
            mean_log_likelihood(xi - delta + delta_k) +
            mean_log_likelihood(xi - delta - delta_k)
        ) / (4 * step[[j]] * step[[k]])
        hessian[k, j] <- hessian[j, k]
      }
    }
  }
  expect_equal(-hessian, information, tolerance = 2e-5)
})

test_that("deterministic derivatives equal the integrated conditional score", {
  xi <- c(2.1, -0.4, 0.8)
  omega <- det_unit(c(0.3, 0.7, -0.2))
  t <- 1.15
  q <- length(xi) - 1L
  kappa <- sqrt(sum(xi^2))
  mu <- xi / kappa
  a <- sum(xi * omega)
  b <- sqrt(max(kappa^2 - a^2, 0))
  threshold <- cos(t)
  score_integral <- vapply(seq_along(xi), function(j) {
    stats::integrate(function(s) {
      projected <- vmf_projected_density_canonical(s, xi, omega)
      C <- projected$one_minus_s2 *
        profile_derivative_bessel_ratio_over_argument(projected$u, q)
      conditional_mean_j <- (s - a * C) * omega[[j]] + C * xi[[j]]
      projected$density *
        (conditional_mean_j - A_q(kappa, q) * mu[[j]])
    }, threshold, 1, rel.tol = 1e-10, abs.tol = 1e-12)$value
  }, numeric(1))
  deterministic <- vmf_profile_and_derivative_xi(
    omega,
    xi,
    t,
    "geodesic",
    16385L
  )$derivative[1L, ]
  expect_equal(deterministic, score_integral, tolerance = 3e-7)

  h_mu <- det_hpoint(0.55, 0.2)
  h_xi <- 7 * h_mu
  h_omega <- det_hpoint(0.8, -0.4)
  h_t <- 1.05
  h_q <- length(h_xi) - 1L
  h_kappa <- 7
  h_a <- hvmf_minkowski_inner_product(h_xi, h_omega)
  h_b <- sqrt(max(h_a^2 - h_kappa^2, 0))
  B_value <- hvmf_mean_resultant_ratio(h_q, h_kappa)
  J <- diag(c(-1, rep.int(1, h_q)))
  h_score_integral <- vapply(seq_along(h_xi), function(j) {
    stats::integrate(function(r) {
      chi <- asinh(h_b / h_kappa)
      density <- hvmf_radial_density(r, h_q, h_kappa, chi)
      u <- h_b * sinh(r)
      D <- sinh(r)^2 *
        profile_derivative_bessel_ratio_over_argument(u, h_q)
      conditional_mean_j <-
        (cosh(r) + h_a * D) * h_omega[[j]] + D * h_xi[[j]]
      sign_j <- if (j == 1L) -1 else 1
      density * sign_j * (conditional_mean_j - B_value * h_mu[[j]])
    }, 0, h_t, rel.tol = 1e-10, abs.tol = 1e-12)$value
  }, numeric(1))
  h_deterministic <- hvmf_profile_and_derivative_xi(
    h_omega,
    h_xi,
    h_t,
    16385L
  )$derivative[1L, ]
  expect_equal(h_deterministic, h_score_integral, tolerance = 4e-7)
})

test_that("vMF b=0, parallel and antiparallel centers are stable", {
  mu <- det_unit(c(1, 2, -1))
  xi <- 5 * mu
  for (omega in list(mu, -mu)) {
    result <- vmf_profile_and_derivative_xi(
      omega,
      xi,
      c(0.4, 1.3),
      "geodesic",
      8193L
    )
    expect_equal(result$table$b, 0, tolerance = 2e-7)
    expect_true(all(is.finite(result$derivative)))
    numerical <- det_componentwise_difference(
      function(value) det_vmf_strict_F(value, omega, 1.3),
      xi
    )
    expect_equal(result$derivative[2L, ], numerical, tolerance = 8e-7)
  }
  expect_equal(
    vmf_profile_and_derivative_xi(mu, xi, c(pi, pi + 1), "geodesic")$derivative,
    matrix(0, 2L, 3L)
  )
  expect_equal(
    vmf_profile_and_derivative_xi(mu, xi, c(2, 2.1), "chordal")$derivative,
    matrix(0, 2L, 3L)
  )
})

test_that("H2 omega=mu has the exact radial profile and radial derivative", {
  kappa <- 5
  mu <- det_hpoint(0.7, -0.2)
  xi <- kappa * mu
  t <- c(0.2, 0.7, 1.2)
  result <- hvmf_profile_and_derivative_xi(mu, xi, t, grid_size = 16385L)
  h <- cosh(t) - 1
  exact_F <- 1 - exp(-kappa * h)
  exact_radial_derivative <- h * exp(-kappa * h)

  expect_equal(result$F, exact_F, tolerance = 8e-7)
  expect_equal(
    drop(result$derivative %*% mu),
    exact_radial_derivative,
    tolerance = 2e-6
  )
  J_mu <- c(-mu[[1L]], mu[-1L])
  orthogonal_residual <- result$derivative -
    outer(drop(result$derivative %*% mu) / sum(J_mu * mu), J_mu)
  expect_lt(max(abs(orthogonal_residual)), 2e-8)
})

test_that("HvMF omega=mu and very small rho use the continuous b=0 limit", {
  mu <- det_hpoint(0.45, 0.1)
  xi <- 20 * mu
  exact <- hvmf_profile_and_derivative_xi(mu, xi, c(0.3, 0.8), 8193L)
  near <- hvmf_profile_and_derivative_xi(
    det_hpoint(0.45 + 1e-8, 0.1),
    xi,
    c(0.3, 0.8),
    8193L
  )
  expect_true(all(is.finite(exact$derivative)))
  expect_true(all(is.finite(near$derivative)))
  expect_lt(max(abs(exact$F - near$F)), 2e-6)
  expect_lt(max(abs(exact$derivative - near$derivative)), 3e-5)
})

test_that("projected densities integrate to one and complete derivatives vanish", {
  xi <- c(2.2, -0.5, 0.9)
  omega <- det_unit(c(0.1, 0.8, 0.4))
  vmf_mass <- stats::integrate(
    function(s) vmf_projected_density_canonical(s, xi, omega)$density,
    -1,
    1,
    rel.tol = 1e-10,
    abs.tol = 1e-12
  )$value
  expect_equal(vmf_mass, 1, tolerance = 2e-9)
  expect_equal(
    unname(vmf_profile_and_derivative_xi(
      omega,
      xi,
      pi,
      "geodesic"
    )$derivative),
    matrix(0, 1L, 3L)
  )

  h_mu <- det_hpoint(0.6, 0.3)
  h_omega <- det_hpoint(0.9, -0.1)
  h_kappa <- 12
  chi <- acosh(-hvmf_minkowski_inner_product(h_mu, h_omega))
  h_mass <- stats::integrate(
    function(r) hvmf_radial_density(r, 2L, h_kappa, chi),
    0,
    Inf,
    rel.tol = 2e-9,
    abs.tol = 1e-11,
    subdivisions = 2000L
  )$value
  expect_equal(h_mass, 1, tolerance = 2e-7)
  expect_equal(
    unname(hvmf_profile_and_derivative_xi(
      h_omega,
      h_kappa * h_mu,
      Inf
    )$derivative),
    matrix(0, 1L, 3L)
  )
})

test_that("continuous A(u)/u extension is accurate at and near zero", {
  for (q in c(2L, 3L, 9L, 10L)) {
    u <- c(0, 1e-10, 1e-6, 1e-3, 1e-2)
    values <- profile_derivative_bessel_ratio_over_argument(u, q)
    expect_equal(values[[1L]], 1 / q, tolerance = 0)
    expect_lt(abs(values[[2L]] - 1 / q), 1e-14)
    direct <- besselI(1e-2, q / 2, expon.scaled = TRUE) /
      besselI(1e-2, q / 2 - 1, expon.scaled = TRUE) / 1e-2
    expect_equal(values[[5L]], direct, tolerance = 1e-12)
    expect_lt(max(abs(diff(values[1:4]))), 1e-7)
  }
})

test_that("vMF derivatives are rotation equivariant", {
  angle <- 0.7
  R <- matrix(c(
    cos(angle), -sin(angle), 0,
    sin(angle), cos(angle), 0,
    0, 0, 1
  ), 3L, 3L, byrow = TRUE)
  xi <- c(2.3, -0.4, 0.8)
  omega <- det_unit(c(0.2, 0.7, -0.1))
  original <- vmf_profile_and_derivative_xi(omega, xi, 1.1, "geodesic", 8193L)
  rotated <- vmf_profile_and_derivative_xi(
    drop(R %*% omega),
    drop(R %*% xi),
    1.1,
    "geodesic",
    8193L
  )
  expect_equal(rotated$F, original$F, tolerance = 2e-10)
  expect_equal(
    rotated$derivative[1L, ],
    drop(R %*% original$derivative[1L, ]),
    tolerance = 4e-8
  )
})

test_that("HvMF derivatives are Lorentz equivariant", {
  phi <- 0.35
  L <- diag(3L)
  L[1:2, 1:2] <- matrix(c(
    cosh(phi), sinh(phi),
    sinh(phi), cosh(phi)
  ), 2L, 2L, byrow = TRUE)
  J <- diag(c(-1, 1, 1))
  expect_equal(t(L) %*% J %*% L, J, tolerance = 2e-14)

  mu <- det_hpoint(0.55, 0.4)
  omega <- det_hpoint(0.8, -0.2)
  xi <- 9 * mu
  original <- hvmf_profile_and_derivative_xi(omega, xi, 1.0, 8193L)
  transformed <- hvmf_profile_and_derivative_xi(
    drop(L %*% omega),
    drop(L %*% xi),
    1.0,
    8193L
  )
  expect_equal(transformed$F, original$F, tolerance = 3e-9)
  expect_equal(
    transformed$derivative[1L, ],
    drop(solve(t(L), original$derivative[1L, ])),
    tolerance = 7e-8
  )
})

test_that("paper-scale concentrations and dimensions remain finite", {
  vmf_cases <- list(
    list(q = 2L, kappa = 0.5),
    list(q = 2L, kappa = 5),
    list(q = 9L, kappa = 10)
  )
  for (case in vmf_cases) {
    mu <- c(1, rep.int(0, case$q))
    omega <- det_unit(c(0.8, 0.6, rep.int(0, case$q - 1L)))
    result <- vmf_profile_and_derivative_xi(
      omega,
      case$kappa * mu,
      c(0.3, 0.8, 1.4),
      "geodesic",
      4097L
    )
    expect_true(all(is.finite(c(result$F, result$derivative))))
  }

  for (case in list(
    list(q = 2L, kappa = 50),
    list(q = 2L, kappa = 200),
    list(q = 10L, kappa = 10)
  )) {
    mu <- c(1, rep.int(0, case$q))
    omega <- det_hpoint(0.6, 0.2, q = case$q)
    result <- hvmf_profile_and_derivative_xi(
      omega,
      case$kappa * mu,
      c(0.3, 0.7, 1.1),
      4097L
    )
    expect_true(all(is.finite(c(result$F, result$derivative))))
    expect_true(all(result$F >= 0 & result$F <= 1))
  }
})

test_that("fast multiplier keeps score_mc and supports quadrature end to end", {
  set.seed(731)
  x <- rotasym::r_vMF(14, mu = c(1, 0, 0), kappa = 3)
  default_prep <- prepare_vmf_fast_multiplier(
    data = x,
    theta_hat = normalize_vmf_theta(compute_mle_xi(x)),
    spec = make_vmf_spec("geodesic", "xi"),
    control = list(),
    distance_type = "geodesic"
  )
  expect_identical(default_prep$derivative_method, "quadrature")

  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 3L,
    seed = 732L,
    bootstrap_method = "fast_multiplier",
    distance_type = "geodesic",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    )
  )
  deterministic <- do.call(multiplier_bootstrap_vmf, c(
    common,
    list(control = list(
      derivative_method = "quadrature",
      vmf_profile_n_u = 513L
    ))
  ))
  monte_carlo <- do.call(multiplier_bootstrap_vmf, c(
    common,
    list(control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 200L,
      derivative_mc_seed = 733L
    ))
  ))
  expect_identical(deterministic$diagnostics$derivative_method, "quadrature")
  expect_identical(monte_carlo$diagnostics$derivative_method, "score_mc")
  legacy <- expect_warning(
    do.call(multiplier_bootstrap_vmf, c(
      common,
      list(control = list(
        derivative_mc_size = 200L,
        derivative_mc_seed = 733L
      ))
    )),
    "legacy.*selecting `quadrature`"
  )
  expect_identical(legacy$diagnostics$derivative_method_effective, "quadrature")
  expect_identical(
    legacy$diagnostics$derivative_method_selection_source,
    "model_default_legacy_mc_controls_ignored"
  )
  expect_identical(
    deterministic$diagnostics$fast_ks_mode,
    "sample_points_unique_distances_streamed"
  )
  expect_true(all(is.finite(c(
    deterministic$inference$ks$p_value,
    deterministic$inference$cvm$p_value,
    monte_carlo$inference$ks$p_value,
    monte_carlo$inference$cvm$p_value
  ))))

  h_mu <- det_hpoint(0.5, 0.1)
  h_data <- rhvmf_h2_polar(14, h_mu, 8)
  h_default_prep <- prepare_hvmf_fast_multiplier(
    spec = make_hvmf_spec(),
    data = h_data,
    theta_hat = hvmf_mle_h2(h_data),
    control = list()
  )
  expect_identical(h_default_prep$derivative_method, "quadrature")
  h_result <- multiplier_bootstrap_hvmf(
    data = h_data,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 3L,
    seed = 734L,
    bootstrap_method = "fast_multiplier",
    fast_multiplier_backend = "r",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "quadrature",
      hvmf_profile_n_y = 513L
    )
  )
  expect_identical(h_result$diagnostics$derivative_method, "quadrature")
  h_legacy <- expect_warning(
    multiplier_bootstrap_hvmf(
      data = h_data,
      null = list(type = "composite"),
      statistics = "ks",
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = 3L,
      seed = 735L,
      bootstrap_method = "fast_multiplier",
      fast_multiplier_backend = "r",
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = list(
        derivative_mc_size = 200L,
        derivative_mc_seed = 736L
      )
    ),
    "legacy.*selecting `quadrature`"
  )
  expect_identical(h_legacy$diagnostics$derivative_method_effective, "quadrature")
  expect_identical(
    h_legacy$diagnostics$derivative_method_selection_source,
    "model_default_legacy_mc_controls_ignored"
  )
  expect_identical(
    h_result$diagnostics$fast_cvm_mode,
    "sample_points_unique_distances_sorted_rows"
  )
})

test_that("HvMF score_mc retains its historical parametrizations exactly", {
  check_case <- function(q, data_seed, derivative_seed) {
    mu <- if (q == 2L) {
      det_hpoint(0.4, 0)
    } else {
      c(sqrt(2), 1, rep.int(0, q - 1L))
    }
    kappa <- if (q == 2L) 3 else 10
    set.seed(data_seed)
    data <- rhvmf_polar(24L, mu = mu, kappa = kappa)
    spec <- make_hvmf_spec(unknown_param = "both")
    theta_hat <- fit_hvmf_theta(
      data,
      null = list(type = "composite")
    )
    prep <- prepare_hvmf_fast_multiplier(
      spec = spec,
      data = data,
      theta_hat = theta_hat,
      control = list(
        derivative_method = "score_mc",
        derivative_mc_size = 80L,
        derivative_mc_seed = derivative_seed
      )
    )

    if (q == 2L) {
      par0 <- c(theta_hat$chi, theta_hat$theta, theta_hat$kappa)
      state_from_par <- function(par) {
        normalize_hvmf_theta(list(
          mu = c(
            cosh(par[[1L]]),
            sinh(par[[1L]]) * cos(par[[2L]]),
            sinh(par[[1L]]) * sin(par[[2L]])
          ),
          kappa = par[[3L]]
        ))
      }
      historical_score <- function(sample, par) {
        state <- state_from_par(par)
        dmu_dchi <- c(
          sinh(state$chi),
          cosh(state$chi) * cos(state$theta),
          cosh(state$chi) * sin(state$theta)
        )
        dmu_dtheta <- c(
          0,
          -sinh(state$chi) * sin(state$theta),
          sinh(state$chi) * cos(state$theta)
        )
        cbind(
          state$kappa * apply(sample, 1L, function(row) {
            minkowski_inner_product(row, dmu_dchi)
          }),
          state$kappa * apply(sample, 1L, function(row) {
            minkowski_inner_product(row, dmu_dtheta)
          }),
          1 / state$kappa + 1 + apply(sample, 1L, function(row) {
            minkowski_inner_product(row, state$mu)
          })
        )
      }
      set.seed(derivative_seed)
      state <- state_from_par(par0)
      auxiliary <- rhvmf_h2_polar(
        80L,
        mu = state$mu,
        kappa = state$kappa
      )
    } else {
      par0 <- c(theta_hat$mu[-1L], log(theta_hat$kappa))
      state_from_par <- function(par) {
        hvmf_hq_theta_from_coordinates(par, q = q)
      }
      historical_score <- function(sample, par) {
        state <- state_from_par(par)
        eta <- state$mu[-1L]
        mu0 <- state$mu[[1L]]
        minkowski_mu <- -sample[, 1L] * mu0 +
          rowSums(
            sample[, -1L, drop = FALSE] *
              rep(eta, each = nrow(sample))
          )
        cbind(
          state$kappa * (
            sample[, -1L, drop = FALSE] -
              tcrossprod(sample[, 1L], eta / mu0)
          ),
          state$kappa * (
            hvmf_mean_resultant_ratio(q, state$kappa) +
              minkowski_mu
          )
        )
      }
      set.seed(derivative_seed)
      state <- state_from_par(par0)
      auxiliary <- rhvmf_polar(
        80L,
        mu = state$mu,
        kappa = state$kappa
      )
    }

    expected_observed_score <- historical_score(data, par0)
    expected_auxiliary_score <- historical_score(auxiliary, par0)
    expect_identical(
      unname(prep$S_obs),
      unname(expected_observed_score)
    )
    expect_identical(
      unname(prep$Psi_aux),
      unname(expected_auxiliary_score)
    )
    expect_identical(
      unname(prep$Vhat),
      unname(
        crossprod(expected_auxiliary_score) /
          nrow(expected_auxiliary_score)
      )
    )
  }

  check_case(q = 2L, data_seed = 88101L, derivative_seed = 88102L)
  check_case(q = 10L, data_seed = 88201L, derivative_seed = 88202L)
})

test_that("the historical HvMF null sampler does not draw unused signs", {
  global_environment <- globalenv()
  had_sample_binding <- exists(
    "sample",
    envir = global_environment,
    inherits = FALSE
  )
  if (had_sample_binding) {
    previous_sample_binding <- get(
      "sample",
      envir = global_environment,
      inherits = FALSE
    )
  }
  on.exit({
    if (had_sample_binding) {
      assign(
        "sample",
        previous_sample_binding,
        envir = global_environment
      )
    } else if (exists(
      "sample",
      envir = global_environment,
      inherits = FALSE
    )) {
      rm("sample", envir = global_environment)
    }
  }, add = TRUE)
  assign(
    "sample",
    function(...) stop("unused sign draw"),
    envir = global_environment
  )

  set.seed(88301L)
  expect_no_error(
    rhvmf_polar(
      8L,
      mu = det_hpoint(0.4, 0),
      kappa = 3
    )
  )
})
