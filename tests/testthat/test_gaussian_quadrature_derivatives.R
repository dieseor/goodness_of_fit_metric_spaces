library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

gaussian_test_reference_F <- function(omega, mu, Sigma, t) {
  eig <- eigen(Sigma, symmetric = TRUE)
  nu <- drop(crossprod(eig$vectors, mu - omega))
  mvnormal_quadform_cdf(
    q = t^2,
    lambda = eig$values,
    h = rep.int(1, length(mu)),
    delta = nu^2 / eig$values,
    control = list(
      mvnormal_quadform_method = "farebrother",
      mvnormal_quadform_eps = 1e-12
    )
  )
}

gaussian_test_theta_to_state <- function(theta, q) {
  list(
    mu = theta[seq_len(q)],
    Sigma = fast_multiplier_ivech(theta[q + seq_len(q * (q + 1L) / 2L)], q)
  )
}

gaussian_test_matrix_gradient <- function(gradient_vech, q) {
  output <- matrix(0, q, q)
  index <- which(lower.tri(output, diag = TRUE), arr.ind = TRUE)
  for (k in seq_len(nrow(index))) {
    i <- index[k, 1L]
    j <- index[k, 2L]
    output[i, j] <- output[j, i] <- gradient_vech[[k]] /
      if (i == j) 1 else 2
  }
  output
}

test_that("univariate Gaussian F and both derivatives match closed forms", {
  mu <- 0.4
  variance <- 1.7
  sigma <- sqrt(variance)
  omega <- -0.2
  thresholds <- c(0, 0.15, 0.8, 2.5, Inf)
  result <- gaussian_ball_profile_quadrature(
    omega, mu, matrix(variance, 1L, 1L), thresholds,
    control = list(gaussian_quadrature_abs_tol = 1e-13)
  )
  finite <- is.finite(thresholds) & thresholds > 0
  upper <- (omega + thresholds[finite] - mu) / sigma
  lower <- (omega - thresholds[finite] - mu) / sigma
  expected_F <- numeric(length(thresholds))
  expected_mu <- expected_variance <- numeric(length(thresholds))
  expected_F[finite] <- pnorm(upper) - pnorm(lower)
  expected_F[is.infinite(thresholds)] <- 1
  expected_mu[finite] <- (dnorm(lower) - dnorm(upper)) / sigma
  expected_variance[finite] <- (
    lower * dnorm(lower) - upper * dnorm(upper)
  ) / (2 * variance)
  expect_equal(result$F, expected_F, tolerance = 2e-12)
  expect_equal(result$gradient_mu[, 1L], expected_mu, tolerance = 2e-12)
  expect_equal(
    result$gradient_vech_sigma[, 1L], expected_variance,
    tolerance = 2e-12
  )
})

test_that("Gaussian dot F matches independent finite differences in vech order", {
  mu <- c(0.35, -0.25)
  Sigma <- matrix(c(1.1, 0.24, 0.24, 0.65), 2L)
  omega <- c(-0.15, 0.3)
  t <- 1.05
  analytic <- gaussian_ball_profile_quadrature(
    omega, mu, Sigma, t,
    control = list(gaussian_quadrature_abs_tol = 1e-6)
  )$derivative[1L, ]
  theta <- c(mu, fast_multiplier_vech(Sigma))
  numerical <- numDeriv::grad(function(value) {
    state <- gaussian_test_theta_to_state(value, 2L)
    gaussian_test_reference_F(omega, state$mu, state$Sigma, t)
  }, theta, method.args = list(eps = 2e-5, d = 1e-4))
  expect_equal(analytic, numerical, tolerance = 2e-7)
})

test_that("q=10 Gaussian dot F matches all independent vech finite differences", {
  q <- 10L
  mu <- seq(-0.2, 0.25, length.out = q)
  Sigma <- toeplitz(0.2^(0:(q - 1L)))
  omega <- rev(mu) / 2
  t <- 3
  analytic <- gaussian_ball_profile_quadrature(
    omega, mu, Sigma, t,
    control = list(gaussian_quadrature_abs_tol = 1e-6)
  )$derivative[1L, ]
  theta <- c(mu, fast_multiplier_vech(Sigma))
  numerical <- vapply(seq_along(theta), function(j) {
    step <- 2e-5 * max(1, abs(theta[[j]]))
    plus <- minus <- theta
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step
    plus_state <- gaussian_test_theta_to_state(plus, q)
    minus_state <- gaussian_test_theta_to_state(minus, q)
    (
      gaussian_test_reference_F(
        omega, plus_state$mu, plus_state$Sigma, t
      ) -
        gaussian_test_reference_F(
          omega, minus_state$mu, minus_state$Sigma, t
        )
    ) / (2 * step)
  }, numeric(1))
  expect_equal(analytic, numerical, tolerance = 8e-7)
})

test_that("Gaussian derivative satisfies score and heat identities", {
  set.seed(8124)
  q <- 2L
  mu <- c(0.3, -0.15)
  Sigma <- matrix(c(0.9, 0.18, 0.18, 0.55), q)
  omega <- c(-0.1, 0.25)
  t <- 0.95
  result <- gaussian_ball_profile_quadrature(
    omega, mu, Sigma, t,
    control = list(gaussian_quadrature_abs_tol = 1e-5)
  )
  n_mc <- 120000L
  x <- mvtnorm::rmvnorm(n_mc, mu, Sigma)
  score <- gaussian_score_matrix_vech(x, mu, Sigma)
  indicator <- as.numeric(sqrt(rowSums((x - rep(omega, each = n_mc))^2)) <= t)
  score_terms <- score * indicator
  estimate <- colMeans(score_terms)
  standard_error <- apply(score_terms, 2L, sd) / sqrt(n_mc)
  expect_true(all(abs(estimate - result$derivative[1L, ]) <=
                    5 * standard_error + 3e-4))

  gradient_matrix <- gaussian_test_matrix_gradient(
    result$gradient_vech_sigma[1L, ], q
  )
  heat_matrix <- 0.5 * result$eigenvectors %*%
    result$hessian_nu[, , 1L] %*% t(result$eigenvectors)
  expect_equal(gradient_matrix, heat_matrix, tolerance = 2e-12)

  numerical_hessian <- numDeriv::hessian(function(value) {
    gaussian_test_reference_F(omega, value, Sigma, t)
  }, mu, method.args = list(eps = 2e-4, d = 1e-3))
  expect_equal(gradient_matrix, 0.5 * numerical_hessian, tolerance = 2e-5)
})

test_that("analytic Gaussian information satisfies all likelihood identities", {
  set.seed(8125)
  q <- 2L
  mu <- c(0.2, -0.3)
  Sigma <- matrix(c(0.8, 0.16, 0.16, 0.5), q)
  information <- gaussian_fisher_information_vech(Sigma)
  n_mc <- 100000L
  x <- mvtnorm::rmvnorm(n_mc, mu, Sigma)
  score <- gaussian_score_matrix_vech(x, mu, Sigma)
  score_outer <- crossprod(score) / n_mc
  expect_lt(
    sqrt(sum((score_outer - information)^2)) /
      sqrt(sum(information^2)),
    0.025
  )

  theta <- c(mu, fast_multiplier_vech(Sigma))
  mean_log_likelihood <- function(value) {
    state <- gaussian_test_theta_to_state(value, q)
    centered <- sweep(x, 2L, state$mu, "-")
    inverse <- solve(state$Sigma)
    -0.5 * determinant(state$Sigma, logarithm = TRUE)$modulus -
      0.5 * mean(rowSums((centered %*% inverse) * centered))
  }
  expected_hessian <- numDeriv::hessian(
    mean_log_likelihood, theta,
    method.args = list(eps = 2e-4, d = 1e-3)
  )
  expect_lt(
    sqrt(sum((-expected_hessian - information)^2)) /
      sqrt(sum(information^2)),
    0.025
  )

  score_jacobian <- numDeriv::jacobian(function(value) {
    state <- gaussian_test_theta_to_state(value, q)
    colMeans(gaussian_score_matrix_vech(x, state$mu, state$Sigma))
  }, theta)
  expect_lt(
    sqrt(sum((-score_jacobian - information)^2)) /
      sqrt(sum(information^2)),
    0.025
  )
})

test_that("Gaussian quadrature is invariant under orthogonal changes and eigensigns", {
  mu <- c(0.4, -0.1, 0.25)
  Sigma <- matrix(c(
    1.0, 0.12, -0.08,
    0.12, 0.7, 0.05,
    -0.08, 0.05, 0.45
  ), 3L)
  omega <- c(-0.2, 0.3, 0.1)
  t <- c(0.4, 1.1, 2.2)
  base <- gaussian_ball_profile_quadrature(
    omega, mu, Sigma, t,
    control = list(gaussian_quadrature_abs_tol = 1e-6)
  )
  Q <- qr.Q(qr(matrix(c(1, 2, -1, 0, 1, 2, 2, -1, 1), 3L)))
  rotated <- gaussian_ball_profile_quadrature(
    drop(Q %*% omega), drop(Q %*% mu), Q %*% Sigma %*% t(Q), t,
    control = list(gaussian_quadrature_abs_tol = 1e-6)
  )
  expect_equal(rotated$F, base$F, tolerance = 2e-11)
  expect_equal(rotated$gradient_mu, base$gradient_mu %*% t(Q), tolerance = 3e-10)
  for (k in seq_along(t)) {
    base_K <- gaussian_test_matrix_gradient(base$gradient_vech_sigma[k, ], 3L)
    rotated_K <- gaussian_test_matrix_gradient(rotated$gradient_vech_sigma[k, ], 3L)
    expect_equal(rotated_K, Q %*% base_K %*% t(Q), tolerance = 4e-10)
  }

  eig <- eigen(Sigma, symmetric = TRUE)
  signed_eig <- eig
  signed_eig$vectors <- sweep(eig$vectors, 2L, c(-1, 1, -1), "*")
  sign_changed <- gaussian_ball_profile_quadrature(
    omega, mu, Sigma, t,
    control = list(gaussian_quadrature_abs_tol = 1e-6),
    spectral = signed_eig
  )
  expect_equal(sign_changed$F, base$F, tolerance = 2e-12)
  expect_equal(sign_changed$derivative, base$derivative, tolerance = 2e-11)
})

test_that("isotropic, repeated, near-repeated and conditioned spectra are stable", {
  cases <- list(
    diag(2L),
    diag(c(1, 1, 0.4)),
    diag(c(1, 1 + 1e-11, 0.4)),
    diag(c(1, 0.2, 0.05))
  )
  for (Sigma in cases) {
    q <- nrow(Sigma)
    mu <- seq_len(q) / (5 * q)
    omega <- rev(mu) / 3
    result <- gaussian_ball_profile_quadrature(
      omega, mu, Sigma, c(0.2, 1, 3),
      control = list(
        gaussian_quadrature_abs_tol = if (q == 2L) 1e-4 else 1e-6,
        gaussian_quadrature_max_terms = 1000000L
      )
    )
    reference <- vapply(c(0.2, 1, 3), function(radius) {
      gaussian_test_reference_F(omega, mu, Sigma, radius)
    }, numeric(1))
    expect_equal(result$F, reference, tolerance = 3e-8)
    expect_true(all(is.finite(result$derivative)))
    expect_lte(
      result$diagnostics$propagated_error_estimate,
      1.01 * if (q == 2L) 1e-4 else 1e-6
    )
  }

  isotropic <- 0.7 * diag(3L)
  mu <- c(0.2, -0.1, 0.3)
  omega <- c(-0.15, 0.25, 0.05)
  canonical <- list(values = rep.int(0.7, 3L), vectors = diag(3L))
  rotation <- qr.Q(qr(matrix(c(1, 2, -1, 0, 1, 2, 2, -1, 1), 3L)))
  alternative <- list(values = rep.int(0.7, 3L), vectors = rotation)
  first <- gaussian_ball_profile_quadrature(
    omega, mu, isotropic, c(0.6, 1.4), spectral = canonical
  )
  second <- gaussian_ball_profile_quadrature(
    omega, mu, isotropic, c(0.6, 1.4), spectral = alternative
  )
  expect_equal(second$F, first$F, tolerance = 2e-12)
  expect_equal(second$derivative, first$derivative, tolerance = 3e-11)
})

test_that("normal and logistic-Gaussian use the same ilr Gaussian engine", {
  set.seed(8126)
  z <- mvtnorm::rmvnorm(
    14L, c(0.25, -0.2), matrix(c(0.7, 0.1, 0.1, 0.45), 2L)
  )
  simplex <- logistic_gaussian_ilr_to_simplex(z, ambient_dim = 3L)
  normal_spec <- make_mvnormal_spec("both")
  lg_spec <- make_logistic_gaussian_spec("both")
  normal_theta <- normal_spec$fit_theta(
    z, NULL, list(type = "composite"), list()
  )
  lg_theta <- lg_spec$fit_theta(
    simplex, NULL, list(type = "composite"), list()
  )
  expect_equal(lg_theta$mu_ilr, normal_theta$mu, tolerance = 1e-13)
  expect_equal(lg_theta$Sigma_ilr, normal_theta$Sigma, tolerance = 1e-13)
  center <- simplex[3L, ]
  center_ilr <- logistic_gaussian_point_to_ilr(center)
  radii <- c(0.3, 0.9, 1.8)
  normal_profile <- gaussian_ball_profile_quadrature(
    center_ilr, normal_theta$mu, normal_theta$Sigma, radii
  )
  lg_profile <- gaussian_ball_profile_quadrature(
    center_ilr, lg_theta$mu_ilr, lg_theta$Sigma_ilr, radii
  )
  expect_equal(lg_profile$F, normal_profile$F, tolerance = 1e-13)
  expect_equal(lg_profile$derivative, normal_profile$derivative, tolerance = 1e-13)
  expect_equal(
    gaussian_score_matrix_vech(z, normal_theta$mu, normal_theta$Sigma),
    gaussian_score_matrix_vech(
      logistic_gaussian_ilr_matrix(simplex), lg_theta$mu_ilr,
      lg_theta$Sigma_ilr
    ),
    tolerance = 1e-13
  )
  expect_equal(
    lg_spec$profile_eval(
      center, radii, lg_theta,
      control = list(derivative_method = "quadrature")
    ),
    normal_profile$F,
    tolerance = 2e-9
  )
  expect_equal(
    gaussian_fisher_information_vech(normal_theta$Sigma),
    gaussian_fisher_information_vech(lg_theta$Sigma_ilr),
    tolerance = 1e-13
  )

  ks_n <- prepare_ks_observed_data(
    z, normal_spec, normal_theta, make_sample_unique_distance_ks_grid(),
    control = list(derivative_method = "quadrature"), light = TRUE
  )
  cvm_n <- prepare_cvm_observed_data_from_sample_ks(
    z, normal_spec, normal_theta, ks_n,
    control = list(derivative_method = "quadrature")
  )
  ks_l <- prepare_ks_observed_data(
    simplex, lg_spec, lg_theta, make_sample_unique_distance_ks_grid(),
    control = list(derivative_method = "quadrature"), light = TRUE
  )
  cvm_l <- prepare_cvm_observed_data_from_sample_ks(
    simplex, lg_spec, lg_theta, ks_l,
    control = list(derivative_method = "quadrature")
  )
  prep_n <- prepare_mvnormal_fast_multiplier(
    normal_spec, z, normal_theta, ks_n, cvm_n,
    list(derivative_method = "quadrature"), "both"
  )
  prep_l <- prepare_logistic_gaussian_fast_multiplier(
    lg_spec, simplex, lg_theta, ks_l, cvm_l,
    list(derivative_method = "quadrature"), "both"
  )
  expect_equal(prep_l$D_ks$derivative_sorted,
               prep_n$D_ks$derivative_sorted, tolerance = 2e-12)
  expect_equal(prep_l$S_obs, prep_n$S_obs, tolerance = 2e-12)
  expect_equal(prep_l$Vhat, prep_n$Vhat, tolerance = 2e-12)
  expect_equal(
    prep_l$D_ks$derivative_sorted %*% t(solve(prep_l$Vhat)),
    prep_n$D_ks$derivative_sorted %*% t(solve(prep_n$Vhat)),
    tolerance = 2e-12
  )
})

test_that("Gaussian quadrature handles q 1, 2 and 10 and boundary probabilities", {
  for (q in c(1L, 2L, 10L)) {
    Sigma <- diag(seq(0.8, 1.2, length.out = q), nrow = q)
    mu <- seq_len(q) / (10 * q)
    result <- gaussian_ball_profile_quadrature(
      rep.int(0, q), mu, Sigma, c(-1, 0, 1e-5, 1, 10, Inf),
      control = list(gaussian_quadrature_abs_tol = if (q == 2L) 1e-4 else 1e-6)
    )
    expect_equal(dim(result$derivative),
                 c(6L, q + q * (q + 1L) / 2L))
    expect_equal(result$F[1:2], c(0, 0))
    expect_equal(result$F[[6L]], 1)
    expect_true(all(result$F >= 0 & result$F <= 1))
    expect_true(all(is.finite(result$derivative)))
  }
})

test_that("auto selects quadrature while explicit and legacy score_mc remain", {
  set.seed(8127)
  x <- mvtnorm::rmvnorm(10L, c(0, 0), diag(2L))
  common <- list(
    data = x,
    null = list(type = "composite"),
    B = 3L,
    seed = 8128L,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    fast_multiplier_backend = "r"
  )
  automatic <- do.call(multiplier_bootstrap_mvnormal, c(
    common, list(control = list(derivative_method = "auto"))
  ))
  expect_identical(automatic$diagnostics$derivative_method_requested, "auto")
  expect_identical(automatic$diagnostics$derivative_method_effective, "quadrature")
  expect_identical(
    automatic$diagnostics$derivative_method_selection_source, "explicit_auto"
  )
  explicit_mc <- do.call(multiplier_bootstrap_mvnormal, c(
    common,
    list(control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 100L,
      derivative_mc_seed = 8129L
    ))
  ))
  expect_identical(explicit_mc$diagnostics$derivative_method_effective, "score_mc")
  legacy <- expect_warning(
    do.call(multiplier_bootstrap_mvnormal, c(
      common,
      list(control = list(
        derivative_mc_size = 100L,
        derivative_mc_seed = 8129L
      ))
    )),
    "legacy.*Selecting `score_mc`"
  )
  expect_identical(legacy$diagnostics$derivative_method_effective, "score_mc")
  expect_identical(
    legacy$diagnostics$derivative_method_selection_source,
    "legacy_mc_controls"
  )
})
