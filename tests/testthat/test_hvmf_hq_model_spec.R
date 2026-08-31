library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

hvmf_hq_test_mu <- function(q) c(sqrt(2), 1, rep.int(0, q - 1L))

hvmf_hq_spatial_rotation <- function(q) {
  Q <- diag(q)
  Q[, 1L] <- 0
  Q[, 2L] <- 0
  Q[1L, 2L] <- 1
  Q[2L, 1L] <- -1
  Q
}

hvmf_hq_rotate <- function(x, Q) {
  if (is.null(dim(x))) x <- matrix(as.numeric(x), nrow = 1L)
  else x <- as.matrix(x)
  cbind(x[, 1L], x[, -1L, drop = FALSE] %*% Q)
}

test_that("general MLE reduces exactly to the established H2 MLE", {
  set.seed(3101)
  x <- rhvmf_polar(80, mu = c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2)), kappa = 3)
  legacy <- hvmf_mle_h2(x)
  general <- hvmf_mle_hq(x)

  expect_equal(general$mu, legacy$mu, tolerance = 0)
  expect_equal(general$kappa, legacy$kappa, tolerance = 0)
  expect_equal(general$R, legacy$R, tolerance = 0)
})

test_that("general HvMF MLE respects integer weights and solves all score equations", {
  q <- 10L
  set.seed(3102)
  x <- rhvmf_polar(32, mu = hvmf_hq_test_mu(q), kappa = 10)
  weights <- rep(c(1, 2), length.out = nrow(x))
  weighted <- hvmf_mle_hq(x, weights = weights)
  replicated <- hvmf_mle_hq(x[rep(seq_len(nrow(x)), weights), , drop = FALSE])

  expect_equal(weighted$mu, replicated$mu, tolerance = 1e-10)
  expect_equal(weighted$kappa, replicated$kappa, tolerance = 1e-9)

  spec <- make_hvmf_spec(unknown_param = "both")
  prep <- prepare_hvmf_hq_fast_multiplier(
    spec = spec,
    data = x,
    theta_hat = hvmf_mle_hq(x),
    control = list(
      derivative_method = "score_mc", derivative_mc_size = 80L,
      derivative_mc_seed = 3103L
    )
  )
  expect_lt(max(abs(colSums(prep$S_obs))), 2e-7)
  expect_true(all(is.finite(prep$Vhat)))
  expect_gt(rcond(prep$Vhat), 1e-12)
})

test_that("Hq distances and fitted parameters are invariant under spatial rotations", {
  q <- 10L
  Q <- hvmf_hq_spatial_rotation(q)
  set.seed(3104)
  x <- rhvmf_polar(70, mu = hvmf_hq_test_mu(q), kappa = 10)
  x_rotated <- hvmf_hq_rotate(x, Q)

  expect_equal(
    hvmf_distance_matrix_hq(x[1:15, , drop = FALSE], x[16:35, , drop = FALSE]),
    hvmf_distance_matrix_hq(x_rotated[1:15, , drop = FALSE], x_rotated[16:35, , drop = FALSE]),
    tolerance = 1e-10
  )
  fit <- hvmf_mle_hq(x)
  fit_rotated <- hvmf_mle_hq(x_rotated)
  expect_equal(fit_rotated$mu, as.numeric(hvmf_hq_rotate(fit$mu, Q)), tolerance = 1e-10)
  expect_equal(fit_rotated$kappa, fit$kappa, tolerance = 1e-10)
})

test_that("Hq geodesic distances have the defining scalar-product form", {
  q <- 10L
  set.seed(31041)
  x <- rhvmf_polar(8, mu = hvmf_hq_test_mu(q), kappa = 10)
  distances <- hvmf_distance_matrix_hq(x, x)
  manual <- outer(seq_len(nrow(x)), seq_len(nrow(x)), Vectorize(function(i, j) {
    acosh(max(-hvmf_minkowski_inner_product(x[i, ], x[j, ]), 1))
  }))

  expect_equal(distances, manual, tolerance = 5e-8)
  expect_equal(diag(distances), rep(0, nrow(x)), tolerance = 5e-8)
  expect_equal(distances, t(distances), tolerance = 5e-8)
})

test_that("general HvMF score agrees with numerical differentiation in H10", {
  q <- 10L
  set.seed(31042)
  x <- rhvmf_polar(1, mu = hvmf_hq_test_mu(q), kappa = 10)
  eta <- hvmf_hq_test_mu(q)[-1L]
  par <- c(eta, log(10))
  log_density <- function(par) {
    theta <- hvmf_hq_theta_from_coordinates(par, q)
    hvmf_log_normalizing_constant(q, theta$kappa) +
      theta$kappa * hvmf_minkowski_inner_product(x[1L, ], theta$mu)
  }
  numerical <- vapply(seq_along(par), function(j) {
    step <- 1e-6 * max(1, abs(par[[j]]))
    plus <- minus <- par
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step
    (log_density(plus) - log_density(minus)) / (2 * step)
  }, numeric(1))
  theta <- hvmf_hq_theta_from_coordinates(par, q)
  analytic <- c(
    theta$kappa * (x[1L, -1L] - x[1L, 1L] * eta / theta$mu[[1L]]),
    theta$kappa * (hvmf_mean_resultant_ratio(q, theta$kappa) +
      hvmf_minkowski_inner_product(x[1L, ], theta$mu))
  )

  expect_equal(analytic, numerical, tolerance = 2e-6)
})

test_that("general distance profiles agree with direct radial integration", {
  radii <- seq(0, 2.5, length.out = 31L)
  for (configuration in list(
    list(q = 2L, kappa = 2),
    list(q = 2L, kappa = 3),
    list(q = 10L, kappa = 10),
    list(q = 10L, kappa = 15)
  )) {
    q <- configuration$q
    mu <- hvmf_hq_test_mu(q)
    omega <- c(1, rep.int(0, q))
    chi <- acosh(-hvmf_minkowski_inner_product(mu, omega))
    profile <- hvmf_distance_profile_hq(
      omega, mu, configuration$kappa, radii, grid_size = 4097L
    )
    direct <- hvmf_radial_cdf(radii, q, configuration$kappa, chi)

    expect_true(all(profile >= 0 & profile <= 1))
    expect_true(all(diff(profile) >= -1e-12))
    expect_lt(max(abs(profile - direct)), 5e-4)
  }
})

test_that("the H5 model adapter uses adaptive radial integration", {
  q <- 5L
  mu <- hvmf_hq_test_mu(q)
  omega <- c(1, rep.int(0, q))
  radii <- seq(0.05, 2.5, length.out = 13L)
  theta <- list(mu = mu, kappa = 5)
  spec <- make_hvmf_spec(unknown_param = "both")

  direct <- hvmf_distance_profile_hq_integral(
    omega = omega, mu = mu, kappa = theta$kappa, t_values = radii
  )
  profile_small_grid <- spec$profile_eval(
    omega, radii, theta, control = list(hvmf_profile_n_y = 17L)
  )
  profile_large_grid <- spec$profile_eval(
    omega, radii, theta, control = list(hvmf_profile_n_y = 16385L)
  )

  expect_equal(profile_small_grid, direct, tolerance = 1e-12)
  expect_equal(profile_large_grid, direct, tolerance = 1e-12)
  expect_null(spec$sample_profile_matrix_eval(
    data = rbind(mu, omega), distance_matrix = matrix(radii[1:4], nrow = 2L),
    theta = theta, control = list(hvmf_profile_n_y = 17L)
  ))
})

test_that("the adaptive profile override is restricted to H5", {
  q <- 10L
  mu <- hvmf_hq_test_mu(q)
  omega <- c(1, rep.int(0, q))
  radii <- seq(0.05, 2.5, length.out = 9L)
  theta <- list(mu = mu, kappa = 10)
  spec <- make_hvmf_spec(unknown_param = "both")
  expected <- hvmf_distance_profile_hq(
    omega, mu, theta$kappa, radii, grid_size = 257L
  )
  actual <- spec$profile_eval(
    omega, radii, theta, control = list(hvmf_profile_n_y = 257L)
  )

  expect_equal(actual, expected, tolerance = 0)
  expect_true(is.matrix(spec$sample_profile_matrix_eval(
    data = rbind(mu, omega), distance_matrix = matrix(radii[1:4], nrow = 2L),
    theta = theta, control = list(hvmf_profile_n_y = 257L)
  )))
})

test_that("the direct Hq profile has the general radial-density formula", {
  for (configuration in list(
    list(q = 2L, kappa = 2), list(q = 5L, kappa = 5)
  )) {
    q <- configuration$q
    mu <- hvmf_hq_test_mu(q)
    omega <- c(1, rep.int(0, q))
    chi <- acosh(-hvmf_minkowski_inner_product(mu, omega))
    radii <- seq(0.1, 2.5, length.out = 7L)
    direct <- hvmf_distance_profile_hq_integral(
      omega, mu, configuration$kappa, radii
    )
    radial <- hvmf_radial_cdf(
      radii, q = q, kappa = configuration$kappa, chi = chi
    )
    expect_equal(direct, radial, tolerance = 1e-12)
    expect_true(all(diff(direct) >= 0))
    expect_true(all(direct > 0 & direct < 1))
  }
})

test_that("the general radial profile density is normalized beyond H2", {
  for (configuration in list(
    list(q = 2L, kappa = 2, chi = asinh(1)),
    list(q = 5L, kappa = 5, chi = asinh(1)),
    list(q = 5L, kappa = 10, chi = 0.4)
  )) {
    integral <- stats::integrate(
      hvmf_radial_density, lower = 0, upper = Inf,
      q = configuration$q, kappa = configuration$kappa,
      chi = configuration$chi, rel.tol = 1e-9, abs.tol = 1e-11
    )
    expect_equal(integral$value, 1, tolerance = 2e-8)
    expect_lt(integral$abs.error, 1e-7)
  }
})

test_that("H5 integral profiles support the quadrature fast-multiplier path", {
  q <- 5L
  set.seed(31045)
  x <- rhvmf_polar(8, mu = hvmf_hq_test_mu(q), kappa = 5)
  result <- multiplier_bootstrap_hvmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 9L,
    n_cores = 1L,
    seed = 31046L,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "quadrature",
      hvmf_profile_n_y = 1025L,
      hvmf_derivative_n_y = 1025L,
      fast_multiplier_cvm_block_size = 8L,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE
    )
  )

  expect_identical(result$diagnostics$derivative_method_effective, "quadrature")
  expect_identical(result$diagnostics$fast_multiplier_backend_effective, "cpp")
  expect_true(result$diagnostics$fast_multiplier_fuse_ks_cvm_effective)
  expect_true(all(c(
    result$inference$ks$p_value,
    result$inference$cvm$p_value
  ) >= 0 & c(
    result$inference$ks$p_value,
    result$inference$cvm$p_value
  ) <= 1))
})

test_that("generic H2 profile agrees with the established exact H2 profile", {
  mu <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))
  omega <- c(1, 0, 0)
  radii <- seq(0.05, 2.5, length.out = 21L)
  generic <- hvmf_distance_profile_hq(
    omega, mu, kappa = 3, t_values = radii, grid_size = 8193L
  )
  established <- theoretical_distance_profile_hvmf(
    omega, mu, kappa = 3, t_values = radii
  )

  expect_equal(generic, established, tolerance = 6e-4)

  direct <- hvmf_distance_profile_hq_integral(
    omega, mu, kappa = 3, t_values = radii
  )
  expect_equal(direct, established, tolerance = 2e-6)
})

test_that("H10 fast multiplier agrees with the reestimated observed statistics", {
  q <- 10L
  set.seed(3105)
  x <- rhvmf_polar(30, mu = hvmf_hq_test_mu(q), kappa = 10)
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 29L,
    alpha = 0.05,
    n_cores = 1L,
    seed = 3106L,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "score_mc", derivative_mc_size = 100L,
      derivative_mc_seed = 3107L, fast_multiplier_cvm_block_size = 25L,
      hvmf_profile_method = "tabulated", hvmf_profile_n_y = 1025L
    ),
    unknown_param = "both"
  )
  fast <- do.call(multiplier_bootstrap_hvmf, c(
    common,
    list(
      bootstrap_method = "fast_multiplier", fast_multiplier_backend = "cpp",
      fuse_ks_cvm = TRUE, cache_block_corrections = "auto"
    )
  ))
  reestimated <- do.call(multiplier_bootstrap_hvmf, c(
    common, list(bootstrap_method = "reestimated")
  ))

  expect_identical(fast$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_identical(fast$diagnostics$fast_multiplier_backend_effective, "cpp")
  expect_true(fast$diagnostics$fast_multiplier_fuse_ks_cvm_effective)
  expect_equal(fast$inference$ks$observed, reestimated$inference$ks$observed, tolerance = 1e-10)
  expect_equal(fast$inference$cvm$observed, reestimated$inference$cvm$observed, tolerance = 1e-10)
  expect_true(all(c(
    fast$inference$ks$p_value, fast$inference$cvm$p_value,
    reestimated$inference$ks$p_value, reestimated$inference$cvm$p_value
  ) >= 0 & c(
    fast$inference$ks$p_value, fast$inference$cvm$p_value,
    reestimated$inference$ks$p_value, reestimated$inference$cvm$p_value
  ) <= 1))
})
