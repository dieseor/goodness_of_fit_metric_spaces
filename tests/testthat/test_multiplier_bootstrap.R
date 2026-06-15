library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

# testthat::test_dir("tests/testthat")

source(file.path("bootstrap", "multiplier_bootstrap.R"))

test_that("multiplier helpers normalize weights correctly", {
  normalized <- normalize_multiplier_weights(c(1, 2, 3))
  expect_equal(normalized, c(0.5, 1, 1.5), tolerance = 1e-12)
  expect_equal(sum(normalized), 3, tolerance = 1e-12)

  spec <- resolve_multiplier_spec(NULL)
  draws <- generate_multiplier_matrix(B = 3, n = 4, multiplier_spec = spec, seed = 123)
  expect_equal(dim(draws), c(3, 4))
  expect_true(all(draws >= 0))
})

test_that("KS cached grid profiles match the naive computation", {
  distance_matrix <- matrix(c(
    0.10, 0.40, 0.70,
    0.20, 0.50, 0.90
  ), nrow = 3, ncol = 2)
  t_grid <- c(0.25, 0.60)
  weights <- c(0.6, 0.9, 1.5)

  order_matrix <- t(vapply(seq_len(ncol(distance_matrix)), function(j) {
    as.integer(order(distance_matrix[, j]))
  }, integer(nrow(distance_matrix))))

  sorted_distance_matrix <- matrix(0, nrow = ncol(distance_matrix), ncol = nrow(distance_matrix))
  for (j in seq_len(ncol(distance_matrix))) {
    sorted_distance_matrix[j, ] <- distance_matrix[order_matrix[j, ], j]
  }

  threshold_index_matrix <- t(vapply(seq_len(ncol(distance_matrix)), function(j) {
    as.integer(findInterval(t_grid, sorted_distance_matrix[j, ]))
  }, integer(length(t_grid))))

  fast_empirical <- compute_grid_empirical_profile(
    distance_matrix,
    t_grid,
    sorted_distance_matrix = sorted_distance_matrix,
    threshold_index_matrix = threshold_index_matrix
  )
  slow_empirical <- compute_grid_empirical_profile(distance_matrix, t_grid)

  fast_weighted <- compute_grid_weighted_profile(
    distance_matrix,
    t_grid,
    weights,
    sorted_distance_matrix = sorted_distance_matrix,
    order_matrix = order_matrix,
    threshold_index_matrix = threshold_index_matrix
  )
  slow_weighted <- compute_grid_weighted_profile(distance_matrix, t_grid, weights)

  expect_equal(fast_empirical, slow_empirical, tolerance = 1e-12)
  expect_equal(fast_weighted, slow_weighted, tolerance = 1e-12)
})

test_that("normal weighted estimators match analytic formulas", {
  x <- c(0, 2)
  weights <- c(1, 3)

  spec_mu <- make_normal_spec(unknown_param = "mu")
  theta_mu <- spec_mu$fit_theta(
    data = x,
    weights = weights,
    null = list(type = "composite", fixed = list(sigma = 1)),
    control = list()
  )
  expect_equal(theta_mu$mu, 1.5, tolerance = 1e-12)
  expect_equal(theta_mu$sigma, 1, tolerance = 1e-12)

  spec_sigma <- make_normal_spec(unknown_param = "sigma")
  theta_sigma <- spec_sigma$fit_theta(
    data = x,
    weights = weights,
    null = list(type = "composite", fixed = list(mu = 1)),
    control = list()
  )
  expect_equal(theta_sigma$mu, 1, tolerance = 1e-12)
  expect_equal(theta_sigma$sigma, 1, tolerance = 1e-12)

  spec_both <- make_normal_spec(unknown_param = "both")
  theta_both <- spec_both$fit_theta(
    data = x,
    weights = weights,
    null = list(type = "composite"),
    control = list()
  )
  expect_equal(theta_both$mu, 1.5, tolerance = 1e-12)
  expect_equal(theta_both$sigma, sqrt(0.75), tolerance = 1e-12)

  theta_simple <- spec_both$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "simple", theta = list(mu = 0, sigma = 2)),
    control = list()
  )
  expect_equal(theta_simple$mu, 0, tolerance = 1e-12)
  expect_equal(theta_simple$sigma, 2, tolerance = 1e-12)
})

test_that("weighted vMF estimator satisfies the resultant equation", {
  x <- rbind(
    c(1, 0, 0),
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  weights <- c(2, 1, 1, 1)

  spec <- make_vmf_spec(distance_type = "chordal", unknown_param = "xi")
  theta_hat <- spec$fit_theta(
    data = x,
    weights = weights,
    null = list(type = "composite"),
    control = list()
  )

  prob_weights <- normalize_probability_weights(weights, n_expected = nrow(x))
  x_normalized <- normalize_vmf_data(x)
  resultant <- colSums(x_normalized * prob_weights)
  r_bar <- sqrt(sum(resultant^2))

  expect_equal(sum(theta_hat$mu^2), 1, tolerance = 1e-10)
  expect_true(is.finite(theta_hat$kappa))
  expect_equal(theta_hat$xi, theta_hat$kappa * theta_hat$mu, tolerance = 1e-10)
  expect_equal(A_q(theta_hat$kappa, q = 2), r_bar, tolerance = 1e-8)

  theta_simple <- spec$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "simple", theta = list(mu = c(1, 0, 0), kappa = 2)),
    control = list()
  )
  expect_equal(theta_simple$kappa, 2, tolerance = 1e-12)
  expect_equal(theta_simple$mu, c(1, 0, 0), tolerance = 1e-12)
})


test_that("JP weighted MLE accepts a warm start", {
  skip_if_not(exists("jp_mle_s2_weighted", mode = "function"))

  x <- rbind(
    c(0.00, 0.00, 1.00),
    c(0.20, 0.00, 0.98),
    c(-0.15, 0.05, 0.9873),
    c(0.05, -0.20, 0.9787),
    c(0.30, 0.10, 0.9487),
    c(-0.10, -0.25, 0.9631)
  )
  x <- normalize_jp_data(x)

  control <- list(
    jp_mle_sign_branches = c(1L),
    jp_mle_psi_abs_starts = c(0.5),
    jp_mle_maxit = 40L,
    jp_mle_reltol = 1e-6
  )

  theta_start <- jp_mle_s2_weighted(
    data = x,
    weights = NULL,
    control = control
  )

  theta_warm <- jp_mle_s2_weighted(
    data = x,
    weights = c(1.2, 0.8, 1.0, 1.1, 0.9, 1.0),
    control = modifyList(control, list(jp_mle_start_theta = theta_start))
  )

  expect_true(is.finite(theta_warm$kappa))
  expect_true(is.finite(theta_warm$psi))
  expect_equal(sum(theta_warm$mu^2), 1, tolerance = 1e-8)
  if (!is.null(theta_warm$loglik)) {
    expect_true(is.finite(theta_warm$loglik))
  }
})

test_that("normal multiplier bootstrap is reproducible and returns valid structure", {
  x <- c(-1.2, -0.5, 0.2, 0.9, 1.4, 2.1)
  ks_grid <- list(
    omega_grid = seq(-2, 2, length.out = 4),
    t_grid = c(0.25, 0.75, 1.5)
  )
  null <- list(type = "simple", theta = list(mu = 0, sigma = 1))

  result_1 <- multiplier_bootstrap_normal(
    data = x,
    null = null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1
  )

  result_2 <- multiplier_bootstrap_normal(
    data = x,
    null = null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(is.list(result_1$observed))
  expect_true(is.list(result_1$inference))
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)
  expect_true(is.finite(result_1$inference$ks$critical_value))
  expect_true(is.finite(result_1$inference$cvm$critical_value))
})

test_that("observed CvM matches the exact double-sum formula", {
  x <- c(-1, 0, 2)
  null <- list(type = "simple", theta = list(mu = 0, sigma = 1))

  result <- multiplier_bootstrap_normal(
    data = x,
    null = null,
    statistics = "cvm",
    B = 4,
    seed = 7,
    n_cores = 1
  )

  n <- length(x)
  distance_matrix <- abs(outer(x, x, FUN = "-"))
  rank_matrix <- t(vapply(seq_len(n), function(i) {
    as.integer(rank(distance_matrix[i, ], ties.method = "max"))
  }, integer(n)))
  empirical_profile <- rank_matrix / n
  theoretical_profile <- t(vapply(seq_len(n), function(i) {
    theoretical_distance_profile_normal(x[i], mu = 0, sigma = 1, t_values = distance_matrix[i, ])
  }, numeric(n)))
  process_matrix <- sqrt(n) * (empirical_profile - theoretical_profile)
  manual_cvm <- mean(process_matrix^2)

  expect_equal(result$observed$cvm$process_matrix, process_matrix, tolerance = 1e-12)
  expect_equal(result$observed$cvm$statistic, manual_cvm, tolerance = 1e-12)
})

test_that("bootstrap CvM uses the exact paper formula with P_n integration", {
  x <- c(-1, 0, 2)
  null <- list(type = "simple", theta = list(mu = 0, sigma = 1))

  result <- multiplier_bootstrap_normal(
    data = x,
    null = null,
    statistics = "cvm",
    B = 1,
    seed = 101,
    n_cores = 1
  )

  n <- length(x)
  multiplier_spec <- resolve_multiplier_spec(NULL)
  raw_draw <- generate_multiplier_matrix(B = 1, n = n, multiplier_spec = multiplier_spec, seed = 101)[1, ]
  normalized_weights <- normalize_multiplier_weights(raw_draw)

  distance_matrix <- abs(outer(x, x, FUN = "-"))
  rank_matrix <- t(vapply(seq_len(n), function(i) {
    as.integer(rank(distance_matrix[i, ], ties.method = "max"))
  }, integer(n)))
  empirical_profile <- rank_matrix / n
  order_list <- lapply(seq_len(n), function(i) {
    order(distance_matrix[i, ])
  })
  f_star <- compute_weighted_sample_profile_matrix(
    order_list = order_list,
    rank_matrix = rank_matrix,
    normalized_weights = normalized_weights
  )
  process_star <- sqrt(n) * (f_star - empirical_profile)
  manual_bootstrap_cvm <- mean(process_star^2)

  expect_equal(result$bootstrap$statistics$cvm[[1]], manual_bootstrap_cvm, tolerance = 1e-12)
})

test_that("normal composite wrappers run for mu, sigma and both", {
  x <- c(-0.8, -0.2, 0.1, 0.5, 1.4, 2.0)
  ks_grid <- list(
    omega_grid = seq(-1.5, 2.5, length.out = 3),
    t_grid = c(0.3, 0.8)
  )

  result_mu <- multiplier_bootstrap_normal(
    data = x,
    null = list(type = "composite", fixed = list(sigma = 1)),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 5,
    seed = 11,
    n_cores = 1,
    unknown_param = "mu"
  )
  expect_true(is.finite(result_mu$observed$theta_hat$mu))
  expect_equal(length(result_mu$bootstrap$statistics$ks), 5)
  expect_equal(length(result_mu$bootstrap$statistics$cvm), 5)

  result_sigma <- multiplier_bootstrap_normal(
    data = x,
    null = list(type = "composite", fixed = list(mu = 0)),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 5,
    seed = 12,
    n_cores = 1,
    unknown_param = "sigma"
  )
  expect_true(is.finite(result_sigma$observed$theta_hat$sigma))

  result_both <- multiplier_bootstrap_normal(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 5,
    seed = 13,
    n_cores = 1,
    unknown_param = "both"
  )
  expect_true(is.finite(result_both$observed$theta_hat$mu))
  expect_true(is.finite(result_both$observed$theta_hat$sigma))
})

test_that("vMF simple and composite wrappers run with both statistics", {
  set.seed(123)
  x <- rotasym::r_vMF(10, mu = c(1, 0, 0), kappa = 3)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(3, dim = 3),
    t_grid = c(0.35, 0.8)
  )

  result_simple <- multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "simple", theta = list(mu = c(1, 0, 0), kappa = 3)),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 4,
    seed = 21,
    n_cores = 1,
    distance_type = "chordal"
  )

  expect_true(is.finite(result_simple$observed$theta_hat$kappa))
  expect_equal(length(result_simple$bootstrap$statistics$ks), 4)
  expect_equal(length(result_simple$bootstrap$statistics$cvm), 4)

  result_composite <- multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 4,
    seed = 22,
    n_cores = 1,
    distance_type = "chordal",
    unknown_param = "xi"
  )

  expect_true(is.finite(result_composite$observed$theta_hat$kappa))
  expect_equal(sum(result_composite$observed$theta_hat$mu^2), 1, tolerance = 1e-8)
  expect_equal(length(result_composite$bootstrap$statistics$ks), 4)
  expect_equal(length(result_composite$bootstrap$statistics$cvm), 4)
})

test_that("tabulated S2 vMF profile stays close to the exact integral", {
  validation <- validate_vmf_s2_grid_profile(
    n_checks = 30,
    kappa_values = c(0.5, 2),
    n_u = 4097L,
    distance_type = "geodesic",
    seed = 1
  )

  expect_true(validation[["max_error"]] < 5e-4)
  expect_true(validation[["q95_error"]] < 2e-4)
})

test_that("tabulated S2 vMF CvM matrix matches the exact row-wise profile closely", {
  set.seed(99)
  x <- rotasym::r_vMF(8, mu = c(1, 0, 0), kappa = 2)
  x <- normalize_vmf_data(x)
  theta <- list(mu = c(1, 0, 0), kappa = 2)
  distance_matrix <- acos(pmax(pmin(x %*% t(x), 1), -1))

  exact_matrix <- compute_theoretical_sample_profile_matrix(
    spec = make_vmf_spec(distance_type = "geodesic", unknown_param = "xi"),
    data = x,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(vmf_profile_method = "exact")
  )
  fast_matrix <- compute_theoretical_sample_profile_matrix(
    spec = make_vmf_spec(distance_type = "geodesic", unknown_param = "xi"),
    data = x,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(vmf_profile_method = "tabulated", vmf_profile_n_u = 4097L)
  )

  expect_true(max(abs(exact_matrix - fast_matrix)) < 5e-4)
})

test_that("vMF fast multiplier bootstrap runs and matches the observed statistic", {
  set.seed(777)
  x <- rotasym::r_vMF(12, mu = c(1, 0, 0), kappa = 2)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(3, dim = 3),
    t_grid = c(0.25, 0.75)
  )

  result_old <- multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 4,
    seed = 901,
    n_cores = 1,
    distance_type = "geodesic",
    unknown_param = "xi",
    bootstrap_method = "reestimated"
  )
  result_fast <- multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 4,
    seed = 901,
    n_cores = 1,
    distance_type = "geodesic",
    unknown_param = "xi",
    bootstrap_method = "fast_multiplier",
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 4000L,
      derivative_mc_seed = 902,
      fast_multiplier_cvm_block_size = 12L
    )
  )

  expect_equal(result_fast$observed$ks$statistic, result_old$observed$ks$statistic, tolerance = 1e-12)
  expect_equal(result_fast$observed$cvm$statistic, result_old$observed$cvm$statistic, tolerance = 1e-12)
  expect_equal(length(result_fast$bootstrap$statistics$ks), 4)
  expect_equal(length(result_fast$bootstrap$statistics$cvm), 4)
  expect_true(all(is.finite(result_fast$bootstrap$statistics$ks)))
  expect_true(all(is.finite(result_fast$bootstrap$statistics$cvm)))
  expect_identical(result_fast$diagnostics$bootstrap_method, "fast_multiplier")
  expect_identical(result_fast$diagnostics$derivative_method, "score_mc")
  expect_equal(result_fast$diagnostics$derivative_mc_size, 4000L)
})

test_that("vMF fast multiplier preparation has the expected dimensions and sign convention", {
  set.seed(778)
  x <- normalize_vmf_data(rotasym::r_vMF(9, mu = c(1, 0, 0), kappa = 2))
  theta_hat <- fit_vmf_theta(x, weights = NULL, null = list(type = "composite"))
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(3, dim = 3),
    t_grid = c(0.4, 0.8)
  )
  ks_prep <- prepare_ks_observed_data(x, make_vmf_spec("geodesic", "xi"), theta_hat, ks_grid)
  cvm_prep <- prepare_cvm_observed_data(x, make_vmf_spec("geodesic", "xi"), theta_hat, control = list(vmf_profile_method = "tabulated"))
  prep <- spec_fast_multiplier_prepare(
    spec = make_vmf_spec("geodesic", "xi"),
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = list(derivative_mc_size = 3000L, derivative_mc_seed = 123L)
  )

  expect_equal(dim(prep$S_obs), c(nrow(x), 3))
  expect_equal(dim(prep$Vhat), c(3, 3))
  expect_equal(dim(prep$D_ks), c(nrow(ks_grid$omega_grid) * length(ks_grid$t_grid), 3))
  expect_equal(dim(prep$D_cvm), c(nrow(x) * nrow(x), 3))

  q <- length(theta_hat$mu) - 1L
  A_q_kappa <- A_q(theta_hat$kappa, q)
  scalar_coef <- 1 - A_q_kappa^2 - ((q + 1) * A_q_kappa / theta_hat$kappa)
  var_X <- (A_q_kappa / theta_hat$kappa) * diag(q + 1) + scalar_coef * outer(theta_hat$mu, theta_hat$mu)
  expect_equal(prep$Vhat, var_X, tolerance = 1e-8)

  correction_v <- prep$S_obs %*% solve(prep$Vhat) %*% t(prep$D_ks)
  correction_a <- -prep$S_obs %*% solve(dot_psi_xi(theta_hat$xi, q)) %*% t(prep$D_ks)
  expect_equal(correction_v, correction_a, tolerance = 1e-8)
})

test_that("vMF score-MC derivative sign agrees with a finite-difference sanity check", {
  theta_hat <- normalize_vmf_theta(list(mu = c(1, 0, 0), kappa = 2))
  ks_grid <- list(
    omega_grid = rbind(c(1, 0, 0)),
    t_grid = c(0.5)
  )
  x <- normalize_vmf_data(rotasym::r_vMF(8, mu = theta_hat$mu, kappa = theta_hat$kappa))
  ks_prep <- prepare_ks_observed_data(x, make_vmf_spec("geodesic", "xi"), theta_hat, ks_grid)
  prep <- spec_fast_multiplier_prepare(
    spec = make_vmf_spec("geodesic", "xi"),
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = NULL,
    control = list(derivative_mc_size = 30000L, derivative_mc_seed = 456L)
  )

  eps <- 1e-4
  for (j in 1:2) {
    xi_plus <- theta_hat$xi
    xi_minus <- theta_hat$xi
    xi_plus[j] <- xi_plus[j] + eps
    xi_minus[j] <- xi_minus[j] - eps
    theta_plus <- normalize_vmf_theta(xi_plus)
    theta_minus <- normalize_vmf_theta(xi_minus)
    fd_value <- (
      theoretical_distance_profile_vmf_s2_fast(ks_grid$omega_grid[1, ], theta_plus$mu, theta_plus$kappa, ks_grid$t_grid[1], "geodesic") -
        theoretical_distance_profile_vmf_s2_fast(ks_grid$omega_grid[1, ], theta_minus$mu, theta_minus$kappa, ks_grid$t_grid[1], "geodesic")
    ) / (2 * eps)
    expect_lt(abs(prep$D_ks[1, j] - fd_value), 0.08)
  }
})

test_that("additional supported models run through the fast multiplier branch", {
  set.seed(204)

  x_hvmf <- rhvmf_h2_gig(12, mu = c(cosh(0.35), sinh(0.35), 0), kappa = 3)
  ks_grid_hvmf <- list(
    omega_grid = x_hvmf[1:3, , drop = FALSE],
    t_grid = c(0.3, 0.8)
  )
  result_hvmf <- multiplier_bootstrap_hvmf(
    data = x_hvmf,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = ks_grid_hvmf,
    B = 3,
    seed = 10,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    control = list(derivative_mc_size = 600L, derivative_mc_seed = 11L)
  )
  expect_identical(result_hvmf$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_true(all(is.finite(result_hvmf$bootstrap$statistics$ks)))

  x_sc <- r_sph_small_circle(12, mu = c(0, 0, 1), kappa = 8, nu = 0.3)
  result_sc <- multiplier_bootstrap_small_circle(
    data = x_sc,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(4, dim = 3), t_grid = c(0.5, 1.0)),
    B = 3,
    seed = 12,
    n_cores = 1,
    distance_type = "geodesic",
    bootstrap_method = "fast_multiplier",
    control = list(derivative_mc_size = 600L, derivative_mc_seed = 13L)
  )
  expect_identical(result_sc$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_true(all(is.finite(result_sc$bootstrap$statistics$ks)))
})

test_that("mixture fast multiplier preparations expose the expected dimensions", {
  set.seed(301)
  theta_bm <- list(
    mu = c(0, 0, 1),
    weight1 = 0.45,
    alpha1 = 2,
    beta1 = 7,
    alpha2 = 7,
    beta2 = 2
  )
  x_bm <- r_sph_beta_mixture2(
    20,
    mu = theta_bm$mu,
    weight1 = theta_bm$weight1,
    alpha1 = theta_bm$alpha1,
    beta1 = theta_bm$beta1,
    alpha2 = theta_bm$alpha2,
    beta2 = theta_bm$beta2
  )
  spec_bm <- make_beta_mixture2_spec(distance_type = "geodesic")
  ks_grid <- list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8))
  ks_prep_bm <- prepare_ks_observed_data(x_bm, spec_bm, theta_bm, ks_grid)
  prep_bm <- spec_fast_multiplier_prepare(
    spec = spec_bm,
    data = x_bm,
    theta_hat = theta_bm,
    ks_prep = ks_prep_bm,
    cvm_prep = NULL,
    control = list(derivative_mc_size = 800L, derivative_mc_seed = 21L)
  )
  expect_equal(dim(prep_bm$S_obs), c(nrow(x_bm), 7))
  expect_equal(dim(prep_bm$Psi_aux), c(800, 7))
  expect_equal(dim(prep_bm$D_ks), c(nrow(ks_grid$omega_grid) * length(ks_grid$t_grid), 7))

  theta_scw <- list(
    mu = c(0, 0, 1),
    pi = 0.55,
    kappa1 = 10,
    nu1 = 0.4,
    kappa2 = 8,
    nu2 = 0.3
  )
  x_scw <- r_sph_small_circle_weighted_mixture2(
    20,
    mu = theta_scw$mu,
    pi = theta_scw$pi,
    kappa1 = theta_scw$kappa1,
    nu1 = theta_scw$nu1,
    kappa2 = theta_scw$kappa2,
    nu2 = theta_scw$nu2
  )
  spec_scw <- make_small_circle_weighted_mixture2_spec(distance_type = "geodesic")
  ks_prep_scw <- prepare_ks_observed_data(x_scw, spec_scw, theta_scw, ks_grid)
  prep_scw <- spec_fast_multiplier_prepare(
    spec = spec_scw,
    data = x_scw,
    theta_hat = theta_scw,
    ks_prep = ks_prep_scw,
    cvm_prep = NULL,
    control = list(derivative_mc_size = 800L, derivative_mc_seed = 22L)
  )
  expect_equal(dim(prep_scw$S_obs), c(nrow(x_scw), 7))
  expect_equal(dim(prep_scw$Psi_aux), c(800, 7))
  expect_equal(dim(prep_scw$D_ks), c(nrow(ks_grid$omega_grid) * length(ks_grid$t_grid), 7))
})

test_that("cardioid rho near one uses the recorded slow fallback", {
  x_card <- r_sph_car(12, mu = c(0, 0, 1), rho = 0.4, k = 2)
  spec_card <- make_cardioid_spec(k = 2, distance_type = "geodesic", unknown_param = "both")
  result <- multiplier_bootstrap_gof(
    data = x_card,
    spec = spec_card,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    B = 3,
    seed = 50,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    observed_theta_hat = list(mu = c(0, 0, 1), rho = 1, k = 2),
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    control = list(derivative_mc_size = 500L, derivative_mc_seed = 51L)
  )
  expect_identical(result$diagnostics$bootstrap_method, "fast_multiplier")
  expect_identical(result$diagnostics$effective_bootstrap_method, "reestimated")
  expect_true(isTRUE(result$diagnostics$fallback_to_reestimated))
  expect_identical(result$diagnostics$fallback_reason, "cardioid_rho_one_boundary")
})

test_that("fast multiplier defaults to Vhat from the auxiliary score outer product", {
  x <- c(-1.4, -0.2, 0.1, 1.0, 1.8, 2.4)
  spec <- make_normal_spec(unknown_param = "both")
  theta_hat <- spec$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "composite"),
    control = list()
  )
  par0 <- c(theta_hat$mu, log(theta_hat$sigma))

  score_matrix_fn <- function(sample, par) {
    mu <- par[[1L]]
    sigma <- exp(par[[2L]])
    cbind(
      (sample - mu) / sigma^2,
      -1 + ((sample - mu)^2) / sigma^2
    )
  }

  sample_fn <- function(n_aux, par) {
    stats::rnorm(n_aux, mean = par[[1L]], sd = exp(par[[2L]]))
  }

  ks_grid <- list(omega_grid = seq(-1, 1, length.out = 2), t_grid = c(0.5, 1))
  ks_prep <- prepare_ks_observed_data(x, spec, theta_hat, ks_grid)
  prep <- prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = NULL,
    control = list(derivative_mc_size = 2000L, derivative_mc_seed = 901L),
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )

  expect_identical(prep$vhat_method, "score_outer_product")
  expect_equal(prep$Vhat, crossprod(prep$Psi_aux) / nrow(prep$Psi_aux), tolerance = 1e-12)
  expect_equal(dim(prep$S_obs), c(length(x), 2))
  expect_equal(dim(prep$Psi_aux), c(2000, 2))
  expect_equal(dim(prep$Vhat), c(2, 2))
  expect_equal(dim(prep$D_ks), c(length(ks_grid$omega_grid) * length(ks_grid$t_grid), 2))
  expect_lt(prep$vhat_diagnostics$score_mean_aux_norm, 0.12)
  expect_gt(min(prep$vhat_diagnostics$Vhat_eigenvalues), 0)
})

test_that("numeric Jacobian Vhat is available only when requested", {
  x <- c(-1.4, -0.2, 0.1, 1.0, 1.8, 2.4)
  spec <- make_normal_spec(unknown_param = "both")
  theta_hat <- spec$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "composite"),
    control = list()
  )
  par0 <- c(theta_hat$mu, log(theta_hat$sigma))

  score_matrix_fn <- function(sample, par) {
    mu <- par[[1L]]
    sigma <- exp(par[[2L]])
    cbind(
      (sample - mu) / sigma^2,
      -1 + ((sample - mu)^2) / sigma^2
    )
  }

  sample_fn <- function(n_aux, par) {
    stats::rnorm(n_aux, mean = par[[1L]], sd = exp(par[[2L]]))
  }

  ks_grid <- list(omega_grid = seq(-1, 1, length.out = 2), t_grid = c(0.5, 1))
  ks_prep <- prepare_ks_observed_data(x, spec, theta_hat, ks_grid)
  prep <- prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = NULL,
    control = list(
      derivative_mc_size = 1000L,
      derivative_mc_seed = 902L,
      fast_multiplier_vhat_method = "numeric_jacobian"
    ),
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )

  expected_vhat <- -fast_multiplier_numeric_jacobian(
    fun = function(par) colMeans(score_matrix_fn(x, par)),
    x0 = par0
  )
  expect_identical(prep$vhat_method, "numeric_jacobian")
  expect_equal(prep$Vhat, expected_vhat, tolerance = 1e-6)
})

test_that("logistic Gaussian fast preparation uses Gaussian MLE influence correction", {
  set.seed(1101)
  x <- rlogistic_gaussian_simplex(
    n = 12,
    mu_ilr = c(0.25, -0.35),
    Sigma_ilr = matrix(c(0.45, 0.08, 0.08, 0.30), nrow = 2L)
  )
  spec <- make_logistic_gaussian_spec(unknown_param = "both")
  theta_hat <- spec$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "composite"),
    control = list()
  )
  ks_grid <- list(
    omega_grid = x[1:4, , drop = FALSE],
    t_grid = c(0.35, 0.8)
  )
  ks_prep <- prepare_ks_observed_data(x, spec, theta_hat, ks_grid)
  cvm_prep <- prepare_cvm_observed_data(x, spec, theta_hat, control = list())
  prep <- spec_fast_multiplier_prepare(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = list(derivative_mc_size = 700L, derivative_mc_seed = 1102L)
  )

  d <- length(theta_hat$mu_ilr)
  p_sigma <- d * (d + 1L) / 2L
  centered_obs <- sweep(logistic_gaussian_ilr_matrix(x), 2L, theta_hat$mu_ilr, FUN = "-")
  expected_if <- cbind(
    centered_obs,
    t(vapply(seq_len(nrow(centered_obs)), function(i) {
      rr <- centered_obs[i, , drop = FALSE]
      fast_multiplier_vech(crossprod(rr) - theta_hat$Sigma_ilr)
    }, numeric(p_sigma)))
  )

  expect_identical(prep$vhat_method, "gaussian_mle_influence_identity")
  expect_equal(dim(prep$S_obs), c(nrow(x), d + p_sigma))
  expect_equal(dim(prep$Psi_aux), c(700, d + p_sigma))
  expect_equal(dim(prep$Vhat), c(d + p_sigma, d + p_sigma))
  expect_equal(dim(prep$D_ks), c(nrow(ks_grid$omega_grid) * length(ks_grid$t_grid), d + p_sigma))
  expect_equal(dim(prep$D_cvm), c(nrow(x) * nrow(x), d + p_sigma))
  expect_equal(prep$S_obs, expected_if, tolerance = 1e-12)
  expect_equal(prep$Vhat, diag(d + p_sigma), tolerance = 1e-12)
  expect_lt(prep$vhat_diagnostics$score_mean_aux_norm, 0.35)
})

test_that("logistic Gaussian fast multiplier runs without Monte Carlo Vhat inversion", {
  set.seed(1103)
  x <- rlogistic_gaussian_simplex(
    n = 16,
    mu_ilr = c(0.2, -0.25),
    Sigma_ilr = matrix(c(0.35, 0.05, 0.05, 0.28), nrow = 2L)
  )
  ks_grid <- list(
    omega_grid = x[1:5, , drop = FALSE],
    t_grid = c(0.3, 0.7)
  )
  result <- multiplier_bootstrap_logistic_gaussian(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 4,
    seed = 1104,
    n_cores = 1,
    unknown_param = "both",
    bootstrap_method = "fast_multiplier",
    control = list(
      derivative_mc_size = 800L,
      derivative_mc_seed = 1105L,
      fast_multiplier_cvm_block_size = 16L
    )
  )

  expect_identical(result$diagnostics$bootstrap_method, "fast_multiplier")
  expect_identical(result$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_identical(result$diagnostics$vhat_method, "gaussian_mle_influence_identity")
  expect_true(all(is.finite(result$bootstrap$statistics$ks)))
  expect_true(all(is.finite(result$bootstrap$statistics$cvm)))
})

test_that("cardioid rho near zero uses the recorded slow fallback", {
  x_card <- r_sph_car(12, mu = c(0, 0, 1), rho = 0.3, k = 2)
  spec_card <- make_cardioid_spec(k = 2, distance_type = "geodesic", unknown_param = "both")
  result <- multiplier_bootstrap_gof(
    data = x_card,
    spec = spec_card,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    B = 3,
    seed = 60,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    observed_theta_hat = list(mu = c(0, 0, 1), rho = 0, k = 2),
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    control = list(derivative_mc_size = 500L, derivative_mc_seed = 61L)
  )
  expect_identical(result$diagnostics$bootstrap_method, "fast_multiplier")
  expect_identical(result$diagnostics$effective_bootstrap_method, "reestimated")
  expect_true(isTRUE(result$diagnostics$fallback_to_reestimated))
  expect_identical(result$diagnostics$fallback_reason, "cardioid_rho_zero_nonidentification")
})

test_that("cardioid interior still uses the fast branch", {
  x_card <- r_sph_car(12, mu = c(0, 0, 1), rho = 0.3, k = 2)
  spec_card <- make_cardioid_spec(k = 2, distance_type = "geodesic", unknown_param = "both")
  result <- multiplier_bootstrap_gof(
    data = x_card,
    spec = spec_card,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    B = 3,
    seed = 62,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    observed_theta_hat = list(mu = c(0, 0, 1), rho = 0.3, k = 2),
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    control = list(derivative_mc_size = 500L, derivative_mc_seed = 63L)
  )
  expect_identical(result$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_false(isTRUE(result$diagnostics$fallback_to_reestimated))
})

test_that("JP near psi zero uses the fast vMF limit branch", {
  x <- rotasym::r_vMF(10, mu = c(0, 0, 1), kappa = 2)
  spec <- make_jp_spec(distance_type = "geodesic")
  theta_hat <- list(mu = c(0, 0, 1), kappa = 2, psi = 0)
  ks_grid <- list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8))
  ks_prep <- prepare_ks_observed_data(x, spec, theta_hat, ks_grid)
  prep <- spec_fast_multiplier_prepare(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = NULL,
    control = list(derivative_mc_size = 700L, derivative_mc_seed = 71L)
  )

  expect_equal(dim(prep$S_obs), c(nrow(x), 3))
  expect_equal(dim(prep$Psi_aux), c(700, 3))
  expect_equal(dim(prep$Vhat), c(3, 3))
})

test_that("JP alpha boundary uses the recorded slow fallback", {
  x <- rotasym::r_vMF(10, mu = c(0, 0, 1), kappa = 2)
  spec <- make_jp_spec(distance_type = "geodesic")
  result <- multiplier_bootstrap_gof(
    data = x,
    spec = spec,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    B = 3,
    seed = 73,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    observed_theta_hat = list(mu = c(0, 0, 1), kappa = 11.2, psi = 0.91),
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    control = list(derivative_mc_size = 700L, derivative_mc_seed = 74L)
  )
  expect_identical(result$diagnostics$effective_bootstrap_method, "reestimated")
  expect_true(isTRUE(result$diagnostics$fallback_to_reestimated))
  expect_identical(result$diagnostics$fallback_reason, "jp_alpha_boundary")
})

test_that("JP interior still uses the fast branch", {
  x <- r_sph_jp(10, mu = c(0, 0, 1), kappa = 2, psi = 0.7)
  spec <- make_jp_spec(distance_type = "geodesic")
  result <- multiplier_bootstrap_gof(
    data = x,
    spec = spec,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    B = 3,
    seed = 75,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    observed_theta_hat = list(mu = c(0, 0, 1), kappa = 2, psi = 0.7),
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    control = list(derivative_mc_size = 700L, derivative_mc_seed = 76L)
  )
  expect_identical(result$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_false(isTRUE(result$diagnostics$fallback_to_reestimated))
})

test_that("spherical Cauchy rho near zero uses the recorded slow fallback", {
  x <- r_sph_spherical_cauchy(12, mu = c(0, 0, 1), rho = 0.35)
  spec <- make_spherical_cauchy_spec(distance_type = "geodesic")
  result <- multiplier_bootstrap_gof(
    data = x,
    spec = spec,
    null = list(type = "composite"),
    statistics = "ks",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    B = 3,
    seed = 77,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    observed_theta_hat = list(mu = c(0, 0, 1), rho = 0),
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE),
    control = list(derivative_mc_size = 700L, derivative_mc_seed = 78L)
  )
  expect_identical(result$diagnostics$effective_bootstrap_method, "reestimated")
  expect_true(isTRUE(result$diagnostics$fallback_to_reestimated))
  expect_identical(result$diagnostics$fallback_reason, "spherical_cauchy_rho_zero_nonidentification")
})

test_that("an interior singular Vhat still errors instead of falling back", {
  x <- c(-1.4, -0.2, 0.1, 1.0, 1.8, 2.4)
  spec <- make_normal_spec(unknown_param = "both")
  theta_hat <- spec$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "composite"),
    control = list()
  )
  par0 <- c(theta_hat$mu, log(theta_hat$sigma))

  score_matrix_fn <- function(sample, par) {
    mu <- par[[1L]]
    sigma <- exp(par[[2L]])
    cbind(
      (sample - mu) / sigma^2,
      -1 + ((sample - mu)^2) / sigma^2
    )
  }

  sample_fn <- function(n_aux, par) {
    stats::rnorm(n_aux, mean = par[[1L]], sd = exp(par[[2L]]))
  }

  expect_error(
    prepare_fast_multiplier_score_model(
      spec = spec,
      data = x,
      theta_hat = theta_hat,
      ks_prep = NULL,
      cvm_prep = NULL,
      control = list(derivative_mc_size = 400L, derivative_mc_seed = 79L),
      par0 = par0,
      score_matrix_fn = score_matrix_fn,
      sample_fn = sample_fn,
      vhat_fn = function(data, par0, S_obs, aux_sample, Psi_aux) {
        matrix(c(1, 1, 1, 1), nrow = 2)
      }
    ),
    "singular or indefinite `Vhat`"
  )
})

test_that("beta mixture weighted MLE respects a user-supplied shape lower bound", {
  x <- r_sph_beta_mixture2(
    n = 40,
    mu = c(0, 0, 1),
    weight1 = 0.4,
    alpha1 = 2,
    beta1 = 8,
    alpha2 = 8,
    beta2 = 2
  )
  theta_hat <- fit_beta_mixture2_theta(
    data = x,
    null = list(type = "composite"),
    control = list(
      beta_mixture2_shape_lower = 1.001,
      beta_mixture2_n_starts = 3L,
      beta_mixture2_optim_control = list(maxit = 80L, reltol = 1e-6)
    )
  )

  expect_gte(min(theta_hat$alpha1, theta_hat$beta1, theta_hat$alpha2, theta_hat$beta2), 1.001 - 1e-10)
})

test_that("fast bootstrap chunking preserves fast multiplier results", {
  x <- rotasym::r_vMF(12, mu = c(1, 0, 0), kappa = 2)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(4, dim = 3),
    t_grid = c(0.4, 0.8)
  )

  result_full <- multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 12,
    seed = 81,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    control = list(
      derivative_mc_size = 700L,
      derivative_mc_seed = 82L
    )
  )

  result_chunked <- multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 12,
    seed = 81,
    n_cores = 1,
    bootstrap_method = "fast_multiplier",
    control = list(
      derivative_mc_size = 700L,
      derivative_mc_seed = 82L,
      fast_bootstrap_chunk_size = 5L
    )
  )

  expect_equal(result_full$bootstrap$statistics$ks, result_chunked$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_full$bootstrap$statistics$cvm, result_chunked$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_equal(result_full$inference$ks$p_value, result_chunked$inference$ks$p_value, tolerance = 1e-12)
  expect_equal(result_full$inference$cvm$p_value, result_chunked$inference$cvm$p_value, tolerance = 1e-12)
})

test_that("cardioid comet runner does not call fast multiplier for simple nulls", {
  source(file.path("scripts", "run_comets_distance_profile_cardioid.R"), local = FALSE)
  model <- make_comet_cardioid_models("Uniform")[[1L]]
  x <- rotasym::r_unif_sphere(n = 12, p = 3)
  result <- run_single_cardioid_comet_model(
    data_matrix = x,
    model = model,
    statistic = "ks",
    B = 3,
    n_cores = 1,
    seed = 72,
    bootstrap_method = "fast_multiplier",
    ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
    control = list()
  )

  expect_identical(result$diagnostics$requested_bootstrap_method, "fast_multiplier")
  expect_identical(result$diagnostics$bootstrap_method, "reestimated")
  expect_identical(result$diagnostics$effective_bootstrap_method, "reestimated")
  expect_identical(result$diagnostics$null_type, "simple")
})

test_that("still-unsupported fast multiplier models fail clearly", {
  x_vmf <- rotasym::r_vMF(6, mu = c(1, 0, 0), kappa = 2)
  expect_error(
    multiplier_bootstrap_logitnormal_mixture2(
      data = x_vmf,
      null = list(type = "composite"),
      statistics = "ks",
      ks_grid = list(omega_grid = generate_canonical_lattice(3, dim = 3), t_grid = c(0.4, 0.8)),
      B = 2,
      seed = 1,
      bootstrap_method = "fast_multiplier"
    ),
    "does not expose the fast multiplier preparation hook"
  )
})
