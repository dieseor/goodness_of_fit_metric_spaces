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
