library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)

model_specs_path <- if (file.exists(file.path("bootstrap", "model_specs.R"))) {
  file.path("bootstrap", "model_specs.R")
} else {
  file.path("..", "..", "bootstrap", "model_specs.R")
}
source(model_specs_path)

small_circle_symmetric_mixture2_model_spec_path <- if (file.exists(file.path("bootstrap", "small_circle_symmetric_mixture2_model_spec.R"))) {
  file.path("bootstrap", "small_circle_symmetric_mixture2_model_spec.R")
} else {
  file.path("..", "..", "bootstrap", "small_circle_symmetric_mixture2_model_spec.R")
}
source(small_circle_symmetric_mixture2_model_spec_path)

test_that("symmetric small-circle-mixture parameters canonicalize deterministically", {
  theta <- small_circle_symmetric_mixture2_canonicalize_theta(list(
    mu = c(-0.2, 0.3, -0.93),
    kappa = 7,
    nu = -0.4
  ))

  expect_gte(theta$mu[which(abs(theta$mu) > 1e-12)[1]], 0)
  expect_equal(theta$nu, 0.4, tolerance = 1e-12)
  expect_equal(sqrt(sum(theta$mu^2)), 1, tolerance = 1e-12)
  expect_equal(
    small_circle_symmetric_mixture2_normalize_theta(theta, ambient_dim = 3L)$ambient_dim,
    3L
  )
})

test_that("symmetric small-circle-mixture nu zero recovers the simple small-circle model", {
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)
  omega <- jp_normalize_unit_vector(c(0.3, 0.4, 0.866), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 17)

  observed <- distance_profile_small_circle_symmetric_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    kappa = 12,
    nu = 0,
    method = "legendre"
  )
  expected <- distance_profile_small_circle(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    kappa = 12,
    nu = 0,
    method = "legendre"
  )

  expect_equal(observed, expected, tolerance = 1e-10)
  expect_equal(
    small_circle_symmetric_mixture2_axis_density(seq(-0.8, 0.8, length.out = 9), kappa = 12, nu = 0),
    small_circle_axis_density(seq(-0.8, 0.8, length.out = 9), kappa = 12, nu = 0),
    tolerance = 1e-12
  )
})

test_that("symmetric small-circle-mixture kappa zero is uniform on S2", {
  z_grid <- seq(-1, 1, length.out = 11)
  expect_equal(
    small_circle_symmetric_mixture2_axis_density(z_grid, kappa = 0, nu = 0.4),
    rep(0.5, length(z_grid)),
    tolerance = 1e-12
  )
  expect_equal(
    small_circle_symmetric_mixture2_axis_cdf(z_grid, kappa = 0, nu = 0.4),
    (z_grid + 1) / 2,
    tolerance = 1e-12
  )

  mu <- c(0, 0, 1)
  omega <- c(1, 0, 0)
  t_grid <- seq(0, pi, length.out = 21)
  observed <- distance_profile_small_circle_symmetric_mixture2(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    kappa = 0,
    nu = 0.5,
    distance_type = "geodesic"
  )
  expect_equal(observed, (1 - cos(t_grid)) / 2, tolerance = 1e-12)
})

test_that("symmetric small-circle-mixture axial and S2 densities are numerically normalized", {
  axial_integral <- integrate(
    f = function(z) small_circle_symmetric_mixture2_axis_density(
      z = z,
      kappa = 12,
      nu = 0.35
    ),
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  s2_integral <- integrate(
    f = function(z) {
      x <- cbind(sqrt(pmax(0, 1 - z^2)), 0, z)
      2 * pi * d_sph_small_circle_symmetric_mixture2_s2(
        x = x,
        mu = c(0, 0, 1),
        kappa = 12,
        nu = 0.35
      )
    },
    lower = -1,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value

  expect_equal(axial_integral, 1, tolerance = 1e-8)
  expect_equal(s2_integral, 1, tolerance = 1e-8)
})

test_that("symmetric small-circle-mixture axis law is even and symmetric", {
  z_grid <- seq(-0.9, 0.9, length.out = 19)
  density_pos <- small_circle_symmetric_mixture2_axis_density(z_grid, kappa = 8, nu = 0.35)
  density_neg <- small_circle_symmetric_mixture2_axis_density(-z_grid, kappa = 8, nu = 0.35)
  cdf_pos <- small_circle_symmetric_mixture2_axis_cdf(z_grid, kappa = 8, nu = 0.35)
  cdf_neg <- small_circle_symmetric_mixture2_axis_cdf(-z_grid, kappa = 8, nu = 0.35)

  expect_equal(density_pos, density_neg, tolerance = 1e-12)
  expect_equal(cdf_pos, 1 - cdf_neg, tolerance = 1e-12)
})

test_that("symmetric small-circle-mixture profile equals the mean of the two components", {
  mu <- jp_normalize_unit_vector(c(0.25, -0.15, 0.956), arg_name = "mu", min_length = 3L)
  omega <- jp_normalize_unit_vector(c(-0.4, 0.2, 0.894), arg_name = "omega", min_length = 3L)
  t_grid <- seq(0.05, pi - 0.05, length.out = 21)

  for (method in c("legendre", "integral")) {
    observed <- distance_profile_small_circle_symmetric_mixture2(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      kappa = 10,
      nu = 0.55,
      method = method,
      l_max = 120L,
      quad_n = 400L
    )
    expected <- 0.5 * (
      distance_profile_small_circle(
        omega = omega,
        t_values = t_grid,
        mu = mu,
        kappa = 10,
        nu = 0.55,
        method = method,
        l_max = 120L,
        quad_n = 400L
      ) +
        distance_profile_small_circle(
          omega = omega,
          t_values = t_grid,
          mu = -mu,
          kappa = 10,
          nu = 0.55,
          method = method,
          l_max = 120L,
          quad_n = 400L
        )
    )
    expect_equal(observed, expected, tolerance = 1e-12)
  }
})

test_that("symmetric small-circle-mixture Legendre and integral profiles agree", {
  mu <- c(0, 0, 1)
  complement <- jp_orthonormal_complement(mu)
  omega_list <- list(
    mu,
    -mu,
    complement[, 1L],
    jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  )
  t_grid <- seq(0.05, pi - 0.05, length.out = 25)
  scenario_grid <- expand.grid(
    kappa = c(0, 5, 20),
    nu = c(0, 0.3, 0.7),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(scenario_grid))) {
    diffs <- vapply(omega_list, function(omega) {
      leg <- distance_profile_small_circle_symmetric_mixture2(
        omega = omega,
        t_values = t_grid,
        mu = mu,
        kappa = scenario_grid$kappa[[i]],
        nu = scenario_grid$nu[[i]],
        method = "legendre",
        l_max = 180L,
        quad_n = 500L
      )
      integ <- distance_profile_small_circle_symmetric_mixture2(
        omega = omega,
        t_values = t_grid,
        mu = mu,
        kappa = scenario_grid$kappa[[i]],
        nu = scenario_grid$nu[[i]],
        method = "integral",
        quad_n = 500L
      )
      max(abs(leg - integ))
    }, numeric(1))

    expect_lte(max(diffs), 5e-4)
  }
})

test_that("symmetric small-circle-mixture grid and CvM helpers match naive evaluation", {
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)
  omega_grid <- rbind(
    mu,
    -mu,
    jp_normalize_unit_vector(c(0.35, 0.4, 0.846), arg_name = "omega", min_length = 3L)
  )
  t_grid <- seq(0.1, pi - 0.1, length.out = 11)
  x <- r_sph_small_circle_symmetric_mixture2(n = 6, mu = mu, kappa = 12, nu = 0.65)
  dot_products <- pmin(pmax(x %*% t(x), -1), 1)

  for (method in c("legendre", "integral")) {
    grid_naive <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
      distance_profile_small_circle_symmetric_mixture2(
        omega = omega_grid[i, ],
        t_values = t_grid,
        mu = mu,
        kappa = 12,
        nu = 0.65,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
    }, numeric(length(t_grid))))
    grid_fast <- distance_profile_small_circle_symmetric_mixture2_grid(
      omega_grid = omega_grid,
      mu = mu,
      kappa = 12,
      nu = 0.65,
      t_grid = t_grid,
      distance_type = "geodesic",
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10
    )
    expect_equal(grid_fast, grid_naive, tolerance = 1e-12)

    cvm_naive <- t(vapply(seq_len(nrow(x)), function(i) {
      distance_profile_small_circle_symmetric_mixture2(
        omega = x[i, ],
        t_values = acos(dot_products[i, ]),
        mu = mu,
        kappa = 12,
        nu = 0.65,
        distance_type = "geodesic",
        method = method,
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
    }, numeric(nrow(x))))
    cvm_fast <- distance_profile_small_circle_symmetric_mixture2_cvm_grid(
      X = x,
      mu = mu,
      kappa = 12,
      nu = 0.65,
      method = method,
      l_max = 80L,
      quad_n = 250L,
      tol = 1e-10
    )
    expect_equal(cvm_fast, cvm_naive, tolerance = 1e-12)
  }
})

test_that("symmetric small-circle-mixture blockwise CvM profile matches dense route", {
  set.seed(20260603)
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)

  for (n in c(50, 100, 300)) {
    x <- r_sph_small_circle_symmetric_mixture2(n = n, mu = mu, kappa = 12, nu = 0.6)
    for (distance_type in c("geodesic", "chordal")) {
      distance_matrix <- {
        dot_products <- pmin(pmax(x %*% t(x), -1), 1)
        if (identical(distance_type, "geodesic")) acos(dot_products) else sqrt(2 * (1 - dot_products))
      }
      dot_matrix <- pmin(pmax(x %*% t(x), -1), 1)
      dense <- distance_profile_small_circle_symmetric_mixture2_cvm_grid(
        X = x,
        mu = mu,
        kappa = 12,
        nu = 0.6,
        distance_matrix = distance_matrix,
        distance_type = distance_type,
        method = "legendre",
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )

      block_rows <- split(seq_len(n), ceiling(seq_len(n) / 37L))
      blockwise <- do.call(rbind, lapply(block_rows, function(idx) {
        small_circle_symmetric_mixture2_cvm_profile_block(
          X_block = x[idx, , drop = FALSE],
          dot_threshold_block = dot_matrix[idx, , drop = FALSE],
          mu = mu,
          kappa = 12,
          nu = 0.6,
          distance_type = distance_type,
          method = "legendre",
          l_max = 80L,
          quad_n = 250L,
          tol = 1e-10
        )
      }))

      expect_equal(blockwise, dense, tolerance = 1e-8)
    }
  }
})

test_that("weighted sample profile block is correctly normalized", {
  dot_matrix <- matrix(c(
    1.0, 0.8, 0.1,
    0.8, 1.0, 0.4,
    0.1, 0.4, 1.0
  ), nrow = 3, byrow = TRUE)
  rank_matrix <- t(vapply(seq_len(nrow(dot_matrix)), function(i) {
    as.integer(rank(-dot_matrix[i, ], ties.method = "max"))
  }, integer(nrow(dot_matrix))))
  order_matrix <- t(vapply(seq_len(nrow(dot_matrix)), function(i) {
    as.integer(order(dot_matrix[i, ], decreasing = TRUE))
  }, integer(nrow(dot_matrix))))
  row_index_matrix <- matrix(rep.int(seq_len(nrow(dot_matrix)), nrow(dot_matrix)), nrow = nrow(dot_matrix), ncol = nrow(dot_matrix))
  rank_linear_index <- row_index_matrix + (rank_matrix - 1L) * nrow(dot_matrix)

  block_profile_ones <- compute_weighted_sample_profile_block(
    order_matrix_block = order_matrix,
    rank_linear_index_block = rank_linear_index,
    normalized_weights = rep(1, nrow(dot_matrix)),
    n_total = nrow(dot_matrix)
  )
  empirical_profile <- rank_matrix / nrow(dot_matrix)
  expect_equal(block_profile_ones, empirical_profile, tolerance = 1e-12)

  weights <- c(2, 5, 3)
  block_profile_weighted <- compute_weighted_sample_profile_block(
    order_matrix_block = order_matrix,
    rank_linear_index_block = rank_linear_index,
    normalized_weights = weights,
    n_total = nrow(dot_matrix)
  )
  naive_profile <- t(vapply(seq_len(nrow(dot_matrix)), function(i) {
    vapply(seq_len(ncol(dot_matrix)), function(j) {
      threshold <- dot_matrix[i, j]
      sum(weights[dot_matrix[i, ] >= threshold]) / sum(weights)
    }, numeric(1))
  }, numeric(ncol(dot_matrix))))

  expect_equal(block_profile_weighted, naive_profile, tolerance = 1e-12)
  expect_true(all(block_profile_weighted >= 0))
  expect_true(all(block_profile_weighted <= 1))
})

test_that("symmetric small-circle-mixture light CvM prep/bootstrap matches dense computation", {
  set.seed(20260604)
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)

  dense_weighted_profile <- function(dot_block, weights) {
    total_weight <- sum(weights)
    t(vapply(seq_len(nrow(dot_block)), function(i) {
      rank_i <- as.integer(rank(-dot_block[i, ], ties.method = "max"))
      vapply(seq_len(ncol(dot_block)), function(j) {
        sum(weights[dot_block[i, ] >= dot_block[i, j]]) / total_weight
      }, numeric(1))
    }, numeric(ncol(dot_block))))
  }

  for (n in c(50, 100, 300)) {
    x <- r_sph_small_circle_symmetric_mixture2(n = n, mu = mu, kappa = 12, nu = 0.6)
    theta_hat <- list(mu = mu, kappa = 12, nu = 0.6)
    theta_star <- list(mu = mu, kappa = 10, nu = 0.5)
    normalized_weights <- stats::rexp(n)
    normalized_weights <- normalized_weights / mean(normalized_weights)

    for (distance_type in c("geodesic", "chordal")) {
      prep <- small_circle_symmetric_mixture2_cvm_prepare(
        data = x,
        theta_hat = theta_hat,
        distance_type = distance_type,
        control = list(
          small_circle_symmetric_mixture2_profile_method = "legendre",
          small_circle_symmetric_mixture2_L_max = 80L,
          small_circle_symmetric_mixture2_quad_n = 250L,
          small_circle_symmetric_mixture2_tol = 1e-10,
          small_circle_symmetric_mixture2_cvm_block_size = 37L
        )
      )

      dot_matrix <- pmin(pmax(x %*% t(x), -1), 1)
      empirical_obs <- dense_weighted_profile(dot_matrix, rep(1, n))
      empirical_star <- dense_weighted_profile(dot_matrix, normalized_weights)
      theoretical_obs <- distance_profile_small_circle_symmetric_mixture2_cvm_grid(
        X = x,
        mu = theta_hat$mu,
        kappa = theta_hat$kappa,
        nu = theta_hat$nu,
        distance_type = distance_type,
        method = "legendre",
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )
      theoretical_star <- distance_profile_small_circle_symmetric_mixture2_cvm_grid(
        X = x,
        mu = theta_star$mu,
        kappa = theta_star$kappa,
        nu = theta_star$nu,
        distance_type = distance_type,
        method = "legendre",
        l_max = 80L,
        quad_n = 250L,
        tol = 1e-10
      )

      observed_dense <- mean((sqrt(n) * (empirical_obs - theoretical_obs))^2)
      bootstrap_dense <- mean((sqrt(n) * ((empirical_star - theoretical_star) - (empirical_obs - theoretical_obs)))^2)

      bootstrap_light <- small_circle_symmetric_mixture2_cvm_bootstrap_stat(
        data = x,
        normalized_weights = normalized_weights,
        theta_star = theta_star,
        cvm_prep = prep,
        null = list(type = "composite"),
        distance_type = distance_type,
        control = list(
          small_circle_symmetric_mixture2_profile_method = "legendre",
          small_circle_symmetric_mixture2_L_max = 80L,
          small_circle_symmetric_mixture2_quad_n = 250L,
          small_circle_symmetric_mixture2_tol = 1e-10
        ),
        scale_factor = 1
      )

      expect_equal(prep$statistic, observed_dense, tolerance = 1e-8)
      expect_equal(bootstrap_light, bootstrap_dense, tolerance = 1e-8)
    }
  }
})

test_that("symmetric small-circle-mixture sampler and weighted MLE behave coherently", {
  set.seed(20260603)
  mu <- jp_normalize_unit_vector(c(0.2, -0.3, 0.93), arg_name = "mu", min_length = 3L)
  sample <- r_sph_small_circle_symmetric_mixture2(n = 180, mu = mu, kappa = 18, nu = 0.55)
  expect_equal(dim(sample), c(180, 3))
  expect_lt(max(abs(sqrt(rowSums(sample^2)) - 1)), 1e-8)

  weights <- rep(c(1, 2, 3), length.out = nrow(sample))
  fit_weighted <- small_circle_symmetric_mixture2_mle_s2_weighted(
    sample,
    weights = weights,
    control = list(
      small_circle_symmetric_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )
  sample_rep <- sample[rep(seq_len(nrow(sample)), times = weights), , drop = FALSE]
  fit_rep <- small_circle_symmetric_mixture2_mle_s2_weighted(
    sample_rep,
    control = list(
      small_circle_symmetric_mixture2_start_theta = fit_weighted,
      small_circle_symmetric_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  expect_gt(abs(sum(fit_weighted$mu * fit_rep$mu)), 0.98)
  expect_equal(fit_weighted$kappa, fit_rep$kappa, tolerance = 1)
  expect_equal(fit_weighted$nu, fit_rep$nu, tolerance = 0.05)
})
