library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

strip_backend_timings <- function(result) {
  timing_names <- grep("_seconds$", names(result$diagnostics), value = TRUE)
  result$diagnostics[c(
    timing_names,
    "distance_profile_backend_requested",
    "distance_profile_backend_effective"
  )] <- NULL
  result
}

capture_profile_outcome <- function(expr) {
  warnings <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(conditionMessage(e), class = "captured_error")
  )
  list(value = value, warnings = warnings)
}

test_that("the default and explicit R backends do not load C++", {
  expect_false(distance_profile_cpp_is_loaded())
  omitted <- theoretical_distance_profile_normal(
    omega = c(-1, 0, 1), mu = 0, sigma = 1, t_values = c(0, 0.5, 2)
  )
  explicit <- theoretical_distance_profile_normal(
    omega = c(-1, 0, 1), mu = 0, sigma = 1, t_values = c(0, 0.5, 2),
    backend = "r"
  )
  expect_identical(omitted, explicit)
  expect_false(distance_profile_cpp_is_loaded())

  bootstrap_args <- list(
    data = c(-1, -0.2, 0.4, 1.1),
    null = list(type = "simple", theta = list(mu = 0, sigma = 1)),
    statistics = c("ks", "cvm"),
    ks_grid = list(omega_grid = c(-1, 0, 1), t_grid = c(0, 0.5, 2)),
    B = 2L,
    seed = 400L,
    n_cores = 1L
  )
  omitted_result <- do.call(multiplier_bootstrap_normal, bootstrap_args)
  explicit_result <- do.call(multiplier_bootstrap_normal, c(
    bootstrap_args, list(distance_profile_backend = "r")
  ))
  expect_identical(
    strip_backend_timings(omitted_result),
    strip_backend_timings(explicit_result)
  )
  expect_false(distance_profile_cpp_is_loaded())
})

test_that("1,000 reproducible normal profile cases are bitwise identical", {
  set.seed(2026072201)
  boundary_t <- c(-Inf, -1, 0, .Machine$double.eps, 1, 1e6, Inf)
  boundary_x <- c(-1e6, -10, -1, 0, 1, 10, 1e6)
  boundary_sigma <- c(1e-12, 0.1, 1, 10, 1e6)

  for (case_index in seq_len(1000L)) {
    n_omega <- sample.int(7L, 1L)
    n_t <- sample.int(7L, 1L)
    omega <- if (case_index <= length(boundary_x)) {
      rep(boundary_x[[case_index]], n_omega)
    } else {
      rnorm(n_omega, sd = sample(c(0.1, 1, 100), 1L))
    }
    t_values <- if (case_index <= length(boundary_t)) {
      rep(boundary_t[[case_index]], n_t)
    } else {
      sample(c(runif(n_t, -1, 5), boundary_t), n_t, replace = TRUE)
    }
    mu <- if (case_index %% 137L == 0L) boundary_x[[1L + case_index %% length(boundary_x)]] else rnorm(1L)
    sigma <- boundary_sigma[[1L + case_index %% length(boundary_sigma)]]

    r_value <- theoretical_distance_profile_normal(
      omega, mu, sigma, t_values, backend = "r"
    )
    cpp_value <- theoretical_distance_profile_normal(
      omega, mu, sigma, t_values, backend = "cpp"
    )
    expect_identical(cpp_value, r_value, info = paste("normal case", case_index))
  }
  expect_true(distance_profile_cpp_is_loaded())
})

test_that("normal matrix kernels preserve values, dimensions, and attributes", {
  set.seed(2026072202)
  spec <- make_normal_spec(unknown_param = "both")
  for (case_index in seq_len(50L)) {
    n <- sample(2:35, 1L)
    x <- rnorm(n, mean = runif(1L, -3, 3), sd = exp(runif(1L, -3, 3)))
    theta <- list(mu = rnorm(1L), sigma = exp(runif(1L, -4, 4)))
    distances <- spec$distance_matrix(x, x, list())
    r_value <- with_distance_profile_backend(
      "r",
      compute_theoretical_sample_profile_matrix(spec, x, distances, theta, list())
    )
    cpp_value <- with_distance_profile_backend(
      "cpp",
      compute_theoretical_sample_profile_matrix(spec, x, distances, theta, list())
    )
    expect_identical(cpp_value, r_value, info = paste("normal matrix case", case_index))
  }
})

test_that("weighted small-circle scalar and grid profiles are bitwise identical", {
  set.seed(2026072203)
  for (case_index in seq_len(30L)) {
    mu <- rnorm(3L)
    mu <- mu / sqrt(sum(mu^2))
    omega <- if (case_index %% 10L == 0L) {
      if (case_index %% 20L == 0L) -mu else mu
    } else {
      value <- rnorm(3L)
      value / sqrt(sum(value^2))
    }
    distance_type <- if (case_index %% 2L) "geodesic" else "chordal"
    boundary <- if (identical(distance_type, "geodesic")) pi else 2
    args <- list(
      omega = omega,
      t_values = c(0, boundary * runif(3L), boundary),
      mu = mu,
      pi = c(0.01, 0.25, 0.5, 0.75, 0.99)[[1L + case_index %% 5L]],
      kappa1 = c(1e-6, 0.2, 3, 25)[[1L + case_index %% 4L]],
      nu1 = c(0, 0.1, 0.4, 0.8, 0.999)[[1L + case_index %% 5L]],
      kappa2 = c(1e-6, 0.5, 8, 40)[[1L + (case_index + 1L) %% 4L]],
      nu2 = c(0, 0.2, 0.5, 0.9, 0.999)[[1L + (case_index + 2L) %% 5L]],
      distance_type = distance_type,
      method = "legendre",
      l_max = 40L,
      quad_n = 121L
    )
    r_value <- do.call(distance_profile_small_circle_weighted_mixture2, c(args, list(backend = "r")))
    cpp_value <- do.call(distance_profile_small_circle_weighted_mixture2, c(args, list(backend = "cpp")))
    expect_identical(cpp_value, r_value, info = paste("weighted scalar case", case_index))
  }

  for (case_index in seq_len(10L)) {
    mu <- c(0, 0, 1)
    omega_grid <- generate_canonical_lattice(3L + case_index %% 5L, dim = 3L)
    args <- list(
      omega_grid = omega_grid,
      mu = mu,
      pi = case_index / 11,
      kappa1 = 2 + case_index,
      nu1 = case_index / 11,
      kappa2 = 1 + case_index / 2,
      nu2 = 1 - case_index / 11,
      t_grid = c(0, seq(0.2, 2.8, length.out = 5L), pi),
      distance_type = "geodesic",
      method = "legendre",
      l_max = 40L,
      quad_n = 121L
    )
    r_value <- do.call(distance_profile_small_circle_weighted_mixture2_grid, c(args, list(backend = "r")))
    cpp_value <- do.call(distance_profile_small_circle_weighted_mixture2_grid, c(args, list(backend = "cpp")))
    expect_identical(cpp_value, r_value, info = paste("weighted grid case", case_index))
  }
})

test_that("invalid inputs have identical errors and warnings", {
  invalid_normal <- list(
    list(omega = 0, mu = 0, sigma = 1, t_values = NA_real_),
    list(omega = NA_real_, mu = 0, sigma = 1, t_values = 1),
    list(omega = 0, mu = NA_real_, sigma = 1, t_values = 1),
    list(omega = 0, mu = 0, sigma = NA_real_, t_values = 1),
    list(omega = Inf, mu = 0, sigma = 1, t_values = Inf),
    list(omega = 1:3, mu = c(0, 1), sigma = 1, t_values = 1:5),
    list(omega = numeric(0), mu = 0, sigma = 1, t_values = numeric(0))
  )
  for (args in invalid_normal) {
    r_outcome <- capture_profile_outcome(do.call(
      theoretical_distance_profile_normal, c(args, list(backend = "r"))
    ))
    cpp_outcome <- capture_profile_outcome(do.call(
      theoretical_distance_profile_normal, c(args, list(backend = "cpp"))
    ))
    expect_identical(cpp_outcome, r_outcome)
  }

  bad_weighted <- list(
    omega = c(1, 0, 0), t_values = 1, mu = c(0, 0, 0), pi = 0.5,
    kappa1 = 2, nu1 = 0.3, kappa2 = 3, nu2 = 0.4
  )
  expect_identical(
    capture_profile_outcome(do.call(
      distance_profile_small_circle_weighted_mixture2,
      c(bad_weighted, list(backend = "cpp"))
    )),
    capture_profile_outcome(do.call(
      distance_profile_small_circle_weighted_mixture2,
      c(bad_weighted, list(backend = "r"))
    ))
  )
  expect_error(
    theoretical_distance_profile_normal(0, 0, 1, 1, backend = "fortran"),
    "either 'r' or 'cpp'",
    fixed = TRUE
  )
})

test_that("weighted empirical profiles preserve ties bitwise", {
  set.seed(2026072208)
  for (case_index in seq_len(25L)) {
    n <- sample(3:40, 1L)
    distance_matrix <- matrix(
      sample(c(0, 0.25, 0.5, 1, 2), n * n, replace = TRUE),
      nrow = n,
      ncol = n
    )
    order_matrix <- t(vapply(seq_len(n), function(i) {
      as.integer(order(distance_matrix[i, ]))
    }, integer(n)))
    rank_matrix <- t(vapply(seq_len(n), function(i) {
      as.integer(rank(distance_matrix[i, ], ties.method = "max"))
    }, integer(n)))
    row_index_matrix <- matrix(rep.int(seq_len(n), n), nrow = n, ncol = n)
    rank_linear_index <- row_index_matrix + (rank_matrix - 1L) * n
    weights <- rexp(n)
    weights <- weights / mean(weights)

    r_value <- with_distance_profile_backend(
      "r",
      compute_weighted_sample_profile_matrix(
        order_matrix = order_matrix,
        rank_linear_index = rank_linear_index,
        normalized_weights = weights
      )
    )
    cpp_value <- with_distance_profile_backend(
      "cpp",
      compute_weighted_sample_profile_matrix(
        order_matrix = order_matrix,
        rank_linear_index = rank_linear_index,
        normalized_weights = weights
      )
    )
    expect_identical(cpp_value, r_value, info = paste("tie case", case_index))
  }
})

test_that("normal complete tests are bitwise identical for simple and composite nulls", {
  set.seed(2026072204)
  x <- rnorm(32L, mean = 0.7, sd = 1.4)
  ks_grid <- list(
    omega_grid = seq(-3, 4, length.out = 9L),
    t_grid = seq(0, 5, length.out = 11L)
  )
  common <- list(
    data = x,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 12L,
    seed = 401L,
    n_cores = 1L,
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    )
  )

  for (null in list(
    list(type = "simple", theta = list(mu = 0.7, sigma = 1.4)),
    list(type = "composite")
  )) {
    unknown <- if (identical(null$type, "simple")) NULL else "both"
    r_result <- do.call(multiplier_bootstrap_normal, c(
      common,
      list(null = null, unknown_param = unknown, distance_profile_backend = "r")
    ))
    cpp_result <- do.call(multiplier_bootstrap_normal, c(
      common,
      list(null = null, unknown_param = unknown, distance_profile_backend = "cpp")
    ))
    expect_identical(strip_backend_timings(cpp_result), strip_backend_timings(r_result))
    expect_identical(cpp_result$diagnostics$distance_profile_backend_requested, "cpp")
    expect_identical(cpp_result$diagnostics$distance_profile_backend_effective, "cpp")
  }
})

test_that("weighted small-circle complete composite tests are bitwise identical", {
  set.seed(2026072205)
  mu <- c(0, 0, 1)
  x <- r_sph_small_circle_weighted_mixture2(24L, mu, 0.6, 10, 0.45, 8, 0.25)
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = list(
      omega_grid = generate_canonical_lattice(5L, dim = 3L),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 5L)
    ),
    B = 3L,
    seed = 402L,
    n_cores = 1L,
    distance_type = "geodesic",
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = list(
      small_circle_weighted_mixture2_L_max = 80L,
      small_circle_weighted_mixture2_quad_n = 201L,
      small_circle_weighted_mixture2_n_starts = 1L,
      small_circle_weighted_mixture2_optim_control = list(maxit = 80L, reltol = 1e-6)
    )
  )
  r_result <- do.call(multiplier_bootstrap_small_circle_weighted_mixture2, c(
    common, list(distance_profile_backend = "r")
  ))
  cpp_result <- do.call(multiplier_bootstrap_small_circle_weighted_mixture2, c(
    common, list(distance_profile_backend = "cpp")
  ))
  expect_identical(strip_backend_timings(cpp_result), strip_backend_timings(r_result))
})

test_that("parallel and sequential C++ execution preserve exact results", {
  skip_on_os("windows")
  set.seed(2026072206)
  x <- rnorm(28L, mean = -0.2, sd = 0.8)
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = list(
      omega_grid = seq(-2, 2, length.out = 7L),
      t_grid = seq(0, 3, length.out = 8L)
    ),
    B = 8L,
    seed = 403L,
    unknown_param = "both",
    keep = list(observed_process = TRUE, bootstrap_statistics = TRUE, bootstrap_thetas = TRUE)
  )
  r_parallel <- do.call(multiplier_bootstrap_normal, c(
    common, list(n_cores = 2L, distance_profile_backend = "r")
  ))
  cpp_parallel <- do.call(multiplier_bootstrap_normal, c(
    common, list(n_cores = 2L, distance_profile_backend = "cpp")
  ))
  cpp_sequential <- do.call(multiplier_bootstrap_normal, c(
    common, list(n_cores = 1L, distance_profile_backend = "cpp")
  ))
  r_parallel$diagnostics$n_cores <- 1L
  cpp_parallel$diagnostics$n_cores <- 1L
  expect_identical(strip_backend_timings(cpp_parallel), strip_backend_timings(r_parallel))
  expect_identical(strip_backend_timings(cpp_parallel), strip_backend_timings(cpp_sequential))
})

test_that("failed families reject C++ and Jones-Pewsey remains untouched", {
  rejected_specs <- c(
    "mvnormal_euclidean",
    "logistic_gaussian_aitchison",
    "vmf_geodesic",
    "hvmf_geodesic",
    "spherical_cauchy_geodesic",
    "small_circle_geodesic",
    "watson_geodesic",
    "cardioid_geodesic",
    "beta_mixture2_geodesic",
    "uniform_beta_mixture_geodesic",
    "logitnormal_mixture2_geodesic",
    "small_circle_symmetric_mixture2_geodesic",
    "axial_truncnorm_mixture2_euclidean"
  )
  for (spec_name in rejected_specs) {
    expect_error(
      multiplier_bootstrap_gof(
        data = 0,
        spec = list(name = spec_name),
        null = list(type = "simple"),
        B = 1L,
        distance_profile_backend = "cpp"
      ),
      "did not pass the exactness and end-to-end performance gates",
      fixed = TRUE,
      info = spec_name
    )
  }
  expect_error(
    theoretical_distance_profile_vmf(
      omega = c(0, 0, 1), t_values = 1, mu = c(0, 0, 1), kappa = 1,
      backend = "cpp"
    ),
    "was not retained",
    fixed = TRUE
  )

  jp_profile_names <- ls(pattern = "^distance_profile_jp")
  expect_false("distance_profile_backend" %in% names(formals(multiplier_bootstrap_jp)))
  expect_true(all(vapply(jp_profile_names, function(name) {
    !"backend" %in% names(formals(get(name))) &&
      !isTRUE(attr(get(name), "distance_profile_backend_wrapper"))
  }, logical(1))))
  set.seed(2026072207)
  jp_mu <- c(0, 0, 1)
  jp_result <- multiplier_bootstrap_jp(
    data = r_sph_jp(8L, mu = jp_mu, kappa = 1, psi = 0.5),
    null = list(type = "simple", theta = list(mu = jp_mu, kappa = 1, psi = 0.5)),
    statistics = "ks",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(3L, dim = 3L),
      t_grid = c(0.3, 1)
    ),
    B = 1L,
    seed = 404L,
    n_cores = 1L,
    distance_type = "geodesic"
  )
  expect_false(any(c(
    "distance_profile_backend_requested",
    "distance_profile_backend_effective"
  ) %in% names(jp_result$diagnostics)))
  expect_error(
    multiplier_bootstrap_gof(
      data = 0,
      spec = list(name = "jp_geodesic"),
      null = list(type = "simple"),
      B = 1L,
      distance_profile_backend = "cpp"
    ),
    "intentionally excluded",
    fixed = TRUE
  )
})
