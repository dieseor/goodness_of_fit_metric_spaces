library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "calibration_study.R"))

make_test_mvnormal_ks_grid <- function(mu, Sigma) {
  theta <- normalize_mvnormal_theta(list(mu = mu, Sigma = Sigma))
  eig <- eigen(theta$Sigma, symmetric = TRUE)
  omega_grid <- rbind(
    theta$mu,
    theta$mu + sqrt(eig$values[[1L]]) * eig$vectors[, 1L],
    theta$mu - sqrt(eig$values[[1L]]) * eig$vectors[, 1L],
    theta$mu + sqrt(eig$values[[2L]]) * eig$vectors[, 2L]
  )

  list(
    omega_grid = omega_grid,
    t_grid = seq(0, 6, length.out = 6L)
  )
}

test_that("mvnormal weighted fit agrees with replicated-data fit", {
  set.seed(1201)
  x <- mvtnorm::rmvnorm(
    n = 10,
    mean = c(0.5, -0.25, 0.75),
    sigma = matrix(c(
      1.2, 0.2, 0.1,
      0.2, 0.8, 0.15,
      0.1, 0.15, 0.6
    ), nrow = 3L, byrow = TRUE)
  )
  weights <- c(1, 2, 1, 3, 2, 1, 2, 1, 2, 1)

  fit_weighted <- fit_mvnormal_theta(
    data = x,
    weights = weights,
    null = list(type = "composite"),
    unknown_param = "both"
  )
  x_replicated <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]
  fit_replicated <- fit_mvnormal_theta(
    data = x_replicated,
    null = list(type = "composite"),
    unknown_param = "both"
  )

  expect_equal(fit_weighted$mu, fit_replicated$mu, tolerance = 1e-12)
  expect_equal(fit_weighted$Sigma, fit_replicated$Sigma, tolerance = 1e-12)
})

test_that("mvnormal bootstrap supports simple and composite nulls", {
  set.seed(1202)
  mu <- c(0.2, -0.4, 0.6)
  Sigma <- matrix(c(
    1.0, 0.25, 0.10,
    0.25, 0.8, 0.20,
    0.10, 0.20, 0.7
  ), nrow = 3L, byrow = TRUE)
  x <- mvtnorm::rmvnorm(18, mean = mu, sigma = Sigma)
  ks_grid <- make_test_mvnormal_ks_grid(mu, Sigma)

  simple_null <- list(type = "simple", theta = list(mu = mu, Sigma = Sigma))
  result_1 <- multiplier_bootstrap_mvnormal(
    data = x,
    null = simple_null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1
  )
  result_2 <- multiplier_bootstrap_mvnormal(
    data = x,
    null = simple_null,
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 8,
    seed = 42,
    n_cores = 1
  )

  expect_equal(result_1$bootstrap$statistics$ks, result_2$bootstrap$statistics$ks, tolerance = 1e-12)
  expect_equal(result_1$bootstrap$statistics$cvm, result_2$bootstrap$statistics$cvm, tolerance = 1e-12)
  expect_true(result_1$inference$ks$p_value >= 0 && result_1$inference$ks$p_value <= 1)
  expect_true(result_1$inference$cvm$p_value >= 0 && result_1$inference$cvm$p_value <= 1)

  composite_result <- multiplier_bootstrap_mvnormal(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = 6,
    seed = 84,
    n_cores = 1,
    unknown_param = "both"
  )

  expect_true(composite_result$inference$ks$p_value >= 0 && composite_result$inference$ks$p_value <= 1)
  expect_true(composite_result$inference$cvm$p_value >= 0 && composite_result$inference$cvm$p_value <= 1)
  expect_equal(length(composite_result$observed$theta_hat$mu), 3L)
  expect_equal(dim(composite_result$observed$theta_hat$Sigma), c(3L, 3L))
})

test_that("mvnormal fast preparation uses Gaussian MLE influence correction", {
  set.seed(1203)
  x <- mvtnorm::rmvnorm(
    n = 14,
    mean = c(0.3, -0.2, 0.1),
    sigma = matrix(c(
      1.1, 0.20, 0.05,
      0.20, 0.9, 0.10,
      0.05, 0.10, 0.7
    ), nrow = 3L, byrow = TRUE)
  )

  spec <- make_mvnormal_spec(unknown_param = "both")
  theta_hat <- spec$fit_theta(
    data = x,
    weights = NULL,
    null = list(type = "composite"),
    control = list()
  )
  ks_grid <- make_test_mvnormal_ks_grid(theta_hat$mu, theta_hat$Sigma)
  ks_prep <- prepare_ks_observed_data(x, spec, theta_hat, ks_grid)
  cvm_prep <- prepare_cvm_observed_data(x, spec, theta_hat, control = list())
  prep <- spec_fast_multiplier_prepare(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = list(
      derivative_mc_size = 600L,
      derivative_mc_seed = 1204L,
      fast_multiplier_store_paper_quantities = TRUE
    )
  )

  d <- length(theta_hat$mu)
  p_sigma <- d * (d + 1L) / 2L
  centered_obs <- sweep(normalize_mvnormal_data(x), 2L, theta_hat$mu, FUN = "-")
  expected_if <- cbind(
    centered_obs,
    t(vapply(seq_len(nrow(centered_obs)), function(i) {
      rr <- centered_obs[i, , drop = FALSE]
      fast_multiplier_vech(crossprod(rr) - theta_hat$Sigma)
    }, numeric(p_sigma)))
  )

  expect_identical(prep$vhat_method, "fitted_gaussian_influence_reparameterization")
  expect_identical(prep$paper_Vhat_method, "analytic_expected_score_jacobian")
  expect_equal(dim(prep$S_obs), c(nrow(x), d + p_sigma))
  expect_equal(dim(prep$Psi_aux), c(600, d + p_sigma))
  expect_equal(dim(prep$Vhat), c(d + p_sigma, d + p_sigma))
  expect_equal(dim(prep$D_ks), c(nrow(ks_grid$omega_grid) * length(ks_grid$t_grid), d + p_sigma))
  expect_equal(dim(prep$D_cvm), c(nrow(x) * nrow(x), d + p_sigma))
  expect_equal(prep$S_obs, expected_if, tolerance = 1e-12)
  expect_equal(prep$Vhat, diag(d + p_sigma), tolerance = 1e-12)
  expect_equal(
    prep$paper_influence_obs,
    -prep$paper_score_obs %*% t(solve(prep$paper_Vhat)),
    tolerance = 1e-10
  )
  expect_lt(prep$vhat_diagnostics$score_mean_aux_norm, 0.35)
})

test_that("mvnormal fast multiplier supports the sample-based KS grid", {
  set.seed(1205)
  x <- mvtnorm::rmvnorm(
    n = 12,
    mean = c(0.15, -0.1, 0.25),
    sigma = matrix(c(
      1.0, 0.15, 0.05,
      0.15, 0.85, 0.10,
      0.05, 0.10, 0.75
    ), nrow = 3L, byrow = TRUE)
  )

  result <- multiplier_bootstrap_mvnormal(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 4,
    seed = 1206,
    n_cores = 1,
    unknown_param = "both",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_mc_size = 250L,
      derivative_mc_seed = 1207L,
      fast_bootstrap_chunk_size = 2L
    ),
    bootstrap_method = "fast_multiplier"
  )

  expect_true(result$inference$ks$p_value >= 0 && result$inference$ks$p_value <= 1)
  expect_true(result$inference$cvm$p_value >= 0 && result$inference$cvm$p_value <= 1)
  expect_equal(length(result$bootstrap$statistics$ks), 4L)
  expect_equal(length(result$bootstrap$statistics$cvm), 4L)
  expect_identical(
    result$diagnostics$fast_multiplier_backend_requested,
    "cpp"
  )
  expect_identical(
    result$diagnostics$fast_multiplier_backend_effective,
    "cpp"
  )
  expect_true(result$diagnostics$fast_multiplier_fuse_ks_cvm_effective)
  expect_true(
    result$diagnostics$fast_multiplier_cache_corrections_effective
  )
})

test_that("mvnormal legacy C++ is bitwise identical to R and contiguous C++ preserves inference", {
  for (d in c(2L, 10L)) {
    set.seed(1300L + d)
    x <- mvtnorm::rmvnorm(
      n = 18L,
      mean = rep(0, d),
      sigma = diag(d)
    )
    common <- list(
      data = x,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      B = 11L,
      seed = 1400L + d,
      n_cores = 1L,
      unknown_param = "both",
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = list(
        derivative_mc_size = 120L,
        derivative_mc_seed = 1500L + d,
        fast_bootstrap_chunk_size = 4L
      ),
      bootstrap_method = "fast_multiplier"
    )

    r_fused <- do.call(
      multiplier_bootstrap_mvnormal,
      c(common, list(
        fast_multiplier_backend = "r",
        fuse_ks_cvm = TRUE,
        cache_block_corrections = "true"
      ))
    )
    cpp_fused <- do.call(
      multiplier_bootstrap_mvnormal,
      c(common, list(
        fast_multiplier_backend = "cpp",
        fuse_ks_cvm = TRUE,
        cache_block_corrections = "true"
      ))
    )
    legacy_common <- common
    legacy_common$control <- utils::modifyList(
      common$control,
      list(fast_multiplier_cpp_kernel = "legacy")
    )
    cpp_legacy <- do.call(
      multiplier_bootstrap_mvnormal,
      c(legacy_common, list(
        fast_multiplier_backend = "cpp",
        fuse_ks_cvm = TRUE,
        cache_block_corrections = "true"
      ))
    )
    cpp_uncached <- do.call(
      multiplier_bootstrap_mvnormal,
      c(common, list(
        fast_multiplier_backend = "cpp",
        fuse_ks_cvm = TRUE,
        cache_block_corrections = "false"
      ))
    )
    r_unfused <- do.call(
      multiplier_bootstrap_mvnormal,
      c(common, list(
        fast_multiplier_backend = "r",
        fuse_ks_cvm = FALSE,
        cache_block_corrections = "true"
      ))
    )

    expect_identical(
      cpp_uncached$bootstrap$statistics,
      cpp_fused$bootstrap$statistics
    )
    expect_identical(
      r_unfused$bootstrap$statistics,
      r_fused$bootstrap$statistics
    )
    expect_identical(cpp_legacy$bootstrap$statistics, r_fused$bootstrap$statistics)
    expect_identical(cpp_legacy$observed, r_fused$observed)
    expect_identical(cpp_legacy$inference, r_fused$inference)
    expect_identical(cpp_fused$observed, r_fused$observed)
    expect_equal(
      cpp_fused$bootstrap$statistics,
      cpp_legacy$bootstrap$statistics,
      tolerance = 1e-12
    )
    expect_identical(cpp_fused$inference, cpp_legacy$inference)
    expect_identical(
      cpp_fused$diagnostics$fast_multiplier_backend_effective,
      "cpp"
    )
    expect_true(
      cpp_fused$diagnostics$fast_multiplier_fuse_ks_cvm_effective
    )
    expect_identical(
      cpp_fused$diagnostics$fast_multiplier_cpp_kernel_effective,
      "contiguous_double"
    )
    expect_identical(
      cpp_legacy$diagnostics$fast_multiplier_cpp_kernel_effective,
      "legacy"
    )
    if (d == 2L) {
      parallel_common <- common
      parallel_common$n_cores <- 2L
      cpp_parallel <- do.call(
        multiplier_bootstrap_mvnormal,
        c(parallel_common, list(
          fast_multiplier_backend = "cpp",
          fuse_ks_cvm = TRUE,
          cache_block_corrections = "true"
        ))
      )
      expect_identical(
        cpp_parallel$bootstrap$statistics,
        cpp_fused$bootstrap$statistics
      )
      expect_identical(cpp_parallel$inference, cpp_fused$inference)
    }
  }
})

test_that("mvnormal fast defaults and cache cutoff are explicit", {
  defaults <- formals(multiplier_bootstrap_mvnormal)
  expect_identical(eval(defaults$fast_multiplier_backend), c("cpp", "r"))
  expect_identical(eval(defaults$fuse_ks_cvm), TRUE)
  expect_identical(
    eval(defaults$cache_block_corrections),
    c("auto", "true", "false")
  )

  at_cutoff <- resolve_fast_sample_correction_cache(
    n_centers = 500L,
    n_thresholds = 500L,
    n_parameters = 3L
  )
  above_cutoff <- resolve_fast_sample_correction_cache(
    n_centers = 501L,
    n_thresholds = 501L,
    n_parameters = 3L
  )
  forced <- resolve_fast_sample_correction_cache(
    n_centers = 501L,
    n_thresholds = 501L,
    n_parameters = 3L,
    control = list(fast_multiplier_cache_corrections = TRUE)
  )

  expect_true(at_cutoff$enabled)
  expect_false(above_cutoff$enabled)
  expect_true(forced$enabled)
  expect_identical(at_cutoff$n_max, 500L)
})
