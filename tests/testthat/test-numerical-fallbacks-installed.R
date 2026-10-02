test_that("weighted sample profiles respect each center's distance order", {
  order <- rbind(1:3, 3:1)
  rank_index <- matrix(c(1L, 2L, 4L, 5L, 7L, 8L), nrow = 2L)
  result <- compute_weighted_sample_profile_block(order, rank_index,
                                                  normalized_weights = c(1, 2, 3),
                                                  n_total = 3L)
  expect_equal(result, rbind(c(1, 3, 6), c(3, 5, 6)) / 6)
  expect_error(compute_weighted_sample_profile_block(order, rank_index,
                                                     c(1, 2), 3L), "dimensions")
  expect_error(compute_weighted_sample_profile_block(order, rank_index,
                                                     c(1, -2, 3), 3L), "nonnegative")
})

test_that("the quadratic-form Monte Carlo fallback preserves the caller's RNG", {
  settings <- list(mc_label = "test", mc_abs_error = 1,
                   mc_batch_size = 1000L, mc_conf_level = 0.95,
                   mc_seed = 41L)
  set.seed(17)
  before <- .Random.seed
  result <- suppressMessages(mvnormal_quadform_mc_cdf(
    q = c(0.5, 1, 2),
    canonical = list(lambda = 1, h = 1, delta = 0),
    settings = settings
  ))
  expect_identical(.Random.seed, before)
  expect_identical(result$n, 1000L)
  expect_true(all(diff(result$probability) >= 0))
  expect_true(all(result$probability >= result$interval[, "lower"] &
                  result$probability <= result$interval[, "upper"]))
})

test_that("the Legendre cutoff respects uniform and capped cases", {
  expect_identical(select_vmf_s2_legendre_l_max(0, min_l = 3L,
                                                max_l = 8L), 3L)
  expect_identical(select_vmf_s2_legendre_l_max(10, min_l = 1L,
                                                max_l = 2L,
                                                tail_tol = 1e-30), 2L)
  expect_error(select_vmf_s2_legendre_l_max(-1), "nonnegative")
})

test_that("chunked parallel bootstrap returns finite inference with progress", {
  skip_on_os("windows")
  x <- c(-1.4, -0.8, -0.4, -0.1, 0.2, 0.5, 0.9, 1.3)
  fit <- suppressMessages(gof_test(
    x, "normal", theta = list(mu = 0, sigma = 1),
    bootstrap = "reestimated", statistics = "ks", B = 4,
    n_cores = 2, seed = 23,
    control = list(progress_bar = TRUE, reestimated_bootstrap_chunk_size = 1L)
  ))
  expect_true(is.finite(fit$inference$ks$p_value))
  expect_identical(fit$bootstrap$used, "reestimated")
})

test_that("fast bootstrap reports progress for chunked multipliers", {
  x <- c(-1.4, -0.8, -0.4, -0.1, 0.2, 0.5, 0.9, 1.3)
  fit <- suppressMessages(gof_test(
    x, "normal", bootstrap = "fast", statistics = "ks", B = 2,
    seed = 24,
    control = list(progress_bar = TRUE, fast_bootstrap_chunk_size = 1L)
  ))
  expect_true(is.finite(fit$inference$ks$p_value))
  expect_identical(fit$bootstrap$used, "fast")
})

test_that("S1 uniform angular CDF and sphere distances match exact values", {
  cdf <- build_vmf_s1_cdf(c(1, 0), kappa = 0, n_grid = 101L)
  expect_equal(cdf$cdf, seq(0, 1, length.out = 101L), tolerance = 1e-12)
  expect_equal(evaluate_vmf_s1_cdf(pi, cdf), 0.5, tolerance = 1e-12)
  expect_equal(sphere_distance(c(1, 0, 0), c(0, 1, 0)), sqrt(2))
  expect_equal(sphere_distance(c(1, 0, 0), c(0, 1, 0), "geodesic"), pi / 2)
  expect_error(sphere_distance(c(1, 0, 0), c(0, 1, 0), "invalid"), "distance_type")
})

test_that("HvMF tabulated projection CDF brackets the supplied grid", {
  tab <- hvmf_tabulate_projection_cdf_h2(alpha = 1.2, kappa = 3,
                                          y_max = 2, grid = c(1.2, 1.6))
  expect_equal(tab$y, c(1, 1.2, 1.6, 2))
  expect_identical(tab$cdf[[1L]], 0)
  expect_true(all(diff(tab$cdf) >= 0))
  expect_true(all(tab$cdf <= 1))
})

test_that("sorted weighted profiles assign one value to each tied distance", {
  result <- compute_sorted_weighted_profile_block(
    order_matrix = rbind(1:3, 3:1),
    sorted_distance_matrix = rbind(c(0, 1, 1), c(0, 0, 2)),
    centered_weights = c(1, 2, 3), row_indices = 1:2
  )
  expect_equal(result, rbind(c(1, 6, 6), c(5, 5, 6)))
})

test_that("memory diagnostics report vector and matrix shapes when enabled", {
  messages <- capture_messages(debug_memory_log(
    control = list(cvm_memory_debug = TRUE), label = "shape check",
    objects = list(vector = 1:3, matrix = matrix(1:4, 2))
  ))
  expect_match(paste(messages, collapse = "\n"), "vector: length=3", fixed = TRUE)
  expect_match(paste(messages, collapse = "\n"), "matrix: dim=2x2", fixed = TRUE)
})
