library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "sunspots", "run_sunspots_cycle23_temporal_beta_composite_gof.R"))
source(file.path("real_data", "sunspots", "run_sunspots_cycle23_joint_spatial_window_kde_plots.R"))

make_joint_theta_for_tests <- function() {
  list(a_N = 0.55, b_N = -0.18, a_S = 0.50, b_S = -0.15, c = 14)
}

test_that("Rgof composite audit checks refit-style phat usage", {
  fake_gof_test <- function(x, vals, pnull, rnull, phat, B, doMethods, maxProcessor = 1L, ...) {
    phat(x)
    for (i in seq_len(B)) {
      x_star <- rnull(phat(x) + c(i * 1e-3, 0))
      phat(x_star)
    }
    list(
      statistics = c(KS = 0.1, CvM = 0.2, AD = 0.3),
      p.values = c(KS = 0.2, CvM = 0.3, AD = 0.4)
    )
  }

  audit <- sunspots_temporal_rgof_composite_refit_audit(
    B = 9L,
    sample_size = 30L,
    seed = 123L,
    gof_test_fun = fake_gof_test
  )
  expect_true(audit$refit_pass)
  expect_gte(audit$phat_calls, audit$required_min_calls)
  expect_true(audit$unique_phat_means > 1L)
})

test_that("spatial windows use the requested non-cumulative rank bins", {
  s <- seq(0.01, 0.99, length.out = 20L)
  windows <- sunspots_joint_spatial_rank_windows(s)
  expect_equal(windows$summary$lower_rank_level, c(0, 0.10, 0.20, 0.40, 0.60))
  expect_equal(windows$summary$upper_rank_level, c(0.10, 0.20, 0.40, 0.60, 0.80))
  expect_equal(windows$summary$center_rank_level, c(0.05, 0.15, 0.30, 0.50, 0.70))
  expect_equal(windows$summary$n, c(2L, 2L, 4L, 4L, 4L))
})

test_that("parametric window density is conditional-only and finite", {
  theta <- make_joint_theta_for_tests()
  x <- jp_normalize_unit_matrix(matrix(c(
    0, 1, 0,
    0.2, 0.3, sqrt(1 - 0.2^2 - 0.3^2),
    -0.4, 0.1, sqrt(1 - 0.4^2 - 0.1^2)
  ), ncol = 3, byrow = TRUE), arg_name = "`x`", min_ncol = 3L)
  center_s <- 0.42
  from_helper <- sunspots_joint_parametric_conditional_density(x, center_s, theta)
  direct <- exp(sunspots_time_varying_log_density(x = x, u = rep(center_s, nrow(x)), theta = theta))
  expect_equal(from_helper, direct, tolerance = 1e-13)
  expect_true(all(is.finite(from_helper)))
})

test_that("HDR thresholds are monotone in credibility level", {
  density_values <- c(0.9, 0.8, 0.2, 0.1)
  area_weights <- c(0.1, 0.2, 0.3, 0.4)
  hdr <- sunspots_joint_hdr_thresholds(density_values, area_weights, levels = c(0.5, 0.8, 0.95))
  expect_true(all(diff(hdr$threshold) <= 0))
})

test_that("DirStats KDE integral is finite and close to one on S2", {
  skip_if_not_installed("DirStats")
  set.seed(1)
  x_window <- matrix(rnorm(120L), ncol = 3L)
  x_window <- x_window / sqrt(rowSums(x_window^2))
  bw <- sunspots_joint_select_bandwidth_lcv_emi(x_window = x_window, bandwidth_seed = 7L)
  integral <- sunspots_joint_lebedev_integral(function(x_eval) {
    DirStats::kde_dir(x = x_eval, data = x_window, h = bw$h, L = NULL)
  })
  expect_true(is.finite(integral))
  expect_equal(integral, 1, tolerance = 0.08)
})
