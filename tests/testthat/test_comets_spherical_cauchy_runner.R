library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("scripts", "run_comets_spherical_cauchy_short_long_fast.R"))

test_that("comets spherical Cauchy fast runner produces short/long KS/CvM outputs", {
  output_root <- file.path(tempdir(), "comets_sc_fast_runner_smoke")

  summary_df <- run_comets_spherical_cauchy_short_long_fast(
    output_root = output_root,
    B = 2L,
    n_cores = 1L,
    ks_t_points = 50L,
    base_seed = 20260708L,
    control = list(
      spherical_cauchy_maxit = 100L,
      spherical_cauchy_reltol = 1e-8,
      spherical_cauchy_optim_method = "BFGS",
      spherical_cauchy_use_gradient = TRUE,
      spherical_cauchy_profile_tol = 1e-8,
      spherical_cauchy_profile_warn = FALSE,
      fast_multiplier_cvm_block_size = 64L,
      fast_multiplier_ks_block_size = 64L
    )
  )

  expect_equal(nrow(summary_df), 4L)
  expect_equal(sort(unique(summary_df$dataset)), c("long_period", "short_period"))
  expect_equal(sort(unique(summary_df$statistic)), c("cvm", "ks"))
  expect_true(all(summary_df$bootstrap_method == "fast_multiplier"))
  expect_true(all(summary_df$effective_bootstrap_method == "fast_multiplier"))
  expect_true(file.exists(file.path(output_root, "comets_spherical_cauchy_short_long_fast_summary.csv")))
})
