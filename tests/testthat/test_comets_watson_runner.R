library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)
source(file.path("scripts", "run_comets_watson_short_long_fast.R"))

test_that("Watson comet runner produces short/long KS/CvM outputs", {
  output_root <- file.path(tempdir(), "comets_watson_fast_runner_smoke")
  summary <- run_comets_watson_short_long_fast(
    output_root = output_root, B = 2L, n_cores = 1L,
    control = list(watson_L_max = 80L, watson_quad_n = 160L, derivative_mc_size = 300L,
                   fast_multiplier_cvm_block_size = 64L, fast_multiplier_ks_block_size = 64L)
  )
  expect_equal(nrow(summary), 4L)
  expect_equal(sort(unique(summary$dataset)), c("long_period", "short_period"))
  expect_equal(sort(unique(summary$statistic)), c("cvm", "ks"))
  expect_true(all(summary$effective_bootstrap_method == "fast_multiplier"))
  expect_true(file.exists(file.path(output_root, "comets_watson_short_long_fast_summary.csv")))
})
