library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
source(file.path("bootstrap", "normal_sigma_Id_bootstrap.R"))

compare_reestimated_sample_fusion <- function(runner, x, extra = list(),
                                               null = list(type = "composite")) {
  common <- c(list(
    data = x,
    null = null,
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = 9L,
    alpha = 0.05,
    n_cores = 1L,
    seed = 20260837L,
    bootstrap_method = "reestimated",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    distance_profile_backend = "r"
  ), extra)
  legacy <- do.call(runner, c(common, list(
    control = list(reestimated_fuse_ks_cvm = FALSE)
  )))
  fused <- do.call(runner, c(common, list(
    control = list(reestimated_fuse_ks_cvm = TRUE)
  )))

  expect_false(legacy$diagnostics$reestimated_fuse_ks_cvm_effective)
  expect_true(fused$diagnostics$reestimated_fuse_ks_cvm_effective)
  expect_true(fused$diagnostics$shared_sample_ks_cvm_cache)
  expect_equal(legacy$observed, fused$observed, tolerance = 1e-13)
  expect_equal(legacy$bootstrap$statistics, fused$bootstrap$statistics, tolerance = 1e-13)
  expect_equal(legacy$inference, fused$inference, tolerance = 1e-13)
}

test_that("fused reestimated sample KS--CvM is identical for restricted single-spiked Normal", {
  set.seed(20260838L)
  x <- rrestricted_spiked_normal(32L, c(1 / sqrt(2), 1 / sqrt(2)), lambda = 2)
  compare_reestimated_sample_fusion(multiplier_bootstrap_restricted_spiked_normal, x)
})

test_that("fused reestimated sample KS--CvM also preserves the simple-null statistic", {
  set.seed(20260841L)
  x <- rnorm(28L, mean = 0.15, sd = 1.2)
  compare_reestimated_sample_fusion(
    multiplier_bootstrap_normal,
    x,
    extra = list(unknown_param = "both"),
    null = list(type = "simple", theta = list(mu = 0.15, sigma = 1.2))
  )
})

test_that("fixed KS grids retain the legacy reestimated route", {
  set.seed(20260842L)
  x <- rnorm(24L)
  fit <- multiplier_bootstrap_normal(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = list(omega_grid = seq(-1, 1, length.out = 3L), t_grid = seq(0, 2, length.out = 4L)),
    B = 5L,
    seed = 20260843L,
    unknown_param = "both",
    bootstrap_method = "reestimated",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE),
    control = list(reestimated_fuse_ks_cvm = TRUE)
  )
  expect_false(fit$diagnostics$reestimated_fuse_ks_cvm_effective)
})

test_that("fused reestimated sample KS--CvM is identical for Normal-sigma-Id", {
  set.seed(20260839L)
  x <- rnormal_sigma_Id(32L, c(0.2, -0.1), sigma = 1.15)
  compare_reestimated_sample_fusion(multiplier_bootstrap_normal_sigma_Id, x)
})

test_that("fused reestimated sample KS--CvM is identical for the general multivariate Normal", {
  set.seed(20260840L)
  x <- mvtnorm::rmvnorm(30L, mean = c(0.15, -0.25), sigma = matrix(c(1.1, 0.2, 0.2, 0.8), 2L))
  compare_reestimated_sample_fusion(
    multiplier_bootstrap_mvnormal,
    x,
    extra = list(unknown_param = "both")
  )
})
