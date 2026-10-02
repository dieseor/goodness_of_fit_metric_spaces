library(gofmetric)

test_that("the public model list is explicit", {
  models <- c("normal", "mvnormal", "restricted_spiked_normal", "normal_sigma_Id",
              "vmf", "vmf_fixed_kappa", "hvmf", "logistic_gaussian",
              "logistic_gaussian_ar1", "cardioid", "small_circle",
              "uniform_beta_mixture", "sunspots_joint_time_space")
  expect_true(all(paste0("multiplier_bootstrap_", models) %in%
                    getNamespaceExports("gofmetric")))
  expect_false(any(grepl("watson|_jp$|spherical_cauchy", getNamespaceExports("gofmetric"))))
})

test_that("the software citation follows DESCRIPTION", {
  entries <- citation("gofmetric")
  expect_length(entries, 2L)
  expect_match(entries[[1L]]$note, as.character(packageVersion("gofmetric")),
               fixed = TRUE)
  expect_length(entries[[2L]]$author, 3L)
})

test_that("normal and one-dimensional mvnormal agree on a simple null", {
  x <- c(-1.3, -0.8, -0.1, 0.2, 0.5, 0.9, 1.4, 1.7)
  normal <- gof_test(x, "normal", theta = list(mu = 0, sigma = 1),
                     bootstrap = "reestimated", statistics = "ks", B = 9, seed = 12)
  mvnormal <- gof_test(matrix(x, ncol = 1), "mvnormal",
                       theta = list(mu = 0, Sigma = matrix(1, 1, 1)),
                       bootstrap = "reestimated", statistics = "ks", B = 9, seed = 12)
  expect_equal(normal$observed$ks$statistic, mvnormal$observed$ks$statistic,
               tolerance = 1e-7)
  expect_true(is.finite(normal$inference$ks$p_value))
  expect_true(is.finite(mvnormal$inference$ks$p_value))
  expect_equal(normal$inference$ks$p_value, mvnormal$inference$ks$p_value)
})

test_that("the joint shared model rejects unavailable re-estimation", {
  x <- cbind(matrix(rep(c(0, 0, 1), 3), ncol = 3, byrow = TRUE),
             c(0.2, 0.5, 0.8))
  expect_error(gof_test(x, "sunspots_joint_time_space", bootstrap = "reestimated"),
               "only.*fast_multiplier")
  expect_error(multiplier_bootstrap_sunspots_joint_time_space(
    x, list(type = "composite"), control = list(hemisphere_regression = "asymmetric")),
    "requires.*shared")
  expect_error(gof_test(x, "sunspots_joint_time_space", theta = list()),
               "requires.*composite")
})

test_that("the installed joint model runs for list and matrix data", {
  set.seed(1)
  x <- matrix(rnorm(36), ncol = 3)
  x <- x / sqrt(rowSums(x * x))
  data <- list(x = x, s = seq(0.1, 0.9, length.out = 12))
  control <- list(profile_l_max = 4L, profile_quad_n = 16L,
                  time_quad_n = 3L, derivative_mc_size = 10L,
                  time_beta_n_starts = 1L,
                  time_beta_nelder_mead_control = list(maxit = 10L),
                  optim_control = list(maxit = 10L))
  a <- gof_test(data, "sunspots_joint_time_space", statistics = "ks",
                B = 1, seed = 7, control = control)
  b <- gof_test(cbind(data$x, data$s), "sunspots_joint_time_space",
                statistics = "ks", B = 1, seed = 7, control = control)
  expect_true(is.finite(a$inference$ks$p_value))
  expect_equal(a$inference$ks$p_value, b$inference$ks$p_value)
  expect_equal(a$bootstrap$used, "fast")
})
