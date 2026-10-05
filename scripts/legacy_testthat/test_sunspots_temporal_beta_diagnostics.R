library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "sunspots", "run_sunspots_cycle23_joint_time_space_gof.R"))

test_that("joint GOF seed defaults are distinct and remain in output metadata", {
  defaults <- formals(run_sunspots_cycle23_joint_time_space_gof)
  seed_names <- c("center_seed", "dequantization_seed", "derivative_mc_seed", "bootstrap_seed")
  seed_values <- vapply(seed_names, function(name) eval(defaults[[name]]), integer(1L))
  expect_equal(length(unique(seed_values)), length(seed_values))
  parsed <- parse_sunspots_joint_time_space_args(c(
    "--center_seed=1", "--dequantization_seed=2", "--derivative_mc_seed=3", "--bootstrap_seed=4"
  ))
  expect_equal(unname(unlist(parsed[seed_names])), 1:4)

  eta <- list(
    weight1 = 0.4, alpha1 = 3, beta1 = 7, alpha2 = 8, beta2 = 3,
    mean1 = 0.3, mean2 = 0.7, opt = list(convergence = 0L),
    boundary_flags = list(weight = FALSE, shape_lower = FALSE, shape_upper = FALSE)
  )
  fit <- list(
    theta_hat = list(a_N = .5, b_N = -.1, a_S = .4, b_S = -.2, c = 12, opt = list(convergence = 0L)),
    eta_hat = eta, n_parameters = 10L, temporal_loglik = -12,
    conditional_loglik = -15, loglik = -27, aic = 74, bic = 88
  )
  settings <- list(
    hemisphere_regression = "asymmetric", B = 2L, n_cores = 1L,
    observed_profile_n_cores = 1L,
    start_date = "1997-06-01", end_date = "2006-01-01",
    center_seed = 1L, dequantization_seed = 2L, derivative_mc_seed = 3L, bootstrap_seed = 4L,
    time_quad_n = 64L, profile_l_max = 100L, spatial_quad_n = 400L,
    distance_profile_backend = "r", allow_boundary_fast = FALSE
  )
  timing <- as.list(stats::setNames(rep(0, 6L), c(
    "temporal_mle_seconds", "spatial_mle_seconds", "observed_profile_seconds",
    "fast_preparation_seconds", "bootstrap_seconds", "total_seconds"
  )))
  summary <- sunspots_joint_summary_row(
    "ks", c(.1, .2), .15, fit,
    list(derivative_mc_size = 10L, diagnostics = list(Vhat_rcond = .1)),
    list(x = matrix(0, nrow = 2L, ncol = 3L)), 1L, settings, timing
  )
  expect_equal(unname(unlist(summary[seed_names])), 1:4)
})