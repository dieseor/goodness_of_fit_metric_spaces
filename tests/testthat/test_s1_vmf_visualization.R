library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("utils.R"))
source(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))
source(file.path("convergence_empirical_process", "gaussian_process_s1_visualization.R"))

test_that("generate_circle_grid returns deterministic unit vectors in angular order", {
  grid <- generate_circle_grid(6)
  expect_true(is.data.frame(grid))
  expect_equal(nrow(grid), 6)
  expect_equal(colnames(grid), c("theta", "x", "y", "label"))
  expect_true(all(diff(grid$theta) > 0))
  expect_equal(grid$theta[1], 0, tolerance = 1e-12)
  norms <- sqrt(grid$x^2 + grid$y^2)
  expect_equal(norms, rep(1, 6), tolerance = 1e-10)
})

test_that("S1 chordal distance profile is monotone with correct endpoints", {
  mu <- c(1, 0)
  omega <- c(1, 0)
  t_values <- seq(0, 2, length.out = 9)
  profile <- theoretical_distance_profile_vmf_s1_chordal(
    omega = omega,
    mu = mu,
    kappa = 2,
    t_values = t_values,
    cdf_grid_size = 4097
  )
  expect_equal(profile[1], 0, tolerance = 1e-12)
  expect_equal(profile[length(profile)], 1, tolerance = 1e-12)
  expect_true(all(diff(profile) >= -1e-8))
})

test_that("S1 chordal inversion is consistent with the deterministic CDF", {
  mu <- c(1, 0)
  omega <- c(cos(pi / 4), sin(pi / 4))
  cdf_object <- build_vmf_s1_cdf(mu, kappa = 2, n_grid = 4097)
  u_values <- seq(0.1, 0.9, by = 0.2)
  t_values <- invert_distance_profile_vmf_s1_chordal(
    omega = omega,
    mu = mu,
    kappa = 2,
    u_values = u_values,
    cdf_object = cdf_object,
    tol = 1e-8
  )
  recovered <- theoretical_distance_profile_vmf_s1_chordal(
    omega = omega,
    mu = mu,
    kappa = 2,
    t_values = t_values,
    cdf_object = cdf_object
  )
  expect_equal(recovered, u_values, tolerance = 1e-4)
})

test_that("angular labels are formatted as fractions of pi", {
  labels <- format_radians_label(c(0, pi / 2, pi, 3 * pi / 2))
  expect_identical(labels, c("0", "π/2", "π", "3π/2"))
})

test_that("S1 Gaussian realization is built from a standard normal draw via Cholesky", {
  sigma <- matrix(c(2, 0.3, 0.3, 1), nrow = 2)
  res <- draw_single_gaussian_realization(sigma, seed = 123)
  expect_equal(length(res$standard_normal_draw), 2)
  expect_equal(
    sigma,
    t(res$chol_upper) %*% res$chol_upper,
    tolerance = 1e-10
  )
  expect_equal(
    res$realization_vec,
    as.numeric(t(res$chol_upper) %*% res$standard_normal_draw),
    tolerance = 1e-12
  )
})

test_that("Exact S1 joint probabilities are symmetric", {
  mu <- c(1, 0)
  omega1 <- c(1, 0)
  omega2 <- c(cos(2 * pi / 3), sin(2 * pi / 3))
  cdf_object <- build_vmf_s1_cdf(mu, kappa = 2, n_grid = 4097)
  p12 <- joint_probability_vmf_s1_chordal_exact(
    omega1 = omega1,
    t1 = 0.75,
    omega2 = omega2,
    t2 = 1.10,
    mu = mu,
    kappa = 2,
    cdf_object = cdf_object
  )
  p21 <- joint_probability_vmf_s1_chordal_exact(
    omega1 = omega2,
    t1 = 1.10,
    omega2 = omega1,
    t2 = 0.75,
    mu = mu,
    kappa = 2,
    cdf_object = cdf_object
  )
  expect_equal(p12, p21, tolerance = 1e-10)
})

test_that("Exact S1 covariance route returns a finite valid matrix on a small grid", {
  mu <- c(1, 0)
  omega_grid <- as.matrix(generate_circle_grid(4)[, c("x", "y")])
  t_grid <- c(0.3, 0.9, 1.5)
  sigma_exact <- cov_vmf(
    omega_grid = omega_grid,
    t_grid = t_grid,
    mu = mu,
    kappa = 2,
    distance_type = "chordal",
    h0 = "simple",
    cov_method = "exact_s1_simple",
    n_cores = 1,
    cdf_grid_size = 4097
  )
  expect_true(all(is.finite(sigma_exact)))
  diag_exact <- validate_covariance_matrix(sigma_exact, stop_on_failure = FALSE)
  expect_true(isTRUE(diag_exact$valid))
})

test_that("S1 visualization works in the exact covariance mode", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("mvtnorm")

  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 4,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 30,
    cdf_grid_size = 4097
  )
  expect_true(!is.null(res$plot))
  expect_true(is.data.frame(res$curve_data))
  expect_true(nrow(res$curve_data) > 0)
})

test_that("S1 visualization can be displayed directly as a function of t", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mvtnorm")

  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 4,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 25,
    plot_domain = "t",
    cdf_grid_size = 4097,
    seed = 123
  )
  expect_true(!is.null(res$plot))
  expect_identical(res$plot_domain, "t")
  expect_equal(range(res$curve_data$plot_x), c(0, 2), tolerance = 1e-12)
})

test_that("u and t displays use the same Gaussian realization when the seed is fixed", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mvtnorm")

  res_u <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 4,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 25,
    plot_domain = "u",
    cdf_grid_size = 4097,
    seed = 123
  )
  res_t <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 4,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 25,
    plot_domain = "t",
    cdf_grid_size = 4097,
    seed = 123
  )
  expect_equal(res_u$realization_vec, res_t$realization_vec, tolerance = 0)
})

test_that("S1 visualization colors are symmetric around the mean direction", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mvtnorm")

  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 6,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 20,
    cdf_grid_size = 4097
  )

  color_lookup <- unique(res$curve_data[, c("theta", "color")])
  color_lookup <- color_lookup[order(color_lookup$theta), , drop = FALSE]
  expect_identical(color_lookup$color[2], color_lookup$color[6])
  expect_identical(color_lookup$color[3], color_lookup$color[5])
})

test_that("legacy symmetric color alias maps to yellow_blue without changing colors", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mvtnorm")

  res_alias <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 6,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 20,
    color_scheme = "symmetric",
    cdf_grid_size = 4097,
    seed = 123
  )
  res_named <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 6,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 20,
    color_scheme = "yellow_blue",
    cdf_grid_size = 4097,
    seed = 123
  )
  expect_identical(res_alias$color_scheme, "yellow_blue")
  expect_identical(res_alias$curve_data$color, res_named$curve_data$color)
})

test_that("S1 visualization supports rainbow colors", {
  skip_if_not_installed("ggplot2")

  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 6,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 20,
    color_scheme = "rainbow",
    cdf_grid_size = 4097,
    seed = 123
  )
  expect_identical(res$color_scheme, "rainbow")
  expect_equal(length(unique(res$curve_data$color)), 6)
})

test_that("S1 visualization can be saved and replotted from stored results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mvtnorm")

  out_rds <- file.path(tempdir(), "limit_gaussian_s1_result.rds")
  out_png <- file.path(tempdir(), "limit_gaussian_s1_replot.png")
  out_base <- sub("\\.png$", "", out_png)
  unlink(c(out_rds, out_png, paste0(out_base, "_left.png"), paste0(out_base, "_right.png")))

  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 6,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 20,
    color_scheme = "rainbow",
    plot_domain = "t",
    cdf_grid_size = 4097,
    seed = 123,
    save_result = out_rds
  )

  replotted <- replot_limit_gaussian_s1_vmf_from_result(
    out_rds,
    color_scheme = "yellow_blue",
    save_plot = out_png
  )

  expect_true(file.exists(out_rds))
  expect_true(file.exists(out_png))
  expect_true(file.exists(paste0(out_base, "_left.png")))
  expect_true(file.exists(paste0(out_base, "_right.png")))
  expect_identical(replotted$plot_domain, "t")
  expect_identical(replotted$realization_vec, res$realization_vec)
  expect_false(identical(replotted$curve_data$color, res$curve_data$color))
})

test_that("S1 visualization save_plot writes combined and separate files", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mvtnorm")

  out_path <- file.path(tempdir(), "limit_gaussian_s1_test.png")
  out_base <- sub("\\.png$", "", out_path)
  unlink(c(out_path, paste0(out_base, "_left.png"), paste0(out_base, "_right.png")))

  visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    n_angles = 4,
    t_grid = c(0.3, 0.9, 1.5),
    cov_method = "exact_s1_simple",
    n_cores = 1,
    curve_points = 20,
    cdf_grid_size = 4097,
    save_plot = out_path
  )

  expect_true(file.exists(out_path))
  expect_true(file.exists(paste0(out_base, "_left.png")))
  expect_true(file.exists(paste0(out_base, "_right.png")))
})

test_that("S1 visualization works in the Monte Carlo covariance mode on a stable small grid", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("mvtnorm")

  omega_grid <- as.matrix(generate_circle_grid(2)[, c("x", "y")])
  res <- visualize_limit_gaussian_s1_vmf(
    mu = c(1, 0),
    kappa = 2,
    omega_grid = omega_grid,
    t_grid = c(0.5, 1.2),
    cov_method = "mc",
    n_mc_samples = 20000,
    n_cores = 1,
    seed = 123,
    curve_points = 25,
    cdf_grid_size = 4097
  )
  expect_true(!is.null(res$plot))
  expect_true(is.data.frame(res$curve_data))
  expect_true(nrow(res$curve_data) > 0)
})
