library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

suppressWarnings(source(file.path("scripts", "run_comets_projected_axis_diagnostics.R")))

test_that("director short axis is the normalized sample mean", {
  comets_data <- load_comets_real_data(finite_normals = "short")
  x <- as.matrix(comets_data$short$normal)
  axis_info <- compute_director_projection_axis(x, dataset_label = "short_period")
  mean_dir <- colMeans(x)
  mean_dir <- mean_dir / sqrt(sum(mean_dir^2))

  expect_equal(axis_info$axis_type, "sample_mean_direction")
  expect_equal(axis_info$axis, mean_dir, tolerance = 1e-12)
})

test_that("director long axis is the leading scatter eigenvector with the stored sign convention", {
  comets_data <- load_comets_real_data(finite_normals = "long")
  x <- as.matrix(comets_data$long$normal)
  axis_info <- compute_director_projection_axis(x, dataset_label = "long_period")
  scatter <- crossprod(x)
  lambda1 <- axis_info$scatter_eigenvalues[[1L]]

  expect_equal(axis_info$axis_type, "scatter_first_eigenvector")
  expect_equal(sum(axis_info$axis^2), 1, tolerance = 1e-12)
  expect_equal(drop(scatter %*% axis_info$axis), lambda1 * axis_info$axis, tolerance = 1e-8)
  expect_gte(axis_info$alignment_with_sample_mean, -1e-12)
})

test_that("small-circle arbitrary-axis helper matches the axial closed form", {
  theta <- list(mu = c(0, 0, 1), kappa = 3.2, nu = 0.15)
  z_grid <- seq(-0.9, 0.9, by = 0.1)

  observed <- small_circle_projected_cdf_on_axis(
    z_grid,
    theta = theta,
    gamma = theta$mu,
    method = "integral",
    quad_n = 512L
  )
  expected <- small_circle_axis_cdf(z_grid, kappa = theta$kappa, nu = theta$nu)

  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("beta-mixture arbitrary-axis helper matches the aligned-axis formula", {
  theta <- list(
    mu = c(0, 0, 1),
    weight1 = 0.35,
    alpha1 = 2.2,
    beta1 = 4.1,
    alpha2 = 7.5,
    beta2 = 1.8
  )
  z_grid <- seq(-0.8, 0.8, by = 0.1)

  observed <- beta_mixture2_projected_cdf_on_axis(
    z_grid,
    theta = theta,
    gamma = theta$mu,
    method = "integral",
    quad_n = 512L
  )
  expected <- beta_mixture2_cdf_y((z_grid + 1) / 2,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )

  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("uniform-beta arbitrary-axis helper matches the aligned-axis formula", {
  theta <- list(
    mu = c(0, 0, 1),
    weight_uniform = 0.2,
    alpha = 3.5,
    beta = 2.1
  )
  z_grid <- seq(-0.8, 0.8, by = 0.1)

  observed <- uniform_beta_mixture_projected_cdf_on_axis(
    z_grid,
    theta = theta,
    gamma = theta$mu,
    method = "integral",
    quad_n = 512L
  )
  expected <- uniform_beta_mixture_cdf_y((z_grid + 1) / 2,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )

  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("runner writes projected diagnostics for a small subset", {
  output_root <- file.path(tempdir(), "comets_projected_axis_smoke")
  summary_df <- run_comets_projected_axis_diagnostics(
    output_root = output_root,
    datasets = "short",
    models = c("small_circle", "uniform_beta_mixture"),
    grid_size = 61L,
    rotational_method = "integral",
    rotational_quad_n = 256L,
    save_plot = TRUE
  )

  expect_equal(nrow(summary_df), 2L)
  expect_true(all(summary_df$dataset == "short_period"))
  expect_true(all(file.exists(file.path(
    output_root,
    c("01_short_period_small_circle", "02_short_period_uniform_beta_mixture"),
    "projected_ecdf_overlay.png"
  ))))
  expect_true(all(file.exists(file.path(
    output_root,
    c("01_short_period_small_circle", "02_short_period_uniform_beta_mixture"),
    "projected_cdf_grid.csv"
  ))))
  expect_true(file.exists(file.path(output_root, "projected_axis_diagnostics_summary.csv")))
})
