normal_sigma_Id_bootstrap_path <- c(
  file.path("bootstrap", "normal_sigma_Id_bootstrap.R"),
  file.path("..", "..", "bootstrap", "normal_sigma_Id_bootstrap.R")
)
normal_sigma_Id_bootstrap_path <- normal_sigma_Id_bootstrap_path[
  file.exists(normal_sigma_Id_bootstrap_path)
][1L]
source(normal_sigma_Id_bootstrap_path)

test_that("normal-sigma-Id MLE has the closed-form value", {
  set.seed(912)
  x <- matrix(rnorm(120), ncol = 3L)
  fit <- fit_normal_sigma_Id_theta(x, null = list(type = "composite"))
  expect_equal(fit$mu, colMeans(x), tolerance = 1e-14)
  expect_equal(fit$sigma^2, sum(rowSums(sweep(x, 2L, colMeans(x))^2)) / (nrow(x) * ncol(x)), tolerance = 1e-14)
})

test_that("normal-sigma-Id score agrees with numerical likelihood derivatives", {
  x <- matrix(c(0.3, -1.2, 0.7), nrow = 1L)
  theta <- list(mu = c(0.2, -0.4, 0.1), sigma = 1.3)
  analytic <- normal_sigma_Id_score_matrix(x, theta)
  par <- c(theta$mu, theta$sigma)
  objective <- function(p) normal_sigma_Id_loglik(x, list(mu = p[1:3], sigma = p[[4L]]))
  h <- 1e-6
  numeric <- vapply(seq_along(par), function(j) {
    plus <- minus <- par; plus[[j]] <- plus[[j]] + h; minus[[j]] <- minus[[j]] - h
    (objective(plus) - objective(minus)) / (2 * h)
  }, numeric(1L))
  expect_equal(as.numeric(analytic), numeric, tolerance = 1e-6)
})

test_that("normal-sigma-Id Fisher equals score covariance and reuses the Gaussian profile", {
  set.seed(913)
  theta <- list(mu = c(0.2, -0.3), sigma = 1.4)
  x <- rnormal_sigma_Id(30000L, theta$mu, theta$sigma)
  expected <- normal_sigma_Id_fisher_information(theta)
  expect_equal(unname(cov(normal_sigma_Id_score_matrix(x, theta))), unname(expected), tolerance = 0.05)
  spherical <- make_normal_sigma_Id_spec()
  general <- make_mvnormal_spec("both")
  t_values <- c(0.2, 0.8, 1.7)
  expect_equal(
    spherical$profile_eval(c(-0.4, 0.5), t_values, theta),
    general$profile_eval(c(-0.4, 0.5), t_values, list(mu = theta$mu, Sigma = diag(theta$sigma^2, 2L))),
    tolerance = 1e-12
  )
})

test_that("normal-sigma-Id fast multiplier uses score MC and the analytic Fisher", {
  set.seed(914)
  x <- rnormal_sigma_Id(25L, c(0.1, -0.2), 1.1)
  fit <- multiplier_bootstrap_normal_sigma_Id(
    x, null = list(type = "composite"), statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(), B = 19L,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE),
    control = list(derivative_mc_size = 1000L, derivative_mc_seed = 22L,
      fast_multiplier_cvm_block_size = 20L),
    fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
    fuse_ks_cvm = TRUE
  )
  expect_identical(fit$diagnostics$derivative_method_effective, "score_mc")
  expect_identical(fit$diagnostics$vhat_method, "normal_sigma_Id_analytic_fisher")
  expect_identical(fit$diagnostics$fast_multiplier_backend_effective, "cpp")
  expect_identical(fit$diagnostics$fast_multiplier_cpp_kernel_effective, "contiguous_double")
  expect_true(fit$diagnostics$fast_multiplier_fuse_ks_cvm_effective)
})
