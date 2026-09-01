library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source("bootstrap/multiplier_bootstrap.R", local = globalenv())
source("bootstrap/vmf_fixed_kappa_model_spec.R", local = globalenv())

test_that("fixed-kappa vMF MLE fits only the resultant direction", {
  set.seed(1401)
  x <- rotasym::r_vMF(80L, mu = c(1, 0, 0), kappa = 2)
  theta <- fit_vmf_fixed_kappa_theta(x, null = list(type = "composite", kappa = 2))
  expect_equal(theta$kappa, 2)
  expect_equal(sqrt(sum(theta$mu^2)), 1, tolerance = 1e-14)
  expect_equal(theta$mu, colSums(x) / sqrt(sum(colSums(x)^2)), tolerance = 1e-12)
  score <- vmf_fixed_kappa_score_matrix(x, theta)
  expect_identical(dim(score), c(80L, 2L))
  expect_lt(max(abs(colSums(score))), 1e-10)
})

test_that("fixed-kappa vMF tangent score has the stated Fisher information", {
  theta <- list(mu = c(1, 0, 0, 0, 0, 0), kappa = 1)
  basis <- vmf_fixed_kappa_tangent_basis(theta$mu)
  expect_equal(crossprod(basis), diag(5L), tolerance = 1e-12)
  expect_equal(drop(crossprod(theta$mu, basis)), rep(0, 5L), tolerance = 1e-12)
  information <- vmf_fixed_kappa_information(theta)
  expect_identical(dim(information), c(5L, 5L))
  expect_true(all(eigen(information, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("fixed-kappa vMF fast multiplier stays on the requested fused C++ path", {
  set.seed(1402)
  x <- rotasym::r_vMF(18L, mu = c(1, 0, 0), kappa = 1)
  fit <- multiplier_bootstrap_vmf_fixed_kappa(
    data = x, kappa = 1, null = list(type = "composite"),
    statistics = c("ks", "cvm"), ks_grid = make_sample_unique_distance_ks_grid(),
    B = 7L, seed = 1403L, n_cores = 1L, bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
    control = list(
      fast_multiplier_cvm_block_size = 5L,
      fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE, vmf_profile_method = "tabulated",
      vmf_profile_n_u = 257L
    ),
    distance_type = "geodesic", distance_profile_backend = "r",
    fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
    fuse_ks_cvm = TRUE
  )
  expect_identical(fit$diagnostics$effective_bootstrap_method, "fast_multiplier")
  expect_identical(fit$diagnostics$derivative_method, "quadrature")
  expect_identical(fit$diagnostics$fast_multiplier_backend_effective, "cpp")
  expect_identical(fit$diagnostics$fast_multiplier_cpp_kernel_effective, "contiguous_double")
  expect_true(isTRUE(fit$diagnostics$fast_multiplier_fuse_ks_cvm_effective))
})
