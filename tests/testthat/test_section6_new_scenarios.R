library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

section6_env <- new.env(parent = globalenv())
sys.source("scripts/run_section6_new_scenarios.R", envir = section6_env)

test_that("Section 6 revised designs have the requested dimensions and beta grid", {
  for (family in c("normal", "lg", "vmf")) {
    design <- section6_env$make_section6_design(family)
    expect_setequal(unique(design$d), c(2L, 10L))
    expect_setequal(unique(design$beta), c(0, 0.5, 1))
    expect_setequal(unique(design$n), c(50L, 100L, 200L))
    expect_length(unique(design$scenario), 2L)
  }
  expect_equal(nrow(section6_env$make_section6_design("hvmf")), 0L)
})

test_that("Section 6 generators implement the new dimension-indexed scenarios", {
  catalog <- section6_env$section6_scenario_catalog()
  for (scenario in names(catalog)) {
    family <- catalog[[scenario]]$family
    design <- section6_env$make_section6_design(
      family = family,
      dimensions = 10L,
      n_values = 25L,
      beta_values = 0.5
    )
    row <- design[design$scenario == scenario, , drop = FALSE][1L, , drop = FALSE]
    x <- section6_env$generate_section6_sample(row)
    expect_equal(nrow(x), 25L)
    if (identical(family, "normal")) {
      expect_equal(ncol(x), 10L)
    } else if (identical(family, "lg")) {
      expect_equal(ncol(x), 11L)
      expect_lt(max(abs(rowSums(x) - 1)), 1e-10)
      expect_true(all(x > 0))
    } else {
      expect_equal(ncol(x), 11L)
      expect_lt(max(abs(rowSums(x^2) - 1)), 1e-10)
    }
  }

  sigma_plus <- section6_env$section6_sigma(10L, "plus")
  sigma_minus <- section6_env$section6_sigma(10L, "minus")
  expect_equal(sigma_plus[1L, 2L], 0.75)
  expect_equal(sigma_minus[1L, 2L], -0.75)
  expect_true(all(eigen(sigma_plus, symmetric = TRUE)$values > 0))
  expect_true(all(eigen(sigma_minus, symmetric = TRUE)$values > 0))
})

test_that("Section 6 fast bootstrap uses the fused C++ sample kernel", {
  ensure_distance_profile_cpp_loaded()
  for (family in c("normal", "lg", "vmf")) {
    design <- section6_env$make_section6_design(
      family = family,
      dimensions = 10L,
      n_values = 50L,
      beta_values = 0
    )
    job <- cbind(design[1L, , drop = FALSE], rep = 1L)
    result <- section6_env$run_section6_job(
      job = job,
      B = 9L,
      base_seed = 20260729L,
      derivative_mc_size = 100L,
      cvm_block_size = 25L
    )
    expect_identical(result$status, "ok")
    expect_identical(result$bootstrap_method_effective, "fast_multiplier")
    expect_identical(result$fast_multiplier_backend_effective, "cpp")
    expect_true(result$fast_multiplier_fuse_ks_cvm_effective)
    expect_identical(result$ks_grid, "sample_unique_distances")
  }
})
