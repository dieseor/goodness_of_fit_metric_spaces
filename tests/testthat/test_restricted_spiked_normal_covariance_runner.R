library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

runner_env <- new.env(parent = globalenv())
sys.source(
  "scripts/run_restricted_spiked_normal_covariance_alternatives.R",
  envir = runner_env
)

test_that("paper mean configurations have exact norms and orthogonal q", {
  expected_norms <- c(
    axis_075 = 0.75, axis_100 = 1, diagonal_150 = 1.5, diagonal_100 = 1
  )
  for (mean_config in names(expected_norms)) {
    for (d in c(2L, 5L)) {
      theta <- runner_env$restricted_spiked_mean_vector(d, mean_config)
      q <- runner_env$restricted_spiked_second_direction(d, mean_config)
      expect_equal(sqrt(sum(theta^2)), expected_norms[[mean_config]],
                   tolerance = 1e-14)
      expect_equal(sqrt(sum(q^2)), 1, tolerance = 1e-14)
      expect_equal(sum(theta * q), 0, tolerance = 1e-14)
    }
  }
  expect_equal(
    runner_env$restricted_spiked_mean_vector(5L, "axis_075"),
    c(0.75, 0, 0, 0, 0)
  )
  expect_equal(
    runner_env$restricted_spiked_mean_vector(5L, "diagonal_150"),
    rep(1.5 / sqrt(5), 5L)
  )
  expect_identical(
    vapply(
      names(expected_norms), runner_env$restricted_spiked_default_seed,
      integer(1)
    ),
    c(
      axis_075 = 20260831L, axis_100 = 20260833L,
      diagonal_150 = 20260831L, diagonal_100 = 20260833L
    )
  )
})

test_that("paper design has 32,000 jobs per mean configuration", {
  runner_defaults <- formals(
    runner_env$run_restricted_spiked_covariance_alternatives
  )
  expect_identical(runner_defaults$M, 1000L)
  expect_identical(runner_defaults$B, 5000L)
  expect_identical(runner_defaults$derivative_mc_size, 10000L)

  all_ids <- integer()
  for (mean_config in names(runner_env$restricted_spiked_mean_catalog())) {
    design <- runner_env$make_restricted_spiked_design(mean_config)
    expect_equal(nrow(design), 32000L)
    expect_setequal(unique(design$d), c(2L, 5L))
    expect_setequal(unique(design$n), c(50L, 100L, 200L, 400L))
    expect_setequal(unique(design$beta), c(0, 0.25, 0.5, 1))
    expect_equal(length(unique(design$design_id)), 32L)
    expect_equal(as.integer(table(design$design_id)), rep(1000L, 32L))
    expect_false(any(design$design_id %in% all_ids))
    all_ids <- c(all_ids, unique(design$design_id))
  }
})

test_that("paper manifest records and enforces every statistical setting", {
  design <- runner_env$make_restricted_spiked_design("axis_075", M = 2L)
  manifest <- runner_env$restricted_spiked_manifest(
    design, M = 2L, B = 19L, base_seed = 123L, lambda = 2,
    derivative_mc_size = 100L, cvm_block_size = 25L
  )
  expect_true(all(manifest$derivative_method == "score_mc"))
  expect_true(all(manifest$vhat_method == "restricted_spiked_analytic_fisher"))
  expect_true(all(manifest$ks_grid == "sample_points_unique_distances"))
  expect_true(all(manifest$cvm_grid == "sample_points_unique_distances"))
  expect_true(all(manifest$ks_cvm_sample_grid_shared))
  expect_true(all(manifest$fast_backend == "cpp"))
  expect_true(all(manifest$fast_kernel == "contiguous_double"))
  expect_true(all(manifest$fast_fused))
  expect_true(all(!manifest$fallback_allowed))

  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  write.csv(manifest, path, row.names = FALSE)
  expect_invisible(runner_env$restricted_spiked_validate_manifest(path, manifest))
  incompatible <- manifest
  incompatible$B <- 20L
  write.csv(incompatible, path, row.names = FALSE)
  expect_error(
    runner_env$restricted_spiked_validate_manifest(path, manifest),
    "B",
    fixed = TRUE
  )
})

test_that("restricted-spiked output lock is exclusive and records Slurm owner", {
  output_dir <- tempfile("restricted-spiked-lock-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  old_job <- Sys.getenv("SLURM_JOB_ID", unset = NA_character_)
  on.exit({
    if (is.na(old_job)) Sys.unsetenv("SLURM_JOB_ID") else Sys.setenv(SLURM_JOB_ID = old_job)
  }, add = TRUE)
  Sys.setenv(SLURM_JOB_ID = "test-restricted-job")
  lock <- runner_env$restricted_spiked_acquire_lock(output_dir)
  owner <- readLines(file.path(lock$path, "owner.txt"), warn = FALSE)
  expect_true("slurm_job_id: test-restricted-job" %in% owner)
  expect_error(runner_env$restricted_spiked_acquire_lock(output_dir), "locked")
  expect_true(runner_env$restricted_spiked_release_lock(lock))
})

test_that("small restricted-spiked run is conforming and resumable", {
  output_dir <- tempfile("restricted-spiked-run-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  settings <- list(
    mean_config = "axis_075", output_dir = output_dir, M = 1L, B = 9L,
    dimensions = 2L, n_values = 50L, beta_values = 0,
    derivative_mc_size = 100L, cvm_block_size = 25L, cores = 1L,
    checkpoint_results = 1L, base_seed = 321L, show_progress = FALSE
  )
  do.call(runner_env$run_restricted_spiked_covariance_alternatives, settings)
  first <- read.csv(file.path(output_dir, "raw_results.csv"),
                    stringsAsFactors = FALSE)
  expect_equal(nrow(first), 1L)
  expect_true(runner_env$restricted_spiked_conforming(first))
  expect_identical(first$sample_grid_mode, "sample_points_unique_distances")
  expect_true(first$shared_sample_ks_cvm_cache)
  expect_equal(first$theta_norm_true, 0.75, tolerance = 1e-14)
  expect_equal(first$beta, 0)

  do.call(runner_env$run_restricted_spiked_covariance_alternatives, settings)
  resumed <- read.csv(file.path(output_dir, "raw_results.csv"),
                      stringsAsFactors = FALSE)
  expect_equal(nrow(resumed), 1L)
  expect_equal(resumed$seed_data, first$seed_data)
  expect_false(dir.exists(runner_env$restricted_spiked_lock_path(output_dir)))
})
