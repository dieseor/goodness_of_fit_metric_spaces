library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

section6_env <- new.env(parent = globalenv())
sys.source("scripts/run_section6_new_scenarios.R", envir = section6_env)

test_that("Section 6 revised designs have the requested dimensions and beta grid", {
  for (family in c("normal", "lg", "vmf", "hvmf")) {
    design <- section6_env$make_section6_design(family)
    expect_setequal(unique(design$d), c(2L, 10L))
    expect_setequal(unique(design$beta), c(0, 0.5, 1))
    expect_setequal(unique(design$n), c(50L, 100L, 200L, 400L))
    expect_length(unique(design$scenario), 2L)
  }
})

test_that("Section 6 generators implement the new dimension-indexed scenarios", {
  catalog <- section6_env$section6_scenario_catalog()
  for (scenario in names(catalog)) {
    family <- catalog[[scenario]]$family
    design <- section6_env$make_section6_design(
      family = family,
      dimensions = 10L,
      n_values = 25L,
      beta_values = 0.5,
      scenarios = scenario
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
    } else if (identical(family, "vmf")) {
      expect_equal(ncol(x), 11L)
      expect_lt(max(abs(rowSums(x^2) - 1)), 1e-10)
    } else {
      expect_equal(ncol(x), 11L)
      expect_true(all(x[, 1L] > 0))
      expect_lt(max(abs(-x[, 1L]^2 + rowSums(x[, -1L, drop = FALSE]^2) + 1)), 1e-8)
    }
  }

  sigma_plus <- section6_env$section6_sigma(10L, "plus")
  sigma_minus <- section6_env$section6_sigma(10L, "minus")
  expect_equal(sigma_plus[1L, 2L], 0.75)
  expect_equal(sigma_minus[1L, 2L], -0.75)
  expect_true(all(eigen(sigma_plus, symmetric = TRUE)$values > 0))
  expect_true(all(eigen(sigma_minus, symmetric = TRUE)$values > 0))
})

test_that("Section 6 transverse prepilot scenarios are opt-in and preserve their fixed norms", {
  expect_false(
    "normal_1_mixture_transverse" %in%
      section6_env$make_section6_design("normal")$scenario
  )

  normal_design <- section6_env$make_section6_design(
    family = "normal", dimensions = c(2L, 10L), n_values = 20L,
    beta_values = c(0, 0.5, 1),
    scenarios = c("normal_1_mixture_transverse", "normal_2_t3")
  )
  lg_design <- section6_env$make_section6_design(
    family = "lg", dimensions = 10L, n_values = 20L,
    beta_values = 0.5, scenarios = "lg_1_mixture_transverse"
  )
  expect_setequal(unique(normal_design$scenario), c("normal_1_mixture_transverse", "normal_2_t3"))
  expect_identical(unique(lg_design$scenario), "lg_1_mixture_transverse")

  for (d in c(2L, 10L)) {
    direction <- section6_env$section6_transverse_direction(d)
    sigma_plus <- section6_env$section6_sigma_transverse(d, "plus")
    sigma_minus <- section6_env$section6_sigma_transverse(d, "minus")
    expect_equal(sum(direction^2), 1, tolerance = 1e-12)
    expect_equal(sqrt(sum((sigma_plus - diag(d))^2)), 0.75, tolerance = 1e-12)
    expect_equal(sqrt(sum((sigma_minus - diag(d))^2)), 0.75, tolerance = 1e-12)
    expect_true(all(eigen(sigma_plus, symmetric = TRUE)$values > 0))
    expect_true(all(eigen(sigma_minus, symmetric = TRUE)$values > 0))
  }
})

test_that("Section 6 fast bootstrap uses the fused C++ sample kernel", {
  ensure_distance_profile_cpp_loaded()
  for (family in c("normal", "lg", "vmf", "hvmf")) {
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
    expect_identical(
      result$fast_multiplier_cpp_kernel_effective,
      "contiguous_double"
    )
    expect_true(result$fast_multiplier_fuse_ks_cvm_effective)
    expect_identical(result$ks_grid, "sample_unique_distances")
  }
})

test_that("Section 6 manifest validation covers every execution setting except cores", {
  design <- section6_env$make_section6_design(
    family = "lg",
    dimensions = 2L,
    n_values = 50L,
    beta_values = 0
  )
  manifest <- section6_env$section6_make_manifest(
    design = design,
    M = 10L,
    B = 99L,
    cores = 2L,
    base_seed = 123L,
    derivative_mc_size = 200L,
    cvm_block_size = 25L
  )
  manifest_path <- tempfile(fileext = ".csv")
  on.exit(unlink(manifest_path, force = TRUE), add = TRUE)
  write.csv(manifest, manifest_path, row.names = FALSE)

  validate <- function() {
    section6_env$section6_validate_manifest_design(
      manifest_path = manifest_path,
      design = design,
      M = 10L,
      B = 99L,
      base_seed = 123L,
      derivative_mc_size = 200L,
      cvm_block_size = 25L
    )
  }
  expect_invisible(validate())

  changed_cores <- manifest
  changed_cores$cores <- 32L
  write.csv(changed_cores, manifest_path, row.names = FALSE)
  expect_invisible(validate())

  mismatches <- list(
    base_seed = 124L,
    derivative_mc_size = 201L,
    cvm_block_size = 26L,
    ks_grid = "different_grid",
    fast_multiplier_backend = "r",
    fast_multiplier_cpp_kernel = "legacy",
    fused_ks_cvm_kernel = FALSE
  )
  for (field in names(mismatches)) {
    incompatible <- manifest
    incompatible[[field]] <- mismatches[[field]]
    write.csv(incompatible, manifest_path, row.names = FALSE)
    expect_error(validate(), field, fixed = TRUE)
  }
})

test_that("Section 6 output locks are exclusive and record their owner", {
  output_dir <- tempfile("section6-lock-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  old_slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = NA_character_)
  on.exit({
    if (is.na(old_slurm_job_id)) {
      Sys.unsetenv("SLURM_JOB_ID")
    } else {
      Sys.setenv(SLURM_JOB_ID = old_slurm_job_id)
    }
  }, add = TRUE)
  Sys.setenv(SLURM_JOB_ID = "test-job-123")

  lock <- section6_env$section6_acquire_output_lock(output_dir)
  owner_path <- file.path(
    section6_env$section6_output_lock_path(output_dir),
    "owner.txt"
  )
  owner <- readLines(owner_path, warn = FALSE)
  expect_true(any(owner == "slurm_job_id: test-job-123"))
  expect_true(any(startsWith(owner, "hostname: ")))
  expect_true(any(startsWith(owner, "pid: ")))
  expect_error(
    section6_env$section6_acquire_output_lock(output_dir),
    "already locked",
    fixed = TRUE
  )
  expect_true(section6_env$section6_release_output_lock(lock))
  expect_false(dir.exists(section6_env$section6_output_lock_path(output_dir)))
})

test_that("Section 6 output locks are released after a controlled error", {
  output_dir <- tempfile("section6-lock-error-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  fail_while_locked <- function() {
    lock <- section6_env$section6_acquire_output_lock(output_dir)
    on.exit(section6_env$section6_release_output_lock(lock), add = TRUE)
    stop("controlled failure")
  }

  expect_error(fail_while_locked(), "controlled failure", fixed = TRUE)
  expect_false(dir.exists(section6_env$section6_output_lock_path(output_dir)))
})

test_that("Section 6 runner releases its lock when manifest validation fails", {
  output_dir <- tempfile("section6-runner-lock-error-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  design <- section6_env$make_section6_design(
    family = "lg",
    dimensions = 2L,
    n_values = 50L,
    beta_values = 0
  )
  manifest <- section6_env$section6_make_manifest(
    design = design,
    M = 1L,
    B = 9L,
    cores = 1L,
    base_seed = 123L,
    derivative_mc_size = 100L,
    cvm_block_size = 25L
  )
  write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE)

  expect_error(
    section6_env$run_section6_family(
      family = "lg",
      output_dir = output_dir,
      M = 1L,
      B = 9L,
      dimensions = 2L,
      n_values = 50L,
      beta_values = 0,
      cores = 1L,
      base_seed = 124L,
      derivative_mc_size = 100L,
      cvm_block_size = 25L,
      show_progress = FALSE
    ),
    "base_seed",
    fixed = TRUE
  )
  expect_false(dir.exists(section6_env$section6_output_lock_path(output_dir)))
})
