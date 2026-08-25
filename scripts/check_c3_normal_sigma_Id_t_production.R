#!/usr/bin/env Rscript

# Read-only preflight and strict final validation for one normal_sigma_Id
# versus standardized multivariate-t production invocation. The only write is
# production_validation_summary.csv after a successful final validation.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source(file.path("scripts", "run_normal_sigma_Id_t_pilot.R"))

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  mode <- tolower(as.character(args$mode %||% "preflight"))
  output_dir <- as.character(args$output_dir %||% "")
  d <- suppressWarnings(as.integer(args$d %||% NA_integer_))
  nu <- suppressWarnings(as.numeric(args$nu %||% NA_real_))
  beta <- suppressWarnings(as.numeric(args$beta %||% NA_real_))
  M <- suppressWarnings(as.integer(args$M %||% 1000L))
  B <- suppressWarnings(as.integer(args$B %||% 5000L))
  base_seed <- suppressWarnings(as.integer(args$seed %||% 20260728L))
  derivative_mc_size <- suppressWarnings(as.integer(
    args$derivative_mc_size %||% 10000L
  ))
  cvm_block_size <- suppressWarnings(as.integer(args$cvm_block_size %||% 50L))
  allow_new <- identical(
    tolower(as.character(args$allow_new %||% "false")), "true"
  )

  if (!mode %in% c("preflight", "final")) {
    stop("`--mode` must be 'preflight' or 'final'.", call. = FALSE)
  }
  if (!nzchar(output_dir)) stop("`--output_dir` is required.", call. = FALSE)
  if (length(d) != 1L || is.na(d) || !d %in% c(2L, 5L)) {
    stop("`--d` must be 2 or 5.", call. = FALSE)
  }
  expected_nu <- if (d == 2L) 3 else 6
  if (length(nu) != 1L || !is.finite(nu) || nu != expected_nu) {
    stop(sprintf("The production pairing requires d=%d with nu=%d.", d, expected_nu),
         call. = FALSE)
  }
  if (length(beta) != 1L || !is.finite(beta) ||
      !beta %in% normal_sigma_Id_t_production_betas) {
    stop("`--beta` must be one of 0, 0.25, 0.5 or 1.", call. = FALSE)
  }
  if (!identical(M, 1000L) || !identical(B, 5000L) ||
      !identical(derivative_mc_size, 10000L) ||
      !identical(cvm_block_size, 50L) || is.na(base_seed)) {
    stop("Production requires M=1000, B=5000, N_deriv=10000 and CvM block size 50.",
         call. = FALSE)
  }

  design <- normal_sigma_Id_t_design(
    dimensions = d,
    n_values = c(50L, 100L, 200L, 400L),
    beta_values = beta,
    M = M
  )
  expected_manifest <- unique(design[c("d", "n", "beta", "design_id")])
  expected_manifest <- transform(
    expected_manifest,
    M = M, B = B, nu = nu, base_seed = base_seed,
    derivative_mc_size = derivative_mc_size,
    cvm_block_size = cvm_block_size,
    null_model = "N_d(mu,sigma^2 I_d)",
    alternative = "standardized_multivariate_t",
    derivative_method = "score_mc", statistics = "ks,cvm",
    ks_grid = "sample_points_unique_distances",
    bootstrap_method = "fast_multiplier", fast_backend = "cpp",
    fast_kernel = "contiguous_double", fast_fused = TRUE
  )
  expected_rows <- nrow(design)
  manifest_path <- file.path(output_dir, "manifest.csv")
  results_path <- file.path(output_dir, "raw_results.csv")
  manifest_exists <- file.exists(manifest_path)
  results_exist <- file.exists(results_path)

  if (!manifest_exists && !results_exist) {
    entries <- if (dir.exists(output_dir)) {
      list.files(output_dir, all.files = TRUE, no.. = TRUE)
    } else {
      character()
    }
    if (identical(mode, "preflight") && isTRUE(allow_new) && !length(entries)) {
      cat("normal_sigma_Id production preflight passed for a new output directory.\n")
      cat(sprintf("d: %d\nnu: %g\n", d, nu))
      cat(sprintf("output_dir: %s\n", output_dir))
      cat(sprintf("completed: 0/%d\npending: %d\n", expected_rows, expected_rows))
      return(invisible(TRUE))
    }
    if (length(entries)) {
      stop("Output directory is nonempty but lacks manifest.csv and raw_results.csv.",
           call. = FALSE)
    }
    stop("Production output has not been initialized.", call. = FALSE)
  }
  if (xor(manifest_exists, results_exist)) {
    stop("manifest.csv and raw_results.csv must coexist.", call. = FALSE)
  }

  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (!isTRUE(all.equal(manifest, expected_manifest, check.attributes = FALSE))) {
    stop("Existing manifest is incompatible with this production invocation.",
         call. = FALSE)
  }
  results <- utils::read.csv(results_path, stringsAsFactors = FALSE)
  required <- names(empty_results())
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop(sprintf("raw_results.csv lacks columns: %s.", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  result_keys <- normal_sigma_Id_t_key(results)
  expected_keys <- normal_sigma_Id_t_key(design)
  if (anyDuplicated(result_keys)) stop("Results contain duplicate keys.", call. = FALSE)
  if (!all(result_keys %in% expected_keys)) {
    stop("Results contain keys outside the requested design.", call. = FALSE)
  }
  expected_design <- unique(design[c("d", "n", "beta", "design_id")])
  design_match <- match(
    paste(results$d, results$n, results$beta, sep = "|"),
    paste(expected_design$d, expected_design$n, expected_design$beta, sep = "|")
  )
  if (anyNA(design_match) || !identical(
      as.integer(results$design_id),
      as.integer(expected_design$design_id[design_match]))) {
    stop("Results contain incompatible design_id values.", call. = FALSE)
  }
  seed_ok <-
    as.integer(results$seed_data) == normal_sigma_Id_t_seed(
      base_seed, results$design_id, results$replication, 0L
    ) &
    as.integer(results$seed_bootstrap) == normal_sigma_Id_t_seed(
      base_seed, results$design_id, results$replication, 1L
    ) &
    as.integer(results$seed_derivative) == normal_sigma_Id_t_seed(
      base_seed, results$design_id, results$replication, 2L
    )
  if (!isTRUE(all(seed_ok))) stop("Results contain incompatible seeds.", call. = FALSE)

  is_conforming <- if (nrow(results)) {
    conforming(results, derivative_mc_size)
  } else {
    logical()
  }
  completed <- sum(is_conforming)
  retryable <- nrow(results) - completed
  if (nrow(results) > expected_rows) stop("Results contain too many rows.", call. = FALSE)
  if (identical(mode, "final") &&
      (nrow(results) != expected_rows || !isTRUE(all(is_conforming)))) {
    stop(sprintf(
      "Final validation found %d/%d conforming rows and %d retryable rows.",
      completed, expected_rows, retryable
    ), call. = FALSE)
  }

  cat(sprintf("normal_sigma_Id %s validation passed.\n", mode))
  cat(sprintf("d: %d\nnu: %g\n", d, nu))
  cat(sprintf("completed: %d/%d\n", completed, expected_rows))
  cat(sprintf("pending: %d\nretryable: %d\n", expected_rows - completed, retryable))

  if (identical(mode, "final")) {
    elapsed <- as.numeric(results$elapsed_seconds)
    if (any(!is.finite(elapsed)) || any(elapsed < 0)) {
      stop("Results contain invalid elapsed times.", call. = FALSE)
    }
    runner_wall_seconds <- suppressWarnings(as.numeric(
      args$runner_wall_seconds %||% NA_real_
    ))
    allocated_cores <- suppressWarnings(as.integer(
      args$allocated_cores %||% NA_integer_
    ))
    effective_workers <- if (is.finite(runner_wall_seconds) &&
                             runner_wall_seconds > 0) {
      sum(elapsed) / runner_wall_seconds
    } else {
      NA_real_
    }
    parallel_efficiency <- if (is.finite(effective_workers) &&
                               is.finite(allocated_cores) &&
                               allocated_cores > 0) {
      effective_workers / allocated_cores
    } else {
      NA_real_
    }
    validation <- data.frame(
      model = "normal_sigma_Id", d = d, nu = nu, status = "ok",
      expected_rows = expected_rows, completed_rows = completed,
      pending_rows = expected_rows - completed, retryable_rows = retryable,
      duplicate_keys = anyDuplicated(result_keys),
      replication_elapsed_total_seconds = sum(elapsed),
      runner_wall_seconds = runner_wall_seconds,
      allocated_cores = allocated_cores,
      effective_workers = effective_workers,
      parallel_efficiency = parallel_efficiency,
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
      validated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z",
                            tz = "Europe/Madrid"),
      output_dir = output_dir, stringsAsFactors = FALSE
    )
    write_atomic_csv(
      validation, file.path(output_dir, "production_validation_summary.csv")
    )
  }
  invisible(TRUE)
}

main()
