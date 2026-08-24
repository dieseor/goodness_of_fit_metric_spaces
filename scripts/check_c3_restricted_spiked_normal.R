#!/usr/bin/env Rscript

# Read-only preflight and strict final validation for one restricted-spiked
# mean configuration. It never initializes or alters simulation results except
# for writing validation_summary.csv after a successful final validation.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source(file.path("scripts", "run_restricted_spiked_normal_covariance_alternatives.R"))

main <- function() {
  args <- restricted_spiked_parse_args(commandArgs(trailingOnly = TRUE))
  mode <- tolower(as.character(args$mode %||% "preflight"))
  mean_config <- as.character(args$mean_config %||% "")
  output_dir <- as.character(args$output_dir %||% "")
  M <- as.integer(args$M %||% 1000L)
  B <- as.integer(args$B %||% 5000L)
  base_seed <- as.integer(args$seed %||% restricted_spiked_default_seed(mean_config))
  derivative_mc_size <- as.integer(args$derivative_mc_size %||% 10000L)
  cvm_block_size <- as.integer(args$cvm_block_size %||% 50L)
  dimensions <- restricted_spiked_parse_csv(args$dimensions, c(2L, 5L), TRUE)
  n_values <- restricted_spiked_parse_csv(
    args$n_values, c(50L, 100L, 200L, 400L), TRUE
  )
  beta_values <- restricted_spiked_parse_csv(
    args$beta_values, c(0, 0.25, 0.5, 1)
  )
  lambda <- as.numeric(args$lambda %||% 2)
  allow_new <- identical(tolower(as.character(args$allow_new %||% "false")), "true")

  if (!mode %in% c("preflight", "final")) stop("Invalid validation mode.")
  if (!mean_config %in% names(restricted_spiked_mean_catalog())) {
    stop("Invalid restricted-spiked mean configuration.")
  }
  if (!nzchar(output_dir)) stop("`--output_dir` is required.")

  design <- make_restricted_spiked_design(
    mean_config, dimensions, n_values, beta_values, M
  )
  expected_manifest <- restricted_spiked_manifest(
    design, M, B, base_seed, lambda, derivative_mc_size, cvm_block_size
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
      cat("Restricted-spiked preflight passed for a new output directory.\n")
      cat(sprintf("mean_config: %s\n", mean_config))
      cat(sprintf("output_dir: %s\n", output_dir))
      cat(sprintf("completed: 0/%d\n", expected_rows))
      cat(sprintf("pending: %d\n", expected_rows))
      return(invisible(TRUE))
    }
    if (length(entries)) {
      stop("Output directory is nonempty but lacks manifest.csv and raw_results.csv.")
    }
    stop("Restricted-spiked output has not been initialized.")
  }
  if (xor(manifest_exists, results_exist)) {
    stop("Output is incomplete: manifest.csv and raw_results.csv must coexist.")
  }

  restricted_spiked_validate_manifest(manifest_path, expected_manifest)
  results <- utils::read.csv(results_path, stringsAsFactors = FALSE)
  required <- names(empty_restricted_spiked_results())
  if (!all(required %in% names(results))) {
    stop("Restricted-spiked results have an incompatible schema.")
  }
  keys <- restricted_spiked_design_key(results, include_replication = TRUE)
  expected_keys <- restricted_spiked_design_key(design, include_replication = TRUE)
  if (anyDuplicated(keys)) stop("Restricted-spiked results contain duplicate keys.")
  if (!all(keys %in% expected_keys)) stop("Results contain keys outside the design.")
  expected_ids <- match(
    restricted_spiked_design_key(results),
    restricted_spiked_design_key(unique(design[c(
      "mean_config", "d", "n", "beta", "design_id"
    )]))
  )
  design_table <- unique(design[c("mean_config", "d", "n", "beta", "design_id")])
  if (anyNA(expected_ids) || !identical(
      as.integer(results$design_id), as.integer(design_table$design_id[expected_ids]))) {
    stop("Results contain incompatible design_id values.")
  }
  seed_ok <-
    as.integer(results$seed_data) == restricted_spiked_seed(
      base_seed, results$design_id, results$replication, 0L
    ) &
    as.integer(results$seed_bootstrap) == restricted_spiked_seed(
      base_seed, results$design_id, results$replication, 1L
    ) &
    as.integer(results$seed_derivative) == restricted_spiked_seed(
      base_seed, results$design_id, results$replication, 2L
    )
  if (!isTRUE(all(seed_ok))) stop("Results contain incompatible seeds.")
  conforming <- restricted_spiked_conforming(results)
  if (!isTRUE(all(conforming))) {
    stop(sprintf("Results contain %d failed or nonconforming rows.", sum(!conforming)))
  }
  if (nrow(results) > expected_rows) stop("Results contain too many rows.")
  if (identical(mode, "final") && nrow(results) != expected_rows) {
    stop(sprintf("Final validation found %d/%d rows.", nrow(results), expected_rows))
  }

  cat(sprintf("Restricted-spiked %s validation passed.\n", mode))
  cat(sprintf("mean_config: %s\n", mean_config))
  cat(sprintf("completed: %d/%d\n", nrow(results), expected_rows))
  cat(sprintf("pending: %d\n", expected_rows - nrow(results)))

  if (identical(mode, "final")) {
    elapsed <- as.numeric(results$elapsed_seconds)
    if (any(!is.finite(elapsed)) || any(elapsed < 0)) {
      stop("Results contain invalid elapsed times.")
    }
    wall_seconds <- as.numeric(args$runner_wall_seconds %||% NA_real_)
    allocated_cores <- as.integer(args$allocated_cores %||% NA_integer_)
    effective_workers <- if (is.finite(wall_seconds) && wall_seconds > 0) {
      sum(elapsed) / wall_seconds
    } else {
      NA_real_
    }
    parallel_efficiency <- if (is.finite(effective_workers) &&
                               is.finite(allocated_cores) && allocated_cores > 0) {
      effective_workers / allocated_cores
    } else {
      NA_real_
    }
    validation <- data.frame(
      mean_config = mean_config, status = "ok", expected_rows = expected_rows,
      completed_rows = nrow(results), pending_rows = expected_rows - nrow(results),
      duplicate_keys = anyDuplicated(keys), nonconforming_rows = sum(!conforming),
      replication_elapsed_total_seconds = sum(elapsed),
      runner_wall_seconds = wall_seconds, allocated_cores = allocated_cores,
      effective_workers = effective_workers, parallel_efficiency = parallel_efficiency,
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
      validated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z",
                            tz = "Europe/Madrid"),
      output_dir = output_dir, stringsAsFactors = FALSE
    )
    restricted_spiked_write_atomic_csv(
      validation, file.path(output_dir, "validation_summary.csv")
    )
  }
  invisible(TRUE)
}

main()
