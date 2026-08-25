#!/usr/bin/env Rscript

# Read-only preflight and strict final validation for the C3 Section 6 normal
# and logistic-Gaussian production runs. Generated simulation results remain
# unversioned; this script only verifies that the required local C3 artifacts
# exist and are compatible before the production runner is allowed to start.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source(file.path("scripts", "run_section6_new_scenarios.R"))

main <- function() {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  family <- tolower(as.character(args$family %||% ""))
  output_dir <- as.character(args$output_dir %||% "")
  base_seed <- suppressWarnings(as.integer(args$seed %||% NA_integer_))
  derivative_method <- tolower(as.character(
    args$derivative_method %||% "score_mc"
  ))
  mode <- tolower(as.character(args$mode %||% "preflight"))
  dimensions <- parse_section6_csv(
    args$dimensions,
    c(2L, 10L),
    "integer"
  )
  allow_new <- identical(
    tolower(as.character(args$allow_new %||% "false")),
    "true"
  )

  if (!family %in% c("normal", "lg")) {
    stop("`--family` must be either 'normal' or 'lg'.", call. = FALSE)
  }
  if (!nzchar(output_dir)) {
    stop("`--output_dir` is required.", call. = FALSE)
  }
  if (is.na(base_seed)) {
    stop("`--seed` must be a finite integer.", call. = FALSE)
  }
  if (!mode %in% c("preflight", "final")) {
    stop("`--mode` must be either 'preflight' or 'final'.", call. = FALSE)
  }
  if (!identical(derivative_method, "score_mc")) {
    stop("Section 6 Normal/LG production requires `--derivative_method=score_mc`.",
         call. = FALSE)
  }
  if (!length(dimensions) || anyNA(dimensions) ||
      any(!is.finite(dimensions)) || any(dimensions < 2L)) {
    stop("`--dimensions` must contain integers greater than or equal to 2.",
         call. = FALSE)
  }
  dimensions <- sort(unique(as.integer(dimensions)))

  production_design <- make_section6_design(
    family = family,
    dimensions = dimensions,
    n_values = c(50L, 100L, 200L, 400L),
    beta_values = c(0, 0.5, 1)
  )
  expected_rows <- nrow(production_design) * 1000L

  manifest_path <- file.path(output_dir, "manifest.csv")
  result_path <- file.path(output_dir, "raw_results.csv")
  manifest_exists <- file.exists(manifest_path)
  results_exist <- file.exists(result_path)

  if (!manifest_exists && !results_exist) {
    if (identical(mode, "preflight") && isTRUE(allow_new)) {
      existing_entries <- if (dir.exists(output_dir)) {
        list.files(
          output_dir,
          all.files = TRUE,
          no.. = TRUE
        )
      } else {
        character()
      }
      if (length(existing_entries)) {
        stop(sprintf(
          paste(
            "The requested new production output directory is not empty:",
            "'%s'. Refusing to initialize it."
          ),
          output_dir
        ), call. = FALSE)
      }
      cat("Production preflight passed for a new output directory.\n")
      cat(sprintf("family: %s\n", family))
      cat(sprintf("dimensions: %s\n", paste(dimensions, collapse = ",")))
      cat(sprintf("output_dir: %s\n", output_dir))
      cat(sprintf("completed: 0/%d\n", expected_rows))
      cat(sprintf("pending: %d\n", expected_rows))
      return(invisible(TRUE))
    }
    if (identical(mode, "preflight") && identical(family, "normal")) {
      stop(paste(
        "The Normal production output has not been initialized.",
        "The 960 conforming Phase B rows are production-compatible and should",
        "be recovered explicitly before submission instead of recomputing them.",
        "Use scripts/run_section6_new_scenarios.R with --recover_from set to",
        "the verified Normal Phase B job directory, --output_dir set to:",
        sprintf("  %s", output_dir),
        "and the final M=1000, B=5000, design and seed settings.",
        "Then re-run this preflight.",
        sep = "\n"
      ), call. = FALSE)
    }
    if (identical(mode, "preflight") && identical(family, "lg")) {
      stop(paste(
        "The logistic-Gaussian production resume files are absent on this",
        "machine. They are intentionally ignored by Git and must be copied",
        "manually to the C3 before submitting the LG production job.",
        "Required files:",
        sprintf("  %s", manifest_path),
        sprintf("  %s", result_path),
        "Files such as summary.csv, progress_status.txt and run.log are useful",
        "but are not required for resume. Re-run the preflight after copying.",
        sep = "\n"
      ), call. = FALSE)
    }
    stop(sprintf(
      "Final production output is absent in '%s'.", output_dir
    ), call. = FALSE)
  }
  if (xor(manifest_exists, results_exist)) {
    missing_path <- if (manifest_exists) result_path else manifest_path
    stop(sprintf(
      paste(
        "Production output is incomplete: '%s' is missing.",
        "Do not submit or resume until manifest.csv and raw_results.csv",
        "come from the same output directory."
      ),
      missing_path
    ), call. = FALSE)
  }

  output_lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(output_lock), add = TRUE)

  section6_validate_manifest_design(
    manifest_path = manifest_path,
    design = production_design,
    M = 1000L,
    B = 5000L,
    base_seed = base_seed,
    derivative_mc_size = 10000L,
    cvm_block_size = 50L,
    derivative_method = derivative_method
  )

  results <- utils::read.csv(result_path, stringsAsFactors = FALSE)
  required <- names(empty_section6_results())
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop(sprintf(
      "Production results lack columns: %s.",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  if (nrow(results) > expected_rows) {
    stop(sprintf(
      "Production results contain more than %d rows.",
      expected_rows
    ), call. = FALSE)
  }

  production_jobs <- merge(
    production_design,
    data.frame(rep = seq_len(1000L)),
    by = NULL
  )
  production_keys <- section6_design_key(production_jobs, include_rep = TRUE)
  result_keys <- section6_design_key(results, include_rep = TRUE)
  if (anyDuplicated(result_keys)) {
    stop("Production results contain duplicate replication keys.",
         call. = FALSE)
  }
  if (!all(result_keys %in% production_keys)) {
    stop("Production results contain keys outside the requested design.",
         call. = FALSE)
  }

  expected_design_id <- match(
    section6_design_key(results),
    section6_design_key(production_design)
  )
  if (anyNA(expected_design_id) ||
      !identical(as.integer(results$design_id), as.integer(expected_design_id))) {
    stop("Production results contain incompatible design_id values.",
         call. = FALSE)
  }

  conforming <- results$status == "ok" &
    results$family == family &
    results$bootstrap_method_requested == "fast_multiplier" &
    results$bootstrap_method_effective == "fast_multiplier" &
    results$fast_multiplier_backend_requested == "cpp" &
    results$fast_multiplier_backend_effective == "cpp" &
    results$fast_multiplier_cpp_kernel_requested == "contiguous_double" &
    results$fast_multiplier_cpp_kernel_effective == "contiguous_double" &
    results$fast_multiplier_fuse_ks_cvm_requested &
    results$fast_multiplier_fuse_ks_cvm_effective &
    results$derivative_method_requested == "score_mc" &
    results$derivative_method_effective == "score_mc" &
    results$derivative_method_selection_source == "explicit" &
    results$ks_grid == "sample_unique_distances" &
    !results$fallback_to_reestimated
  if (!isTRUE(all(conforming))) {
    failed <- which(is.na(conforming) | !conforming)
    stop(sprintf(
      "Production results contain %d failed or nonconforming rows.",
      length(failed)
    ), call. = FALSE)
  }

  seed_ok <-
    as.integer(results$seed_data) ==
      section6_seed(base_seed, results$design_id, results$rep, 0L) &
    as.integer(results$seed_bootstrap) ==
      section6_seed(base_seed, results$design_id, results$rep, 1L) &
    as.integer(results$seed_derivative) ==
      section6_seed(base_seed, results$design_id, results$rep, 2L)
  if (!isTRUE(all(seed_ok))) {
    stop("Production results contain incompatible seeds.", call. = FALSE)
  }

  completed <- nrow(results)
  if (identical(mode, "final") && completed != expected_rows) {
    stop(sprintf(
      "Final production validation found %d/%d rows.",
      completed,
      expected_rows
    ), call. = FALSE)
  }

  cat(sprintf("Production %s validation passed.\n", mode))
  cat(sprintf("family: %s\n", family))
  cat(sprintf("dimensions: %s\n", paste(dimensions, collapse = ",")))
  cat(sprintf("output_dir: %s\n", output_dir))
  cat(sprintf("completed: %d/%d\n", completed, expected_rows))
  cat(sprintf("pending: %d\n", expected_rows - completed))

  if (identical(mode, "final")) {
    elapsed <- as.numeric(results$elapsed_seconds)
    if (any(!is.finite(elapsed)) || any(elapsed < 0)) {
      stop("Final production results contain invalid elapsed times.",
           call. = FALSE)
    }
    summary <- data.frame(
      family = family,
      dimensions = paste(dimensions, collapse = ","),
      status = "ok",
      completed_rows = completed,
      pending_rows = expected_rows - completed,
      duplicate_keys = anyDuplicated(result_keys),
      nonconforming_rows = sum(!conforming),
      replication_elapsed_total_seconds = sum(elapsed),
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
      validated_at = format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S %Z",
        tz = "Europe/Madrid"
      ),
      output_dir = output_dir,
      stringsAsFactors = FALSE
    )
    summary_path <- file.path(
      output_dir,
      "production_validation_summary.csv"
    )
    section6_write_atomic_csv(summary, summary_path)
    cat(sprintf("summary: %s\n", summary_path))
  }

  invisible(TRUE)
}

main()
