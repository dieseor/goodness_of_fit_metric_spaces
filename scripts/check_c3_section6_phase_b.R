#!/usr/bin/env Rscript

# Validate a completed C3 Section 6 Phase B benchmark and estimate production
# cost from its cell-specific timings. This script does not run or alter any
# simulation; it only reads results and writes benchmark summaries.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source(file.path("scripts", "run_section6_new_scenarios.R"))

main <- function() {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  family <- tolower(as.character(args$family %||% ""))
  output_dir <- as.character(args$output_dir %||% "")
  base_seed <- suppressWarnings(as.integer(args$seed %||% NA_integer_))
  cores <- suppressWarnings(as.integer(args$cores %||% NA_integer_))
  runner_wall_seconds <- suppressWarnings(as.numeric(
    args$runner_wall_seconds %||% NA_real_
  ))
  completed_results_path <- as.character(args$completed_results %||% "")

  if (!family %in% c("normal", "lg")) {
    stop("`--family` must be either 'normal' or 'lg'.", call. = FALSE)
  }
  if (!nzchar(output_dir)) {
    stop("`--output_dir` is required.", call. = FALSE)
  }
  if (is.na(base_seed)) {
    stop("`--seed` must be a finite integer.", call. = FALSE)
  }
  if (is.na(cores) || cores != 32L) {
    stop("Phase B validation requires `--cores=32`.", call. = FALSE)
  }
  if (!is.finite(runner_wall_seconds) || runner_wall_seconds <= 0) {
    stop("`--runner_wall_seconds` must be positive.", call. = FALSE)
  }

  manifest_path <- file.path(output_dir, "manifest.csv")
  result_path <- file.path(output_dir, "raw_results.csv")
  if (!file.exists(manifest_path) || !file.exists(result_path)) {
    stop(
      sprintf("Phase B output is incomplete in '%s'.", output_dir),
      call. = FALSE
    )
  }

  output_lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(output_lock), add = TRUE)

  phase_b_design <- make_section6_design(
    family = family,
    dimensions = c(2L, 10L),
    n_values = c(50L, 100L, 200L, 400L),
    beta_values = c(0, 0.5, 1)
  )
  section6_validate_manifest_design(
    manifest_path = manifest_path,
    design = phase_b_design,
    M = 20L,
    B = 5000L,
    base_seed = base_seed,
    derivative_mc_size = 1000L,
    cvm_block_size = 50L
  )
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (!all(as.integer(manifest$cores) == cores)) {
    stop("Phase B manifest does not record 32 outer workers.", call. = FALSE)
  }

  results <- utils::read.csv(result_path, stringsAsFactors = FALSE)
  required <- names(empty_section6_results())
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop(
      sprintf("Phase B results lack columns: %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  expected_rows <- nrow(phase_b_design) * 20L
  expected_jobs <- merge(
    phase_b_design,
    data.frame(rep = seq_len(20L)),
    by = NULL
  )
  result_keys <- section6_design_key(results, include_rep = TRUE)
  expected_keys <- section6_design_key(expected_jobs, include_rep = TRUE)
  if (nrow(results) != expected_rows) {
    stop(sprintf(
      "Phase B produced %d rows; exactly %d were expected.",
      nrow(results), expected_rows
    ), call. = FALSE)
  }
  if (anyDuplicated(result_keys)) {
    stop("Phase B results contain duplicate replication keys.", call. = FALSE)
  }
  if (!setequal(result_keys, expected_keys)) {
    stop("Phase B result keys do not match the requested design.", call. = FALSE)
  }

  conforming <- results$status == "ok" &
    results$bootstrap_method_requested == "fast_multiplier" &
    results$bootstrap_method_effective == "fast_multiplier" &
    results$fast_multiplier_backend_requested == "cpp" &
    results$fast_multiplier_backend_effective == "cpp" &
    results$fast_multiplier_cpp_kernel_requested == "contiguous_double" &
    results$fast_multiplier_cpp_kernel_effective == "contiguous_double" &
    results$fast_multiplier_fuse_ks_cvm_requested &
    results$fast_multiplier_fuse_ks_cvm_effective &
    results$ks_grid == "sample_unique_distances" &
    !results$fallback_to_reestimated
  if (!isTRUE(all(conforming))) {
    failed <- which(is.na(conforming) | !conforming)
    stop(sprintf(
      "Phase B has %d nonconforming or failed rows (row indices: %s).",
      length(failed), paste(failed, collapse = ",")
    ), call. = FALSE)
  }

  expected_seed_data <- section6_seed(
    base_seed, results$design_id, results$rep, 0L
  )
  expected_seed_bootstrap <- section6_seed(
    base_seed, results$design_id, results$rep, 1L
  )
  expected_seed_derivative <- section6_seed(
    base_seed, results$design_id, results$rep, 2L
  )
  if (!identical(as.integer(results$seed_data), expected_seed_data) ||
      !identical(as.integer(results$seed_bootstrap), expected_seed_bootstrap) ||
      !identical(as.integer(results$seed_derivative), expected_seed_derivative)) {
    stop("Phase B result seeds do not match the requested design.", call. = FALSE)
  }

  elapsed <- as.numeric(results$elapsed_seconds)
  if (any(!is.finite(elapsed)) || any(elapsed < 0)) {
    stop("Phase B contains invalid elapsed times.", call. = FALSE)
  }
  results$cell_key <- section6_design_key(results)

  timing_groups <- split(results, results$cell_key)
  cell_timings <- do.call(rbind, lapply(timing_groups, function(x) {
    data.frame(
      scenario = x$scenario[[1L]],
      family = x$family[[1L]],
      d = x$d[[1L]],
      n = x$n[[1L]],
      beta = x$beta[[1L]],
      design_id = x$design_id[[1L]],
      benchmark_rows = nrow(x),
      elapsed_mean_seconds = mean(x$elapsed_seconds),
      elapsed_median_seconds = stats::median(x$elapsed_seconds),
      elapsed_p95_seconds = unname(stats::quantile(x$elapsed_seconds, 0.95)),
      elapsed_max_seconds = max(x$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  }))
  rownames(cell_timings) <- NULL
  cell_timings <- cell_timings[order(cell_timings$design_id), , drop = FALSE]
  if (nrow(cell_timings) != 48L ||
      !all(cell_timings$benchmark_rows == 20L)) {
    stop("Phase B does not contain exactly 20 timings in each of 48 cells.",
         call. = FALSE)
  }

  completed_counts <- integer(nrow(cell_timings))
  completed_rows <- 0L
  if (nzchar(completed_results_path)) {
    if (!file.exists(completed_results_path)) {
      stop(sprintf(
        "Completed production results do not exist: %s.",
        completed_results_path
      ), call. = FALSE)
    }
    completed <- utils::read.csv(
      completed_results_path,
      stringsAsFactors = FALSE
    )
    completed_missing <- setdiff(required, names(completed))
    if (length(completed_missing)) {
      stop(sprintf(
        "Completed production results lack columns: %s.",
        paste(completed_missing, collapse = ", ")
      ), call. = FALSE)
    }
    completed_conforming <- completed$status == "ok" &
      completed$family == family &
      completed$bootstrap_method_effective == "fast_multiplier" &
      completed$fast_multiplier_backend_effective == "cpp" &
      completed$fast_multiplier_fuse_ks_cvm_effective &
      completed$ks_grid == "sample_unique_distances" &
      !completed$fallback_to_reestimated
    if (!isTRUE(all(completed_conforming))) {
      stop("Completed production results contain nonconforming rows.",
           call. = FALSE)
    }
    completed_keys <- section6_design_key(completed, include_rep = TRUE)
    production_jobs <- merge(
      phase_b_design,
      data.frame(rep = seq_len(1000L)),
      by = NULL
    )
    production_keys <- section6_design_key(production_jobs, include_rep = TRUE)
    if (anyDuplicated(completed_keys)) {
      stop("Completed production results contain duplicate keys.",
           call. = FALSE)
    }
    if (!all(completed_keys %in% production_keys)) {
      stop("Completed production results contain keys outside production.",
           call. = FALSE)
    }
    completed_seed_ok <-
      as.integer(completed$seed_data) ==
        section6_seed(base_seed, completed$design_id, completed$rep, 0L) &
      as.integer(completed$seed_bootstrap) ==
        section6_seed(base_seed, completed$design_id, completed$rep, 1L) &
      as.integer(completed$seed_derivative) ==
        section6_seed(base_seed, completed$design_id, completed$rep, 2L)
    if (!isTRUE(all(completed_seed_ok))) {
      stop("Completed production results have incompatible seeds.",
           call. = FALSE)
    }
    completed_cell_counts <- table(section6_design_key(completed))
    cell_keys <- section6_design_key(cell_timings)
    matched_counts <- completed_cell_counts[cell_keys]
    matched_counts[is.na(matched_counts)] <- 0L
    completed_counts <- as.integer(matched_counts)
    completed_rows <- nrow(completed)
  }

  if (any(completed_counts < 0L) || any(completed_counts > 1000L)) {
    stop("Invalid completed-row counts in the production design.",
         call. = FALSE)
  }
  cell_timings$production_completed_rows <- completed_counts
  cell_timings$production_pending_rows <- 1000L - completed_counts
  cell_timings$estimated_pending_cpu_seconds <-
    cell_timings$production_pending_rows *
      cell_timings$elapsed_mean_seconds

  replication_elapsed_total <- sum(elapsed)
  effective_parallel_workers <-
    replication_elapsed_total / runner_wall_seconds
  parallel_efficiency <- effective_parallel_workers / cores
  if (!is.finite(parallel_efficiency) || parallel_efficiency <= 0) {
    stop("Could not compute a valid parallel efficiency.", call. = FALSE)
  }

  estimated_pending_cpu_seconds <-
    sum(cell_timings$estimated_pending_cpu_seconds)
  estimated_pending_wall_seconds <-
    estimated_pending_cpu_seconds / (cores * parallel_efficiency)
  pending_rows <- sum(cell_timings$production_pending_rows)

  summary <- data.frame(
    family = family,
    status = "ok",
    expected_rows = expected_rows,
    produced_rows = nrow(results),
    ok_rows = sum(results$status == "ok"),
    M = 20L,
    B = 5000L,
    dimensions = "2,10",
    n_values = "50,100,200,400",
    beta_values = "0,0.5,1",
    derivative_mc_size = 1000L,
    cvm_block_size = 50L,
    base_seed = base_seed,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    cores = cores,
    runner_wall_seconds = runner_wall_seconds,
    replication_elapsed_total_seconds = replication_elapsed_total,
    effective_parallel_workers = effective_parallel_workers,
    parallel_efficiency = parallel_efficiency,
    production_total_rows = 48000L,
    production_completed_rows = completed_rows,
    production_pending_rows = pending_rows,
    estimated_pending_cpu_hours = estimated_pending_cpu_seconds / 3600,
    estimated_pending_wall_hours = estimated_pending_wall_seconds / 3600,
    estimated_pending_wall_hours_with_25pct_margin =
      1.25 * estimated_pending_wall_seconds / 3600,
    estimate_basis =
      "48 cell means from Phase B, weighted by pending production rows",
    completed_results = if (nzchar(completed_results_path)) {
      completed_results_path
    } else {
      NA_character_
    },
    validated_at = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S %Z",
      tz = "Europe/Madrid"
    ),
    output_dir = output_dir,
    stringsAsFactors = FALSE
  )

  summary_path <- file.path(output_dir, "phase_b_validation_summary.csv")
  cell_path <- file.path(output_dir, "phase_b_cell_timings.csv")
  section6_write_atomic_csv(summary, summary_path)
  section6_write_atomic_csv(cell_timings, cell_path)

  cat("Phase B validation passed.\n")
  cat(sprintf("family: %s\n", family))
  cat(sprintf("results: %d/%d conforming\n", nrow(results), expected_rows))
  cat(sprintf("runner wall seconds: %.3f\n", runner_wall_seconds))
  cat(sprintf(
    "parallel efficiency: %.4f (%0.2f effective workers of %d)\n",
    parallel_efficiency, effective_parallel_workers, cores
  ))
  cat(sprintf(
    "production rows: completed=%d pending=%d total=48000\n",
    completed_rows, pending_rows
  ))
  cat(sprintf(
    "estimated pending production: %.3f CPU-hours, %.3f wall-hours",
    estimated_pending_cpu_seconds / 3600,
    estimated_pending_wall_seconds / 3600
  ))
  cat(sprintf(
    " (%.3f wall-hours with 25%% margin)\n",
    1.25 * estimated_pending_wall_seconds / 3600
  ))
  cat(sprintf("summary: %s\n", summary_path))
  cat(sprintf("cell timings: %s\n", cell_path))
}

main()
