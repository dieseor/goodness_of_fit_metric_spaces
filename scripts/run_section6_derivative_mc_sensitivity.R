#!/usr/bin/env Rscript

# Paired derivative-Monte-Carlo sensitivity experiment for the Section 6
# Normal and logistic-Gaussian d = 10 null calibrations.  The complete design
# is built before filtering, preserving design_id and hence the data and
# multiplier seeds of the corresponding final campaign.  Only the auxiliary
# derivative sample size is changed.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_section6_new_scenarios.R")

run_section6_derivative_sensitivity <- function(family,
                                                output_dir,
                                                M = 1000L,
                                                B = 5000L,
                                                cores = 6L,
                                                base_seed = 20260727L,
                                                derivative_mc_size = 10000L,
                                                cvm_block_size = 50L,
                                                checkpoint_results = 16L,
                                                show_progress = TRUE) {
  family <- tolower(as.character(family))
  if (!family %in% c("normal", "lg")) {
    stop("`family` must be either 'normal' or 'lg'.", call. = FALSE)
  }
  M <- as.integer(M)
  B <- as.integer(B)
  cores <- as.integer(cores)
  derivative_mc_size <- as.integer(derivative_mc_size)
  cvm_block_size <- as.integer(cvm_block_size)
  checkpoint_results <- as.integer(checkpoint_results)
  if (any(!is.finite(c(M, B, cores, derivative_mc_size, cvm_block_size, checkpoint_results))) ||
      any(c(M, B, cores, derivative_mc_size, cvm_block_size, checkpoint_results) < 1L)) {
    stop("All integer controls must be strictly positive.", call. = FALSE)
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Outer parallelism requires a Unix platform.", call. = FALSE)
  }

  # Do not construct the filtered design directly: design_id determines all
  # random-number streams in the final runner.  Keeping these identifiers
  # makes this a paired comparison with the existing d = 10 campaign.
  full_design <- make_section6_design(
    family = family,
    dimensions = c(2L, 10L),
    n_values = c(50L, 100L, 200L, 400L),
    beta_values = c(0, 0.5, 1)
  )
  design <- full_design[full_design$d == 10L & full_design$beta == 0, , drop = FALSE]
  if (nrow(design) != 8L) stop("Internal error: expected eight d = 10 null cells.")

  ensure_distance_profile_cpp_loaded()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(output_lock), add = TRUE)

  result_path <- file.path(output_dir, "raw_results.csv")
  manifest_path <- file.path(output_dir, "manifest.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  log_path <- file.path(output_dir, "run.log")
  section6_validate_manifest_design(
    manifest_path, design, M = M, B = B, base_seed = base_seed,
    derivative_mc_size = derivative_mc_size, cvm_block_size = cvm_block_size
  )
  if (!file.exists(manifest_path)) {
    section6_write_atomic_csv(
      section6_make_manifest(design, M, B, cores, base_seed, derivative_mc_size, cvm_block_size),
      manifest_path
    )
    writeLines(c(
      "purpose: paired d=10 null-calibration sensitivity to derivative_mc_size",
      "paired_reference_design: Section 6 final design with dimensions=2,10; n=50,100,200,400; beta=0,0.5,1",
      "retained_cells: d=10; beta=0",
      "unchanged: data seed, multiplier seed, M, B, sample KS, C++ fused fast bootstrap",
      sprintf("changed: derivative_mc_size=%d", derivative_mc_size)
    ), file.path(output_dir, "sensitivity_manifest.txt"))
  }

  existing <- if (file.exists(result_path)) {
    utils::read.csv(result_path, stringsAsFactors = FALSE)
  } else {
    empty_section6_results()
  }
  jobs <- merge(design, data.frame(rep = seq_len(M)), by = NULL)
  done <- if (nrow(existing)) {
    section6_design_key(existing[existing$status == "ok", , drop = FALSE], include_rep = TRUE)
  } else character()
  pending <- jobs[!section6_design_key(jobs, include_rep = TRUE) %in% done, , drop = FALSE]
  started <- Sys.time()
  completed_before <- nrow(jobs) - nrow(pending)
  cat(sprintf(
    "%s derivative-MC sensitivity family=%s d=10 beta=0 M=%d B=%d Nderiv=%d cores=%d pending=%d\n",
    format(started, tz = "Europe/Madrid"), family, M, B, derivative_mc_size, cores, nrow(pending)
  ), file = log_path, append = TRUE)
  section6_write_status(status_path, family, nrow(jobs), completed_before, existing, started, cores)
  if (isTRUE(show_progress)) section6_progress(completed_before, nrow(jobs), started, existing, cores)
  if (!nrow(pending)) return(invisible(list(results = existing, summary = summarize_section6_results(existing))))

  section6_run_dynamic_queue(
    pending = pending, B = B, base_seed = base_seed,
    derivative_mc_size = derivative_mc_size, cvm_block_size = cvm_block_size,
    cores = cores, checkpoint_results = checkpoint_results,
    on_checkpoint = function(rows, finished) {
      existing <<- rbind(existing, rows)
      existing <- existing[order(existing$design_id, existing$rep), , drop = FALSE]
      section6_write_atomic_csv(existing, result_path)
      section6_write_atomic_csv(summarize_section6_results(existing), summary_path)
      completed <- completed_before + finished
      section6_write_status(status_path, family, nrow(jobs), completed, existing, started, cores)
      cat(sprintf("%s completed=%d/%d\n", format(Sys.time(), tz = "Europe/Madrid"), completed, nrow(jobs)), file = log_path, append = TRUE)
    },
    on_progress = function(finished) {
      if (isTRUE(show_progress)) {
        section6_progress(completed_before + finished, nrow(jobs), started, existing, cores)
      }
    }
  )
  invisible(list(results = existing, summary = summarize_section6_results(existing)))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  family <- tolower(as.character(args$family %||% "normal"))
  M <- as.integer(args$M %||% 1000L)
  B <- as.integer(args$B %||% 5000L)
  output_dir <- as.character(args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios",
    sprintf("sensitivity_%s_d10_beta0_M%d_B%d_Nderiv%s", family, M, B,
      args$derivative_mc_size %||% "10000")
  ))
  run_section6_derivative_sensitivity(
    family = family, output_dir = output_dir, M = M, B = B,
    cores = as.integer(args$cores %||% 6L),
    base_seed = as.integer(args$seed %||% 20260727L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 10000L),
    cvm_block_size = as.integer(args$cvm_block_size %||% 50L),
    checkpoint_results = as.integer(args$checkpoint_results %||% 16L),
    show_progress = !identical(tolower(as.character(args$show_progress %||% "true")), "false")
  )
}
