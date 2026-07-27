#!/usr/bin/env Rscript

# Validate a completed C3 Section 6 Phase A benchmark without changing the
# statistical procedure. The validator requires exactly the 12 conforming
# results implied by the Phase A design and writes a compact timing summary.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source(file.path("scripts", "run_section6_new_scenarios.R"))

main <- function() {
args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
family <- tolower(as.character(args$family %||% ""))
output_dir <- as.character(args$output_dir %||% "")
base_seed <- suppressWarnings(as.integer(args$seed %||% NA_integer_))

if (!family %in% c("normal", "lg")) {
  stop("`--family` must be either 'normal' or 'lg'.", call. = FALSE)
}
if (!nzchar(output_dir)) {
  stop("`--output_dir` is required.", call. = FALSE)
}
if (is.na(base_seed)) {
  stop("`--seed` must be a finite integer.", call. = FALSE)
}

manifest_path <- file.path(output_dir, "manifest.csv")
result_path <- file.path(output_dir, "raw_results.csv")
if (!file.exists(manifest_path) || !file.exists(result_path)) {
  stop(
    sprintf("Phase A output is incomplete in '%s'.", output_dir),
    call. = FALSE
  )
}

output_lock <- section6_acquire_output_lock(output_dir)
on.exit(section6_release_output_lock(output_lock), add = TRUE)

phase_a_design <- make_section6_design(
  family = family,
  dimensions = 10L,
  n_values = 400L,
  beta_values = c(0, 0.5, 1)
)
section6_validate_manifest_design(
  manifest_path = manifest_path,
  design = phase_a_design,
  M = 2L,
  B = 5000L,
  base_seed = base_seed,
  derivative_mc_size = 1000L,
  cvm_block_size = 50L
)

results <- utils::read.csv(result_path, stringsAsFactors = FALSE)
required <- names(empty_section6_results())
missing <- setdiff(required, names(results))
if (length(missing)) {
  stop(
    sprintf("Phase A results lack columns: %s.", paste(missing, collapse = ", ")),
    call. = FALSE
  )
}

expected_rows <- nrow(phase_a_design) * 2L
result_keys <- section6_design_key(results, include_rep = TRUE)
expected_jobs <- merge(
  phase_a_design,
  data.frame(rep = seq_len(2L)),
  by = NULL
)
expected_keys <- section6_design_key(expected_jobs, include_rep = TRUE)

if (nrow(results) != expected_rows) {
  stop(sprintf(
    "Phase A produced %d rows; exactly %d were expected.",
    nrow(results), expected_rows
  ), call. = FALSE)
}
if (anyDuplicated(result_keys)) {
  stop("Phase A results contain duplicate replication keys.", call. = FALSE)
}
if (!setequal(result_keys, expected_keys)) {
  stop("Phase A result keys do not match the requested design.", call. = FALSE)
}

conforming <- results$status == "ok" &
  results$bootstrap_method_effective == "fast_multiplier" &
  results$fast_multiplier_backend_effective == "cpp" &
  results$fast_multiplier_fuse_ks_cvm_effective &
  results$ks_grid == "sample_unique_distances" &
  !results$fallback_to_reestimated
if (!isTRUE(all(conforming))) {
  failed <- which(is.na(conforming) | !conforming)
  stop(sprintf(
    "Phase A has %d nonconforming or failed rows (row indices: %s).",
    length(failed), paste(failed, collapse = ",")
  ), call. = FALSE)
}

elapsed <- as.numeric(results$elapsed_seconds)
if (any(!is.finite(elapsed)) || any(elapsed < 0)) {
  stop("Phase A contains invalid elapsed times.", call. = FALSE)
}

slurm_cpus <- suppressWarnings(as.integer(Sys.getenv(
  "SLURM_CPUS_PER_TASK",
  unset = NA_character_
)))
summary <- data.frame(
  family = family,
  status = "ok",
  expected_rows = expected_rows,
  produced_rows = nrow(results),
  ok_rows = sum(results$status == "ok"),
  M = 2L,
  B = 5000L,
  dimensions = "10",
  n_values = "400",
  beta_values = "0,0.5,1",
  derivative_mc_size = 1000L,
  cvm_block_size = 50L,
  base_seed = base_seed,
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  slurm_cpus_per_task = slurm_cpus,
  replication_elapsed_total_seconds = sum(elapsed),
  replication_elapsed_mean_seconds = mean(elapsed),
  replication_elapsed_median_seconds = stats::median(elapsed),
  replication_elapsed_p95_seconds = unname(stats::quantile(elapsed, 0.95)),
  replication_elapsed_max_seconds = max(elapsed),
  validated_at = format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S %Z",
    tz = "Europe/Madrid"
  ),
  output_dir = output_dir,
  stringsAsFactors = FALSE
)

summary_path <- file.path(output_dir, "phase_a_validation_summary.csv")
section6_write_atomic_csv(summary, summary_path)

cat("Phase A validation passed.\n")
cat(sprintf("family: %s\n", family))
cat(sprintf("results: %d/%d conforming\n", nrow(results), expected_rows))
cat(sprintf(
  "replication elapsed seconds: total=%.3f mean=%.3f median=%.3f p95=%.3f max=%.3f\n",
  sum(elapsed), mean(elapsed), stats::median(elapsed),
  unname(stats::quantile(elapsed, 0.95)), max(elapsed)
))
cat(sprintf("summary: %s\n", summary_path))
}

main()
