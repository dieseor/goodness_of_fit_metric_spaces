#!/usr/bin/env Rscript

# Resumable null-calibration comparison between the fast multiplier and the
# reestimated bootstrap for the two closed Normal scenarios of Section 6:
#   1. restricted mean-aligned single-spiked Normal;
#   2. Normal with covariance sigma^2 I_d.
# This diagnostic is intentionally isolated from the existing simulation
# runners and does not alter either model implementation.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    fields <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[fields[[1L]]]] <- if (length(fields) == 1L) "TRUE" else paste(fields[-1L], collapse = "=")
  }
  out
}

parse_csv_integer <- function(value, default, name, lower = 1L) {
  if (is.null(value) || !nzchar(value)) return(as.integer(default))
  out <- as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!is.finite(out)) || any(out < lower)) {
    stop("`--", name, "` must contain integers at least ", lower, ".")
  }
  sort(unique(out))
}

parse_csv_character <- function(value, default, name, choices) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  aliases <- c("1.1" = "normal_1_single_spiked", "1.2" = "normal_2_sigma_Id")
  alias_values <- unname(aliases[out])
  out[!is.na(alias_values)] <- alias_values[!is.na(alias_values)]
  if (!length(out) || any(!out %in% choices)) {
    stop("`--", name, "` must be a comma-separated subset of ",
         paste(c(choices, names(aliases)), collapse = ", "), ".")
  }
  unique(out)
}

atomic_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not atomically write '", path, "'.")
}

result_key <- function(x) {
  paste(x$scenario, x$d, x$n, x$replication, sep = "|")
}

empty_results <- function() {
  data.frame(
    scenario = character(), d = integer(), n = integer(), replication = integer(),
    seed_data = integer(), seed_bootstrap = integer(), seed_derivative = integer(),
    status = character(), error_message = character(), warning_message = character(),
    fast_ks_statistic = numeric(), fast_cvm_statistic = numeric(),
    slow_ks_statistic = numeric(), slow_cvm_statistic = numeric(),
    fast_ks_pvalue = numeric(), fast_cvm_pvalue = numeric(),
    slow_ks_pvalue = numeric(), slow_cvm_pvalue = numeric(),
    fast_ks_reject = logical(), fast_cvm_reject = logical(),
    slow_ks_reject = logical(), slow_cvm_reject = logical(),
    fast_method = character(), fast_backend = character(), fast_kernel = character(),
    fast_fused = logical(), fast_derivative_method = character(),
    fast_derivative_mc_size = integer(), slow_method = character(),
    slow_fused_ks_cvm = logical(),
    elapsed_fast_seconds = numeric(), elapsed_slow_seconds = numeric(),
    elapsed_pair_seconds = numeric(), stringsAsFactors = FALSE
  )
}

summarize_results <- function(x) {
  if (!nrow(x)) return(data.frame())
  groups <- split(x, interaction(x$scenario, x$d, x$n, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    ok <- z[z$status == "ok", , drop = FALSE]
    data.frame(
      scenario = z$scenario[[1L]], d = z$d[[1L]], n = z$n[[1L]],
      completed_pairs = nrow(ok), failed_pairs = sum(z$status != "ok"),
      fast_ks_rejection_percent = if (nrow(ok)) 100 * mean(ok$fast_ks_reject) else NA_real_,
      fast_cvm_rejection_percent = if (nrow(ok)) 100 * mean(ok$fast_cvm_reject) else NA_real_,
      slow_ks_rejection_percent = if (nrow(ok)) 100 * mean(ok$slow_ks_reject) else NA_real_,
      slow_cvm_rejection_percent = if (nrow(ok)) 100 * mean(ok$slow_cvm_reject) else NA_real_,
      mean_abs_pvalue_difference_ks = if (nrow(ok)) mean(abs(ok$fast_ks_pvalue - ok$slow_ks_pvalue)) else NA_real_,
      mean_abs_pvalue_difference_cvm = if (nrow(ok)) mean(abs(ok$fast_cvm_pvalue - ok$slow_cvm_pvalue)) else NA_real_,
      max_abs_observed_difference_ks = if (nrow(ok)) max(abs(ok$fast_ks_statistic - ok$slow_ks_statistic)) else NA_real_,
      max_abs_observed_difference_cvm = if (nrow(ok)) max(abs(ok$fast_cvm_statistic - ok$slow_cvm_statistic)) else NA_real_,
      mean_fast_seconds = if (nrow(ok)) mean(ok$elapsed_fast_seconds) else NA_real_,
      mean_slow_seconds = if (nrow(ok)) mean(ok$elapsed_slow_seconds) else NA_real_,
      all_slow_fused_ks_cvm = if (nrow(ok)) all(ok$slow_fused_ks_cvm) else NA,
      stringsAsFactors = FALSE
    )
  }))
}

write_progress <- function(path, total, results, started, cores) {
  completed <- sum(results$status == "ok")
  pair_time <- if (completed) mean(results$elapsed_pair_seconds[results$status == "ok"]) else NA_real_
  eta <- if (is.finite(pair_time)) (total - completed) * pair_time / cores else NA_real_
  writeLines(c(
    sprintf("updated: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("completed_pairs: %d/%d", completed, total),
    sprintf("failed_pairs: %d", sum(results$status != "ok")),
    sprintf("elapsed_seconds: %.1f", as.numeric(difftime(Sys.time(), started, units = "secs"))),
    sprintf("mean_seconds_per_completed_pair: %s", if (is.finite(pair_time)) format(round(pair_time, 2), nsmall = 2) else "NA"),
    sprintf("eta_seconds: %s", if (is.finite(eta)) as.character(round(eta)) else "NA"),
    sprintf("outer_cores: %d; inner_cores: 1", cores)
  ), path)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
M <- as.integer(args$M %||% 1000L)
B <- as.integer(args$B %||% 5000L)
cores <- as.integer(args$cores %||% 8L)
dimensions <- parse_csv_integer(args$dimensions, c(2L, 5L), "dimensions", lower = 2L)
n_values <- parse_csv_integer(args$n_values %||% args$n, c(200L, 1000L), "n", lower = 2L)
derivative_mc_size <- as.integer(args$derivative_mc_size %||% 10000L)
cvm_block_size <- as.integer(args$cvm_block_size %||% 50L)
checkpoint_pairs <- as.integer(args$checkpoint_pairs %||% cores)
seed <- as.integer(args$seed %||% 20260826L)
output_dir <- args$output_dir %||% NULL

if (any(!is.finite(c(M, B, cores, derivative_mc_size, cvm_block_size, checkpoint_pairs, seed))) ||
    M < 1L || B < 1L || cores < 1L || derivative_mc_size < 1L ||
    cvm_block_size < 1L || checkpoint_pairs < 1L) {
  stop("M, B, cores, auxiliary size, block size, checkpoint size and seed must be positive.")
}
if (.Platform$OS.type != "unix" && cores > 1L) stop("Outer parallelism requires a Unix platform.")

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
source(file.path("bootstrap", "normal_sigma_Id_bootstrap.R"))

scenario_catalog <- data.frame(
  scenario = c("normal_1_single_spiked", "normal_2_sigma_Id"),
  scenario_code = c(1L, 2L),
  stringsAsFactors = FALSE
)
scenarios <- parse_csv_character(
  args$scenarios %||% args$scenario,
  scenario_catalog$scenario,
  "scenarios",
  scenario_catalog$scenario
)
scenario_catalog <- scenario_catalog[scenario_catalog$scenario %in% scenarios, , drop = FALSE]
if (is.null(output_dir)) {
  output_dir <- file.path(
    "simulation_results", "section6_new_scenarios",
    sprintf("closed_normal_fast_vs_reestimated_%s_d%s_n%s_M%d_B%d_Nderiv%d",
            paste(scenarios, collapse = "_"), paste(dimensions, collapse = "_"),
            paste(n_values, collapse = "_"), M, B, derivative_mc_size)
  )
}
design <- merge(
  merge(scenario_catalog, expand.grid(d = dimensions, n = n_values), by = NULL),
  data.frame(replication = seq_len(M)), by = NULL
)
design <- design[order(design$scenario, design$d, design$n, design$replication), , drop = FALSE]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
lock_path <- file.path(output_dir, ".closed_normal_fast_vs_reestimated.lock")
if (!dir.create(lock_path, showWarnings = FALSE)) {
  stop("Output directory is locked by another process: ", output_dir)
}
on.exit(unlink(lock_path, recursive = TRUE, force = TRUE), add = TRUE)

manifest_path <- file.path(output_dir, "manifest.csv")
results_path <- file.path(output_dir, "matched_results.csv")
summary_path <- file.path(output_dir, "summary.csv")
progress_path <- file.path(output_dir, "progress_status.txt")
if (file.exists(manifest_path)) {
  old <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  expected <- transform(design, M = M, B = B, derivative_mc_size = derivative_mc_size,
                        cvm_block_size = cvm_block_size, base_seed = seed,
                        same_multiplier_seed = TRUE, reestimated_fuse_ks_cvm = TRUE,
                        fast_cache_corrections = "true")
  if (!isTRUE(all.equal(old, expected, check.attributes = FALSE))) {
    stop("Existing manifest is incompatible. Use a new output directory.")
  }
} else {
  manifest <- transform(design, M = M, B = B, derivative_mc_size = derivative_mc_size,
                        cvm_block_size = cvm_block_size, base_seed = seed,
                        same_multiplier_seed = TRUE, reestimated_fuse_ks_cvm = TRUE,
                        fast_cache_corrections = "true")
  atomic_write_csv(manifest, manifest_path)
}

results <- if (file.exists(results_path)) utils::read.csv(results_path, stringsAsFactors = FALSE) else empty_results()
if (!all(names(empty_results()) %in% names(results))) stop("Existing result file has an incompatible schema.")
done <- if (nrow(results)) result_key(results[results$status == "ok", , drop = FALSE]) else character()
pending <- design[!result_key(design) %in% done, , drop = FALSE]

generate_null_data <- function(row, data_seed) {
  set.seed(data_seed)
  d <- as.integer(row$d)
  n <- as.integer(row$n)
  if (identical(as.character(row$scenario), "normal_1_single_spiked")) {
    theta <- rep(1 / sqrt(d), d)
    return(rrestricted_spiked_normal(n, theta = theta, lambda = 2))
  }
  matrix(stats::rnorm(n * d), nrow = n, ncol = d)
}

run_one <- function(row) {
  started <- proc.time()[["elapsed"]]
  seed_for <- function(stream) {
    as.integer((as.numeric(seed) + 1000003 * as.integer(row$scenario_code) +
      10007 * as.integer(row$d) + 10000019 * as.integer(row$n) +
      1009 * as.integer(row$replication) + 97 * as.integer(stream)) %% 2147483647) + 1L
  }
  data_seed <- seed_for(0L)
  bootstrap_seed <- seed_for(1L)
  derivative_seed <- seed_for(2L)
  base <- data.frame(
    scenario = as.character(row$scenario), d = as.integer(row$d), n = as.integer(row$n),
    replication = as.integer(row$replication), seed_data = data_seed,
    seed_bootstrap = bootstrap_seed, seed_derivative = derivative_seed,
    stringsAsFactors = FALSE
  )
  warnings <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  tryCatch({
    x <- generate_null_data(row, data_seed)
    common <- list(
      data = x, null = list(type = "composite"), statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = 0.05,
      n_cores = 1L, keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                                bootstrap_thetas = FALSE),
      control = list(derivative_method = "score_mc", derivative_mc_size = derivative_mc_size,
                     derivative_mc_seed = derivative_seed, fast_multiplier_cvm_block_size = cvm_block_size,
                     fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
                     fast_multiplier_fuse_ks_cvm = TRUE, fast_multiplier_cache_corrections = "true",
                     fast_multiplier_stream_chunk_size = 100L,
                     reestimated_fuse_ks_cvm = TRUE),
      distance_profile_backend = "r"
    )
    runner <- if (identical(base$scenario[[1L]], "normal_1_single_spiked")) {
      multiplier_bootstrap_restricted_spiked_normal
    } else multiplier_bootstrap_normal_sigma_Id
    fast_started <- proc.time()[["elapsed"]]
    fast <- capture(do.call(runner, c(common, list(
      seed = bootstrap_seed, bootstrap_method = "fast_multiplier",
      fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
      fuse_ks_cvm = TRUE, cache_block_corrections = "true"
    ))))
    elapsed_fast <- proc.time()[["elapsed"]] - fast_started
    slow_started <- proc.time()[["elapsed"]]
    slow <- capture(do.call(runner, c(common, list(
      seed = bootstrap_seed, bootstrap_method = "reestimated"
    ))))
    elapsed_slow <- proc.time()[["elapsed"]] - slow_started
    diagnostics <- fast$diagnostics
    is_conforming_fast <- identical(diagnostics$effective_bootstrap_method, "fast_multiplier") &&
      !isTRUE(diagnostics$fallback_to_reestimated) &&
      identical(diagnostics$fast_multiplier_backend_effective, "cpp") &&
      identical(diagnostics$fast_multiplier_cpp_kernel_effective, "contiguous_double") &&
      isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective) &&
      identical(diagnostics$derivative_method, "score_mc") &&
      identical(as.integer(diagnostics$derivative_mc_size), derivative_mc_size)
    if (!is_conforming_fast) stop("The requested conforming fast configuration was not effective.")
    cbind(base, data.frame(
      status = "ok", error_message = NA_character_,
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      fast_ks_statistic = fast$inference$ks$observed, fast_cvm_statistic = fast$inference$cvm$observed,
      slow_ks_statistic = slow$inference$ks$observed, slow_cvm_statistic = slow$inference$cvm$observed,
      fast_ks_pvalue = fast$inference$ks$p_value, fast_cvm_pvalue = fast$inference$cvm$p_value,
      slow_ks_pvalue = slow$inference$ks$p_value, slow_cvm_pvalue = slow$inference$cvm$p_value,
      fast_ks_reject = fast$inference$ks$reject, fast_cvm_reject = fast$inference$cvm$reject,
      slow_ks_reject = slow$inference$ks$reject, slow_cvm_reject = slow$inference$cvm$reject,
      fast_method = diagnostics$effective_bootstrap_method,
      fast_backend = diagnostics$fast_multiplier_backend_effective,
      fast_kernel = diagnostics$fast_multiplier_cpp_kernel_effective,
      fast_fused = isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective),
      fast_derivative_method = diagnostics$derivative_method,
      fast_derivative_mc_size = as.integer(diagnostics$derivative_mc_size),
      slow_method = slow$diagnostics$effective_bootstrap_method %||% "reestimated",
      slow_fused_ks_cvm = isTRUE(slow$diagnostics$reestimated_fuse_ks_cvm_effective),
      elapsed_fast_seconds = elapsed_fast, elapsed_slow_seconds = elapsed_slow,
      elapsed_pair_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cbind(base, data.frame(
      status = "error", error_message = conditionMessage(e),
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      fast_ks_statistic = NA_real_, fast_cvm_statistic = NA_real_, slow_ks_statistic = NA_real_, slow_cvm_statistic = NA_real_,
      fast_ks_pvalue = NA_real_, fast_cvm_pvalue = NA_real_, slow_ks_pvalue = NA_real_, slow_cvm_pvalue = NA_real_,
      fast_ks_reject = NA, fast_cvm_reject = NA, slow_ks_reject = NA, slow_cvm_reject = NA,
      fast_method = NA_character_, fast_backend = NA_character_, fast_kernel = NA_character_, fast_fused = NA,
      fast_derivative_method = NA_character_, fast_derivative_mc_size = NA_integer_, slow_method = NA_character_,
      slow_fused_ks_cvm = NA,
      elapsed_fast_seconds = NA_real_, elapsed_slow_seconds = NA_real_,
      elapsed_pair_seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE
    ))
  })
}

message(sprintf("Closed Normal paired fast/reestimated comparison: %d pending / %d total pairs; M=%d, B=%d, cores=%d.",
                nrow(pending), nrow(design), M, B, cores))
started_all <- Sys.time()
write_progress(progress_path, nrow(design), results, started_all, cores)
while (nrow(pending)) {
  take <- min(checkpoint_pairs, nrow(pending))
  batch <- pending[seq_len(take), , drop = FALSE]
  pending <- pending[-seq_len(take), , drop = FALSE]
  rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) run_one(batch[i, , drop = FALSE]),
                             mc.cores = min(cores, nrow(batch)), mc.preschedule = FALSE, mc.set.seed = FALSE)
  rows <- do.call(rbind, rows)
  keys <- result_key(rows)
  if (nrow(results)) results <- results[!result_key(results) %in% keys, , drop = FALSE]
  results <- rbind(results, rows)
  results <- results[order(results$scenario, results$d, results$n, results$replication), , drop = FALSE]
  atomic_write_csv(results, results_path)
  atomic_write_csv(summarize_results(results), summary_path)
  write_progress(progress_path, nrow(design), results, started_all, cores)
  message(sprintf("completed %d/%d pairs", sum(results$status == "ok"), nrow(design)))
}

message("Results: ", normalizePath(results_path, mustWork = FALSE))
message("Summary: ", normalizePath(summary_path, mustWork = FALSE))
