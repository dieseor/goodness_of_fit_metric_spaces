#!/usr/bin/env Rscript

# Small, resumable matched fast-versus-reestimated diagnostic for the
# restricted mean-aligned single-spiked Gaussian model.  It is deliberately
# separate from the general Gaussian and the main restricted-model runners.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  value <- args[startsWith(args, key)]
  if (!length(value)) return(default)
  sub(key, "", value[[1L]], fixed = TRUE)
}
parse_csv_integer <- function(value, name) {
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (!length(out) || any(!is.finite(out)) || any(out < 2L)) {
    stop("`--", name, "` must contain integers at least two.")
  }
  out
}

B <- as.integer(arg_value("B", "39"))
cores <- as.integer(arg_value("cores", "8"))
n_values <- parse_csv_integer(arg_value("n", "50,200"), "n")
lambda <- as.numeric(arg_value("lambda", "2"))
axis_norm <- as.numeric(arg_value("axis_norm", "0.75"))
diagonal_norm <- as.numeric(arg_value("diagonal_norm", "1.5"))
derivative_mc_size <- as.integer(arg_value("derivative_mc_size", "1000"))
cvm_block_size <- as.integer(arg_value("cvm_block_size", "50"))
seed <- as.integer(arg_value("seed", "20260832"))
output_dir <- arg_value(
  "output_dir",
  "simulation_results/restricted_spiked_normal_fast_vs_reestimated_lambda2_B39"
)
if (any(!is.finite(c(B, cores, lambda, axis_norm, diagonal_norm,
                     derivative_mc_size, cvm_block_size, seed))) ||
    B < 1L || cores < 1L || lambda <= 0 || axis_norm <= 0 ||
    diagonal_norm <= 0 || derivative_mc_size < 1L || cvm_block_size < 1L) {
  stop("B, cores, lambda, norms, and auxiliary sizes must be strictly positive.")
}

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theta_for <- function(d, theta_name) {
  if (identical(theta_name, "axis_moderate")) {
    c(axis_norm, rep(0, d - 1L))
  } else {
    rep(diagonal_norm / sqrt(d), d)
  }
}

design <- expand.grid(
  d = c(2L, 5L),
  theta_name = c("axis_moderate", "diagonal_large"),
  n = n_values,
  stringsAsFactors = FALSE
)
design$job_id <- seq_len(nrow(design))
utils::write.csv(design, file.path(output_dir, "design.csv"), row.names = FALSE)
writeLines(c(
  "matched fast versus reestimated diagnostic: restricted single-spiked normal",
  sprintf("lambda=%g; B=%d; outer cores=%d; seed=%d", lambda, B, cores, seed),
  sprintf("d=2,5; n=%s", paste(n_values, collapse = ",")),
  sprintf("theta designs: %g e1; %g (1,...,1)/sqrt(d)", axis_norm, diagonal_norm),
  "Each pair uses the identical dataset and the same bootstrap seed.",
  "KS: sample unique-distance grid; fast: C++ contiguous_double fused KS-CvM.",
  sprintf("derivative=score_mc, auxiliary size=%d; CvM block size=%d.", derivative_mc_size, cvm_block_size)
), file.path(output_dir, "manifest.txt"))

run_one <- function(row) {
  started <- proc.time()[["elapsed"]]
  job_seed <- seed + 100003L * as.integer(row$job_id)
  set.seed(job_seed)
  d <- as.integer(row$d)
  theta <- theta_for(d, as.character(row$theta_name))
  x <- rrestricted_spiked_normal(as.integer(row$n), theta, lambda)
  base <- data.frame(
    job_id = as.integer(row$job_id), d = d, theta_name = as.character(row$theta_name),
    theta_true = paste(theta, collapse = ";"), theta_norm_true = sqrt(sum(theta^2)),
    n = as.integer(row$n), lambda = lambda, seed = job_seed,
    stringsAsFactors = FALSE
  )
  warnings <- character()
  tryCatch({
    common <- list(
      data = x, null = list(type = "composite"), statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = 0.05,
      n_cores = 1L, seed = job_seed + 1L,
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                  bootstrap_thetas = FALSE),
      control = list(
        derivative_method = "score_mc", derivative_mc_size = derivative_mc_size,
        derivative_mc_seed = job_seed + 2L,
        fast_multiplier_cvm_block_size = cvm_block_size,
        fast_multiplier_backend = "cpp",
        fast_multiplier_cpp_kernel = "contiguous_double",
        fast_multiplier_fuse_ks_cvm = TRUE,
        fast_multiplier_cache_corrections = "auto",
        progress_bar = FALSE
      ),
      distance_profile_backend = "r"
    )
    capture <- function(expression) withCallingHandlers(expression, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
    fast <- capture(do.call(multiplier_bootstrap_restricted_spiked_normal, c(
      common,
      list(bootstrap_method = "fast_multiplier", fast_multiplier_backend = "cpp",
           fast_multiplier_cpp_kernel = "contiguous_double", fuse_ks_cvm = TRUE,
           cache_block_corrections = "auto")
    )))
    slow <- capture(do.call(multiplier_bootstrap_restricted_spiked_normal, c(
      common, list(bootstrap_method = "reestimated")
    )))
    cbind(base, data.frame(
      status = "ok", error_message = NA_character_,
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      fast_ks_pvalue = fast$inference$ks$p_value,
      fast_cvm_pvalue = fast$inference$cvm$p_value,
      reestimated_ks_pvalue = slow$inference$ks$p_value,
      reestimated_cvm_pvalue = slow$inference$cvm$p_value,
      difference_ks = fast$inference$ks$p_value - slow$inference$ks$p_value,
      difference_cvm = fast$inference$cvm$p_value - slow$inference$cvm$p_value,
      fast_backend = fast$diagnostics$fast_multiplier_backend_effective,
      fast_kernel = fast$diagnostics$fast_multiplier_cpp_kernel_effective,
      fast_fused = fast$diagnostics$fast_multiplier_fuse_ks_cvm_effective,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    ))
  }, error = function(error) {
    cbind(base, data.frame(
      status = "error", error_message = conditionMessage(error),
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      fast_ks_pvalue = NA_real_, fast_cvm_pvalue = NA_real_,
      reestimated_ks_pvalue = NA_real_, reestimated_cvm_pvalue = NA_real_,
      difference_ks = NA_real_, difference_cvm = NA_real_,
      fast_backend = NA_character_, fast_kernel = NA_character_, fast_fused = NA,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    ))
  })
}

results_path <- file.path(output_dir, "matched_results.csv")
existing <- if (file.exists(results_path)) utils::read.csv(results_path, stringsAsFactors = FALSE) else data.frame()
done <- if (nrow(existing)) as.integer(existing$job_id) else integer()
pending <- design[!design$job_id %in% done, , drop = FALSE]
message(sprintf("Restricted-spiked matched fast/reestimated: %d pending / %d pairs; B=%d; outer cores=%d.",
                nrow(pending), nrow(design), B, cores))
started_all <- proc.time()[["elapsed"]]
if (nrow(pending)) {
  for (begin in seq.int(1L, nrow(pending), by = cores)) {
    end <- min(nrow(pending), begin + cores - 1L)
    batch <- pending[begin:end, , drop = FALSE]
    rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) run_one(batch[i, , drop = FALSE]),
                               mc.cores = min(cores, nrow(batch)), mc.preschedule = FALSE)
    existing <- rbind(existing, do.call(rbind, rows))
    utils::write.csv(existing, results_path, row.names = FALSE)
    elapsed <- proc.time()[["elapsed"]] - started_all
    message(sprintf("completed %d/%d pairs (this run %d/%d); elapsed %.1fs; ETA %.1fs",
                    nrow(existing), nrow(design), end, nrow(pending), elapsed,
                    elapsed / end * (nrow(pending) - end)))
  }
}

ok <- existing[existing$status == "ok", , drop = FALSE]
summary <- if (nrow(ok)) {
  data.frame(
    completed_pairs = nrow(ok), failed_pairs = sum(existing$status != "ok"),
    mean_abs_difference_ks = mean(abs(ok$difference_ks)),
    mean_abs_difference_cvm = mean(abs(ok$difference_cvm)),
    median_abs_difference_ks = stats::median(abs(ok$difference_ks)),
    median_abs_difference_cvm = stats::median(abs(ok$difference_cvm)),
    mean_elapsed_seconds = mean(ok$elapsed_seconds), stringsAsFactors = FALSE
  )
} else data.frame()
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
message("Results: ", normalizePath(results_path, mustWork = FALSE))
message("Summary: ", normalizePath(file.path(output_dir, "summary.csv"), mustWork = FALSE))
