#!/usr/bin/env Rscript

# Matched fast-versus-reestimated validation for the two Section 6 HvMF null
# settings in H^2 and H^10. Each row uses the same simulated data for both
# calibrations, so their rejection-rate difference is assessed with a paired
# Monte Carlo standard error. The full run is deliberately separate from the
# paper power campaign.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_section6_new_scenarios.R")

hvmf_validation_args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
hvmf_validation_csv <- function(name, default, mode = c("integer", "numeric")) {
  mode <- match.arg(mode)
  parse_section6_csv(hvmf_validation_args[[name]], default, mode)
}

hvmf_validation_empty_results <- function() {
  data.frame(
    scenario = character(), family = character(), alternative = character(),
    d = integer(), n = integer(), rep = integer(), method = character(),
    ks_pvalue = numeric(), cvm_pvalue = numeric(),
    ks_reject = logical(), cvm_reject = logical(),
    ks_observed = numeric(), cvm_observed = numeric(),
    effective_bootstrap_method = character(), fast_backend = character(),
    fast_fused = logical(), fallback_to_reestimated = logical(),
    elapsed_seconds = numeric(), status = character(), error_message = character(),
    stringsAsFactors = FALSE
  )
}

hvmf_validation_key <- function(x) {
  paste(x$scenario, x$d, x$n, x$rep, x$method, sep = "|")
}

hvmf_validation_pair_key <- function(x) {
  paste(x$scenario, x$d, x$n, x$rep, sep = "|")
}

hvmf_validation_one <- function(job, B, base_seed, derivative_mc_size, cvm_block_size) {
  started <- proc.time()[["elapsed"]]
  data_seed <- section6_seed(base_seed, job$design_id, job$rep, stream = 0L)
  x <- tryCatch({
    set.seed(data_seed)
    generate_section6_sample(job)
  }, error = function(e) e)
  methods <- c("fast_multiplier", "reestimated")
  rows <- lapply(seq_along(methods), function(index) {
    method <- methods[[index]]
    out <- data.frame(
      scenario = as.character(job$scenario), family = "hvmf",
      alternative = as.character(job$alternative), d = as.integer(job$d),
      n = as.integer(job$n), rep = as.integer(job$rep), method = method,
      ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA,
      cvm_reject = NA, ks_observed = NA_real_, cvm_observed = NA_real_,
      effective_bootstrap_method = NA_character_, fast_backend = NA_character_,
      fast_fused = NA, fallback_to_reestimated = NA,
      elapsed_seconds = NA_real_, status = "ok", error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    if (inherits(x, "error")) {
      out$status <- "error"
      out$error_message <- conditionMessage(x)
      return(out)
    }
    fit <- tryCatch(
      run_section6_bootstrap(
        design_row = job, x = x, B = B,
        # The two calibrations use the same simulated sample and multiplier
        # seed.  This removes avoidable Monte Carlo noise from their paired
        # comparison; the derivative seed is relevant only to the fast route.
        seed = section6_seed(base_seed, job$design_id, job$rep, stream = 1L),
        derivative_seed = section6_seed(base_seed, job$design_id, job$rep, stream = 11L),
        derivative_mc_size = derivative_mc_size, cvm_block_size = cvm_block_size,
        bootstrap_method = method
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      out$status <- "error"
      out$error_message <- conditionMessage(fit)
      return(out)
    }
    diagnostics <- fit$diagnostics
    out$ks_pvalue <- fit$inference$ks$p_value
    out$cvm_pvalue <- fit$inference$cvm$p_value
    out$ks_reject <- fit$inference$ks$reject
    out$cvm_reject <- fit$inference$cvm$reject
    out$ks_observed <- fit$inference$ks$observed
    out$cvm_observed <- fit$inference$cvm$observed
    out$effective_bootstrap_method <- diagnostics$effective_bootstrap_method %||% NA_character_
    out$fast_backend <- diagnostics$fast_multiplier_backend_effective %||% NA_character_
    out$fast_fused <- isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective)
    out$fallback_to_reestimated <- isTRUE(diagnostics$fallback_to_reestimated)
    if (identical(method, "fast_multiplier") &&
        (!identical(out$effective_bootstrap_method, "fast_multiplier") ||
         !identical(out$fast_backend, "cpp") || !isTRUE(out$fast_fused) ||
         isTRUE(out$fallback_to_reestimated))) {
      out$status <- "nonconforming"
      out$error_message <- "Requested fused C++ fast bootstrap was not effective."
    }
    out
  })
  output <- do.call(rbind, rows)
  output$elapsed_seconds <- proc.time()[["elapsed"]] - started
  output
}

hvmf_validation_summary <- function(results, alpha = 0.05) {
  ok <- results[results$status == "ok", , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  cell_key <- interaction(ok$scenario, ok$d, ok$n, drop = TRUE)
  cells <- lapply(split(ok, cell_key), function(cell) {
    fast <- cell[cell$method == "fast_multiplier", , drop = FALSE]
    slow <- cell[cell$method == "reestimated", , drop = FALSE]
    merge_key <- c("scenario", "d", "n", "rep")
    paired <- merge(
      fast[, c(merge_key, "ks_reject", "cvm_reject", "ks_observed", "cvm_observed")],
      slow[, c(merge_key, "ks_reject", "cvm_reject", "ks_observed", "cvm_observed")],
      by = merge_key, suffixes = c("_fast", "_reestimated")
    )
    summarize_statistic <- function(statistic) {
      delta <- paired[[paste0(statistic, "_reject_fast")]] -
        paired[[paste0(statistic, "_reject_reestimated")]]
      data.frame(
        statistic = statistic, M = nrow(paired),
        size_fast = mean(paired[[paste0(statistic, "_reject_fast")]]),
        size_reestimated = mean(paired[[paste0(statistic, "_reject_reestimated")]]),
        paired_difference = mean(delta),
        paired_mc_se = if (length(delta) > 1L) stats::sd(delta) / sqrt(length(delta)) else NA_real_,
        max_abs_observed_difference = max(abs(
          paired[[paste0(statistic, "_observed_fast")]] -
            paired[[paste0(statistic, "_observed_reestimated")]]
        )),
        stringsAsFactors = FALSE
      )
    }
    base <- cell[1L, c("scenario", "d", "n")]
    statistic_summary <- rbind(summarize_statistic("ks"), summarize_statistic("cvm"))
    cbind(base[rep.int(1L, nrow(statistic_summary)), , drop = FALSE], statistic_summary)
  })
  do.call(rbind, cells)
}

if (sys.nframe() == 0L) {
  M <- as.integer(hvmf_validation_args$M %||% 100L)
  B <- as.integer(hvmf_validation_args$B %||% 299L)
  cores <- as.integer(hvmf_validation_args$cores %||% 1L)
  if (!is.finite(M) || M < 1L || !is.finite(B) || B < 1L || !is.finite(cores) || cores < 1L) {
    stop("`M`, `B`, and `cores` must be strictly positive integers.")
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Parallel matched validation requires a Unix platform.")
  }

  output_dir <- as.character(hvmf_validation_args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios", "validation_hvmf_hq_fast_vs_reestimated_M100_B299"
  ))
  design <- make_section6_design(
    family = "hvmf", dimensions = hvmf_validation_csv("dimensions", c(2L, 10L), "integer"),
    n_values = hvmf_validation_csv("n_values", c(50L, 100L, 200L, 400L), "integer"), beta_values = 0
  )
  jobs <- merge(design, data.frame(rep = seq_len(M)), by = NULL)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else hvmf_validation_empty_results()
  # A completed matched pair is the atomic checkpoint. Rows from an interrupted
  # pair are recomputed, which prevents a one-method-only replication from
  # entering the paired comparison.
  completed_pair_key <- character()
  if (nrow(existing)) {
    ok <- existing[existing$status == "ok", , drop = FALSE]
    by_pair <- split(ok$method, hvmf_validation_pair_key(ok))
    completed_pair_key <- names(by_pair)[vapply(
      by_pair,
      function(methods) setequal(methods, c("fast_multiplier", "reestimated")),
      logical(1)
    )]
  }
  pending <- jobs[!hvmf_validation_pair_key(jobs) %in% completed_pair_key, , drop = FALSE]

  started <- Sys.time()
  cat(sprintf(
    "HvMF H^q matched validation: %d pending pairs (M=%d, B=%d, cores=%d)\n",
    nrow(pending), M, B, cores
  ))
  base_seed <- as.integer(hvmf_validation_args$seed %||% 20260730L)
  derivative_mc_size <- as.integer(hvmf_validation_args$derivative_mc_size %||% 1000L)
  cvm_block_size <- as.integer(hvmf_validation_args$cvm_block_size %||% 50L)
  completed_now <- 0L
  batch_starts <- seq.int(1L, nrow(pending), by = cores)
  for (first in batch_starts) {
    last <- min(first + cores - 1L, nrow(pending))
    batch_indices <- seq.int(first, last)
    evaluate_pair <- function(index) {
      hvmf_validation_one(
        pending[index, , drop = FALSE], B = B,
        base_seed = base_seed, derivative_mc_size = derivative_mc_size,
        cvm_block_size = cvm_block_size
      )
    }
    batch_rows <- if (length(batch_indices) == 1L) {
      list(evaluate_pair(batch_indices))
    } else {
      parallel::mclapply(
        batch_indices, evaluate_pair,
        mc.cores = min(cores, length(batch_indices)), mc.preschedule = FALSE
      )
    }
    existing <- rbind(existing, do.call(rbind, batch_rows))
    section6_write_atomic_csv(existing, result_path)
    section6_write_atomic_csv(hvmf_validation_summary(existing), summary_path)
    completed_now <- completed_now + length(batch_indices)
    cat(sprintf("\rcompleted pairs: %d/%d", completed_now, nrow(pending)))
    flush.console()
  }
  if (nrow(pending)) cat("\n")
  summary <- hvmf_validation_summary(existing)
  section6_write_atomic_csv(summary, summary_path)
  writeLines(c(
    sprintf("M: %d", M), sprintf("B: %d", B), sprintf("outer_cores: %d", cores),
    "families: hvmf", "dimensions: 2,10", "n: 50,100,200,400",
    "scenarios: hvmf_1_mixture,hvmf_2_angular", "beta: 0",
    "KS: sample_unique_distances", "fast: fused C++ multiplier bootstrap",
    "comparison: paired fast-versus-reestimated on identical simulated datasets"
  ), file.path(output_dir, "manifest.txt"))
  print(summary)
}
