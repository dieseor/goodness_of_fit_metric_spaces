#!/usr/bin/env Rscript

# Matched fast-versus-reestimated null-calibration check for one Section 6
# family.  Every completed unit contains both methods for the same simulated
# sample and multiplier seed, which makes the comparison paired and resumable.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_section6_new_scenarios.R")

validation_args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
validation_csv <- function(name, default, mode = c("integer", "numeric")) {
  mode <- match.arg(mode)
  parse_section6_csv(validation_args[[name]], default, mode)
}

empty_validation_results <- function() {
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

validation_pair_key <- function(x) {
  paste(x$scenario, x$d, x$n, x$rep, sep = "|")
}

run_matched_pair <- function(job, B, base_seed, derivative_mc_size,
                             cvm_block_size, bootstrap_cores = 1L) {
  started <- proc.time()[["elapsed"]]
  data_seed <- section6_seed(base_seed, job$design_id, job$rep, stream = 0L)
  x <- tryCatch({
    set.seed(data_seed)
    generate_section6_sample(job)
  }, error = function(e) e)

  rows <- lapply(c("fast_multiplier", "reestimated"), function(method) {
    out <- data.frame(
      scenario = as.character(job$scenario), family = as.character(job$family),
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
        # The multiplier seed is shared.  The auxiliary seed is used only by
        # the fast calculation.
        seed = section6_seed(base_seed, job$design_id, job$rep, stream = 1L),
        derivative_seed = section6_seed(base_seed, job$design_id, job$rep, stream = 11L),
        derivative_mc_size = derivative_mc_size, cvm_block_size = cvm_block_size,
        bootstrap_method = method, n_cores = bootstrap_cores
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
  result <- do.call(rbind, rows)
  result$elapsed_seconds <- proc.time()[["elapsed"]] - started
  result
}

summarize_matched_validation <- function(results) {
  ok <- results[results$status == "ok", , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  groups <- split(ok, interaction(ok$scenario, ok$d, ok$n, drop = TRUE))
  summaries <- lapply(groups, function(cell) {
    merge_key <- c("scenario", "family", "alternative", "d", "n", "rep")
    fast <- cell[cell$method == "fast_multiplier", , drop = FALSE]
    slow <- cell[cell$method == "reestimated", , drop = FALSE]
    paired <- merge(
      fast[, c(merge_key, "ks_pvalue", "cvm_pvalue", "ks_reject", "cvm_reject", "ks_observed", "cvm_observed")],
      slow[, c(merge_key, "ks_pvalue", "cvm_pvalue", "ks_reject", "cvm_reject", "ks_observed", "cvm_observed")],
      by = merge_key, suffixes = c("_fast", "_reestimated")
    )
    one_stat <- function(statistic) {
      fast_reject <- paired[[paste0(statistic, "_reject_fast")]]
      slow_reject <- paired[[paste0(statistic, "_reject_reestimated")]]
      delta <- fast_reject - slow_reject
      pvalue_delta <- paired[[paste0(statistic, "_pvalue_fast")]] -
        paired[[paste0(statistic, "_pvalue_reestimated")]]
      data.frame(
        statistic = statistic, M = nrow(paired),
        size_fast = mean(fast_reject), size_reestimated = mean(slow_reject),
        paired_difference = mean(delta),
        paired_mc_se = if (nrow(paired) > 1L) stats::sd(delta) / sqrt(nrow(paired)) else NA_real_,
        mean_pvalue_fast_minus_reestimated = mean(pvalue_delta),
        median_pvalue_fast_minus_reestimated = stats::median(pvalue_delta),
        pvalue_difference_mc_se = if (nrow(paired) > 1L) stats::sd(pvalue_delta) / sqrt(nrow(paired)) else NA_real_,
        min_pvalue_fast_minus_reestimated = min(pvalue_delta),
        max_pvalue_fast_minus_reestimated = max(pvalue_delta),
        n_pvalue_fast_larger = sum(pvalue_delta > 0),
        max_abs_observed_difference = max(abs(
          paired[[paste0(statistic, "_observed_fast")]] -
            paired[[paste0(statistic, "_observed_reestimated")]]
        )),
        stringsAsFactors = FALSE
      )
    }
    base <- cell[1L, c("scenario", "family", "alternative", "d", "n"), drop = FALSE]
    rownames(base) <- NULL
    statistic_summary <- rbind(one_stat("ks"), one_stat("cvm"))
    rownames(statistic_summary) <- NULL
    cbind(base[rep.int(1L, nrow(statistic_summary)), , drop = FALSE], statistic_summary)
  })
  do.call(rbind, summaries)
}

if (sys.nframe() == 0L) {
  family <- tolower(as.character(validation_args$family %||% "vmf"))
  if (!family %in% section6_families) {
    stop("`family` must be one of normal, lg, vmf, or hvmf.")
  }
  M <- as.integer(validation_args$M %||% 30L)
  B <- as.integer(validation_args$B %||% 199L)
  cores <- as.integer(validation_args$cores %||% 1L)
  outer_cores <- as.integer(validation_args$outer_cores %||% cores)
  bootstrap_cores <- as.integer(validation_args$bootstrap_cores %||% 1L)
  if (any(!is.finite(c(M, B, cores, outer_cores, bootstrap_cores))) ||
      any(c(M, B, cores, outer_cores, bootstrap_cores) < 1L)) {
    stop("`M`, `B`, `cores`, `outer_cores`, and `bootstrap_cores` must be strictly positive integers.")
  }
  if (outer_cores * bootstrap_cores > cores) {
    stop("`outer_cores * bootstrap_cores` cannot exceed `cores`.")
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Parallel matched validation requires a Unix platform.")
  }

  dimensions <- validation_csv("dimensions", 10L, "integer")
  n_values <- validation_csv("n_values", 200L, "integer")
  output_dir <- as.character(validation_args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios",
    sprintf("validation_%s_d%s_n%s_M%d_B%d_fast_vs_reestimated",
      family, paste(dimensions, collapse = "_"), paste(n_values, collapse = "_"), M, B)
  ))
  derivative_mc_size <- as.integer(validation_args$derivative_mc_size %||% 10000L)
  cvm_block_size <- as.integer(validation_args$cvm_block_size %||% 50L)
  base_seed <- as.integer(validation_args$seed %||% 20260802L)

  design <- make_section6_design(
    family = family, dimensions = dimensions, n_values = n_values, beta_values = 0
  )
  jobs <- merge(design, data.frame(rep = seq_len(M)), by = NULL)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  existing <- if (file.exists(result_path)) {
    utils::read.csv(result_path, stringsAsFactors = FALSE)
  } else {
    empty_validation_results()
  }
  completed <- character()
  if (nrow(existing)) {
    ok <- existing[existing$status == "ok", , drop = FALSE]
    methods_by_pair <- split(ok$method, validation_pair_key(ok))
    completed <- names(methods_by_pair)[vapply(
      methods_by_pair,
      function(x) setequal(x, c("fast_multiplier", "reestimated")), logical(1)
    )]
  }
  pending <- jobs[!validation_pair_key(jobs) %in% completed, , drop = FALSE]
  cat(sprintf(
    "%s matched fast-versus-reestimated validation: %d pending pairs (M=%d, B=%d, cores=%d)\n",
    family, nrow(pending), M, B, cores
  ))
  completed_now <- 0L
  for (first in seq.int(1L, nrow(pending), by = outer_cores)) {
    indices <- seq.int(first, min(first + outer_cores - 1L, nrow(pending)))
    evaluate <- function(index) run_matched_pair(
      pending[index, , drop = FALSE], B = B, base_seed = base_seed,
      derivative_mc_size = derivative_mc_size, cvm_block_size = cvm_block_size,
      bootstrap_cores = bootstrap_cores
    )
    rows <- if (length(indices) == 1L) list(evaluate(indices)) else parallel::mclapply(
      indices, evaluate, mc.cores = min(outer_cores, length(indices)), mc.preschedule = FALSE
    )
    existing <- rbind(existing, do.call(rbind, rows))
    section6_write_atomic_csv(existing, result_path)
    section6_write_atomic_csv(summarize_matched_validation(existing), summary_path)
    completed_now <- completed_now + length(indices)
    cat(sprintf("\rcompleted pairs: %d/%d", completed_now, nrow(pending)))
    flush.console()
  }
  if (nrow(pending)) cat("\n")
  writeLines(c(
    sprintf("family: %s", family), sprintf("M: %d", M), sprintf("B: %d", B),
    sprintf("total_core_budget: %d", cores), sprintf("outer_cores: %d", outer_cores),
    sprintf("bootstrap_cores: %d", bootstrap_cores),
    sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
    sprintf("n: %s", paste(n_values, collapse = ",")), "beta: 0",
    "KS: sample_unique_distances", "fast: requested fused C++ multiplier bootstrap",
    "comparison: paired fast-versus-reestimated on identical simulated samples and multipliers"
  ), file.path(output_dir, "manifest.txt"))
  print(summarize_matched_validation(existing))
}
