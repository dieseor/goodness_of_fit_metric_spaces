#!/usr/bin/env Rscript

# Simple-null control for the Section 6 Normal d=10 experiment.  No parameter
# is fitted in either the observed statistic or the multiplier replicates.
# `bootstrap_method = reestimated` is merely the shared engine name: under a
# simple null it performs no re-estimation.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_section6_new_scenarios.R")

simple_theta0_section6 <- function(job) {
  d <- as.integer(job$d)
  if (identical(as.character(job$scenario), "normal_1_mixture")) {
    return(list(mu = 0.5 * section6_e(d), Sigma = section6_sigma(d, "plus")))
  }
  if (identical(as.character(job$scenario), "normal_2_t3")) {
    return(list(mu = rep(0, d), Sigma = diag(d)))
  }
  stop("Unsupported Normal scenario in simple-null diagnostic.")
}

run_one_simple <- function(job, rep, B, base_seed, cvm_block_size) {
  started <- proc.time()[["elapsed"]]
  data_seed <- section6_seed(base_seed, job$design_id, rep, 0L)
  bootstrap_seed <- section6_seed(base_seed, job$design_id, rep, 1L)
  tryCatch({
    set.seed(data_seed)
    x <- generate_section6_sample(job)
    spec <- make_mvnormal_spec("both")
    fit <- multiplier_bootstrap_gof(
      data = x, spec = spec,
      null = list(type = "simple", theta = simple_theta0_section6(job)),
      statistics = c("ks", "cvm"), ks_grid = make_sample_unique_distance_ks_grid(),
      B = as.integer(B), alpha = .05, n_cores = 1L, seed = as.integer(bootstrap_seed),
      bootstrap_method = "reestimated",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                  bootstrap_thetas = FALSE),
      control = list(
        cvm_block_size = as.integer(cvm_block_size),
        fast_multiplier_backend = "cpp"
      )
    )
    data.frame(
      scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
      rep = as.integer(rep), status = "ok", error_message = NA_character_,
      ks_pvalue = fit$inference$ks$p_value, cvm_pvalue = fit$inference$cvm$p_value,
      ks_reject = fit$inference$ks$reject, cvm_reject = fit$inference$cvm$reject,
      method = fit$diagnostics$effective_bootstrap_method,
      seed_data = data_seed, seed_bootstrap = bootstrap_seed,
      elapsed_seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE
    )
  }, error = function(e) data.frame(
    scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
    rep = as.integer(rep), status = "error", error_message = conditionMessage(e),
    ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
    method = NA_character_, seed_data = data_seed, seed_bootstrap = bootstrap_seed,
    elapsed_seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE
  ))
}

summarize_simple <- function(x) {
  ok <- x[x$status == "ok", , drop = FALSE]
  do.call(rbind, lapply(split(ok, ok$scenario), function(z) data.frame(
    scenario = z$scenario[[1L]], M = nrow(z),
    ks_size = mean(z$ks_reject), cvm_size = mean(z$cvm_reject),
    median_ks_pvalue = stats::median(z$ks_pvalue),
    median_cvm_pvalue = stats::median(z$cvm_pvalue),
    stringsAsFactors = FALSE
  )))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  M <- as.integer(args$M %||% 100L)
  B <- as.integer(args$B %||% 299L)
  cores <- as.integer(args$cores %||% 3L)
  cvm_block_size <- as.integer(args$cvm_block_size %||% 50L)
  base_seed <- as.integer(args$seed %||% 20260809L)
  output_dir <- as.character(args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios",
    "diagnostic_normal_simple_d10_n100_M100_B299"
  ))
  if (M < 2L || B < 99L || cores < 1L) stop("Invalid controls.")
  ensure_distance_profile_cpp_loaded()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(lock), add = TRUE)
  result_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  design <- make_section6_design("normal", dimensions = 10L, n_values = 100L, beta_values = 0)
  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else data.frame()
  tasks <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) cbind(
    design[rep(i, M), , drop = FALSE], rep = seq_len(M)
  )))
  done <- if (nrow(existing)) with(existing[existing$status == "ok", ], paste(scenario, rep, sep = "|")) else character()
  pending <- tasks[!paste(tasks$scenario, tasks$rep, sep = "|") %in% done, , drop = FALSE]
  writeLines(c(
    "purpose: simple-null control; not a composite-null simulation",
    "no parameter fitting or derivative/V correction is used",
    "engine label reestimated means no re-estimation under a simple null",
    sprintf("M per scenario: %d; B: %d; outer cores: %d", M, B, cores)
  ), file.path(output_dir, "README.txt"))
  total <- nrow(tasks)
  message(sprintf("Simple-null diagnostic: %d/%d completed; %d pending.", total - nrow(pending), total, nrow(pending)))
  while (nrow(pending)) {
    batch <- pending[seq_len(min(cores, nrow(pending))), , drop = FALSE]
    rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) {
      run_one_simple(batch[i, , drop = FALSE], batch$rep[[i]], B, base_seed, cvm_block_size)
    }, mc.cores = min(cores, nrow(batch)), mc.set.seed = FALSE)
    existing <- rbind(existing, do.call(rbind, rows))
    section6_write_atomic_csv(existing, result_path)
    section6_write_atomic_csv(summarize_simple(existing), summary_path)
    pending <- pending[-seq_len(nrow(batch)), , drop = FALSE]
    message(sprintf("completed replications: %d/%d", total - nrow(pending), total))
  }
}
