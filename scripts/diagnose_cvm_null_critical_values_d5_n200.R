#!/usr/bin/env Rscript

# Diagnostic for the finite-sample CvM calibration discrepancy between the
# two revised Normal scenarios under H0, using exactly d = 5 and n = 200 by
# default.  For each Monte Carlo sample it stores the observed CvM statistic
# and the conditional 0.95 multiplier-bootstrap critical value.  The final
# summary compares those critical values with the Monte Carlo 0.95 quantile
# of the observed null statistics.
#
# Scenario 1.1:
#   X ~ N_d(theta_d, I_d + 2 theta_d theta_d^T),
#   theta_d = d^{-1/2} 1_d,
# fitted with the restricted mean-aligned single-spike model.
#
# Scenario 1.2 under beta = 0:
#   X ~ N_d(0, I_d),
# fitted with N_d(mu, sigma^2 I_d).

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

args <- parse_args(commandArgs(trailingOnly = TRUE))
M <- as.integer(args$M %||% 1000L)
B <- as.integer(args$B %||% 5000L)
cores <- as.integer(args$cores %||% 10L)
d <- as.integer(args$d %||% 5L)
n <- as.integer(args$n %||% 200L)
derivative_mc_size <- as.integer(args$derivative_mc_size %||% 10000L)
cvm_block_size <- as.integer(args$cvm_block_size %||% 50L)
seed <- as.integer(args$seed %||% 20260826L)
output_dir <- as.character(args$output_dir %||% file.path(
  "simulation_results", "cvm_null_critical_value_diagnostic",
  sprintf("d%d_n%d_M%d_B%d_Nderiv%d", d, n, M, B, derivative_mc_size)
))

if (any(!is.finite(c(M, B, cores, d, n, derivative_mc_size, cvm_block_size, seed))) ||
    M < 1L || B < 2L || cores < 1L || d < 2L || n < 2L ||
    derivative_mc_size < 1L || cvm_block_size < 1L) {
  stop("Invalid diagnostic settings.")
}
if (.Platform$OS.type != "unix" && cores > 1L) {
  stop("Outer parallelism with cores > 1 requires a Unix-like platform.")
}

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
source(file.path("bootstrap", "normal_sigma_Id_bootstrap.R"))

write_atomic_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not update ", path)
}

seed_for <- function(base_seed, scenario_id, replication, stream) {
  as.integer((as.numeric(base_seed) +
    1000003 * scenario_id +
    1009 * replication +
    10000019 * stream) %% 2147483647) + 1L
}

scenario_names <- c("normal_1_1_single_spike", "normal_1_2_sigma_Id")
design <- expand.grid(
  scenario = scenario_names,
  replication = seq_len(M),
  stringsAsFactors = FALSE
)
design$scenario_id <- match(design$scenario, scenario_names)
design <- design[order(design$scenario_id, design$replication), , drop = FALSE]

empty_results <- function() {
  data.frame(
    scenario = character(), scenario_id = integer(), replication = integer(),
    d = integer(), n = integer(), status = character(), error_message = character(),
    seed_data = integer(), seed_bootstrap = integer(), seed_derivative = integer(),
    cvm_observed = numeric(), cvm_boot_q90 = numeric(), cvm_boot_q95 = numeric(),
    cvm_boot_q975 = numeric(), cvm_boot_q99 = numeric(), cvm_boot_mean = numeric(),
    cvm_boot_sd = numeric(), cvm_pvalue = numeric(), cvm_reject = logical(),
    effective_bootstrap_method = character(), fallback_to_reestimated = logical(),
    fast_backend = character(), fast_kernel = character(), derivative_method = character(),
    derivative_mc_size = integer(), vhat_method = character(), elapsed_seconds = numeric(),
    stringsAsFactors = FALSE
  )
}

result_key <- function(x) paste(x$scenario, x$replication, sep = "|")

run_one <- function(row) {
  scenario <- as.character(row$scenario)
  scenario_id <- as.integer(row$scenario_id)
  replication <- as.integer(row$replication)
  seed_data <- seed_for(seed, scenario_id, replication, 0L)
  seed_bootstrap <- seed_for(seed, scenario_id, replication, 1L)
  seed_derivative <- seed_for(seed, scenario_id, replication, 2L)
  began <- proc.time()[["elapsed"]]

  base <- data.frame(
    scenario = scenario,
    scenario_id = scenario_id,
    replication = replication,
    d = d,
    n = n,
    status = "ok",
    error_message = NA_character_,
    seed_data = seed_data,
    seed_bootstrap = seed_bootstrap,
    seed_derivative = seed_derivative,
    stringsAsFactors = FALSE
  )

  tryCatch({
    set.seed(seed_data)

    if (identical(scenario, "normal_1_1_single_spike")) {
      theta0 <- rep(1 / sqrt(d), d)
      sigma0 <- diag(d) + 2 * tcrossprod(theta0)
      x <- mvtnorm::rmvnorm(n, mean = theta0, sigma = sigma0)
      fit <- multiplier_bootstrap_restricted_spiked_normal(
        data = x,
        null = list(type = "composite"),
        statistics = "cvm",
        B = B,
        alpha = 0.05,
        n_cores = 1L,
        seed = seed_bootstrap,
        bootstrap_method = "fast_multiplier",
        keep = list(
          observed_process = FALSE,
          bootstrap_statistics = TRUE,
          bootstrap_thetas = FALSE
        ),
        control = list(
          derivative_method = "score_mc",
          derivative_mc_size = derivative_mc_size,
          derivative_mc_seed = seed_derivative,
          fast_multiplier_cvm_block_size = cvm_block_size,
          fast_multiplier_backend = "cpp",
          fast_multiplier_cpp_kernel = "contiguous_double",
          fast_multiplier_cache_corrections = "auto",
          fast_multiplier_stream_chunk_size = 100L,
          progress_bar = FALSE
        ),
        distance_profile_backend = "r",
        fast_multiplier_backend = "cpp",
        fast_multiplier_cpp_kernel = "contiguous_double",
        cache_block_corrections = "auto"
      )
    } else if (identical(scenario, "normal_1_2_sigma_Id")) {
      x <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
      fit <- multiplier_bootstrap_normal_sigma_Id(
        data = x,
        null = list(type = "composite"),
        statistics = "cvm",
        B = B,
        alpha = 0.05,
        n_cores = 1L,
        seed = seed_bootstrap,
        bootstrap_method = "fast_multiplier",
        keep = list(
          observed_process = FALSE,
          bootstrap_statistics = TRUE,
          bootstrap_thetas = FALSE
        ),
        control = list(
          derivative_method = "score_mc",
          derivative_mc_size = derivative_mc_size,
          derivative_mc_seed = seed_derivative,
          fast_multiplier_cvm_block_size = cvm_block_size,
          fast_multiplier_backend = "cpp",
          fast_multiplier_cpp_kernel = "contiguous_double",
          fast_multiplier_cache_corrections = "auto",
          fast_multiplier_stream_chunk_size = 100L,
          progress_bar = FALSE
        ),
        distance_profile_backend = "r",
        fast_multiplier_backend = "cpp",
        fast_multiplier_cpp_kernel = "contiguous_double",
        cache_block_corrections = "auto"
      )
    } else {
      stop("Unknown scenario: ", scenario)
    }

    cvm_star <- as.numeric(fit$bootstrap$statistics$cvm)
    if (length(cvm_star) != B || any(!is.finite(cvm_star))) {
      stop("Invalid stored CvM bootstrap statistics.")
    }

    q <- as.numeric(stats::quantile(
      cvm_star,
      probs = c(0.90, 0.95, 0.975, 0.99),
      names = FALSE,
      type = 8
    ))
    diagnostics <- fit$diagnostics

    out <- cbind(base, data.frame(
      cvm_observed = as.numeric(fit$inference$cvm$observed),
      cvm_boot_q90 = q[[1L]],
      cvm_boot_q95 = q[[2L]],
      cvm_boot_q975 = q[[3L]],
      cvm_boot_q99 = q[[4L]],
      cvm_boot_mean = mean(cvm_star),
      cvm_boot_sd = stats::sd(cvm_star),
      cvm_pvalue = as.numeric(fit$inference$cvm$p_value),
      cvm_reject = isTRUE(fit$inference$cvm$reject),
      effective_bootstrap_method = diagnostics$effective_bootstrap_method %||% NA_character_,
      fallback_to_reestimated = isTRUE(diagnostics$fallback_to_reestimated),
      fast_backend = diagnostics$fast_multiplier_backend_effective %||% NA_character_,
      fast_kernel = diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_,
      derivative_method = diagnostics$derivative_method_effective %||%
        diagnostics$derivative_method %||% NA_character_,
      derivative_mc_size = diagnostics$derivative_mc_size %||% NA_integer_,
      vhat_method = diagnostics$vhat_method %||% NA_character_,
      elapsed_seconds = proc.time()[["elapsed"]] - began,
      stringsAsFactors = FALSE
    ))

    conforming <- identical(out$effective_bootstrap_method[[1L]], "fast_multiplier") &&
      !isTRUE(out$fallback_to_reestimated[[1L]]) &&
      identical(out$fast_backend[[1L]], "cpp") &&
      identical(out$fast_kernel[[1L]], "contiguous_double") &&
      identical(out$derivative_method[[1L]], "score_mc") &&
      identical(as.integer(out$derivative_mc_size[[1L]]), derivative_mc_size)

    if (!conforming) {
      out$status <- "nonconforming"
      out$error_message <- "Requested fast C++ / score_mc configuration was not effective."
    }
    out
  }, error = function(e) {
    cbind(base, data.frame(
      cvm_observed = NA_real_, cvm_boot_q90 = NA_real_, cvm_boot_q95 = NA_real_,
      cvm_boot_q975 = NA_real_, cvm_boot_q99 = NA_real_, cvm_boot_mean = NA_real_,
      cvm_boot_sd = NA_real_, cvm_pvalue = NA_real_, cvm_reject = NA,
      effective_bootstrap_method = NA_character_, fallback_to_reestimated = NA,
      fast_backend = NA_character_, fast_kernel = NA_character_,
      derivative_method = NA_character_, derivative_mc_size = NA_integer_,
      vhat_method = NA_character_, elapsed_seconds = proc.time()[["elapsed"]] - began,
      stringsAsFactors = FALSE
    ), status = "error", error_message = conditionMessage(e))
  })
}

summarize_results <- function(results) {
  ok <- results[results$status == "ok", , drop = FALSE]
  if (!nrow(ok)) return(data.frame())

  rows <- lapply(split(ok, ok$scenario), function(z) {
    q_mc <- as.numeric(stats::quantile(z$cvm_observed, 0.95, names = FALSE, type = 8))
    gaps <- z$cvm_boot_q95 - q_mc
    data.frame(
      scenario = z$scenario[[1L]],
      d = z$d[[1L]],
      n = z$n[[1L]],
      M_ok = nrow(z),
      B = B,
      mc_q95_observed_cvm = q_mc,
      mean_boot_q95 = mean(z$cvm_boot_q95),
      median_boot_q95 = stats::median(z$cvm_boot_q95),
      boot_q95_q25 = as.numeric(stats::quantile(z$cvm_boot_q95, 0.25, names = FALSE, type = 8)),
      boot_q95_q75 = as.numeric(stats::quantile(z$cvm_boot_q95, 0.75, names = FALSE, type = 8)),
      mean_boot_minus_mc_q95 = mean(gaps),
      median_boot_minus_mc_q95 = stats::median(gaps),
      median_boot_to_mc_q95_ratio = stats::median(z$cvm_boot_q95) / q_mc,
      bootstrap_rejection_percent = 100 * mean(z$cvm_reject),
      rejection_at_mc_q95_percent = 100 * mean(z$cvm_observed > q_mc),
      mean_bootstrap_pvalue = mean(z$cvm_pvalue),
      mean_elapsed_seconds = mean(z$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(output_dir, "raw_results.csv")
summary_path <- file.path(output_dir, "summary.csv")
manifest_path <- file.path(output_dir, "manifest.txt")

writeLines(c(
  "CvM null critical-value diagnostic",
  sprintf("d=%d; n=%d; M=%d; B=%d; outer cores=%d; inner bootstrap cores=1", d, n, M, B, cores),
  sprintf("derivative=score_mc; N_deriv=%d; fast backend=C++ contiguous_double", derivative_mc_size),
  "Scenario 1.1 H0: N_d(theta_d, I_d + 2 theta_d theta_d^T), theta_d=d^{-1/2}1_d; restricted mean-aligned single spike",
  "Scenario 1.2 H0: N_d(0,I_d); fitted null N_d(mu,sigma^2 I_d)",
  "Per replication: observed CvM, bootstrap q90/q95/q97.5/q99, bootstrap mean/sd, p-value and reject",
  "Final comparison: Monte Carlo q95 of observed CvM versus distribution of conditional bootstrap q95 values",
  sprintf("base seed=%d", seed)
), manifest_path)

results <- if (file.exists(raw_path)) {
  utils::read.csv(raw_path, stringsAsFactors = FALSE)
} else {
  empty_results()
}
if (!all(names(empty_results()) %in% names(results))) {
  stop("Existing raw_results.csv has an incompatible schema; use a new output_dir.")
}

done <- if (nrow(results)) result_key(results[results$status == "ok", , drop = FALSE]) else character()
pending <- design[!result_key(design) %in% done, , drop = FALSE]

message(sprintf("CvM diagnostic: %d pending / %d total jobs; %d outer cores.",
                nrow(pending), nrow(design), cores))

while (nrow(pending)) {
  batch_n <- min(cores, nrow(pending))
  batch <- pending[seq_len(batch_n), , drop = FALSE]
  pending <- pending[-seq_len(batch_n), , drop = FALSE]

  rows <- parallel::mclapply(
    seq_len(nrow(batch)),
    function(i) run_one(batch[i, , drop = FALSE]),
    mc.cores = min(cores, nrow(batch)),
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  rows <- do.call(rbind, rows)

  new_keys <- result_key(rows)
  if (nrow(results)) {
    results <- results[!result_key(results) %in% new_keys, , drop = FALSE]
  }
  results <- rbind(results, rows)
  results <- results[order(results$scenario_id, results$replication), , drop = FALSE]
  write_atomic_csv(results, raw_path)
  write_atomic_csv(summarize_results(results), summary_path)

  message(sprintf(
    "completed %d/%d (ok=%d, errors=%d, nonconforming=%d)",
    nrow(results), nrow(design),
    sum(results$status == "ok"),
    sum(results$status == "error"),
    sum(results$status == "nonconforming")
  ))
}

summary <- summarize_results(results)
write_atomic_csv(summary, summary_path)
print(summary, row.names = FALSE)
message("Results: ", output_dir)
