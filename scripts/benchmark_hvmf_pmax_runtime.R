#!/usr/bin/env Rscript

# Runtime-only diagnostic for the HvMF polar sampler regularisation.
# It does not change the sampler defaults or write simulation results.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

source("scripts/run_section6_new_scenarios.R", encoding = "UTF-8")

benchmark_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  values <- values[startsWith(values, prefix)]
  if (!length(values)) return(default)
  substring(values[[length(values)]], nchar(prefix) + 1L)
}

benchmark_csv <- function(name, default, storage.mode = c("integer", "numeric")) {
  storage.mode <- match.arg(storage.mode)
  raw <- benchmark_arg(name, NULL)
  if (is.null(raw) || !nzchar(raw)) return(default)
  values <- trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])
  switch(storage.mode, integer = as.integer(values), numeric = as.numeric(values))
}

benchmark_hvmf_design <- function(dimensions, n_values, p_max_values, repetitions) {
  base <- expand.grid(
    scenario = c("hvmf_1_mixture", "hvmf_2_angular"),
    d = sort(unique(as.integer(dimensions))),
    n = sort(unique(as.integer(n_values))),
    p_max = sort(unique(as.numeric(p_max_values))),
    repetition = seq_len(as.integer(repetitions)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(base$d)) || any(base$d < 2L) ||
      any(!is.finite(base$n)) || any(base$n < 1L) ||
      any(!is.finite(base$p_max)) || any(base$p_max <= 0 | base$p_max >= 1)) {
    stop("Invalid dimension, sample-size, or `p_max` grid.", call. = FALSE)
  }
  if (length(repetitions) != 1L || !is.finite(repetitions) || repetitions < 1L ||
      repetitions != as.integer(repetitions)) {
    stop("`repetitions` must be a positive integer.", call. = FALSE)
  }
  ## Pair p_max values through identical random-number streams.  The seed key
  ## deliberately excludes p_max, whereas design_id records the run itself.
  base$seed_id <- match(
    interaction(base$scenario, base$d, base$n, base$repetition, drop = TRUE),
    unique(interaction(base$scenario, base$d, base$n, base$repetition, drop = TRUE))
  )
  ## Randomising the execution order avoids systematically charging one p_max
  ## value for residual cache or worker-start costs.
  base <- base[sample.int(nrow(base)), , drop = FALSE]
  base$design_id <- seq_len(nrow(base))
  base
}

benchmark_hvmf_null_sample <- function(scenario, d, n, p_max, probability_step) {
  mu0 <- c(sqrt(2), section6_e(d))
  kappa <- switch(
    scenario,
    hvmf_1_mixture = as.numeric(d),
    hvmf_2_angular = 1.5 * as.numeric(d),
    stop("Unsupported scenario.", call. = FALSE)
  )
  rhvmf_polar(
    n = n, mu = mu0, kappa = kappa,
    p_max = p_max, probability_step = probability_step
  )
}

benchmark_hvmf_pmax_job <- function(job, B, seed, probability_step, cvm_block_size) {
  data_seed <- section6_seed(seed, job$seed_id, 1L, 0L)
  bootstrap_seed <- section6_seed(seed, job$seed_id, 1L, 1L)
  output <- data.frame(
    scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
    repetition = as.integer(job$repetition),
    p_max = as.numeric(job$p_max), probability_step = as.numeric(probability_step),
    seed_id = as.integer(job$seed_id),
    B = as.integer(B), status = "ok", error_message = NA_character_,
    generation_seconds = NA_real_, bootstrap_seconds = NA_real_, total_seconds = NA_real_,
    ks_pvalue = NA_real_, cvm_pvalue = NA_real_, stringsAsFactors = FALSE
  )
  tryCatch({
    set.seed(data_seed)
    generated_at <- proc.time()[["elapsed"]]
    x <- benchmark_hvmf_null_sample(
      scenario = job$scenario, d = job$d, n = job$n, p_max = job$p_max,
      probability_step = probability_step
    )
    output$generation_seconds <- proc.time()[["elapsed"]] - generated_at

    row <- data.frame(
      scenario = job$scenario, family = "hvmf", alternative = "null_runtime_benchmark",
      d = job$d, n = job$n, beta = 0, design_id = job$design_id,
      stringsAsFactors = FALSE
    )
    bootstrapped_at <- proc.time()[["elapsed"]]
    fit <- run_section6_bootstrap(
      design_row = row, x = x, B = B, seed = bootstrap_seed,
      derivative_seed = NA_integer_, derivative_mc_size = 10000L,
      cvm_block_size = cvm_block_size, derivative_method = "quadrature", n_cores = 1L
    )
    output$bootstrap_seconds <- proc.time()[["elapsed"]] - bootstrapped_at
    output$total_seconds <- output$generation_seconds + output$bootstrap_seconds
    output$ks_pvalue <- fit$inference$ks$p_value
    output$cvm_pvalue <- fit$inference$cvm$p_value
    output
  }, error = function(error) {
    output$status <- "error"
    output$error_message <- conditionMessage(error)
    output
  })
}

run_hvmf_pmax_runtime_benchmark <- function(output_dir,
                                             dimensions = c(2L, 5L),
                                             n_values = c(50L, 400L, 800L),
                                             p_max_values = c(0.99, 0.999),
                                             probability_step = 0.01,
                                             B = 1000L,
                                             repetitions = 3L,
                                             cores = 8L,
                                             seed = 20260906L,
                                             cvm_block_size = 50L) {
  if (length(B) != 1L || !is.finite(B) || B < 1L || B != as.integer(B)) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }
  if (length(cores) != 1L || !is.finite(cores) || cores < 1L) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  if (length(probability_step) != 1L || !is.finite(probability_step) || probability_step <= 0) {
    stop("`probability_step` must be strictly positive.", call. = FALSE)
  }

  design <- benchmark_hvmf_design(dimensions, n_values, p_max_values, repetitions)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(design, file.path(output_dir, "manifest.csv"), row.names = FALSE)

  ## Initialise the C++/quadrature path once before timing forked workers.
  ## This warm-up is deliberately excluded from timings.csv.
  warmup_job <- data.frame(
    scenario = "hvmf_1_mixture", d = min(dimensions), n = 10L,
    p_max = min(p_max_values), repetition = 0L, seed_id = 0L, design_id = 0L,
    stringsAsFactors = FALSE
  )
  invisible(benchmark_hvmf_pmax_job(
    warmup_job, B = min(3L, as.integer(B)), seed = seed,
    probability_step = probability_step, cvm_block_size = cvm_block_size
  ))

  started <- Sys.time()
  rows <- parallel::mclapply(
    seq_len(nrow(design)),
    function(i) benchmark_hvmf_pmax_job(
      design[i, , drop = FALSE], B = B, seed = seed,
      probability_step = probability_step, cvm_block_size = cvm_block_size
    ),
    mc.cores = min(as.integer(cores), nrow(design)), mc.preschedule = FALSE
  )
  results <- do.call(rbind, rows)
  utils::write.csv(results, file.path(output_dir, "timings.csv"), row.names = FALSE)

  successful <- results[results$status == "ok", , drop = FALSE]
  summary_by_pmax <- aggregate(
    cbind(generation_seconds, bootstrap_seconds, total_seconds) ~
      scenario + d + n + p_max + B + probability_step,
    data = successful,
    FUN = median
  )
  baseline <- summary_by_pmax[summary_by_pmax$p_max == min(p_max_values), , drop = FALSE]
  candidate <- summary_by_pmax[summary_by_pmax$p_max == max(p_max_values), , drop = FALSE]
  key <- c("scenario", "d", "n", "B", "probability_step")
  comparison <- merge(
    baseline[, c(key, "generation_seconds", "bootstrap_seconds", "total_seconds")],
    candidate[, c(key, "generation_seconds", "bootstrap_seconds", "total_seconds")],
    by = key, suffixes = c("_p99", "_p999"), all = TRUE
  )
  comparison$generation_ratio <- comparison$generation_seconds_p999 / comparison$generation_seconds_p99
  comparison$bootstrap_ratio <- comparison$bootstrap_seconds_p999 / comparison$bootstrap_seconds_p99
  comparison$total_ratio <- comparison$total_seconds_p999 / comparison$total_seconds_p99
  utils::write.csv(comparison, file.path(output_dir, "comparison.csv"), row.names = FALSE)
  writeLines(c(
    sprintf("started: %s", format(started, tz = "Europe/Madrid")),
    sprintf("finished: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("elapsed_seconds: %.3f", as.numeric(difftime(Sys.time(), started, units = "secs"))),
    "purpose: runtime comparison only; no sampler default was changed",
    sprintf("p_max_values: %s", paste(p_max_values, collapse = ",")),
    sprintf("probability_step: %.17g", probability_step),
    sprintf("B: %d; repetitions: %d; outer_cores: %d; inner_cores: 1",
            B, as.integer(repetitions), as.integer(cores))
  ), file.path(output_dir, "README.txt"))

  print(comparison)
  invisible(list(results = results, comparison = comparison))
}

if (sys.nframe() == 0L) {
  run_hvmf_pmax_runtime_benchmark(
    output_dir = benchmark_arg(
      "output_dir",
      "simulation_results/section6_new_scenarios/benchmark_hvmf_pmax_099_vs_0999"
    ),
    dimensions = benchmark_csv("dimensions", c(2L, 5L), "integer"),
    n_values = benchmark_csv("n_values", c(50L, 400L, 800L), "integer"),
    p_max_values = benchmark_csv("p_max_values", c(0.99, 0.999), "numeric"),
    probability_step = as.numeric(benchmark_arg("probability_step", "0.01")),
    B = as.integer(benchmark_arg("B", "1000")),
    repetitions = as.integer(benchmark_arg("repetitions", "3")),
    cores = as.integer(benchmark_arg("cores", "8")),
    seed = as.integer(benchmark_arg("seed", "20260906")),
    cvm_block_size = as.integer(benchmark_arg("cvm_block_size", "50"))
  )
}
