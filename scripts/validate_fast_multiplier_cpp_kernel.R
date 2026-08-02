#!/usr/bin/env Rscript

# Paired validation and timing of the production cache-friendly C++ fast
# multiplier kernel against the explicit historical `legacy` fallback.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
source("bootstrap/multiplier_bootstrap.R")

option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  value <- commandArgs(trailingOnly = TRUE)
  value <- value[startsWith(value, prefix)]
  if (!length(value)) default else substring(value[[length(value)]], nchar(prefix) + 1L)
}

parse_integer_vector <- function(name, default) {
  value <- strsplit(option(name, default), ",", fixed = TRUE)[[1L]]
  value <- suppressWarnings(as.integer(value))
  if (!length(value) || any(!is.finite(value)) || any(value < 1L)) {
    stop(sprintf("`--%s` must be a nonempty vector of positive integers.", name))
  }
  unique(value)
}

sizes <- parse_integer_vector("sizes", "50,100,200,400")
dimensions <- parse_integer_vector("dimensions", "2,10")
B_chunk <- as.integer(option("B-chunk", "100"))
repetitions <- as.integer(option("reps", "5"))
timing_calls <- as.integer(option("timing-calls", "5"))
cores <- as.integer(option("cores", "1"))
seed <- as.integer(option("seed", "20260820"))
output <- option(
  "output",
  file.path("benchmarks", "fast_multiplier_cpp_kernel_validation.csv")
)
if (!is.finite(B_chunk) || B_chunk < 1L ||
    !is.finite(repetitions) || repetitions < 2L ||
    !is.finite(timing_calls) || timing_calls < 1L ||
    !is.finite(cores) || cores < 1L || !is.finite(seed)) {
  stop("Invalid scalar validation arguments.")
}

make_tie_end <- function(n, mode) {
  if (identical(mode, "none")) return(seq_len(n))
  groups <- integer()
  remaining <- n
  while (remaining > 0L) {
    size <- min(remaining, sample.int(5L, 1L))
    groups <- c(groups, size)
    remaining <- remaining - size
  }
  rep.int(cumsum(groups), groups)
}

median_elapsed <- function(fun, repetitions, timing_calls) {
  invisible(fun())
  timings <- numeric(repetitions)
  value <- NULL
  for (replicate_index in seq_len(repetitions)) {
    elapsed <- system.time({
      for (call_index in seq_len(timing_calls)) value <- fun()
    })[["elapsed"]]
    timings[[replicate_index]] <- elapsed / timing_calls
  }
  list(value = value, median_seconds = median(timings),
       minimum_seconds = min(timings))
}

run_case <- function(item) {
  n <- item$n
  q <- item$q
  p <- q + 1L
  n_centers <- n
  set.seed(seed + 100000L * q + 1000L * n + item$case_index)
  centered_weights <- matrix(
    stats::rexp(B_chunk * n) - 1,
    nrow = B_chunk,
    ncol = n
  )
  if (identical(item$weight_regime, "extreme")) {
    centered_weights[1L, ] <- seq(-1e8, 1e8, length.out = n)
    centered_weights[2L, ] <- rep(c(1e-8, -1e-8), length.out = n)
  }
  score_block <- matrix(stats::rnorm(B_chunk * p), B_chunk, p)
  obs_order_matrix <- t(replicate(n_centers, sample.int(n)))
  tie_end_row <- make_tie_end(n, item$tie_regime)
  tie_end_matrix <- matrix(
    rep.int(tie_end_row, n_centers), nrow = n_centers, byrow = TRUE
  )
  correction_matrix <- matrix(
    stats::rnorm(n_centers * n * p), nrow = n_centers * n, ncol = p
  )
  arguments <- list(
    centered_weights = centered_weights,
    score_block = score_block,
    obs_order_matrix = obs_order_matrix,
    tie_end_matrix = tie_end_matrix,
    correction_matrix = correction_matrix,
    scale_factor = 1,
    compute_ks = TRUE,
    compute_cvm = TRUE
  )
  legacy_call <- function() do.call(
    distance_profile_cpp_call,
    c(list("cpp_fast_sample_ks_cvm_stats"), arguments)
  )
  candidate_call <- function() do.call(
    distance_profile_cpp_call,
    c(list("cpp_fast_sample_ks_cvm_stats_contiguous_double"), arguments)
  )
  legacy <- median_elapsed(legacy_call, repetitions, timing_calls)
  candidate <- median_elapsed(candidate_call, repetitions, timing_calls)
  ks_difference <- max(abs(legacy$value$ks - candidate$value$ks))
  cvm_difference <- max(abs(legacy$value$cvm - candidate$value$cvm))
  data.frame(
    q = q,
    n = n,
    p = p,
    B_chunk = B_chunk,
    repetitions = repetitions,
    timing_calls = timing_calls,
    tie_regime = item$tie_regime,
    weight_regime = item$weight_regime,
    legacy_median_seconds = legacy$median_seconds,
    contiguous_double_median_seconds = candidate$median_seconds,
    speedup = legacy$median_seconds / candidate$median_seconds,
    legacy_minimum_seconds = legacy$minimum_seconds,
    contiguous_double_minimum_seconds = candidate$minimum_seconds,
    max_abs_ks_difference = ks_difference,
    max_abs_cvm_difference = cvm_difference,
    max_relative_ks_difference = ks_difference / max(1, max(abs(legacy$value$ks))),
    max_relative_cvm_difference = cvm_difference / max(1, max(abs(legacy$value$cvm))),
    correction_matrix_MiB = as.numeric(object.size(correction_matrix)) / 1024^2,
    stringsAsFactors = FALSE
  )
}

case_grid <- expand.grid(
  q = dimensions,
  n = sizes,
  tie_regime = c("none", "random"),
  weight_regime = c("ordinary", "extreme"),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
case_grid$case_index <- seq_len(nrow(case_grid))

ensure_distance_profile_cpp_loaded()
items <- lapply(seq_len(nrow(case_grid)), function(index) case_grid[index, ])
results <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    items,
    run_case,
    mc.cores = min(cores, length(items)),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  lapply(items, run_case)
}
results <- do.call(rbind, results)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output, row.names = FALSE)
print(results, row.names = FALSE)
cat(sprintf("Wrote %s\n", output))

if (any(results$max_relative_ks_difference > 1e-12) ||
    any(results$max_relative_cvm_difference > 1e-12)) {
  stop("The contiguous-double kernel exceeded relative agreement 1e-12.")
}
