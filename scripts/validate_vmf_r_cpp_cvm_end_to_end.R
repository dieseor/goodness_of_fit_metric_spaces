#!/usr/bin/env Rscript

# Paired end-to-end validation of the vMF fast-multiplier implementation.
# For every cell, the R and C++ routes use the identical data, fitted MLE,
# deterministic profile derivative and multiplier draws.  The only changed
# control is `fast_multiplier_backend`.

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
  values <- commandArgs(trailingOnly = TRUE)
  hit <- values[startsWith(values, prefix)]
  if (!length(hit)) default else substring(hit[[length(hit)]], nchar(prefix) + 1L)
}

integer_vector <- function(name, default) {
  value <- suppressWarnings(as.integer(strsplit(option(name, default), ",", fixed = TRUE)[[1L]]))
  if (!length(value) || any(!is.finite(value)) || any(value < 2L)) {
    stop(sprintf("`--%s` must be a nonempty vector of integers at least two.", name))
  }
  unique(value)
}

dimensions <- integer_vector("dimensions", "5,10")
sizes <- integer_vector("sizes", "50,200,400")
B <- as.integer(option("B", "299"))
cores <- as.integer(option("cores", "6"))
seed <- as.integer(option("seed", "20260824"))
output_dir <- option(
  "output-dir",
  "simulation_results/section6_new_scenarios/validation_vmf_r_cpp_cvm"
)
if (!is.finite(B) || B < 2L || !is.finite(cores) || cores < 1L || !is.finite(seed)) {
  stop("`B`, `cores`, and `seed` must be finite, with B >= 2 and cores >= 1.")
}

run_backend <- function(x, B, bootstrap_seed, backend) {
  multiplier_bootstrap_vmf(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    seed = bootstrap_seed,
    n_cores = 1L,
    distance_type = "geodesic",
    unknown_param = "xi",
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "quadrature",
      fast_multiplier_backend = backend,
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_bootstrap_chunk_size = 100L
    )
  )
}

maximum <- function(x) if (length(x)) max(abs(x)) else NA_real_

run_case <- function(item) {
  data_seed <- seed + 100000L * item$d + 1000L * item$n + item$index
  bootstrap_seed <- data_seed + 500L
  set.seed(data_seed)
  x <- normalize_vmf_data(rotasym::r_vMF(
    item$n,
    mu = c(1, rep.int(0, item$d)),
    kappa = item$d
  ))
  r_result <- run_backend(x, B, bootstrap_seed, "r")
  cpp_result <- run_backend(x, B, bootstrap_seed, "cpp")
  if (!identical(r_result$diagnostics$fast_multiplier_backend_effective, "r") ||
      !identical(cpp_result$diagnostics$fast_multiplier_backend_effective, "cpp") ||
      !identical(cpp_result$diagnostics$fast_multiplier_cpp_kernel_effective, "contiguous_double")) {
    stop("The intended R/C++ backend was not selected.")
  }
  data.frame(
    d = item$d,
    n = item$n,
    B = B,
    data_seed = data_seed,
    bootstrap_seed = bootstrap_seed,
    observed_ks_abs_difference = abs(r_result$inference$ks$observed - cpp_result$inference$ks$observed),
    observed_cvm_abs_difference = abs(r_result$inference$cvm$observed - cpp_result$inference$cvm$observed),
    bootstrap_ks_max_abs_difference = maximum(r_result$bootstrap$statistics$ks - cpp_result$bootstrap$statistics$ks),
    bootstrap_cvm_max_abs_difference = maximum(r_result$bootstrap$statistics$cvm - cpp_result$bootstrap$statistics$cvm),
    critical_ks_abs_difference = abs(r_result$inference$ks$critical_value - cpp_result$inference$ks$critical_value),
    critical_cvm_abs_difference = abs(r_result$inference$cvm$critical_value - cpp_result$inference$cvm$critical_value),
    p_value_ks_abs_difference = abs(r_result$inference$ks$p_value - cpp_result$inference$ks$p_value),
    p_value_cvm_abs_difference = abs(r_result$inference$cvm$p_value - cpp_result$inference$cvm$p_value),
    r_total_seconds = r_result$diagnostics$fast_total_seconds,
    cpp_total_seconds = cpp_result$diagnostics$fast_total_seconds,
    r_backend = r_result$diagnostics$fast_multiplier_backend_effective,
    cpp_backend = cpp_result$diagnostics$fast_multiplier_backend_effective,
    cpp_kernel = cpp_result$diagnostics$fast_multiplier_cpp_kernel_effective,
    stringsAsFactors = FALSE
  )
}

grid <- expand.grid(d = dimensions, n = sizes, KEEP.OUT.ATTRS = FALSE)
grid$index <- seq_len(nrow(grid))
items <- lapply(seq_len(nrow(grid)), function(i) grid[i, , drop = FALSE])
message(sprintf("Paired vMF R/C++ validation: %d cells (B=%d, cores=%d)", length(items), B, cores))
rows <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(items, run_case, mc.cores = min(cores, length(items)), mc.set.seed = FALSE)
} else {
  lapply(items, run_case)
}
summary <- do.call(rbind, rows)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
if (any(as.matrix(summary[, c(
  "observed_ks_abs_difference", "observed_cvm_abs_difference",
  "bootstrap_ks_max_abs_difference", "bootstrap_cvm_max_abs_difference",
  "critical_ks_abs_difference", "critical_cvm_abs_difference",
  "p_value_ks_abs_difference", "p_value_cvm_abs_difference"
)]) > 1e-10)) {
  stop("The paired R/C++ comparison exceeded the 1e-10 tolerance.")
}
message("R/C++ validation passed.")
