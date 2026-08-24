#!/usr/bin/env Rscript

# Small, resumable-in-the-output sense pilot runner for the restricted-spiked
# Gaussian model.  It has no paper defaults: M, B, lambdas and n values must
# be supplied explicitly by the caller.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
arg_value <- function(name) {
  key <- paste0("--", name, "=")
  value <- args[startsWith(args, key)]
  if (!length(value)) return(NULL)
  sub(key, "", value[[1L]], fixed = TRUE)
}
required <- c("M", "B", "n", "lambdas", "cores", "output_dir")
missing <- required[vapply(required, function(x) is.null(arg_value(x)), logical(1L))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "))

M <- as.integer(arg_value("M")); B <- as.integer(arg_value("B")); cores <- as.integer(arg_value("cores"))
n_values <- as.integer(strsplit(arg_value("n"), ",", fixed = TRUE)[[1L]])
lambda_values <- as.numeric(strsplit(arg_value("lambdas"), ",", fixed = TRUE)[[1L]])
output_dir <- arg_value("output_dir")
seed <- as.integer(arg_value("seed") %||% "20260826")
if (any(!is.finite(c(M, B, cores, n_values, lambda_values, seed))) ||
    M < 1L || B < 1L || cores < 1L || any(n_values < 2L) || any(lambda_values <= 0)) {
  stop("Invalid pilot arguments.")
}

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

design <- expand.grid(d = c(2L, 5L), n = n_values, lambda = lambda_values,
                      replication = seq_len(M), stringsAsFactors = FALSE)
design$job_id <- seq_len(nrow(design))
result_path <- file.path(output_dir, "results.csv")
existing <- if (file.exists(result_path)) read.csv(result_path) else data.frame()
completed <- if (nrow(existing)) existing$job_id else integer()

run_one <- function(row) {
  set.seed(seed + 100003L * row$job_id)
  theta <- c(0.8, rep(0, row$d - 1L))
  x <- rrestricted_spiked_normal(row$n, theta, row$lambda)
  result <- multiplier_bootstrap_restricted_spiked_normal(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = cores,
    seed = seed + 100003L * row$job_id + 1L,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                bootstrap_thetas = FALSE),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 1000L,
      derivative_mc_seed = seed + 100003L * row$job_id + 2L,
      fast_multiplier_cvm_block_size = 50L,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto",
      progress_bar = FALSE
    ),
    distance_profile_backend = "r",
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fuse_ks_cvm = TRUE,
    cache_block_corrections = "auto"
  )
  data.frame(
    job_id = row$job_id, d = row$d, n = row$n, lambda = row$lambda,
    ks_pvalue = result$inference$ks$p_value,
    cvm_pvalue = result$inference$cvm$p_value,
    fast_backend = result$diagnostics$fast_multiplier_backend_effective,
    fast_kernel = result$diagnostics$fast_multiplier_cpp_kernel_effective,
    fast_fused = result$diagnostics$fast_multiplier_fuse_ks_cvm_effective,
    derivative_method = result$diagnostics$derivative_method_effective,
    vhat_method = result$diagnostics$vhat_method,
    elapsed_seconds = result$diagnostics$elapsed_seconds
  )
}

for (i in seq_len(nrow(design))) {
  if (design$job_id[[i]] %in% completed) next
  row_result <- run_one(design[i, ])
  existing <- rbind(existing, row_result)
  utils::write.csv(existing, result_path, row.names = FALSE)
  message(sprintf("completed %d/%d pilot jobs", nrow(existing), nrow(design)))
}

summary <- aggregate(cbind(ks_pvalue, cvm_pvalue) ~ d + n + lambda, existing, mean)
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
