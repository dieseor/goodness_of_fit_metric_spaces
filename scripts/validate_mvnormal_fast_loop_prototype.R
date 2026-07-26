#!/usr/bin/env Rscript

source(file.path("bootstrap", "multiplier_bootstrap.R"))

ensure_distance_profile_cpp_loaded()

args <- commandArgs(trailingOnly = TRUE)
cases_arg <- grep("^--cases=", args, value = TRUE)
n_cases <- if (length(cases_arg)) {
  as.integer(sub("^--cases=", "", cases_arg[[1L]]))
} else {
  500L
}
if (!is.finite(n_cases) || n_cases <= 0L) stop("`--cases` must be positive.")

set.seed(20260723L)

for (case_id in seq_len(n_cases)) {
  d <- sample(c(2L, 10L), 1L)
  p <- d + d * (d + 1L) / 2L
  n <- sample(2:50, 1L)
  B <- sample.int(30L, 1L)
  scale_factor <- runif(1L, min = 0.25, max = 2.5)

  raw_weights <- matrix(rexp(B * n), nrow = B, ncol = n)
  centered_weights <- raw_weights / rowMeans(raw_weights) - 1
  S_obs <- matrix(rnorm(n * p), nrow = n, ncol = p)
  obs_order_matrix <- t(vapply(seq_len(n), function(i) sample.int(n), integer(n)))

  # Rounded values deliberately introduce ties, including complete ties in a
  # subset of cases.  The order matrix remains a valid permutation because the
  # production cache stores ordering and sorted distances separately.
  obs_sorted_distance_matrix <- t(vapply(seq_len(n), function(i) {
    if (case_id %% 25L == 0L) {
      rep.int(0, n)
    } else {
      sort(round(runif(n), digits = sample(0:3, 1L)))
    }
  }, numeric(n)))
  correction_values <- matrix(rnorm(n * n * p), nrow = n * n, ncol = p)
  correction_cache <- list(values = correction_values, bytes = as.double(length(correction_values)) * 8)
  tie_end_matrix <- build_sorted_tie_end_matrix(
    obs_sorted_distance_matrix
  )

  ks_stream <- list(
    mode = "sample_points_unique_distances_streamed",
    S_obs = S_obs,
    obs_order_matrix = obs_order_matrix,
    obs_sorted_distance_matrix = obs_sorted_distance_matrix,
    tie_end_matrix = tie_end_matrix,
    n = n,
    correction_cache = correction_cache
  )
  cvm_stream <- list(
    mode = "sample_points_unique_distances_sorted_rows",
    S_obs = S_obs,
    obs_order_matrix = obs_order_matrix,
    obs_sorted_distance_matrix = obs_sorted_distance_matrix,
    tie_end_matrix = tie_end_matrix,
    n = n,
    block_size = sample.int(n, 1L),
    correction_cache = correction_cache
  )

  expected_ks <- compute_fast_ks_sample_stats_streamed(
    centered_weight_block = centered_weights,
    ks_sample_stream_prep = ks_stream,
    scale_factor = scale_factor
  )
  expected_cvm <- compute_fast_cvm_stats_streamed(
    centered_weight_block = centered_weights,
    cvm_stream_prep = cvm_stream,
    scale_factor = scale_factor
  )
  observed <- compute_fast_sample_ks_cvm_stats_cpp(
    centered_weight_block = centered_weights,
    stream_prep = ks_stream,
    scale_factor = scale_factor,
    compute_ks = TRUE,
    compute_cvm = TRUE,
    fuse_ks_cvm = TRUE
  )

  if (!identical(observed$ks, expected_ks)) {
    stop(sprintf("KS differs in differential case %d (d=%d, n=%d, B=%d).", case_id, d, n, B))
  }
  if (!identical(observed$cvm, expected_cvm)) {
    stop(sprintf("CvM differs in differential case %d (d=%d, n=%d, B=%d).", case_id, d, n, B))
  }

  chunk_size <- sample.int(B, 1L)
  blocks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
  chunked_ks <- numeric(B)
  chunked_cvm <- numeric(B)
  for (idx in blocks) {
    block <- centered_weights[idx, , drop = FALSE]
    chunked <- compute_fast_sample_ks_cvm_stats_cpp(
      centered_weight_block = block,
      stream_prep = ks_stream,
      scale_factor = scale_factor,
      compute_ks = TRUE,
      compute_cvm = TRUE,
      fuse_ks_cvm = TRUE
    )
    chunked_ks[idx] <- chunked$ks
    chunked_cvm[idx] <- chunked$cvm
  }
  if (!identical(chunked_ks, expected_ks) || !identical(chunked_cvm, expected_cvm)) {
    stop(sprintf("Chunked output differs in case %d (chunk size %d).", case_id, chunk_size))
  }
}

cat(sprintf(
  "Validated %d differential cases: KS and CvM are bitwise identical for d in {2, 10}.\n",
  n_cases
))
