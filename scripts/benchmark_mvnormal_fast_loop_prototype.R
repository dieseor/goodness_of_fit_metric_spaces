#!/usr/bin/env Rscript

source(file.path("bootstrap", "multiplier_bootstrap.R"))

ensure_distance_profile_cpp_loaded()

parse_integer_vector <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 0L) return(default)
  as.integer(strsplit(sub(prefix, "", hit[[1L]], fixed = TRUE), ",", fixed = TRUE)[[1L]])
}

dimensions <- parse_integer_vector("dimensions", c(2L, 10L))
sample_sizes <- parse_integer_vector("sizes", c(50L, 100L, 200L))
B <- parse_integer_vector("B", 500L)[[1L]]
repetitions <- parse_integer_vector("repetitions", 5L)[[1L]]
chunk_sizes <- parse_integer_vector("chunk", 100L)
derivative_mc_size <- parse_integer_vector("derivative", 1000L)[[1L]]

make_case <- function(d, n) {
  set.seed(71000L + 100L * d + n)
  Sigma <- diag(d)
  if (d >= 2L) Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.35
  mu <- seq(-0.2, 0.2, length.out = d)
  data <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
  null <- list(type = "composite")
  spec <- make_mvnormal_spec(unknown_param = "both")
  control <- list(
    derivative_mc_size = derivative_mc_size,
    derivative_mc_seed = 81000L + 100L * d + n,
    fast_multiplier_stream_chunk_size = chunk_sizes[[1L]],
    fast_multiplier_cache_corrections = "true",
    logistic_gaussian_quadform_method = "farebrother"
  )
  theta_hat <- spec$fit_theta(data, weights = NULL, null = null, control = control)
  # Build only the ordering structures used by the fast loop.  The observed
  # theoretical profile is intentionally excluded: it is identical for both
  # loop implementations and is dominated by the generalized chi-square CDF.
  distance_matrix <- spec$distance_matrix(data, data, control)
  sorted_rows <- sort_distance_matrix_rows(distance_matrix)
  ks_prep <- list(
    ks_grid_mode = "sample_points_unique_distances",
    omega_grid = data,
    distance_matrix = distance_matrix,
    order_matrix = sorted_rows$order_matrix,
    sorted_distance_matrix = sorted_rows$sorted_distance_matrix,
    light = TRUE,
    statistic = NA_real_
  )
  cvm_prep <- list(
    order_matrix = sorted_rows$order_matrix,
    sorted_distance_matrix = sorted_rows$sorted_distance_matrix,
    statistic = NA_real_,
    light = TRUE,
    shared_with_ks = TRUE
  )
  fast_prep <- spec_fast_multiplier_prepare(
    spec = spec,
    data = data,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control
  )
  ks_stream <- prepare_fast_ks_sample_stream_prep(
    S_obs = fast_prep$S_obs,
    Vhat = fast_prep$Vhat,
    Psi_aux = fast_prep$Psi_aux,
    ks_prep = ks_prep,
    D_ks_info = fast_prep$D_ks,
    control = control
  )
  cvm_stream <- prepare_fast_cvm_stream_prep(
    S_obs = fast_prep$S_obs,
    Vhat = fast_prep$Vhat,
    D_cvm = fast_prep$D_cvm,
    observed_distance_matrix = cvm_prep$distance_matrix %||% NULL,
    Psi_aux = fast_prep$Psi_aux,
    cvm_prep = cvm_prep,
    correction_cache = ks_stream$correction_cache,
    control = control
  )
  no_cache_control <- utils::modifyList(
    control,
    list(fast_multiplier_cache_corrections = "false")
  )
  ks_stream_no_cache <- prepare_fast_ks_sample_stream_prep(
    S_obs = fast_prep$S_obs,
    Vhat = fast_prep$Vhat,
    Psi_aux = fast_prep$Psi_aux,
    ks_prep = ks_prep,
    D_ks_info = fast_prep$D_ks,
    control = no_cache_control
  )
  cvm_stream_no_cache <- prepare_fast_cvm_stream_prep(
    S_obs = fast_prep$S_obs,
    Vhat = fast_prep$Vhat,
    D_cvm = fast_prep$D_cvm,
    observed_distance_matrix = cvm_prep$distance_matrix %||% NULL,
    Psi_aux = fast_prep$Psi_aux,
    cvm_prep = cvm_prep,
    correction_cache = NULL,
    control = no_cache_control
  )
  set.seed(91000L + 100L * d + n)
  raw <- matrix(rexp(B * n), nrow = B, ncol = n)
  weights <- raw / rowMeans(raw)
  centered <- weights - 1
  tie_end_matrix <- t(vapply(seq_len(n), function(j) {
    sorted_tie_end_positions(ks_stream$obs_sorted_distance_matrix[j, ])
  }, integer(n)))

  list(
    d = d,
    n = n,
    centered = centered,
    ks_stream = ks_stream,
    cvm_stream = cvm_stream,
    ks_stream_no_cache = ks_stream_no_cache,
    cvm_stream_no_cache = cvm_stream_no_cache,
    tie_end_matrix = tie_end_matrix
  )
}

run_r <- function(case, chunk_size, cached = TRUE) {
  ks_stream <- if (cached) case$ks_stream else case$ks_stream_no_cache
  cvm_stream <- if (cached) case$cvm_stream else case$cvm_stream_no_cache
  blocks <- split(seq_len(nrow(case$centered)), ceiling(seq_len(nrow(case$centered)) / chunk_size))
  ks <- numeric(nrow(case$centered))
  cvm <- numeric(nrow(case$centered))
  for (idx in blocks) {
    block <- case$centered[idx, , drop = FALSE]
    ks[idx] <- compute_fast_ks_sample_stats_streamed(block, ks_stream, 1)
    cvm[idx] <- compute_fast_cvm_stats_streamed(block, cvm_stream, 1)
  }
  list(ks = ks, cvm = cvm)
}

run_r_fused <- function(case, chunk_size) {
  blocks <- split(seq_len(nrow(case$centered)), ceiling(seq_len(nrow(case$centered)) / chunk_size))
  ks <- numeric(nrow(case$centered))
  cvm <- numeric(nrow(case$centered))
  for (idx in blocks) {
    block <- case$centered[idx, , drop = FALSE]
    value <- compute_fast_sample_ks_cvm_stats_fused_r(
      centered_weight_block = block,
      stream_prep = case$ks_stream,
      scale_factor = 1,
      compute_ks = TRUE,
      compute_cvm = TRUE
    )
    ks[idx] <- value$ks
    cvm[idx] <- value$cvm
  }
  list(ks = ks, cvm = cvm)
}

run_cpp <- function(case, chunk_size) {
  blocks <- split(seq_len(nrow(case$centered)), ceiling(seq_len(nrow(case$centered)) / chunk_size))
  ks <- numeric(nrow(case$centered))
  cvm <- numeric(nrow(case$centered))
  for (idx in blocks) {
    block <- case$centered[idx, , drop = FALSE]
    value <- compute_fast_sample_ks_cvm_stats_cpp(
      centered_weight_block = block,
      stream_prep = case$ks_stream,
      scale_factor = 1,
      compute_ks = TRUE,
      compute_cvm = TRUE,
      fuse_ks_cvm = TRUE
    )
    ks[idx] <- value$ks
    cvm[idx] <- value$cvm
  }
  list(ks = ks, cvm = cvm)
}

rows <- list()
row_id <- 1L
for (d in dimensions) {
  for (n in sample_sizes) {
    message(sprintf("Preparing d=%d, n=%d", d, n))
    case <- make_case(d, n)
    for (chunk_size in chunk_sizes) {
      reference <- run_r(case, chunk_size = chunk_size, cached = TRUE)
      fused_reference <- run_r_fused(case, chunk_size = chunk_size)
      benchmark_uncached <- identical(as.integer(chunk_size), 100L)
      uncached_reference <- if (benchmark_uncached) {
        run_r(case, chunk_size = chunk_size, cached = FALSE)
      } else {
        reference
      }
      candidate <- run_cpp(case, chunk_size = chunk_size)
      identical_ks <- identical(candidate$ks, reference$ks) &&
        identical(fused_reference$ks, reference$ks) &&
        identical(uncached_reference$ks, reference$ks)
      identical_cvm <- identical(candidate$cvm, reference$cvm) &&
        identical(fused_reference$cvm, reference$cvm) &&
        identical(uncached_reference$cvm, reference$cvm)
      max_abs_ks <- max(abs(candidate$ks - reference$ks))
      max_abs_cvm <- max(abs(candidate$cvm - reference$cvm))
      if (!identical(fused_reference$ks, reference$ks) ||
          !identical(fused_reference$cvm, reference$cvm)) {
        message(sprintf(
          "R-fused mismatch d=%d n=%d chunk=%d: values KS=%s CvM=%s; attrs KS=%s CvM=%s",
          d, n, chunk_size,
          identical(as.numeric(fused_reference$ks), as.numeric(reference$ks)),
          identical(as.numeric(fused_reference$cvm), as.numeric(reference$cvm)),
          identical(attributes(fused_reference$ks), attributes(reference$ks)),
          identical(attributes(fused_reference$cvm), attributes(reference$cvm))
        ))
      }

      for (rep_id in seq_len(repetitions)) {
        backend_order <- switch(
          as.character((rep_id - 1L) %% 4L + 1L),
          `1` = c("r_uncached", "r_cached", "r_fused", "cpp"),
          `2` = c("cpp", "r_uncached", "r_cached", "r_fused"),
          `3` = c("r_fused", "r_cached", "cpp", "r_uncached"),
          `4` = c("r_cached", "cpp", "r_uncached", "r_fused")
        )
        if (!benchmark_uncached) backend_order <- backend_order[backend_order != "r_uncached"]
        elapsed <- numeric(4L)
        names(elapsed) <- c("r_uncached", "r_cached", "r_fused", "cpp")
        elapsed[["r_uncached"]] <- NA_real_
        for (backend in backend_order) {
          elapsed[[backend]] <- unname(system.time(
            switch(
              backend,
              r_uncached = run_r(case, chunk_size = chunk_size, cached = FALSE),
              r_cached = run_r(case, chunk_size = chunk_size, cached = TRUE),
              r_fused = run_r_fused(case, chunk_size = chunk_size),
              cpp = run_cpp(case, chunk_size = chunk_size)
            )
          )[["elapsed"]])
        }
        rows[[row_id]] <- data.frame(
          d = d,
          n = n,
          B = B,
          chunk_size = chunk_size,
          repetition = rep_id,
          r_uncached_seconds = elapsed[["r_uncached"]],
          r_cached_seconds = elapsed[["r_cached"]],
          r_fused_seconds = elapsed[["r_fused"]],
          cpp_seconds = elapsed[["cpp"]],
          cache_speedup = elapsed[["r_uncached"]] / elapsed[["r_cached"]],
          cpp_speedup_vs_uncached = elapsed[["r_uncached"]] / elapsed[["cpp"]],
          cpp_speedup_vs_cached = elapsed[["r_cached"]] / elapsed[["cpp"]],
          cpp_speedup_vs_r_fused = elapsed[["r_fused"]] / elapsed[["cpp"]],
          identical_ks = identical_ks,
          identical_cvm = identical_cvm,
          max_abs_ks = max_abs_ks,
          max_abs_cvm = max_abs_cvm
        )
        row_id <- row_id + 1L
      }
    }
  }
}

raw_results <- do.call(rbind, rows)
summary_results <- do.call(rbind, lapply(
  split(raw_results, interaction(raw_results$d, raw_results$n, raw_results$chunk_size, drop = TRUE)),
  function(x) data.frame(
    d = x$d[[1L]],
    n = x$n[[1L]],
    B = x$B[[1L]],
    chunk_size = x$chunk_size[[1L]],
    repetitions = nrow(x),
    median_r_uncached_seconds = if (all(is.na(x$r_uncached_seconds))) NA_real_ else median(x$r_uncached_seconds, na.rm = TRUE),
    median_r_cached_seconds = median(x$r_cached_seconds),
    median_r_fused_seconds = median(x$r_fused_seconds),
    median_cpp_seconds = median(x$cpp_seconds),
    median_cache_speedup = if (all(is.na(x$cache_speedup))) NA_real_ else median(x$cache_speedup, na.rm = TRUE),
    median_cpp_speedup_vs_uncached = if (all(is.na(x$cpp_speedup_vs_uncached))) NA_real_ else median(x$cpp_speedup_vs_uncached, na.rm = TRUE),
    median_cpp_speedup_vs_cached = median(x$cpp_speedup_vs_cached),
    median_cpp_speedup_vs_r_fused = median(x$cpp_speedup_vs_r_fused),
    identical_ks = all(x$identical_ks),
    identical_cvm = all(x$identical_cvm),
    max_abs_ks = max(x$max_abs_ks),
    max_abs_cvm = max(x$max_abs_cvm)
  )
))
row.names(summary_results) <- NULL

print(summary_results)

output_dir <- file.path("output", "mvnormal_fast_loop_prototype")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(raw_results, file.path(output_dir, "raw.csv"), row.names = FALSE)
write.csv(summary_results, file.path(output_dir, "summary.csv"), row.names = FALSE)
