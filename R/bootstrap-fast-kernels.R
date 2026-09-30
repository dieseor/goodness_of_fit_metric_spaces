build_fast_ks_indicator_matrix <- function(ks_prep) {
  n <- nrow(ks_prep$distance_matrix)
  threshold_matrix <- matrix(
    rep(as.numeric(ks_prep$t_grid), each = n),
    nrow = n,
    ncol = length(ks_prep$t_grid)
  )
  indicator_blocks <- lapply(seq_len(ncol(ks_prep$distance_matrix)), function(k) {
    distance_block <- matrix(
      rep.int(ks_prep$distance_matrix[, k], length(ks_prep$t_grid)),
      nrow = n,
      ncol = length(ks_prep$t_grid),
      byrow = FALSE
    )
    distance_block <= threshold_matrix
  })

  do.call(cbind, indicator_blocks) * 1
}
resolve_fast_sample_correction_cache <- function(n_centers,
                                                 n_thresholds,
                                                 n_parameters,
                                                 control = list()) {
  n_centers <- as.integer(n_centers)
  n_thresholds <- as.integer(n_thresholds)
  n_parameters <- as.integer(n_parameters)
  if (any(!is.finite(c(n_centers, n_thresholds, n_parameters))) ||
      any(c(n_centers, n_thresholds, n_parameters) <= 0L)) {
    stop("Sample-correction cache dimensions must be strictly positive integers.")
  }

  requested <- normalize_fast_multiplier_cache(
    control$fast_multiplier_cache_corrections %||% "auto"
  )
  n_max <- as.integer(
    control$fast_multiplier_correction_cache_n_max %||% 500L
  )
  if (length(n_max) != 1L || !is.finite(n_max) || n_max <= 0L) {
    stop(
      "`control$fast_multiplier_correction_cache_n_max` must be a strictly positive integer."
    )
  }
  max_bytes <- as.numeric(control$fast_multiplier_correction_cache_max_bytes %||% (128 * 1024^2))
  if (!is.finite(max_bytes) || max_bytes <= 0) {
    stop("`control$fast_multiplier_correction_cache_max_bytes` must be a positive finite number.")
  }

  bytes <- as.double(n_centers) * as.double(n_thresholds) *
    as.double(n_parameters) * 8
  enabled <- if (identical(requested, "true")) {
    TRUE
  } else if (identical(requested, "false")) {
    FALSE
  } else {
    n_centers <= n_max && bytes <= max_bytes
  }

  list(
    enabled = enabled,
    bytes = bytes,
    requested = requested,
    max_bytes = max_bytes,
    n_max = n_max
  )
}
build_fast_sample_correction_cache <- function(Psi_aux_solved,
                                               aux_order_matrix,
                                               aux_sorted_distance_matrix,
                                               obs_sorted_distance_matrix,
                                               control = list()) {
  Psi_aux_solved <- as.matrix(Psi_aux_solved)
  aux_order_matrix <- as.matrix(aux_order_matrix)
  aux_sorted_distance_matrix <- as.matrix(aux_sorted_distance_matrix)
  obs_sorted_distance_matrix <- as.matrix(obs_sorted_distance_matrix)
  n_centers <- nrow(obs_sorted_distance_matrix)
  n_thresholds <- ncol(obs_sorted_distance_matrix)
  n_aux <- nrow(Psi_aux_solved)
  n_parameters <- ncol(Psi_aux_solved)

  if (!identical(dim(aux_order_matrix), c(n_centers, n_aux)) ||
      !identical(dim(aux_sorted_distance_matrix), c(n_centers, n_aux))) {
    stop("The auxiliary ordering data are incompatible with the observed sample centers.")
  }

  decision <- resolve_fast_sample_correction_cache(
    n_centers = n_centers,
    n_thresholds = n_thresholds,
    n_parameters = n_parameters,
    control = control
  )
  if (!isTRUE(decision$enabled)) {
    return(NULL)
  }

  values <- matrix(0, nrow = n_centers * n_thresholds, ncol = n_parameters)
  for (center_idx in seq_len(n_centers)) {
    aux_cumpsi_solved <- col_cumsums_base(
      Psi_aux_solved[aux_order_matrix[center_idx, ], , drop = FALSE]
    ) / n_aux
    aux_basis_full <- rbind(0, aux_cumpsi_solved)
    selected_counts <- findInterval(
      obs_sorted_distance_matrix[center_idx, ],
      aux_sorted_distance_matrix[center_idx, ]
    )
    idx <- ((center_idx - 1L) * n_thresholds + 1L):(center_idx * n_thresholds)
    values[idx, ] <- aux_basis_full[selected_counts + 1L, , drop = FALSE]
  }

  list(
    values = values,
    bytes = decision$bytes,
    requested = decision$requested,
    n_max = decision$n_max
  )
}
get_fast_sample_correction <- function(stream_prep, center_idx) {
  correction_cache <- stream_prep$correction_cache
  n_thresholds <- ncol(stream_prep$obs_sorted_distance_matrix)
  if (!is.null(correction_cache)) {
    idx <- ((center_idx - 1L) * n_thresholds + 1L):(center_idx * n_thresholds)
    return(correction_cache$values[idx, , drop = FALSE])
  }

  aux_order <- stream_prep$aux_order_matrix[center_idx, ]
  aux_cumpsi_solved <- col_cumsums_base(
    stream_prep$Psi_aux_solved[aux_order, , drop = FALSE]
  ) / stream_prep$n_aux
  aux_basis_full <- rbind(0, aux_cumpsi_solved)
  selected_counts <- findInterval(
    stream_prep$obs_sorted_distance_matrix[center_idx, ],
    stream_prep$aux_sorted_distance_matrix[center_idx, ]
  )
  aux_basis_full[selected_counts + 1L, , drop = FALSE]
}
prepare_fast_ks_sample_cache <- function(S_obs,
                                         Vhat,
                                         Psi_aux,
                                         ks_prep,
                                         D_ks_info) {
  if (!identical(D_ks_info$mode %||% "", "sample_points_unique_distances")) {
    stop("The sample-based KS cache requires `D_ks_info$mode = 'sample_points_unique_distances'`.")
  }

  n <- nrow(S_obs)
  n_aux <- nrow(Psi_aux)
  Vhat_inv <- solve(Vhat)
  n_omega <- ncol(ks_prep$distance_matrix)
  omega_cache <- vector("list", n_omega)

  for (j in seq_len(n_omega)) {
    aux_order <- D_ks_info$aux_order_matrix[j, ]
    aux_cumpsi <- col_cumsums_base(Psi_aux[aux_order, , drop = FALSE]) / n_aux
    aux_basis_full <- rbind(0, aux_cumpsi)
    aux_sorted_distances <- D_ks_info$aux_sorted_distance_matrix[j, ]
    obs_sorted_distances <- ks_prep$sorted_distance_matrix[j, ]
    aux_counts_sorted <- findInterval(obs_sorted_distances, aux_sorted_distances)
    correction_selected <- aux_basis_full[aux_counts_sorted + 1L, , drop = FALSE]

    omega_cache[[j]] <- list(
      obs_order = ks_prep$order_matrix[j, ],
      correction_basis = Vhat_inv %*% t(correction_selected)
    )
  }

  list(
    mode = "sample_points_unique_distances",
    omega_cache = omega_cache,
    n = n
  )
}
prepare_fast_ks_sample_stream_prep <- function(S_obs,
                                               Vhat,
                                               Psi_aux,
                                               ks_prep,
                                               D_ks_info,
                                               control = list()) {
  if (!identical(D_ks_info$mode %||% "", "sample_points_unique_distances")) {
    stop("The sample-based KS stream prep requires `D_ks_info$mode = 'sample_points_unique_distances'`.")
  }

  vhat_inverse <- fast_multiplier_solve_vhat(
    Vhat,
    diag(ncol(Vhat)),
    label = "the fast sample-KS inverse"
  )
  if (!is.null(D_ks_info$derivative_sorted)) {
    derivative_sorted <- as.matrix(D_ks_info$derivative_sorted)
    n_centers <- nrow(ks_prep$sorted_distance_matrix)
    n_thresholds <- ncol(ks_prep$sorted_distance_matrix)
    if (!identical(
      dim(derivative_sorted),
      c(n_centers * n_thresholds, ncol(S_obs))
    )) {
      stop("The deterministic sample-KS derivative table has incompatible dimensions.")
    }
    correction_values <- derivative_sorted %*% t(vhat_inverse)
    return(list(
      mode = "sample_points_unique_distances_streamed",
      S_obs = as.matrix(S_obs),
      Psi_aux_solved = NULL,
      aux_order_matrix = NULL,
      aux_sorted_distance_matrix = NULL,
      obs_order_matrix = ks_prep$order_matrix,
      obs_sorted_distance_matrix = ks_prep$sorted_distance_matrix,
      tie_end_matrix = build_sorted_tie_end_matrix(
        ks_prep$sorted_distance_matrix
      ),
      n = nrow(S_obs),
      n_aux = 0L,
      correction_cache = list(
        values = correction_values,
        bytes = as.double(length(correction_values)) * 8,
        requested = "deterministic",
        n_max = n_centers
      )
    ))
  }
  Psi_aux_solved <- as.matrix(Psi_aux) %*% t(vhat_inverse)
  list(
    mode = "sample_points_unique_distances_streamed",
    S_obs = as.matrix(S_obs),
    Psi_aux_solved = Psi_aux_solved,
    aux_order_matrix = D_ks_info$aux_order_matrix,
    aux_sorted_distance_matrix = D_ks_info$aux_sorted_distance_matrix,
    obs_order_matrix = ks_prep$order_matrix,
    obs_sorted_distance_matrix = ks_prep$sorted_distance_matrix,
    tie_end_matrix = build_sorted_tie_end_matrix(
      ks_prep$sorted_distance_matrix
    ),
    n = nrow(S_obs),
    n_aux = nrow(Psi_aux),
    correction_cache = build_fast_sample_correction_cache(
      Psi_aux_solved = Psi_aux_solved,
      aux_order_matrix = D_ks_info$aux_order_matrix,
      aux_sorted_distance_matrix = D_ks_info$aux_sorted_distance_matrix,
      obs_sorted_distance_matrix = ks_prep$sorted_distance_matrix,
      control = control
    )
  )
}
compute_fast_ks_sample_stats_reference <- function(centered_weight_block,
                                                   S_obs,
                                                   H_ks_sample_cache,
                                                   scale_factor) {
  score_block <- centered_weight_block %*% S_obs
  block_max <- rep.int(0, nrow(centered_weight_block))

  for (omega_info in H_ks_sample_cache$omega_cache) {
    ordered_weights <- centered_weight_block[, omega_info$obs_order, drop = FALSE]
    empirical_selected <- row_cumsums_base(ordered_weights)
    correction_selected <- score_block %*% omega_info$correction_basis
    process_selected <- scale_factor * (empirical_selected - correction_selected) /
      sqrt(H_ks_sample_cache$n)
    block_max <- pmax(block_max, apply(abs(process_selected), 1L, max))
  }

  block_max
}
compute_fast_ks_sample_stats_streamed <- function(centered_weight_block,
                                                  ks_sample_stream_prep,
                                                  scale_factor,
                                                  control = list()) {
  n_omega <- nrow(ks_sample_stream_prep$obs_order_matrix)
  omega_block_size <- normalize_ks_block_size(
    block_size = control$fast_multiplier_ks_block_size %||%
      control$ks_block_size %||% NULL,
    n_rows = n_omega,
    arg_name = "`control$fast_multiplier_ks_block_size`"
  )
  score_block <- centered_weight_block %*% ks_sample_stream_prep$S_obs
  block_max <- rep.int(0, nrow(centered_weight_block))

  for (block_start in seq.int(1L, n_omega, by = omega_block_size)) {
    for (j in block_start:min(block_start + omega_block_size - 1L, n_omega)) {
      obs_order <- ks_sample_stream_prep$obs_order_matrix[j, ]
      correction_selected_solved <- get_fast_sample_correction(ks_sample_stream_prep, j)
      ordered_weights <- centered_weight_block[, obs_order, drop = FALSE]
      empirical_selected <- row_cumsums_base(ordered_weights)
      empirical_selected <- empirical_selected[
        , ks_sample_stream_prep$tie_end_matrix[j, ], drop = FALSE
      ]
      correction_selected <- score_block %*% t(correction_selected_solved)
      process_selected <- scale_factor * (empirical_selected - correction_selected) /
        sqrt(ks_sample_stream_prep$n)
      block_max <- pmax(block_max, apply(abs(process_selected), 1L, max))
    }
  }

  block_max
}
compute_fast_sample_ks_cvm_stats_fused_r <- function(
    centered_weight_block,
    stream_prep,
    scale_factor,
    compute_ks = TRUE,
    compute_cvm = TRUE) {
  if (!isTRUE(compute_ks) && !isTRUE(compute_cvm)) {
    stop("At least one fast statistic must be requested.")
  }

  score_block <- centered_weight_block %*% stream_prep$S_obs
  n_reps <- nrow(centered_weight_block)
  n <- stream_prep$n
  n_centers <- nrow(stream_prep$obs_order_matrix)
  ks <- if (compute_ks) rep.int(0, n_reps) else NULL
  cvm_sum <- if (compute_cvm) rep.int(0, n_reps) else NULL

  for (center_idx in seq_len(n_centers)) {
    ordered_weights <- centered_weight_block[
      , stream_prep$obs_order_matrix[center_idx, ], drop = FALSE
    ]
    empirical_selected <- row_cumsums_base(ordered_weights)
    empirical_selected <- empirical_selected[
      , stream_prep$tie_end_matrix[center_idx, ], drop = FALSE
    ]
    correction_selected <- score_block %*% t(
      get_fast_sample_correction(stream_prep, center_idx)
    )
    process_selected <- scale_factor * (
      empirical_selected - correction_selected
    ) / sqrt(n)
    if (compute_ks) {
      ks <- pmax(ks, apply(abs(process_selected), 1L, max))
    }
    if (compute_cvm) {
      cvm_sum <- cvm_sum + rowSums(process_selected^2)
    }
  }

  list(
    ks = ks,
    cvm = if (compute_cvm) cvm_sum / (n * n) else NULL
  )
}
compute_fast_sample_ks_cvm_stats_cpp <- function(
    centered_weight_block,
    stream_prep,
    scale_factor,
    compute_ks = TRUE,
    compute_cvm = TRUE,
    fuse_ks_cvm = TRUE,
    cpp_kernel = "contiguous_double") {
  if (!isTRUE(compute_ks) && !isTRUE(compute_cvm)) {
    stop("At least one fast statistic must be requested.")
  }
  cpp_kernel <- normalize_fast_multiplier_cpp_kernel(cpp_kernel)
  ensure_distance_profile_cpp_loaded()
  score_block <- centered_weight_block %*% stream_prep$S_obs
  n_reps <- nrow(centered_weight_block)
  ks <- if (compute_ks) rep.int(0, n_reps) else NULL
  cvm_sum <- if (compute_cvm) rep.int(0, n_reps) else NULL

  call_kernel <- function(obs_order_matrix,
                          tie_end_matrix,
                          correction_matrix,
                          want_ks,
                          want_cvm) {
    distance_profile_cpp_call(
      if (identical(cpp_kernel, "contiguous_double")) {
        "cpp_fast_sample_ks_cvm_stats_contiguous_double"
      } else {
        "cpp_fast_sample_ks_cvm_stats"
      },
      centered_weights = centered_weight_block,
      score_block = score_block,
      obs_order_matrix = obs_order_matrix,
      tie_end_matrix = tie_end_matrix,
      correction_matrix = correction_matrix,
      scale_factor = scale_factor,
      compute_ks = want_ks,
      compute_cvm = want_cvm
    )
  }
  update_results <- function(value, want_ks, want_cvm) {
    if (want_ks) {
      ks <<- pmax(ks, value$ks)
    }
    if (want_cvm) {
      cvm_sum <<- cvm_sum + value$cvm_sum
    }
  }
  evaluate <- function(want_ks, want_cvm) {
    if (!is.null(stream_prep$correction_cache)) {
      update_results(
        call_kernel(
          obs_order_matrix = stream_prep$obs_order_matrix,
          tie_end_matrix = stream_prep$tie_end_matrix,
          correction_matrix = stream_prep$correction_cache$values,
          want_ks = want_ks,
          want_cvm = want_cvm
        ),
        want_ks = want_ks,
        want_cvm = want_cvm
      )
      return(invisible(NULL))
    }

    for (center_idx in seq_len(nrow(stream_prep$obs_order_matrix))) {
      update_results(
        call_kernel(
          obs_order_matrix = stream_prep$obs_order_matrix[
            center_idx, , drop = FALSE
          ],
          tie_end_matrix = stream_prep$tie_end_matrix[
            center_idx, , drop = FALSE
          ],
          correction_matrix = get_fast_sample_correction(
            stream_prep, center_idx
          ),
          want_ks = want_ks,
          want_cvm = want_cvm
        ),
        want_ks = want_ks,
        want_cvm = want_cvm
      )
    }
    invisible(NULL)
  }

  if (isTRUE(fuse_ks_cvm) || !isTRUE(compute_ks) || !isTRUE(compute_cvm)) {
    evaluate(compute_ks, compute_cvm)
  } else {
    evaluate(TRUE, FALSE)
    evaluate(FALSE, TRUE)
  }

  list(
    ks = ks,
    cvm = if (compute_cvm) {
      cvm_sum / (stream_prep$n * stream_prep$n)
    } else {
      NULL
    }
  )
}
compute_fast_ks_sample_stats_blocked <- function(centered_weight_block,
                                                 S_obs,
                                                 H_ks_sample_cache,
                                                 scale_factor,
                                                 control = list()) {
  n_omega <- length(H_ks_sample_cache$omega_cache)
  omega_block_size <- normalize_ks_block_size(
    block_size = control$fast_multiplier_ks_block_size %||%
      control$ks_block_size %||% NULL,
    n_rows = n_omega,
    arg_name = "`control$fast_multiplier_ks_block_size`"
  )
  score_block <- centered_weight_block %*% S_obs
  block_max <- rep.int(0, nrow(centered_weight_block))

  for (block_start in seq.int(1L, n_omega, by = omega_block_size)) {
    block_end <- min(block_start + omega_block_size - 1L, n_omega)
    omega_block <- H_ks_sample_cache$omega_cache[block_start:block_end]

    for (omega_info in omega_block) {
      ordered_weights <- centered_weight_block[, omega_info$obs_order, drop = FALSE]
      empirical_selected <- row_cumsums_base(ordered_weights)
      correction_selected <- score_block %*% omega_info$correction_basis
      process_selected <- scale_factor * (empirical_selected - correction_selected) /
        sqrt(H_ks_sample_cache$n)
      block_max <- pmax(block_max, apply(abs(process_selected), 1L, max))
    }
  }

  block_max
}
build_fast_cvm_H_block <- function(block_start,
                                   correction_obs,
                                   D_cvm,
                                   observed_distance_matrix,
                                   block_size) {
  n <- nrow(correction_obs)
  if (nrow(D_cvm) != n * n) {
    stop("The fast multiplier CvM derivative matrix has incompatible dimensions.")
  }
  observed_distance_matrix <- as.matrix(observed_distance_matrix)
  block_end <- min(block_start + block_size - 1L, n)
  block_rows <- block_end - block_start + 1L
  idx_flat <- ((block_start - 1L) * n + 1L):(block_end * n)
  D_block <- D_cvm[idx_flat, , drop = FALSE]
  Y_block <- matrix(0, nrow = n, ncol = block_rows * n)

  for (offset in seq_len(block_rows)) {
    center_idx <- block_start + offset - 1L
    thresholds <- observed_distance_matrix[center_idx, ]
    distance_to_center <- matrix(
      rep.int(observed_distance_matrix[center_idx, ], n),
      nrow = n,
      ncol = n,
      byrow = FALSE
    )
    Y_block[, ((offset - 1L) * n + 1L):(offset * n)] <- (distance_to_center <=
      matrix(thresholds, nrow = n, ncol = n, byrow = TRUE)) * 1
  }

  Y_block - correction_obs %*% t(D_block)
}
prepare_fast_cvm_H_blocks <- function(S_obs,
                                      Vhat,
                                      D_cvm,
                                      observed_distance_matrix,
                                      control = list()) {
  n <- nrow(S_obs)
  block_size <- as.integer(
    control$fast_multiplier_cvm_block_size %||%
      default_fast_multiplier_cvm_block_size(n)
  )
  if (!is.finite(block_size) || block_size <= 0L) {
    stop("`control$fast_multiplier_cvm_block_size` must be a strictly positive integer.")
  }
  correction_obs <- t(fast_multiplier_solve_vhat(
    t(Vhat),
    t(S_obs),
    label = "the fast CvM correction"
  ))

  H_blocks <- vector("list", ceiling(n / block_size))
  block_id <- 1L
  for (block_start in seq.int(1L, n, by = block_size)) {
    H_blocks[[block_id]] <- build_fast_cvm_H_block(
      block_start = block_start,
      correction_obs = correction_obs,
      D_cvm = D_cvm,
      observed_distance_matrix = observed_distance_matrix,
      block_size = block_size
    )
    block_id <- block_id + 1L
  }
  H_blocks
}
resolve_fast_cvm_h_cache <- function(n, control = list()) {
  n <- as.integer(n)
  if (!is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer for the fast CvM H cache.")
  }
  requested <- tolower(as.character(control$fast_multiplier_cache_cvm_h %||% "auto"))
  if (length(requested) != 1L || !requested %in% c("auto", "true", "false")) {
    stop("`control$fast_multiplier_cache_cvm_h` must be TRUE, FALSE, or 'auto'.")
  }
  max_bytes <- as.numeric(control$fast_multiplier_cvm_h_cache_max_bytes %||% (128 * 1024^2))
  if (!is.finite(max_bytes) || max_bytes <= 0) {
    stop("`control$fast_multiplier_cvm_h_cache_max_bytes` must be a positive finite number.")
  }
  bytes <- as.double(n) * as.double(n) * 8
  list(
    enabled = if (identical(requested, "true")) TRUE else if (identical(requested, "false")) FALSE else bytes <= max_bytes,
    bytes = bytes
  )
}
compute_fast_cvm_stat_chunked <- function(centered_weights,
                                          H_blocks,
                                          scale_factor) {
  if (!is.list(H_blocks) || length(H_blocks) == 0L) {
    stop("`H_blocks` must be a non-empty list of precomputed fast CvM blocks.")
  }
  cvm_sum <- 0

  for (H_block in H_blocks) {
    n <- nrow(H_block)
    process_block <- scale_factor * drop(crossprod(centered_weights, H_block)) / sqrt(n)
    cvm_sum <- cvm_sum + sum(process_block^2)
  }

  n <- nrow(H_blocks[[1L]])
  cvm_sum / (n * n)
}
prepare_fast_cvm_stream_prep <- function(S_obs,
                                         Vhat,
                                         D_cvm,
                                         observed_distance_matrix,
                                         Psi_aux = NULL,
                                         cvm_prep = NULL,
                                         correction_cache = NULL,
                                         control = list()) {
  n <- nrow(S_obs)
  block_size <- as.integer(
    control$fast_multiplier_cvm_block_size %||%
      default_fast_multiplier_cvm_block_size(n)
  )
  if (!is.finite(block_size) || block_size <= 0L) {
    stop("`control$fast_multiplier_cvm_block_size` must be a strictly positive integer.")
  }

  if (is.list(D_cvm) &&
      identical(D_cvm$mode %||% "", "sample_points_unique_distances_sorted_rows")) {
    if (is.null(cvm_prep)) {
      stop("The streamed fast CvM prep requires `cvm_prep` in lightweight mode.")
    }
    vhat_inverse <- fast_multiplier_solve_vhat(
      Vhat,
      diag(ncol(Vhat)),
      label = "the fast CvM inverse"
    )
    if (!is.null(D_cvm$derivative_sorted)) {
      derivative_sorted <- as.matrix(D_cvm$derivative_sorted)
      n_centers <- nrow(cvm_prep$sorted_distance_matrix)
      n_thresholds <- ncol(cvm_prep$sorted_distance_matrix)
      if (!identical(
        dim(derivative_sorted),
        c(n_centers * n_thresholds, ncol(S_obs))
      )) {
        stop("The deterministic streamed CvM derivative table has incompatible dimensions.")
      }
      if (is.null(correction_cache)) {
        correction_values <- derivative_sorted %*% t(vhat_inverse)
        correction_cache <- list(
          values = correction_values,
          bytes = as.double(length(correction_values)) * 8,
          requested = "deterministic",
          n_max = n_centers
        )
      }
      return(list(
        mode = "sample_points_unique_distances_sorted_rows",
        S_obs = as.matrix(S_obs),
        Psi_aux_solved = NULL,
        aux_order_matrix = NULL,
        aux_sorted_distance_matrix = NULL,
        obs_order_matrix = cvm_prep$order_matrix,
        obs_sorted_distance_matrix = cvm_prep$sorted_distance_matrix,
        tie_end_matrix = build_sorted_tie_end_matrix(
          cvm_prep$sorted_distance_matrix
        ),
        n = n,
        n_aux = 0L,
        block_size = block_size,
        correction_cache = correction_cache
      ))
    }
    if (is.null(Psi_aux)) {
      stop("The Monte Carlo streamed CvM prep requires `Psi_aux`.")
    }
    Psi_aux_solved <- as.matrix(Psi_aux) %*% t(vhat_inverse)
    if (is.null(correction_cache)) {
      correction_cache <- build_fast_sample_correction_cache(
        Psi_aux_solved = Psi_aux_solved,
        aux_order_matrix = D_cvm$aux_order_matrix,
        aux_sorted_distance_matrix = D_cvm$aux_sorted_distance_matrix,
        obs_sorted_distance_matrix = cvm_prep$sorted_distance_matrix,
        control = control
      )
    }
    return(list(
      mode = "sample_points_unique_distances_sorted_rows",
      S_obs = as.matrix(S_obs),
      Psi_aux_solved = Psi_aux_solved,
      aux_order_matrix = D_cvm$aux_order_matrix,
      aux_sorted_distance_matrix = D_cvm$aux_sorted_distance_matrix,
      obs_order_matrix = cvm_prep$order_matrix,
      obs_sorted_distance_matrix = cvm_prep$sorted_distance_matrix,
      tie_end_matrix = build_sorted_tie_end_matrix(
        cvm_prep$sorted_distance_matrix
      ),
      n = n,
      n_aux = nrow(Psi_aux),
      block_size = block_size,
      correction_cache = correction_cache
    ))
  }

  if (nrow(D_cvm) != n * n) {
    stop("The fast multiplier CvM derivative matrix has incompatible dimensions.")
  }
  correction_obs <- t(fast_multiplier_solve_vhat(
    t(Vhat),
    t(S_obs),
    label = "the fast CvM correction"
  ))

  h_cache <- resolve_fast_cvm_h_cache(n, control)
  H_blocks <- if (isTRUE(h_cache$enabled)) {
    prepare_fast_cvm_H_blocks(
      S_obs = S_obs,
      Vhat = Vhat,
      D_cvm = D_cvm,
      observed_distance_matrix = observed_distance_matrix,
      control = control
    )
  } else {
    NULL
  }

  list(
    mode = "dense_matrix",
    correction_obs = correction_obs,
    D_cvm = D_cvm,
    observed_distance_matrix = as.matrix(observed_distance_matrix),
    n = n,
    block_size = block_size,
    H_blocks = H_blocks,
    H_cache_bytes = if (isTRUE(h_cache$enabled)) h_cache$bytes else 0
  )
}
compute_fast_cvm_stats_streamed <- function(centered_weight_block,
                                            cvm_stream_prep,
                                            scale_factor) {
  if (identical(cvm_stream_prep$mode %||% "dense_matrix", "sample_points_unique_distances_sorted_rows")) {
    score_block <- centered_weight_block %*% cvm_stream_prep$S_obs
    cvm_sum <- rep.int(0, nrow(centered_weight_block))
    n <- cvm_stream_prep$n

    for (block_start in seq.int(1L, n, by = cvm_stream_prep$block_size)) {
      block_end <- min(block_start + cvm_stream_prep$block_size - 1L, n)
      for (center_idx in block_start:block_end) {
        obs_order <- cvm_stream_prep$obs_order_matrix[center_idx, ]
        correction_selected_solved <- get_fast_sample_correction(cvm_stream_prep, center_idx)
        ordered_weights <- centered_weight_block[, obs_order, drop = FALSE]
        empirical_selected <- row_cumsums_base(ordered_weights)
        empirical_selected <- empirical_selected[
          , cvm_stream_prep$tie_end_matrix[center_idx, ], drop = FALSE
        ]
        process_center <- scale_factor * (
          empirical_selected -
            score_block %*% t(correction_selected_solved)
        ) / sqrt(n)
        cvm_sum <- cvm_sum + rowSums(process_center^2)
      }
    }

    return(cvm_sum / (n * n))
  }

  n <- cvm_stream_prep$n
  cvm_sum <- rep.int(0, nrow(centered_weight_block))

  if (!is.null(cvm_stream_prep$H_blocks)) {
    for (H_block in cvm_stream_prep$H_blocks) {
      process_block <- scale_factor * centered_weight_block %*% H_block / sqrt(n)
      cvm_sum <- cvm_sum + rowSums(process_block^2)
    }
    return(cvm_sum / (n * n))
  }

  for (block_start in seq.int(1L, n, by = cvm_stream_prep$block_size)) {
    H_block <- build_fast_cvm_H_block(
      block_start = block_start,
      correction_obs = cvm_stream_prep$correction_obs,
      D_cvm = cvm_stream_prep$D_cvm,
      observed_distance_matrix = cvm_stream_prep$observed_distance_matrix,
      block_size = cvm_stream_prep$block_size
    )
    process_block <- scale_factor * centered_weight_block %*% H_block / sqrt(n)
    cvm_sum <- cvm_sum + rowSums(process_block^2)
  }

  cvm_sum / (n * n)
}
