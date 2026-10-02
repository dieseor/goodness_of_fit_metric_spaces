rebuild_ks_prep_for_reestimated_fallback <- function(ks_prep,
                                                     data,
                                                     spec,
                                                     theta_hat,
                                                     control = list()) {
  if (is.null(ks_prep) || !isTRUE(ks_prep$light)) {
    return(ks_prep)
  }

  ks_grid <- if (identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
    make_sample_unique_distance_ks_grid()
  } else {
    list(
      omega_grid = ks_prep$omega_grid,
      t_grid = ks_prep$t_grid
    )
  }

  prepare_ks_observed_data(
    data = data,
    spec = spec,
    theta_hat = theta_hat,
    ks_grid = ks_grid,
    control = control,
    light = FALSE
  )
}
rebuild_cvm_prep_for_reestimated_fallback <- function(cvm_prep,
                                                      data,
                                                      spec,
                                                      theta_hat,
                                                      control = list()) {
  if (is.null(cvm_prep) || !isTRUE(cvm_prep$light)) {
    return(cvm_prep)
  }

  prepare_cvm_observed_data(
    data = data,
    spec = spec,
    theta_hat = theta_hat,
    control = control,
    light = FALSE
  )
}
run_fast_multiplier_bootstrap <- function(weight_matrix,
                                          spec,
                                          data,
                                          null,
                                          control,
                                          scale_factor,
                                          ks_prep = NULL,
                                          cvm_prep = NULL,
                                          want_ks = FALSE,
                                          want_cvm = FALSE,
                                          theta_hat,
                                          keep_bootstrap_thetas = FALSE,
                                          n_cores = 1L) {
  if (!identical(null$type, "composite")) {
    stop("The fast multiplier branch is only implemented for composite nulls.")
  }
  fast_backend_requested <- normalize_fast_multiplier_backend(
    control$fast_multiplier_backend %||% "cpp"
  )
  cpp_kernel_requested <- normalize_fast_multiplier_cpp_kernel(
    control$fast_multiplier_cpp_kernel %||% "contiguous_double"
  )
  fusion_requested <- normalize_fast_multiplier_fusion(
    control$fast_multiplier_fuse_ks_cvm %||% TRUE
  )
  cache_requested <- normalize_fast_multiplier_cache(
    control$fast_multiplier_cache_corrections %||% "auto"
  )
  fast_prep <- spec_fast_multiplier_prepare(
    spec = spec,
    data = data,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control
  )
  if (is.null(fast_prep)) {
    stop(sprintf(
      "Model '%s' does not expose the fast multiplier preparation hook required for `bootstrap_method = 'fast_multiplier'`.",
      spec$name
    ))
  }
  if (isTRUE(fast_prep$fallback_to_reestimated)) {
    ks_prep_fallback <- rebuild_ks_prep_for_reestimated_fallback(
      ks_prep = ks_prep,
      data = data,
      spec = spec,
      theta_hat = theta_hat,
      control = control
    )
    cvm_prep_fallback <- rebuild_cvm_prep_for_reestimated_fallback(
      cvm_prep = cvm_prep,
      data = data,
      spec = spec,
      theta_hat = theta_hat,
      control = control
    )
    fallback_chunks <- run_reestimated_bootstrap_chunks(
      weight_matrix = weight_matrix,
      spec = spec,
      data = data,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep_fallback,
      cvm_prep = cvm_prep_fallback,
      want_ks = want_ks,
      want_cvm = want_cvm,
      keep_bootstrap_thetas = keep_bootstrap_thetas,
      theta_hat = theta_hat,
      n_cores = n_cores
    )
    fallback_result <- list(
      ks = if (want_ks) {
        unlist(lapply(fallback_chunks, `[[`, "ks"), use.names = FALSE)
      } else {
        NULL
      },
      cvm = if (want_cvm) {
        unlist(lapply(fallback_chunks, `[[`, "cvm"), use.names = FALSE)
      } else {
        NULL
      },
      theta = if (keep_bootstrap_thetas) {
        unlist(lapply(fallback_chunks, `[[`, "theta"), recursive = FALSE, use.names = FALSE)
      } else {
        NULL
      },
      prep_seconds = sum(vapply(fallback_chunks, function(x) as.numeric(x$prep_seconds %||% 0), numeric(1))),
      loop_seconds = sum(vapply(fallback_chunks, function(x) as.numeric(x$loop_seconds %||% 0), numeric(1))),
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_,
      vhat_method = NA_character_,
      fast_ks_mode = NA_character_,
      fast_cvm_mode = NA_character_,
      S_obs_dim = NA_integer_,
      Psi_aux_dim = NA_integer_,
      D_ks_dim = NA_integer_,
      D_cvm_dim = NA_integer_,
      Vhat_dim = NA_integer_,
      score_mean_aux = NA_real_,
      score_mean_aux_norm = NA_real_,
      Vhat_eigenvalues = NA_real_,
      Vhat_rcond = NA_real_,
      Vhat_condition_number = NA_real_,
      fast_parameter_summary = NA_real_,
      fast_multiplier_backend_requested = fast_backend_requested,
      fast_multiplier_backend_effective = "r",
      fast_multiplier_cpp_kernel_requested = cpp_kernel_requested,
      fast_multiplier_cpp_kernel_effective = "not_used",
      fast_multiplier_fuse_ks_cvm_requested = fusion_requested,
      fast_multiplier_fuse_ks_cvm_effective = FALSE,
      fast_multiplier_cache_corrections_requested = cache_requested,
      fast_multiplier_cache_corrections_effective = FALSE,
      fallback_to_reestimated = TRUE,
      fallback_reason = fast_prep$fallback_reason %||% NA_character_,
      effective_bootstrap_method = "reestimated"
    )
    return(fallback_result)
  }

  prep_start <- proc.time()[["elapsed"]]
  S_obs <- as.matrix(fast_prep$S_obs)
  Vhat <- as.matrix(fast_prep$Vhat)
  if (nrow(S_obs) != spec_n_obs(spec, data, control)) {
    stop("The fast multiplier observed score matrix has incompatible dimensions.")
  }

  H_ks <- NULL
  ks_sample_stream_prep <- NULL
  if (want_ks) {
    if (is.list(fast_prep$D_ks) &&
        identical(fast_prep$D_ks$mode %||% "", "sample_points_unique_distances")) {
      ks_sample_stream_prep <- prepare_fast_ks_sample_stream_prep(
        S_obs = S_obs,
        Vhat = Vhat,
        Psi_aux = as.matrix(fast_prep$Psi_aux),
        ks_prep = ks_prep,
        D_ks_info = fast_prep$D_ks,
        control = control
      )
    } else {
      D_ks <- as.matrix(fast_prep$D_ks)
      Y_ks <- build_fast_ks_indicator_matrix(ks_prep)
      H_ks <- Y_ks - S_obs %*% fast_multiplier_solve_vhat(
        Vhat,
        t(D_ks),
        label = "the fast KS correction"
      )
    }
  }
  cvm_stream_prep <- NULL
  if (want_cvm) {
    cvm_stream_prep <- prepare_fast_cvm_stream_prep(
      S_obs = S_obs,
      Vhat = Vhat,
      D_cvm = fast_prep$D_cvm,
      observed_distance_matrix = fast_prep$observed_cvm_distance_matrix %||% cvm_prep$distance_matrix %||% NULL,
      Psi_aux = fast_prep$Psi_aux,
      cvm_prep = cvm_prep,
      correction_cache = if (is.list(fast_prep$D_cvm) &&
        isTRUE(fast_prep$D_cvm$shared_with_ks)) {
        ks_sample_stream_prep$correction_cache
      } else {
        NULL
      },
      control = control
    )
  }
  ks_sample_eligible <- !want_ks || (
    !is.null(ks_sample_stream_prep) &&
      identical(
        ks_sample_stream_prep$mode,
        "sample_points_unique_distances_streamed"
      )
  )
  cvm_sample_eligible <- !want_cvm || (
    !is.null(cvm_stream_prep) &&
      identical(
        cvm_stream_prep$mode,
        "sample_points_unique_distances_sorted_rows"
      )
  )
  shared_sample_stream <- if (want_ks && want_cvm &&
      ks_sample_eligible && cvm_sample_eligible) {
    identical(
      ks_sample_stream_prep$obs_order_matrix,
      cvm_stream_prep$obs_order_matrix
    ) &&
      identical(
        ks_sample_stream_prep$tie_end_matrix,
        cvm_stream_prep$tie_end_matrix
      ) &&
      identical(
        ks_sample_stream_prep$aux_order_matrix,
        cvm_stream_prep$aux_order_matrix
      )
  } else {
    TRUE
  }
  sample_backend_eligible <- ks_sample_eligible &&
    cvm_sample_eligible && shared_sample_stream
  joint_stream_prep <- if (sample_backend_eligible && want_ks) {
    ks_sample_stream_prep
  } else if (sample_backend_eligible && want_cvm) {
    cvm_stream_prep
  } else {
    NULL
  }
  fast_backend_effective <- if (
    identical(fast_backend_requested, "cpp") &&
      sample_backend_eligible &&
      !is.null(joint_stream_prep)
  ) {
    ensure_distance_profile_cpp_loaded()
    "cpp"
  } else {
    "r"
  }
  fusion_effective <- isTRUE(fusion_requested) &&
    want_ks && want_cvm && sample_backend_eligible
  if (want_ks && want_cvm &&
      (!isTRUE(cvm_prep$shared_with_ks) ||
       !isTRUE(shared_sample_stream) ||
       !isTRUE(fusion_effective))) {
    stop(
      paste(
        "Joint fast-multiplier KS and CvM were requested, but their shared",
        "sample preparation and fused evaluation are not active.",
        "Use the sample-points unique-distance KS grid, set",
        "`keep$observed_process = FALSE`, and leave",
        "`control$fast_multiplier_fuse_ks_cvm = TRUE`."
      ),
      call. = FALSE
    )
  }
  prep_seconds <- proc.time()[["elapsed"]] - prep_start

  n_reps <- nrow(weight_matrix)
  ks_values <- if (want_ks) numeric(n_reps) else NULL
  cvm_values <- if (want_cvm) numeric(n_reps) else NULL
  chunk_size <- control$fast_bootstrap_chunk_size %||% NULL
  if (is.null(chunk_size) &&
      ((want_ks &&
        !is.null(ks_sample_stream_prep) &&
        identical(ks_sample_stream_prep$mode, "sample_points_unique_distances_streamed")) ||
       (want_cvm &&
        !is.null(cvm_stream_prep) &&
        identical(cvm_stream_prep$mode, "sample_points_unique_distances_sorted_rows")))) {
    chunk_size <- as.integer(control$fast_multiplier_stream_chunk_size %||% 100L)
  }
  if (!is.null(chunk_size)) {
    chunk_size <- as.integer(chunk_size)
    if (!is.finite(chunk_size) || chunk_size <= 0L) {
      stop("`control$fast_bootstrap_chunk_size` must be a strictly positive integer when supplied.")
    }
  }
  n_cores <- max(1L, as.integer(n_cores))
  replicate_blocks <- if (is.null(chunk_size)) {
    split(seq_len(n_reps), rep(seq_len(min(n_cores, n_reps)), length.out = n_reps))
  } else {
    split(seq_len(n_reps), ceiling(seq_len(n_reps) / chunk_size))
  }
  run_fast_block <- function(block_indices) {
    out_ks <- if (want_ks) numeric(length(block_indices)) else NULL
    out_cvm <- if (want_cvm) numeric(length(block_indices)) else NULL
    centered_weight_block <- weight_matrix[block_indices, , drop = FALSE] - 1
    joint_stats <- NULL

    if (!is.null(joint_stream_prep) &&
        (identical(fast_backend_effective, "cpp") ||
          isTRUE(fusion_effective))) {
      joint_stats <- if (identical(fast_backend_effective, "cpp")) {
        compute_fast_sample_ks_cvm_stats_cpp(
          centered_weight_block = centered_weight_block,
          stream_prep = joint_stream_prep,
          scale_factor = scale_factor,
          compute_ks = want_ks,
          compute_cvm = want_cvm,
          fuse_ks_cvm = fusion_effective,
          cpp_kernel = cpp_kernel_requested
        )
      } else {
        compute_fast_sample_ks_cvm_stats_fused_r(
          centered_weight_block = centered_weight_block,
          stream_prep = joint_stream_prep,
          scale_factor = scale_factor,
          compute_ks = want_ks,
          compute_cvm = want_cvm
        )
      }
      if (want_ks) {
        out_ks[] <- joint_stats$ks
      }
      if (want_cvm) {
        out_cvm[] <- joint_stats$cvm
      }
    }

    if (is.null(joint_stats) && want_ks &&
        !is.null(ks_sample_stream_prep)) {
      out_ks[] <- compute_fast_ks_sample_stats_streamed(
        centered_weight_block = centered_weight_block,
        ks_sample_stream_prep = ks_sample_stream_prep,
        scale_factor = scale_factor,
        control = control
      )
    }

    for (j in seq_along(block_indices)) {
      b <- block_indices[[j]]
      centered_weights <- as.numeric(centered_weight_block[j, ])
      if (is.null(joint_stats) && want_ks &&
          is.null(ks_sample_stream_prep)) {
        process_ks <- scale_factor * drop(crossprod(centered_weights, H_ks)) / sqrt(nrow(S_obs))
        out_ks[[j]] <- max(abs(process_ks))
      }
    }
    if (is.null(joint_stats) && want_cvm) {
      out_cvm[] <- compute_fast_cvm_stats_streamed(
        centered_weight_block = centered_weight_block,
        cvm_stream_prep = cvm_stream_prep,
        scale_factor = scale_factor
      )
    }
    list(indices = block_indices, ks = out_ks, cvm = out_cvm)
  }
  show_progress <- isTRUE(control$progress_bar %||% FALSE)
  progress_label <- as.character(control$progress_label %||% "bootstrap")
  run_fast_blocks_with_progress <- function() {
    pb <- utils::txtProgressBar(
      min = 0,
      max = n_reps,
      style = 3,
      file = stderr()
    )
    cat(sprintf("\n[%s] ", progress_label), file = stderr())
    on.exit({
      close(pb)
      cat("\n", file = stderr())
    }, add = TRUE)

    results <- vector("list", length(replicate_blocks))
    completed_reps <- 0L
    update_progress <- function(reps) {
      completed_reps <<- completed_reps + as.integer(reps)
      utils::setTxtProgressBar(pb, completed_reps)
    }

    if (.Platform$OS.type != "unix") {
      for (i in seq_along(replicate_blocks)) {
        results[[i]] <- run_fast_block(replicate_blocks[[i]])
        update_progress(length(replicate_blocks[[i]]))
      }
      return(results)
    }

    active_jobs <- list()
    next_idx <- 1L
    launch_job <- function(i) {
      job <- parallel::mcparallel(
        run_fast_block(replicate_blocks[[i]]),
        detached = FALSE,
        silent = TRUE
      )
      active_jobs[[as.character(job$pid)]] <<- list(
        job = job,
        index = i,
        reps = length(replicate_blocks[[i]])
      )
    }

    while (next_idx <= length(replicate_blocks) &&
        length(active_jobs) < min(n_cores, length(replicate_blocks))) {
      launch_job(next_idx)
      next_idx <- next_idx + 1L
    }
    while (length(active_jobs) > 0L) {
      collected <- parallel::mccollect(
        jobs = lapply(active_jobs, `[[`, "job"),
        wait = TRUE,
        timeout = 0.5
      )
      if (is.null(collected) || length(collected) == 0L) {
        next
      }
      for (pid in names(collected)) {
        meta <- active_jobs[[pid]]
        value <- collected[[pid]]
        if (inherits(value, "try-error")) {
          stop(sprintf("Parallel fast bootstrap block failed: %s", as.character(value)))
        }
        results[[meta$index]] <- value
        update_progress(meta$reps)
        active_jobs[[pid]] <- NULL
        if (next_idx <= length(replicate_blocks)) {
          launch_job(next_idx)
          next_idx <- next_idx + 1L
        }
      }
    }
    results
  }
  loop_start <- proc.time()[["elapsed"]]
  block_results <- if (length(replicate_blocks) == 1L) {
    if (show_progress) {
      run_fast_blocks_with_progress()
    } else {
      list(run_fast_block(replicate_blocks[[1L]]))
    }
  } else if (show_progress) {
    run_fast_blocks_with_progress()
  } else if (.Platform$OS.type == "unix") {
    parallel::mclapply(
      replicate_blocks,
      run_fast_block,
      mc.cores = min(n_cores, length(replicate_blocks)),
      mc.preschedule = TRUE
    )
  } else {
    lapply(replicate_blocks, run_fast_block)
  }
  for (block_result in block_results) {
    idx <- block_result$indices
    if (want_ks) {
      ks_values[idx] <- block_result$ks
    }
    if (want_cvm) {
      cvm_values[idx] <- block_result$cvm
    }
  }
  loop_seconds <- proc.time()[["elapsed"]] - loop_start
  shared_correction_cache <- is.list(fast_prep$D_cvm) &&
    isTRUE(fast_prep$D_cvm$shared_with_ks) &&
    !is.null(ks_sample_stream_prep$correction_cache)
  correction_cache_bytes <- c(
    if (!is.null(ks_sample_stream_prep$correction_cache)) {
      ks_sample_stream_prep$correction_cache$bytes
    } else {
      0
    },
    if (!is.null(cvm_stream_prep$correction_cache) && !shared_correction_cache) {
      cvm_stream_prep$correction_cache$bytes
    } else {
      0
    }
  )

  list(
    ks = ks_values,
    cvm = cvm_values,
    theta = if (keep_bootstrap_thetas) vector("list", 0L) else NULL,
    prep_seconds = prep_seconds,
    loop_seconds = loop_seconds,
    derivative_method = fast_prep$derivative_method,
    derivative_method_requested = fast_prep$derivative_method_requested %||%
      fast_prep$derivative_method,
    derivative_method_effective = fast_prep$derivative_method_effective %||%
      fast_prep$derivative_method,
    derivative_method_selection_source =
      fast_prep$derivative_method_selection_source %||% NA_character_,
    derivative_mc_size = fast_prep$derivative_mc_size,
    derivative_mc_seed = fast_prep$derivative_mc_seed,
    quadrature_algorithm = paste(
      fast_prep$quadrature_diagnostics$algorithms_effective %||%
        fast_prep$quadrature_settings$algorithm %||% NA_character_,
      collapse = "+"
    ),
    quadrature_abs_tol = fast_prep$quadrature_settings$abs_tol %||% NA_real_,
    quadrature_max_terms = fast_prep$quadrature_settings$max_terms %||% NA_integer_,
    quadrature_initial_upper =
      fast_prep$quadrature_settings$initial_upper %||% NA_real_,
    quadrature_max_upper =
      fast_prep$quadrature_settings$max_upper %||% NA_real_,
    quadrature_tail_consecutive =
      fast_prep$quadrature_settings$tail_consecutive %||% NA_integer_,
    quadrature_eigen_rel_tol =
      fast_prep$quadrature_settings$eigen_rel_tol %||% NA_real_,
    quadrature_clip_tol = fast_prep$quadrature_settings$clip_tol %||% NA_real_,
    quadrature_center_evaluations =
      fast_prep$quadrature_diagnostics$center_evaluations %||% NA_integer_,
    quadrature_max_terms_used =
      fast_prep$quadrature_diagnostics$max_terms_used %||% NA_integer_,
    quadrature_max_residual_error_estimate =
      fast_prep$quadrature_diagnostics$max_residual_error_estimate %||% NA_real_,
    quadrature_max_condition_number =
      fast_prep$quadrature_diagnostics$max_condition_number %||% NA_real_,
    quadrature_max_propagated_error_estimate =
      fast_prep$quadrature_diagnostics$max_propagated_error_estimate %||% NA_real_,
    quadrature_max_upper_used =
      fast_prep$quadrature_diagnostics$max_upper_limit %||% NA_real_,
    quadrature_max_evaluations =
      fast_prep$quadrature_diagnostics$max_evaluations %||% NA_integer_,
    vhat_method = fast_prep$vhat_method %||% NA_character_,
    correction_representation = fast_prep$correction_representation %||% "score",
    paper_Vhat_method = fast_prep$paper_Vhat_method %||% NA_character_,
    paper_Vhat_eigenvalues = fast_prep$paper_Vhat_diagnostics$eigenvalues %||% NA_real_,
    paper_Vhat_rcond = fast_prep$paper_Vhat_diagnostics$rcond %||% NA_real_,
    paper_Vhat_condition_number =
      fast_prep$paper_Vhat_diagnostics$condition_number %||% NA_real_,
    fast_ks_mode = if (!is.null(ks_sample_stream_prep)) ks_sample_stream_prep$mode else if (!is.null(H_ks)) "dense_matrix" else NA_character_,
    fast_cvm_mode = cvm_stream_prep$mode %||% NA_character_,
    sample_correction_cache_bytes = sum(correction_cache_bytes),
    shared_sample_correction_cache = shared_correction_cache,
    S_obs_dim = fast_prep$vhat_diagnostics$S_obs_dim %||% NA_integer_,
    Psi_aux_dim = fast_prep$vhat_diagnostics$Psi_aux_dim %||% NA_integer_,
    D_ks_dim = if (!is.null(fast_prep$D_ks)) dim(fast_prep$D_ks) else NULL,
    D_cvm_dim = if (!is.null(fast_prep$D_cvm)) dim(fast_prep$D_cvm) else NULL,
    Vhat_dim = fast_prep$vhat_diagnostics$Vhat_dim %||% NA_integer_,
    score_mean_aux = fast_prep$vhat_diagnostics$score_mean_aux %||% NA_real_,
    score_mean_aux_norm = fast_prep$vhat_diagnostics$score_mean_aux_norm %||% NA_real_,
    Vhat_eigenvalues = fast_prep$vhat_diagnostics$Vhat_eigenvalues %||% NA_real_,
    Vhat_rcond = fast_prep$vhat_diagnostics$Vhat_rcond %||% NA_real_,
    Vhat_condition_number = fast_prep$vhat_diagnostics$Vhat_condition_number %||% NA_real_,
    fast_parameter_summary = fast_prep$vhat_diagnostics$par0 %||% NA_real_,
    fast_multiplier_backend_requested = fast_backend_requested,
    fast_multiplier_backend_effective = fast_backend_effective,
    fast_multiplier_cpp_kernel_requested = cpp_kernel_requested,
    fast_multiplier_cpp_kernel_effective = if (
      identical(fast_backend_effective, "cpp")
    ) cpp_kernel_requested else "not_used",
    fast_multiplier_fuse_ks_cvm_requested = fusion_requested,
    fast_multiplier_fuse_ks_cvm_effective = fusion_effective,
    fast_multiplier_cache_corrections_requested = cache_requested,
    fast_multiplier_cache_corrections_effective =
      any(correction_cache_bytes > 0),
    fallback_to_reestimated = FALSE,
    fallback_reason = NA_character_,
    effective_bootstrap_method = "fast_multiplier"
  )
}
run_reestimated_bootstrap_chunks <- function(weight_matrix,
                                             spec,
                                             data,
                                             null,
                                             control,
                                             scale_factor,
                                             ks_prep,
                                             cvm_prep,
                                             want_ks,
                                             want_cvm,
                                             fuse_sample_ks_cvm = FALSE,
                                             keep_bootstrap_thetas,
                                             theta_hat,
                                             n_cores = 1L) {
  n_reps <- nrow(weight_matrix)
  n_cores_effective <- min(max(1L, as.integer(n_cores)), n_reps)
  show_progress <- isTRUE(control$progress_bar %||% FALSE)
  progress_label <- as.character(control$progress_label %||% "bootstrap")
  chunk_size <- control$reestimated_bootstrap_chunk_size %||% NULL
  if (is.null(chunk_size)) {
    if (show_progress) {
      target_chunks <- max(n_cores_effective, 8L * n_cores_effective)
      chunk_size <- max(1L, ceiling(n_reps / target_chunks))
    } else {
      chunk_size <- ceiling(n_reps / n_cores_effective)
    }
  }
  chunk_size <- as.integer(chunk_size)
  if (!is.finite(chunk_size) || chunk_size <= 0L) {
    stop("`control$reestimated_bootstrap_chunk_size` must be a strictly positive integer when supplied.")
  }
  chunk_starts <- seq.int(1L, n_reps, by = chunk_size)
  chunk_ids <- lapply(chunk_starts, function(start_idx) {
    seq.int(start_idx, min(start_idx + chunk_size - 1L, n_reps))
  })
  weight_chunks <- lapply(chunk_ids, function(indices) {
    weight_matrix[indices, , drop = FALSE]
  })
  replicate_index_chunks <- unname(chunk_ids)
  task_fun <- function(i) {
    run_bootstrap_chunk(
      weight_chunk = weight_chunks[[i]],
      spec = spec,
      data = data,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      fuse_sample_ks_cvm = fuse_sample_ks_cvm,
      keep_bootstrap_thetas = keep_bootstrap_thetas,
      theta_start = theta_hat,
      replicate_indices = replicate_index_chunks[[i]]
    )
  }

  update_progress_bar <- function(pb, completed_reps) {
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, completed_reps)
    }
  }

  if (n_cores_effective == 1L) {
    pb <- if (show_progress) utils::txtProgressBar(min = 0, max = n_reps, style = 3,
                                                   file = stderr()) else NULL
    if (!is.null(pb)) {
      cat(sprintf("\n[%s] ", progress_label), file = stderr())
    }
    on.exit(if (!is.null(pb)) {
      close(pb)
      cat("\n", file = stderr())
    }, add = TRUE)
    completed_reps <- 0L
    out <- vector("list", length(weight_chunks))
    for (i in seq_along(weight_chunks)) {
      out[[i]] <- task_fun(i)
      completed_reps <- completed_reps + nrow(weight_chunks[[i]])
      update_progress_bar(pb, completed_reps)
    }
    return(out)
  }

  if (.Platform$OS.type == "unix") {
    if (!show_progress) {
      return(parallel::mclapply(
        seq_along(weight_chunks),
        task_fun,
        mc.cores = n_cores_effective,
        mc.preschedule = TRUE
      ))
    }

    pb <- utils::txtProgressBar(min = 0, max = n_reps, style = 3, file = stderr())
    cat(sprintf("\n[%s] ", progress_label), file = stderr())
    on.exit({
      close(pb)
      cat("\n", file = stderr())
    }, add = TRUE)

    results <- vector("list", length(weight_chunks))
    active_jobs <- list()
    next_idx <- 1L
    completed_reps <- 0L

    launch_job <- function(i) {
      job <- parallel::mcparallel(task_fun(i), detached = FALSE, silent = TRUE)
      active_jobs[[as.character(job$pid)]] <<- list(
        job = job,
        index = i,
        reps = nrow(weight_chunks[[i]])
      )
    }

    while (next_idx <= length(weight_chunks) && length(active_jobs) < n_cores_effective) {
      launch_job(next_idx)
      next_idx <- next_idx + 1L
    }

    while (length(active_jobs) > 0L) {
      collected <- parallel::mccollect(
        jobs = lapply(active_jobs, `[[`, "job"),
        wait = TRUE,
        timeout = 0.5
      )
      if (is.null(collected) || length(collected) == 0L) {
        next
      }

      for (pid in names(collected)) {
        meta <- active_jobs[[pid]]
        value <- collected[[pid]]
        if (inherits(value, "try-error")) {
          stop(sprintf("Parallel bootstrap chunk failed: %s", as.character(value)))
        }
        results[[meta$index]] <- value
        completed_reps <- completed_reps + meta$reps
        update_progress_bar(pb, completed_reps)
        active_jobs[[pid]] <- NULL

        if (next_idx <= length(weight_chunks)) {
          launch_job(next_idx)
          next_idx <- next_idx + 1L
        }
      }
    }

    return(results)
  }

  cl <- parallel::makeCluster(n_cores_effective)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, { library(gofmetric); NULL })

  worker_symbols <- c(
    "spec",
    "data",
    "control",
    "ks_prep",
    "cvm_prep",
    "want_ks",
    "want_cvm",
    "fuse_sample_ks_cvm",
    "scale_factor",
    "run_bootstrap_chunk",
    "compute_grid_weighted_profile",
    "compute_theoretical_profile_matrix",
    "compute_theoretical_sample_profile_matrix",
    "compute_weighted_sample_profile_matrix",
    "spec_observation_at",
    "grid_n_points",
    "grid_point_at",
    "ensure_profile_matrix",
    "null",
    "keep_bootstrap_thetas",
    "theta_hat",
    "clip_cardioid_dot_products",
    "normalize_cardioid_data",
    "normalize_cardioid_theta",
    "weighted_cardioid_resultant",
    "normalize_cardioid_mle_weights",
    "cardioid_distance_threshold",
    "theoretical_distance_profile_cardioid",
    "mle_sph_car_weighted",
    "fit_cardioid_theta",
    "normalize_small_circle_data",
    "normalize_small_circle_theta",
    "fit_small_circle_theta",
    "make_small_circle_spec"
  )

  parallel::clusterExport(cl, worker_symbols, envir = environment())
  parallel::clusterExport(cl, c("replicate_index_chunks", "weight_chunks"), envir = environment())

  parallel::parLapply(cl, seq_along(weight_chunks), function(i) {
    run_bootstrap_chunk(
      weight_chunk = weight_chunks[[i]],
      spec = spec,
      data = data,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      fuse_sample_ks_cvm = fuse_sample_ks_cvm,
      keep_bootstrap_thetas = keep_bootstrap_thetas,
      theta_start = theta_hat,
      replicate_indices = replicate_index_chunks[[i]]
    )
  })
}
