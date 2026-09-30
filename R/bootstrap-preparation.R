
resolve_multiplier_bootstrap_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}
normalize_requested_statistics <- function(statistics) {
  statistics <- unique(tolower(as.character(statistics)))
  valid_statistics <- c("ks", "cvm")

  if (length(statistics) == 0) {
    stop("`statistics` cannot be empty.")
  }
  if (!all(statistics %in% valid_statistics)) {
    stop("`statistics` must be a subset of c('ks', 'cvm').")
  }

  statistics
}
normalize_fast_multiplier_backend <- function(backend = c("cpp", "r")) {
  backend <- tolower(as.character(backend))
  if (length(backend) > 1L) {
    backend <- backend[[1L]]
  }
  if (length(backend) != 1L || is.na(backend) ||
      !backend %in% c("cpp", "r")) {
    stop("`fast_multiplier_backend` must be either 'cpp' or 'r'.")
  }
  backend
}
normalize_fast_multiplier_cpp_kernel <- function(kernel = c(
                                                "contiguous_double",
                                                "legacy")) {
  kernel <- tolower(as.character(kernel))
  if (length(kernel) > 1L) {
    kernel <- kernel[[1L]]
  }
  if (length(kernel) != 1L || is.na(kernel) ||
      !kernel %in% c("legacy", "contiguous_double")) {
    stop(
      "`fast_multiplier_cpp_kernel` must be either 'legacy' or 'contiguous_double'."
    )
  }
  kernel
}
normalize_fast_multiplier_fusion <- function(fuse_ks_cvm = TRUE) {
  if (length(fuse_ks_cvm) != 1L || is.na(fuse_ks_cvm) ||
      !is.logical(fuse_ks_cvm)) {
    stop("`fuse_ks_cvm` must be TRUE or FALSE.")
  }
  isTRUE(fuse_ks_cvm)
}
normalize_fast_multiplier_cache <- function(cache_block_corrections = c(
                                               "auto", "true", "false")) {
  cache_block_corrections <- tolower(as.character(cache_block_corrections))
  if (length(cache_block_corrections) > 1L) {
    cache_block_corrections <- cache_block_corrections[[1L]]
  }
  if (length(cache_block_corrections) != 1L ||
      is.na(cache_block_corrections) ||
      !cache_block_corrections %in% c("auto", "true", "false")) {
    stop(
      "`cache_block_corrections` must be TRUE, FALSE, or one of 'auto', 'true', and 'false'."
    )
  }
  cache_block_corrections
}
normalize_keep_options <- function(keep) {
  defaults <- list(
    observed_process = TRUE,
    bootstrap_statistics = TRUE,
    bootstrap_thetas = FALSE
  )

  if (is.null(keep)) {
    return(defaults)
  }

  keep <- utils::modifyList(defaults, keep)
  keep$observed_process <- isTRUE(keep$observed_process)
  keep$bootstrap_statistics <- isTRUE(keep$bootstrap_statistics)
  keep$bootstrap_thetas <- isTRUE(keep$bootstrap_thetas)
  keep
}
default_fast_multiplier_cvm_block_size <- function(n) {
  n <- as.integer(n)
  if (!is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }

  max(1L, floor(n / 20L))
}
validate_null_object <- function(null) {
  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list with a field `type`.")
  }

  null$type <- tolower(as.character(null$type))
  if (!null$type %in% c("simple", "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }
  if (identical(null$type, "simple") && is.null(null$theta)) {
    stop("Simple nulls require `null$theta`.")
  }

  null
}
resolve_multiplier_spec <- function(multipliers = NULL) {
  if (is.null(multipliers)) {
    return(list(
      name = "Exp(1)",
      generator = function(n) stats::rexp(n, rate = 1),
      mean = 1,
      sd = 1
    ))
  }

  if (!is.list(multipliers) || !is.function(multipliers$generator)) {
    stop("`multipliers` must be NULL or a list with a `generator` function.")
  }

  mean_value <- as.numeric(multipliers$mean)
  sd_value <- as.numeric(multipliers$sd)

  if (length(mean_value) != 1L || !is.finite(mean_value) || mean_value <= 0) {
    stop("`multipliers$mean` must be a strictly positive finite scalar.")
  }
  if (length(sd_value) != 1L || !is.finite(sd_value) || sd_value <= 0) {
    stop("`multipliers$sd` must be a strictly positive finite scalar.")
  }

  list(
    name = multipliers$name %||% "custom",
    generator = multipliers$generator,
    mean = mean_value,
    sd = sd_value
  )
}
generate_multiplier_matrix <- function(B, n, multiplier_spec, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  output <- matrix(0, nrow = B, ncol = n)
  for (b in seq_len(B)) {
    raw_draw <- as.numeric(multiplier_spec$generator(n))
    if (length(raw_draw) != n) {
      stop("Multiplier generator returned an object of incompatible length.")
    }
    if (any(!is.finite(raw_draw))) {
      stop("Multiplier generator returned non-finite values.")
    }
    if (any(raw_draw < 0)) {
      stop("Multiplier generator returned negative values.")
    }
    output[b, ] <- raw_draw
  }

  output
}
normalize_multiplier_weights <- function(raw_multipliers) {
  raw_multipliers <- as.numeric(raw_multipliers)
  if (length(raw_multipliers) == 0) {
    stop("`raw_multipliers` cannot be empty.")
  }
  if (any(!is.finite(raw_multipliers))) {
    stop("`raw_multipliers` must be finite.")
  }
  if (any(raw_multipliers < 0)) {
    stop("`raw_multipliers` must be nonnegative.")
  }

  multiplier_mean <- mean(raw_multipliers)
  if (multiplier_mean <= 0) {
    stop("The sampled multiplier mean must be strictly positive.")
  }

  raw_multipliers / multiplier_mean
}
ensure_profile_matrix <- function(values, n_rows, n_cols) {
  matrix(as.numeric(values), nrow = n_rows, ncol = n_cols)
}
grid_n_points <- function(omega_grid) {
  if (is.matrix(omega_grid) || is.data.frame(omega_grid)) {
    return(nrow(omega_grid))
  }
  if (is.list(omega_grid) && is.null(dim(omega_grid))) {
    return(length(omega_grid))
  }
  length(omega_grid)
}
grid_point_at <- function(omega_grid, idx) {
  if (is.matrix(omega_grid) || is.data.frame(omega_grid)) {
    return(as.numeric(omega_grid[idx, , drop = TRUE]))
  }
  if (is.list(omega_grid) && is.null(dim(omega_grid))) {
    return(omega_grid[[idx]])
  }
  omega_grid[[idx]]
}
make_sample_unique_distance_ks_grid <- function() {
  list(mode = "sample_points_unique_distances")
}
is_sample_unique_distance_ks_grid <- function(ks_grid) {
  is.list(ks_grid) &&
    identical(
      tolower(as.character(ks_grid$mode %||% "")),
      "sample_points_unique_distances"
    )
}
build_sample_omega_grid <- function(spec, data, control = list()) {
  normalized_data <- spec_normalize_data(spec, data, control)

  if (is.matrix(normalized_data) || is.data.frame(normalized_data)) {
    return(as.matrix(normalized_data))
  }

  n <- spec_n_obs_normalized(spec, normalized_data, control)
  do.call(rbind, lapply(seq_len(n), function(i) {
    as.numeric(spec_observation_at_normalized(spec, normalized_data, i, control))
  }))
}
derive_sample_ks_grid_components <- function(data, spec, control = list()) {
  normalized_data <- spec_normalize_data(spec, data, control)
  omega_grid <- build_sample_omega_grid(spec, data = normalized_data, control = control)
  distance_matrix <- spec$distance_matrix(normalized_data, omega_grid, control)

  list(
    omega_grid = omega_grid,
    distance_matrix = distance_matrix,
    normalized_data = normalized_data
  )
}
row_cumsums_base <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  t(apply(mat, 1L, cumsum))
}
col_cumsums_base <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  matrix(
    as.numeric(apply(mat, 2L, cumsum)),
    nrow = nrow(mat),
    ncol = ncol(mat)
  )
}
sorted_tie_end_positions <- function(sorted_values) {
  sorted_values <- as.numeric(sorted_values)
  if (length(sorted_values) == 0L) {
    integer(0)
  } else {
    as.integer(findInterval(sorted_values, sorted_values))
  }
}
build_sorted_tie_end_matrix <- function(sorted_distance_matrix) {
  sorted_distance_matrix <- as.matrix(sorted_distance_matrix)
  t(vapply(seq_len(nrow(sorted_distance_matrix)), function(i) {
    sorted_tie_end_positions(sorted_distance_matrix[i, ])
  }, integer(ncol(sorted_distance_matrix))))
}
sort_distance_matrix_rows <- function(distance_matrix) {
  distance_matrix <- as.matrix(distance_matrix)
  n_rows <- nrow(distance_matrix)
  n_cols <- ncol(distance_matrix)

  order_matrix <- t(vapply(seq_len(n_rows), function(i) {
    as.integer(order(distance_matrix[i, ]))
  }, integer(n_cols)))

  sorted_distance_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(n_rows)) {
    sorted_distance_matrix[i, ] <- distance_matrix[i, order_matrix[i, ]]
  }

  list(
    order_matrix = order_matrix,
    sorted_distance_matrix = sorted_distance_matrix
  )
}
compute_sorted_empirical_profile_block <- function(sorted_distance_matrix,
                                                   row_indices) {
  sorted_block <- sorted_distance_matrix[row_indices, , drop = FALSE]
  n_total <- ncol(sorted_block)
  output <- matrix(0, nrow = nrow(sorted_block), ncol = n_total)

  for (k in seq_len(nrow(sorted_block))) {
    output[k, ] <- sorted_tie_end_positions(sorted_block[k, ]) / n_total
  }

  output
}
compute_sorted_weighted_profile_block <- function(order_matrix,
                                                  sorted_distance_matrix,
                                                  centered_weights,
                                                  row_indices) {
  order_block <- order_matrix[row_indices, , drop = FALSE]
  sorted_block <- sorted_distance_matrix[row_indices, , drop = FALSE]
  block_n <- nrow(order_block)
  n_total <- ncol(order_block)
  ordered_weight_block <- matrix(
    centered_weights[order_block],
    nrow = block_n,
    ncol = n_total
  )
  cumulative_block <- row_cumsums_base(ordered_weight_block)
  output <- matrix(0, nrow = block_n, ncol = n_total)

  for (k in seq_len(block_n)) {
    tie_end <- sorted_tie_end_positions(sorted_block[k, ])
    output[k, ] <- cumulative_block[k, tie_end]
  }

  output
}
compute_grid_empirical_profile <- function(distance_matrix,
                                           t_grid,
                                           sorted_distance_matrix = NULL,
                                           threshold_index_matrix = NULL) {
  n <- nrow(distance_matrix)
  n_omega <- ncol(distance_matrix)

  if (!is.null(threshold_index_matrix)) {
    profile_values <- threshold_index_matrix / n
    return(ensure_profile_matrix(profile_values, n_rows = n_omega, n_cols = length(t_grid)))
  }

  profile_values <- vapply(t_grid, function(t_value) {
    colMeans(distance_matrix <= t_value)
  }, numeric(n_omega))

  ensure_profile_matrix(profile_values, n_rows = n_omega, n_cols = length(t_grid))
}
compute_grid_weighted_profile <- function(distance_matrix,
                                          t_grid,
                                          normalized_weights,
                                          sorted_distance_matrix = NULL,
                                          order_matrix = NULL,
                                          threshold_index_matrix = NULL) {
  n <- nrow(distance_matrix)
  n_omega <- ncol(distance_matrix)

  if (is.null(sorted_distance_matrix) || is.null(order_matrix) || is.null(threshold_index_matrix)) {
    profile_values <- vapply(t_grid, function(t_value) {
      colSums((distance_matrix <= t_value) * normalized_weights) / n
    }, numeric(n_omega))

    return(ensure_profile_matrix(profile_values, n_rows = n_omega, n_cols = length(t_grid)))
  }

  profile_values <- matrix(0, nrow = n_omega, ncol = length(t_grid))

  for (j in seq_len(n_omega)) {
    ordered_weights <- normalized_weights[order_matrix[j, ]]
    cumulative_weights <- cumsum(ordered_weights)
    threshold_indices <- threshold_index_matrix[j, ]
    row_values <- numeric(length(threshold_indices))
    positive <- threshold_indices > 0L
    if (any(positive)) {
      row_values[positive] <- cumulative_weights[threshold_indices[positive]] / n
    }
    profile_values[j, ] <- row_values
  }

  ensure_profile_matrix(profile_values, n_rows = n_omega, n_cols = length(t_grid))
}
precompute_ks_grid_cache <- function(distance_matrix, t_grid) {
  n <- nrow(distance_matrix)
  n_omega <- ncol(distance_matrix)

  order_matrix <- t(vapply(seq_len(n_omega), function(j) {
    as.integer(order(distance_matrix[, j]))
  }, integer(n)))

  sorted_distance_matrix <- matrix(0, nrow = n_omega, ncol = n)
  for (j in seq_len(n_omega)) {
    sorted_distance_matrix[j, ] <- distance_matrix[order_matrix[j, ], j]
  }

  threshold_index_matrix <- t(vapply(seq_len(n_omega), function(j) {
    as.integer(findInterval(t_grid, sorted_distance_matrix[j, ]))
  }, integer(length(t_grid))))

  list(
    order_matrix = order_matrix,
    sorted_distance_matrix = sorted_distance_matrix,
    threshold_index_matrix = threshold_index_matrix
  )
}
compute_theoretical_profile_matrix <- function(spec, omega_grid, t_grid, theta, control = list()) {
  fast_output <- spec_profile_matrix_eval(
    spec = spec,
    omega_grid = omega_grid,
    t_grid = t_grid,
    theta = theta,
    control = control
  )
  if (!is.null(fast_output)) {
    return(ensure_profile_matrix(
      fast_output,
      n_rows = grid_n_points(omega_grid),
      n_cols = length(t_grid)
    ))
  }

  n_omega <- grid_n_points(omega_grid)
  n_t <- length(t_grid)

  output <- matrix(0, nrow = n_omega, ncol = n_t)
  for (i in seq_len(n_omega)) {
    omega_i <- grid_point_at(omega_grid, i)
    output[i, ] <- as.numeric(spec$profile_eval(omega_i, t_grid, theta, control))
  }

  output
}
prepare_ks_observed_data <- function(data,
                                     spec,
                                     theta_hat,
                                     ks_grid,
                                     control = list(),
                                     light = FALSE,
                                     share_cvm_statistic = FALSE) {
  if (!is.list(ks_grid)) {
    stop("KS requires `ks_grid = list(omega_grid = ..., t_grid = ...)`.")
  }

  if (is_sample_unique_distance_ks_grid(ks_grid)) {
    derived_grid <- derive_sample_ks_grid_components(data = data, spec = spec, control = control)
    omega_grid <- derived_grid$omega_grid
    distance_matrix <- derived_grid$distance_matrix
    ks_grid_mode <- "sample_points_unique_distances"
    sorted_rows <- sort_distance_matrix_rows(distance_matrix)
    order_matrix <- sorted_rows$order_matrix
    sorted_distance_matrix <- sorted_rows$sorted_distance_matrix
    n <- nrow(sorted_distance_matrix)

    if (isTRUE(light)) {
      shared_statistics <- if (isTRUE(share_cvm_statistic)) {
        compute_sample_ks_cvm_observed_stats_light(
          spec = spec,
          normalized_data = derived_grid$normalized_data,
          sorted_distance_matrix = sorted_distance_matrix,
          theta = theta_hat,
          control = control
        )
      } else {
        list(
          ks = compute_sample_ks_observed_stat_light(
            spec = spec,
            normalized_data = derived_grid$normalized_data,
            sorted_distance_matrix = sorted_distance_matrix,
            theta = theta_hat,
            control = control
          ),
          cvm = NULL
        )
      }

      return(list(
        ks_grid_mode = ks_grid_mode,
        omega_grid = omega_grid,
        t_grid = NULL,
        order_matrix = order_matrix,
        sorted_distance_matrix = sorted_distance_matrix,
        statistic = shared_statistics$ks,
        shared_cvm_statistic = shared_statistics$cvm,
        light = TRUE
      ))
    }

    rank_matrix <- t(vapply(seq_len(n), function(i) {
      as.integer(rank(distance_matrix[i, ], ties.method = "max"))
    }, integer(n)))
    row_index_matrix <- matrix(rep.int(seq_len(n), n), nrow = n, ncol = n)
    rank_linear_index <- row_index_matrix + (rank_matrix - 1L) * n
    empirical_profile <- rank_matrix / n
    theoretical_profile <- compute_theoretical_sample_profile_matrix(
      spec = spec,
      data = data,
      distance_matrix = distance_matrix,
      theta = theta_hat,
      control = control
    )
    process_matrix <- sqrt(n) * (empirical_profile - theoretical_profile)

    return(list(
      ks_grid_mode = ks_grid_mode,
      omega_grid = omega_grid,
      t_grid = NULL,
      distance_matrix = distance_matrix,
      rank_matrix = rank_matrix,
      order_matrix = order_matrix,
      sorted_distance_matrix = sorted_distance_matrix,
      rank_linear_index = rank_linear_index,
      empirical_profile = empirical_profile,
      theoretical_profile = theoretical_profile,
      process_matrix = process_matrix,
      statistic = max(abs(process_matrix)),
      light = FALSE
    ))
  } else {
    if (is.null(ks_grid$omega_grid) || is.null(ks_grid$t_grid)) {
      stop("KS requires `ks_grid = list(omega_grid = ..., t_grid = ...)`.")
    }
    omega_grid <- ks_grid$omega_grid
    t_grid <- as.numeric(ks_grid$t_grid)
    if (length(t_grid) == 0) {
      stop("`ks_grid$t_grid` cannot be empty.")
    }
    distance_matrix <- spec$distance_matrix(data, omega_grid, control)
    ks_grid_mode <- "fixed"
  }

  ks_cache <- precompute_ks_grid_cache(distance_matrix, t_grid)
  empirical_profile <- compute_grid_empirical_profile(
    distance_matrix,
    t_grid,
    sorted_distance_matrix = ks_cache$sorted_distance_matrix,
    threshold_index_matrix = ks_cache$threshold_index_matrix
  )
  theoretical_profile <- compute_theoretical_profile_matrix(
    spec = spec,
    omega_grid = omega_grid,
    t_grid = t_grid,
    theta = theta_hat,
    control = control
  )

  n <- spec_n_obs(spec, data, control)
  process_matrix <- sqrt(n) * (empirical_profile - theoretical_profile)

  list(
    ks_grid_mode = ks_grid_mode,
    omega_grid = omega_grid,
    t_grid = t_grid,
    distance_matrix = distance_matrix,
    order_matrix = ks_cache$order_matrix,
    sorted_distance_matrix = ks_cache$sorted_distance_matrix,
    threshold_index_matrix = ks_cache$threshold_index_matrix,
    empirical_profile = empirical_profile,
    theoretical_profile = theoretical_profile,
    process_matrix = process_matrix,
    statistic = max(abs(process_matrix)),
    light = FALSE
  )
}
compute_theoretical_sample_profile_matrix <- function(spec,
                                                     data,
                                                     distance_matrix,
                                                     theta,
                                                     control = list()) {
  debug_memory_log(
    control,
    sprintf("compute_theoretical_sample_profile_matrix: enter spec=%s", spec$name),
    list(distance_matrix = distance_matrix)
  )
  fast_output <- spec_sample_profile_matrix_eval(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta,
    control = control
  )
  if (!is.null(fast_output)) {
    n <- nrow(distance_matrix)
    debug_memory_log(control, "compute_theoretical_sample_profile_matrix: fast_output", list(fast_output = fast_output))
    return(ensure_profile_matrix(fast_output, n_rows = n, n_cols = n))
  }

  n <- nrow(distance_matrix)
  output <- matrix(0, nrow = n, ncol = n)
  normalized_data <- spec_normalize_data(spec, data, control)

  for (i in seq_len(n)) {
    omega_i <- spec_observation_at_normalized(spec, normalized_data, i, control)
    output[i, ] <- as.numeric(spec$profile_eval(omega_i, distance_matrix[i, ], theta, control))
  }

  output
}
compute_theoretical_sample_profile_sorted_block <- function(spec,
                                                            normalized_data,
                                                            sorted_distance_matrix,
                                                            theta,
                                                            row_indices,
                                                            prepared = NULL,
                                                            control = list()) {
  fast_output <- spec_sample_profile_sorted_block_eval(
    spec = spec,
    data = normalized_data,
    sorted_distance_matrix = sorted_distance_matrix,
    theta = theta,
    row_indices = row_indices,
    prepared = prepared,
    control = control
  )
  if (!is.null(fast_output)) {
    return(ensure_profile_matrix(
      fast_output,
      n_rows = length(row_indices),
      n_cols = ncol(sorted_distance_matrix)
    ))
  }

  output <- matrix(
    0,
    nrow = length(row_indices),
    ncol = ncol(sorted_distance_matrix)
  )

  for (k in seq_along(row_indices)) {
    i <- row_indices[[k]]
    omega_i <- spec_observation_at_normalized(
      spec,
      normalized_data,
      i,
      control
    )
    output[k, ] <- as.numeric(
      spec$profile_eval(
        omega_i,
        sorted_distance_matrix[i, ],
        theta,
        control
      )
    )
  }

  output
}
normalize_observed_profile_n_cores <- function(n_cores = 1L,
                                               n_blocks = 1L) {
  n_cores <- as.integer(n_cores)
  n_blocks <- as.integer(n_blocks)

  if (length(n_cores) != 1L ||
      !is.finite(n_cores) ||
      n_cores < 1L) {
    stop(
      "`control$observed_profile_n_cores` must be a positive integer."
    )
  }
  if (length(n_blocks) != 1L ||
      !is.finite(n_blocks) ||
      n_blocks < 1L) {
    stop("`n_blocks` must be a positive integer.")
  }

  min(n_cores, n_blocks)
}
make_observed_profile_row_blocks <- function(n_rows, block_size) {
  starts <- seq.int(1L, n_rows, by = block_size)
  lapply(starts, function(block_start) {
    block_end <- min(block_start + block_size - 1L, n_rows)
    block_start:block_end
  })
}
map_observed_profile_blocks <- function(row_blocks,
                                        worker,
                                        n_cores = 1L) {
  if (!is.list(row_blocks) || length(row_blocks) == 0L) {
    stop("`row_blocks` must be a non-empty list.")
  }
  if (!is.function(worker)) {
    stop("`worker` must be a function.")
  }

  n_cores <- normalize_observed_profile_n_cores(
    n_cores = n_cores,
    n_blocks = length(row_blocks)
  )

  results <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(
      row_blocks,
      worker,
      mc.cores = n_cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(row_blocks, worker)
  }

  failed <- vapply(
    results,
    inherits,
    logical(1L),
    what = "try-error"
  )
  if (any(failed)) {
    stop(
      sprintf(
        "Observed-profile worker failed: %s",
        as.character(results[[which(failed)[[1L]]]])
      ),
      call. = FALSE
    )
  }

  results
}
compute_sample_ks_observed_stat_light <- function(spec,
                                                  normalized_data,
                                                  sorted_distance_matrix,
                                                  theta,
                                                  control = list()) {
  n <- nrow(sorted_distance_matrix)
  block_size <- normalize_ks_block_size(
    block_size = control$ks_block_size %||% NULL,
    n_rows = n
  )
  prepared <- spec_sample_profile_sorted_prepare(
    spec = spec,
    data = normalized_data,
    sorted_distance_matrix = sorted_distance_matrix,
    theta = theta,
    control = control
  )
  row_blocks <- make_observed_profile_row_blocks(n, block_size)

  block_results <- map_observed_profile_blocks(
    row_blocks = row_blocks,
    n_cores = control$observed_profile_n_cores %||% 1L,
    worker = function(row_indices) {
      empirical_block <- compute_sorted_empirical_profile_block(
        sorted_distance_matrix = sorted_distance_matrix,
        row_indices = row_indices
      )
      theoretical_block <- compute_theoretical_sample_profile_sorted_block(
        spec = spec,
        normalized_data = normalized_data,
        sorted_distance_matrix = sorted_distance_matrix,
        theta = theta,
        row_indices = row_indices,
        prepared = prepared,
        control = control
      )
      process_block <- sqrt(n) * (empirical_block - theoretical_block)
      max(abs(process_block))
    }
  )

  max(unlist(block_results, use.names = FALSE))
}
compute_sample_ks_cvm_observed_stats_light <- function(
    spec,
    normalized_data,
    sorted_distance_matrix,
    theta,
    control = list()) {
  n <- nrow(sorted_distance_matrix)
  ks_block_size <- normalize_ks_block_size(
    block_size = control$ks_block_size %||% NULL,
    n_rows = n
  )
  cvm_block_size <- normalize_ks_block_size(
    block_size = control$cvm_block_size %||% control$ks_block_size %||% NULL,
    n_rows = n,
    arg_name = "`control$cvm_block_size`"
  )
  block_size <- min(ks_block_size, cvm_block_size)
  prepared <- spec_sample_profile_sorted_prepare(
    spec = spec,
    data = normalized_data,
    sorted_distance_matrix = sorted_distance_matrix,
    theta = theta,
    control = control
  )
  row_blocks <- make_observed_profile_row_blocks(n, block_size)

  block_results <- map_observed_profile_blocks(
    row_blocks = row_blocks,
    n_cores = control$observed_profile_n_cores %||% 1L,
    worker = function(row_indices) {
      empirical_block <- compute_sorted_empirical_profile_block(
        sorted_distance_matrix = sorted_distance_matrix,
        row_indices = row_indices
      )
      theoretical_block <- compute_theoretical_sample_profile_sorted_block(
        spec = spec,
        normalized_data = normalized_data,
        sorted_distance_matrix = sorted_distance_matrix,
        theta = theta,
        row_indices = row_indices,
        prepared = prepared,
        control = control
      )
      process_block <- sqrt(n) * (empirical_block - theoretical_block)
      list(
        ks = max(abs(process_block)),
        cvm_sum = sum(process_block^2)
      )
    }
  )

  list(
    ks = max(vapply(block_results, `[[`, numeric(1L), "ks")),
    cvm = sum(vapply(block_results, `[[`, numeric(1L), "cvm_sum")) /
      (n * n)
  )
}
compute_cvm_observed_stat_light <- function(spec,
                                            normalized_data,
                                            sorted_distance_matrix,
                                            theta,
                                            control = list()) {
  n <- nrow(sorted_distance_matrix)
  block_size <- normalize_ks_block_size(
    block_size = control$cvm_block_size %||% control$ks_block_size %||% NULL,
    n_rows = n,
    arg_name = "`control$cvm_block_size`"
  )
  prepared <- spec_sample_profile_sorted_prepare(
    spec = spec,
    data = normalized_data,
    sorted_distance_matrix = sorted_distance_matrix,
    theta = theta,
    control = control
  )
  row_blocks <- make_observed_profile_row_blocks(n, block_size)

  block_results <- map_observed_profile_blocks(
    row_blocks = row_blocks,
    n_cores = control$observed_profile_n_cores %||% 1L,
    worker = function(row_indices) {
      empirical_block <- compute_sorted_empirical_profile_block(
        sorted_distance_matrix = sorted_distance_matrix,
        row_indices = row_indices
      )
      theoretical_block <- compute_theoretical_sample_profile_sorted_block(
        spec = spec,
        normalized_data = normalized_data,
        sorted_distance_matrix = sorted_distance_matrix,
        theta = theta,
        row_indices = row_indices,
        prepared = prepared,
        control = control
      )
      process_block <- sqrt(n) * (empirical_block - theoretical_block)
      sum(process_block^2)
    }
  )

  sum(unlist(block_results, use.names = FALSE)) / (n * n)
}
prepare_cvm_observed_data <- function(data,
                                      spec,
                                      theta_hat,
                                      control = list(),
                                      light = FALSE) {
  # A specialised prep may discard the ordering data required by the streamed
  # fast multiplier path.  When only the statistic is requested, use the
  # generic lightweight representation instead.
  fast_prep <- if (isTRUE(light)) NULL else {
    spec_cvm_prepare(spec, data = data, theta_hat = theta_hat, control = control)
  }
  if (!is.null(fast_prep)) {
    return(fast_prep)
  }

  debug_memory_log(control, "prepare_cvm_observed_data: before distance_matrix")
  distance_matrix <- spec$distance_matrix(data, data, control)
  debug_memory_log(control, "prepare_cvm_observed_data: after distance_matrix", list(distance_matrix = distance_matrix))
  sorted_rows <- sort_distance_matrix_rows(distance_matrix)
  order_matrix <- sorted_rows$order_matrix
  sorted_distance_matrix <- sorted_rows$sorted_distance_matrix
  n <- nrow(distance_matrix)

  if (isTRUE(light)) {
    statistic <- compute_cvm_observed_stat_light(
      spec = spec,
      normalized_data = spec_normalize_data(spec, data, control),
      sorted_distance_matrix = sorted_distance_matrix,
      theta = theta_hat,
      control = control
    )

    return(list(
      order_matrix = order_matrix,
      sorted_distance_matrix = sorted_distance_matrix,
      statistic = statistic,
      light = TRUE
    ))
  }

  rank_matrix <- t(vapply(seq_len(n), function(i) {
    as.integer(rank(distance_matrix[i, ], ties.method = "max"))
  }, integer(n)))
  debug_memory_log(control, "prepare_cvm_observed_data: after rank_matrix", list(rank_matrix = rank_matrix))

  order_list <- lapply(seq_len(n), function(i) {
    order(distance_matrix[i, ])
  })
  row_index_matrix <- matrix(rep.int(seq_len(n), n), nrow = n, ncol = n)
  rank_linear_index <- row_index_matrix + (rank_matrix - 1L) * n
  debug_memory_log(
    control,
    "prepare_cvm_observed_data: after ordering structures",
    list(
      order_matrix = order_matrix,
      rank_linear_index = rank_linear_index
    )
  )

  empirical_profile <- rank_matrix / n
  debug_memory_log(control, "prepare_cvm_observed_data: after empirical_profile", list(empirical_profile = empirical_profile))
  theoretical_profile <- compute_theoretical_sample_profile_matrix(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta_hat,
    control = control
  )
  debug_memory_log(control, "prepare_cvm_observed_data: after theoretical_profile", list(theoretical_profile = theoretical_profile))

  process_matrix <- sqrt(n) * (empirical_profile - theoretical_profile)
  debug_memory_log(control, "prepare_cvm_observed_data: after process_matrix", list(process_matrix = process_matrix))

  list(
    distance_matrix = distance_matrix,
    rank_matrix = rank_matrix,
    order_list = order_list,
    order_matrix = order_matrix,
    sorted_distance_matrix = sorted_distance_matrix,
    rank_linear_index = rank_linear_index,
    empirical_profile = empirical_profile,
    theoretical_profile = theoretical_profile,
    process_matrix = process_matrix,
    statistic = mean(process_matrix^2),
    light = FALSE
  )
}
prepare_cvm_observed_data_from_sample_ks <- function(data,
                                                      spec,
                                                      theta_hat,
                                                      ks_prep,
                                                      control = list()) {
  if (!isTRUE(ks_prep$light) ||
      !identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
    stop("A lightweight sample-based KS preparation is required to share its CvM cache.")
  }

  normalized_data <- spec_normalize_data(spec, data, control)
  n <- spec_n_obs_normalized(spec, normalized_data, control)
  if (!identical(dim(ks_prep$order_matrix), c(n, n)) ||
      !identical(dim(ks_prep$sorted_distance_matrix), c(n, n))) {
    stop("The KS ordering cache has incompatible dimensions for the CvM statistic.")
  }

  list(
    order_matrix = ks_prep$order_matrix,
    sorted_distance_matrix = ks_prep$sorted_distance_matrix,
    statistic = if (!is.null(ks_prep$shared_cvm_statistic)) {
      ks_prep$shared_cvm_statistic
    } else {
      compute_cvm_observed_stat_light(
        spec = spec,
        normalized_data = normalized_data,
        sorted_distance_matrix = ks_prep$sorted_distance_matrix,
        theta = theta_hat,
        control = control
      )
    },
    light = TRUE,
    shared_with_ks = TRUE
  )
}
prepare_cvm_observed_data_from_full_sample_ks <- function(ks_prep) {
  required <- c(
    "distance_matrix", "rank_matrix", "order_matrix", "sorted_distance_matrix",
    "rank_linear_index", "empirical_profile", "theoretical_profile",
    "process_matrix", "statistic", "ks_grid_mode", "light"
  )
  if (!all(required %in% names(ks_prep)) || isTRUE(ks_prep$light) ||
      !identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
    stop("A full sample-based KS preparation is required to share its CvM cache.")
  }

  list(
    distance_matrix = ks_prep$distance_matrix,
    rank_matrix = ks_prep$rank_matrix,
    order_matrix = ks_prep$order_matrix,
    sorted_distance_matrix = ks_prep$sorted_distance_matrix,
    rank_linear_index = ks_prep$rank_linear_index,
    empirical_profile = ks_prep$empirical_profile,
    theoretical_profile = ks_prep$theoretical_profile,
    process_matrix = ks_prep$process_matrix,
    statistic = mean(ks_prep$process_matrix^2),
    light = FALSE,
    shared_with_ks = TRUE
  )
}
reestimated_sample_ks_cvm_fusion_eligible <- function(spec,
                                                       ks_prep,
                                                       cvm_prep,
                                                       want_ks,
                                                       want_cvm,
                                                       control = list()) {
  requested <- control$reestimated_fuse_ks_cvm %||% TRUE
  if (length(requested) != 1L || is.na(requested) || !is.logical(requested)) {
    stop("`control$reestimated_fuse_ks_cvm` must be TRUE or FALSE when supplied.")
  }
  if (!isTRUE(requested) || !isTRUE(want_ks) || !isTRUE(want_cvm) ||
      is.null(ks_prep) || is.null(cvm_prep) ||
      !identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances") ||
      !isTRUE(cvm_prep$shared_with_ks) || isTRUE(ks_prep$light) ||
      isTRUE(cvm_prep$light) || is.function(spec$cvm_prepare) ||
      is.function(spec$cvm_bootstrap_stat)) {
    return(FALSE)
  }
  identical(ks_prep$distance_matrix, cvm_prep$distance_matrix) &&
    identical(ks_prep$rank_matrix, cvm_prep$rank_matrix) &&
    identical(ks_prep$rank_linear_index, cvm_prep$rank_linear_index) &&
    identical(ks_prep$empirical_profile, cvm_prep$empirical_profile) &&
    identical(ks_prep$theoretical_profile, cvm_prep$theoretical_profile)
}
compute_weighted_sample_profile_matrix <- function(order_matrix = NULL,
                                                   rank_linear_index = NULL,
                                                   normalized_weights,
                                                   order_list = NULL,
                                                   rank_matrix = NULL) {
  if (is.null(order_matrix) && !is.null(order_list)) {
    n_from_list <- length(order_list)
    order_matrix <- t(vapply(seq_len(n_from_list), function(i) {
      as.integer(order_list[[i]])
    }, integer(length(order_list[[1]]))))
  }
  if (is.null(rank_linear_index) && !is.null(rank_matrix)) {
    n_from_rank <- nrow(rank_matrix)
    row_index_matrix <- matrix(rep.int(seq_len(n_from_rank), n_from_rank), nrow = n_from_rank, ncol = n_from_rank)
    rank_linear_index <- row_index_matrix + (rank_matrix - 1L) * n_from_rank
  }
  if (is.null(order_matrix) || is.null(rank_linear_index)) {
    stop("Weighted sample profile requires either `(order_matrix, rank_linear_index)` or `(order_list, rank_matrix)`.")
  }

  if (identical(distance_profile_backend_current(), "cpp")) {
    return(distance_profile_cpp_call(
      "cpp_dp_weighted_sample_profile_linear",
      order_matrix,
      as.integer(rank_linear_index),
      normalized_weights
    ))
  }

  n <- nrow(order_matrix)
  ordered_weights_matrix <- matrix(
    normalized_weights[order_matrix],
    nrow = n,
    ncol = n
  )
  cumulative_weights_matrix <- ordered_weights_matrix

  if (n >= 2L) {
    for (j in 2:n) {
      cumulative_weights_matrix[, j] <- cumulative_weights_matrix[, j] +
        cumulative_weights_matrix[, j - 1L]
    }
  }

  matrix(cumulative_weights_matrix[rank_linear_index] / n, nrow = n, ncol = n)
}
normalize_ks_block_size <- function(block_size,
                                    n_rows,
                                    default = min(as.integer(n_rows), 64L),
                                    arg_name = "`control$ks_block_size`") {
  resolved <- block_size %||% default
  resolved <- as.integer(resolved)
  if (!is.finite(resolved) || resolved <= 0L) {
    stop(sprintf("%s must be a strictly positive integer.", arg_name))
  }
  min(resolved, as.integer(n_rows))
}
compute_weighted_sample_profile_rows <- function(order_matrix,
                                                 rank_matrix,
                                                 normalized_weights,
                                                 row_indices) {
  order_block <- order_matrix[row_indices, , drop = FALSE]
  rank_block <- rank_matrix[row_indices, , drop = FALSE]
  if (identical(distance_profile_backend_current(), "cpp")) {
    return(distance_profile_cpp_call(
      "cpp_dp_weighted_sample_profile_rows",
      order_block,
      rank_block,
      normalized_weights
    ))
  }
  block_n <- nrow(order_block)
  n_total <- ncol(order_block)

  ordered_weights_matrix <- matrix(
    normalized_weights[order_block],
    nrow = block_n,
    ncol = n_total
  )
  cumulative_weights_matrix <- ordered_weights_matrix
  if (n_total >= 2L) {
    for (j in 2:n_total) {
      cumulative_weights_matrix[, j] <- cumulative_weights_matrix[, j] +
        cumulative_weights_matrix[, j - 1L]
    }
  }

  block_row_index <- matrix(rep.int(seq_len(block_n), n_total), nrow = block_n, ncol = n_total)
  block_linear_index <- block_row_index + (rank_block - 1L) * block_n
  matrix(cumulative_weights_matrix[block_linear_index] / n_total, nrow = block_n, ncol = n_total)
}
compute_theoretical_sample_profile_block <- function(spec,
                                                     normalized_data,
                                                     distance_matrix,
                                                     theta,
                                                     row_indices,
                                                     control = list()) {
  output <- matrix(0, nrow = length(row_indices), ncol = ncol(distance_matrix))

  for (k in seq_along(row_indices)) {
    i <- row_indices[[k]]
    omega_i <- spec_observation_at_normalized(spec, normalized_data, i, control)
    output[k, ] <- as.numeric(spec$profile_eval(omega_i, distance_matrix[i, ], theta, control))
  }

  output
}
compute_ks_sample_stat_blocked <- function(spec,
                                           normalized_data,
                                           ks_prep,
                                           normalized_weights,
                                           scale_factor = 1,
                                           theta_star = NULL,
                                           null_type = c("simple", "composite"),
                                           control = list()) {
  null_type <- match.arg(null_type)
  n <- nrow(ks_prep$distance_matrix)
  block_size <- normalize_ks_block_size(
    block_size = control$ks_block_size %||% NULL,
    n_rows = n
  )
  block_max <- 0

  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    row_indices <- block_start:block_end
    empirical_star_block <- compute_weighted_sample_profile_rows(
      order_matrix = ks_prep$order_matrix,
      rank_matrix = ks_prep$rank_matrix,
      normalized_weights = normalized_weights,
      row_indices = row_indices
    )

    if (identical(null_type, "simple")) {
      process_block <- scale_factor * sqrt(n) * (
        empirical_star_block - ks_prep$empirical_profile[row_indices, , drop = FALSE]
      )
    } else {
      theoretical_star_block <- compute_theoretical_sample_profile_block(
        spec = spec,
        normalized_data = normalized_data,
        distance_matrix = ks_prep$distance_matrix,
        theta = theta_star,
        row_indices = row_indices,
        control = control
      )
      process_block <- scale_factor * sqrt(n) * (
        (empirical_star_block - theoretical_star_block) -
          (ks_prep$empirical_profile[row_indices, , drop = FALSE] -
             ks_prep$theoretical_profile[row_indices, , drop = FALSE])
      )
    }

    block_max <- max(block_max, max(abs(process_block)))
  }

  block_max
}
