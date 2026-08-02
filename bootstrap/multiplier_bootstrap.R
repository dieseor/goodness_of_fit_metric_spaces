# Multiplier bootstrap engine for goodness-of-fit tests in metric spaces

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

model_specs_path <- resolve_multiplier_bootstrap_path("bootstrap", "model_specs.R")
if (!exists("make_normal_spec", mode = "function") ||
    !exists("make_small_circle_weighted_mixture2_spec", mode = "function")) {
  source(model_specs_path)
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

# The contiguous-double kernel is the production default.  It implements the
# same statistic as the legacy C++ kernel, but traverses the bootstrap arrays
# in their contiguous layout.  `legacy` remains an explicit reproducibility
# fallback for historical checks.
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
                                                            control = list()) {
  output <- matrix(0, nrow = length(row_indices), ncol = ncol(sorted_distance_matrix))

  for (k in seq_along(row_indices)) {
    i <- row_indices[[k]]
    omega_i <- spec_observation_at_normalized(spec, normalized_data, i, control)
    output[k, ] <- as.numeric(spec$profile_eval(omega_i, sorted_distance_matrix[i, ], theta, control))
  }

  output
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
  block_max <- 0

  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    row_indices <- block_start:block_end
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
      control = control
    )
    process_block <- sqrt(n) * (empirical_block - theoretical_block)
    block_max <- max(block_max, max(abs(process_block)))
  }

  block_max
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
  block_max <- 0
  cvm_sum <- 0

  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    row_indices <- block_start:block_end
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
      control = control
    )
    process_block <- sqrt(n) * (empirical_block - theoretical_block)
    block_max <- max(block_max, max(abs(process_block)))
    cvm_sum <- cvm_sum + sum(process_block^2)
  }

  list(
    ks = block_max,
    cvm = cvm_sum / (n * n)
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
  cvm_sum <- 0

  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    row_indices <- block_start:block_end
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
      control = control
    )
    process_block <- sqrt(n) * (empirical_block - theoretical_block)
    cvm_sum <- cvm_sum + sum(process_block^2)
  }

  cvm_sum / (n * n)
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

run_bootstrap_chunk <- function(weight_chunk,
                                spec,
                                data,
                                null,
                                control,
                                scale_factor,
                                ks_prep = NULL,
                                cvm_prep = NULL,
                                want_ks = FALSE,
                                want_cvm = FALSE,
                                keep_bootstrap_thetas = FALSE,
                                theta_start = NULL,
                                replicate_indices = NULL) {
  scoped_backend <- distance_profile_backend_from_control(control)
  previous_backend <- distance_profile_backend_current()
  if (!identical(scoped_backend, previous_backend)) {
    if (identical(scoped_backend, "cpp")) ensure_distance_profile_cpp_loaded()
    .distance_profile_cpp_state$active_backend <- scoped_backend
    on.exit({
      .distance_profile_cpp_state$active_backend <- previous_backend
    }, add = TRUE)
  }
  n_reps <- nrow(weight_chunk)
  ks_values <- if (want_ks) numeric(n_reps) else NULL
  cvm_values <- if (want_cvm) numeric(n_reps) else NULL
  prep_seconds_total <- 0
  loop_seconds_total <- 0
  theta_values <- if (keep_bootstrap_thetas && identical(null$type, "composite")) {
    vector("list", n_reps)
  } else {
    NULL
  }
  n <- spec_n_obs(spec, data, control)

  for (b in seq_len(n_reps)) {
    replicate_start <- proc.time()[["elapsed"]]
    normalized_weights <- weight_chunk[b, ]
    replicate_index <- if (is.null(replicate_indices)) b else as.integer(replicate_indices[[b]])
    debug_memory_log(
      control,
      sprintf("run_bootstrap_chunk: start replicate %d/%d", replicate_index, n_reps),
      list(
        weight_chunk = weight_chunk,
        normalized_weights = normalized_weights
      )
    )

    theta_star <- NULL
    bootstrap_fit_warnings <- character()
    theta_star_loglik <- NA_real_
    theta_star_convergence <- NA_integer_
    if (identical(null$type, "composite")) {
      bootstrap_control <- control
      if (!is.null(theta_start) && grepl("^jp_", spec$name)) {
        # JP composite bootstrap refits use a warm-started local re-optimization.
        # Together with the logic in jp_mle_s2_weighted(), this keeps the refit
        # on the observed sign branch of psi unless the caller explicitly
        # overrides it. This is a stabilization device for the JP optimizer, not
        # the fully unconstrained composite re-fit.
        bootstrap_control$jp_mle_start_theta <- theta_start
        bootstrap_control$jp_mle_warm_start_only <- TRUE
        bootstrap_control$jp_mle_bootstrap_refit <- TRUE
      } else if (!is.null(theta_start) && grepl("^beta_mixture2_", spec$name)) {
        bootstrap_control$beta_mixture2_start_theta <- theta_start
        bootstrap_control$beta_mixture2_warm_start_only <- TRUE
        bootstrap_control$beta_mixture2_n_starts <- bootstrap_control$beta_mixture2_bootstrap_n_starts %||% 1L
        bootstrap_control$beta_mixture2_optim_control <- bootstrap_control$beta_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^uniform_beta_mixture_", spec$name)) {
        bootstrap_control$uniform_beta_mixture_start_theta <- theta_start
        bootstrap_control$uniform_beta_mixture_warm_start_only <- TRUE
        bootstrap_control$uniform_beta_mixture_n_starts <- bootstrap_control$uniform_beta_mixture_bootstrap_n_starts %||% 1L
        bootstrap_control$uniform_beta_mixture_optim_control <- bootstrap_control$uniform_beta_mixture_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^small_circle_symmetric_mixture2_", spec$name)) {
        bootstrap_control$small_circle_symmetric_mixture2_start_theta <- theta_start
        bootstrap_control$small_circle_symmetric_mixture2_warm_start_only <- TRUE
        bootstrap_control$small_circle_symmetric_mixture2_n_starts <-
          bootstrap_control$small_circle_symmetric_mixture2_bootstrap_n_starts %||% 1L
        bootstrap_control$small_circle_symmetric_mixture2_optim_control <-
          bootstrap_control$small_circle_symmetric_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^small_circle_weighted_mixture2_", spec$name)) {
        bootstrap_control$small_circle_weighted_mixture2_start_theta <- theta_start
        bootstrap_control$small_circle_weighted_mixture2_warm_start_only <- TRUE
        bootstrap_control$small_circle_weighted_mixture2_n_starts <-
          bootstrap_control$small_circle_weighted_mixture2_bootstrap_n_starts %||% 1L
        bootstrap_control$small_circle_weighted_mixture2_optim_control <-
          bootstrap_control$small_circle_weighted_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^axial_truncnorm_mixture2_", spec$name)) {
        bootstrap_control$axial_truncnorm_mixture2_start_theta <- theta_start
        bootstrap_control$axial_truncnorm_mixture2_optim_control <-
          bootstrap_control$axial_truncnorm_mixture2_bootstrap_optim_control %||%
          list(maxit = 80L, reltol = 1e-6)
      } else if (!is.null(theta_start) && grepl("^logitnormal_mixture2_", spec$name)) {
        bootstrap_control$logitnormal_mixture2_start_theta <- theta_start
        bootstrap_control$logitnormal_mixture2_warm_start_only <- TRUE
      } else if (!is.null(theta_start) && is.null(bootstrap_control$jp_mle_start_theta)) {
        bootstrap_control$jp_mle_start_theta <- theta_start
      }
      theta_star <- tryCatch(
        withCallingHandlers(
          spec$fit_theta(
            data = data,
            weights = normalized_weights,
            null = null,
            control = bootstrap_control
          ),
          warning = function(w) {
            bootstrap_fit_warnings <<- c(bootstrap_fit_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          bootstrap_fit_warnings <<- c(
            bootstrap_fit_warnings,
            paste0("Bootstrap MLE error: ", conditionMessage(e))
          )
          NULL
        }
      )
      theta_star_loglik <- as.numeric(theta_star$loglik %||% NA_real_)
      theta_star_convergence <- as.integer(theta_star$opt$convergence %||% NA_integer_)
      if (grepl("^small_circle_symmetric_mixture2_", spec$name) &&
          (is.null(theta_star) ||
             any(!is.finite(as.numeric(c(theta_star$mu, theta_star$kappa, theta_star$nu)))))) {
        bootstrap_fit_warnings <- c(
          bootstrap_fit_warnings,
          "Bootstrap MLE returned non-finite theta_star; falling back to observed theta_hat."
        )
        theta_star <- theta_start
        theta_star_loglik <- as.numeric(theta_start$loglik %||% NA_real_)
        theta_star_convergence <- as.integer(theta_start$opt$convergence %||% NA_integer_)
      }
      if (grepl("^axial_truncnorm_mixture2_", spec$name) &&
          (is.null(theta_star) ||
             any(!is.finite(as.numeric(c(
               theta_star$pi,
               theta_star$kappa1,
               theta_star$nu1,
               theta_star$kappa2,
               theta_star$nu2
             )))))) {
        bootstrap_fit_warnings <- c(
          bootstrap_fit_warnings,
          "Bootstrap axial MLE returned an invalid theta_star; falling back to observed theta_hat."
        )
        theta_star <- theta_start
        theta_star_loglik <- as.numeric(theta_start$loglik %||% NA_real_)
        theta_star_convergence <- as.integer(theta_start$opt$convergence %||% NA_integer_)
      }
      if (!is.null(theta_values)) {
        theta_values[[b]] <- theta_star
      }
    }

    if (want_ks) {
      if (identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
        ks_values[b] <- compute_ks_sample_stat_blocked(
          spec = spec,
          normalized_data = data,
          ks_prep = ks_prep,
          normalized_weights = normalized_weights,
          scale_factor = scale_factor,
          theta_star = if (identical(null$type, "composite")) theta_star else NULL,
          null_type = null$type,
          control = control
        )
      } else {
        f_star_grid <- compute_grid_weighted_profile(
          distance_matrix = ks_prep$distance_matrix,
          t_grid = ks_prep$t_grid,
          normalized_weights = normalized_weights,
          sorted_distance_matrix = ks_prep$sorted_distance_matrix,
          order_matrix = ks_prep$order_matrix,
          threshold_index_matrix = ks_prep$threshold_index_matrix
        )

        if (identical(null$type, "simple")) {
          process_star_grid <- scale_factor * sqrt(n) * (f_star_grid - ks_prep$empirical_profile)
        } else {
          f_theta_star <- compute_theoretical_profile_matrix(
            spec = spec,
            omega_grid = ks_prep$omega_grid,
            t_grid = ks_prep$t_grid,
            theta = theta_star,
            control = control
          )
          process_star_grid <- scale_factor * sqrt(n) * (
            (f_star_grid - f_theta_star) -
              (ks_prep$empirical_profile - ks_prep$theoretical_profile)
          )
        }

        ks_values[b] <- max(abs(process_star_grid))
      }
    }

    if (want_cvm) {
      cvm_control <- utils::modifyList(
        control,
        list(
          small_circle_symmetric_mixture2_bootstrap_replicate_index = replicate_index,
          small_circle_symmetric_mixture2_bootstrap_warnings = bootstrap_fit_warnings,
          small_circle_symmetric_mixture2_bootstrap_loglik = theta_star_loglik,
          small_circle_symmetric_mixture2_bootstrap_convergence = theta_star_convergence
        )
      )
      cvm_stat_fast <- spec_cvm_bootstrap_stat(
        spec = spec,
        data = data,
        normalized_weights = normalized_weights,
        theta_star = theta_star,
        cvm_prep = cvm_prep,
        null = null,
        control = cvm_control,
        scale_factor = scale_factor
      )
      if (!is.null(cvm_stat_fast)) {
        cvm_values[b] <- cvm_stat_fast
      } else {
        f_star_sample <- compute_weighted_sample_profile_matrix(
          order_matrix = cvm_prep$order_matrix,
          rank_linear_index = cvm_prep$rank_linear_index,
          normalized_weights = normalized_weights
        )
        debug_memory_log(control, sprintf("run_bootstrap_chunk: replicate %d after f_star_sample", b), list(f_star_sample = f_star_sample))

        if (identical(null$type, "simple")) {
          process_star_sample <- scale_factor * sqrt(n) * (f_star_sample - cvm_prep$empirical_profile)
        } else {
          f_theta_star_sample <- compute_theoretical_sample_profile_matrix(
            spec = spec,
            data = data,
            distance_matrix = cvm_prep$distance_matrix,
            theta = theta_star,
            control = cvm_control
          )
          debug_memory_log(
            control,
            sprintf("run_bootstrap_chunk: replicate %d after f_theta_star_sample", b),
            list(f_theta_star_sample = f_theta_star_sample)
          )
          process_star_sample <- scale_factor * sqrt(n) * (
            (f_star_sample - f_theta_star_sample) -
              (cvm_prep$empirical_profile - cvm_prep$theoretical_profile)
          )
        }
        debug_memory_log(
          control,
          sprintf("run_bootstrap_chunk: replicate %d after process_star_sample", b),
          list(process_star_sample = process_star_sample)
        )

        cvm_values[b] <- mean(process_star_sample^2)
      }
    }
    loop_seconds_total <- loop_seconds_total + (proc.time()[["elapsed"]] - replicate_start)
  }

  list(
    ks = ks_values,
    cvm = cvm_values,
    theta = theta_values,
    prep_seconds = prep_seconds_total,
    loop_seconds = loop_seconds_total
  )
}

compute_inference_summary <- function(observed_statistics, bootstrap_statistics, alpha) {
  output <- list()

  for (stat_name in names(observed_statistics)) {
    observed_value <- observed_statistics[[stat_name]]
    bootstrap_values <- bootstrap_statistics[[stat_name]]

    critical_value <- as.numeric(stats::quantile(
      bootstrap_values,
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    ))
    p_value <- (1 + sum(bootstrap_values >= observed_value)) / (length(bootstrap_values) + 1)

    output[[stat_name]] <- list(
      observed = observed_value,
      critical_value = critical_value,
      p_value = p_value,
      reject = isTRUE(p_value <= alpha)
    )
  }

  output
}

build_observed_output <- function(theta_hat, ks_prep, cvm_prep, keep_options) {
  output <- list(theta_hat = theta_hat)

  if (!is.null(ks_prep)) {
    output$ks <- list(statistic = ks_prep$statistic)
    if (keep_options$observed_process) {
      output$ks$process_matrix <- ks_prep$process_matrix
      output$ks$empirical_profile <- ks_prep$empirical_profile
      output$ks$theoretical_profile <- ks_prep$theoretical_profile
    }
  }

  if (!is.null(cvm_prep)) {
    output$cvm <- list(statistic = cvm_prep$statistic)
    if (keep_options$observed_process && !is.null(cvm_prep$process_matrix)) {
      output$cvm$process_matrix <- cvm_prep$process_matrix
      output$cvm$distance_matrix <- cvm_prep$distance_matrix
      output$cvm$empirical_profile <- cvm_prep$empirical_profile
      output$cvm$theoretical_profile <- cvm_prep$theoretical_profile
    }
  }

  output
}

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
  loop_start <- proc.time()[["elapsed"]]
  block_results <- if (length(replicate_blocks) == 1L) {
    list(run_fast_block(replicate_blocks[[1L]]))
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

  utils_path_worker <- normalizePath(resolve_multiplier_bootstrap_path("utils.R"), winslash = "/", mustWork = TRUE)
  model_specs_path_worker <- normalizePath(resolve_multiplier_bootstrap_path("bootstrap", "model_specs.R"), winslash = "/", mustWork = TRUE)

  cl <- parallel::makeCluster(n_cores_effective)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    c(
      "utils_path_worker",
      "model_specs_path_worker"
    ),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    source(utils_path_worker)
    source(model_specs_path_worker)
    NULL
  })

  worker_symbols <- c(
    "spec",
    "data",
    "control",
    "ks_prep",
    "cvm_prep",
    "want_ks",
    "want_cvm",
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
    "make_small_circle_spec",
    "normalize_watson_data",
    "normalize_watson_theta",
    "fit_watson_theta",
    "make_watson_spec",
    "normalize_small_circle_symmetric_mixture2_data",
    "normalize_small_circle_symmetric_mixture2_theta",
    "fit_small_circle_symmetric_mixture2_theta",
    "make_small_circle_symmetric_mixture2_spec",
    "normalize_small_circle_weighted_mixture2_data",
    "normalize_small_circle_weighted_mixture2_theta",
    "fit_small_circle_weighted_mixture2_theta",
    "make_small_circle_weighted_mixture2_spec",
    "normalize_axial_truncnorm_mixture2_data",
    "normalize_axial_truncnorm_mixture2_theta",
    "fit_axial_truncnorm_mixture2_theta",
    "make_axial_truncnorm_mixture2_spec"
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
      keep_bootstrap_thetas = keep_bootstrap_thetas,
      theta_start = theta_hat,
      replicate_indices = replicate_index_chunks[[i]]
    )
  })
}

multiplier_bootstrap_gof <- function(data,
                                     spec,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 5000,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     observed_theta_hat = NULL,
                                     bootstrap_method = c("reestimated", "fast_multiplier"),
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list()) {
  validate_model_spec(spec)
  null <- validate_null_object(null)
  statistics <- normalize_requested_statistics(statistics)
  bootstrap_method <- match.arg(bootstrap_method)
  keep <- normalize_keep_options(keep)

  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  if (!is.finite(B) || B <= 0) {
    stop("`B` must be a strictly positive integer.")
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must belong to (0, 1).")
  }
  if (!is.finite(n_cores) || n_cores <= 0) {
    stop("`n_cores` must be a strictly positive integer.")
  }

  want_ks <- "ks" %in% statistics
  want_cvm <- "cvm" %in% statistics

  data_normalized <- spec_normalize_data(spec, data, control)
  n <- spec_n_obs(spec, data_normalized, control)
  use_lightweight_sample_ks_prep <- identical(bootstrap_method, "fast_multiplier") &&
    want_ks &&
    is_sample_unique_distance_ks_grid(ks_grid %||% list()) &&
    !keep$observed_process
  use_lightweight_cvm_prep <- identical(bootstrap_method, "fast_multiplier") &&
    want_cvm &&
    !keep$observed_process

  if (want_ks && is.null(ks_grid)) {
    ks_grid <- make_sample_unique_distance_ks_grid()
    use_lightweight_sample_ks_prep <- identical(bootstrap_method, "fast_multiplier") &&
      is_sample_unique_distance_ks_grid(ks_grid) &&
      !keep$observed_process
  }

  multiplier_spec <- resolve_multiplier_spec(multipliers)
  scale_factor <- multiplier_spec$mean / multiplier_spec$sd

  start_time <- Sys.time()
  common_observed_start <- proc.time()[["elapsed"]]

  theta_hat <- if (is.null(observed_theta_hat)) {
    spec$fit_theta(
      data = data_normalized,
      weights = NULL,
      null = null,
      control = control
    )
  } else {
    observed_theta_hat
  }

  ks_prep <- if (want_ks) {
    prepare_ks_observed_data(
      data = data_normalized,
      spec = spec,
      theta_hat = theta_hat,
      ks_grid = ks_grid,
      control = control,
      light = use_lightweight_sample_ks_prep,
      share_cvm_statistic = use_lightweight_sample_ks_prep &&
        use_lightweight_cvm_prep && want_cvm
    )
  } else {
    NULL
  }

  cvm_prep <- if (want_cvm) {
    if (isTRUE(use_lightweight_cvm_prep) &&
        isTRUE(ks_prep$light) &&
        identical(ks_prep$ks_grid_mode %||% "", "sample_points_unique_distances")) {
      prepare_cvm_observed_data_from_sample_ks(
        data = data_normalized,
        spec = spec,
        theta_hat = theta_hat,
        ks_prep = ks_prep,
        control = control
      )
    } else {
      prepare_cvm_observed_data(
        data = data_normalized,
        spec = spec,
        theta_hat = theta_hat,
        control = control,
        light = use_lightweight_cvm_prep
      )
    }
  } else {
    NULL
  }

  raw_multiplier_matrix <- generate_multiplier_matrix(
    B = B,
    n = n,
    multiplier_spec = multiplier_spec,
    seed = seed
  )
  normalized_multiplier_matrix <- raw_multiplier_matrix / rowMeans(raw_multiplier_matrix)
  common_observed_seconds <- proc.time()[["elapsed"]] - common_observed_start

  n_cores_effective <- min(n_cores, B)

  if (identical(bootstrap_method, "fast_multiplier")) {
    chunk_results <- list(run_fast_multiplier_bootstrap(
      weight_matrix = normalized_multiplier_matrix,
      spec = spec,
      data = data_normalized,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      theta_hat = theta_hat,
      keep_bootstrap_thetas = keep$bootstrap_thetas,
      n_cores = n_cores_effective
    ))
  } else {
    chunk_results <- run_reestimated_bootstrap_chunks(
      weight_matrix = normalized_multiplier_matrix,
      spec = spec,
      data = data_normalized,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      keep_bootstrap_thetas = keep$bootstrap_thetas,
      theta_hat = theta_hat,
      n_cores = n_cores_effective
    )
  }

  bootstrap_statistics_internal <- list()
  if (want_ks) {
    bootstrap_statistics_internal$ks <- unlist(lapply(chunk_results, `[[`, "ks"), use.names = FALSE)
  }
  if (want_cvm) {
    bootstrap_statistics_internal$cvm <- unlist(lapply(chunk_results, `[[`, "cvm"), use.names = FALSE)
  }
  bootstrap_theta_internal <- if (keep$bootstrap_thetas &&
    identical(null$type, "composite") &&
    identical(bootstrap_method, "reestimated")) {
    unlist(lapply(chunk_results, `[[`, "theta"), recursive = FALSE, use.names = FALSE)
  } else {
    NULL
  }

  observed_statistics <- list()
  if (want_ks) {
    observed_statistics$ks <- ks_prep$statistic
  }
  if (want_cvm) {
    observed_statistics$cvm <- cvm_prep$statistic
  }

  inference <- compute_inference_summary(
    observed_statistics = observed_statistics,
    bootstrap_statistics = bootstrap_statistics_internal,
    alpha = alpha
  )

  end_time <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  branch_prep_seconds <- sum(vapply(chunk_results, function(x) as.numeric(x$prep_seconds %||% 0), numeric(1)))
  branch_loop_seconds <- sum(vapply(chunk_results, function(x) as.numeric(x$loop_seconds %||% 0), numeric(1)))
  fallback_to_reestimated <- any(vapply(chunk_results, function(x) isTRUE(x$fallback_to_reestimated), logical(1)))
  effective_bootstrap_method <- chunk_results[[1L]]$effective_bootstrap_method %||% bootstrap_method
  fallback_reason <- chunk_results[[1L]]$fallback_reason %||% NA_character_

  result <- list(
    observed = build_observed_output(
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      keep_options = keep
    ),
    bootstrap = list(
      statistics = if (keep$bootstrap_statistics) bootstrap_statistics_internal else NULL,
      theta_star = bootstrap_theta_internal,
      multiplier = list(
        name = multiplier_spec$name,
        mean = multiplier_spec$mean,
        sd = multiplier_spec$sd
      ),
      B = B
    ),
    inference = inference,
    grid = if (want_ks) ks_grid else NULL,
    diagnostics = list(
      n = n,
      B = B,
      alpha = alpha,
      seed = seed,
      n_cores = n_cores_effective,
      null_type = null$type,
      spec_name = spec$name,
      bootstrap_method = bootstrap_method,
      effective_bootstrap_method = effective_bootstrap_method,
      fallback_to_reestimated = fallback_to_reestimated,
      fallback_reason = fallback_reason,
      engine = "multiplier_bootstrap_gof",
      method = "distance_profiles",
      weighted_mle = isTRUE(spec$weighted_mle),
      lightweight_ks_prep = isTRUE(ks_prep$light),
      lightweight_cvm_prep = isTRUE(cvm_prep$light),
      shared_sample_ks_cvm_cache = isTRUE(cvm_prep$shared_with_ks),
      ks_prep_bytes = if (!is.null(ks_prep)) as.numeric(object.size(ks_prep)) else NA_real_,
      cvm_prep_bytes = if (!is.null(cvm_prep)) as.numeric(object.size(cvm_prep)) else NA_real_,
      derivative_method = chunk_results[[1L]]$derivative_method %||% NA_character_,
      derivative_method_requested =
        chunk_results[[1L]]$derivative_method_requested %||% NA_character_,
      derivative_method_effective =
        chunk_results[[1L]]$derivative_method_effective %||%
          chunk_results[[1L]]$derivative_method %||% NA_character_,
      derivative_method_selection_source =
        chunk_results[[1L]]$derivative_method_selection_source %||% NA_character_,
      derivative_mc_size = chunk_results[[1L]]$derivative_mc_size %||% NA_integer_,
      derivative_mc_seed = chunk_results[[1L]]$derivative_mc_seed %||% NA_integer_,
      quadrature_algorithm =
        chunk_results[[1L]]$quadrature_algorithm %||% NA_character_,
      quadrature_abs_tol =
        chunk_results[[1L]]$quadrature_abs_tol %||% NA_real_,
      quadrature_max_terms =
        chunk_results[[1L]]$quadrature_max_terms %||% NA_integer_,
      quadrature_initial_upper =
        chunk_results[[1L]]$quadrature_initial_upper %||% NA_real_,
      quadrature_max_upper =
        chunk_results[[1L]]$quadrature_max_upper %||% NA_real_,
      quadrature_tail_consecutive =
        chunk_results[[1L]]$quadrature_tail_consecutive %||% NA_integer_,
      quadrature_eigen_rel_tol =
        chunk_results[[1L]]$quadrature_eigen_rel_tol %||% NA_real_,
      quadrature_clip_tol =
        chunk_results[[1L]]$quadrature_clip_tol %||% NA_real_,
      quadrature_center_evaluations =
        chunk_results[[1L]]$quadrature_center_evaluations %||% NA_integer_,
      quadrature_max_terms_used =
        chunk_results[[1L]]$quadrature_max_terms_used %||% NA_integer_,
      quadrature_max_residual_error_estimate =
        chunk_results[[1L]]$quadrature_max_residual_error_estimate %||% NA_real_,
      quadrature_max_condition_number =
        chunk_results[[1L]]$quadrature_max_condition_number %||% NA_real_,
      quadrature_max_propagated_error_estimate =
        chunk_results[[1L]]$quadrature_max_propagated_error_estimate %||% NA_real_,
      quadrature_max_upper_used =
        chunk_results[[1L]]$quadrature_max_upper_used %||% NA_real_,
      quadrature_max_evaluations =
        chunk_results[[1L]]$quadrature_max_evaluations %||% NA_integer_,
      vhat_method = chunk_results[[1L]]$vhat_method %||% NA_character_,
      S_obs_dim = chunk_results[[1L]]$S_obs_dim %||% NA_integer_,
      Psi_aux_dim = chunk_results[[1L]]$Psi_aux_dim %||% NA_integer_,
      D_ks_dim = chunk_results[[1L]]$D_ks_dim %||% NA_integer_,
      D_cvm_dim = chunk_results[[1L]]$D_cvm_dim %||% NA_integer_,
      fast_ks_mode = chunk_results[[1L]]$fast_ks_mode %||% NA_character_,
      fast_cvm_mode = chunk_results[[1L]]$fast_cvm_mode %||% NA_character_,
      sample_correction_cache_bytes = chunk_results[[1L]]$sample_correction_cache_bytes %||% NA_real_,
      shared_sample_correction_cache = chunk_results[[1L]]$shared_sample_correction_cache %||% NA,
      fast_multiplier_backend_requested =
        chunk_results[[1L]]$fast_multiplier_backend_requested %||% NA_character_,
      fast_multiplier_backend_effective =
        chunk_results[[1L]]$fast_multiplier_backend_effective %||% NA_character_,
      fast_multiplier_cpp_kernel_requested =
        chunk_results[[1L]]$fast_multiplier_cpp_kernel_requested %||% NA_character_,
      fast_multiplier_cpp_kernel_effective =
        chunk_results[[1L]]$fast_multiplier_cpp_kernel_effective %||% NA_character_,
      fast_multiplier_fuse_ks_cvm_requested =
        chunk_results[[1L]]$fast_multiplier_fuse_ks_cvm_requested %||% NA,
      fast_multiplier_fuse_ks_cvm_effective =
        chunk_results[[1L]]$fast_multiplier_fuse_ks_cvm_effective %||% NA,
      fast_multiplier_cache_corrections_requested =
        chunk_results[[1L]]$fast_multiplier_cache_corrections_requested %||% NA_character_,
      fast_multiplier_cache_corrections_effective =
        chunk_results[[1L]]$fast_multiplier_cache_corrections_effective %||% NA,
      Vhat_dim = chunk_results[[1L]]$Vhat_dim %||% NA_integer_,
      score_mean_aux = chunk_results[[1L]]$score_mean_aux %||% NA_real_,
      score_mean_aux_norm = chunk_results[[1L]]$score_mean_aux_norm %||% NA_real_,
      Vhat_eigenvalues = chunk_results[[1L]]$Vhat_eigenvalues %||% NA_real_,
      Vhat_rcond = chunk_results[[1L]]$Vhat_rcond %||% NA_real_,
      Vhat_condition_number = chunk_results[[1L]]$Vhat_condition_number %||% NA_real_,
      fast_parameter_summary = chunk_results[[1L]]$fast_parameter_summary %||% NA_real_,
      common_observed_seconds = common_observed_seconds,
      old_prep_seconds = if (identical(bootstrap_method, "reestimated")) branch_prep_seconds else NA_real_,
      old_loop_seconds = if (identical(bootstrap_method, "reestimated")) branch_loop_seconds else NA_real_,
      old_total_seconds = if (identical(bootstrap_method, "reestimated")) common_observed_seconds + branch_prep_seconds + branch_loop_seconds else NA_real_,
      fast_prep_seconds = if (identical(bootstrap_method, "fast_multiplier")) branch_prep_seconds else NA_real_,
      fast_loop_seconds = if (identical(bootstrap_method, "fast_multiplier")) branch_loop_seconds else NA_real_,
      fast_total_seconds = if (identical(bootstrap_method, "fast_multiplier")) common_observed_seconds + branch_prep_seconds + branch_loop_seconds else NA_real_,
      elapsed_seconds = elapsed_seconds
    )
  )

  class(result) <- c("multiplier_bootstrap_gof_result", "list")
  result
}

.multiplier_bootstrap_gof_backend_implementation <- multiplier_bootstrap_gof

multiplier_bootstrap_gof <- function(data,
                                     spec,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 5000,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     observed_theta_hat = NULL,
                                     bootstrap_method = c("reestimated", "fast_multiplier"),
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list(),
                                     distance_profile_backend = c("r", "cpp")) {
  backend <- normalize_distance_profile_backend(distance_profile_backend)
  spec_name <- as.character(spec$name)
  is_jones_pewsey <- length(spec_name) == 1L && !is.na(spec_name) && grepl("^jp_", spec_name)
  if (is_jones_pewsey && identical(backend, "r")) {
    return(.multiplier_bootstrap_gof_backend_implementation(
      data = data,
      spec = spec,
      null = null,
      statistics = statistics,
      ks_grid = ks_grid,
      B = B,
      alpha = alpha,
      multipliers = multipliers,
      n_cores = n_cores,
      seed = seed,
      observed_theta_hat = observed_theta_hat,
      bootstrap_method = bootstrap_method,
      keep = keep,
      control = control
    ))
  }
  if (identical(backend, "cpp")) {
    assert_distance_profile_cpp_spec_available(spec$name)
  }
  control$distance_profile_backend <- backend
  result <- with_distance_profile_backend(
    backend,
    .multiplier_bootstrap_gof_backend_implementation(
      data = data,
      spec = spec,
      null = null,
      statistics = statistics,
      ks_grid = ks_grid,
      B = B,
      alpha = alpha,
      multipliers = multipliers,
      n_cores = n_cores,
      seed = seed,
      observed_theta_hat = observed_theta_hat,
      bootstrap_method = bootstrap_method,
      keep = keep,
      control = control
    )
  )
  result$diagnostics$distance_profile_backend_requested <- backend
  result$diagnostics$distance_profile_backend_effective <- backend
  result
}

multiplier_bootstrap_normal <- function(data,
                                        null,
                                        statistics = c("ks", "cvm"),
                                        ks_grid = NULL,
                                        B = 5000,
                                        alpha = 0.05,
                                        multipliers = NULL,
                                        n_cores = 1,
                                        seed = NULL,
                                        bootstrap_method = c("reestimated", "fast_multiplier"),
                                        keep = list(
                                          observed_process = TRUE,
                                          bootstrap_statistics = TRUE,
                                          bootstrap_thetas = FALSE
                                        ),
                                        control = list(),
                                        unknown_param = NULL,
                                        distance_profile_backend = c("r", "cpp")) {
  spec <- make_normal_spec(unknown_param = unknown_param)
  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_mvnormal <- function(data,
                                          null,
                                          statistics = c("ks", "cvm"),
                                          ks_grid = NULL,
                                          B = 5000,
                                          alpha = 0.05,
                                          multipliers = NULL,
                                          n_cores = 1,
                                          seed = NULL,
                                          bootstrap_method = c("reestimated", "fast_multiplier"),
                                          keep = list(
                                            observed_process = TRUE,
                                            bootstrap_statistics = TRUE,
                                            bootstrap_thetas = FALSE
                                          ),
                                          control = list(),
                                          unknown_param = "both",
                                          distance_profile_backend = c("r", "cpp"),
                                          fast_multiplier_backend = c("cpp", "r"),
                                          fuse_ks_cvm = TRUE,
                                          cache_block_corrections = c(
                                            "auto", "true", "false"
                                          )) {
  fast_multiplier_backend <- normalize_fast_multiplier_backend(
    fast_multiplier_backend
  )
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(fuse_ks_cvm)
  cache_block_corrections <- normalize_fast_multiplier_cache(
    cache_block_corrections
  )
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  legacy_mc_control <- is.null(control$derivative_method) &&
    (!is.null(control$derivative_mc_size) ||
       !is.null(control$derivative_mc_seed))
  if (legacy_mc_control) {
    warning(
      paste(
        "Multivariate-normal fast multiplier: legacy derivative MC controls",
        "were supplied without `derivative_method`. Selecting `score_mc`."
      ),
      call. = FALSE
    )
  }
  requested_method <- tolower(as.character(
    control$derivative_method %||% if (legacy_mc_control) "score_mc" else "auto"
  ))
  effective_method <- if (requested_method %in% c("auto", "deterministic")) {
    "quadrature"
  } else {
    requested_method
  }
  selection_source <- if (!is.null(control$derivative_method)) {
    if (identical(requested_method, "auto")) "explicit_auto" else "explicit"
  } else if (legacy_mc_control) {
    "legacy_mc_controls"
  } else {
    "model_default"
  }
  control$derivative_method <- effective_method
  spec <- make_mvnormal_spec(unknown_param = unknown_param)
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
  result$diagnostics$fast_multiplier_backend_requested <-
    fast_multiplier_backend
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <-
    fuse_ks_cvm
  result$diagnostics$fast_multiplier_cache_corrections_requested <-
    cache_block_corrections
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$derivative_method_selection_source <- selection_source
  if (!identical(
      result$diagnostics$effective_bootstrap_method,
      "fast_multiplier"
  )) {
    result$diagnostics$fast_multiplier_backend_effective <- "r"
    result$diagnostics$fast_multiplier_fuse_ks_cvm_effective <- FALSE
    result$diagnostics$fast_multiplier_cache_corrections_effective <- FALSE
  }
  result
}

multiplier_bootstrap_vmf <- function(data,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 5000,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     bootstrap_method = c("reestimated", "fast_multiplier"),
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list(),
                                     distance_type = c("chordal", "geodesic"),
                                     unknown_param = "xi",
                                     distance_profile_backend = c("r", "cpp"),
                                     fast_multiplier_backend = c("cpp", "r"),
                                     fast_multiplier_cpp_kernel = c(
                                       "contiguous_double", "legacy"
                                     ),
                                     fuse_ks_cvm = TRUE,
                                     cache_block_corrections = c(
                                       "auto", "true", "false"
                                     )) {
  distance_type <- match.arg(distance_type)
  fast_multiplier_backend <- normalize_fast_multiplier_backend(
    if (missing(fast_multiplier_backend)) {
      control$fast_multiplier_backend %||% "cpp"
    } else {
      fast_multiplier_backend
    }
  )
  fast_multiplier_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(
    if (missing(fast_multiplier_cpp_kernel)) {
      control$fast_multiplier_cpp_kernel %||% "contiguous_double"
    } else {
      fast_multiplier_cpp_kernel
    }
  )
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(
    if (missing(fuse_ks_cvm)) {
      control$fast_multiplier_fuse_ks_cvm %||% TRUE
    } else {
      fuse_ks_cvm
    }
  )
  cache_block_corrections <- normalize_fast_multiplier_cache(
    if (missing(cache_block_corrections)) {
      control$fast_multiplier_cache_corrections %||% "auto"
    } else {
      cache_block_corrections
    }
  )
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_cpp_kernel <- fast_multiplier_cpp_kernel
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  spec <- make_vmf_spec(
    distance_type = distance_type,
    unknown_param = unknown_param
  )

  legacy_mc_control <- is.null(control$derivative_method) &&
    (!is.null(control$derivative_mc_size) ||
       !is.null(control$derivative_mc_seed))
  requested_method <- control$derivative_method %||% if (legacy_mc_control) {
    "score_mc"
  } else {
    "quadrature"
  }
  selection_source <- if (!is.null(control$derivative_method)) {
    "explicit"
  } else if (legacy_mc_control) {
    "legacy_mc_controls"
  } else {
    "model_default"
  }
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_selection_source <- selection_source
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$fast_multiplier_backend_requested <-
    fast_multiplier_backend
  result$diagnostics$fast_multiplier_cpp_kernel_requested <-
    fast_multiplier_cpp_kernel
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <-
    fuse_ks_cvm
  result$diagnostics$fast_multiplier_cache_corrections_requested <-
    cache_block_corrections
  if (!identical(
      result$diagnostics$effective_bootstrap_method,
      "fast_multiplier"
  )) {
    result$diagnostics$fast_multiplier_backend_effective <- "r"
    result$diagnostics$fast_multiplier_cpp_kernel_effective <- "not_used"
    result$diagnostics$fast_multiplier_fuse_ks_cvm_effective <- FALSE
    result$diagnostics$fast_multiplier_cache_corrections_effective <- FALSE
  }
  result
}

multiplier_bootstrap_jp <- function(data,
                                    null,
                                    statistics = c("ks", "cvm"),
                                    ks_grid = NULL,
                                    B = 5000,
                                    alpha = 0.05,
                                    multipliers = NULL,
                                    n_cores = 1,
                                    seed = NULL,
                                    bootstrap_method = c("reestimated", "fast_multiplier"),
                                    keep = list(
                                      observed_process = TRUE,
                                      bootstrap_statistics = TRUE,
                                      bootstrap_thetas = FALSE
                                    ),
                                    control = list(),
                                    distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)
  spec <- make_jp_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control
  )
}

multiplier_bootstrap_hvmf <- function(data,
                                      null,
                                      statistics = c("ks", "cvm"),
                                      ks_grid = NULL,
                                      B = 5000,
                                      alpha = 0.05,
                                      multipliers = NULL,
                                      n_cores = 1,
                                      seed = NULL,
                                      bootstrap_method = c("reestimated", "fast_multiplier"),
                                      keep = list(
                                        observed_process = TRUE,
                                        bootstrap_statistics = TRUE,
                                        bootstrap_thetas = FALSE
                                      ),
                                      control = list(),
                                      unknown_param = "both",
                                      distance_profile_backend = c("r", "cpp"),
                                      fast_multiplier_backend = c("cpp", "r"),
                                      fast_multiplier_cpp_kernel = c(
                                        "contiguous_double", "legacy"
                                      ),
                                      fuse_ks_cvm = TRUE,
                                      cache_block_corrections = c(
                                        "auto", "true", "false"
                                      )) {
  fast_multiplier_backend <- normalize_fast_multiplier_backend(
    if (missing(fast_multiplier_backend)) {
      control$fast_multiplier_backend %||% "cpp"
    } else {
      fast_multiplier_backend
    }
  )
  fast_multiplier_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(
    if (missing(fast_multiplier_cpp_kernel)) {
      control$fast_multiplier_cpp_kernel %||% "contiguous_double"
    } else {
      fast_multiplier_cpp_kernel
    }
  )
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(
    if (missing(fuse_ks_cvm)) {
      control$fast_multiplier_fuse_ks_cvm %||% TRUE
    } else {
      fuse_ks_cvm
    }
  )
  cache_block_corrections <- normalize_fast_multiplier_cache(
    if (missing(cache_block_corrections)) {
      control$fast_multiplier_cache_corrections %||% "auto"
    } else {
      cache_block_corrections
    }
  )
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_cpp_kernel <- fast_multiplier_cpp_kernel
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  spec <- make_hvmf_spec(unknown_param = unknown_param)

  legacy_mc_control <- is.null(control$derivative_method) &&
    (!is.null(control$derivative_mc_size) ||
       !is.null(control$derivative_mc_seed))
  requested_method <- control$derivative_method %||% if (legacy_mc_control) {
    "score_mc"
  } else {
    "quadrature"
  }
  selection_source <- if (!is.null(control$derivative_method)) {
    "explicit"
  } else if (legacy_mc_control) {
    "legacy_mc_controls"
  } else {
    "model_default"
  }
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
  result$diagnostics$fast_multiplier_backend_requested <-
    fast_multiplier_backend
  result$diagnostics$fast_multiplier_cpp_kernel_requested <-
    fast_multiplier_cpp_kernel
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <-
    fuse_ks_cvm
  result$diagnostics$fast_multiplier_cache_corrections_requested <-
    cache_block_corrections
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_selection_source <- selection_source
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  if (!identical(
      result$diagnostics$effective_bootstrap_method,
      "fast_multiplier"
  )) {
    result$diagnostics$fast_multiplier_backend_effective <- "r"
    result$diagnostics$fast_multiplier_fuse_ks_cvm_effective <- FALSE
    result$diagnostics$fast_multiplier_cache_corrections_effective <- FALSE
  }
  result
}

multiplier_bootstrap_logistic_gaussian <- function(data,
                                                   null,
                                                   statistics = c("ks", "cvm"),
                                                   ks_grid = NULL,
                                                   B = 5000,
                                                   alpha = 0.05,
                                                   multipliers = NULL,
                                                   n_cores = 1,
                                                   seed = NULL,
                                                   bootstrap_method = c("reestimated", "fast_multiplier"),
                                                   keep = list(
                                                     observed_process = TRUE,
                                                     bootstrap_statistics = TRUE,
                                                     bootstrap_thetas = FALSE
                                                   ),
                                                   control = list(),
                                                   unknown_param = "both",
                                                   distance_profile_backend = c("r", "cpp")) {
  legacy_mc_control <- is.null(control$derivative_method) &&
    (!is.null(control$derivative_mc_size) ||
       !is.null(control$derivative_mc_seed))
  if (legacy_mc_control) {
    warning(
      paste(
        "Logistic-Gaussian fast multiplier: legacy derivative MC controls",
        "were supplied without `derivative_method`. Selecting `score_mc`."
      ),
      call. = FALSE
    )
  }
  requested_method <- tolower(as.character(
    control$derivative_method %||% if (legacy_mc_control) "score_mc" else "auto"
  ))
  effective_method <- if (requested_method %in% c("auto", "deterministic")) {
    "quadrature"
  } else {
    requested_method
  }
  selection_source <- if (!is.null(control$derivative_method)) {
    if (identical(requested_method, "auto")) "explicit_auto" else "explicit"
  } else if (legacy_mc_control) {
    "legacy_mc_controls"
  } else {
    "model_default"
  }
  control$derivative_method <- effective_method
  spec <- make_logistic_gaussian_spec(unknown_param = unknown_param)

  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$derivative_method_selection_source <- selection_source
  result
}

multiplier_bootstrap_beta_mixture2 <- function(data,
                                                          null,
                                                          statistics = c("ks", "cvm"),
                                                          ks_grid = NULL,
                                                          B = 5000,
                                                          alpha = 0.05,
                                                          multipliers = NULL,
                                                          n_cores = 1,
                                                          seed = NULL,
                                                          bootstrap_method = c("reestimated", "fast_multiplier"),
                                                          keep = list(
                                                            observed_process = TRUE,
                                                            bootstrap_statistics = TRUE,
                                                            bootstrap_thetas = FALSE
                                                          ),
                                                          control = list(),
                                                          distance_type = c("chordal", "geodesic"),
                                                          distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_beta_mixture2_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_uniform_beta_mixture <- function(data,
                                                      null,
                                                      statistics = c("ks", "cvm"),
                                                      ks_grid = NULL,
                                                      B = 5000,
                                                      alpha = 0.05,
                                                      multipliers = NULL,
                                                      n_cores = 1,
                                                      seed = NULL,
                                                      bootstrap_method = c("reestimated", "fast_multiplier"),
                                                      keep = list(
                                                        observed_process = TRUE,
                                                        bootstrap_statistics = TRUE,
                                                        bootstrap_thetas = FALSE
                                                      ),
                                                      control = list(),
                                                      distance_type = c("chordal", "geodesic"),
                                                      distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_uniform_beta_mixture_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_logitnormal_mixture2 <- function(data,
                                                                 null,
                                                                 statistics = c("ks", "cvm"),
                                                                 ks_grid = NULL,
                                                                 B = 5000,
                                                                 alpha = 0.05,
                                                                 multipliers = NULL,
                                                                 n_cores = 1,
                                                                 seed = NULL,
                                                                 bootstrap_method = c("reestimated", "fast_multiplier"),
                                                                 keep = list(
                                                                   observed_process = TRUE,
                                                                   bootstrap_statistics = TRUE,
                                                                   bootstrap_thetas = FALSE
                                                                   ),
                                                                   control = list(),
                                                                   distance_type = c("chordal", "geodesic"),
                                                                   distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_logitnormal_mixture2_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_cardioid <- function(data,
                                          null,
                                          k,
                                          statistics = c("ks", "cvm"),
                                          ks_grid = NULL,
                                          B = 5000,
                                          alpha = 0.05,
                                          multipliers = NULL,
                                          n_cores = 1,
                                          seed = NULL,
                                          bootstrap_method = c("reestimated", "fast_multiplier"),
                                          keep = list(
                                            observed_process = TRUE,
                                            bootstrap_statistics = TRUE,
                                            bootstrap_thetas = FALSE
                                          ),
                                          control = list(),
                                          distance_type = c("chordal", "geodesic"),
                                          unknown_param = "both",
                                          distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_cardioid_spec(
    k = as.integer(k),
    distance_type = distance_type,
    unknown_param = unknown_param
  )

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_spherical_cauchy <- function(data,
                                                  null,
                                                  statistics = c("ks", "cvm"),
                                                  ks_grid = NULL,
                                                  B = 5000,
                                                  alpha = 0.05,
                                                  multipliers = NULL,
                                                  n_cores = 1,
                                                  seed = NULL,
                                                  bootstrap_method = c("reestimated", "fast_multiplier"),
                                                  keep = list(
                                                    observed_process = TRUE,
                                                    bootstrap_statistics = TRUE,
                                                    bootstrap_thetas = FALSE
                                                    ),
                                                    control = list(),
                                                    distance_type = c("chordal", "geodesic"),
                                                    distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_spherical_cauchy_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_small_circle <- function(data,
                                              null,
                                              statistics = c("ks", "cvm"),
                                              ks_grid = NULL,
                                              B = 5000,
                                              alpha = 0.05,
                                              multipliers = NULL,
                                              n_cores = 1,
                                              seed = NULL,
                                              bootstrap_method = c("reestimated", "fast_multiplier"),
                                              keep = list(
                                                observed_process = TRUE,
                                                bootstrap_statistics = TRUE,
                                                bootstrap_thetas = FALSE
                                                ),
                                                control = list(),
                                                distance_type = c("chordal", "geodesic"),
                                                distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_small_circle_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_watson <- function(data,
                                        null,
                                        statistics = c("ks", "cvm"),
                                        ks_grid = NULL,
                                        B = 5000,
                                        alpha = 0.05,
                                        multipliers = NULL,
                                        n_cores = 1,
                                        seed = NULL,
                                        bootstrap_method = c("reestimated", "fast_multiplier"),
                                        keep = list(
                                          observed_process = TRUE,
                                          bootstrap_statistics = TRUE,
                                          bootstrap_thetas = FALSE
                                          ),
                                          control = list(),
                                          distance_type = c("chordal", "geodesic"),
                                          distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  multiplier_bootstrap_gof(
    data = data,
    spec = make_watson_spec(distance_type = distance_type),
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_small_circle_symmetric_mixture2 <- function(data,
                                                                 null,
                                                                 statistics = c("ks", "cvm"),
                                                                 ks_grid = NULL,
                                                                 B = 5000,
                                                                 alpha = 0.05,
                                                                 multipliers = NULL,
                                                                 n_cores = 1,
                                                                 seed = NULL,
                                                                 bootstrap_method = c("reestimated", "fast_multiplier"),
                                                                 keep = list(
                                                                   observed_process = TRUE,
                                                                   bootstrap_statistics = TRUE,
                                                                   bootstrap_thetas = FALSE
                                                                   ),
                                                                   control = list(),
                                                                   distance_type = c("chordal", "geodesic"),
                                                                   distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_small_circle_symmetric_mixture2_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_small_circle_weighted_mixture2 <- function(data,
                                                                 null,
                                                                 statistics = c("ks", "cvm"),
                                                                ks_grid = NULL,
                                                                B = 5000,
                                                                alpha = 0.05,
                                                                multipliers = NULL,
                                                                n_cores = 1,
                                                                seed = NULL,
                                                                bootstrap_method = c("reestimated", "fast_multiplier"),
                                                                keep = list(
                                                                  observed_process = TRUE,
                                                                  bootstrap_statistics = TRUE,
                                                                  bootstrap_thetas = FALSE
                                                                  ),
                                                                  control = list(),
                                                                  distance_type = c("chordal", "geodesic"),
                                                                  distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_small_circle_weighted_mixture2_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}

multiplier_bootstrap_axial_truncnorm_mixture2 <- function(data,
                                                           null,
                                                           statistics = c("ks", "cvm"),
                                                           ks_grid = NULL,
                                                           B = 5000,
                                                           alpha = 0.05,
                                                           multipliers = NULL,
                                                           n_cores = 1,
                                                           seed = NULL,
                                                           keep = list(
                                                             observed_process = TRUE,
                                                             bootstrap_statistics = TRUE,
                                                             bootstrap_thetas = FALSE
                                                           ),
                                                           control = list(),
                                                           distance_profile_backend = c("r", "cpp")) {
  spec <- make_axial_truncnorm_mixture2_spec(distance_type = "euclidean")
  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}
