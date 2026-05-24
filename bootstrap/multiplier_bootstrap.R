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
if (!exists("make_normal_spec", mode = "function")) {
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

prepare_ks_observed_data <- function(data, spec, theta_hat, ks_grid, control = list()) {
  if (!is.list(ks_grid) || is.null(ks_grid$omega_grid) || is.null(ks_grid$t_grid)) {
    stop("KS requires `ks_grid = list(omega_grid = ..., t_grid = ...)`.")
  }

  omega_grid <- ks_grid$omega_grid
  t_grid <- as.numeric(ks_grid$t_grid)
  if (length(t_grid) == 0) {
    stop("`ks_grid$t_grid` cannot be empty.")
  }

  distance_matrix <- spec$distance_matrix(data, omega_grid, control)
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
    omega_grid = omega_grid,
    t_grid = t_grid,
    distance_matrix = distance_matrix,
    order_matrix = ks_cache$order_matrix,
    sorted_distance_matrix = ks_cache$sorted_distance_matrix,
    threshold_index_matrix = ks_cache$threshold_index_matrix,
    empirical_profile = empirical_profile,
    theoretical_profile = theoretical_profile,
    process_matrix = process_matrix,
    statistic = max(abs(process_matrix))
  )
}

compute_theoretical_sample_profile_matrix <- function(spec,
                                                     data,
                                                     distance_matrix,
                                                     theta,
                                                     control = list()) {
  fast_output <- spec_sample_profile_matrix_eval(
    spec = spec,
    data = data,
    distance_matrix = distance_matrix,
    theta = theta,
    control = control
  )
  if (!is.null(fast_output)) {
    n <- nrow(distance_matrix)
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

prepare_cvm_observed_data <- function(data, spec, theta_hat, control = list()) {
  distance_matrix <- spec$distance_matrix(data, data, control)
  n <- nrow(distance_matrix)

  rank_matrix <- t(vapply(seq_len(n), function(i) {
    as.integer(rank(distance_matrix[i, ], ties.method = "max"))
  }, integer(n)))

  order_list <- lapply(seq_len(n), function(i) {
    order(distance_matrix[i, ])
  })
  order_matrix <- t(vapply(seq_len(n), function(i) {
    as.integer(order(distance_matrix[i, ]))
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

  list(
    distance_matrix = distance_matrix,
    rank_matrix = rank_matrix,
    order_list = order_list,
    order_matrix = order_matrix,
    rank_linear_index = rank_linear_index,
    empirical_profile = empirical_profile,
    theoretical_profile = theoretical_profile,
    process_matrix = process_matrix,
    statistic = mean(process_matrix^2)
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
                                keep_bootstrap_thetas = FALSE) {
  n_reps <- nrow(weight_chunk)
  ks_values <- if (want_ks) numeric(n_reps) else NULL
  cvm_values <- if (want_cvm) numeric(n_reps) else NULL
  theta_values <- if (keep_bootstrap_thetas && identical(null$type, "composite")) {
    vector("list", n_reps)
  } else {
    NULL
  }
  n <- spec_n_obs(spec, data, control)

  for (b in seq_len(n_reps)) {
    normalized_weights <- weight_chunk[b, ]

    theta_star <- NULL
    if (identical(null$type, "composite")) {
      theta_star <- spec$fit_theta(
        data = data,
        weights = normalized_weights,
        null = null,
        control = control
      )
      if (!is.null(theta_values)) {
        theta_values[[b]] <- theta_star
      }
    }

    if (want_ks) {
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

    if (want_cvm) {
      f_star_sample <- compute_weighted_sample_profile_matrix(
        order_matrix = cvm_prep$order_matrix,
        rank_linear_index = cvm_prep$rank_linear_index,
        normalized_weights = normalized_weights
      )

      if (identical(null$type, "simple")) {
        process_star_sample <- scale_factor * sqrt(n) * (f_star_sample - cvm_prep$empirical_profile)
      } else {
        f_theta_star_sample <- compute_theoretical_sample_profile_matrix(
          spec = spec,
          data = data,
          distance_matrix = cvm_prep$distance_matrix,
          theta = theta_star,
          control = control
        )
        process_star_sample <- scale_factor * sqrt(n) * (
          (f_star_sample - f_theta_star_sample) -
            (cvm_prep$empirical_profile - cvm_prep$theoretical_profile)
        )
      }

      cvm_values[b] <- mean(process_star_sample^2)
    }
  }

  list(ks = ks_values, cvm = cvm_values, theta = theta_values)
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
    if (keep_options$observed_process) {
      output$cvm$process_matrix <- cvm_prep$process_matrix
      output$cvm$distance_matrix <- cvm_prep$distance_matrix
      output$cvm$empirical_profile <- cvm_prep$empirical_profile
      output$cvm$theoretical_profile <- cvm_prep$theoretical_profile
    }
  }

  output
}

multiplier_bootstrap_gof <- function(data,
                                     spec,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 999,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list()) {
  validate_model_spec(spec)
  null <- validate_null_object(null)
  statistics <- normalize_requested_statistics(statistics)
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

  if (want_ks && is.null(ks_grid)) {
    stop("KS requires a user-provided `ks_grid`.")
  }

  multiplier_spec <- resolve_multiplier_spec(multipliers)
  scale_factor <- multiplier_spec$mean / multiplier_spec$sd

  start_time <- Sys.time()

  theta_hat <- spec$fit_theta(
    data = data_normalized,
    weights = NULL,
    null = null,
    control = control
  )

  ks_prep <- if (want_ks) {
    prepare_ks_observed_data(
      data = data_normalized,
      spec = spec,
      theta_hat = theta_hat,
      ks_grid = ks_grid,
      control = control
    )
  } else {
    NULL
  }

  cvm_prep <- if (want_cvm) {
    prepare_cvm_observed_data(
      data = data_normalized,
      spec = spec,
      theta_hat = theta_hat,
      control = control
    )
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

  n_cores_effective <- min(n_cores, B)
  chunk_ids <- split(seq_len(B), rep(seq_len(n_cores_effective), length.out = B))
  weight_chunks <- lapply(chunk_ids, function(indices) {
    normalized_multiplier_matrix[indices, , drop = FALSE]
  })

  if (n_cores_effective == 1L) {
    chunk_results <- lapply(weight_chunks, run_bootstrap_chunk,
      spec = spec,
      data = data_normalized,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      want_ks = want_ks,
      want_cvm = want_cvm,
      keep_bootstrap_thetas = keep$bootstrap_thetas
    )
  } else {
    utils_path_worker <- normalizePath(resolve_multiplier_bootstrap_path("utils.R"), winslash = "/", mustWork = TRUE)
    model_specs_path_worker <- normalizePath(resolve_multiplier_bootstrap_path("bootstrap", "model_specs.R"), winslash = "/", mustWork = TRUE)
    cardioid_model_spec_path_worker <- normalizePath(resolve_multiplier_bootstrap_path("bootstrap", "cardioid_model_spec.R"), winslash = "/", mustWork = TRUE)

    cl <- parallel::makeCluster(n_cores_effective)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      c("utils_path_worker", "model_specs_path_worker", "cardioid_model_spec_path_worker"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      source(utils_path_worker)
      source(model_specs_path_worker)
      source(cardioid_model_spec_path_worker)
      NULL
    })

    worker_symbols <- c(
      "spec",
      "data_normalized",
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
      "keep"
    )

    parallel::clusterExport(cl, worker_symbols, envir = environment())

    chunk_results <- parallel::parLapply(cl, weight_chunks, function(chunk) {
      run_bootstrap_chunk(
        weight_chunk = chunk,
        spec = spec,
        data = data_normalized,
        null = null,
        control = control,
        scale_factor = scale_factor,
        ks_prep = ks_prep,
        cvm_prep = cvm_prep,
        want_ks = want_ks,
        want_cvm = want_cvm,
        keep_bootstrap_thetas = keep$bootstrap_thetas
      )
    })
  }

  bootstrap_statistics_internal <- list()
  if (want_ks) {
    bootstrap_statistics_internal$ks <- unlist(lapply(chunk_results, `[[`, "ks"), use.names = FALSE)
  }
  if (want_cvm) {
    bootstrap_statistics_internal$cvm <- unlist(lapply(chunk_results, `[[`, "cvm"), use.names = FALSE)
  }
  bootstrap_theta_internal <- if (keep$bootstrap_thetas && identical(null$type, "composite")) {
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
      engine = "multiplier_bootstrap_gof",
      method = "distance_profiles",
      weighted_mle = isTRUE(spec$weighted_mle),
      elapsed_seconds = elapsed_seconds
    )
  )

  class(result) <- c("multiplier_bootstrap_gof_result", "list")
  result
}

multiplier_bootstrap_normal <- function(data,
                                        null,
                                        statistics = c("ks", "cvm"),
                                        ks_grid = NULL,
                                        B = 999,
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
                                        unknown_param = NULL) {
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
    keep = keep,
    control = control
  )
}

multiplier_bootstrap_vmf <- function(data,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 999,
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
                                     distance_type = c("chordal", "geodesic"),
                                     unknown_param = "xi") {
  distance_type <- match.arg(distance_type)
  spec <- make_vmf_spec(
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
    keep = keep,
    control = control
  )
}

multiplier_bootstrap_hvmf <- function(data,
                                      null,
                                      statistics = c("ks", "cvm"),
                                      ks_grid = NULL,
                                      B = 999,
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
                                      unknown_param = "both") {
  spec <- make_hvmf_spec(unknown_param = unknown_param)

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
    control = control
  )
}
