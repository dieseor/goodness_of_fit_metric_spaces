# Infrastructure for model-specific multiplier bootstrap adapters

resolve_bootstrap_path <- function(...) {
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

safe_source_text <- function(path, envir = parent.frame()) {
  expr <- parse(text = readLines(path, warn = FALSE))
  for (node in expr) {
    eval(node, envir = envir)
  }
  invisible(TRUE)
}

utils_path_model_specs <- resolve_bootstrap_path("utils.R")
if (!exists("theoretical_distance_profile_normal", mode = "function")) {
  safe_source_text(utils_path_model_specs, envir = environment())
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) {
    rhs
  } else {
    lhs
  }
}

new_model_spec <- function(name,
                           fit_theta,
                           distance_matrix,
                           profile_eval,
                           normalize_data = NULL,
                           n_obs = NULL,
                           observation_at = NULL,
                           extras = list()) {
  spec <- c(
    list(
      name = name,
      fit_theta = fit_theta,
      distance_matrix = distance_matrix,
      profile_eval = profile_eval,
      normalize_data = normalize_data,
      n_obs = n_obs,
      observation_at = observation_at
    ),
    extras
  )

  class(spec) <- c("gof_model_spec", "list")
  validate_model_spec(spec)
  spec
}

validate_model_spec <- function(spec) {
  if (!is.list(spec)) {
    stop("`spec` must be a list.")
  }

  required_fields <- c("name", "fit_theta", "distance_matrix", "profile_eval")
  missing_fields <- required_fields[!vapply(required_fields, function(field) {
    !is.null(spec[[field]])
  }, logical(1))]

  if (length(missing_fields) > 0) {
    stop(sprintf(
      "Model spec is missing required fields: %s",
      paste(missing_fields, collapse = ", ")
    ))
  }

  function_fields <- c("fit_theta", "distance_matrix", "profile_eval")
  invalid_functions <- function_fields[!vapply(function_fields, function(field) {
    is.function(spec[[field]])
  }, logical(1))]

  if (length(invalid_functions) > 0) {
    stop(sprintf(
      "Model spec fields must be functions: %s",
      paste(invalid_functions, collapse = ", ")
    ))
  }

  invisible(spec)
}

spec_normalize_data <- function(spec, data, control = list()) {
  if (is.function(spec$normalize_data)) {
    return(spec$normalize_data(data, control))
  }
  data
}

spec_n_obs_normalized <- function(spec, normalized_data, control = list()) {
  if (is.function(spec$n_obs)) {
    return(spec$n_obs(normalized_data, control))
  }

  if (is.matrix(normalized_data) || is.data.frame(normalized_data)) {
    return(nrow(normalized_data))
  }

  length(normalized_data)
}

spec_n_obs <- function(spec, data, control = list()) {
  normalized_data <- spec_normalize_data(spec, data, control)
  spec_n_obs_normalized(spec, normalized_data, control)
}

spec_observation_at_normalized <- function(spec, normalized_data, idx, control = list()) {
  if (is.function(spec$observation_at)) {
    return(spec$observation_at(normalized_data, idx, control))
  }

  if (is.matrix(normalized_data) || is.data.frame(normalized_data)) {
    return(as.numeric(normalized_data[idx, , drop = TRUE]))
  }

  normalized_data[[idx]]
}

spec_observation_at <- function(spec, data, idx, control = list()) {
  normalized_data <- spec_normalize_data(spec, data, control)
  spec_observation_at_normalized(spec, normalized_data, idx, control)
}

spec_profile_matrix_eval <- function(spec,
                                     omega_grid,
                                     t_grid,
                                     theta,
                                     control = list()) {
  if (!is.function(spec$profile_matrix_eval)) {
    return(NULL)
  }

  spec$profile_matrix_eval(omega_grid, t_grid, theta, control)
}

spec_sample_profile_matrix_eval <- function(spec,
                                            data,
                                            distance_matrix,
                                            theta,
                                            control = list()) {
  if (!is.function(spec$sample_profile_matrix_eval)) {
    return(NULL)
  }

  spec$sample_profile_matrix_eval(data, distance_matrix, theta, control)
}

spec_cvm_prepare <- function(spec,
                             data,
                             theta_hat,
                             control = list()) {
  if (!is.function(spec$cvm_prepare)) {
    return(NULL)
  }

  spec$cvm_prepare(data, theta_hat, control)
}

spec_fast_multiplier_prepare <- function(spec,
                                         data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
  if (!is.function(spec$fast_multiplier_prepare)) {
    return(NULL)
  }

  spec$fast_multiplier_prepare(
    data = data,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control
  )
}

fast_multiplier_numeric_jacobian <- function(fun,
                                             x0,
                                             step = NULL) {
  x0 <- as.numeric(x0)
  p <- length(x0)
  f0 <- as.numeric(fun(x0))
  m <- length(f0)

  if (is.null(step)) {
    step <- pmax(1e-6, abs(x0) * 1e-4)
  }
  step <- as.numeric(step)
  if (length(step) == 1L) {
    step <- rep.int(step, p)
  }
  if (length(step) != p || any(!is.finite(step)) || any(step <= 0)) {
    stop("`step` must be a strictly positive finite scalar or vector of length `length(x0)`.")
  }

  jacobian <- matrix(0, nrow = m, ncol = p)
  for (j in seq_len(p)) {
    delta <- rep.int(0, p)
    delta[[j]] <- step[[j]]
    f_plus <- as.numeric(fun(x0 + delta))
    f_minus <- as.numeric(fun(x0 - delta))
    if (length(f_plus) != m || length(f_minus) != m) {
      stop("The Jacobian target returned values of inconsistent length.")
    }
    jacobian[, j] <- (f_plus - f_minus) / (2 * step[[j]])
  }

  jacobian
}

fast_multiplier_parse_derivative_control <- function(control = list()) {
  derivative_method <- tolower(as.character(control$derivative_method %||% "score_mc"))
  if (!identical(derivative_method, "score_mc")) {
    stop("The fast multiplier bootstrap currently supports only `control$derivative_method = 'score_mc'`.")
  }

  derivative_mc_size <- as.integer(control$derivative_mc_size %||% 1000L)
  if (!is.finite(derivative_mc_size) || derivative_mc_size <= 0L) {
    stop("`control$derivative_mc_size` must be a strictly positive integer.")
  }

  derivative_mc_seed <- control$derivative_mc_seed %||% control$seed %||% NULL

  list(
    derivative_method = derivative_method,
    derivative_mc_size = derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_mc_seed)) NULL else as.integer(derivative_mc_seed)
  )
}

fast_multiplier_parse_vhat_control <- function(control = list()) {
  vhat_method <- tolower(as.character(control$fast_multiplier_vhat_method %||% "score_outer_product"))
  supported_methods <- c("score_outer_product", "numeric_jacobian")
  if (!vhat_method %in% supported_methods) {
    stop(sprintf(
      "`control$fast_multiplier_vhat_method` must be one of %s.",
      paste(sprintf("'%s'", supported_methods), collapse = ", ")
    ))
  }

  rcond_tol <- as.numeric(control$fast_multiplier_vhat_rcond_tol %||% 1e-12)
  if (!is.finite(rcond_tol) || rcond_tol <= 0) {
    stop("`control$fast_multiplier_vhat_rcond_tol` must be a strictly positive finite number.")
  }

  list(
    vhat_method = vhat_method,
    rcond_tol = rcond_tol
  )
}

fast_multiplier_vhat_diagnostics <- function(S_obs,
                                             Psi_aux,
                                             Vhat,
                                             par0) {
  S_obs <- as.matrix(S_obs)
  Psi_aux <- as.matrix(Psi_aux)
  Vhat <- as.matrix(Vhat)
  score_mean_aux <- colMeans(Psi_aux)
  eigenvalues <- tryCatch(
    eigen((Vhat + t(Vhat)) / 2, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep.int(NA_real_, ncol(Vhat))
  )
  rcond_value <- tryCatch(rcond(Vhat), error = function(e) NA_real_)
  condition_number <- if (is.finite(rcond_value) && rcond_value > 0) {
    1 / rcond_value
  } else {
    Inf
  }

  list(
    S_obs_dim = dim(S_obs),
    Psi_aux_dim = dim(Psi_aux),
    Vhat_dim = dim(Vhat),
    score_mean_aux = score_mean_aux,
    score_mean_aux_norm = sqrt(sum(score_mean_aux^2)),
    Vhat_eigenvalues = as.numeric(eigenvalues),
    Vhat_rcond = as.numeric(rcond_value),
    Vhat_condition_number = as.numeric(condition_number),
    par0 = as.numeric(par0)
  )
}

fast_multiplier_validate_vhat <- function(Vhat,
                                          diagnostics,
                                          label = "fast multiplier preparation",
                                          rcond_tol = 1e-12) {
  if (any(!is.finite(Vhat))) {
    stop(sprintf(
      "The %s produced a non-finite `Vhat`. par0 = %s",
      label,
      paste(signif(diagnostics$par0, 8), collapse = ", ")
    ))
  }

  if (any(!is.finite(diagnostics$Vhat_eigenvalues))) {
    stop(sprintf(
      "The %s produced non-finite eigenvalues for `Vhat`. par0 = %s",
      label,
      paste(signif(diagnostics$par0, 8), collapse = ", ")
    ))
  }

  min_eig <- min(diagnostics$Vhat_eigenvalues)
  if (min_eig <= 0) {
    stop(sprintf(
      paste(
        "The %s produced a singular or indefinite `Vhat`.",
        "min_eigenvalue = %.6e, rcond = %.6e, condition_number = %.6e,",
        "score_mean_norm = %.6e, par0 = %s"
      ),
      label,
      min_eig,
      diagnostics$Vhat_rcond,
      diagnostics$Vhat_condition_number,
      diagnostics$score_mean_aux_norm,
      paste(signif(diagnostics$par0, 8), collapse = ", ")
    ))
  }

  if (!is.finite(diagnostics$Vhat_rcond) || diagnostics$Vhat_rcond <= rcond_tol) {
    stop(sprintf(
      paste(
        "The %s produced an ill-conditioned `Vhat`.",
        "rcond = %.6e <= %.6e, min_eigenvalue = %.6e, condition_number = %.6e,",
        "score_mean_norm = %.6e, par0 = %s"
      ),
      label,
      diagnostics$Vhat_rcond,
      rcond_tol,
      min_eig,
      diagnostics$Vhat_condition_number,
      diagnostics$score_mean_aux_norm,
      paste(signif(diagnostics$par0, 8), collapse = ", ")
    ))
  }
}

fast_multiplier_solve_vhat <- function(Vhat,
                                       rhs,
                                       label = "fast multiplier correction") {
  tryCatch(
    solve(Vhat, rhs),
    error = function(e) {
      stop(sprintf(
        "Could not invert `Vhat` while computing %s: %s",
        label,
        conditionMessage(e)
      ))
    }
  )
}

fast_multiplier_vech <- function(mat) {
  mat <- as.matrix(mat)
  mat[lower.tri(mat, diag = TRUE)]
}

fast_multiplier_ivech <- function(values, dim_size) {
  dim_size <- as.integer(dim_size)
  values <- as.numeric(values)
  expected <- dim_size * (dim_size + 1L) / 2L
  if (length(values) != expected) {
    stop("`values` has incompatible length for the requested symmetric dimension.")
  }

  out <- matrix(0, nrow = dim_size, ncol = dim_size)
  out[lower.tri(out, diag = TRUE)] <- values
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  out
}

fast_multiplier_sym_score_to_vech <- function(mat) {
  mat <- as.matrix(mat)
  idx <- which(lower.tri(mat, diag = TRUE), arr.ind = TRUE)
  out <- numeric(nrow(idx))
  for (i in seq_len(nrow(idx))) {
    row_idx <- idx[i, 1L]
    col_idx <- idx[i, 2L]
    factor <- if (row_idx == col_idx) 1 else 2
    out[[i]] <- factor * mat[row_idx, col_idx]
  }
  out
}

fast_multiplier_sphere_chart <- function(mu) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  list(
    mu0 = mu,
    basis = jp_orthonormal_complement(mu),
    ambient_dim = length(mu)
  )
}

fast_multiplier_sphere_chart_map <- function(chart,
                                             a) {
  a <- as.numeric(a)
  if (length(a) != ncol(chart$basis)) {
    stop("Sphere chart coordinates have incompatible dimension.")
  }

  u <- chart$mu0 + drop(chart$basis %*% a)
  norm_u <- sqrt(sum(u^2))
  if (!is.finite(norm_u) || norm_u <= 0) {
    stop("The local spherical chart produced an invalid nonpositive norm.")
  }

  list(
    mu = u / norm_u,
    u = u,
    norm_u = norm_u
  )
}

fast_multiplier_sphere_chart_jacobian <- function(chart,
                                                  a) {
  mapped <- fast_multiplier_sphere_chart_map(chart, a)
  mu <- mapped$mu
  projector <- diag(length(mu)) - tcrossprod(mu)
  (projector %*% chart$basis) / mapped$norm_u
}

fast_multiplier_compute_D_ks <- function(spec,
                                         aux_sample,
                                         Psi_aux,
                                         ks_prep,
                                         control = list()) {
  if (is.null(ks_prep)) {
    return(NULL)
  }

  aux_distance_matrix <- spec$distance_matrix(aux_sample, ks_prep$omega_grid, control)
  if (identical(ks_prep$ks_grid_mode %||% "fixed", "sample_points_unique_distances")) {
    n_aux <- nrow(aux_distance_matrix)
    n_omega <- ncol(aux_distance_matrix)
    aux_order_matrix <- t(vapply(seq_len(n_omega), function(j) {
      as.integer(order(aux_distance_matrix[, j]))
    }, integer(n_aux)))
    aux_sorted_distance_matrix <- matrix(0, nrow = n_omega, ncol = n_aux)
    for (j in seq_len(n_omega)) {
      aux_sorted_distance_matrix[j, ] <- aux_distance_matrix[aux_order_matrix[j, ], j]
    }
    return(list(
      mode = "sample_points_unique_distances",
      aux_order_matrix = aux_order_matrix,
      aux_sorted_distance_matrix = aux_sorted_distance_matrix
    ))
  }

  n_aux <- nrow(aux_distance_matrix)
  threshold_matrix <- matrix(
    rep(as.numeric(ks_prep$t_grid), each = n_aux),
    nrow = n_aux,
    ncol = length(ks_prep$t_grid)
  )
  indicator_blocks <- lapply(seq_len(ncol(aux_distance_matrix)), function(k) {
    distance_block <- matrix(
      rep.int(aux_distance_matrix[, k], length(ks_prep$t_grid)),
      nrow = n_aux,
      ncol = length(ks_prep$t_grid),
      byrow = FALSE
    )
    (distance_block <= threshold_matrix) * 1
  })

  crossprod(do.call(cbind, indicator_blocks), Psi_aux) / nrow(Psi_aux)
}

fast_multiplier_compute_D_cvm <- function(spec,
                                          aux_sample,
                                          Psi_aux,
                                          data,
                                          cvm_prep,
                                          control = list()) {
  if (is.null(cvm_prep)) {
    return(NULL)
  }

  data_normalized <- spec_normalize_data(spec, data, control)
  aux_distance_matrix <- spec$distance_matrix(aux_sample, data_normalized, control)
  if (isTRUE(cvm_prep$light) &&
      !is.null(cvm_prep$order_matrix) &&
      !is.null(cvm_prep$sorted_distance_matrix)) {
    n_aux <- nrow(aux_distance_matrix)
    n_centers <- ncol(aux_distance_matrix)
    aux_order_matrix <- t(vapply(seq_len(n_centers), function(j) {
      as.integer(order(aux_distance_matrix[, j]))
    }, integer(n_aux)))
    aux_sorted_distance_matrix <- matrix(0, nrow = n_centers, ncol = n_aux)
    for (j in seq_len(n_centers)) {
      aux_sorted_distance_matrix[j, ] <- aux_distance_matrix[aux_order_matrix[j, ], j]
    }
    return(list(
      mode = "sample_points_unique_distances_sorted_rows",
      aux_order_matrix = aux_order_matrix,
      aux_sorted_distance_matrix = aux_sorted_distance_matrix
    ))
  }

  observed_distance_matrix <- cvm_prep$distance_matrix
  if (is.null(observed_distance_matrix)) {
    observed_distance_matrix <- spec$distance_matrix(data_normalized, data_normalized, control)
  }
  observed_distance_matrix <- as.matrix(observed_distance_matrix)
  n <- nrow(observed_distance_matrix)
  p <- ncol(Psi_aux)
  default_block_size <- max(1L, floor(as.integer(n) / 20L))
  block_size <- as.integer(control$fast_multiplier_cvm_block_size %||% default_block_size)
  if (!is.finite(block_size) || block_size <= 0L) {
    stop("`control$fast_multiplier_cvm_block_size` must be a strictly positive integer.")
  }

  D_cvm <- matrix(0, nrow = n * n, ncol = p)
  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    block_rows <- block_end - block_start + 1L
    D_block <- matrix(0, nrow = block_rows * n, ncol = p)

    for (offset in seq_len(block_rows)) {
      center_idx <- block_start + offset - 1L
      indicator_block <- (
        matrix(
          rep.int(aux_distance_matrix[, center_idx], n),
          nrow = nrow(aux_distance_matrix),
          ncol = n,
          byrow = FALSE
        ) <= matrix(
          observed_distance_matrix[center_idx, ],
          nrow = nrow(aux_distance_matrix),
          ncol = n,
          byrow = TRUE
        )
      ) * 1
      D_block[((offset - 1L) * n + 1L):(offset * n), ] <- crossprod(indicator_block, Psi_aux) / nrow(Psi_aux)
    }

    idx_flat <- ((block_start - 1L) * n + 1L):(block_end * n)
    D_cvm[idx_flat, ] <- D_block
  }

  D_cvm
}

fast_multiplier_reuse_sample_ks_derivative_for_cvm <- function(D_ks,
                                                                cvm_prep) {
  if (!isTRUE(cvm_prep$shared_with_ks) ||
      !is.list(D_ks) ||
      !identical(D_ks$mode %||% "", "sample_points_unique_distances")) {
    return(NULL)
  }

  list(
    mode = "sample_points_unique_distances_sorted_rows",
    aux_order_matrix = D_ks$aux_order_matrix,
    aux_sorted_distance_matrix = D_ks$aux_sorted_distance_matrix,
    shared_with_ks = TRUE
  )
}

prepare_fast_multiplier_score_model <- function(spec,
                                                data,
                                                theta_hat,
                                                ks_prep = NULL,
                                                cvm_prep = NULL,
                                                control = list(),
                                                par0,
                                                score_matrix_fn,
                                                sample_fn,
                                                vhat_fn = NULL) {
  derivative_control <- fast_multiplier_parse_derivative_control(control)
  vhat_control <- fast_multiplier_parse_vhat_control(control)
  allow_invalid_vhat_diagnostics <- isTRUE(control$fast_multiplier_allow_singular_vhat_diagnostics %||% FALSE)
  data_normalized <- spec_normalize_data(spec, data, control)
  par0 <- as.numeric(par0)

  S_obs <- as.matrix(score_matrix_fn(data_normalized, par0))
  if (nrow(S_obs) != spec_n_obs(spec, data_normalized, control)) {
    stop("The fast multiplier observed score matrix has incompatible dimensions.")
  }

  if (!is.null(derivative_control$derivative_mc_seed)) {
    set.seed(derivative_control$derivative_mc_seed)
  }
  aux_sample <- sample_fn(derivative_control$derivative_mc_size, par0)
  Psi_aux <- as.matrix(score_matrix_fn(aux_sample, par0))
  if (ncol(S_obs) != ncol(Psi_aux)) {
    stop("The fast multiplier observed and auxiliary score matrices have incompatible column dimensions.")
  }

  Vhat <- if (is.function(vhat_fn)) {
    as.matrix(vhat_fn(
      data = data_normalized,
      par0 = par0,
      S_obs = S_obs,
      aux_sample = aux_sample,
      Psi_aux = Psi_aux
    ))
  } else if (identical(vhat_control$vhat_method, "numeric_jacobian")) {
    -fast_multiplier_numeric_jacobian(
      fun = function(par) {
        colMeans(score_matrix_fn(data_normalized, par))
      },
      x0 = par0
    )
  } else {
    crossprod(Psi_aux) / nrow(Psi_aux)
  }
  if (!identical(dim(Vhat), c(ncol(S_obs), ncol(S_obs)))) {
    stop("The fast multiplier preparation produced a `Vhat` with incompatible dimensions.")
  }

  vhat_diagnostics <- fast_multiplier_vhat_diagnostics(
    S_obs = S_obs,
    Psi_aux = Psi_aux,
    Vhat = Vhat,
    par0 = par0
  )
  if (!allow_invalid_vhat_diagnostics) {
    fast_multiplier_validate_vhat(
      Vhat = Vhat,
      diagnostics = vhat_diagnostics,
      label = sprintf("fast multiplier preparation for model '%s'", spec$name),
      rcond_tol = vhat_control$rcond_tol
    )
  }
  observed_cvm_distance_matrix <- NULL
  if (!is.null(cvm_prep) && !isTRUE(cvm_prep$light)) {
    observed_cvm_distance_matrix <- cvm_prep$distance_matrix
    if (is.null(observed_cvm_distance_matrix)) {
      observed_cvm_distance_matrix <- spec$distance_matrix(data_normalized, data_normalized, control)
    }
    observed_cvm_distance_matrix <- as.matrix(observed_cvm_distance_matrix)
  }

  D_ks <- fast_multiplier_compute_D_ks(
    spec = spec,
    aux_sample = aux_sample,
    Psi_aux = Psi_aux,
    ks_prep = ks_prep,
    control = control
  )
  D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
  if (is.null(D_cvm)) {
    D_cvm <- fast_multiplier_compute_D_cvm(
      spec = spec,
      aux_sample = aux_sample,
      Psi_aux = Psi_aux,
      data = data_normalized,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = Psi_aux,
    vhat_method = if (is.function(vhat_fn)) "custom" else vhat_control$vhat_method,
    vhat_diagnostics = vhat_diagnostics,
    observed_cvm_distance_matrix = observed_cvm_distance_matrix,
    derivative_method = derivative_control$derivative_method,
    derivative_mc_size = derivative_control$derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_control$derivative_mc_seed)) NA_integer_ else derivative_control$derivative_mc_seed,
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

spec_cvm_bootstrap_stat <- function(spec,
                                    data,
                                    normalized_weights,
                                    theta_star,
                                    cvm_prep,
                                    null,
                                    control = list(),
                                    scale_factor = 1) {
  if (!is.function(spec$cvm_bootstrap_stat)) {
    return(NULL)
  }

  spec$cvm_bootstrap_stat(
    data = data,
    normalized_weights = normalized_weights,
    theta_star = theta_star,
    cvm_prep = cvm_prep,
    null = null,
    control = control,
    scale_factor = scale_factor
  )
}

normalize_probability_weights <- function(weights, n_expected = NULL) {
  weights <- as.numeric(weights)

  if (!is.null(n_expected) && length(weights) != n_expected) {
    stop("`weights` has incompatible length.")
  }
  if (length(weights) == 0) {
    stop("`weights` cannot be empty.")
  }
  if (any(!is.finite(weights))) {
    stop("`weights` must be finite.")
  }
  if (any(weights < 0)) {
    stop("`weights` must be nonnegative.")
  }

  total_weight <- sum(weights)
  if (total_weight <= 0) {
    stop("`weights` must have strictly positive sum.")
  }

  weights / total_weight
}

normalize_normal_data <- function(data, control = list()) {
  if (is.matrix(data) || is.data.frame(data)) {
    if (ncol(as.matrix(data)) != 1L) {
      stop("Normal data must be a numeric vector or a one-column matrix/data frame.")
    }
    data <- as.matrix(data)[, 1]
  }

  data <- as.numeric(data)
  if (length(data) == 0) {
    stop("`data` cannot be empty.")
  }
  if (any(!is.finite(data))) {
    stop("Normal data must be finite.")
  }

  data
}

normalize_normal_theta <- function(theta) {
  if (is.numeric(theta) && length(theta) == 2L && is.null(names(theta))) {
    theta <- list(mu = theta[[1]], sigma = theta[[2]])
  }

  if (!is.list(theta)) {
    stop("Normal theta must be a list with entries `mu` and `sigma`.")
  }

  mu <- as.numeric(theta$mu)
  sigma <- as.numeric(theta$sigma)

  if (length(mu) != 1L || !is.finite(mu)) {
    stop("Normal theta requires a finite scalar `mu`.")
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("Normal theta requires a strictly positive finite scalar `sigma`.")
  }

  list(mu = mu, sigma = sigma)
}

fit_normal_theta <- function(data,
                             weights = NULL,
                             null,
                             unknown_param = NULL,
                             control = list()) {
  x <- normalize_normal_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_normal_theta(null$theta))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (is.null(unknown_param)) {
    stop("Normal composite null requires `unknown_param` in the spec.")
  }

  fixed <- null$fixed %||% list()
  prob_weights <- if (is.null(weights)) {
    rep.int(1 / length(x), length(x))
  } else {
    normalize_probability_weights(weights, length(x))
  }

  weighted_mean <- function(values) {
    sum(prob_weights * values)
  }

  weighted_sigma <- function(values, center) {
    sqrt(sum(prob_weights * (values - center)^2))
  }

  if (identical(unknown_param, "mu")) {
    sigma_fixed <- as.numeric(fixed$sigma)
    if (length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0) {
      stop("Composite normal null with unknown `mu` requires `null$fixed$sigma`.")
    }
    return(list(mu = weighted_mean(x), sigma = sigma_fixed))
  }

  if (identical(unknown_param, "sigma")) {
    mu_fixed <- as.numeric(fixed$mu)
    if (length(mu_fixed) != 1L || !is.finite(mu_fixed)) {
      stop("Composite normal null with unknown `sigma` requires `null$fixed$mu`.")
    }
    return(list(mu = mu_fixed, sigma = weighted_sigma(x, mu_fixed)))
  }

  if (identical(unknown_param, "both")) {
    mu_hat <- weighted_mean(x)
    return(list(mu = mu_hat, sigma = weighted_sigma(x, mu_hat)))
  }

  stop("Unsupported `unknown_param` for the normal model.")
}

prepare_normal_fast_multiplier <- function(spec,
                                           data,
                                           theta_hat,
                                           ks_prep = NULL,
                                           cvm_prep = NULL,
                                           control = list(),
                                           unknown_param = "both") {
  x <- normalize_normal_data(data, control)
  theta_hat <- normalize_normal_theta(theta_hat)
  unknown_param <- tolower(as.character(unknown_param %||% "both"))

  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the fast normal multiplier bootstrap.")
  }

  par0 <- switch(
    unknown_param,
    mu = c(theta_hat$mu),
    sigma = c(log(theta_hat$sigma)),
    both = c(theta_hat$mu, log(theta_hat$sigma))
  )

  state_from_par <- function(par) {
    if (identical(unknown_param, "mu")) {
      normalize_normal_theta(list(mu = par[[1L]], sigma = theta_hat$sigma))
    } else if (identical(unknown_param, "sigma")) {
      normalize_normal_theta(list(mu = theta_hat$mu, sigma = exp(par[[1L]])))
    } else {
      normalize_normal_theta(list(mu = par[[1L]], sigma = exp(par[[2L]])))
    }
  }

  score_matrix_fn <- function(sample, par) {
    theta <- state_from_par(par)
    sample <- normalize_normal_data(sample)
    centered <- sample - theta$mu
    sigma_sq <- theta$sigma^2

    if (identical(unknown_param, "mu")) {
      return(matrix(centered / sigma_sq, ncol = 1L))
    }
    if (identical(unknown_param, "sigma")) {
      return(matrix(-1 + (centered^2) / sigma_sq, ncol = 1L))
    }

    cbind(
      centered / sigma_sq,
      -1 + (centered^2) / sigma_sq
    )
  }

  sample_fn <- function(n_aux, par) {
    theta <- state_from_par(par)
    stats::rnorm(n_aux, mean = theta$mu, sd = theta$sigma)
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

make_normal_spec <- function(unknown_param = NULL) {
  new_model_spec(
    name = "normal",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_normal_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_normal_data(data, control)
      omega_vec <- as.numeric(omega)
      if (length(omega_vec) == 0) {
        stop("`omega` cannot be empty.")
      }
      abs(outer(x, omega_vec, FUN = "-"))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_normal_theta(theta)
      theoretical_distance_profile_normal(
        omega = as.numeric(omega),
        mu = theta$mu,
        sigma = theta$sigma,
        t_values = as.numeric(t)
      )
    },
    normalize_data = normalize_normal_data,
    n_obs = function(data, control = list()) {
      length(normalize_normal_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_normal_data(data, control)[[idx]]
    },
    extras = list(
      unknown_param = unknown_param,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        if (is.null(unknown_param)) {
          stop("The fast normal multiplier bootstrap requires `unknown_param` in the model spec.")
        }
        prepare_normal_fast_multiplier(
          spec = make_normal_spec(unknown_param = unknown_param),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          unknown_param = unknown_param
        )
      }
    )
  )
}

normalize_mvnormal_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (!is.matrix(data) || nrow(data) == 0L || ncol(data) == 0L) {
    stop("Multivariate normal data must be a non-empty matrix.")
  }

  data <- matrix(
    as.numeric(data),
    nrow = nrow(data),
    ncol = ncol(data)
  )

  if (any(!is.finite(data))) {
    stop("Multivariate normal data must be finite.")
  }

  data
}

normalize_mvnormal_theta <- function(theta,
                                     ambient_dim = NULL,
                                     control = list()) {
  if (!is.list(theta)) {
    stop("Multivariate normal theta must be a list containing `mu` and `Sigma`.")
  }

  mu <- as.numeric(theta$mu)
  Sigma <- theta$Sigma

  if (length(mu) == 0L || any(!is.finite(mu))) {
    stop("Multivariate normal theta requires a finite vector `mu`.")
  }

  ambient_dim <- ambient_dim %||% theta$ambient_dim %||% length(mu)
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 1L) {
    stop("`ambient_dim` must be a strictly positive integer.")
  }
  if (length(mu) != ambient_dim) {
    stop("Multivariate normal theta has incompatible dimension.")
  }

  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(ambient_dim, ambient_dim))) {
    stop("Multivariate normal theta requires a square covariance matrix `Sigma`.")
  }
  if (any(!is.finite(Sigma))) {
    stop("Multivariate normal covariance must be finite.")
  }

  Sigma <- 0.5 * (Sigma + t(Sigma))
  tol <- as.numeric(control$mvnormal_cov_tol %||% control$logistic_gaussian_cov_tol %||% 1e-10)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`control$mvnormal_cov_tol` must be a strictly positive finite scalar.")
  }

  eigen_sigma <- eigen(Sigma, symmetric = TRUE)
  if (any(eigen_sigma$values < -tol)) {
    stop("Multivariate normal covariance must be positive semidefinite.")
  }

  eigenvalues_full <- pmax(eigen_sigma$values, 0)
  positive_idx <- eigenvalues_full > tol

  list(
    mu = mu,
    Sigma = Sigma,
    ambient_dim = ambient_dim,
    eigenvalues_full = eigenvalues_full,
    eigenvectors_full = eigen_sigma$vectors,
    positive_idx = positive_idx,
    rank = sum(positive_idx)
  )
}

mvnormal_weighted_covariance <- function(x, prob_weights, center) {
  centered <- sweep(x, 2L, center, FUN = "-")
  crossprod(centered * sqrt(prob_weights))
}

rmvnormal_euclidean <- function(n, mu, Sigma, control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- normalize_mvnormal_theta(
    list(mu = mu, Sigma = Sigma),
    control = control
  )

  mvtnorm::rmvnorm(n = n, mean = theta$mu, sigma = theta$Sigma)
}

fit_mvnormal_theta <- function(data,
                               weights = NULL,
                               null,
                               unknown_param = "both",
                               control = list()) {
  x <- normalize_mvnormal_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_mvnormal_theta(
      null$theta,
      ambient_dim = ncol(x),
      control = control
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the multivariate normal model.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }

  fixed <- null$fixed %||% list()
  fixed_mu_raw <- fixed$mu
  fixed_sigma_raw <- fixed$Sigma
  weighted_mean <- colSums(x * prob_weights)

  if (identical(unknown_param, "mu")) {
    if (is.null(fixed_sigma_raw)) {
      stop("Composite multivariate normal null with unknown `mu` requires `null$fixed$Sigma`.")
    }
    theta_out <- list(
      mu = weighted_mean,
      Sigma = normalize_mvnormal_theta(
        list(mu = rep.int(0, ncol(x)), Sigma = fixed_sigma_raw),
        ambient_dim = ncol(x),
        control = control
      )$Sigma
    )
    return(normalize_mvnormal_theta(theta_out, ambient_dim = ncol(x), control = control))
  }

  if (identical(unknown_param, "sigma")) {
    if (is.null(fixed_mu_raw)) {
      stop("Composite multivariate normal null with unknown `Sigma` requires `null$fixed$mu`.")
    }
    fixed_mu <- normalize_mvnormal_theta(
      list(mu = fixed_mu_raw, Sigma = diag(ncol(x))),
      ambient_dim = ncol(x),
      control = control
    )$mu
    theta_out <- list(
      mu = fixed_mu,
      Sigma = mvnormal_weighted_covariance(x, prob_weights, fixed_mu)
    )
    return(normalize_mvnormal_theta(theta_out, ambient_dim = ncol(x), control = control))
  }

  theta_out <- list(
    mu = weighted_mean,
    Sigma = mvnormal_weighted_covariance(x, prob_weights, weighted_mean)
  )
  normalize_mvnormal_theta(theta_out, ambient_dim = ncol(x), control = control)
}

prepare_mvnormal_fast_multiplier <- function(spec,
                                             data,
                                             theta_hat,
                                             ks_prep = NULL,
                                             cvm_prep = NULL,
                                             control = list(),
                                             unknown_param = "both") {
  normalized <- normalize_mvnormal_data(data, control)
  theta_hat <- normalize_mvnormal_theta(theta_hat, ambient_dim = ncol(normalized), control = control)
  d <- theta_hat$ambient_dim
  sigma_vech_hat <- fast_multiplier_vech(theta_hat$Sigma)
  unknown_param <- tolower(as.character(unknown_param %||% "both"))

  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the fast multivariate normal multiplier bootstrap.")
  }

  par0 <- switch(
    unknown_param,
    mu = theta_hat$mu,
    sigma = sigma_vech_hat,
    both = c(theta_hat$mu, sigma_vech_hat)
  )

  gaussian_score_matrix_from_matrix <- function(x, theta) {
    Sigma_inv <- solve(theta$Sigma)
    centered <- sweep(x, 2L, theta$mu, FUN = "-")
    score_mu <- centered %*% t(Sigma_inv)
    score_sigma <- t(vapply(seq_len(nrow(x)), function(i) {
      rr <- centered[i, , drop = FALSE]
      matrix_score <- 0.5 * (Sigma_inv %*% crossprod(rr) %*% Sigma_inv - Sigma_inv)
      fast_multiplier_sym_score_to_vech(matrix_score)
    }, numeric(d * (d + 1L) / 2L)))

    if (identical(unknown_param, "mu")) {
      return(score_mu)
    }
    if (identical(unknown_param, "sigma")) {
      return(score_sigma)
    }
    cbind(score_mu, score_sigma)
  }

  gaussian_influence_matrix_from_matrix <- function(x, theta) {
    centered <- sweep(x, 2L, theta$mu, FUN = "-")
    if_mu <- centered
    if_sigma <- t(vapply(seq_len(nrow(x)), function(i) {
      rr <- centered[i, , drop = FALSE]
      fast_multiplier_vech(crossprod(rr) - theta$Sigma)
    }, numeric(d * (d + 1L) / 2L)))

    if (identical(unknown_param, "mu")) {
      return(if_mu)
    }
    if (identical(unknown_param, "sigma")) {
      return(if_sigma)
    }
    cbind(if_mu, if_sigma)
  }

  derivative_control <- fast_multiplier_parse_derivative_control(control)
  S_obs <- gaussian_influence_matrix_from_matrix(normalized, theta_hat)

  if (!is.null(derivative_control$derivative_mc_seed)) {
    set.seed(derivative_control$derivative_mc_seed)
  }
  aux_sample <- rmvnormal_euclidean(
    n = derivative_control$derivative_mc_size,
    mu = theta_hat$mu,
    Sigma = theta_hat$Sigma,
    control = control
  )
  Psi_aux <- gaussian_score_matrix_from_matrix(aux_sample, theta_hat)
  Vhat <- diag(ncol(S_obs))
  vhat_diagnostics <- fast_multiplier_vhat_diagnostics(
    S_obs = S_obs,
    Psi_aux = Psi_aux,
    Vhat = Vhat,
    par0 = par0
  )

  observed_cvm_distance_matrix <- NULL
  if (!is.null(cvm_prep) && !isTRUE(cvm_prep$light)) {
    observed_cvm_distance_matrix <- cvm_prep$distance_matrix
    if (is.null(observed_cvm_distance_matrix)) {
      observed_cvm_distance_matrix <- spec$distance_matrix(normalized, normalized, control)
    }
    observed_cvm_distance_matrix <- as.matrix(observed_cvm_distance_matrix)
  }

  D_ks <- fast_multiplier_compute_D_ks(
    spec = spec,
    aux_sample = aux_sample,
    Psi_aux = Psi_aux,
    ks_prep = ks_prep,
    control = control
  )
  D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
  if (is.null(D_cvm)) {
    D_cvm <- fast_multiplier_compute_D_cvm(
      spec = spec,
      aux_sample = aux_sample,
      Psi_aux = Psi_aux,
      data = normalized,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = Psi_aux,
    vhat_method = "gaussian_mle_influence_identity",
    vhat_diagnostics = vhat_diagnostics,
    observed_cvm_distance_matrix = observed_cvm_distance_matrix,
    derivative_method = derivative_control$derivative_method,
    derivative_mc_size = derivative_control$derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_control$derivative_mc_seed)) NA_integer_ else derivative_control$derivative_mc_seed,
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

make_mvnormal_spec <- function(unknown_param = "both") {
  new_model_spec(
    name = "mvnormal_euclidean",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_mvnormal_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_mvnormal_data(data, control)
      omega_matrix <- normalize_mvnormal_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      x_sq <- rowSums(x^2)
      omega_sq <- rowSums(omega_matrix^2)
      sq_distances <- outer(x_sq, omega_sq, FUN = "+") - 2 * (x %*% t(omega_matrix))
      sqrt(pmax(sq_distances, 0))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_mvnormal_theta(theta, control = control)
      omega_vec <- as.numeric(normalize_mvnormal_data(omega, control)[1L, , drop = TRUE])
      evaluate_mvnorm_distance_profile(
        shift = theta$mu - omega_vec,
        t_values = as.numeric(t),
        eigenvalues_full = theta$eigenvalues_full,
        eigenvectors_full = theta$eigenvectors_full,
        positive_idx = theta$positive_idx,
        control = control
      )
    },
    normalize_data = normalize_mvnormal_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_mvnormal_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_mvnormal_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_mvnormal_theta(theta, control = control)
        omega_matrix <- normalize_mvnormal_data(omega_grid, control)
        if (ncol(omega_matrix) != theta$ambient_dim) {
          stop("`omega_grid` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu, each = nrow(omega_matrix)),
          nrow = nrow(omega_matrix),
          ncol = length(theta$mu)
        ) - omega_matrix
        t_matrix <- matrix(
          rep(as.numeric(t_grid), each = nrow(shift_matrix)),
          nrow = nrow(shift_matrix),
          ncol = length(t_grid)
        )
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = t_matrix,
          eigenvalues_full = theta$eigenvalues_full,
          eigenvectors_full = theta$eigenvectors_full,
          positive_idx = theta$positive_idx,
          control = control
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_mvnormal_theta(theta, control = control)
        x <- normalize_mvnormal_data(data, control)
        if (ncol(x) != theta$ambient_dim) {
          stop("`data` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu, each = nrow(x)),
          nrow = nrow(x),
          ncol = length(theta$mu)
        ) - x
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = distance_matrix,
          eigenvalues_full = theta$eigenvalues_full,
          eigenvectors_full = theta$eigenvectors_full,
          positive_idx = theta$positive_idx,
          control = control
        )
      },
      unknown_param = unknown_param,
      distance_type = "euclidean",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_mvnormal_fast_multiplier(
          spec = make_mvnormal_spec(unknown_param = unknown_param),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          unknown_param = unknown_param
        )
      }
    )
  )
}

logistic_gaussian_ilr_basis <- function(ambient_dim) {
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 2L) {
    stop("`ambient_dim` must be an integer greater than or equal to 2.")
  }

  basis <- matrix(0, nrow = ambient_dim, ncol = ambient_dim - 1L)
  for (j in seq_len(ambient_dim - 1L)) {
    basis[seq_len(j), j] <- 1 / sqrt(j * (j + 1))
    basis[j + 1L, j] <- -sqrt(j / (j + 1))
  }

  basis
}

row_softmax_matrix <- function(values) {
  values <- as.matrix(values)
  row_max <- apply(values, 1L, max)
  shifted <- sweep(values, 1L, row_max, FUN = "-")
  exp_shifted <- exp(shifted)
  exp_shifted / rowSums(exp_shifted)
}

coerce_logistic_gaussian_simplex_matrix <- function(data) {
  if (inherits(data, "logistic_gaussian_simplex_data")) {
    return(data$simplex)
  }

  if (is.list(data) && is.null(dim(data))) {
    if (length(data) == 0L) {
      stop("`data` cannot be empty.")
    }
    data <- do.call(rbind, lapply(data, as.numeric))
  } else if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (!is.matrix(data) || nrow(data) == 0L || ncol(data) < 2L) {
    stop("Logistic Gaussian data must be a non-empty matrix with at least two columns.")
  }

  data
}

normalize_logistic_gaussian_data <- function(data, control = list()) {
  if (inherits(data, "logistic_gaussian_simplex_data")) {
    return(data)
  }

  simplex <- coerce_logistic_gaussian_simplex_matrix(data)
  simplex <- matrix(
    as.numeric(simplex),
    nrow = nrow(simplex),
    ncol = ncol(simplex)
  )

  if (any(!is.finite(simplex))) {
    stop("Logistic Gaussian data must be finite.")
  }
  if (any(simplex <= 0)) {
    stop("Logistic Gaussian data must lie in the interior of the simplex.")
  }

  row_sums <- rowSums(simplex)
  tol <- as.numeric(control$logistic_gaussian_simplex_tol %||% 1e-8)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`control$logistic_gaussian_simplex_tol` must be a strictly positive finite scalar.")
  }
  if (any(abs(row_sums - 1) > tol)) {
    stop("Logistic Gaussian data rows must sum to 1 up to tolerance.")
  }

  simplex <- simplex / row_sums
  log_simplex <- log(simplex)
  clr_matrix <- log_simplex - rowMeans(log_simplex)
  basis <- logistic_gaussian_ilr_basis(ncol(simplex))
  ilr_matrix <- clr_matrix %*% basis

  structure(
    list(
      simplex = simplex,
      ilr = ilr_matrix,
      basis = basis,
      ambient_dim = ncol(simplex),
      ilr_dim = ncol(simplex) - 1L
    ),
    class = c("logistic_gaussian_simplex_data", "list")
  )
}

logistic_gaussian_simplex_matrix <- function(data, control = list()) {
  normalize_logistic_gaussian_data(data, control)$simplex
}

logistic_gaussian_ilr_matrix <- function(data, control = list()) {
  normalize_logistic_gaussian_data(data, control)$ilr
}

logistic_gaussian_point_to_ilr <- function(point, ambient_dim = NULL, control = list()) {
  normalized <- normalize_logistic_gaussian_data(point, control)
  if (!is.null(ambient_dim) && normalized$ambient_dim != as.integer(ambient_dim)) {
    stop("Logistic Gaussian point has incompatible ambient dimension.")
  }
  as.numeric(normalized$ilr[1L, , drop = TRUE])
}

logistic_gaussian_ilr_to_simplex <- function(z, ambient_dim = NULL) {
  if (is.vector(z)) {
    z <- matrix(as.numeric(z), nrow = 1L)
  } else {
    z <- as.matrix(z)
  }

  if (nrow(z) == 0L || ncol(z) == 0L || any(!is.finite(z))) {
    stop("`z` must be a non-empty finite matrix.")
  }

  ilr_dim <- ncol(z)
  ambient_dim <- ambient_dim %||% (ilr_dim + 1L)
  ambient_dim <- as.integer(ambient_dim)
  if (ambient_dim != ilr_dim + 1L) {
    stop("`ambient_dim` is incompatible with the number of ilr coordinates.")
  }

  basis <- logistic_gaussian_ilr_basis(ambient_dim)
  clr_matrix <- z %*% t(basis)
  row_softmax_matrix(clr_matrix)
}

rlogistic_gaussian_simplex <- function(n, mu_ilr, Sigma_ilr, control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- normalize_logistic_gaussian_theta(
    list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    control = control
  )
  z <- mvtnorm::rmvnorm(n = n, mean = theta$mu_ilr, sigma = theta$Sigma_ilr)
  logistic_gaussian_ilr_to_simplex(z, ambient_dim = theta$ambient_dim)
}

normalize_logistic_gaussian_theta <- function(theta,
                                              ambient_dim = NULL,
                                              control = list()) {
  if (!is.list(theta)) {
    stop("Logistic Gaussian theta must be a list.")
  }

  mu_ilr <- as.numeric(theta$mu_ilr %||% theta$mu)
  Sigma_ilr <- theta$Sigma_ilr %||% theta$Sigma

  if (length(mu_ilr) == 0L || any(!is.finite(mu_ilr))) {
    stop("Logistic Gaussian theta requires a finite vector `mu_ilr`.")
  }

  ambient_dim <- ambient_dim %||% theta$ambient_dim %||% (length(mu_ilr) + 1L)
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 2L) {
    stop("`ambient_dim` must be an integer greater than or equal to 2.")
  }
  if (length(mu_ilr) != ambient_dim - 1L) {
    stop("Logistic Gaussian theta has incompatible ilr dimension.")
  }

  Sigma_ilr <- as.matrix(Sigma_ilr)
  if (!all(dim(Sigma_ilr) == c(length(mu_ilr), length(mu_ilr)))) {
    stop("Logistic Gaussian theta requires a square covariance matrix `Sigma_ilr`.")
  }
  if (any(!is.finite(Sigma_ilr))) {
    stop("Logistic Gaussian covariance must be finite.")
  }

  Sigma_ilr <- 0.5 * (Sigma_ilr + t(Sigma_ilr))
  tol <- as.numeric(control$logistic_gaussian_cov_tol %||% 1e-10)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`control$logistic_gaussian_cov_tol` must be a strictly positive finite scalar.")
  }

  eigen_sigma <- eigen(Sigma_ilr, symmetric = TRUE)
  if (any(eigen_sigma$values < -tol)) {
    stop("Logistic Gaussian covariance must be positive semidefinite.")
  }

  eigenvalues_full <- pmax(eigen_sigma$values, 0)
  positive_idx <- eigenvalues_full > tol
  basis <- theta$basis %||% logistic_gaussian_ilr_basis(ambient_dim)

  list(
    mu_ilr = mu_ilr,
    Sigma_ilr = Sigma_ilr,
    basis = basis,
    ambient_dim = ambient_dim,
    ilr_dim = ambient_dim - 1L,
    eigenvalues_full = eigenvalues_full,
    eigenvectors_full = eigen_sigma$vectors,
    positive_idx = positive_idx,
    rank = sum(positive_idx)
  )
}

evaluate_mvnorm_distance_profile <- function(shift,
                                             t_values,
                                             eigenvalues_full,
                                             eigenvectors_full,
                                             positive_idx,
                                             tol = 1e-12,
                                             control = list()) {
  shift <- as.numeric(shift)
  t_values <- as.numeric(t_values)
  if (any(!is.finite(shift))) {
    stop("`shift` must be finite.")
  }
  if (any(!is.finite(t_values))) {
    stop("`t_values` must be finite.")
  }

  nu <- as.vector(crossprod(eigenvectors_full, shift))
  null_const <- sum(nu[!positive_idx]^2)
  threshold_sq <- t_values^2 - null_const
  output <- numeric(length(t_values))
  positive_threshold <- threshold_sq >= -tol

  if (!any(positive_idx)) {
    output[positive_threshold] <- 1
    return(pmax(pmin(output, 1), 0))
  }

  threshold_sq[threshold_sq < 0] <- 0
  if (any(positive_threshold)) {
    positive_values <- threshold_sq[positive_threshold]
    lambda_pos <- eigenvalues_full[positive_idx]
    nu_pos <- nu[positive_idx]
    delta <- nu_pos^2 / lambda_pos
    method <- tolower(as.character(control$logistic_gaussian_quadform_method %||% "auto"))
    if (identical(method, "hbe")) {
      output[positive_threshold] <- sphunif::p_wschisq(
        x = positive_values,
        weights = lambda_pos,
        dfs = rep.int(1, length(lambda_pos)),
        ncps = delta,
        method = "HBE"
      )
    } else {
      output[positive_threshold] <- vapply(
        positive_values,
        function(threshold) {
          if (!is.finite(threshold)) {
            return(NA_real_)
          }
          if (threshold <= 0) {
            return(0)
          }

          logistic_gaussian_quadform_tail_probability(
            q = threshold,
            lambda = lambda_pos,
            h = rep.int(1, length(lambda_pos)),
            delta = delta,
            control = control
          )
        },
        numeric(1)
      )
    }
  }

  pmax(pmin(output, 1), 0)
}

logistic_gaussian_quadform_tail_probability <- function(q,
                                                        lambda,
                                                        h,
                                                        delta,
                                                        control = list()) {
  # Empirical dispatcher for the weighted noncentral chi-squared CDFs that
  # arise in the logistic Gaussian distance-profile calculations.
  #
  # Benchmark evidence stored under `tests/benchmark_outputs/logistic_gaussian_quadform`
  # showed that:
  # - `farebrother` is typically the fastest accurate exact route in the
  #   regular regime;
  # - highly ill-conditioned spectra and strongly noncentral cases are the main
  #   regimes where `farebrother` becomes
  #   unreliable or expensive;
  # - in those ill-conditioned regimes `imhof` preserved accuracy well;
  # - if `imhof` also fails, `HBE` is the most accurate approximation among the
  #   fast alternatives we benchmarked.
  #
  # The ill-conditioning rule is empirical, not theoretical:
  #   cond(lambda) > 1e4, or
  #   length(lambda) >= 5, max(delta) > 10, and q / E(Q) > 0.1, or
  #   max(delta) > 1e3 and q / E(Q) > 1
  # with E(Q) = sum_j lambda_j (h_j + delta_j).
  method <- tolower(as.character(control$logistic_gaussian_quadform_method %||% "auto"))
  if (!method %in% c("auto", "farebrother", "imhof", "hbe", "davies")) {
    stop("`control$logistic_gaussian_quadform_method` must be one of 'auto', 'farebrother', 'imhof', 'hbe', or 'davies'.")
  }

  cond_threshold <- as.numeric(control$logistic_gaussian_ill_conditioned_cond_threshold %||% 1e4)
  dim_threshold <- as.integer(control$logistic_gaussian_ill_conditioned_dim_threshold %||% 5L)
  delta_threshold <- as.numeric(control$logistic_gaussian_ill_conditioned_delta_threshold %||% 10)
  q_ratio_threshold <- as.numeric(control$logistic_gaussian_ill_conditioned_q_ratio_threshold %||% 0.1)
  extreme_delta_threshold <- as.numeric(control$logistic_gaussian_ill_conditioned_extreme_delta_threshold %||% 1e3)
  extreme_q_ratio_threshold <- as.numeric(control$logistic_gaussian_ill_conditioned_extreme_q_ratio_threshold %||% 1)
  q_mean <- sum(lambda * (h + delta))
  q_ratio <- if (is.finite(q_mean) && q_mean > 0) q / q_mean else Inf
  cond_number <- max(lambda) / min(lambda)

  is_ill_conditioned <- is.finite(cond_number) &&
    (
      cond_number > cond_threshold ||
        (length(lambda) >= dim_threshold && max(delta) > delta_threshold && q_ratio > q_ratio_threshold) ||
        (max(delta) > extreme_delta_threshold && q_ratio > extreme_q_ratio_threshold)
    )

  auto_selected_method <- if (identical(method, "auto")) {
    if (is_ill_conditioned) "imhof" else "farebrother"
  } else {
    NA_character_
  }
  if (identical(method, "auto")) {
    method <- auto_selected_method
  }

  compute_farebrother <- function() {
    eps <- as.numeric(control$logistic_gaussian_quadform_eps %||% 1e-8)
    maxit <- as.integer(control$logistic_gaussian_quadform_maxit %||% 100000L)
    if (!is.finite(eps) || eps <= 0) {
      eps <- 1e-8
    }
    if (!is.finite(maxit) || maxit < 1000L) {
      maxit <- 100000L
    }
    res <- CompQuadForm::farebrother(
      q = q,
      lambda = lambda,
      h = h,
      delta = delta,
      maxit = maxit,
      eps = eps
    )
    list(prob = 1 - res$Qq, ifault = as.integer(res$ifault %||% NA_integer_))
  }

  compute_imhof <- function() {
    epsabs <- as.numeric(control$logistic_gaussian_quadform_imhof_epsabs %||% 1e-8)
    epsrel <- as.numeric(control$logistic_gaussian_quadform_imhof_epsrel %||% 1e-8)
    limit <- as.integer(control$logistic_gaussian_quadform_imhof_limit %||% 20000L)
    if (!is.finite(epsabs) || epsabs <= 0) {
      epsabs <- 1e-8
    }
    if (!is.finite(epsrel) || epsrel <= 0) {
      epsrel <- 1e-8
    }
    if (!is.finite(limit) || limit < 1000L) {
      limit <- 20000L
    }
    res <- CompQuadForm::imhof(
      q = q,
      lambda = lambda,
      h = h,
      delta = delta,
      epsabs = epsabs,
      epsrel = epsrel,
      limit = limit
    )
    list(prob = 1 - res$Qq)
  }

  compute_hbe <- function() {
    prob <- sphunif::p_wschisq(
      x = q,
      weights = lambda,
      dfs = h,
      ncps = delta,
      method = "HBE"
    )
    list(prob = prob)
  }

  is_valid_prob <- function(prob) {
    is.numeric(prob) && length(prob) == 1L && is.finite(prob) && prob >= 0 && prob <= 1
  }

  if (identical(method, "davies")) {
    acc <- as.numeric(control$logistic_gaussian_quadform_acc %||% 1e-8)
    lim <- as.integer(control$logistic_gaussian_quadform_lim %||% 20000L)
    if (!is.finite(acc) || acc <= 0) {
      acc <- 1e-8
    }
    if (!is.finite(lim) || lim < 1000L) {
      lim <- 20000L
    }
    res <- CompQuadForm::davies(
      q = q,
      lambda = lambda,
      h = h,
      delta = delta,
      acc = acc,
      lim = lim
    )
    return(pmin(pmax(1 - res$Qq, 0), 1))
  }

  if (identical(method, "farebrother")) {
    fare <- compute_farebrother()
    fare_prob <- pmin(pmax(fare$prob, 0), 1)
    if (identical(auto_selected_method, "farebrother")) {
      fare_valid <- is_valid_prob(fare$prob) && (is.na(fare$ifault) || fare$ifault == 0L)
      if (!fare_valid) {
        imh <- compute_imhof()
        if (is_valid_prob(imh$prob)) {
          return(imh$prob)
        }
        hbe <- compute_hbe()
        return(pmin(pmax(hbe$prob, 0), 1))
      }
    }
    return(fare_prob)
  }

  if (identical(method, "imhof")) {
    imh <- compute_imhof()
    if (is_valid_prob(imh$prob)) {
      return(imh$prob)
    }
    hbe <- compute_hbe()
    return(pmin(pmax(hbe$prob, 0), 1))
  }

  if (identical(method, "hbe")) {
    hbe <- compute_hbe()
    return(pmin(pmax(hbe$prob, 0), 1))
  }

  stop("Unsupported logistic Gaussian quadratic-form method.")
}

evaluate_mvnorm_distance_profile_matrix <- function(shift_matrix,
                                                    t_matrix,
                                                    eigenvalues_full,
                                                    eigenvectors_full,
                                                    positive_idx,
                                                    tol = 1e-12,
                                                    control = list()) {
  shift_matrix <- as.matrix(shift_matrix)
  t_matrix <- as.matrix(t_matrix)

  if (nrow(shift_matrix) == 0L || ncol(shift_matrix) == 0L) {
    stop("`shift_matrix` must be a non-empty matrix.")
  }
  if (nrow(t_matrix) != nrow(shift_matrix)) {
    stop("`t_matrix` and `shift_matrix` must have the same number of rows.")
  }
  if (any(!is.finite(shift_matrix))) {
    stop("`shift_matrix` must be finite.")
  }
  if (any(!is.finite(t_matrix))) {
    stop("`t_matrix` must be finite.")
  }

  nu_matrix <- shift_matrix %*% eigenvectors_full
  if (any(!positive_idx)) {
    null_const <- rowSums(nu_matrix[, !positive_idx, drop = FALSE]^2)
  } else {
    null_const <- rep.int(0, nrow(shift_matrix))
  }

  threshold_sq <- t_matrix^2 - matrix(
    null_const,
    nrow = nrow(t_matrix),
    ncol = ncol(t_matrix)
  )
  positive_threshold <- threshold_sq >= -tol
  output <- matrix(0, nrow = nrow(t_matrix), ncol = ncol(t_matrix))

  if (!any(positive_idx)) {
    output[positive_threshold] <- 1
    return(pmax(pmin(output, 1), 0))
  }

  threshold_sq[threshold_sq < 0] <- 0
  lambda_pos <- eigenvalues_full[positive_idx]
  nu_pos_matrix <- nu_matrix[, positive_idx, drop = FALSE]
  delta_matrix <- sweep(nu_pos_matrix^2, 2L, lambda_pos, FUN = "/")
  df_vec <- rep.int(1, length(lambda_pos))

  for (i in seq_len(nrow(t_matrix))) {
    positive_i <- positive_threshold[i, ]
    if (!any(positive_i)) {
      next
    }

    positive_values <- threshold_sq[i, positive_i]
    method <- tolower(as.character(control$logistic_gaussian_quadform_method %||% "auto"))
    if (identical(method, "hbe")) {
      output[i, positive_i] <- sphunif::p_wschisq(
        x = positive_values,
        weights = lambda_pos,
        dfs = df_vec,
        ncps = delta_matrix[i, ],
        method = "HBE"
      )
    } else {
      output[i, positive_i] <- vapply(
        positive_values,
        function(threshold) {
          if (!is.finite(threshold)) {
            return(NA_real_)
          }
          if (threshold <= 0) {
            return(0)
          }

          logistic_gaussian_quadform_tail_probability(
            q = threshold,
            lambda = lambda_pos,
            h = df_vec,
            delta = delta_matrix[i, ],
            control = control
          )
        },
        numeric(1)
      )
    }
  }
  pmax(pmin(output, 1), 0)
}

logistic_gaussian_weighted_covariance <- function(z, prob_weights, center) {
  centered <- sweep(z, 2L, center, FUN = "-")
  crossprod(centered * sqrt(prob_weights))
}

fit_logistic_gaussian_theta <- function(data,
                                        weights = NULL,
                                        null,
                                        unknown_param = "both",
                                        control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_logistic_gaussian_theta(
      null$theta,
      ambient_dim = normalized$ambient_dim,
      control = control
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the Logistic Gaussian model.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(z), nrow(z))
  } else {
    normalize_probability_weights(weights, nrow(z))
  }

  fixed <- null$fixed %||% list()
  fixed_mu_raw <- fixed$mu_ilr %||% fixed$mu
  fixed_sigma_raw <- fixed$Sigma_ilr %||% fixed$Sigma
  weighted_mean <- colSums(z * prob_weights)

  if (identical(unknown_param, "mu")) {
    if (is.null(fixed_sigma_raw)) {
      stop("Composite Logistic Gaussian null with unknown `mu` requires `null$fixed$Sigma_ilr`.")
    }
    theta_out <- list(
      mu_ilr = weighted_mean,
      Sigma_ilr = normalize_logistic_gaussian_theta(
        list(mu_ilr = rep.int(0, normalized$ilr_dim), Sigma_ilr = fixed_sigma_raw),
        ambient_dim = normalized$ambient_dim,
        control = control
      )$Sigma_ilr
    )
    return(normalize_logistic_gaussian_theta(theta_out, ambient_dim = normalized$ambient_dim, control = control))
  }

  if (identical(unknown_param, "sigma")) {
    if (is.null(fixed_mu_raw)) {
      stop("Composite Logistic Gaussian null with unknown `Sigma` requires `null$fixed$mu_ilr`.")
    }
    fixed_mu <- normalize_logistic_gaussian_theta(
      list(mu_ilr = fixed_mu_raw, Sigma_ilr = diag(normalized$ilr_dim)),
      ambient_dim = normalized$ambient_dim,
      control = control
    )$mu_ilr
    theta_out <- list(
      mu_ilr = fixed_mu,
      Sigma_ilr = logistic_gaussian_weighted_covariance(z, prob_weights, fixed_mu)
    )
    return(normalize_logistic_gaussian_theta(theta_out, ambient_dim = normalized$ambient_dim, control = control))
  }

  theta_out <- list(
    mu_ilr = weighted_mean,
    Sigma_ilr = logistic_gaussian_weighted_covariance(z, prob_weights, weighted_mean)
  )
  normalize_logistic_gaussian_theta(theta_out, ambient_dim = normalized$ambient_dim, control = control)
}

prepare_logistic_gaussian_fast_multiplier <- function(spec,
                                                      data,
                                                      theta_hat,
                                                      ks_prep = NULL,
                                                      cvm_prep = NULL,
                                                      control = list(),
                                                      unknown_param = "both") {
  normalized <- normalize_logistic_gaussian_data(data, control)
  theta_hat <- normalize_logistic_gaussian_theta(theta_hat, ambient_dim = normalized$ambient_dim, control = control)
  d <- theta_hat$ilr_dim
  sigma_vech_hat <- fast_multiplier_vech(theta_hat$Sigma_ilr)
  unknown_param <- tolower(as.character(unknown_param %||% "both"))

  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the fast Logistic Gaussian multiplier bootstrap.")
  }

  par0 <- switch(
    unknown_param,
    mu = theta_hat$mu_ilr,
    sigma = sigma_vech_hat,
    both = c(theta_hat$mu_ilr, sigma_vech_hat)
  )

  gaussian_score_matrix_from_ilr <- function(z, theta) {
    Sigma_inv <- solve(theta$Sigma_ilr)
    centered <- sweep(z, 2L, theta$mu_ilr, FUN = "-")
    score_mu <- centered %*% t(Sigma_inv)
    score_sigma <- t(vapply(seq_len(nrow(z)), function(i) {
      rr <- centered[i, , drop = FALSE]
      matrix_score <- 0.5 * (Sigma_inv %*% crossprod(rr) %*% Sigma_inv - Sigma_inv)
      fast_multiplier_sym_score_to_vech(matrix_score)
    }, numeric(d * (d + 1L) / 2L)))

    if (identical(unknown_param, "mu")) {
      return(score_mu)
    }
    if (identical(unknown_param, "sigma")) {
      return(score_sigma)
    }
    cbind(score_mu, score_sigma)
  }

  gaussian_influence_matrix_from_ilr <- function(z, theta) {
    centered <- sweep(z, 2L, theta$mu_ilr, FUN = "-")
    if_mu <- centered
    if_sigma <- t(vapply(seq_len(nrow(z)), function(i) {
      rr <- centered[i, , drop = FALSE]
      fast_multiplier_vech(crossprod(rr) - theta$Sigma_ilr)
    }, numeric(d * (d + 1L) / 2L)))

    if (identical(unknown_param, "mu")) {
      return(if_mu)
    }
    if (identical(unknown_param, "sigma")) {
      return(if_sigma)
    }
    cbind(if_mu, if_sigma)
  }

  derivative_control <- fast_multiplier_parse_derivative_control(control)
  z_obs <- normalized$ilr
  S_obs <- gaussian_influence_matrix_from_ilr(z_obs, theta_hat)

  if (!is.null(derivative_control$derivative_mc_seed)) {
    set.seed(derivative_control$derivative_mc_seed)
  }
  aux_sample <- rlogistic_gaussian_simplex(
    n = derivative_control$derivative_mc_size,
    mu_ilr = theta_hat$mu_ilr,
    Sigma_ilr = theta_hat$Sigma_ilr,
    control = control
  )
  z_aux <- logistic_gaussian_ilr_matrix(aux_sample, control)
  Psi_aux <- gaussian_score_matrix_from_ilr(z_aux, theta_hat)
  Vhat <- diag(ncol(S_obs))
  vhat_diagnostics <- fast_multiplier_vhat_diagnostics(
    S_obs = S_obs,
    Psi_aux = Psi_aux,
    Vhat = Vhat,
    par0 = par0
  )

  observed_cvm_distance_matrix <- NULL
  if (!is.null(cvm_prep) && !isTRUE(cvm_prep$light)) {
    observed_cvm_distance_matrix <- cvm_prep$distance_matrix
    if (is.null(observed_cvm_distance_matrix)) {
      observed_cvm_distance_matrix <- spec$distance_matrix(normalized, normalized, control)
    }
    observed_cvm_distance_matrix <- as.matrix(observed_cvm_distance_matrix)
  }

  D_ks <- fast_multiplier_compute_D_ks(
    spec = spec,
    aux_sample = aux_sample,
    Psi_aux = Psi_aux,
    ks_prep = ks_prep,
    control = control
  )
  D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
  if (is.null(D_cvm)) {
    D_cvm <- fast_multiplier_compute_D_cvm(
      spec = spec,
      aux_sample = aux_sample,
      Psi_aux = Psi_aux,
      data = normalized,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = Psi_aux,
    vhat_method = "gaussian_mle_influence_identity",
    vhat_diagnostics = vhat_diagnostics,
    observed_cvm_distance_matrix = observed_cvm_distance_matrix,
    derivative_method = derivative_control$derivative_method,
    derivative_mc_size = derivative_control$derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_control$derivative_mc_seed)) NA_integer_ else derivative_control$derivative_mc_seed,
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

make_logistic_gaussian_spec <- function(unknown_param = "both") {
  new_model_spec(
    name = "logistic_gaussian_aitchison",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_logistic_gaussian_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_logistic_gaussian_data(data, control)
      omega_normalized <- normalize_logistic_gaussian_data(omega, control)

      if (x$ambient_dim != omega_normalized$ambient_dim) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      x_sq <- rowSums(x$ilr^2)
      omega_sq <- rowSums(omega_normalized$ilr^2)
      sq_distances <- outer(x_sq, omega_sq, FUN = "+") - 2 * (x$ilr %*% t(omega_normalized$ilr))
      sqrt(pmax(sq_distances, 0))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_logistic_gaussian_theta(theta, control = control)
      omega_ilr <- logistic_gaussian_point_to_ilr(
        omega,
        ambient_dim = theta$ambient_dim,
        control = control
      )
      evaluate_mvnorm_distance_profile(
        shift = theta$mu_ilr - omega_ilr,
        t_values = as.numeric(t),
        eigenvalues_full = theta$eigenvalues_full,
        eigenvectors_full = theta$eigenvectors_full,
        positive_idx = theta$positive_idx,
        control = control
      )
    },
    normalize_data = normalize_logistic_gaussian_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_logistic_gaussian_data(data, control)$simplex)
    },
    observation_at = function(data, idx, control = list()) {
      normalize_logistic_gaussian_data(data, control)$simplex[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_logistic_gaussian_theta(theta, control = control)
        omega_normalized <- normalize_logistic_gaussian_data(omega_grid, control)
        if (omega_normalized$ambient_dim != theta$ambient_dim) {
          stop("`omega_grid` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu_ilr, each = nrow(omega_normalized$ilr)),
          nrow = nrow(omega_normalized$ilr),
          ncol = length(theta$mu_ilr)
        ) - omega_normalized$ilr
        t_matrix <- matrix(
          rep(as.numeric(t_grid), each = nrow(shift_matrix)),
          nrow = nrow(shift_matrix),
          ncol = length(t_grid)
        )
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = t_matrix,
          eigenvalues_full = theta$eigenvalues_full,
          eigenvectors_full = theta$eigenvectors_full,
          positive_idx = theta$positive_idx,
          control = control
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_logistic_gaussian_theta(theta, control = control)
        data_normalized <- normalize_logistic_gaussian_data(data, control)
        if (data_normalized$ambient_dim != theta$ambient_dim) {
          stop("`data` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu_ilr, each = nrow(data_normalized$ilr)),
          nrow = nrow(data_normalized$ilr),
          ncol = length(theta$mu_ilr)
        ) - data_normalized$ilr
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = distance_matrix,
          eigenvalues_full = theta$eigenvalues_full,
          eigenvectors_full = theta$eigenvectors_full,
          positive_idx = theta$positive_idx,
          control = control
        )
      },
      unknown_param = unknown_param,
      distance_type = "aitchison",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_logistic_gaussian_fast_multiplier(
          spec = make_logistic_gaussian_spec(unknown_param = unknown_param),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          unknown_param = unknown_param
        )
      }
    )
  )
}

normalize_vmf_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0 || ncol(data) < 2L) {
    stop("vMF data must be a non-empty matrix with at least two columns.")
  }
  if (any(!is.finite(data))) {
    stop("vMF data must be finite.")
  }

  norms <- sqrt(rowSums(data^2))
  if (any(norms <= 0)) {
    stop("vMF data rows must have strictly positive norm.")
  }

  data / norms
}

default_unit_vector <- function(dim) {
  output <- rep(0, dim)
  output[[1]] <- 1
  output
}

normalize_vmf_theta <- function(theta, ambient_dim = NULL) {
  if (is.numeric(theta) && !is.list(theta)) {
    xi <- as.numeric(theta)
    if (length(xi) < 2L) {
      stop("vMF xi must have length at least 2.")
    }
    if (any(!is.finite(xi))) {
      stop("vMF xi must be finite.")
    }

    kappa <- sqrt(sum(xi^2))
    mu <- if (kappa > 0) {
      xi / kappa
    } else {
      default_unit_vector(length(xi))
    }

    return(list(
      xi = xi,
      mu = mu,
      kappa = kappa,
      q = length(xi) - 1L
    ))
  }

  if (!is.list(theta)) {
    stop("vMF theta must be either a xi vector or a list containing `mu` and `kappa`.")
  }

  if (!is.null(theta$xi)) {
    return(normalize_vmf_theta(theta$xi, ambient_dim = ambient_dim))
  }

  mu <- as.numeric(theta$mu)
  kappa <- as.numeric(theta$kappa)

  if (length(mu) < 2L || any(!is.finite(mu))) {
    stop("vMF theta requires a finite vector `mu` of length at least 2.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("vMF theta requires a nonnegative finite scalar `kappa`.")
  }

  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 0) {
    stop("vMF theta requires `mu` with strictly positive norm.")
  }

  mu <- mu / mu_norm
  xi <- if (kappa > 0) kappa * mu else rep(0, length(mu))
  if (!is.null(ambient_dim) && length(mu) != ambient_dim) {
    stop("vMF theta has incompatible dimension.")
  }

  list(
    xi = xi,
    mu = mu,
    kappa = kappa,
    q = length(mu) - 1L
  )
}

solve_vmf_kappa_from_rbar <- function(r_bar,
                                      q,
                                      tol = 1e-10,
                                      max_kappa = 1e6) {
  r_bar <- as.numeric(r_bar)
  if (length(r_bar) != 1L || !is.finite(r_bar) || r_bar < 0 || r_bar > 1) {
    stop("`r_bar` must be a scalar in [0, 1].")
  }
  if (length(q) != 1L || !is.finite(q) || q < 1) {
    stop("`q` must be a positive scalar.")
  }

  if (r_bar <= tol) {
    return(0)
  }
  if (r_bar >= 1 - tol) {
    return(max_kappa)
  }

  objective <- function(kappa) A_q(kappa, q) - r_bar
  upper <- 1

  while (objective(upper) < 0 && upper < max_kappa) {
    upper <- upper * 2
  }

  if (upper >= max_kappa && objective(upper) < 0) {
    return(max_kappa)
  }

  stats::uniroot(objective, interval = c(0, upper), tol = tol)$root
}

fit_vmf_theta <- function(data,
                          weights = NULL,
                          null,
                          unknown_param = "xi",
                          control = list()) {
  x <- normalize_vmf_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_vmf_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (!identical(unknown_param, "xi")) {
    stop("The vMF adapter currently supports only `unknown_param = 'xi'`.")
  }

  if (is.null(weights)) {
    xi_hat <- compute_mle_xi(x)
    return(normalize_vmf_theta(xi_hat, ambient_dim = ncol(x)))
  }

  prob_weights <- normalize_probability_weights(weights, nrow(x))
  resultant <- colSums(x * prob_weights)
  r_bar <- sqrt(sum(resultant^2))

  if (r_bar <= 1e-12) {
    mu_hat <- default_unit_vector(ncol(x))
    kappa_hat <- 0
  } else {
    mu_hat <- resultant / r_bar
    kappa_hat <- solve_vmf_kappa_from_rbar(r_bar, q = ncol(x) - 1L)
  }

  normalize_vmf_theta(
    list(mu = mu_hat, kappa = kappa_hat),
    ambient_dim = ncol(x)
  )
}

prepare_vmf_fast_multiplier <- function(data,
                                        theta_hat,
                                        spec,
                                        ks_prep = NULL,
                                        cvm_prep = NULL,
                                        control = list(),
                                        distance_type = "geodesic") {
  x <- normalize_vmf_data(data, control)
  theta_hat <- normalize_vmf_theta(theta_hat, ambient_dim = ncol(x))
  p <- length(theta_hat$xi)
  q <- theta_hat$q
  derivative_control <- fast_multiplier_parse_derivative_control(control)

  S_obs <- t(vapply(seq_len(nrow(x)), function(i) {
    psi_xi(x[i, ], theta_hat$xi, q)
  }, numeric(p)))
  Vhat <- -dot_psi_xi(theta_hat$xi, q)

  if (isTRUE(any(!is.finite(Vhat)))) {
    stop("The vMF fast multiplier preparation produced a non-finite `Vhat`.")
  }

  if (!is.null(derivative_control$derivative_mc_seed)) {
    set.seed(derivative_control$derivative_mc_seed)
  }
  aux_sample <- rotasym::r_vMF(
    derivative_control$derivative_mc_size,
    mu = theta_hat$mu,
    kappa = theta_hat$kappa
  )
  aux_sample <- normalize_vmf_data(aux_sample, control)
  Psi_aux <- t(vapply(seq_len(nrow(aux_sample)), function(i) {
    psi_xi(aux_sample[i, ], theta_hat$xi, q)
  }, numeric(p)))

  D_ks <- fast_multiplier_compute_D_ks(
    spec = spec,
    aux_sample = aux_sample,
    Psi_aux = Psi_aux,
    ks_prep = ks_prep,
    control = control
  )
  D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
  if (is.null(D_cvm)) {
    D_cvm <- fast_multiplier_compute_D_cvm(
      spec = spec,
      aux_sample = aux_sample,
      Psi_aux = Psi_aux,
      data = x,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = Psi_aux,
    derivative_method = derivative_control$derivative_method,
    derivative_mc_size = derivative_control$derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_control$derivative_mc_seed)) NA_integer_ else derivative_control$derivative_mc_seed,
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

make_vmf_spec <- function(distance_type = c("chordal", "geodesic"),
                          unknown_param = "xi") {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("vmf_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_vmf_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_vmf_data(data, control)
      omega_matrix <- normalize_vmf_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- x %*% t(omega_matrix)
      dot_products <- pmax(pmin(dot_products, 1), -1)

      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_vmf_theta(theta)
      ambient_dim <- length(theta$mu)
      profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
      l_max <- control$vmf_profile_l_max %||% NULL
      tail_tol <- as.numeric(control$vmf_profile_legendre_tail_tol %||% 1e-10)

      if (ambient_dim == 3L) {
        if (identical(profile_method, "legendre")) {
          return(distance_profile_vmf_s2_legendre(
            omega = omega,
            mu = theta$mu,
            kappa = theta$kappa,
            t_values = as.numeric(t),
            distance_type = distance_type,
            l_max = l_max,
            tail_tol = tail_tol
          ))
        }

        return(theoretical_distance_profile_vmf_s2_fast(
          omega = omega,
          mu = theta$mu,
          kappa = theta$kappa,
          t_values = as.numeric(t),
          distance_type = distance_type
        ))
      }

      theoretical_distance_profile_vmf(
        omega = omega,
        mu = theta$mu,
        kappa = theta$kappa,
        t_values = as.numeric(t),
        distance_type = distance_type
      )
    },
    normalize_data = normalize_vmf_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_vmf_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_vmf_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_vmf_theta(theta)
        ambient_dim <- length(theta$mu)
        profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
        n_u <- as.integer(control$vmf_profile_n_u %||% 4097L)
        l_max <- control$vmf_profile_l_max %||% NULL
        tail_tol <- as.numeric(control$vmf_profile_legendre_tail_tol %||% 1e-10)

        if (ambient_dim == 3L && identical(profile_method, "tabulated")) {
          return(distance_profile_vmf_s2_grid(
            omega_grid = omega_grid,
            mu = theta$mu,
            kappa = theta$kappa,
            t_grid = t_grid,
            distance_type = distance_type,
            n_u = n_u
          ))
        }
        if (ambient_dim == 3L && identical(profile_method, "legendre")) {
          return(distance_profile_vmf_s2_legendre_grid(
            omega_grid = omega_grid,
            mu = theta$mu,
            kappa = theta$kappa,
            t_grid = t_grid,
            distance_type = distance_type,
            l_max = l_max,
            tail_tol = tail_tol
          ))
        }

        NULL
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_vmf_theta(theta)
        ambient_dim <- length(theta$mu)
        profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
        n_u <- as.integer(control$vmf_profile_n_u %||% 4097L)
        l_max <- control$vmf_profile_l_max %||% NULL
        tail_tol <- as.numeric(control$vmf_profile_legendre_tail_tol %||% 1e-10)

        if (ambient_dim == 3L && identical(profile_method, "tabulated")) {
          return(distance_profile_vmf_s2_cvm_grid(
            X = data,
            mu = theta$mu,
            kappa = theta$kappa,
            n_u = n_u
          ))
        }
        if (ambient_dim == 3L && identical(profile_method, "legendre")) {
          return(distance_profile_vmf_s2_legendre_cvm_grid(
            X = data,
            mu = theta$mu,
            kappa = theta$kappa,
            l_max = l_max,
            tail_tol = tail_tol
          ))
        }

        NULL
      },
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        if (!identical(unknown_param, "xi")) {
          stop("The fast vMF multiplier bootstrap currently supports only `unknown_param = 'xi'`.")
        }

        prepare_vmf_fast_multiplier(
          data = data,
          theta_hat = theta_hat,
          spec = make_vmf_spec(distance_type = distance_type, unknown_param = unknown_param),
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          distance_type = distance_type
        )
      },
      distance_type = distance_type,
      unknown_param = unknown_param
    )
  )
}

normalize_jp_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) < 3L) {
    stop("JP data must be a non-empty matrix with at least three columns.")
  }
  if (any(!is.finite(data))) {
    stop("JP data must be finite.")
  }

  norms <- sqrt(rowSums(data^2))
  if (any(norms <= 0)) {
    stop("JP data rows must have strictly positive norm.")
  }

  data / norms
}

normalize_jp_theta <- function(theta, ambient_dim = NULL) {
  if (!is.list(theta)) {
    stop("JP theta must be a list containing `mu`, `kappa`, and `psi`.")
  }

  mu <- as.numeric(theta$mu)
  kappa <- as.numeric(theta$kappa)
  psi <- as.numeric(theta$psi)

  if (length(mu) < 3L || any(!is.finite(mu))) {
    stop("JP theta requires a finite vector `mu` of length at least 3.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("JP theta requires a nonnegative finite scalar `kappa`.")
  }
  if (length(psi) != 1L || !is.finite(psi)) {
    stop("JP theta requires a finite scalar `psi`.")
  }

  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 0) {
    stop("JP theta requires `mu` with strictly positive norm.")
  }
  mu <- mu / mu_norm

  if (!is.null(ambient_dim) && length(mu) != ambient_dim) {
    stop("JP theta has incompatible ambient dimension.")
  }

  q <- length(mu) - 1L
  if (q < 2L) {
    stop("JP theta currently supports only q >= 2.")
  }

  list(
    mu = mu,
    kappa = kappa,
    psi = psi,
    alpha = if (psi == 0) 0 else tanh(kappa * psi),
    beta = if (psi == 0) NA_real_ else 1 / psi,
    q = q
  )
}

fast_multiplier_small_circle_norm_constant <- function(kappa, nu) {
  root_kappa <- sqrt(as.numeric(kappa))
  (sqrt(pi) / (2 * root_kappa)) * (
    small_circle_erf(root_kappa * (1 - nu)) +
      small_circle_erf(root_kappa * (1 + nu))
  )
}

fast_multiplier_small_circle_dlogc_dkappa <- function(kappa, nu) {
  c_value <- fast_multiplier_small_circle_norm_constant(kappa, nu)
  -1 / (2 * kappa) +
    (
      (1 - nu) * exp(-kappa * (1 - nu)^2) +
        (1 + nu) * exp(-kappa * (1 + nu)^2)
    ) / (2 * kappa * c_value)
}

fast_multiplier_small_circle_dlogc_dnu <- function(kappa, nu) {
  c_value <- fast_multiplier_small_circle_norm_constant(kappa, nu)
  (
    exp(-kappa * (1 + nu)^2) -
      exp(-kappa * (1 - nu)^2)
  ) / c_value
}

fast_multiplier_small_circle_component_scores_natural <- function(s,
                                                                  kappa,
                                                                  nu) {
  cbind(
    -fast_multiplier_small_circle_dlogc_dkappa(kappa, nu) - (s - nu)^2,
    -fast_multiplier_small_circle_dlogc_dnu(kappa, nu) + 2 * kappa * (s - nu)
  )
}

fast_multiplier_small_circle_component_log_density <- function(s,
                                                               kappa,
                                                               nu) {
  -log(fast_multiplier_small_circle_norm_constant(kappa, nu)) -
    kappa * (s - nu)^2
}

fast_multiplier_cardioid_legendre <- function(z, k) {
  z <- as.numeric(z)
  switch(
    as.character(as.integer(k)),
    `1` = z,
    `2` = (3 * z^2 - 1) / 2,
    `3` = (5 * z^3 - 3 * z) / 2,
    `4` = (35 * z^4 - 30 * z^2 + 3) / 8,
    stop("Fast cardioid support currently covers only k = 1, 2, 3, 4.")
  )
}

fast_multiplier_cardioid_legendre_prime <- function(z, k) {
  z <- as.numeric(z)
  switch(
    as.character(as.integer(k)),
    `1` = rep.int(1, length(z)),
    `2` = 3 * z,
    `3` = (15 * z^2 - 3) / 2,
    `4` = (35 * z^3 - 15 * z) / 2,
    stop("Fast cardioid support currently covers only k = 1, 2, 3, 4.")
  )
}

fast_multiplier_jp_axial_expectations <- function(alpha,
                                                  beta,
                                                  quad_n = 1024L) {
  quad <- rotational_gauss_legendre(as.integer(quad_n))
  term <- 1 + alpha * quad$nodes
  if (any(!is.finite(term)) || any(term <= 0)) {
    stop("JP axial expectation quadrature encountered nonpositive support.")
  }

  weights <- quad$weights * term^beta
  denom <- sum(weights)
  if (!is.finite(denom) || denom <= 0) {
    stop("JP axial expectation quadrature produced a nonpositive normalizing constant.")
  }

  list(
    e_t_over_one_plus_alpha_t = sum(weights * (quad$nodes / term)) / denom,
    e_log_one_plus_alpha_t = sum(weights * log(term)) / denom
  )
}

prepare_jp_fast_multiplier <- function(spec,
                                       data,
                                       theta_hat,
                                       ks_prep = NULL,
                                       cvm_prep = NULL,
                                       control = list(),
                                       distance_type = "geodesic") {
  x <- normalize_jp_data(data, control)
  theta_hat <- normalize_jp_theta(theta_hat, ambient_dim = ncol(x))

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = length(theta_hat$mu),
    kappa = theta_hat$kappa,
    psi = theta_hat$psi,
    abs_kappa_psi_tol = as.numeric(control$jp_vmf_switch_abs_kappa_psi %||% jp_vmf_near_zero_abs_kappa_psi_default)
  )) {
    theta_vmf <- normalize_vmf_theta(list(mu = theta_hat$mu, kappa = theta_hat$kappa))
    return(prepare_vmf_fast_multiplier(
      data = x,
      theta_hat = theta_vmf,
      spec = make_vmf_spec(distance_type = distance_type, unknown_param = "xi"),
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control,
      distance_type = distance_type
    ))
  }

  jp_alpha_boundary_eps <- as.numeric(control$jp_fast_alpha_boundary_eps %||% 1e-8)
  if (abs(theta_hat$alpha) >= 1 - jp_alpha_boundary_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "jp_alpha_boundary",
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_
    ))
  }

  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(0, 0, theta_hat$kappa, theta_hat$psi)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_jp_theta(
      list(mu = mapped$mu, kappa = par[[3L]], psi = par[[4L]]),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_jp_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    alpha <- theta_state$alpha
    beta <- theta_state$beta
    term <- 1 + alpha * z
    axial <- fast_multiplier_jp_axial_expectations(
      alpha = alpha,
      beta = beta,
      quad_n = as.integer(control$jp_fast_axial_quad_n %||% 1024L)
    )

    coeff_mu <- alpha * beta / term
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))
    psi_alpha <- beta * (z / term - axial$e_t_over_one_plus_alpha_t)
    psi_beta <- log(term) - axial$e_log_one_plus_alpha_t

    cbind(
      score_mu,
      theta_state$psi * (1 - alpha^2) * psi_alpha,
      theta_state$kappa * (1 - alpha^2) * psi_alpha - (theta_state$psi^-2) * psi_beta
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_jp(
      n = n_aux,
      mu = theta_state$mu,
      kappa = theta_state$kappa,
      psi = theta_state$psi
    )
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

prepare_beta_mixture2_fast_multiplier <- function(spec,
                                                  data,
                                                  theta_hat,
                                                  ks_prep = NULL,
                                                  cvm_prep = NULL,
                                                  control = list()) {
  x <- normalize_beta_mixture2_data(data, control)
  theta_hat <- normalize_beta_mixture2_theta(theta_hat, ambient_dim = ncol(x))
  fast_shape_regular_eps <- as.numeric(control$beta_mixture2_fast_shape_regular_eps %||% 0)
  shape_values <- c(theta_hat$alpha1, theta_hat$beta1, theta_hat$alpha2, theta_hat$beta2)
  if (any(!is.finite(shape_values))) {
    stop("The fitted beta_mixture2 shape parameters are not finite.")
  }
  if (min(shape_values) <= 1 + fast_shape_regular_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "beta_mixture2_shape_nonregular"
    ))
  }
  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(
    0,
    0,
    stats::qlogis(theta_hat$weight1),
    log(theta_hat$alpha1),
    log(theta_hat$beta1),
    log(theta_hat$alpha2),
    log(theta_hat$beta2)
  )
  weight_eps <- as.numeric(control$beta_mixture2_weight_eps %||% 0.01)
  shape_lower <- as.numeric(control$beta_mixture2_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$beta_mixture2_shape_upper %||% 1e3)
  y_eps <- as.numeric(control$beta_mixture2_eps %||% 1e-12)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_beta_mixture2_theta(
      list(
        mu = mapped$mu,
        weight1 = rotational_bounded_weight(par[[3L]], weight_eps = weight_eps),
        alpha1 = rotational_positive_parameter(par[[4L]], lower = shape_lower, upper = shape_upper),
        beta1 = rotational_positive_parameter(par[[5L]], lower = shape_lower, upper = shape_upper),
        alpha2 = rotational_positive_parameter(par[[6L]], lower = shape_lower, upper = shape_upper),
        beta2 = rotational_positive_parameter(par[[7L]], lower = shape_lower, upper = shape_upper)
      ),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_beta_mixture2_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    y <- rotational_clamp_unit_interval((z + 1) / 2, eps = y_eps)
    log_b1 <- stats::dbeta(y, theta_state$alpha1, theta_state$beta1, log = TRUE)
    log_b2 <- stats::dbeta(y, theta_state$alpha2, theta_state$beta2, log = TRUE)
    log_m <- rotational_logsumexp2(
      log(theta_state$weight1) + log_b1,
      log1p(-theta_state$weight1) + log_b2
    )
    r1 <- exp(log(theta_state$weight1) + log_b1 - log_m)
    r2 <- 1 - r1

    g1 <- (theta_state$alpha1 - 1) / y - (theta_state$beta1 - 1) / (1 - y)
    g2 <- (theta_state$alpha2 - 1) / y - (theta_state$beta2 - 1) / (1 - y)
    coeff_mu <- 0.5 * (r1 * g1 + r2 * g2)
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))

    cbind(
      score_mu,
      r1 - theta_state$weight1,
      theta_state$alpha1 * r1 * (
        log(y) - digamma(theta_state$alpha1) + digamma(theta_state$alpha1 + theta_state$beta1)
      ),
      theta_state$beta1 * r1 * (
        log1p(-y) - digamma(theta_state$beta1) + digamma(theta_state$alpha1 + theta_state$beta1)
      ),
      theta_state$alpha2 * r2 * (
        log(y) - digamma(theta_state$alpha2) + digamma(theta_state$alpha2 + theta_state$beta2)
      ),
      theta_state$beta2 * r2 * (
        log1p(-y) - digamma(theta_state$beta2) + digamma(theta_state$alpha2 + theta_state$beta2)
      )
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_beta_mixture2(
      n = n_aux,
      mu = theta_state$mu,
      weight1 = theta_state$weight1,
      alpha1 = theta_state$alpha1,
      beta1 = theta_state$beta1,
      alpha2 = theta_state$alpha2,
      beta2 = theta_state$beta2
    )
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

prepare_uniform_beta_mixture_fast_multiplier <- function(spec,
                                                         data,
                                                         theta_hat,
                                                         ks_prep = NULL,
                                                         cvm_prep = NULL,
                                                         control = list()) {
  x <- normalize_uniform_beta_mixture_data(data, control)
  theta_hat <- normalize_uniform_beta_mixture_theta(theta_hat, ambient_dim = ncol(x))
  fast_shape_regular_eps <- as.numeric(control$uniform_beta_mixture_fast_shape_regular_eps %||% 0)
  shape_values <- c(theta_hat$alpha, theta_hat$beta)
  if (any(!is.finite(shape_values))) {
    stop("The fitted uniform_beta_mixture shape parameters are not finite.")
  }
  if (min(shape_values) <= 1 + fast_shape_regular_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "uniform_beta_mixture_shape_nonregular"
    ))
  }

  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(
    0,
    0,
    stats::qlogis(theta_hat$weight_uniform),
    log(theta_hat$alpha),
    log(theta_hat$beta)
  )
  weight_eps <- as.numeric(control$uniform_beta_mixture_weight_eps %||% 0.01)
  shape_lower <- as.numeric(control$uniform_beta_mixture_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$uniform_beta_mixture_shape_upper %||% 1e3)
  y_eps <- as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_uniform_beta_mixture_theta(
      list(
        mu = mapped$mu,
        weight_uniform = rotational_bounded_weight(par[[3L]], weight_eps = weight_eps),
        alpha = rotational_positive_parameter(par[[4L]], lower = shape_lower, upper = shape_upper),
        beta = rotational_positive_parameter(par[[5L]], lower = shape_lower, upper = shape_upper)
      ),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_uniform_beta_mixture_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    y <- rotational_clamp_unit_interval((z + 1) / 2, eps = y_eps)
    log_b <- stats::dbeta(y, theta_state$alpha, theta_state$beta, log = TRUE)
    log_m <- rotational_logsumexp2(
      log(theta_state$weight_uniform),
      log1p(-theta_state$weight_uniform) + log_b
    )
    r_beta <- exp(log1p(-theta_state$weight_uniform) + log_b - log_m)
    r_uniform <- 1 - r_beta

    g_beta <- (theta_state$alpha - 1) / y - (theta_state$beta - 1) / (1 - y)
    coeff_mu <- 0.5 * r_beta * g_beta
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))

    cbind(
      score_mu,
      r_uniform - theta_state$weight_uniform,
      theta_state$alpha * r_beta * (
        log(y) - digamma(theta_state$alpha) + digamma(theta_state$alpha + theta_state$beta)
      ),
      theta_state$beta * r_beta * (
        log1p(-y) - digamma(theta_state$beta) + digamma(theta_state$alpha + theta_state$beta)
      )
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_uniform_beta_mixture(
      n = n_aux,
      mu = theta_state$mu,
      weight_uniform = theta_state$weight_uniform,
      alpha = theta_state$alpha,
      beta = theta_state$beta
    )
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

fit_jp_theta <- function(data,
                         weights = NULL,
                         null,
                         control = list()) {
  x <- normalize_jp_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_jp_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (ncol(x) != 3L) {
    stop("The JP composite adapter currently supports only S^2, i.e. data with three columns.")
  }

  jp_control <- control
  jp_control$jp_data_already_normalized <- TRUE

  normalize_jp_theta(
    jp_mle_s2_weighted(
      data = x,
      weights = weights,
      control = jp_control
    ),
    ambient_dim = ncol(x)
  )
}

make_jp_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("jp_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_jp_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_jp_data(data, control)
      omega_matrix <- normalize_jp_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- x %*% t(omega_matrix)
      dot_products <- pmax(pmin(dot_products, 1), -1)

      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_jp_theta(theta)
      distance_profile_jp(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        kappa = theta$kappa,
        psi = theta$psi,
        distance_type = distance_type
      )
    },
    normalize_data = normalize_jp_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_jp_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_jp_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_jp_theta(theta)
        profile_method <- tolower(as.character(control$jp_profile_method %||% "tabulated"))
        n_u <- as.integer(control$jp_profile_n_u %||% 1025L)
        n_delta <- as.integer(control$jp_profile_n_delta %||% 257L)

        if (!identical(profile_method, "tabulated")) {
          return(NULL)
        }
        if (theta$psi == 0 && length(theta$mu) != 3L) {
          return(NULL)
        }

        distance_profile_jp_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          kappa = theta$kappa,
          psi = theta$psi,
          t_grid = t_grid,
          distance_type = distance_type,
          n_u = n_u,
          n_delta = n_delta
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_jp_theta(theta)
        profile_method <- tolower(as.character(control$jp_profile_method %||% "tabulated"))
        n_u <- as.integer(control$jp_profile_n_u %||% 1025L)
        n_delta <- as.integer(control$jp_profile_n_delta %||% 257L)

        if (!identical(profile_method, "tabulated")) {
          return(NULL)
        }
        if (theta$psi == 0 && length(theta$mu) != 3L) {
          return(NULL)
        }

        distance_profile_jp_cvm_grid(
          X = data,
          mu = theta$mu,
          kappa = theta$kappa,
          psi = theta$psi,
          n_u = n_u,
          n_delta = n_delta
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_jp_fast_multiplier(
          spec = make_jp_spec(distance_type = distance_type),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          distance_type = distance_type
        )
      }
    )
  )
}

normalize_hvmf_data <- function(data, control = list()) {
  tol <- as.numeric(control$hvmf_tol %||% 1e-10)
  normalize_hvmf_h2_data(data, tol = tol)
}

normalize_hvmf_theta <- function(theta, control = list()) {
  tol <- as.numeric(control$hvmf_tol %||% 1e-10)

  if (!is.list(theta)) {
    stop("HvMF theta must be a list containing `mu` and `kappa`.")
  }

  mu <- theta$mu
  kappa <- as.numeric(theta$kappa)

  mu_matrix <- normalize_hvmf_h2_data(mu, tol = tol)
  mu <- as.numeric(mu_matrix[1L, , drop = TRUE])

  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("HvMF theta requires a strictly positive finite scalar `kappa`.")
  }

  sinh_chi <- sqrt(sum(mu[-1L]^2))
  chi <- asinh(sinh_chi)
  theta_angle <- atan2(mu[[3L]], mu[[2L]])
  theta_deg <- (theta_angle * 180 / pi) %% 360

  list(
    mu = mu,
    kappa = kappa,
    chi = chi,
    sinh_chi = sinh_chi,
    theta = theta_angle,
    theta_deg = theta_deg,
    q = 2L
  )
}

fit_hvmf_theta <- function(data,
                           weights = NULL,
                           null,
                           unknown_param = "both",
                           control = list()) {
  x <- normalize_hvmf_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_hvmf_theta(null$theta, control = control))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (!is.null(unknown_param) && !identical(unknown_param, "both")) {
    stop("The HvMF adapter currently supports only `unknown_param = 'both'`.")
  }

  tol <- as.numeric(control$hvmf_tol %||% 1e-10)
  fit <- hvmf_mle_h2(x, weights = weights, tol = tol)

  normalize_hvmf_theta(
    list(mu = fit$mu, kappa = fit$kappa),
    control = control
  )
}

prepare_hvmf_fast_multiplier <- function(spec,
                                         data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list(),
                                         unknown_param = "both") {
  if (!identical(unknown_param, "both")) {
    stop("The fast HvMF multiplier bootstrap currently supports only `unknown_param = 'both'`.")
  }

  x <- normalize_hvmf_data(data, control)
  theta_hat <- normalize_hvmf_theta(theta_hat, control = control)
  par0 <- c(theta_hat$chi, theta_hat$theta, theta_hat$kappa)

  state_from_par <- function(par) {
    chi <- as.numeric(par[[1L]])
    theta_angle <- as.numeric(par[[2L]])
    kappa <- as.numeric(par[[3L]])
    mu <- c(
      cosh(chi),
      sinh(chi) * cos(theta_angle),
      sinh(chi) * sin(theta_angle)
    )
    normalize_hvmf_theta(list(mu = mu, kappa = kappa), control = control)
  }

  dmu_dchi <- function(chi, theta_angle) {
    c(
      sinh(chi),
      cosh(chi) * cos(theta_angle),
      cosh(chi) * sin(theta_angle)
    )
  }

  dmu_dtheta <- function(chi, theta_angle) {
    c(
      0,
      -sinh(chi) * sin(theta_angle),
      sinh(chi) * cos(theta_angle)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_hvmf_data(sample, control)
    chi <- theta_state$chi
    theta_angle <- theta_state$theta
    mu <- theta_state$mu
    kappa <- theta_state$kappa
    dchi_vec <- dmu_dchi(chi, theta_angle)
    dtheta_vec <- dmu_dtheta(chi, theta_angle)
    mink_mu <- apply(sample, 1L, function(row) minkowski_inner_product(row, mu))
    score_loc_chi <- kappa * apply(sample, 1L, function(row) minkowski_inner_product(row, dchi_vec))
    score_loc_theta <- kappa * apply(sample, 1L, function(row) minkowski_inner_product(row, dtheta_vec))
    cbind(
      score_loc_chi,
      score_loc_theta,
      1 / kappa + 1 + mink_mu
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    rhvmf_h2_polar(n_aux, mu = theta_state$mu, kappa = theta_state$kappa)
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

make_hvmf_spec <- function(unknown_param = "both") {
  new_model_spec(
    name = "hvmf_geodesic",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_hvmf_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_hvmf_data(data, control)
      omega_matrix <- normalize_hvmf_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      t(hvmf_distance_matrix(omega_matrix, x))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_hvmf_theta(theta, control = control)
      theoretical_distance_profile_hvmf(
        omega = omega,
        mu = theta$mu,
        kappa = theta$kappa,
        t_values = as.numeric(t)
      )
    },
    normalize_data = normalize_hvmf_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_hvmf_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_hvmf_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_hvmf_theta(theta, control = control)
        profile_method <- tolower(as.character(control$hvmf_profile_method %||% "tabulated"))
        n_y <- as.integer(control$hvmf_profile_n_y %||% control$hvmf_profile_grid_size %||% 4097L)

        if (!profile_method %in% c("exact", "tabulated")) {
          stop("`control$hvmf_profile_method` must be either 'exact' or 'tabulated'.")
        }

        if (identical(profile_method, "tabulated")) {
          return(hvmf_cvm_profile_matrix_tabulated(
            data = data,
            theta = theta,
            grid_size = n_y,
            distance_matrix = distance_matrix,
            tol = as.numeric(control$hvmf_tol %||% 1e-10)
          ))
        }

        NULL
      },
      distance_type = "geodesic",
      unknown_param = unknown_param,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_hvmf_fast_multiplier(
          spec = make_hvmf_spec(unknown_param = unknown_param),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          unknown_param = unknown_param
        )
      }
    )
  )
}

normalize_spherical_cauchy_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Spherical Cauchy data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_spherical_cauchy_theta <- function(theta, ambient_dim = 3L) {
  theta <- spherical_cauchy_normalize_theta(theta, ambient_dim = ambient_dim)
  if (theta$ambient_dim != 3L) {
    stop("Spherical Cauchy theta currently supports only S^2.")
  }
  theta
}

fit_spherical_cauchy_theta <- function(data,
                                       weights = NULL,
                                       null,
                                       control = list()) {
  x <- normalize_spherical_cauchy_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_spherical_cauchy_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- spherical_cauchy_mle_s2_weighted(
    data = x,
    weights = weights,
    control = control
  )

  normalize_spherical_cauchy_theta(fit, ambient_dim = ncol(x))
}

prepare_spherical_cauchy_fast_multiplier <- function(spec,
                                                     data,
                                                     theta_hat,
                                                     ks_prep = NULL,
                                                     cvm_prep = NULL,
                                                     control = list()) {
  x <- normalize_spherical_cauchy_data(data, control)
  theta_hat <- normalize_spherical_cauchy_theta(theta_hat, ambient_dim = ncol(x))
  rho_zero_eps <- as.numeric(control$spherical_cauchy_fast_boundary_eps %||%
    control$spherical_cauchy_fast_zero_eps %||% 1e-8)
  if (theta_hat$rho <= rho_zero_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "spherical_cauchy_rho_zero_nonidentification",
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_
    ))
  }
  par0 <- as.numeric(theta_hat$phi)

  state_from_par <- function(par) {
    normalize_spherical_cauchy_theta(list(phi = par), ambient_dim = ncol(x))
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_spherical_cauchy_data(sample, control)
    phi <- as.numeric(theta_state$phi)
    phi_norm_sq <- sum(phi^2)
    denom <- as.numeric(1 - 2 * (sample %*% phi) + phi_norm_sq)
    const_term <- matrix(
      rep(-2 * phi / (1 - phi_norm_sq), each = nrow(sample)),
      nrow = nrow(sample),
      ncol = length(phi)
    )
    phi_matrix <- matrix(
      rep(phi, each = nrow(sample)),
      nrow = nrow(sample),
      ncol = length(phi)
    )
    const_term + 3 * sample / denom - sweep(phi_matrix, 1L, 3 / denom, FUN = "*")
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_spherical_cauchy(
      n = n_aux,
      mu = theta_state$mu,
      rho = theta_state$rho
    )
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

make_spherical_cauchy_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("spherical_cauchy_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_spherical_cauchy_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_spherical_cauchy_data(data, control)
      omega_matrix <- normalize_spherical_cauchy_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- x %*% t(omega_matrix)
      dot_products <- pmin(pmax(dot_products, -1), 1)

      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_spherical_cauchy_theta(theta)
      l_max <- control$spherical_cauchy_profile_l_max %||% NULL
      max_l_max <- control$spherical_cauchy_profile_max_l_max %||% NULL
      if (!is.null(l_max)) {
        l_max <- as.integer(l_max)
      }
      if (!is.null(max_l_max)) {
        max_l_max <- as.integer(max_l_max)
      }
      distance_profile_spherical_cauchy(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        rho = theta$rho,
        distance_type = distance_type,
        tail_tol = as.numeric(control$spherical_cauchy_profile_tol %||% 1e-10),
        l_max = l_max,
        max_l_max = max_l_max,
        warn = isTRUE(control$spherical_cauchy_profile_warn %||% FALSE)
      )
    },
    normalize_data = normalize_spherical_cauchy_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_spherical_cauchy_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_spherical_cauchy_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_spherical_cauchy_theta(theta)

        l_max <- control$spherical_cauchy_profile_l_max %||% NULL
        max_l_max <- control$spherical_cauchy_profile_max_l_max %||% NULL
        if (!is.null(l_max)) {
          l_max <- as.integer(l_max)
        }
        if (!is.null(max_l_max)) {
          max_l_max <- as.integer(max_l_max)
        }

        distance_profile_spherical_cauchy_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          rho = theta$rho,
          t_grid = t_grid,
          distance_type = distance_type,
          tail_tol = as.numeric(control$spherical_cauchy_profile_tol %||% 1e-10),
          l_max = l_max,
          max_l_max = max_l_max,
          warn = isTRUE(control$spherical_cauchy_profile_warn %||% FALSE)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_spherical_cauchy_theta(theta)

        l_max <- control$spherical_cauchy_profile_l_max %||% NULL
        max_l_max <- control$spherical_cauchy_profile_max_l_max %||% NULL
        if (!is.null(l_max)) {
          l_max <- as.integer(l_max)
        }
        if (!is.null(max_l_max)) {
          max_l_max <- as.integer(max_l_max)
        }

        distance_profile_spherical_cauchy_cvm_grid(
          data = data,
          mu = theta$mu,
          rho = theta$rho,
          distance_matrix = distance_matrix,
          distance_type = distance_type,
          tail_tol = as.numeric(control$spherical_cauchy_profile_tol %||% 1e-10),
          l_max = l_max,
          max_l_max = max_l_max,
          warn = isTRUE(control$spherical_cauchy_profile_warn %||% FALSE)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_spherical_cauchy_fast_multiplier(
          spec = make_spherical_cauchy_spec(distance_type = distance_type),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control
        )
      }
    )
  )
}

normalize_beta_mixture2_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Rotational beta-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_beta_mixture2_theta <- function(theta, ambient_dim = 3L) {
  beta_mixture2_normalize_theta(theta, ambient_dim = ambient_dim)
}

fit_beta_mixture2_theta <- function(data,
                                               weights = NULL,
                                               null,
                                               control = list()) {
  x <- normalize_beta_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_beta_mixture2_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- beta_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_beta_mixture2_theta(fit, ambient_dim = ncol(x))
}

make_beta_mixture2_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("beta_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_beta_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_beta_mixture2_data(data, control)
      omega_matrix <- normalize_beta_mixture2_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- pmin(pmax(x %*% t(omega_matrix), -1), 1)
      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_beta_mixture2_theta(theta)
      distance_profile_beta_mixture2(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
        distance_type = distance_type,
        method = control$beta_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$beta_mixture2_L_max %||% 150L),
        quad_n = as.integer(control$beta_mixture2_quad_n %||% 1000L),
        tol = as.numeric(control$beta_mixture2_tol %||% 1e-6),
        validate_against_integral = isTRUE(control$beta_mixture2_validate_against_integral),
        validation_tol = as.numeric(control$beta_mixture2_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_beta_mixture2_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_beta_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_beta_mixture2_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_beta_mixture2_theta(theta)
        distance_profile_beta_mixture2_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$beta_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$beta_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$beta_mixture2_quad_n %||% 1000L),
          tol = as.numeric(control$beta_mixture2_tol %||% 1e-6)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_beta_mixture2_theta(theta)
        distance_profile_beta_mixture2_cvm_grid(
          X = data,
          mu = theta$mu,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2,
          method = control$beta_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$beta_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$beta_mixture2_quad_n %||% 1000L),
          tol = as.numeric(control$beta_mixture2_tol %||% 1e-6)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_beta_mixture2_fast_multiplier(
          spec = make_beta_mixture2_spec(distance_type = distance_type),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control
        )
      }
    )
  )
}

normalize_uniform_beta_mixture_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Rotational uniform-beta-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_uniform_beta_mixture_theta <- function(theta, ambient_dim = 3L) {
  uniform_beta_mixture_normalize_theta(theta, ambient_dim = ambient_dim)
}

fit_uniform_beta_mixture_theta <- function(data,
                                           weights = NULL,
                                           null,
                                           control = list()) {
  x <- normalize_uniform_beta_mixture_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_uniform_beta_mixture_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- uniform_beta_mixture_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_uniform_beta_mixture_theta(fit, ambient_dim = ncol(x))
}

make_uniform_beta_mixture_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("uniform_beta_mixture_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_uniform_beta_mixture_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_uniform_beta_mixture_data(data, control)
      omega_matrix <- normalize_uniform_beta_mixture_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- pmin(pmax(x %*% t(omega_matrix), -1), 1)
      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_uniform_beta_mixture_theta(theta)
      distance_profile_uniform_beta_mixture(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta,
        distance_type = distance_type,
        method = control$uniform_beta_mixture_profile_method %||% "legendre",
        l_max = as.integer(control$uniform_beta_mixture_L_max %||% 150L),
        quad_n = as.integer(control$uniform_beta_mixture_quad_n %||% 1000L),
        tol = as.numeric(control$uniform_beta_mixture_tol %||% 1e-6),
        eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12),
        validate_against_integral = isTRUE(control$uniform_beta_mixture_validate_against_integral),
        validation_tol = as.numeric(control$uniform_beta_mixture_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_uniform_beta_mixture_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_uniform_beta_mixture_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_uniform_beta_mixture_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_uniform_beta_mixture_theta(theta)
        distance_profile_uniform_beta_mixture_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$uniform_beta_mixture_profile_method %||% "legendre",
          l_max = as.integer(control$uniform_beta_mixture_L_max %||% 150L),
          quad_n = as.integer(control$uniform_beta_mixture_quad_n %||% 1000L),
          tol = as.numeric(control$uniform_beta_mixture_tol %||% 1e-6),
          eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_uniform_beta_mixture_theta(theta)
        distance_profile_uniform_beta_mixture_cvm_grid(
          X = data,
          mu = theta$mu,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          distance_matrix = distance_matrix,
          method = control$uniform_beta_mixture_profile_method %||% "legendre",
          l_max = as.integer(control$uniform_beta_mixture_L_max %||% 150L),
          quad_n = as.integer(control$uniform_beta_mixture_quad_n %||% 1000L),
          tol = as.numeric(control$uniform_beta_mixture_tol %||% 1e-6),
          eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_uniform_beta_mixture_fast_multiplier(
          spec = make_uniform_beta_mixture_spec(distance_type = distance_type),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control
        )
      }
    )
  )
}

normalize_logitnormal_mixture2_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Rotational logit-normal-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_logitnormal_mixture2_theta <- function(theta, ambient_dim = 3L) {
  logitnormal_mixture2_normalize_theta(theta, ambient_dim = ambient_dim)
}

fit_logitnormal_mixture2_theta <- function(data,
                                                      weights = NULL,
                                                      null,
                                                      control = list()) {
  x <- normalize_logitnormal_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_logitnormal_mixture2_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- logitnormal_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_logitnormal_mixture2_theta(fit, ambient_dim = ncol(x))
}

make_logitnormal_mixture2_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("logitnormal_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_logitnormal_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_logitnormal_mixture2_data(data, control)
      omega_matrix <- normalize_logitnormal_mixture2_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- pmin(pmax(x %*% t(omega_matrix), -1), 1)
      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_logitnormal_mixture2_theta(theta)
      distance_profile_logitnormal_mixture2(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        distance_type = distance_type,
        method = control$logitnormal_mixture2_profile_method %||% "integral",
        l_max = as.integer(control$logitnormal_mixture2_L_max %||% 150L),
        quad_n = as.integer(control$logitnormal_mixture2_quad_n %||% 1000L),
        tol = as.numeric(control$logitnormal_mixture2_tol %||% 1e-6),
        eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12),
        validate_against_integral = isTRUE(control$logitnormal_mixture2_validate_against_integral),
        validation_tol = as.numeric(control$logitnormal_mixture2_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_logitnormal_mixture2_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_logitnormal_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_logitnormal_mixture2_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_logitnormal_mixture2_theta(theta)
        distance_profile_logitnormal_mixture2_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$logitnormal_mixture2_profile_method %||% "integral",
          l_max = as.integer(control$logitnormal_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$logitnormal_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$logitnormal_mixture2_tol %||% 1e-6),
          eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_logitnormal_mixture2_theta(theta)
        distance_profile_logitnormal_mixture2_cvm_grid(
          X = data,
          mu = theta$mu,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          method = control$logitnormal_mixture2_profile_method %||% "integral",
          l_max = as.integer(control$logitnormal_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$logitnormal_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$logitnormal_mixture2_tol %||% 1e-6),
          eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE
    )
  )
}

normalize_axial_truncnorm_mixture2_data <- function(data, control = list()) {
  if (is.matrix(data) || is.data.frame(data)) {
    if (ncol(as.matrix(data)) != 1L) {
      stop("Axial truncated-normal-mixture data must be a numeric vector or one-column matrix/data frame.")
    }
    data <- as.matrix(data)[, 1L]
  }

  z <- as.numeric(data)
  if (length(z) == 0L) {
    stop("`data` cannot be empty.")
  }
  if (any(!is.finite(z))) {
    stop("Axial truncated-normal-mixture data must be finite.")
  }
  if (any(z < -1 - 1e-12 | z > 1 + 1e-12)) {
    stop("Axial truncated-normal-mixture data must lie in [-1, 1].")
  }

  pmin(1, pmax(-1, z))
}

clip_axial_truncnorm_prob <- function(x, eps = 1e-10) {
  pmin(1 - eps, pmax(eps, as.numeric(x)))
}

clip_axial_truncnorm_kappa <- function(x, min_value = 1e-8, max_value = 1e8) {
  pmin(max_value, pmax(min_value, as.numeric(x)))
}

normalize_axial_truncnorm_mixture2_theta <- function(theta) {
  if (!is.list(theta)) {
    stop("Axial truncated-normal-mixture theta must be a list with entries `pi`, `kappa1`, `nu1`, `kappa2`, `nu2`.")
  }

  pi <- clip_axial_truncnorm_prob(theta$pi)
  kappa1 <- as.numeric(theta$kappa1)
  kappa2 <- as.numeric(theta$kappa2)
  nu1 <- clip_axial_truncnorm_prob(theta$nu1, eps = 1e-10)
  nu2 <- clip_axial_truncnorm_prob(theta$nu2, eps = 1e-10)

  if (length(kappa1) != 1L || !is.finite(kappa1) || kappa1 <= 0) {
    stop("Axial truncated-normal-mixture theta requires a strictly positive finite scalar `kappa1`.")
  }
  if (length(kappa2) != 1L || !is.finite(kappa2) || kappa2 <= 0) {
    stop("Axial truncated-normal-mixture theta requires a strictly positive finite scalar `kappa2`.")
  }

  list(
    pi = pi,
    kappa1 = clip_axial_truncnorm_kappa(kappa1),
    nu1 = nu1,
    kappa2 = clip_axial_truncnorm_kappa(kappa2),
    nu2 = nu2
  )
}

axial_truncnorm_interval_mass <- function(kappa, mean_value, lower, upper) {
  if (!is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be strictly positive and finite.")
  }

  sqrt_kappa <- sqrt(kappa)
  upper_std <- sqrt(2) * sqrt_kappa * (upper - mean_value)
  lower_std <- sqrt(2) * sqrt_kappa * (lower - mean_value)

  stats::pnorm(upper_std) - stats::pnorm(lower_std)
}

axial_truncnorm_logdiffexp <- function(log_x, log_y) {
  if (any(log_y > log_x, na.rm = TRUE)) {
    stop("`logdiffexp` requires `log_x >= log_y`.")
  }
  log_x + log1p(-exp(log_y - log_x))
}

axial_truncnorm_log_interval_mass <- function(kappa, mean_value, lower, upper) {
  if (!is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be strictly positive and finite.")
  }

  sqrt_kappa <- sqrt(kappa)
  upper_std <- sqrt(2) * sqrt_kappa * (upper - mean_value)
  lower_std <- sqrt(2) * sqrt_kappa * (lower - mean_value)

  if (lower_std >= 0) {
    log_tail_lower <- stats::pnorm(lower_std, lower.tail = FALSE, log.p = TRUE)
    log_tail_upper <- stats::pnorm(upper_std, lower.tail = FALSE, log.p = TRUE)

    if (!is.finite(log_tail_lower)) {
      return(log_tail_lower)
    }
    if (!is.finite(log_tail_upper)) {
      return(log_tail_lower)
    }

    return(axial_truncnorm_logdiffexp(log_tail_lower, log_tail_upper))
  }

  if (upper_std <= 0) {
    log_cdf_upper <- stats::pnorm(upper_std, log.p = TRUE)
    log_cdf_lower <- stats::pnorm(lower_std, log.p = TRUE)

    if (!is.finite(log_cdf_upper)) {
      return(log_cdf_upper)
    }
    if (!is.finite(log_cdf_lower)) {
      return(log_cdf_upper)
    }

    return(axial_truncnorm_logdiffexp(log_cdf_upper, log_cdf_lower))
  }

  left_tail <- stats::pnorm(lower_std)
  right_tail <- stats::pnorm(upper_std, lower.tail = FALSE)
  mass <- 1 - left_tail - right_tail

  if (!is.finite(mass) || mass <= 0) {
    mass <- axial_truncnorm_interval_mass(kappa, mean_value, lower, upper)
  }
  if (!is.finite(mass) || mass <= 0) {
    stop("Failed to compute a positive interval mass for the axial truncated-normal component.")
  }

  log(mass)
}

axial_truncnorm_log_normconst <- function(kappa, mean_value) {
  log_mass <- axial_truncnorm_log_interval_mass(kappa, mean_value, lower = -1, upper = 1)
  if (!is.finite(log_mass)) {
    stop("Failed to compute a positive normalizing constant for the axial truncated-normal component.")
  }

  0.5 * log(pi / kappa) + log_mass
}

axial_truncnorm_component_log_density <- function(z, kappa, mean_value) {
  -(kappa * (z - mean_value)^2) - axial_truncnorm_log_normconst(kappa, mean_value)
}

axial_truncnorm_component_cdf <- function(z, kappa, mean_value) {
  z_clipped <- pmin(1, pmax(-1, as.numeric(z)))
  log_denom <- axial_truncnorm_log_interval_mass(kappa, mean_value, lower = -1, upper = 1)
  log_numer <- vapply(z_clipped, function(one_z) {
    if (one_z <= -1) {
      return(-Inf)
    }
    axial_truncnorm_log_interval_mass(kappa, mean_value, lower = -1, upper = one_z)
  }, numeric(1))
  out <- exp(log_numer - log_denom)
  out[z <= -1] <- 0
  out[z >= 1] <- 1
  pmin(1, pmax(0, out))
}

axial_truncnorm_logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

axial_truncnorm_mixture_log_density <- function(z, theta) {
  theta <- normalize_axial_truncnorm_mixture2_theta(theta)
  log_left <- log(theta$pi) + axial_truncnorm_component_log_density(z, theta$kappa1, theta$nu1)
  log_right <- log1p(-theta$pi) + axial_truncnorm_component_log_density(z, theta$kappa2, -theta$nu2)
  axial_truncnorm_logsumexp2(log_left, log_right)
}

axial_truncnorm_mixture_cdf <- function(z, theta) {
  theta <- normalize_axial_truncnorm_mixture2_theta(theta)
  theta$pi * axial_truncnorm_component_cdf(z, theta$kappa1, theta$nu1) +
    (1 - theta$pi) * axial_truncnorm_component_cdf(z, theta$kappa2, -theta$nu2)
}

axial_truncnorm_mixture_distance_profile <- function(omega, t_values, theta) {
  theta <- normalize_axial_truncnorm_mixture2_theta(theta)
  omega <- as.numeric(omega)
  if (length(omega) != 1L || !is.finite(omega) || omega < -1 || omega > 1) {
    stop("`omega` must be a finite scalar in [-1, 1].")
  }
  t_values <- pmax(0, as.numeric(t_values))
  upper <- pmin(1, omega + t_values)
  lower <- pmax(-1, omega - t_values)
  axial_truncnorm_mixture_cdf(upper, theta) - axial_truncnorm_mixture_cdf(lower, theta)
}

axial_truncnorm_mixture_distance_profile_grid <- function(omega_grid, t_grid, theta) {
  omega_grid <- as.numeric(omega_grid)
  t_grid <- pmax(0, as.numeric(t_grid))

  profile_values <- vapply(omega_grid, function(omega) {
    axial_truncnorm_mixture_distance_profile(omega = omega, t_values = t_grid, theta = theta)
  }, numeric(length(t_grid)))

  t(profile_values)
}

axial_truncnorm_encode_theta <- function(theta) {
  theta <- normalize_axial_truncnorm_mixture2_theta(theta)
  c(
    qlogis(theta$pi),
    log(theta$kappa1),
    qlogis(theta$nu1),
    log(theta$kappa2),
    qlogis(theta$nu2)
  )
}

axial_truncnorm_decode_theta <- function(par) {
  if (length(par) != 5L) {
    stop("Internal axial truncated-normal-mixture parameter vector must have length 5.")
  }

  nu_eps <- 1e-10
  log_kappa_min <- log(1e-8)
  log_kappa_max <- log(1e8)
  list(
    pi = clip_axial_truncnorm_prob(stats::plogis(par[[1L]])),
    kappa1 = clip_axial_truncnorm_kappa(exp(pmin(pmax(par[[2L]], log_kappa_min), log_kappa_max))),
    nu1 = clip_axial_truncnorm_prob(stats::plogis(par[[3L]]), eps = nu_eps),
    kappa2 = clip_axial_truncnorm_kappa(exp(pmin(pmax(par[[4L]], log_kappa_min), log_kappa_max))),
    nu2 = clip_axial_truncnorm_prob(stats::plogis(par[[5L]]), eps = nu_eps)
  )
}

axial_truncnorm_default_theta_start <- function(data,
                                                weights = NULL,
                                                prob_eps = 1e-4,
                                                kappa_min = 1e-8,
                                                kappa_max = 1e8) {
  z <- normalize_axial_truncnorm_mixture2_data(data)
  obs_weights <- if (is.null(weights)) {
    rep.int(1 / length(z), length(z))
  } else {
    normalize_probability_weights(weights, length(z))
  }

  weighted_mean <- function(values, value_weights, fallback) {
    if (length(values) == 0L || sum(value_weights) <= 0) {
      return(fallback)
    }
    sum(value_weights * values) / sum(value_weights)
  }

  weighted_variance <- function(values, value_weights, center, fallback) {
    if (length(values) <= 1L || sum(value_weights) <= 0) {
      return(fallback)
    }
    sum(value_weights * (values - center)^2) / sum(value_weights)
  }

  moment_kappa <- function(values, value_weights, center, fallback_var = 0.05^2) {
    variance_value <- weighted_variance(values, value_weights, center = center, fallback = fallback_var)
    clip_axial_truncnorm_kappa(1 / (2 * max(variance_value, 1 / (2 * kappa_max))), min_value = kappa_min, max_value = kappa_max)
  }

  north_idx <- z >= 0
  south_idx <- !north_idx

  north_mass <- sum(obs_weights[north_idx])
  south_mass <- sum(obs_weights[south_idx])
  total_mass <- north_mass + south_mass
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Failed to construct a valid weighted start for the axial truncated-normal mixture.")
  }

  pi0 <- clip_axial_truncnorm_prob(north_mass / total_mass, eps = prob_eps)

  z_pos <- z[north_idx]
  w_pos <- obs_weights[north_idx]
  z_neg <- -z[south_idx]
  w_neg <- obs_weights[south_idx]

  fallback_center <- clip_axial_truncnorm_prob(sum(obs_weights * abs(z)), eps = prob_eps)
  nu1 <- clip_axial_truncnorm_prob(
    weighted_mean(z_pos, w_pos, fallback = fallback_center),
    eps = prob_eps
  )
  nu2 <- clip_axial_truncnorm_prob(
    weighted_mean(z_neg, w_neg, fallback = fallback_center),
    eps = prob_eps
  )

  list(
    pi = pi0,
    kappa1 = moment_kappa(z_pos, w_pos, center = nu1),
    nu1 = nu1,
    kappa2 = moment_kappa(z_neg, w_neg, center = nu2),
    nu2 = nu2
  )
}

fit_axial_truncnorm_mixture2_theta <- function(data,
                                                weights = NULL,
                                                null,
                                                control = list()) {
  z <- normalize_axial_truncnorm_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_axial_truncnorm_mixture2_theta(null$theta))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  obs_weights <- if (is.null(weights)) rep.int(1, length(z)) else as.numeric(weights)
  if (length(obs_weights) != length(z) || any(!is.finite(obs_weights)) || any(obs_weights < 0) || sum(obs_weights) <= 0) {
    stop("`weights` must be NULL or a nonnegative finite vector with positive sum and length equal to `data`.")
  }

  theta_start <- control$axial_truncnorm_mixture2_start_theta %||% axial_truncnorm_default_theta_start(
    data = z,
    weights = obs_weights
  )
  theta_start <- normalize_axial_truncnorm_mixture2_theta(theta_start)
  start_par <- axial_truncnorm_encode_theta(theta_start)

  objective <- function(par) {
    tryCatch({
      theta <- axial_truncnorm_decode_theta(par)
      log_density <- axial_truncnorm_mixture_log_density(z, theta)
      if (any(!is.finite(log_density))) {
        return(Inf)
      }
      -sum(obs_weights * log_density)
    }, error = function(e) {
      Inf
    })
  }

  optim_control <- modifyList(list(maxit = 300L, reltol = 1e-8), control$axial_truncnorm_mixture2_optim_control %||% list())
  fit <- optim(
    par = start_par,
    fn = objective,
    method = "BFGS",
    control = optim_control
  )

  theta_hat <- axial_truncnorm_decode_theta(fit$par)
  c(
    theta_hat,
    list(
      loglik = -fit$value,
      opt = fit
    )
  )
}

make_axial_truncnorm_mixture2_spec <- function(distance_type = c("euclidean")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("axial_truncnorm_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_axial_truncnorm_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      z <- normalize_axial_truncnorm_mixture2_data(data, control)
      omega_vec <- normalize_axial_truncnorm_mixture2_data(omega, control)
      abs(outer(z, omega_vec, FUN = "-"))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      axial_truncnorm_mixture_distance_profile(
        omega = omega,
        t_values = as.numeric(t),
        theta = theta
      )
    },
    normalize_data = normalize_axial_truncnorm_mixture2_data,
    n_obs = function(data, control = list()) {
      length(normalize_axial_truncnorm_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_axial_truncnorm_mixture2_data(data, control)[[idx]]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        axial_truncnorm_mixture_distance_profile_grid(
          omega_grid = omega_grid,
          t_grid = t_grid,
          theta = theta
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE
    )
  )
}

small_circle_model_spec_path <- resolve_bootstrap_path("bootstrap", "small_circle_model_spec.R")
if (!exists("make_small_circle_spec", mode = "function") && file.exists(small_circle_model_spec_path)) {
  source(small_circle_model_spec_path)
}

watson_model_spec_path <- resolve_bootstrap_path("bootstrap", "watson_model_spec.R")
if (!exists("make_watson_spec", mode = "function") && file.exists(watson_model_spec_path)) {
  source(watson_model_spec_path)
}

small_circle_symmetric_mixture2_model_spec_path <- resolve_bootstrap_path("bootstrap", "small_circle_symmetric_mixture2_model_spec.R")
if (!exists("make_small_circle_symmetric_mixture2_spec", mode = "function") &&
    file.exists(small_circle_symmetric_mixture2_model_spec_path)) {
  source(small_circle_symmetric_mixture2_model_spec_path)
}
