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

deterministic_derivatives_path <- resolve_bootstrap_path(
  "bootstrap",
  "deterministic_profile_derivatives.R"
)
safe_source_text(deterministic_derivatives_path, envir = environment())

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


spec_sample_profile_sorted_prepare <- function(spec,
                                               data,
                                               sorted_distance_matrix,
                                               theta,
                                               control = list()) {
  if (!is.function(spec$sample_profile_sorted_prepare)) {
    return(NULL)
  }

  spec$sample_profile_sorted_prepare(
    data = data,
    sorted_distance_matrix = sorted_distance_matrix,
    theta = theta,
    control = control
  )
}

spec_sample_profile_sorted_block_eval <- function(spec,
                                                  data,
                                                  sorted_distance_matrix,
                                                  theta,
                                                  row_indices,
                                                  prepared = NULL,
                                                  control = list()) {
  if (!is.function(spec$sample_profile_sorted_block_eval)) {
    return(NULL)
  }

  spec$sample_profile_sorted_block_eval(
    data = data,
    sorted_distance_matrix = sorted_distance_matrix,
    theta = theta,
    row_indices = row_indices,
    prepared = prepared,
    control = control
  )
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

fast_multiplier_parse_derivative_control <- function(control = list(),
                                                     default_method = "score_mc") {
  default_method <- tolower(as.character(default_method))
  if (length(default_method) != 1L ||
      !default_method %in% c("score_mc", "quadrature")) {
    stop("`default_method` must be either 'score_mc' or 'quadrature'.")
  }
  method_was_supplied <- !is.null(control$derivative_method)
  derivative_method_requested <- tolower(as.character(
    control$derivative_method %||% "auto"
  ))
  derivative_method <- derivative_method_requested
  if (identical(derivative_method, "deterministic")) {
    derivative_method <- "quadrature"
    derivative_method_requested <- "quadrature"
  }
  supported_methods <- c("auto", "score_mc", "quadrature")
  if (length(derivative_method) != 1L ||
      !derivative_method %in% supported_methods) {
    stop(sprintf(
      "`control$derivative_method` must be one of %s.",
      paste(sprintf("'%s'", supported_methods), collapse = ", ")
    ))
  }
  if (identical(derivative_method, "auto")) {
    derivative_method <- default_method
  }

  derivative_mc_size <- as.integer(control$derivative_mc_size %||% 1000L)
  if (!is.finite(derivative_mc_size) || derivative_mc_size <= 0L) {
    stop("`control$derivative_mc_size` must be a strictly positive integer.")
  }

  derivative_mc_seed <- control$derivative_mc_seed %||% control$seed %||% NULL

  list(
    derivative_method_requested = derivative_method_requested,
    derivative_method = derivative_method,
    derivative_method_effective = derivative_method,
    derivative_method_selection_source = if (method_was_supplied) {
      if (identical(derivative_method_requested, "auto")) "explicit_auto" else "explicit"
    } else {
      "model_default"
    },
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

fast_multiplier_duplication_matrix <- function(dim_size) {
  dim_size <- as.integer(dim_size)
  if (length(dim_size) != 1L || !is.finite(dim_size) || dim_size < 1L) {
    stop("`dim_size` must be a positive integer.")
  }
  template <- matrix(0, dim_size, dim_size)
  index <- which(lower.tri(template, diag = TRUE), arr.ind = TRUE)
  D <- matrix(0, nrow = dim_size^2, ncol = nrow(index))
  for (column in seq_len(nrow(index))) {
    i <- index[column, 1L]
    j <- index[column, 2L]
    D[i + (j - 1L) * dim_size, column] <- 1
    if (i != j) D[j + (i - 1L) * dim_size, column] <- 1
  }
  D
}

gaussian_score_matrix_vech <- function(x,
                                       mu,
                                       Sigma,
                                       unknown_param = "both") {
  x <- as.matrix(x)
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  q <- length(mu)
  if (ncol(x) != q || !identical(dim(Sigma), c(q, q)) ||
      any(!is.finite(c(x, mu, Sigma)))) {
    stop("Gaussian score inputs have incompatible dimensions or non-finite values.")
  }
  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported Gaussian unknown-parameter specification.")
  }
  Sigma_inv <- solve(0.5 * (Sigma + t(Sigma)))
  centered <- sweep(x, 2L, mu, FUN = "-")
  score_mu <- centered %*% Sigma_inv
  transformed <- score_mu
  index <- which(lower.tri(Sigma, diag = TRUE), arr.ind = TRUE)
  score_sigma <- vapply(seq_len(nrow(index)), function(k) {
    i <- index[k, 1L]
    j <- index[k, 2L]
    symmetry_factor <- if (i == j) 1 else 2
    0.5 * symmetry_factor * (
      transformed[, i] * transformed[, j] - Sigma_inv[i, j]
    )
  }, numeric(nrow(x)))
  if (is.null(dim(score_sigma))) {
    score_sigma <- matrix(score_sigma, nrow = nrow(x))
  }
  switch(
    unknown_param,
    mu = score_mu,
    sigma = score_sigma,
    both = cbind(score_mu, score_sigma)
  )
}

gaussian_fisher_information_vech <- function(Sigma,
                                              unknown_param = "both") {
  Sigma <- as.matrix(Sigma)
  q <- nrow(Sigma)
  if (!identical(dim(Sigma), c(q, q)) || any(!is.finite(Sigma))) {
    stop("`Sigma` must be a finite square matrix.")
  }
  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported Gaussian unknown-parameter specification.")
  }
  Sigma_inv <- solve(0.5 * (Sigma + t(Sigma)))
  D <- fast_multiplier_duplication_matrix(q)
  information_sigma <- 0.5 * crossprod(
    D,
    kronecker(Sigma_inv, Sigma_inv) %*% D
  )
  switch(
    unknown_param,
    mu = Sigma_inv,
    sigma = information_sigma,
    both = {
      output <- matrix(0, q + ncol(D), q + ncol(D))
      output[seq_len(q), seq_len(q)] <- Sigma_inv
      sigma_index <- q + seq_len(ncol(D))
      output[sigma_index, sigma_index] <- information_sigma
      output
    }
  )
}

gaussian_mle_influence_matrix_vech <- function(x,
                                                mu,
                                                Sigma,
                                                unknown_param = "both") {
  score <- gaussian_score_matrix_vech(x, mu, Sigma, unknown_param)
  information <- gaussian_fisher_information_vech(Sigma, unknown_param)
  score %*% solve(information)
}

# For a Gaussian likelihood parametrised by (mu, vech(Sigma)), this returns
# the covariance of the MLE influence function for one observation.  Hence
# the matrix V in the paper's notation is -solve(output), because V is the
# derivative of theta -> E{psi_theta(X)} and psi is the likelihood score.
fast_multiplier_gaussian_influence_covariance <- function(Sigma,
                                                           unknown_param = "both") {
  Sigma <- as.matrix(Sigma)
  d <- nrow(Sigma)
  if (!identical(dim(Sigma), c(d, d)) || any(!is.finite(Sigma))) {
    stop("`Sigma` must be a finite square matrix.")
  }
  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported Gaussian unknown-parameter specification.")
  }
  idx <- which(lower.tri(Sigma, diag = TRUE), arr.ind = TRUE)
  p_sigma <- nrow(idx)
  covariance_sigma <- matrix(0, p_sigma, p_sigma)
  for (a in seq_len(p_sigma)) {
    i <- idx[a, 1L]
    j <- idx[a, 2L]
    for (b in seq_len(p_sigma)) {
      k <- idx[b, 1L]
      ell <- idx[b, 2L]
      covariance_sigma[a, b] <-
        Sigma[i, k] * Sigma[j, ell] + Sigma[i, ell] * Sigma[j, k]
    }
  }
  switch(
    unknown_param,
    mu = Sigma,
    sigma = covariance_sigma,
    both = {
      out <- matrix(0, d + p_sigma, d + p_sigma)
      out[seq_len(d), seq_len(d)] <- Sigma
      out[d + seq_len(p_sigma), d + seq_len(p_sigma)] <- covariance_sigma
      out
    }
  )
}

fast_multiplier_gaussian_paper_vhat <- function(Sigma,
                                                 unknown_param = "both") {
  -solve(fast_multiplier_gaussian_influence_covariance(
    Sigma = Sigma,
    unknown_param = unknown_param
  ))
}

fast_multiplier_matrix_condition_diagnostics <- function(matrix_value) {
  matrix_value <- as.matrix(matrix_value)
  rcond_value <- tryCatch(rcond(matrix_value), error = function(e) NA_real_)
  list(
    eigenvalues = as.numeric(eigen((matrix_value + t(matrix_value)) / 2,
                                   symmetric = TRUE, only.values = TRUE)$values),
    rcond = as.numeric(rcond_value),
    condition_number = if (is.finite(rcond_value) && rcond_value > 0) {
      as.numeric(1 / rcond_value)
    } else {
      Inf
    }
  )
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
  if (!identical(derivative_control$derivative_method, "score_mc")) {
    stop(sprintf(
      "Model '%s' does not implement deterministic profile derivatives.",
      spec$name
    ))
  }
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

# Load every model adapter from this compatibility entry point.  Keeping the
# common interface here lets scripts that source `model_specs.R` continue to
# work, while each family owns its own implementation file.
# EN DUDA (2026-07-26): MVN quadratic-form infrastructure is deliberately
# loaded before both the MVN and logistic-Gaussian adapters.  Both adapters
# then evaluate exactly the same weighted noncentral chi-square CDF.
mvnormal_quadform_path <- resolve_bootstrap_path("bootstrap", "mvnormal_quadform.R")
source(mvnormal_quadform_path)

model_spec_adapter_files <- c(
  "normal_model_spec.R",
  "mvnormal_model_spec.R",
  "logistic_gaussian_model_spec.R",
  "vmf_model_spec.R",
  "jp_model_spec.R",
  "hvmf_model_spec.R",
  "spherical_cauchy_model_spec.R",
  "beta_mixture2_model_spec.R",
  "uniform_beta_mixture_model_spec.R",
  "logitnormal_mixture2_model_spec.R",
  "axial_truncnorm_mixture2_model_spec.R",
  "cardioid_model_spec.R",
  "small_circle_model_spec.R",
  "watson_model_spec.R",
  "small_circle_symmetric_mixture2_model_spec.R",
  "small_circle_weighted_mixture2_model_spec.R",
  "sunspots_joint_time_space_model_spec.R"
)

for (adapter_file in model_spec_adapter_files) {
  adapter_path <- resolve_bootstrap_path("bootstrap", adapter_file)
  source(adapter_path)
}
rm(adapter_file, adapter_path, mvnormal_quadform_path)

install_distance_profile_backend_wrappers(
  c(
    "evaluate_mvnorm_distance_profile",
    "evaluate_mvnorm_distance_profile_matrix",
    "axial_truncnorm_mixture_distance_profile",
    "axial_truncnorm_mixture_distance_profile_grid"
  ),
  envir = environment(),
  cpp_supported = FALSE
)
