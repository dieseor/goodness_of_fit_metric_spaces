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
    extras = list(unknown_param = unknown_param)
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

  pmax(pmin(output, 1), 0)
}

logistic_gaussian_quadform_tail_probability <- function(q,
                                                        lambda,
                                                        h,
                                                        delta,
                                                        control = list()) {
  method <- tolower(as.character(control$logistic_gaussian_quadform_method %||% "farebrother"))
  if (!method %in% c("farebrother", "davies")) {
    stop("`control$logistic_gaussian_quadform_method` must be 'farebrother' or 'davies'.")
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
  } else {
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
  }

  pmin(pmax(1 - res$Qq, 0), 1)
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
      weighted_mle = TRUE
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

      if (ambient_dim == 3L) {
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

        NULL
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_vmf_theta(theta)
        ambient_dim <- length(theta$mu)
        profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
        n_u <- as.integer(control$vmf_profile_n_u %||% 4097L)

        if (ambient_dim == 3L && identical(profile_method, "tabulated")) {
          return(distance_profile_vmf_s2_cvm_grid(
            X = data,
            mu = theta$mu,
            kappa = theta$kappa,
            n_u = n_u
          ))
        }

        NULL
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
      unknown_param = "theta"
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
      unknown_param = unknown_param
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
      unknown_param = "theta"
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
      weighted_mle = TRUE
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

small_circle_model_spec_path <- resolve_bootstrap_path("bootstrap", "small_circle_model_spec.R")
if (!exists("make_small_circle_spec", mode = "function") && file.exists(small_circle_model_spec_path)) {
  source(small_circle_model_spec_path)
}

small_circle_symmetric_mixture2_model_spec_path <- resolve_bootstrap_path("bootstrap", "small_circle_symmetric_mixture2_model_spec.R")
if (!exists("make_small_circle_symmetric_mixture2_spec", mode = "function") &&
    file.exists(small_circle_symmetric_mixture2_model_spec_path)) {
  source(small_circle_symmetric_mixture2_model_spec_path)
}
