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

utils_path_model_specs <- resolve_bootstrap_path("utils.R")
if (!exists("theoretical_distance_profile_normal", mode = "function")) {
  source(utils_path_model_specs)
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
