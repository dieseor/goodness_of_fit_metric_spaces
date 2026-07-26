# Axial truncated-normal-mixture model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_axial_truncnorm_mixture2 <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_axial_truncnorm_mixture2 <- model_specs_candidates_axial_truncnorm_mixture2[file.exists(model_specs_candidates_axial_truncnorm_mixture2)][1L]
  if (is.na(model_specs_path_axial_truncnorm_mixture2)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_axial_truncnorm_mixture2)
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
