# Cardioid model adapter for multiplier bootstrap GOF tests

resolve_cardioid_model_spec_path <- function(...) {
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

model_specs_path_cardioid <- resolve_cardioid_model_spec_path("bootstrap", "model_specs.R")
if (!exists("new_model_spec", mode = "function")) {
  source(model_specs_path_cardioid)
}

cardioid_source_path_model_spec <- resolve_cardioid_model_spec_path(
  "real_data",
  "comets",
  "cardioid",
  "legacy_materials",
  "Comets",
  "unregalitonavideno",
  "cardioid-source.R"
)
if (!exists("d_sph_car", mode = "function") || !exists("p_proj_car_gamma", mode = "function")) {
  source(cardioid_source_path_model_spec)
}

clip_cardioid_dot_products <- function(dot_products) {
  pmax(pmin(dot_products, 1), -1)
}

normalize_cardioid_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) < 2L) {
    stop("Cardioid data must be a non-empty matrix with at least two columns.")
  }
  if (any(!is.finite(data))) {
    stop("Cardioid data must be finite.")
  }

  norms <- sqrt(rowSums(data^2))
  if (any(norms <= 0)) {
    stop("Cardioid data rows must have strictly positive norm.")
  }

  data / norms
}

normalize_cardioid_theta <- function(theta,
                                     k = NULL,
                                     ambient_dim = NULL) {
  if (!is.list(theta)) {
    stop("Cardioid theta must be a list containing `mu` and `rho`.")
  }

  rho <- as.numeric(theta$rho)
  mu <- as.numeric(theta$mu)
  theta_k <- theta$k %||% k

  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho > 1) {
    stop("Cardioid theta requires `rho` in [0, 1].")
  }
  if (is.null(theta_k)) {
    stop("Cardioid theta requires `k`.")
  }
  theta_k <- as.integer(theta_k)
  if (length(theta_k) != 1L || !is.finite(theta_k) || theta_k < 1L) {
    stop("Cardioid theta requires a positive integer `k`.")
  }
  if (length(mu) < 2L || any(!is.finite(mu))) {
    stop("Cardioid theta requires a finite vector `mu` of length at least 2.")
  }
  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 0) {
    stop("Cardioid theta requires `mu` with strictly positive norm.")
  }
  mu <- mu / mu_norm

  if (!is.null(ambient_dim) && length(mu) != ambient_dim) {
    stop("Cardioid theta has incompatible ambient dimension.")
  }

  list(
    mu = mu,
    rho = rho,
    k = theta_k,
    p = length(mu)
  )
}

weighted_cardioid_resultant <- function(X, weights = NULL) {
  if (is.null(weights)) {
    return(colMeans(X))
  }

  prob_weights <- normalize_probability_weights(weights, n_expected = nrow(X))
  colSums(X * prob_weights)
}

normalize_cardioid_mle_weights <- function(weights, n_expected) {
  if (is.null(weights)) {
    return(rep.int(1, n_expected))
  }

  weights <- as.numeric(weights)
  if (length(weights) != n_expected) {
    stop("`weights` has incompatible length.")
  }
  if (any(!is.finite(weights))) {
    stop("`weights` must be finite.")
  }
  if (any(weights < 0)) {
    stop("`weights` must be nonnegative.")
  }

  weight_mean <- mean(weights)
  if (weight_mean <= 0) {
    stop("`weights` must have strictly positive mean.")
  }

  weights / weight_mean
}

cardioid_distance_threshold <- function(t_values, distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  t_values <- as.numeric(t_values)

  if (identical(distance_type, "geodesic")) {
    return(cos(t_values))
  }

  1 - (t_values^2) / 2
}

theoretical_distance_profile_cardioid <- function(omega,
                                                  mu,
                                                  rho,
                                                  k,
                                                  t_values,
                                                  distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega <- normalize_cardioid_data(omega)
  theta <- normalize_cardioid_theta(
    list(mu = mu, rho = rho, k = k),
    ambient_dim = ncol(omega)
  )

  omega_vec <- as.numeric(omega[1L, , drop = TRUE])
  t_values <- as.numeric(t_values)
  output <- numeric(length(t_values))

  if (identical(distance_type, "geodesic")) {
    output[t_values <= 0] <- 0
    output[t_values >= pi] <- 1
    active <- which(t_values > 0 & t_values < pi)
  } else {
    output[t_values <= 0] <- 0
    output[t_values >= 2] <- 1
    active <- which(t_values > 0 & t_values < 2)
  }

  if (length(active) == 0L) {
    return(output)
  }

  thresholds <- cardioid_distance_threshold(t_values[active], distance_type = distance_type)
  thresholds <- clip_cardioid_dot_products(thresholds)

  output[active] <- 1 - p_proj_car_gamma(
    x = thresholds,
    rho = theta$rho,
    k = theta$k,
    p = theta$p,
    mu = theta$mu,
    gamma = omega_vec
  )

  output
}

mle_sph_car_weighted <- function(X,
                                 k,
                                 weights = NULL,
                                 mu0 = NULL,
                                 rho0 = NULL,
                                 theta_start = NULL,
                                 control = list()) {
  X <- normalize_cardioid_data(X, control)
  p <- ncol(X)
  n <- nrow(X)
  scaled_weights <- normalize_cardioid_mle_weights(weights, n_expected = n)

  if (!is.null(theta_start)) {
    theta_start <- normalize_cardioid_theta(theta_start, k = k, ambient_dim = p)
    if (is.null(mu0)) {
      mu0 <- theta_start$mu
    }
    if (is.null(rho0)) {
      rho0 <- theta_start$rho
    }
  }

  resultant <- weighted_cardioid_resultant(X, weights = weights)
  resultant_norm <- sqrt(sum(resultant^2))
  if (is.null(mu0)) {
    if (resultant_norm <= 1e-12) {
      mu0 <- rep(0, p)
      mu0[[1L]] <- 1
    } else {
      mu0 <- resultant / resultant_norm
    }
  } else {
    mu0 <- as.numeric(mu0)
    if (length(mu0) != p || any(!is.finite(mu0))) {
      stop("`mu0` has incompatible dimension or contains non-finite values.")
    }
    mu0_norm <- sqrt(sum(mu0^2))
    if (mu0_norm <= 0) {
      stop("`mu0` must have strictly positive norm.")
    }
    mu0 <- mu0 / mu0_norm
  }

  if (is.null(rho0)) {
    rho0 <- if (resultant_norm <= 1e-12) 0 else min(max(resultant_norm, 0), 1)
  } else {
    rho0 <- as.numeric(rho0)
    if (length(rho0) != 1L || !is.finite(rho0)) {
      stop("`rho0` must be a finite scalar.")
    }
    rho0 <- min(max(rho0, 0), 1)
  }

  objective <- function(par) {
    rho <- min(max(par[[1L]], 0), 1)
    mu_raw <- par[-1L]
    mu_norm <- sqrt(sum(mu_raw^2))
    if (!is.finite(mu_norm) || mu_norm <= 0) {
      return(.Machine$double.xmax / 100)
    }

    mu <- mu_raw / mu_norm
    logdens <- d_sph_car(x = X, mu = mu, rho = rho, k = k, log = TRUE)
    if (any(!is.finite(logdens))) {
      return(.Machine$double.xmax / 100)
    }

    -sum(scaled_weights * logdens)
  }

  optim_control <- control$cardioid_optim_control %||% list()
  opt <- stats::optim(
    par = c(rho0, mu0),
    fn = objective,
    method = "L-BFGS-B",
    lower = c(0, rep(-Inf, p)),
    upper = c(1, rep(Inf, p)),
    control = optim_control
  )
  mu_hat <- opt$par[-1L]
  mu_hat <- mu_hat / sqrt(sum(mu_hat^2))
  rho_hat <- min(max(opt$par[[1L]], 0), 1)

  list(
    mu = mu_hat,
    rho = rho_hat,
    k = as.integer(k),
    p = p,
    ll = -opt$value,
    opt = opt
  )
}

fit_cardioid_theta <- function(data,
                               k,
                               weights = NULL,
                               null,
                               unknown_param = "both",
                               control = list()) {
  X <- normalize_cardioid_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    theta_simple <- normalize_cardioid_theta(null$theta, k = k, ambient_dim = ncol(X))
    return(c(theta_simple, list(weighted_mle = FALSE)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (!identical(unknown_param, "both")) {
    stop("The cardioid adapter currently supports only `unknown_param = 'both'`.")
  }
  if (!is.null(null$fixed) && length(null$fixed) > 0L) {
    stop("The cardioid adapter does not currently support fixed sub-parameters.")
  }

  theta_start <- control$theta_start %||% null$theta_start %||% NULL
  fit <- mle_sph_car_weighted(
    X = X,
    k = k,
    weights = weights,
    theta_start = theta_start,
    control = control
  )

  c(fit, list(weighted_mle = TRUE))
}

prepare_cardioid_fast_multiplier <- function(spec,
                                             data,
                                             theta_hat,
                                             ks_prep = NULL,
                                             cvm_prep = NULL,
                                             control = list(),
                                             k) {
  x <- normalize_cardioid_data(data, control)
  theta_hat <- normalize_cardioid_theta(theta_hat, k = k, ambient_dim = ncol(x))
  rho_boundary_eps <- as.numeric(control$cardioid_fast_boundary_eps %||% 1e-8)
  if (abs(theta_hat$rho) <= rho_boundary_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "cardioid_rho_zero_nonidentification",
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_
    ))
  }
  if (theta_hat$rho >= 1 - rho_boundary_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "cardioid_rho_one_boundary",
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_
    ))
  }

  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(0, 0, theta_hat$rho)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_cardioid_theta(
      list(
        mu = mapped$mu,
        rho = min(max(par[[3L]], 0), 1 - 1e-10),
        k = k
      ),
      k = k,
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_cardioid_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- clip_cardioid_dot_products(as.numeric(sample %*% theta_state$mu))
    qk <- fast_multiplier_cardioid_legendre(z, k = k)
    qk_prime <- fast_multiplier_cardioid_legendre_prime(z, k = k)
    denom <- 1 + theta_state$rho * qk
    coeff_mu <- theta_state$rho * qk_prime / denom
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))
    cbind(score_mu, qk / denom)
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_car(
      n = n_aux,
      mu = theta_state$mu,
      rho = theta_state$rho,
      k = theta_state$k
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

make_cardioid_spec <- function(k,
                               distance_type = c("geodesic", "chordal"),
                               unknown_param = "both") {
  distance_type <- match.arg(distance_type)
  k <- as.integer(k)
  if (length(k) != 1L || !is.finite(k) || k < 1L) {
    stop("`k` must be a positive integer.")
  }

  new_model_spec(
    name = sprintf("cardioid_k%d_%s", k, distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_cardioid_theta(
        data = data,
        k = k,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_cardioid_data(data, control)
      omega_matrix <- normalize_cardioid_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- clip_cardioid_dot_products(x %*% t(omega_matrix))
      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_cardioid_theta(theta, k = k)
      theoretical_distance_profile_cardioid(
        omega = omega,
        mu = theta$mu,
        rho = theta$rho,
        k = theta$k,
        t_values = as.numeric(t),
        distance_type = distance_type
      )
    },
    normalize_data = normalize_cardioid_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_cardioid_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_cardioid_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      distance_type = distance_type,
      unknown_param = unknown_param,
      k = k,
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_cardioid_fast_multiplier(
          spec = make_cardioid_spec(
            k = k,
            distance_type = distance_type,
            unknown_param = unknown_param
          ),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          k = k
        )
      }
    )
  )
}
