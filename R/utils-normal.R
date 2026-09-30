# Mathematical and numerical helpers: normal.

simulate_limit_gaussian <- function(cov_matrix, M = 10000, seed = NULL, tol = 1e-10) {
  cov_matrix <- as.matrix(cov_matrix)
  n_total <- nrow(cov_matrix)

  validate_covariance_matrix(
    cov_matrix,
    symmetry_tol = tol,
    psd_tol = tol,
    stop_on_failure = TRUE
  )

  cat("Generating", M, "multivariate normal samples from", n_total, "dimensional process...\n")
  if (!is.null(seed)) set.seed(seed)
  gaussian_samples <- tryCatch(
    mvtnorm::rmvnorm(M, mean = rep(0, n_total), sigma = cov_matrix),
    error = function(e) {
      stop(sprintf("simulate_limit_gaussian failed: %s", e$message))
    }
  )

  supremum_values <- apply(gaussian_samples, 1, function(row) max(abs(row)))
  return(supremum_values)
}
theoretical_distance_profile_normal <- function(omega, mu, sigma, t_values) {
  n <- max(length(omega), length(t_values))
  omega <- rep(omega, length.out = n)
  t_values <- rep(t_values, length.out = n)
  result <- numeric(n)
  valid_t <- t_values > 0
  if (any(valid_t)) {
    if (exists("distance_profile_backend_current", mode = "function") &&
        identical(distance_profile_backend_current(), "cpp") &&
        length(mu) == 1L && length(sigma) == 1L) {
      return(distance_profile_cpp_call(
        "cpp_dp_normal_profile",
        as.numeric(omega),
        as.numeric(t_values),
        as.numeric(mu),
        as.numeric(sigma)
      ))
    }
    upper_bound <- (omega[valid_t] + t_values[valid_t] - mu) / sigma
    lower_bound <- (omega[valid_t] - t_values[valid_t] - mu) / sigma
    result[valid_t] <- pnorm(upper_bound) - pnorm(lower_bound)
  }
  return(result)
}
dotF_mu_normal <- function(omega, mu, sigma, t) {
  -(1/sigma) * dnorm((omega - mu + t)/sigma) +
    (1/sigma) * dnorm((omega - mu - t)/sigma)
}
dotF_sigma_normal <- function(omega, mu, sigma, t) {
  -(omega - mu + t)/sigma^2 * dnorm((omega - mu + t)/sigma) +
    (omega - mu - t)/sigma^2 * dnorm((omega - mu - t)/sigma)
}
beta_mixture2_normalize_theta <- function(theta,
                                                     ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Beta-mixture theta must be a list.")
  }

  params <- beta_mixture2_canonicalize_theta(theta)
  if (params$ambient_dim != ambient_dim) {
    stop("Beta-mixture theta has incompatible ambient dimension.")
  }
  params
}
beta_mixture2_gauss_jacobi <- function(n, alpha, beta) {
  n <- as.integer(n)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= -1 ||
      length(beta) != 1L || !is.finite(beta) || beta <= -1) {
    stop("Gauss-Jacobi parameters must be finite and greater than -1.")
  }

  ab <- alpha + beta
  total_mass <- exp((ab + 1) * log(2) + lgamma(alpha + 1) +
    lgamma(beta + 1) - lgamma(ab + 2))

  diagonal <- numeric(n)
  diagonal[[1L]] <- (beta - alpha) / (ab + 2)
  if (n > 1L) {
    k_diag <- seq_len(n - 1L)
    diagonal[-1L] <- (beta^2 - alpha^2) /
      ((2 * k_diag + ab) * (2 * k_diag + ab + 2))

    k_off <- seq_len(n - 1L)
    offdiag <- sqrt(
      4 * k_off * (k_off + alpha) * (k_off + beta) * (k_off + ab) /
        ((2 * k_off + ab)^2 * ((2 * k_off + ab)^2 - 1))
    )
  } else {
    offdiag <- numeric(0L)
  }

  jacobi_matrix <- matrix(0, nrow = n, ncol = n)
  diag(jacobi_matrix) <- diagonal
  if (n > 1L) {
    jacobi_matrix[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- offdiag
    jacobi_matrix[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- offdiag
  }

  eig <- eigen(jacobi_matrix, symmetric = TRUE)
  order_idx <- order(eig$values)
  nodes <- eig$values[order_idx]
  weights <- total_mass * eig$vectors[1L, order_idx]^2

  list(nodes = nodes, weights = weights, total_mass = total_mass)
}
uniform_beta_mixture_normalize_theta <- function(theta,
                                                 ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Uniform-beta-mixture theta must be a list.")
  }

  params <- uniform_beta_mixture_validate_parameters(
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
  if (params$ambient_dim != ambient_dim) {
    stop("Uniform-beta-mixture theta has incompatible ambient dimension.")
  }
  params
}
logitnormal_mixture2_validate_parameters <- function(mu,
                                                                weight1,
                                                                mean1,
                                                                sd1,
                                                                mean2,
                                                                sd2) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Rotational logit-normal-mixture utilities currently support only S^2.")
  }

  weight1 <- as.numeric(weight1)
  mean1 <- as.numeric(mean1)
  mean2 <- as.numeric(mean2)
  sd1 <- as.numeric(sd1)
  sd2 <- as.numeric(sd2)

  if (length(weight1) != 1L || !is.finite(weight1) || weight1 <= 0 || weight1 >= 1) {
    stop("`weight1` must be a finite scalar in (0, 1).")
  }
  if (length(mean1) != 1L || !is.finite(mean1) || length(mean2) != 1L || !is.finite(mean2)) {
    stop("Logit-normal mixture means must be finite scalars.")
  }
  if (length(sd1) != 1L || !is.finite(sd1) || sd1 <= 0 || length(sd2) != 1L || !is.finite(sd2) || sd2 <= 0) {
    stop("Logit-normal mixture standard deviations must be strictly positive finite scalars.")
  }

  list(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2,
    ambient_dim = 3L
  )
}
logitnormal_mixture2_canonicalize_theta <- function(theta) {
  params <- logitnormal_mixture2_validate_parameters(
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2
  )

  if (params$mean1 <= params$mean2) {
    return(params)
  }

  list(
    mu = params$mu,
    weight1 = 1 - params$weight1,
    mean1 = params$mean2,
    sd1 = params$sd2,
    mean2 = params$mean1,
    sd2 = params$sd1,
    ambient_dim = params$ambient_dim
  )
}
logitnormal_mixture2_normalize_theta <- function(theta,
                                                            ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Logit-normal-mixture theta must be a list.")
  }
  params <- logitnormal_mixture2_canonicalize_theta(theta)
  if (params$ambient_dim != ambient_dim) {
    stop("Logit-normal-mixture theta has incompatible ambient dimension.")
  }
  params
}
logitnormal_mixture2_density_y <- function(y,
                                                      weight1,
                                                      mean1,
                                                      sd1,
                                                      mean2,
                                                      sd2,
                                                      log = FALSE,
                                                      eps = 1e-12) {
  y <- as.numeric(y)
  out <- rep(if (log) -Inf else 0, length(y))
  active <- y > 0 & y < 1
  if (!any(active)) {
    return(out)
  }

  y_active <- rotational_clamp_unit_interval(y[active], eps = eps)
  x_active <- stats::qlogis(y_active)
  log_jacobian <- -log(y_active) - log1p(-y_active)
  log_density <- rotational_logsumexp2(
    log(weight1) + stats::dnorm(x_active, mean = mean1, sd = sd1, log = TRUE) + log_jacobian,
    log1p(-weight1) + stats::dnorm(x_active, mean = mean2, sd = sd2, log = TRUE) + log_jacobian
  )
  out[active] <- if (log) log_density else exp(log_density)
  out
}
logitnormal_mixture2_cdf_y <- function(y,
                                                  weight1,
                                                  mean1,
                                                  sd1,
                                                  mean2,
                                                  sd2,
                                                  eps = 1e-12) {
  y <- as.numeric(y)
  out <- numeric(length(y))
  out[y <= 0] <- 0
  out[y >= 1] <- 1
  active <- which(y > 0 & y < 1)
  if (length(active) == 0L) {
    return(out)
  }

  y_active <- rotational_clamp_unit_interval(y[active], eps = eps)
  x_active <- stats::qlogis(y_active)
  out[active] <- weight1 * stats::pnorm((x_active - mean1) / sd1) +
    (1 - weight1) * stats::pnorm((x_active - mean2) / sd2)
  pmin(pmax(out, 0), 1)
}
logitnormal_mixture2_density_h <- function(z,
                                                      weight1,
                                                      mean1,
                                                      sd1,
                                                      mean2,
                                                      sd2,
                                                      eps = 1e-12) {
  z <- as.numeric(z)
  y <- (z + 1) / 2
  out <- numeric(length(z))
  valid <- z > -1 & z < 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- logitnormal_mixture2_density_y(
    y = y[valid],
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2,
    log = FALSE,
    eps = eps
  )
  out
}
logitnormal_mixture2_density_gz <- function(z,
                                                       weight1,
                                                       mean1,
                                                       sd1,
                                                       mean2,
                                                       sd2,
                                                       eps = 1e-12) {
  0.5 * logitnormal_mixture2_density_h(
    z = z,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2,
    eps = eps
  )
}
d_sph_logitnormal_mixture2_s2 <- function(x,
                                                     mu,
                                                     weight1,
                                                     mean1,
                                                     sd1,
                                                     mean2,
                                                     sd2,
                                                     log = FALSE,
                                                     eps = 1e-12) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )
  z <- pmin(pmax(as.numeric(x %*% params$mu), -1 + eps, 1 - eps), 1 - eps)
  y <- (z + 1) / 2
  log_density <- -log(4 * pi) + logitnormal_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    mean1 = params$mean1,
    sd1 = params$sd1,
    mean2 = params$mean2,
    sd2 = params$sd2,
    log = TRUE,
    eps = eps
  )
  if (log) log_density else exp(log_density)
}
logitnormal_mixture2_weighted_loglik_s2 <- function(mu,
                                                               weight1,
                                                               mean1,
                                                               sd1,
                                                               mean2,
                                                               sd2,
                                                               x,
                                                               prob_weights = NULL,
                                                               eps = 1e-12) {
  params <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  z <- pmin(pmax(as.numeric(x %*% params$mu), -1 + eps, 1 - eps), 1 - eps)
  y <- (z + 1) / 2
  sum(prob_weights * logitnormal_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    mean1 = params$mean1,
    sd1 = params$sd1,
    mean2 = params$mean2,
    sd2 = params$sd2,
    log = TRUE,
    eps = eps
  ))
}
logitnormal_mixture2_start_thetas_s2 <- function(x,
                                                            weights = NULL,
                                                            control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  warm_start <- control$logitnormal_mixture2_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- logitnormal_mixture2_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$logitnormal_mixture2_warm_start_only %||% FALSE)) {
      return(list(warm_start))
    }
  }

  resultant <- colSums(x * prob_weights)
  second_moment <- rotational_weighted_covariance_matrix(x, prob_weights)
  eig <- eigen(second_moment, symmetric = TRUE)
  mu_candidates <- rotational_unique_mu_candidates(list(
    resultant,
    -resultant,
    eig$vectors[, which.max(eig$values)],
    -eig$vectors[, which.max(eig$values)],
    eig$vectors[, which.min(eig$values)],
    -eig$vectors[, which.min(eig$values)],
    c(0, 0, 1),
    c(0, 0, -1)
  ))

  split_probs <- c(0.35, 0.5, 0.65)
  for (mu0 in mu_candidates) {
    y <- rotational_clamp_unit_interval((pmin(pmax(as.numeric(x %*% mu0), -1 + 1e-10), 1 - 1e-10) + 1) / 2, eps = 1e-8)
    x_logit <- stats::qlogis(y)
    for (split_prob in split_probs) {
      threshold <- rotational_weighted_quantile(x_logit, prob_weights, prob = split_prob)
      group1 <- x_logit <= threshold
      if (all(group1) || !any(group1)) {
        next
      }

      w1 <- prob_weights[group1]
      w2 <- prob_weights[!group1]
      p1 <- sum(w1)
      if (p1 <= 1e-6 || p1 >= 1 - 1e-6) {
        next
      }

      mean1 <- rotational_weighted_mean(x_logit[group1], w1 / sum(w1))
      mean2 <- rotational_weighted_mean(x_logit[!group1], w2 / sum(w2))
      sd1 <- sqrt(max(rotational_weighted_variance(x_logit[group1], w1 / sum(w1), center = mean1), 0.05^2))
      sd2 <- sqrt(max(rotational_weighted_variance(x_logit[!group1], w2 / sum(w2), center = mean2), 0.05^2))
      out[[length(out) + 1L]] <- logitnormal_mixture2_canonicalize_theta(list(
        mu = mu0,
        weight1 = p1,
        mean1 = mean1,
        sd1 = sd1,
        mean2 = mean2,
        sd2 = sd2
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- logitnormal_mixture2_canonicalize_theta(list(
      mu = c(0, 0, 1),
      weight1 = 0.5,
      mean1 = -1,
      sd1 = 0.7,
      mean2 = 1,
      sd2 = 0.7
    ))
  }

  out
}
logitnormal_mixture2_pack_par <- function(theta,
                                                     control = list()) {
  theta <- logitnormal_mixture2_normalize_theta(theta, ambient_dim = 3L)
  list(
    par = c(
      stats::qlogis(theta$weight1),
      theta$mean1,
      log(theta$sd1),
      theta$mean2,
      log(theta$sd2),
      theta$mu
    ),
    theta = theta
  )
}
logitnormal_mixture2_unpack_par <- function(par,
                                                       control = list()) {
  clip_means <- isTRUE(control$logitnormal_mixture2_clip_means %||% FALSE)
  mean_lower <- as.numeric(control$logitnormal_mixture2_mean_lower %||% -8)
  mean_upper <- as.numeric(control$logitnormal_mixture2_mean_upper %||% 8)
  sd_lower <- as.numeric(control$logitnormal_mixture2_sd_lower %||% 0.05)
  sd_upper <- as.numeric(control$logitnormal_mixture2_sd_upper %||% 5)
  weight_eps <- as.numeric(control$logitnormal_mixture2_weight_eps %||% 0.01)

  mu_raw <- par[6:8]
  mu_hat <- rotational_unit_vector_fallback(mu_raw)
  mean1_hat <- par[[2L]]
  mean2_hat <- par[[4L]]
  if (clip_means) {
    mean1_hat <- pmin(pmax(mean1_hat, mean_lower), mean_upper)
    mean2_hat <- pmin(pmax(mean2_hat, mean_lower), mean_upper)
  }
  logitnormal_mixture2_canonicalize_theta(list(
    mu = mu_hat,
    weight1 = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    mean1 = mean1_hat,
    sd1 = rotational_positive_parameter(par[[3L]], lower = sd_lower, upper = sd_upper),
    mean2 = mean2_hat,
    sd2 = rotational_positive_parameter(par[[5L]], lower = sd_lower, upper = sd_upper)
  ))
}
logitnormal_mixture2_mle_s2_weighted <- function(x,
                                                            weights = NULL,
                                                            control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- logitnormal_mixture2_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(length(candidate_thetas), as.integer(control$logitnormal_mixture2_n_starts %||% 12L)))]

  objective <- function(par) {
    theta <- logitnormal_mixture2_unpack_par(par, control = control)
    value <- -logitnormal_mixture2_weighted_loglik_s2(
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      x = x,
      prob_weights = prob_weights,
      eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12)
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$logitnormal_mixture2_optim_method %||% "BFGS"
  optim_control <- control$logitnormal_mixture2_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- logitnormal_mixture2_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- logitnormal_mixture2_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$logitnormal_mixture2_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$logitnormal_mixture2_warm_start_only <- FALSE
    return(logitnormal_mixture2_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Logit-normal-mixture weighted MLE failed for all starting values.")
  }

  c(
    best$theta,
    list(
      loglik = -best$opt$value,
      opt = best$opt,
      weighted_mle = TRUE,
      start_theta = best$start_theta
    )
  )
}
logitnormal_mixture2_legendre_coefficients <- function(theta,
                                                                  l_max = 150L,
                                                                  quad_n = 1000L,
                                                                  tol = 1e-10,
                                                                  eps = 1e-12) {
  theta <- logitnormal_mixture2_normalize_theta(theta, ambient_dim = 3L)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  quad <- rotational_gauss_hermite(quad_n)
  x1 <- theta$mean1 + theta$sd1 * sqrt(2) * quad$nodes
  x2 <- theta$mean2 + theta$sd2 * sqrt(2) * quad$nodes
  z1 <- tanh(x1 / 2)
  z2 <- tanh(x2 / 2)

  legendre_1 <- rotational_legendre_matrix(z1, l_max = l_max)
  legendre_2 <- rotational_legendre_matrix(z2, l_max = l_max)
  expectation_1 <- as.numeric(crossprod(legendre_1, quad$weights / sqrt(pi)))
  expectation_2 <- as.numeric(crossprod(legendre_2, quad$weights / sqrt(pi)))
  expectations <- theta$weight1 * expectation_1 + (1 - theta$weight1) * expectation_2

  ell <- 0:l_max
  coeffs <- (2 * ell + 1) * expectations
  coeffs[[1L]] <- 1
  a0_error <- abs(expectations[[1L]] - 1)
  if (a0_error > tol) {
    stop(sprintf("Rotational logit-normal Legendre coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error
  )
}
distance_profile_logitnormal_mixture2 <- function(omega,
                                                             t_values,
                                                             mu,
                                                             weight1,
                                                             mean1,
                                                             sd1,
                                                             mean2,
                                                             sd2,
                                                             distance_type = c("geodesic", "chordal"),
                                                             method = c("legendre", "integral"),
                                                             l_max = 150L,
                                                             quad_n = 1000L,
                                                             tol = 1e-10,
                                                             eps = 1e-12,
                                                             validate_against_integral = FALSE,
                                                             validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }
    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_logitnormal_mixture2(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        eps = eps,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  r_value <- sum(omega * theta$mu)
  if (abs(r_value - 1) <= 1e-12) {
    out[active] <- 1 - logitnormal_mixture2_cdf_y(
      y = (thresholds + 1) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- logitnormal_mixture2_cdf_y(
      y = (1 - thresholds) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }

  if (identical(method, "integral")) {
    return(rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        logitnormal_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  coeffs <- logitnormal_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    eps = eps
  )$coefficients
  out[active] <- 1 - rotational_projection_cdf_legendre(
    x = thresholds,
    r = r_value,
    coefficients = coeffs
  )
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral)) {
    out_integral <- rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        logitnormal_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Logit-normal-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}
distance_profile_logitnormal_mixture2_grid <- function(omega_grid,
                                                                  mu,
                                                                  weight1,
                                                                  mean1,
                                                                  sd1,
                                                                  mean2,
                                                                  sd2,
                                                                  t_grid,
                                                                  distance_type = c("geodesic", "chordal"),
                                                                  method = c("legendre", "integral"),
                                                                  l_max = 150L,
                                                                  quad_n = 1000L,
                                                                  tol = 1e-10,
                                                                  eps = 1e-12) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  if (identical(method, "integral")) {
    omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * logitnormal_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
    r_values <- as.numeric(omega_grid %*% theta$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - logitnormal_mixture2_cdf_y(
        y = (thresholds + 1) / 2,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        eps = eps
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- logitnormal_mixture2_cdf_y(
        y = (1 - thresholds) / 2,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        eps = eps
      )
    }
    return(out)
  }

  coeffs <- logitnormal_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    eps = eps
  )$coefficients
  rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = theta$mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}
distance_profile_logitnormal_mixture2_cvm_grid <- function(X,
                                                                      mu,
                                                                      weight1,
                                                                      mean1,
                                                                      sd1,
                                                                      mean2,
                                                                      sd2,
                                                                      method = c("legendre", "integral"),
                                                                      l_max = 150L,
                                                                      quad_n = 1000L,
                                                                      tol = 1e-10,
                                                                      eps = 1e-12) {
  method <- match.arg(method)
  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)

  if (identical(method, "integral")) {
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * logitnormal_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- logitnormal_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    eps = eps
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  out <- 1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )

  pos_idx <- which(abs(r_values - 1) <= 1e-12)
  if (length(pos_idx) > 0L) {
    out[pos_idx, ] <- 1 - logitnormal_mixture2_cdf_y(
      y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
  }

  neg_idx <- which(abs(r_values + 1) <= 1e-12)
  if (length(neg_idx) > 0L) {
    out[neg_idx, ] <- logitnormal_mixture2_cdf_y(
      y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
  }

  out <- pmin(pmax(out, 0), 1)
  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = acos(dot_products[i, ]),
      values = out[i, ],
      upper_bound = pi
    )
  }
  out
}
r_sph_logitnormal_mixture2 <- function(n,
                                                  mu,
                                                  weight1,
                                                  mean1,
                                                  sd1,
                                                  mean2,
                                                  sd2,
                                                  check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  component1 <- stats::runif(n) <= theta$weight1
  t_values <- numeric(n)
  n1 <- sum(component1)
  n2 <- n - n1
  if (n1 > 0L) {
    t_values[component1] <- stats::rnorm(n1, mean = theta$mean1, sd = theta$sd1)
  }
  if (n2 > 0L) {
    t_values[!component1] <- stats::rnorm(n2, mean = theta$mean2, sd = theta$sd2)
  }

  y <- stats::plogis(t_values)
  z <- pmin(pmax(2 * y - 1, -1), 1)
  phi <- stats::runif(n, min = 0, max = 2 * pi)
  basis <- jp_orthonormal_complement(theta$mu)
  tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
  radial <- sqrt(pmax(0, 1 - z^2))
  x <- tcrossprod(z, theta$mu) + sweep(tangent, 1L, radial, FUN = "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Logit-normal-mixture sampler returned non-unit vectors.")
    }
  }
  x
}
