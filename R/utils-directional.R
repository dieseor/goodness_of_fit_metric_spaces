# Mathematical and numerical helpers: directional.

sphere_distance <- function(omega1, omega2, distance_type = "chordal") {
  omega1 <- omega1 / sqrt(sum(omega1^2))
  omega2 <- omega2 / sqrt(sum(omega2^2))
  dot_product <- sum(omega1 * omega2)
  dot_product <- pmax(pmin(dot_product, 1), -1)
  if (distance_type == "chordal") {
    return(sqrt(2 * (1 - dot_product)))
  } else if (distance_type == "geodesic") {
    return(acos(dot_product))
  } else {
    stop("distance_type must be 'chordal' or 'geodesic'")
  }
}
normalize_hvmf_h2_data <- function(data, tol = 1e-10) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) != 3L) {
    stop("`data` must be a non-empty n x 3 matrix with rows in H^2.")
  }
  if (any(!is.finite(data))) {
    stop("HvMF data must be finite.")
  }
  if (any(data[, 1L] <= 0)) {
    stop("HvMF data rows must satisfy x[1] > 0.")
  }

  minkowski_norms <- -data[, 1L]^2 + rowSums(data[, -1L, drop = FALSE]^2)
  if (any(abs(minkowski_norms + 1) > tol)) {
    stop("HvMF data rows must satisfy <x, x>_M = -1 up to tolerance.")
  }

  data
}
normalize_hvmf_hq_data <- function(data, q = NULL, tol = 1e-10) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) < 3L) {
    stop("`data` must be a non-empty n x (q + 1) matrix with q >= 2.")
  }

  inferred_q <- ncol(data) - 1L
  if (is.null(q)) {
    q <- inferred_q
  } else {
    q_numeric <- as.numeric(q)
    if (length(q_numeric) != 1L || !is.finite(q_numeric) ||
        q_numeric < 2 || q_numeric != floor(q_numeric)) {
      stop("`q` must be an integer greater than or equal to 2.")
    }
    q <- as.integer(q_numeric)
    if (inferred_q != q) {
      stop("`data` must have exactly q + 1 columns.")
    }
  }

  if (any(!is.finite(data))) {
    stop("HvMF data must be finite.")
  }
  if (any(data[, 1L] <= 0)) {
    stop("HvMF data rows must satisfy x[1] > 0.")
  }

  minkowski_norms <- -data[, 1L]^2 + rowSums(data[, -1L, drop = FALSE]^2)
  if (any(abs(minkowski_norms + 1) > tol)) {
    stop("HvMF data rows must satisfy <x, x>_M = -1 up to tolerance.")
  }

  data
}
hvmf_validate_dimension <- function(q) {
  q_numeric <- as.numeric(q)
  if (length(q_numeric) != 1L || !is.finite(q_numeric) ||
      q_numeric < 2 || q_numeric != floor(q_numeric)) {
    stop("`q` must be an integer greater than or equal to 2.")
  }
  as.integer(q_numeric)
}
hvmf_validate_kappa <- function(kappa) {
  kappa <- as.numeric(kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }
  kappa
}
hvmf_log_normalizing_constant <- function(q, kappa) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  order <- (q - 1) / 2
  scaled_bessel_k <- besselK(kappa, nu = order, expon.scaled = TRUE)
  if (!is.finite(scaled_bessel_k) || scaled_bessel_k <= 0) {
    stop("Could not evaluate the scaled Bessel K function for the HvMF normalising constant.")
  }

  kappa - log(2) - order * log(2 * pi / kappa) - log(scaled_bessel_k)
}
hvmf_radial_density <- function(u, q, kappa, chi, log = FALSE) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  chi <- as.numeric(chi)
  if (length(chi) != 1L || !is.finite(chi) || chi < 0) {
    stop("`chi` must be a finite nonnegative scalar.")
  }
  if (!requireNamespace("rotasym", quietly = TRUE)) {
    stop("Package 'rotasym' is required for the HvMF radial density.")
  }

  u <- as.numeric(u)
  log_density <- rep.int(-Inf, length(u))
  valid <- is.finite(u) & u > 0 & u < log(.Machine$double.xmax) - 2
  if (any(valid)) {
    u_valid <- u[valid]
    lambda <- kappa * sinh(chi) * sinh(u_valid)
    log_radial_jacobian <- (q - 1) * log(sinh(u_valid))
    log_angular_integral <- -rotasym::c_vMF(p = q, kappa = lambda, log = TRUE)
    base_exponent <- -kappa * cosh(chi) * cosh(u_valid)

    stable <- is.finite(lambda) & lambda > 1e-6
    if (any(stable)) {
      order <- q / 2 - 1
      scaled_bessel_i <- besselI(
        lambda[stable],
        nu = order,
        expon.scaled = TRUE
      )
      log_angular_integral[stable] <-
        (q / 2) * log(2 * pi) +
        log(scaled_bessel_i) -
        order * log(lambda[stable])
      base_exponent[stable] <- -kappa * cosh(u_valid[stable] - chi)
    }

    log_density[valid] <-
      hvmf_log_normalizing_constant(q = q, kappa = kappa) +
      log_radial_jacobian + base_exponent + log_angular_integral
  }

  if (isTRUE(log)) log_density else exp(log_density)
}
hvmf_radial_cdf <- function(u,
                            q,
                            kappa,
                            chi,
                            rel.tol = 1e-8,
                            abs.tol = 1e-10) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  chi <- as.numeric(chi)
  if (length(chi) != 1L || !is.finite(chi) || chi < 0) {
    stop("`chi` must be a finite nonnegative scalar.")
  }

  u <- as.numeric(u)
  output <- numeric(length(u))
  output[is.infinite(u) & u > 0] <- 1
  finite_positive <- is.finite(u) & u > 0
  if (any(finite_positive)) {
    output[finite_positive] <- vapply(u[finite_positive], function(upper_bound) {
      integral <- stats::integrate(
        hvmf_radial_density,
        lower = 0,
        upper = upper_bound,
        q = q,
        kappa = kappa,
        chi = chi,
        rel.tol = rel.tol,
        abs.tol = abs.tol,
        stop.on.error = FALSE
      )
      if (!identical(integral$message, "OK") || !is.finite(integral$value)) {
        stop(sprintf("Could not integrate the HvMF radial density up to u = %.8g.", upper_bound))
      }
      integral$value
    }, numeric(1))
  }
  pmin(pmax(output, 0), 1)
}
hvmf_build_radial_quantile_table <- function(q,
                                             kappa,
                                             chi,
                                             p_max = 0.999,
                                             probability_step = 0.01,
                                             upper = 3) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  chi <- as.numeric(chi)
  p_max <- as.numeric(p_max)
  probability_step <- as.numeric(probability_step)
  upper <- as.numeric(upper)

  if (length(chi) != 1L || !is.finite(chi) || chi < 0) {
    stop("`chi` must be a finite nonnegative scalar.")
  }
  if (length(p_max) != 1L || !is.finite(p_max) || p_max <= 0 || p_max >= 1) {
    stop("`p_max` must be a finite scalar in (0, 1).")
  }
  if (length(probability_step) != 1L || !is.finite(probability_step) ||
      probability_step <= 0 || probability_step > p_max) {
    stop("`probability_step` must be a finite scalar in (0, p_max].")
  }
  if (length(upper) != 1L || !is.finite(upper) || upper <= 0) {
    stop("`upper` must be a strictly positive finite scalar.")
  }
  if (!requireNamespace("gbutils", quietly = TRUE)) {
    stop("Package 'gbutils' is required for the polar HvMF sampler.")
  }

  cdf_u <- function(upper_bounds) {
    hvmf_radial_cdf(
      u = upper_bounds,
      q = q,
      kappa = kappa,
      chi = chi,
      rel.tol = 1e-8,
      abs.tol = 1e-10
    )
  }

  effective_upper <- upper
  upper_cdf <- cdf_u(effective_upper)
  while (upper_cdf < p_max) {
    effective_upper <- 2 * effective_upper
    if (!is.finite(effective_upper) || effective_upper > 96) {
      stop("Could not bracket the requested HvMF radial quantiles.")
    }
    upper_cdf <- cdf_u(effective_upper)
  }

  probabilities <- seq(0, p_max, by = probability_step)
  if (tail(probabilities, 1L) < p_max) {
    probabilities <- c(probabilities, p_max)
  }
  probabilities <- unique(pmin(probabilities, p_max))

  quantile_values <- vapply(probabilities, function(probability) {
    if (probability == 0) {
      return(0)
    }
    gbutils::cdf2quantile(
      p = probability,
      lower = 0,
      upper = effective_upper,
      cdf = cdf_u
    )
  }, numeric(1))

  output <- data.frame(
    base = probabilities,
    quantile_values = quantile_values
  )
  attr(output, "upper") <- effective_upper
  output
}
hvmf_mean_projection_density <- function(y, q, kappa, log = FALSE) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  y <- as.numeric(y)

  log_density <- rep.int(-Inf, length(y))
  valid <- is.finite(y) & y >= 1
  if (any(valid)) {
    y_valid <- y[valid]
    power <- (q - 2) / 2
    log_polynomial <- if (power == 0) {
      numeric(length(y_valid))
    } else {
      polynomial_log_base <- rep.int(-Inf, length(y_valid))
      above_boundary <- y_valid > 1
      polynomial_log_base[above_boundary] <-
        2 * log(y_valid[above_boundary]) +
        log1p(-1 / y_valid[above_boundary]^2)
      power * polynomial_log_base
    }
    log_sphere_area <- log(2) + (q / 2) * log(pi) - lgamma(q / 2)
    log_density[valid] <-
      hvmf_log_normalizing_constant(q = q, kappa = kappa) +
      log_sphere_area + log_polynomial - kappa * y_valid
  }

  if (isTRUE(log)) log_density else exp(log_density)
}
hvmf_polar_sample_core <- local({
  quantile_cache <- new.env(parent = emptyenv())

  function(n,
           mu,
           kappa,
           delta = 0,
           tangent = NULL,
           check = TRUE,
           p_max = 0.999,
           probability_step = 0.01,
           upper = 3) {
    n_numeric <- as.numeric(n)
    if (length(n_numeric) != 1L || !is.finite(n_numeric) ||
        n_numeric < 1 || n_numeric != floor(n_numeric)) {
      stop("`n` must be a strictly positive integer.")
    }
    n <- as.integer(n_numeric)

    mu <- as.numeric(mu)
    if (length(mu) < 3L || any(!is.finite(mu))) {
      stop("`mu` must be a finite numeric vector of length q + 1 with q >= 2.")
    }
    q <- length(mu) - 1L
    if (isTRUE(check)) {
      mu <- as.numeric(normalize_hvmf_hq_data(mu, q = q, tol = 1e-10)[1L, , drop = TRUE])
    }

    kappa <- hvmf_validate_kappa(kappa)
    delta <- as.numeric(delta)
    if (length(delta) != 1L || !is.finite(delta) || delta < 0 || delta > pi) {
      stop("`delta` must be a finite scalar in [0, pi].")
    }
    p_max <- as.numeric(p_max)
    probability_step <- as.numeric(probability_step)
    upper <- as.numeric(upper)
    if (length(p_max) != 1L || !is.finite(p_max) || p_max <= 0 || p_max >= 1) {
      stop("`p_max` must be a finite scalar in (0, 1).")
    }
    if (length(probability_step) != 1L || !is.finite(probability_step) ||
        probability_step <= 0 || probability_step > p_max) {
      stop("`probability_step` must be a finite scalar in (0, p_max].")
    }
    if (length(upper) != 1L || !is.finite(upper) || upper <= 0) {
      stop("`upper` must be a strictly positive finite scalar.")
    }
    if (!requireNamespace("rotasym", quietly = TRUE)) {
      stop("Package 'rotasym' is required for the polar HvMF sampler.")
    }

    spatial_norm <- sqrt(sum(mu[-1L]^2))
    chi <- asinh(spatial_norm)
    mean_direction <- if (spatial_norm > sqrt(.Machine$double.eps)) {
      mu[-1L] / spatial_norm
    } else {
      c(1, rep.int(0, q - 1L))
    }

    tangent_direction <- NULL
    if (delta > 0) {
      if (is.null(tangent)) {
        if (q != 2L) {
          stop("`tangent` is required for a positive angular displacement when q > 2.")
        }
        tangent_direction <- c(-mean_direction[[2L]], mean_direction[[1L]])
      } else {
        tangent <- as.numeric(tangent)
        if (length(tangent) != q || any(!is.finite(tangent))) {
          stop("`tangent` must be a finite numeric vector of length q.")
        }
        tangent_residual <- tangent - sum(tangent * mean_direction) * mean_direction
        tangent_norm <- sqrt(sum(tangent_residual^2))
        if (!is.finite(tangent_norm) || tangent_norm <= 1e-10) {
          stop("`tangent` must have a nonzero component orthogonal to the spatial mean direction.")
        }
        tangent_direction <- tangent_residual / tangent_norm
      }
    }

    cache_key <- paste(
      q,
      formatC(kappa, digits = 17, format = "fg"),
      formatC(chi, digits = 17, format = "fg"),
      formatC(p_max, digits = 17, format = "fg"),
      formatC(probability_step, digits = 17, format = "fg"),
      formatC(upper, digits = 17, format = "fg"),
      sep = "_"
    )
    if (exists(cache_key, envir = quantile_cache, inherits = FALSE)) {
      quantile_table <- get(cache_key, envir = quantile_cache, inherits = FALSE)
    } else {
      quantile_table <- hvmf_build_radial_quantile_table(
        q = q,
        kappa = kappa,
        chi = chi,
        p_max = p_max,
        probability_step = probability_step,
        upper = upper
      )
      assign(cache_key, quantile_table, envir = quantile_cache)
    }

    radial_probability <- stats::runif(n, min = 0, max = p_max)
    u <- stats::approx(
      x = quantile_table$base,
      y = quantile_table$quantile_values,
      xout = radial_probability,
      rule = 2
    )$y
    lambda <- kappa * sinh(chi) * sinh(u)
    signs <- if (delta == 0) {
      rep.int(0, n)
    } else {
      sample(c(-1, 1), size = n, replace = TRUE)
    }
    angular_part <- matrix(NA_real_, nrow = n, ncol = q)
    zero_lambda <- !is.finite(lambda) | lambda <= .Machine$double.eps

    if (any(zero_lambda)) {
      if (q == 2L) {
        uniform_angles <- stats::runif(sum(zero_lambda), min = -pi, max = pi)
        angular_part[zero_lambda, ] <- cbind(cos(uniform_angles), sin(uniform_angles))
      } else {
        angular_part[zero_lambda, ] <- rotasym::r_unif_sphere(
          n = sum(zero_lambda),
          p = q
        )
      }
    }

    for (index in which(!zero_lambda)) {
      conditional_mean <- if (delta == 0) {
        mean_direction
      } else {
        cos(delta) * mean_direction +
          signs[[index]] * sin(delta) * tangent_direction
      }
      angular_part[index, ] <- rotasym::r_vMF(
        n = 1L,
        mu = conditional_mean,
        kappa = lambda[[index]]
      )
    }

    x <- cbind(cosh(u), sinh(u) * angular_part)

    if (isTRUE(check)) {
      minkowski_norms <- -x[, 1L]^2 + rowSums(x[, -1L, drop = FALSE]^2)
      if (any(!is.finite(x)) || any(abs(minkowski_norms + 1) > 1e-8)) {
        stop("Sampler output does not satisfy <x, x>_M = -1 up to tolerance.")
      }
      if (any(x[, 1L] <= 0)) {
        stop("Sampler output must satisfy x[1] > 0.")
      }
    }

    x
  }
})
rhvmf_polar <- function(n,
                        mu,
                        kappa,
                        check = TRUE,
                        p_max = 0.999,
                        probability_step = 0.01,
                        upper = 3) {
  hvmf_polar_sample_core(
    n = n,
    mu = mu,
    kappa = kappa,
    delta = 0,
    check = check,
    p_max = p_max,
    probability_step = probability_step,
    upper = upper
  )
}
rhvmf_angular_mixture <- function(n,
                                  mu,
                                  kappa,
                                  delta,
                                  tangent = NULL,
                                  check = TRUE,
                                  p_max = 0.999,
                                  probability_step = 0.01,
                                  upper = 3) {
  hvmf_polar_sample_core(
    n = n,
    mu = mu,
    kappa = kappa,
    delta = delta,
    tangent = tangent,
    check = check,
    p_max = p_max,
    probability_step = probability_step,
    upper = upper
  )
}
rhvmf_h2_polar <- function(n, mu, kappa, check = TRUE) {
  rhvmf_polar(n = n, mu = mu, kappa = kappa, check = check)
}
rhvmf_h2_angular_mixture <- function(n, mu, kappa, delta, check = TRUE) {
  rhvmf_angular_mixture(
    n = n,
    mu = mu,
    kappa = kappa,
    delta = delta,
    tangent = NULL,
    check = check
  )
}
hvmf_mle_h2 <- function(data, weights = NULL, tol = 1e-10) {
  data <- normalize_hvmf_h2_data(data, tol = tol)
  n <- nrow(data)

  if (is.null(weights)) {
    weights <- rep.int(1, n)
  }

  weights <- as.numeric(weights)
  if (length(weights) != n) {
    stop("`weights` must have length nrow(data).")
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and nonnegative.")
  }

  W <- sum(weights)
  if (!is.finite(W) || W <= 0) {
    stop("sum(weights) must be strictly positive.")
  }

  S <- colSums(data * weights)
  S_inner <- minkowski_inner_product(S, S)
  R_sq <- -S_inner
  if (!is.finite(R_sq) || R_sq <= 0) {
    stop("The weighted resultant has non-positive or non-finite hyperbolic norm.")
  }

  R <- sqrt(R_sq)
  if (!is.finite(R)) {
    stop("The weighted resultant has non-finite hyperbolic length.")
  }
  if (R <= W + tol) {
    stop("Degenerate or near-degenerate HvMF MLE: R/W must be greater than 1.")
  }

  xi_hat <- S / R
  xi_inner <- minkowski_inner_product(xi_hat, xi_hat)
  if (!is.finite(xi_inner) || abs(xi_inner + 1) > 100 * tol) {
    stop("Internal HvMF MLE check failed: the estimated `xi` is not on H^2.")
  }
  if (xi_hat[[1]] <= 0) {
    stop("Internal HvMF MLE check failed: the estimated `xi` does not satisfy x[1] > 0.")
  }

  kappa_hat <- W / (R - W)
  if (!is.finite(kappa_hat) || kappa_hat <= 0) {
    stop("Internal HvMF MLE check failed: the estimated `kappa` is not positive and finite.")
  }

  sinh_chi_hat <- sqrt(sum(xi_hat[-1L]^2))
  chi_hat <- asinh(sinh_chi_hat)
  theta_hat <- atan2(xi_hat[[3L]], xi_hat[[2L]])
  theta_hat_deg <- (theta_hat * 180 / pi) %% 360

  list(
    xi = xi_hat,
    mu = xi_hat,
    kappa = kappa_hat,
    chi = chi_hat,
    sinh_chi = sinh_chi_hat,
    theta = theta_hat,
    theta_deg = theta_hat_deg,
    R = R,
    W = W,
    xi_inner = xi_inner,
    resultant = S
  )
}
hvmf_minkowski_inner_product <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 3L || length(x) != length(y) ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("`x` and `y` must be finite vectors of the same length q + 1, with q >= 2.")
  }
  -x[[1L]] * y[[1L]] + sum(x[-1L] * y[-1L])
}
hvmf_mean_resultant_ratio <- function(q, kappa) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  nu <- (q - 1) / 2
  numerator <- besselK(kappa, nu = nu + 1, expon.scaled = TRUE)
  denominator <- besselK(kappa, nu = nu, expon.scaled = TRUE)
  value <- numerator / denominator
  if (!is.finite(value) || value <= 1) {
    stop("Could not evaluate the HvMF mean-resultant ratio.")
  }
  value
}
hvmf_kappa_from_mean_resultant_ratio <- function(q, ratio) {
  q <- hvmf_validate_dimension(q)
  ratio <- as.numeric(ratio)
  if (length(ratio) != 1L || !is.finite(ratio) || ratio <= 1) {
    stop("`ratio` must be a strictly greater than one finite scalar.")
  }

  objective <- function(kappa) hvmf_mean_resultant_ratio(q, kappa) - ratio
  lower <- 1e-8
  upper <- max(1, 1 / (ratio - 1))
  while (objective(upper) > 0) {
    upper <- 2 * upper
    if (!is.finite(upper) || upper > 1e12) {
      stop("Could not bracket the HvMF concentration MLE.")
    }
  }
  stats::uniroot(objective, lower = lower, upper = upper, tol = 1e-10)$root
}
hvmf_mle_hq <- function(data, weights = NULL, tol = 1e-10) {
  x <- normalize_hvmf_hq_data(data, tol = tol)
  q <- ncol(x) - 1L
  if (q == 2L) return(hvmf_mle_h2(x, weights = weights, tol = tol))

  n <- nrow(x)
  if (is.null(weights)) weights <- rep.int(1, n)
  weights <- as.numeric(weights)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be a finite nonnegative vector of length nrow(data).")
  }
  W <- sum(weights)
  if (!is.finite(W) || W <= 0) stop("sum(weights) must be strictly positive.")

  resultant <- colSums(x * weights)
  resultant_norm_sq <- -hvmf_minkowski_inner_product(resultant, resultant)
  if (!is.finite(resultant_norm_sq) || resultant_norm_sq <= 0) {
    stop("The weighted resultant has non-positive or non-finite hyperbolic norm.")
  }
  R <- sqrt(resultant_norm_sq)
  ratio <- R / W
  if (!is.finite(ratio) || ratio <= 1 + tol) {
    stop("Degenerate or near-degenerate HvMF MLE: R/W must be greater than one.")
  }
  mu <- resultant / R
  if (mu[[1L]] <= 0 || abs(hvmf_minkowski_inner_product(mu, mu) + 1) > 100 * tol) {
    stop("Internal HvMF MLE check failed for the estimated mean direction.")
  }
  kappa <- hvmf_kappa_from_mean_resultant_ratio(q, ratio)
  list(
    xi = mu, mu = mu, kappa = kappa, q = q, R = R, W = W,
    resultant = resultant, resultant_ratio = ratio,
    xi_inner = hvmf_minkowski_inner_product(mu, mu)
  )
}
hvmf_distance_matrix_hq <- function(omega_grid, data, tol = 1e-12) {
  omega_matrix <- normalize_hvmf_hq_data(omega_grid, tol = tol)
  data_matrix <- normalize_hvmf_hq_data(data, q = ncol(omega_matrix) - 1L, tol = tol)
  minkowski_data <- data_matrix
  minkowski_data[, 1L] <- -minkowski_data[, 1L]
  acosh(pmax(-(omega_matrix %*% t(minkowski_data)), 1))
}
hvmf_distance_profile_tabulated <- function(t_values,
                                            q,
                                            kappa,
                                            chi,
                                            grid_size = 4097L) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  chi <- as.numeric(chi)
  t_values <- as.numeric(t_values)
  grid_size <- as.integer(grid_size)
  if (length(chi) != 1L || !is.finite(chi) || chi < 0) {
    stop("`chi` must be a finite nonnegative scalar.")
  }
  if (any(!is.finite(t_values)) || any(t_values < 0)) {
    stop("`t_values` must be finite and nonnegative.")
  }
  if (!is.finite(grid_size) || grid_size < 3L) {
    stop("`grid_size` must be an integer of at least three.")
  }
  upper <- max(t_values)
  if (upper == 0) return(rep.int(0, length(t_values)))

  grid <- seq(0, upper, length.out = grid_size)
  density <- hvmf_radial_density(grid, q = q, kappa = kappa, chi = chi)
  increments <- diff(grid) * (density[-1L] + density[-length(density)]) / 2
  cdf <- c(0, cumsum(increments))
  values <- stats::approx(grid, cdf, xout = t_values, rule = 2)$y
  pmin(pmax(values, 0), 1)
}
hvmf_distance_profile_hq <- function(omega, mu, kappa, t_values,
                                     grid_size = 4097L) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  q <- length(mu) - 1L
  normalize_hvmf_hq_data(omega, q = q)
  normalize_hvmf_hq_data(mu, q = q)
  chi <- acosh(pmax(-hvmf_minkowski_inner_product(mu, omega), 1))
  hvmf_distance_profile_tabulated(
    t_values = t_values, q = q, kappa = kappa, chi = chi,
    grid_size = grid_size
  )
}
hvmf_distance_profile_hq_integral <- function(omega, mu, kappa, t_values,
                                              rel.tol = 1e-8,
                                              abs.tol = 1e-10) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  q <- length(mu) - 1L
  normalize_hvmf_hq_data(omega, q = q)
  normalize_hvmf_hq_data(mu, q = q)
  chi <- acosh(pmax(-hvmf_minkowski_inner_product(mu, omega), 1))
  hvmf_radial_cdf(
    u = t_values, q = q, kappa = kappa, chi = chi,
    rel.tol = rel.tol, abs.tol = abs.tol
  )
}
hvmf_distance_matrix <- function(omega_grid, data, tol = 1e-12) {
  omega_matrix <- normalize_hvmf_h2_data(omega_grid, tol = tol)
  data_matrix <- normalize_hvmf_h2_data(data, tol = tol)

  minkowski_weighted_data <- data_matrix
  minkowski_weighted_data[, 1L] <- -minkowski_weighted_data[, 1L]

  inner_products <- omega_matrix %*% t(minkowski_weighted_data)
  acosh(pmax(-inner_products, 1))
}
theoretical_distance_profile_hvmf <- function(omega, mu, kappa, t_values) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  t_values <- as.numeric(t_values)

  if (length(omega) != 3L || length(mu) != 3L) {
    stop("`omega` and `mu` must both have length 3.")
  }

  mu_omega <- minkowski_inner_product(mu, omega)
  q <- length(mu) - 1L

  sapply(t_values, function(t) {
    density_R <- function(r) {
      nu <- (q - 1) / 2
      log_bessel_k <- log(besselK(kappa, nu = nu, expon.scaled = TRUE)) - kappa

      log_sinh_power_term <- (q - 1) * log(sinh(r))
      log_exp_term <- kappa * mu_omega * cosh(r)
      log_numerator <- nu * log(kappa) - nu * log(2 * pi) - log(2) - log_bessel_k +
        log_sinh_power_term + log_exp_term

      mu_omega_sq_minus_1 <- mu_omega^2 - 1
      if (mu_omega_sq_minus_1 < 0) {
        if (abs(mu_omega_sq_minus_1) < 1e-12) {
          mu_omega_sq_term <- 0
        } else {
          stop(paste(
            "Mathematical error: mu_omega^2 - 1 =",
            mu_omega_sq_minus_1,
            "unexpected negative value inside square root."
          ))
        }
      } else {
        mu_omega_sq_term <- mu_omega_sq_minus_1
      }

      kappa_vmf <- kappa * sinh(r) * sqrt(mu_omega_sq_term)
      log_denominator <- rotasym::c_vMF(p = q, kappa = kappa_vmf, log = TRUE)

      exp(log_numerator - log_denominator)
    }

    cdf_result <- integrate(
      density_R,
      lower = 0,
      upper = t,
      rel.tol = 1e-6,
      abs.tol = 1e-8,
      stop.on.error = FALSE
    )

    cdf_result$value
  })
}
hvmf_projection_density_h2 <- function(y, alpha, kappa) {
  y <- as.numeric(y)
  alpha <- as.numeric(alpha)
  kappa <- as.numeric(kappa)

  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 1) {
    stop("`alpha` must be a finite scalar with alpha >= 1.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }

  output <- numeric(length(y))
  valid <- is.finite(y) & (y >= 1)
  if (!any(valid)) {
    return(output)
  }

  beta_sq <- max(alpha^2 - 1, 0)
  beta <- sqrt(beta_sq)
  y_valid <- y[valid]
  radial_term <- sqrt(pmax(y_valid^2 - 1, 0))

  if (beta == 0) {
    log_density <- log(kappa) + kappa - kappa * y_valid
  } else {
    z <- kappa * beta * radial_term
    log_i0 <- log(besselI(z, nu = 0, expon.scaled = TRUE)) + z
    log_density <- log(kappa) + kappa - kappa * alpha * y_valid + log_i0
  }

  output[valid] <- exp(log_density)
  output
}
hvmf_tabulate_projection_cdf_h2 <- function(alpha,
                                            kappa,
                                            y_max,
                                            grid_size = 4097L,
                                            grid = NULL) {
  alpha <- as.numeric(alpha)
  kappa <- as.numeric(kappa)
  y_max <- as.numeric(y_max)

  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 1) {
    stop("`alpha` must be a finite scalar with alpha >= 1.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }
  if (length(y_max) != 1L || !is.finite(y_max) || y_max < 1) {
    stop("`y_max` must be a finite scalar with y_max >= 1.")
  }

  if (is.null(grid)) {
    grid_size <- as.integer(grid_size)
    if (!is.finite(grid_size) || grid_size < 2L) {
      stop("`grid_size` must be an integer >= 2.")
    }
    y_grid <- seq(1, y_max, length.out = grid_size)
  } else {
    y_grid <- sort(as.numeric(grid))
    if (length(y_grid) < 2L) {
      stop("`grid` must contain at least two points.")
    }
    if (any(!is.finite(y_grid))) {
      stop("`grid` must be finite.")
    }
    if (min(y_grid) < 1 || max(y_grid) > y_max) {
      stop("`grid` must be contained in [1, y_max].")
    }
    if (y_grid[[1L]] > 1) {
      y_grid <- c(1, y_grid)
    }
    if (y_grid[[length(y_grid)]] < y_max) {
      y_grid <- c(y_grid, y_max)
    }
    y_grid <- unique(y_grid)
  }

  density_values <- hvmf_projection_density_h2(
    y = y_grid,
    alpha = alpha,
    kappa = kappa
  )

  cdf_values <- numeric(length(y_grid))
  if (length(y_grid) >= 2L) {
    increments <- diff(y_grid) * (density_values[-1L] + density_values[-length(y_grid)]) / 2
    cdf_values[-1L] <- cumsum(increments)
  }

  cdf_values <- pmin(pmax(cdf_values, 0), 1)
  cdf_values[[1L]] <- 0

  data.frame(
    y = y_grid,
    cdf = cdf_values
  )
}
hvmf_eval_projection_cdf_tabulated <- function(y, cdf_table) {
  hvmf_eval_projection_cdf_tabulated_spline(y = y, cdf_table = cdf_table)
}
hvmf_prepare_projection_cdf_interpolator <- function(cdf_table) {
  if (!is.data.frame(cdf_table) || !all(c("y", "cdf") %in% names(cdf_table))) {
    stop("`cdf_table` must be a data frame with columns `y` and `cdf`.")
  }

  y_grid <- as.numeric(cdf_table$y)
  cdf_grid <- as.numeric(cdf_table$cdf)
  if (length(y_grid) == 0L) {
    stop("`cdf_table` cannot be empty.")
  }

  if (length(y_grid) == 1L) {
    return(list(
      type = "constant",
      x_min = y_grid[[1L]],
      x_max = y_grid[[1L]],
      y_min = cdf_grid[[1L]],
      y_max = cdf_grid[[1L]]
    ))
  }

  list(
    type = "spline",
    x_min = y_grid[[1L]],
    x_max = y_grid[[length(y_grid)]],
    y_min = min(cdf_grid),
    y_max = max(cdf_grid),
    fn = stats::splinefun(
      x = y_grid,
      y = cdf_grid,
      method = "monoH.FC"
    )
  )
}
hvmf_eval_projection_cdf_tabulated_linear <- function(y, cdf_table) {
  y <- as.numeric(y)

  if (!is.data.frame(cdf_table) || !all(c("y", "cdf") %in% names(cdf_table))) {
    stop("`cdf_table` must be a data frame with columns `y` and `cdf`.")
  }

  y_grid <- as.numeric(cdf_table$y)
  cdf_grid <- as.numeric(cdf_table$cdf)
  if (length(y_grid) == 0L) {
    stop("`cdf_table` cannot be empty.")
  }
  if (length(y_grid) == 1L) {
    output <- rep(cdf_grid[[1L]], length(y))
    output[y <= y_grid[[1L]]] <- 0
    return(pmin(pmax(output, 0), 1))
  }

  output <- approx(
    x = y_grid,
    y = cdf_grid,
    xout = y,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y

  output[y <= 1] <- 0
  pmin(pmax(output, 0), 1)
}
hvmf_eval_projection_cdf_tabulated_spline <- function(y,
                                                      cdf_table,
                                                      interpolator = NULL) {
  y <- as.numeric(y)
  interpolator <- interpolator %||% hvmf_prepare_projection_cdf_interpolator(cdf_table)

  if (identical(interpolator$type, "constant")) {
    output <- rep(interpolator$y_min, length(y))
    output[y <= interpolator$x_min] <- 0
    return(pmin(pmax(output, 0), 1))
  }

  y_eval <- pmin(pmax(y, interpolator$x_min), interpolator$x_max)
  output <- interpolator$fn(y_eval)
  output[y <= 1] <- 0
  pmin(pmax(output, 0), 1)
}
hvmf_cvm_profile_matrix_tabulated <- function(data,
                                              theta,
                                              grid_size = 4097L,
                                              distance_matrix = NULL,
                                              tol = 1e-10) {
  x <- normalize_hvmf_h2_data(data, tol = tol)
  if (!is.list(theta)) {
    stop("`theta` must be a list containing `mu` and `kappa`.")
  }
  mu_matrix <- normalize_hvmf_h2_data(theta$mu, tol = tol)
  mu <- as.numeric(mu_matrix[1L, , drop = TRUE])
  kappa <- as.numeric(theta$kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("HvMF theta requires a strictly positive finite scalar `kappa`.")
  }
  n <- nrow(x)

  if (is.null(distance_matrix)) {
    y_matrix <- pmax(cosh(t(hvmf_distance_matrix(x, x))), 1)
  } else {
    distance_matrix <- as.matrix(distance_matrix)
    if (!all(dim(distance_matrix) == c(n, n))) {
      stop("`distance_matrix` must be an n x n matrix compatible with `data`.")
    }
    y_matrix <- pmax(cosh(distance_matrix), 1)
  }

  alpha_values <- as.numeric(x %*% c(mu[[1L]], -mu[[2L]], -mu[[3L]]))
  alpha_values <- pmax(alpha_values, 1)

  output <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    row_y <- pmax(y_matrix[i, ], 1)
    cdf_table <- hvmf_tabulate_projection_cdf_h2(
      alpha = alpha_values[[i]],
      kappa = kappa,
      y_max = max(row_y),
      grid_size = grid_size
    )
    cdf_interpolator <- hvmf_prepare_projection_cdf_interpolator(cdf_table)
    output[i, ] <- hvmf_eval_projection_cdf_tabulated_spline(
      y = row_y,
      cdf_table = cdf_table,
      interpolator = cdf_interpolator
    )
  }

  output
}
hvmf_cvm_profile_matrix_tabulated_linear <- function(data,
                                                     theta,
                                                     grid_size = 4097L,
                                                     distance_matrix = NULL,
                                                     tol = 1e-10) {
  x <- normalize_hvmf_h2_data(data, tol = tol)
  if (!is.list(theta)) {
    stop("`theta` must be a list containing `mu` and `kappa`.")
  }
  mu_matrix <- normalize_hvmf_h2_data(theta$mu, tol = tol)
  mu <- as.numeric(mu_matrix[1L, , drop = TRUE])
  kappa <- as.numeric(theta$kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("HvMF theta requires a strictly positive finite scalar `kappa`.")
  }
  n <- nrow(x)

  if (is.null(distance_matrix)) {
    y_matrix <- pmax(cosh(t(hvmf_distance_matrix(x, x))), 1)
  } else {
    distance_matrix <- as.matrix(distance_matrix)
    if (!all(dim(distance_matrix) == c(n, n))) {
      stop("`distance_matrix` must be an n x n matrix compatible with `data`.")
    }
    y_matrix <- pmax(cosh(distance_matrix), 1)
  }

  alpha_values <- as.numeric(x %*% c(mu[[1L]], -mu[[2L]], -mu[[3L]]))
  alpha_values <- pmax(alpha_values, 1)

  output <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    row_y <- pmax(y_matrix[i, ], 1)
    cdf_table <- hvmf_tabulate_projection_cdf_h2(
      alpha = alpha_values[[i]],
      kappa = kappa,
      y_max = max(row_y),
      grid_size = grid_size
    )
    output[i, ] <- hvmf_eval_projection_cdf_tabulated_linear(row_y, cdf_table)
  }

  output
}
vmf_s1_angle_density <- function(phi, mu, kappa, log = FALSE) {
  mu <- as.numeric(mu)
  if (length(mu) != 2) {
    stop("`mu` must have length 2 for S^1 computations.")
  }
  mu <- mu / sqrt(sum(mu^2))
  mu_angle <- circle_angle_from_point(mu)
  log_i0 <- log(besselI(kappa, nu = 0, expon.scaled = TRUE)) + kappa
  log_density <- kappa * cos(phi - mu_angle) - log(2 * pi) - log_i0
  if (log) {
    return(log_density)
  }
  exp(log_density)
}
build_vmf_s1_cdf <- function(mu, kappa, n_grid = 16385) {
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid < 3) {
    stop("`n_grid` must be an integer >= 3.")
  }
  phi <- seq(0, 2 * pi, length.out = n_grid)
  density <- vmf_s1_angle_density(phi, mu, kappa)
  dx <- diff(phi)
  cdf <- numeric(n_grid)
  cdf[-1] <- cumsum((density[-n_grid] + density[-1]) * dx / 2)
  total_mass <- cdf[n_grid]
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Failed to build a deterministic angular CDF for vMF on S^1.")
  }
  cdf <- cdf / total_mass
  cdf[n_grid] <- 1
  list(
    mu = as.numeric(mu) / sqrt(sum(mu^2)),
    kappa = kappa,
    phi = phi,
    density = density,
    cdf = cdf
  )
}
evaluate_vmf_s1_cdf <- function(phi, cdf_object, tol = 1e-12) {
  phi <- as.numeric(phi)
  phi_wrapped <- phi %% (2 * pi)
  is_full_turn <- abs(phi_wrapped) < tol & phi > 0
  phi_wrapped[is_full_turn] <- 2 * pi
  approx(
    x = cdf_object$phi,
    y = cdf_object$cdf,
    xout = phi_wrapped,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}
vmf_s1_segments_probability <- function(segments, cdf_object) {
  if (length(segments) == 0 || nrow(segments) == 0) {
    return(0)
  }
  sum(vapply(seq_len(nrow(segments)), function(i) {
    evaluate_vmf_s1_cdf(segments[i, 2], cdf_object) -
      evaluate_vmf_s1_cdf(segments[i, 1], cdf_object)
  }, numeric(1)))
}
theoretical_distance_profile_vmf_s1_chordal <- function(omega,
                                                        mu,
                                                        kappa,
                                                        t_values,
                                                        cdf_object = NULL,
                                                        cdf_grid_size = 16385) {
  if (is.matrix(omega)) {
    n <- nrow(omega)
    if (length(t_values) == 1) t_values <- rep(t_values, n)
    stopifnot(length(t_values) == n)
    return(vapply(seq_len(n), function(i) {
      theoretical_distance_profile_vmf_s1_chordal(
        omega = omega[i, ],
        mu = mu,
        kappa = kappa,
        t_values = t_values[i],
        cdf_object = cdf_object,
        cdf_grid_size = cdf_grid_size
      )
    }, numeric(1)))
  }

  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }

  t_values <- as.numeric(t_values)
  vapply(t_values, function(t) {
    if (t <= 0) return(0)
    if (t >= 2) return(1)
    segments <- s1_event_segments_chordal(omega, t)
    vmf_s1_segments_probability(segments, cdf_object)
  }, numeric(1))
}
sphere_distance_to_dot_threshold <- function(t, distance_type = "chordal") {
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  input_is_matrix <- is.matrix(t)
  t_dims <- if (input_is_matrix) dim(t) else NULL
  t <- as.numeric(t)
  if (distance_type == "chordal") {
    threshold <- 1 - (t^2) / 2
  } else {
    threshold <- cos(t)
  }
  threshold <- pmin(pmax(threshold, -1), 1)
  if (isTRUE(input_is_matrix)) {
    matrix(threshold, nrow = t_dims[[1L]], ncol = t_dims[[2L]])
  } else {
    threshold
  }
}
vmf_s2_projected_density <- function(u, omega, mu, kappa, log = FALSE) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3 || length(mu) != 3) {
    stop("`omega` and `mu` must both have length 3 for S^2 projected densities.")
  }
  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))

  u <- as.numeric(u)
  out <- rep(if (log) -Inf else 0, length(u))
  valid <- u >= -1 & u <= 1
  if (!any(valid)) {
    return(out)
  }

  m1 <- sum(mu * omega)
  m1 <- pmin(pmax(m1, -1), 1)
  radial_sq <- pmax(0, 1 - u[valid]^2)
  tangential_sq <- max(0, 1 - m1^2)
  lambda <- kappa * sqrt(radial_sq * tangential_sq)

  log_density <- rotasym::c_vMF(p = 3, kappa = kappa, log = TRUE) +
    kappa * m1 * u[valid] -
    rotasym::c_vMF(p = 2, kappa = lambda, log = TRUE)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}
distance_profile_vmf_s2_integral <- function(omega,
                                             mu,
                                             kappa,
                                             t_values,
                                             distance_type = "chordal",
                                             rel.tol = 1e-8,
                                             abs.tol = 1e-10,
                                             subdivisions = 200L) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3 || length(mu) != 3) {
    stop("`omega` and `mu` must both have length 3 for S^2 distance profiles.")
  }

  if (is.matrix(omega)) {
    n <- nrow(omega)
    if (length(t_values) == 1) t_values <- rep(t_values, n)
    stopifnot(length(t_values) == n)
    return(vapply(seq_len(n), function(i) {
      distance_profile_vmf_s2_integral(
        omega = omega[i, ],
        mu = mu,
        kappa = kappa,
        t_values = t_values[i],
        distance_type = distance_type,
        rel.tol = rel.tol,
        abs.tol = abs.tol,
        subdivisions = subdivisions
      )
    }, numeric(1)))
  }

  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t_values <- as.numeric(t_values)

  vapply(t_values, function(t) {
    a <- sphere_distance_to_dot_threshold(t, distance_type)
    if (a <= -1) return(1)
    if (a >= 1) return(0)
    integrate(
      f = function(u) vmf_s2_projected_density(u, omega = omega, mu = mu, kappa = kappa),
      lower = a,
      upper = 1,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      subdivisions = subdivisions
    )$value
  }, numeric(1))
}
vmf_log_normalizing_constant_intrinsic <- function(q, kappa) {
  q <- as.integer(q)
  kappa <- as.numeric(kappa)

  if (length(q) != 1L || !is.finite(q) || q < 1L ||
      length(kappa) == 0L || any(!is.finite(kappa)) || any(kappa < 0)) {
    stop("Invalid vMF dimension or concentration.")
  }

  output <- numeric(length(kappa))
  zero <- kappa == 0

  if (any(zero)) {
    output[zero] <-
      lgamma((q + 1) / 2) -
      log(2) -
      ((q + 1) / 2) * log(pi)
  }

  positive <- !zero
  if (any(positive)) {
    nu <- (q - 1) / 2
    scaled_bessel <- besselI(
      kappa[positive],
      nu = nu,
      expon.scaled = TRUE
    )

    if (any(!is.finite(scaled_bessel)) || any(scaled_bessel <= 0)) {
      stop("Could not evaluate the scaled Bessel function for the vMF normalising constant.")
    }

    output[positive] <-
      nu * log(kappa[positive]) -
      ((q + 1) / 2) * log(2 * pi) -
      log(scaled_bessel) -
      kappa[positive]
  }

  output
}
vmf_projected_density_canonical_state <- function(xi, omega) {
  xi <- as.numeric(xi)
  omega <- as.numeric(omega)
  q <- length(xi) - 1L

  if (q < 2L || length(omega) != length(xi) ||
      any(!is.finite(c(xi, omega)))) {
    stop("The vMF projected density requires finite vectors on S^q with q >= 2.")
  }

  omega_norm <- sqrt(sum(omega^2))
  if (!is.finite(omega_norm) || omega_norm <= 0) {
    stop("`omega` must have strictly positive norm.")
  }
  omega <- omega / omega_norm

  kappa <- sqrt(sum(xi^2))
  a <- sum(xi * omega)
  b_sq <- profile_derivative_nonnegative_square(
    kappa^2 - a^2,
    scale = max(kappa^2, a^2),
    label = "vMF b^2"
  )

  list(
    a = a,
    b = sqrt(b_sq),
    omega = omega,
    kappa = kappa,
    q = q,
    log_normalizing_constant =
      vmf_log_normalizing_constant_intrinsic(q, kappa)
  )
}
vmf_projected_density_canonical <- function(s,
                                            xi = NULL,
                                            omega = NULL,
                                            state = NULL) {
  s <- as.numeric(s)
  if (any(!is.finite(s))) {
    stop("The vMF projected density requires finite projection values.")
  }

  if (is.null(state)) {
    if (is.null(xi) || is.null(omega)) {
      stop("Supply either `state` or both `xi` and `omega`.")
    }
    state <- vmf_projected_density_canonical_state(xi = xi, omega = omega)
  }

  q <- state$q
  one_minus_s2 <- pmax(0, 1 - s^2)
  u <- state$b * sqrt(one_minus_s2)

  log_density <- rep.int(-Inf, length(s))
  interior_or_q2 <- one_minus_s2 > 0 | q == 2L

  if (any(interior_or_q2)) {
    power_term <- if (q == 2L) {
      rep.int(0, sum(interior_or_q2))
    } else {
      ((q - 2) / 2) * log(one_minus_s2[interior_or_q2])
    }

    log_density[interior_or_q2] <-
      state$log_normalizing_constant -
      vmf_log_normalizing_constant_intrinsic(
        q = q - 1L,
        kappa = u[interior_or_q2]
      ) +
      state$a * s[interior_or_q2] +
      power_term
  }

  list(
    density = exp(log_density),
    a = state$a,
    b = state$b,
    u = u,
    one_minus_s2 = one_minus_s2,
    omega = state$omega,
    kappa = state$kappa,
    q = q
  )
}
theoretical_distance_profile_vmf <- function(omega, mu, kappa, t_values, distance_type = "chordal") {
  # If omega is a matrix, vectorize over rows
  if (is.matrix(omega)) {
    n <- nrow(omega)
    # If t_values is scalar, repeat
    if (length(t_values) == 1) t_values <- rep(t_values, n)
    stopifnot(length(t_values) == n)
    return(sapply(1:n, function(i) {
      theoretical_distance_profile_vmf(omega[i, ], mu, kappa, t_values[i], distance_type)
    }))
  }
  if (length(mu) == 2) {
    if (distance_type == "chordal") {
      return(theoretical_distance_profile_vmf_s1_chordal(
        omega, mu, kappa, t_values
      ))
    }

    chordal_t <- sqrt(pmax(
      0,
      2 * (1 - cos(as.numeric(t_values)))
    ))
    return(theoretical_distance_profile_vmf_s1_chordal(
      omega, mu, kappa, chordal_t
    ))
  }
  rho <- sum(mu * omega)
  p <- length(mu)
  log_normalizing_constant <- rotasym::c_vMF(
    p = p,
    kappa = kappa,
    log = TRUE
  )
  sapply(t_values, function(t) {
    threshold <- ifelse(distance_type == "chordal", 1 - (t^2)/2, cos(t))
    density_T <- function(s) {
      exp(
        log_normalizing_constant +
          kappa * rho * s +
          ((p - 3) / 2) * log(1 - s^2) -
          rotasym::c_vMF(
            p = p - 1L,
            kappa = kappa * sqrt((1 - s^2) * (1 - rho^2)),
            log = TRUE
          )
      )
    }
    cdf_at_threshold <- integrate(density_T, lower = -1 + 1e-8, upper = threshold,
                                 rel.tol = 1e-8, abs.tol = 1e-10)$value
    return(1 - cdf_at_threshold)
  })
}
theoretical_distance_profile_vmf_s2_fast <- function(omega,
                                                     mu,
                                                     kappa,
                                                     t_values,
                                                     distance_type = "chordal",
                                                     rel.tol = 1e-8,
                                                     abs.tol = 1e-10) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3L || length(mu) != 3L) {
    stop("`omega` and `mu` must both have length 3 for S^2 vMF distance profiles.")
  }

  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t_values <- as.numeric(t_values)

  rho <- sum(mu * omega)
  rho <- pmin(pmax(rho, -1), 1)
  tangential_norm <- sqrt(pmax(0, 1 - rho^2))
  log_c3 <- rotasym::c_vMF(p = 3, kappa = kappa, log = TRUE)

  vapply(t_values, function(t) {
    threshold <- if (identical(distance_type, "chordal")) {
      1 - (t^2) / 2
    } else {
      cos(t)
    }

    if (!is.finite(threshold) || threshold <= -1) {
      return(1)
    }
    if (threshold >= 1) {
      return(0)
    }

    density_u <- function(u) {
      lambda <- kappa * tangential_norm * sqrt(pmax(0, 1 - u^2))
      log_i0 <- ifelse(
        lambda <= 0,
        0,
        log(besselI(lambda, nu = 0, expon.scaled = TRUE)) + lambda
      )
      exp(log_c3 + kappa * rho * u + log(2 * pi) + log_i0)
    }

    cdf_at_threshold <- integrate(
      f = density_u,
      lower = -1 + 1e-8,
      upper = threshold,
      rel.tol = rel.tol,
      abs.tol = abs.tol
    )$value

    1 - cdf_at_threshold
  }, numeric(1))
}
build_vmf_s2_u_grid <- function(n_u = 4097L) {
  n_u <- as.integer(n_u)
  if (!is.finite(n_u) || n_u < 17L || n_u %% 2L != 1L) {
    stop("`n_u` must be an odd integer >= 17.")
  }

  u_grid <- seq(-1, 1, length.out = n_u)
  list(
    u = u_grid,
    du = u_grid[[2]] - u_grid[[1]],
    sqrt_one_minus_u2 = sqrt(pmax(0, 1 - u_grid^2))
  )
}
vmf_s2_projected_density_matrix <- function(rho,
                                            kappa,
                                            u_grid_object) {
  rho <- as.numeric(rho)
  rho <- pmin(pmax(rho, -1), 1)
  n_u <- length(u_grid_object$u)
  n_rho <- length(rho)

  if (n_rho == 0L) {
    return(matrix(0, nrow = n_u, ncol = 0L))
  }

  if (kappa <= 1e-12) {
    return(matrix(0.5, nrow = n_u, ncol = n_rho))
  }

  sqrt_one_minus_rho2 <- sqrt(pmax(0, 1 - rho^2))
  lambda <- kappa * outer(u_grid_object$sqrt_one_minus_u2, sqrt_one_minus_rho2)
  log_i0 <- matrix(0, nrow = n_u, ncol = n_rho)
  positive_lambda <- lambda > 1e-14
  log_i0[positive_lambda] <- log(besselI(
    lambda[positive_lambda],
    nu = 0,
    expon.scaled = TRUE
  )) + lambda[positive_lambda]

  log_sinh_kappa <- kappa + log1p(-exp(-2 * kappa)) - log(2)
  log_const <- log(kappa) - log(2) - log_sinh_kappa
  log_density <- log_const +
    kappa * outer(u_grid_object$u, rho) +
    log_i0

  exp(log_density)
}
build_vmf_s2_projected_cdf_matrix <- function(rho,
                                              kappa,
                                              n_u = 4097L) {
  u_grid_object <- build_vmf_s2_u_grid(n_u)
  density_matrix <- vmf_s2_projected_density_matrix(
    rho = rho,
    kappa = kappa,
    u_grid_object = u_grid_object
  )
  n_u <- nrow(density_matrix)
  n_rho <- ncol(density_matrix)

  if (n_rho == 0L) {
    return(list(
      u = u_grid_object$u,
      cdf = matrix(0, nrow = n_u, ncol = 0L)
    ))
  }

  increments <- 0.5 * (
    density_matrix[-1, , drop = FALSE] +
      density_matrix[-n_u, , drop = FALSE]
  ) * u_grid_object$du
  cdf_matrix <- rbind(rep(0, n_rho), apply(increments, 2, cumsum))
  cdf_matrix <- sweep(cdf_matrix, 2, cdf_matrix[n_u, ], "/", check.margin = FALSE)
  cdf_matrix[n_u, ] <- 1

  list(
    u = u_grid_object$u,
    cdf = cdf_matrix
  )
}
interpolate_vmf_s2_upper_tail <- function(cdf_values,
                                          u_grid,
                                          a_values) {
  a_values <- pmin(pmax(as.numeric(a_values), -1), 1)
  n_u <- length(u_grid)
  if (n_u < 2L) {
    stop("`u_grid` must have length at least 2.")
  }

  du <- u_grid[[2]] - u_grid[[1]]
  scaled_position <- (a_values - u_grid[[1]]) / du
  left_index <- floor(scaled_position) + 1L
  left_index <- pmin(pmax(left_index, 1L), n_u - 1L)
  lambda <- scaled_position - (left_index - 1L)
  lambda <- pmin(pmax(lambda, 0), 1)

  cdf_left <- cdf_values[left_index]
  cdf_right <- cdf_values[left_index + 1L]
  cdf_at_a <- (1 - lambda) * cdf_left + lambda * cdf_right

  cdf_at_a[a_values <= -1] <- 0
  cdf_at_a[a_values >= 1] <- 1
  1 - cdf_at_a
}
distance_profile_vmf_s2_grid <- function(omega_grid,
                                         mu,
                                         kappa,
                                         t_grid,
                                         distance_type = "geodesic",
                                         n_u = 4097L) {
  omega_grid <- as.matrix(omega_grid)
  mu <- as.numeric(mu)

  if (ncol(omega_grid) != 3L || length(mu) != 3L) {
    stop("`omega_grid` must have 3 columns and `mu` must have length 3.")
  }

  omega_norms <- sqrt(rowSums(omega_grid^2))
  if (any(omega_norms <= 0)) {
    stop("`omega_grid` must contain nonzero rows.")
  }
  omega_grid <- omega_grid / omega_norms
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t_grid <- as.numeric(t_grid)

  rho <- as.numeric(omega_grid %*% mu)
  cdf_object <- build_vmf_s2_projected_cdf_matrix(
    rho = rho,
    kappa = kappa,
    n_u = n_u
  )
  a_values <- if (identical(distance_type, "geodesic")) {
    cos(t_grid)
  } else {
    1 - (t_grid^2) / 2
  }
  a_values <- pmin(pmax(a_values, -1), 1)

  output <- t(vapply(seq_along(rho), function(i) {
    interpolate_vmf_s2_upper_tail(
      cdf_values = cdf_object$cdf[, i],
      u_grid = cdf_object$u,
      a_values = a_values
    )
  }, numeric(length(t_grid))))

  output[, t_grid <= 0] <- 0
  if (identical(distance_type, "geodesic")) {
    output[, t_grid >= pi] <- 1
  } else {
    output[, t_grid >= 2] <- 1
  }

  output
}
vmf_s2_legendre_coefficients <- function(kappa, l_max) {
  kappa <- as.numeric(kappa)
  l_max <- as.integer(l_max)

  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite nonnegative scalar.")
  }
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  key <- paste0(formatC(kappa, digits = 17, format = "fg"), "::", l_max)
  if (exists(key, envir = vmf_s2_legendre_cache, inherits = FALSE)) {
    return(get(key, envir = vmf_s2_legendre_cache, inherits = FALSE))
  }

  coeffs <- numeric(l_max + 1L)
  coeffs[[1L]] <- 1
  if (l_max > 0L && kappa > 0) {
    ell <- 0:l_max
    ive <- besselI(kappa, nu = ell + 0.5, expon.scaled = TRUE)
    i_ell <- sqrt(pi / (2 * kappa)) * exp(kappa) * ive
    coeffs <- (kappa / sinh(kappa)) * (2 * ell + 1) * i_ell
    coeffs[[1L]] <- 1
  }

  assign(key, coeffs, envir = vmf_s2_legendre_cache)
  coeffs
}
select_vmf_s2_legendre_l_max <- function(kappa,
                                         tail_tol = 1e-10,
                                         min_l = 10L,
                                         max_l = 200L,
                                         tail_length = 5L) {
  kappa <- as.numeric(kappa)
  tail_tol <- as.numeric(tail_tol)
  min_l <- as.integer(min_l)
  max_l <- as.integer(max_l)
  tail_length <- as.integer(tail_length)

  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite nonnegative scalar.")
  }
  if (length(tail_tol) != 1L || !is.finite(tail_tol) || tail_tol <= 0) {
    stop("`tail_tol` must be a strictly positive scalar.")
  }
  if (length(min_l) != 1L || !is.finite(min_l) || min_l < 0L) {
    stop("`min_l` must be a nonnegative integer.")
  }
  if (length(max_l) != 1L || !is.finite(max_l) || max_l < min_l) {
    stop("`max_l` must be an integer not smaller than `min_l`.")
  }
  if (length(tail_length) != 1L || !is.finite(tail_length) || tail_length < 1L) {
    stop("`tail_length` must be a strictly positive integer.")
  }

  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = max_l + tail_length)
  abs_coeffs <- abs(coeffs)

  for (l_max_candidate in min_l:max_l) {
    tail_idx <- (l_max_candidate + 2L):(l_max_candidate + tail_length + 1L)
    if (max(abs_coeffs[tail_idx]) <= tail_tol) {
      return(l_max_candidate)
    }
  }

  max_l
}
distance_profile_vmf_s2_legendre <- function(omega,
                                             mu,
                                             kappa,
                                             t_values,
                                             distance_type = "geodesic",
                                             l_max = NULL,
                                             tail_tol = 1e-10) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3L || length(mu) != 3L) {
    stop("`omega` and `mu` must both have length 3 for S^2 Legendre vMF distance profiles.")
  }

  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  if (is.null(l_max)) {
    l_max <- select_vmf_s2_legendre_l_max(kappa = kappa, tail_tol = tail_tol)
  }
  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = l_max)

  rotational_profile_legendre(
    t = t_values,
    omega = omega,
    mu = mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}
distance_profile_vmf_s2_legendre_grid <- function(omega_grid,
                                                  mu,
                                                  kappa,
                                                  t_grid,
                                                  distance_type = "geodesic",
                                                  l_max = NULL,
                                                  tail_tol = 1e-10) {
  omega_grid <- as.matrix(omega_grid)
  mu <- as.numeric(mu)
  if (ncol(omega_grid) != 3L || length(mu) != 3L) {
    stop("`omega_grid` must have 3 columns and `mu` must have length 3.")
  }

  omega_norms <- sqrt(rowSums(omega_grid^2))
  if (any(omega_norms <= 0)) {
    stop("`omega_grid` must contain nonzero rows.")
  }
  omega_grid <- omega_grid / omega_norms
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  if (is.null(l_max)) {
    l_max <- select_vmf_s2_legendre_l_max(kappa = kappa, tail_tol = tail_tol)
  }
  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = l_max)

  rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}
distance_profile_vmf_s2_legendre_cvm_grid <- function(X,
                                                      mu,
                                                      kappa,
                                                      l_max = NULL,
                                                      tail_tol = 1e-10) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)

  if (ncol(X) != 3L || length(mu) != 3L) {
    stop("`X` must have 3 columns and `mu` must have length 3.")
  }

  X_norms <- sqrt(rowSums(X^2))
  if (any(X_norms <= 0)) {
    stop("`X` must contain nonzero rows.")
  }
  X <- X / X_norms
  mu <- mu / sqrt(sum(mu^2))
  if (is.null(l_max)) {
    l_max <- select_vmf_s2_legendre_l_max(kappa = kappa, tail_tol = tail_tol)
  }
  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = l_max)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  r_values <- as.numeric(X %*% mu)

  1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )
}
distance_profile_vmf_s2_cvm_grid <- function(X,
                                             mu,
                                             kappa,
                                             n_u = 4097L) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)

  if (ncol(X) != 3L || length(mu) != 3L) {
    stop("`X` must have 3 columns and `mu` must have length 3.")
  }

  X_norms <- sqrt(rowSums(X^2))
  if (any(X_norms <= 0)) {
    stop("`X` must contain nonzero rows.")
  }
  X <- X / X_norms
  mu <- mu / sqrt(sum(mu^2))

  rho <- as.numeric(X %*% mu)
  a_matrix <- X %*% t(X)
  a_matrix <- pmin(pmax(a_matrix, -1), 1)

  cdf_object <- build_vmf_s2_projected_cdf_matrix(
    rho = rho,
    kappa = kappa,
    n_u = n_u
  )

  t(vapply(seq_along(rho), function(i) {
    interpolate_vmf_s2_upper_tail(
      cdf_values = cdf_object$cdf[, i],
      u_grid = cdf_object$u,
      a_values = a_matrix[i, ]
    )
  }, numeric(nrow(X))))
}
validate_vmf_s2_grid_profile <- function(n_checks = 200,
                                         kappa_values = c(0.5, 2, 5),
                                         n_u = 4097L,
                                         distance_type = "geodesic",
                                         seed = 1) {
  set.seed(seed)
  errors <- numeric(0)

  for (kappa in kappa_values) {
    mu <- stats::rnorm(3)
    mu <- mu / sqrt(sum(mu^2))
    omega_grid <- matrix(stats::rnorm(n_checks * 3), ncol = 3)
    omega_grid <- omega_grid / sqrt(rowSums(omega_grid^2))
    t_max <- if (identical(distance_type, "geodesic")) pi else 2
    t_grid <- stats::runif(n_checks, 0, t_max)

    fast_matrix <- distance_profile_vmf_s2_grid(
      omega_grid = omega_grid,
      mu = mu,
      kappa = kappa,
      t_grid = t_grid,
      distance_type = distance_type,
      n_u = n_u
    )
    fast_values <- diag(fast_matrix)
    exact_values <- vapply(seq_len(n_checks), function(i) {
      distance_profile_vmf_s2_integral(
        omega = omega_grid[i, ],
        mu = mu,
        kappa = kappa,
        t_values = t_grid[i],
        distance_type = distance_type
      )
    }, numeric(1))

    errors <- c(errors, abs(fast_values - exact_values))
  }

  c(
    max_error = max(errors),
    mean_error = mean(errors),
    q95_error = unname(stats::quantile(errors, probs = 0.95, names = FALSE, type = 8))
  )
}
sphere_surface_area <- function(q) {
  q <- as.integer(q)
  if (length(q) != 1L || !is.finite(q) || q < 0L) {
    stop("`q` must be a nonnegative integer.")
  }

  2 * pi^((q + 1) / 2) / gamma((q + 1) / 2)
}
jp_normalize_unit_vector <- function(x, arg_name = "`x`", min_length = 3L) {
  x <- as.numeric(x)
  if (length(x) < min_length) {
    stop(sprintf("%s must have length at least %d.", arg_name, min_length))
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s must be finite.", arg_name))
  }

  x_norm <- sqrt(sum(x^2))
  if (!is.finite(x_norm) || x_norm <= 0) {
    stop(sprintf("%s must have strictly positive norm.", arg_name))
  }

  x / x_norm
}
jp_normalize_unit_matrix <- function(x, arg_name = "`x`", min_ncol = 3L) {
  if (is.vector(x)) {
    x <- matrix(as.numeric(x), nrow = 1L)
  } else {
    x <- as.matrix(x)
  }

  if (nrow(x) == 0L || ncol(x) < min_ncol) {
    stop(sprintf("%s must be a non-empty matrix with at least %d columns.", arg_name, min_ncol))
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s must be finite.", arg_name))
  }

  norms <- sqrt(rowSums(x^2))
  if (any(!is.finite(norms)) || any(norms <= 0)) {
    stop(sprintf("%s must contain rows with strictly positive norm.", arg_name))
  }

  x / norms
}
jp_uniform_sphere <- function(n, ambient_dim) {
  z <- matrix(stats::rnorm(n * ambient_dim), nrow = n, ncol = ambient_dim)
  z / sqrt(rowSums(z^2))
}
jp_orthonormal_complement <- function(mu) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  q_complete <- qr.Q(qr(cbind(mu, diag(length(mu)))), complete = TRUE)
  q_complete[, -1L, drop = FALSE]
}
jp_normalize_probability_weights <- function(weights, n_expected) {
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

  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("`weights` must have strictly positive sum.")
  }

  weights / total_weight
}
small_circle_validate_parameters <- function(mu, kappa, nu, allow_negative_nu = FALSE) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("The Small Circle implementation currently supports only S^2, i.e. `mu` of length 3.")
  }

  kappa <- as.numeric(kappa)
  nu <- as.numeric(nu)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite scalar in [0, Inf).")
  }
  if (length(nu) != 1L || !is.finite(nu)) {
    stop("`nu` must be a finite scalar.")
  }

  if (!allow_negative_nu && (nu < 0 || nu >= 1)) {
    stop("`nu` must lie in [0, 1).")
  }
  if (allow_negative_nu && abs(nu) >= 1) {
    stop("`nu` must lie in (-1, 1) before canonicalization.")
  }

  list(mu = mu, kappa = kappa, nu = nu)
}
small_circle_canonicalize_theta <- function(mu, nu, tol = 1e-12) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  nu <- as.numeric(nu)
  if (length(nu) != 1L || !is.finite(nu) || abs(nu) >= 1) {
    stop("`nu` must be a finite scalar in (-1, 1).")
  }

  if (nu < 0) {
    mu <- -mu
    nu <- -nu
  }
  if (abs(nu) <= tol) {
    nu <- 0
  }

  list(mu = mu, nu = nu)
}
small_circle_erf <- function(x) {
  2 * stats::pnorm(sqrt(2) * as.numeric(x)) - 1
}
small_circle_log_norm_constant <- function(kappa, nu) {
  # Keep the scalar-kappa small-circle model while allowing batched axial offsets.
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = 0)
  nu <- as.numeric(nu)
  if (length(nu) == 0L || any(!is.finite(nu)) || any(nu < 0 | nu >= 1)) {
    stop("`nu` must contain finite values in [0, 1).")
  }
  if (params$kappa <= 0) {
    return(rep(0, length(nu)))
  }

  root_kappa <- sqrt(params$kappa)
  erf_sum <- small_circle_erf(root_kappa * (1 - nu)) +
    small_circle_erf(root_kappa * (1 + nu))
  log(sqrt(pi) / (4 * root_kappa)) + log(erf_sum)
}
small_circle_norm_constant <- function(kappa, nu) {
  exp(small_circle_log_norm_constant(kappa = kappa, nu = nu))
}
small_circle_axis_density <- function(z, kappa, nu, log = FALSE) {
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = nu)
  z <- as.numeric(z)
  out <- rep(if (log) -Inf else 0, length(z))
  valid <- is.finite(z) & z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  if (params$kappa <= 0) {
    log_density <- rep(log(0.5), sum(valid))
  } else {
    log_density <- -log(2) - small_circle_log_norm_constant(params$kappa, params$nu) -
      params$kappa * (z[valid] - params$nu)^2
  }

  out[valid] <- if (log) log_density else exp(log_density)
  out
}
small_circle_axis_cdf <- function(z, kappa, nu) {
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = nu)
  z <- as.numeric(z)
  out <- numeric(length(z))
  out[z <= -1] <- 0
  out[z >= 1] <- 1
  active <- which(is.finite(z) & z > -1 & z < 1)
  if (length(active) == 0L) {
    return(out)
  }

  if (params$kappa <= 0) {
    out[active] <- (z[active] + 1) / 2
    return(out)
  }

  root_kappa <- sqrt(params$kappa)
  denominator <- small_circle_erf(root_kappa * (1 - params$nu)) +
    small_circle_erf(root_kappa * (1 + params$nu))
  numerator <- small_circle_erf(root_kappa * (z[active] - params$nu)) +
    small_circle_erf(root_kappa * (1 + params$nu))
  out[active] <- numerator / denominator
  pmin(pmax(out, 0), 1)
}
small_circle_gauss_legendre <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("`n` must be an integer >= 2.")
  }

  key <- as.character(n)
  if (exists(key, envir = small_circle_gauss_legendre_cache, inherits = FALSE)) {
    return(get(key, envir = small_circle_gauss_legendre_cache, inherits = FALSE))
  }

  beta <- seq_len(n - 1L) / sqrt(4 * seq_len(n - 1L)^2 - 1)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- beta
  jacobi[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- beta
  eig <- eigen(jacobi, symmetric = TRUE)
  order_idx <- order(eig$values)
  nodes <- eig$values[order_idx]
  weights <- 2 * (eig$vectors[1L, order_idx]^2)
  out <- list(nodes = nodes, weights = weights)
  assign(key, out, envir = small_circle_gauss_legendre_cache)
  out
}
small_circle_legendre_matrix <- function(x, l_max) {
  x <- pmin(pmax(as.numeric(x), -1), 1)
  l_max <- as.integer(l_max)
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  if (exists("distance_profile_backend_current", mode = "function") &&
      identical(distance_profile_backend_current(), "cpp")) {
    return(distance_profile_cpp_call("cpp_dp_legendre_matrix", x, l_max))
  }

  out <- matrix(0, nrow = length(x), ncol = l_max + 1L)
  out[, 1L] <- 1
  if (l_max == 0L) {
    return(out)
  }

  out[, 2L] <- x
  if (l_max == 1L) {
    return(out)
  }

  for (ell in seq_len(l_max - 1L)) {
    out[, ell + 2L] <- ((2 * ell + 1) * x * out[, ell + 1L] - ell * out[, ell]) / (ell + 1)
  }
  out
}
small_circle_legendre_coefficients <- function(kappa,
                                               nu,
                                               l_max = 150L,
                                               quad_n = 1000L,
                                               tol = 1e-10) {
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = nu)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (params$kappa <= 0) {
    coeffs <- c(1, rep(0, l_max))
    return(list(coefficients = coeffs, a0_error = 0, max_abs_nonzero = 0))
  }

  quad <- small_circle_gauss_legendre(quad_n)
  legendre_matrix <- small_circle_legendre_matrix(quad$nodes, l_max = l_max)
  h_values <- exp(-params$kappa * (quad$nodes - params$nu)^2 -
    small_circle_log_norm_constant(params$kappa, params$nu))
  raw_moments <- as.numeric(crossprod(legendre_matrix, quad$weights * h_values))
  ell <- 0:l_max
  coeffs <- ((2 * ell + 1) / 2) * raw_moments
  coeffs[[1L]] <- 1

  a0_error <- abs(((1 / 2) * raw_moments[[1L]]) - 1)
  if (a0_error > tol) {
    stop(sprintf("Small Circle Legendre coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error,
    max_abs_nonzero = if (l_max >= 1L) max(abs(coeffs[-1L])) else 0
  )
}
small_circle_projection_cdf_legendre <- function(x,
                                                 r,
                                                 coefficients,
                                                 enforce_bounds = TRUE) {
  x <- pmin(pmax(as.numeric(x), -1), 1)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != 1L || !is.finite(r)) {
    stop("`r` must be a finite scalar.")
  }

  out <- (x + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- as.numeric(small_circle_legendre_matrix(r, l_max = l_max))
  p_x <- small_circle_legendre_matrix(x, l_max = l_max + 1L)
  basis_matrix <- matrix(0, nrow = length(x), ncol = l_max)
  for (ell in seq_len(l_max)) {
    basis_matrix[, ell] <- (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1))
  }

  out <- out + as.numeric(basis_matrix %*% (coefficients[-1L] * p_r[-1L]))
  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}
small_circle_projection_cdf_legendre_matrix <- function(x_matrix,
                                                        r,
                                                        coefficients,
                                                        enforce_bounds = TRUE) {
  x_matrix <- as.matrix(x_matrix)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != nrow(x_matrix) || any(!is.finite(r))) {
    stop("`r` must be a finite vector of length nrow(`x_matrix`).")
  }

  if (exists("distance_profile_backend_current", mode = "function") &&
      identical(distance_profile_backend_current(), "cpp")) {
    return(distance_profile_cpp_call(
      "cpp_dp_projection_cdf_legendre_matrix",
      x_matrix,
      r,
      coefficients,
      isTRUE(enforce_bounds)
    ))
  }

  x_matrix <- pmin(pmax(x_matrix, -1), 1)
  out <- (x_matrix + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- small_circle_legendre_matrix(r, l_max = l_max)
  p_x <- small_circle_legendre_matrix(as.numeric(x_matrix), l_max = l_max + 1L)

  for (ell in seq_len(l_max)) {
    basis_ell <- matrix(
      (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1)),
      nrow = nrow(x_matrix),
      ncol = ncol(x_matrix)
    )
    out <- out + sweep(
      basis_ell,
      1L,
      coefficients[[ell + 1L]] * p_r[, ell + 1L],
      FUN = "*"
    )
  }

  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}
small_circle_projection_kernel <- function(c_thresholds, z_nodes, r) {
  c_thresholds <- pmin(pmax(as.numeric(c_thresholds), -1), 1)
  z_nodes <- pmin(pmax(as.numeric(z_nodes), -1), 1)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  one_minus_r2 <- pmax(0, 1 - r^2)
  a_values <- r * z_nodes
  b_values <- sqrt(one_minus_r2) * sqrt(pmax(0, 1 - z_nodes^2))
  kernel <- matrix(0, nrow = length(c_thresholds), ncol = length(z_nodes))

  for (j in seq_along(z_nodes)) {
    if (b_values[[j]] <= 1e-15) {
      kernel[, j] <- as.numeric(a_values[[j]] >= c_thresholds)
    } else {
      ratio <- (c_thresholds - a_values[[j]]) / b_values[[j]]
      kernel[, j] <- ifelse(
        ratio <= -1,
        1,
        ifelse(ratio >= 1, 0, acos(pmin(pmax(ratio, -1), 1)) / pi)
      )
    }
  }

  kernel
}
small_circle_distance_profile_integral <- function(omega,
                                                   t_values,
                                                   mu,
                                                   kappa,
                                                   nu,
                                                   distance_type = c("geodesic", "chordal"),
                                                   quad_n = 1000L) {
  distance_type <- match.arg(distance_type)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  if (length(omega) != 3L) {
    stop("The Small Circle implementation currently supports only S^2.")
  }

  t_values <- as.numeric(t_values)
  out <- numeric(length(t_values))
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  if (params$kappa <= 0) {
    out[active] <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_values[active])) / 2
    } else {
      (t_values[active]^2) / 4
    }
    return(out)
  }

  r_value <- sum(omega * params$mu)
  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  if (abs(r_value - 1) <= 1e-12) {
    out[active] <- 1 - small_circle_axis_cdf(thresholds, kappa = params$kappa, nu = params$nu)
    return(pmin(pmax(out, 0), 1))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- small_circle_axis_cdf(-thresholds, kappa = params$kappa, nu = params$nu)
    return(pmin(pmax(out, 0), 1))
  }

  quad <- small_circle_gauss_legendre(as.integer(quad_n))
  density_z <- small_circle_axis_density(quad$nodes, kappa = params$kappa, nu = params$nu)
  kernel <- small_circle_projection_kernel(
    c_thresholds = thresholds,
    z_nodes = quad$nodes,
    r = r_value
  )
  out[active] <- as.numeric(kernel %*% (quad$weights * density_z))
  pmin(pmax(out, 0), 1)
}
small_circle_monotone_clip <- function(t_values, values, upper_bound) {
  if (length(values) <= 1L) {
    return(pmin(pmax(values, 0), 1))
  }

  order_idx <- order(t_values)
  sorted_values <- pmin(pmax(values[order_idx], 0), 1)
  sorted_values <- cummax(sorted_values)
  sorted_values[t_values[order_idx] <= 0] <- 0
  sorted_values[t_values[order_idx] >= upper_bound] <- 1
  out <- numeric(length(sorted_values))
  out[order_idx] <- sorted_values
  out
}
small_circle_profile_matrix_legendre <- function(t_grid,
                                                 omega_grid,
                                                 mu,
                                                 coeffs,
                                                 l_max = length(coeffs) - 1L,
                                                 distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(omega_grid))
  coeffs <- as.numeric(coeffs)
  l_max <- as.integer(l_max)
  if (l_max != length(coeffs) - 1L) {
    stop("`l_max` must match `length(coeffs) - 1`.")
  }

  t_grid <- as.numeric(t_grid)
  n_omega <- nrow(omega_grid)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- matrix(0, nrow = n_omega, ncol = length(t_grid))

  if (length(t_grid) == 0L) {
    return(out)
  }

  active <- which(is.finite(t_grid) & t_grid > 0 & t_grid < upper_bound)
  out[, t_grid >= upper_bound] <- 1
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_grid[active], distance_type = distance_type)
  threshold_matrix <- matrix(thresholds, nrow = n_omega, ncol = length(active), byrow = TRUE)
  out[, active] <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = threshold_matrix,
    r = as.numeric(omega_grid %*% mu),
    coefficients = coeffs
  )
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(n_omega)) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = t_grid,
      values = out[i, ],
      upper_bound = upper_bound
    )
  }
  out
}
small_circle_sample_profile_matrix_legendre <- function(X,
                                                        mu,
                                                        coeffs,
                                                        l_max = length(coeffs) - 1L) {
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(X))
  coeffs <- as.numeric(coeffs)
  l_max <- as.integer(l_max)
  if (l_max != length(coeffs) - 1L) {
    stop("`l_max` must match `length(coeffs) - 1`.")
  }

  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  out <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = as.numeric(X %*% mu),
    coefficients = coeffs
  )
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
distance_profile_small_circle <- function(omega,
                                          t_values,
                                          mu,
                                          kappa,
                                          nu,
                                          distance_type = c("geodesic", "chordal"),
                                          method = c("legendre", "integral"),
                                          l_max = 150L,
                                          quad_n = 1000L,
                                          tol = 1e-10,
                                          validate_against_integral = FALSE,
                                          validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  t_values <- as.numeric(t_values)

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (ncol(omega) != 3L) {
      stop("`omega` must have three columns for Small Circle profiles on S^2.")
    }
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }

    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_small_circle(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = params$mu,
        kappa = params$kappa,
        nu = params$nu,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  if (params$kappa <= 0) {
    out[active] <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_values[active])) / 2
    } else {
      (t_values[active]^2) / 4
    }
    return(out)
  }

  if (identical(method, "integral")) {
    return(small_circle_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = params$mu,
      kappa = params$kappa,
      nu = params$nu,
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  coeffs <- small_circle_legendre_coefficients(
    kappa = params$kappa,
    nu = params$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  cdf_values <- small_circle_projection_cdf_legendre(
    x = thresholds,
    r = sum(omega * params$mu),
    coefficients = coeffs
  )
  out[active] <- 1 - cdf_values
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral)) {
    integral_values <- small_circle_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = params$mu,
      kappa = params$kappa,
      nu = params$nu,
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - integral_values))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Small Circle Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}
distance_profile_small_circle_grid <- function(omega_grid,
                                               mu,
                                               kappa,
                                               nu,
                                               t_grid,
                                               distance_type = c("geodesic", "chordal"),
                                               method = c("legendre", "integral"),
                                               l_max = 150L,
                                               quad_n = 1000L,
                                               tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  t_grid <- as.numeric(t_grid)

  if (params$kappa <= 0) {
    base_profile <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_grid)) / 2
    } else {
      (t_grid^2) / 4
    }
    return(matrix(base_profile, nrow = nrow(omega_grid), ncol = length(t_grid), byrow = TRUE))
  }

  if (identical(method, "integral")) {
    quad <- small_circle_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * small_circle_axis_density(
      quad$nodes,
      kappa = params$kappa,
      nu = params$nu
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = params$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
    r_values <- as.numeric(omega_grid %*% params$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - small_circle_axis_cdf(
        thresholds,
        kappa = params$kappa,
        nu = params$nu
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- small_circle_axis_cdf(
        -thresholds,
        kappa = params$kappa,
        nu = params$nu
      )
    }
    return(out)
  }

  coeffs <- small_circle_legendre_coefficients(
    kappa = params$kappa,
    nu = params$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  small_circle_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = params$mu,
    coeffs = coeffs,
    l_max = l_max,
    distance_type = distance_type
  )
}
distance_profile_small_circle_cvm_grid <- function(X,
                                                   mu,
                                                   kappa,
                                                   nu,
                                                   method = c("legendre", "integral"),
                                                   l_max = 150L,
                                                   quad_n = 1000L,
                                                   tol = 1e-10) {
  method <- match.arg(method)
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)

  if (params$kappa <= 0) {
    dot_products <- pmin(pmax(X %*% t(X), -1), 1)
    return((1 - dot_products) / 2)
  }

  if (identical(method, "integral")) {
    quad <- small_circle_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * small_circle_axis_density(
      quad$nodes,
      kappa = params$kappa,
      nu = params$nu
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = params$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- small_circle_legendre_coefficients(
    kappa = params$kappa,
    nu = params$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  small_circle_sample_profile_matrix_legendre(
    X = X,
    mu = params$mu,
    coeffs = coeffs,
    l_max = l_max
  )
}
r_sph_small_circle <- function(n, mu, kappa, nu, check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  if (params$kappa <= 0) {
    x <- jp_uniform_sphere(n = n, ambient_dim = 3L)
    if (isTRUE(check)) {
      expect_norms <- sqrt(rowSums(x^2))
      if (max(abs(expect_norms - 1)) > 1e-8) {
        stop("Small Circle uniform sampler returned non-unit vectors.")
      }
    }
    return(x)
  }

  sigma <- 1 / sqrt(2 * params$kappa)
  lower_prob <- stats::pnorm((-1 - params$nu) / sigma)
  upper_prob <- stats::pnorm((1 - params$nu) / sigma)
  uniforms <- stats::runif(n)
  z <- params$nu + sigma * stats::qnorm(lower_prob + uniforms * (upper_prob - lower_prob))
  z <- pmin(pmax(z, -1), 1)

  phi <- stats::runif(n, min = 0, max = 2 * pi)
  basis <- jp_orthonormal_complement(params$mu)
  tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
  radial <- sqrt(pmax(0, 1 - z^2))
  x <- tcrossprod(z, params$mu) + sweep(tangent, 1, radial, "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Small Circle sampler returned non-unit vectors.")
    }
  }

  x
}
small_circle_weighted_loglik_s2 <- function(mu,
                                            kappa,
                                            nu,
                                            x,
                                            prob_weights = NULL) {
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  projections <- as.numeric(x %*% params$mu)
  -small_circle_log_norm_constant(params$kappa, params$nu) -
    params$kappa * sum(prob_weights * (projections - params$nu)^2)
}
small_circle_start_theta_s2 <- function(x,
                                        weights = NULL,
                                        nu_min = 1e-6,
                                        nu_eps = 1e-6,
                                        kappa_min = 1e-8,
                                        kappa_max = 1e6) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  xbar <- colSums(x * prob_weights)
  weighted_x <- sweep(x, 1, sqrt(prob_weights), "*")
  s_matrix <- crossprod(weighted_x)
  sigma_matrix <- s_matrix - tcrossprod(xbar)
  eig <- eigen(sigma_matrix, symmetric = TRUE)
  mu0 <- eig$vectors[, which.min(eig$values)]
  nu0 <- sum(mu0 * xbar)
  canonical <- small_circle_canonicalize_theta(mu0, nu0)
  mu0 <- canonical$mu
  nu0 <- canonical$nu
  nu0 <- min(max(nu0, nu_min), 1 - nu_eps)

  t3 <- max(min(eig$values), .Machine$double.eps)
  kappa0 <- min(max(1 / (2 * t3), kappa_min), kappa_max)

  list(mu = mu0, kappa = kappa0, nu = nu0)
}
small_circle_mle_s2_weighted <- function(x,
                                         weights = NULL,
                                         control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  nu_eps <- as.numeric(control$small_circle_nu_eps %||% 1e-6)
  kappa_min <- as.numeric(control$small_circle_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_kappa_max %||% 1e6)
  nu_optim <- match.arg(
    as.character(control$small_circle_nu_optim %||% "logistic"),
    choices = c("logistic", "box")
  )
  theta_start <- control$small_circle_mle_start_theta %||% control$theta_start %||% control$jp_mle_start_theta %||% NULL

  if (is.null(theta_start)) {
    theta_start <- small_circle_start_theta_s2(
      x = x,
      weights = prob_weights,
      nu_min = nu_eps,
      nu_eps = nu_eps,
      kappa_min = kappa_min,
      kappa_max = kappa_max
    )
  } else {
    theta_start <- small_circle_validate_parameters(
      mu = theta_start$mu,
      kappa = theta_start$kappa,
      nu = theta_start$nu,
      allow_negative_nu = TRUE
    )
    canonical <- small_circle_canonicalize_theta(theta_start$mu, theta_start$nu)
    theta_start$mu <- canonical$mu
    theta_start$nu <- min(max(canonical$nu, nu_eps), 1 - nu_eps)
    theta_start$kappa <- min(max(theta_start$kappa, kappa_min), kappa_max)
  }

  objective <- function(par) {
    if (identical(nu_optim, "logistic")) {
      kappa_value <- min(max(log1p(exp(par[[1L]])), kappa_min), kappa_max)
      nu_value <- small_circle_logistic_bounded(par[[2L]], upper = 1 - nu_eps)
    } else {
      kappa_value <- min(max(par[[1L]], kappa_min), kappa_max)
      nu_value <- min(max(par[[2L]], 0), 1 - nu_eps)
    }

    mu_raw <- par[3:5]
    mu_norm <- sqrt(sum(mu_raw^2))
    if (!is.finite(mu_norm) || mu_norm <= 0) {
      return(.Machine$double.xmax / 100)
    }

    mu_value <- mu_raw / mu_norm
    value <- -small_circle_weighted_loglik_s2(
      mu = mu_value,
      kappa = kappa_value,
      nu = nu_value,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      return(.Machine$double.xmax / 100)
    }
    value
  }

  optim_control <- control$small_circle_optim_control %||% list(maxit = 500L, reltol = 1e-10)
  if (identical(nu_optim, "logistic")) {
    opt <- stats::optim(
      par = c(
        log(expm1(theta_start$kappa)),
        small_circle_inverse_logistic_bounded(theta_start$nu, upper = 1 - nu_eps),
        theta_start$mu
      ),
      fn = objective,
      method = control$small_circle_optim_method %||% "BFGS",
      control = optim_control
    )
    kappa_hat <- min(max(log1p(exp(opt$par[[1L]])), kappa_min), kappa_max)
    nu_hat <- small_circle_logistic_bounded(opt$par[[2L]], upper = 1 - nu_eps)
  } else {
    optim_method_box <- control$small_circle_optim_method %||% "L-BFGS-B"
    if (identical(optim_method_box, "L-BFGS-B") && !is.null(optim_control$reltol)) {
      reltol_value <- as.numeric(optim_control$reltol)
      if (length(reltol_value) == 1L && is.finite(reltol_value) && reltol_value > 0) {
        optim_control$factr <- reltol_value / .Machine$double.eps
      }
      optim_control$reltol <- NULL
    }
    lower_bounds <- c(kappa_min, 0, rep(-Inf, 3L))
    upper_bounds <- c(kappa_max, 1 - nu_eps, rep(Inf, 3L))
    opt <- stats::optim(
      par = c(theta_start$kappa, theta_start$nu, theta_start$mu),
      fn = objective,
      method = optim_method_box,
      lower = lower_bounds,
      upper = upper_bounds,
      control = optim_control
    )
    kappa_hat <- min(max(opt$par[[1L]], kappa_min), kappa_max)
    nu_hat <- min(max(opt$par[[2L]], 0), 1 - nu_eps)
  }

  mu_hat <- opt$par[3:5] / sqrt(sum(opt$par[3:5]^2))
  canonical <- small_circle_canonicalize_theta(mu_hat, nu_hat)

  list(
    mu = canonical$mu,
    kappa = kappa_hat,
    nu = canonical$nu,
    loglik = -opt$value,
    opt = opt,
    weighted_mle = TRUE,
    start_theta = theta_start,
    nu_optim = nu_optim
  )
}
rotational_logsumexp2 <- function(log_x, log_y) {
  m <- pmax(log_x, log_y)
  m + log(exp(log_x - m) + exp(log_y - m))
}
rotational_clamp_unit_interval <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}
rotational_bounded_weight <- function(eta, weight_eps = 1e-6) {
  pmin(pmax(stats::plogis(eta), weight_eps), 1 - weight_eps)
}
rotational_positive_parameter <- function(log_value,
                                          lower = 1e-6,
                                          upper = 1e6) {
  pmin(pmax(exp(log_value), lower), upper)
}
rotational_weighted_mean <- function(x, weights) {
  sum(as.numeric(weights) * as.numeric(x))
}
rotational_weighted_variance <- function(x, weights, center = NULL) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  if (is.null(center)) {
    center <- rotational_weighted_mean(x, weights)
  }
  sum(weights * (x - center)^2)
}
rotational_weighted_quantile <- function(x, weights, prob) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  prob <- as.numeric(prob)

  if (length(prob) != 1L || !is.finite(prob) || prob < 0 || prob > 1) {
    stop("`prob` must be a finite scalar in [0, 1].")
  }

  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- weights[ord]
  cum_weights <- cumsum(w_sorted)
  total <- cum_weights[[length(cum_weights)]]
  x_sorted[[which(cum_weights >= prob * total)[[1L]]]]
}
rotational_weighted_covariance_matrix <- function(x, weights) {
  x <- as.matrix(x)
  weights <- as.numeric(weights)
  center <- colSums(x * weights)
  centered <- sweep(x, 2L, center, FUN = "-")
  crossprod(sweep(centered, 1L, sqrt(weights), FUN = "*"))
}
rotational_unit_vector_fallback <- function(x,
                                            fallback = c(0, 0, 1),
                                            tol = 1e-12) {
  x <- as.numeric(x)
  norm_x <- sqrt(sum(x^2))
  if (!is.finite(norm_x) || norm_x <= tol) {
    fallback <- as.numeric(fallback)
    fallback / sqrt(sum(fallback^2))
  } else {
    x / norm_x
  }
}
rotational_unique_mu_candidates <- function(mu_candidates,
                                            tol = 1e-8) {
  out <- list()
  for (mu in mu_candidates) {
    mu_vec <- rotational_unit_vector_fallback(mu)
    keep <- TRUE
    if (length(out) > 0L) {
      for (existing in out) {
        if (max(abs(mu_vec - existing)) <= tol || max(abs(mu_vec + existing)) <= tol) {
          keep <- FALSE
          break
        }
      }
    }
    if (keep) {
      out[[length(out) + 1L]] <- mu_vec
    }
  }
  out
}
rotational_gauss_legendre <- function(n) {
  small_circle_gauss_legendre(n)
}
rotational_legendre_matrix <- function(x, l_max) {
  small_circle_legendre_matrix(x, l_max)
}
rotational_projection_cdf_legendre <- function(x,
                                               r,
                                               coefficients,
                                               enforce_bounds = TRUE) {
  x <- pmin(pmax(as.numeric(x), -1), 1)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != 1L || !is.finite(r)) {
    stop("`r` must be a finite scalar.")
  }

  out <- (x + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- as.numeric(rotational_legendre_matrix(r, l_max = l_max))
  p_x <- rotational_legendre_matrix(x, l_max = l_max + 1L)
  basis_matrix <- matrix(0, nrow = length(x), ncol = l_max)
  for (ell in seq_len(l_max)) {
    basis_matrix[, ell] <- (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1))
  }

  out <- out + as.numeric(basis_matrix %*% (coefficients[-1L] * p_r[-1L]))
  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}
rotational_projection_cdf_legendre_matrix <- function(x_matrix,
                                                      r,
                                                      coefficients,
                                                      enforce_bounds = TRUE) {
  x_matrix <- as.matrix(x_matrix)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != nrow(x_matrix) || any(!is.finite(r))) {
    stop("`r` must be a finite vector of length nrow(`x_matrix`).")
  }

  x_matrix <- pmin(pmax(x_matrix, -1), 1)
  out <- (x_matrix + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- rotational_legendre_matrix(r, l_max = l_max)
  p_x <- rotational_legendre_matrix(as.numeric(x_matrix), l_max = l_max + 1L)

  for (ell in seq_len(l_max)) {
    basis_ell <- matrix(
      (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1)),
      nrow = nrow(x_matrix),
      ncol = ncol(x_matrix)
    )
    out <- out + sweep(
      basis_ell,
      1L,
      coefficients[[ell + 1L]] * p_r[, ell + 1L],
      FUN = "*"
    )
  }

  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}
rotational_profile_legendre <- function(t,
                                        omega,
                                        mu,
                                        coeffs,
                                        Lmax = length(coeffs) - 1L,
                                        distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = length(omega))
  coeffs <- as.numeric(coeffs)
  Lmax <- as.integer(Lmax)
  if (Lmax != length(coeffs) - 1L) {
    stop("`Lmax` must match `length(coeffs) - 1`.")
  }

  t <- as.numeric(t)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t))
  out[t <= 0] <- 0
  out[t >= upper_bound] <- 1
  active <- which(is.finite(t) & t > 0 & t < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t[active], distance_type = distance_type)
  out[active] <- 1 - rotational_projection_cdf_legendre(
    x = thresholds,
    r = sum(omega * mu),
    coefficients = coeffs
  )
  small_circle_monotone_clip(t_values = t, values = out, upper_bound = upper_bound)
}
rotational_profile_matrix_legendre <- function(t_grid,
                                               omega_grid,
                                               mu,
                                               coeffs,
                                               Lmax = length(coeffs) - 1L,
                                               distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(omega_grid))
  coeffs <- as.numeric(coeffs)
  Lmax <- as.integer(Lmax)
  if (Lmax != length(coeffs) - 1L) {
    stop("`Lmax` must match `length(coeffs) - 1`.")
  }

  t_grid <- as.numeric(t_grid)
  n_omega <- nrow(omega_grid)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- matrix(0, nrow = n_omega, ncol = length(t_grid))

  if (length(t_grid) == 0L) {
    return(out)
  }

  active <- which(is.finite(t_grid) & t_grid > 0 & t_grid < upper_bound)
  out[, t_grid >= upper_bound] <- 1
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_grid[active], distance_type = distance_type)
  p_r <- rotational_legendre_matrix(as.numeric(omega_grid %*% mu), l_max = Lmax)
  if (Lmax == 0L) {
    out[, active] <- matrix((1 - thresholds) / 2, nrow = n_omega, ncol = length(active), byrow = TRUE)
    return(out)
  }

  p_x <- rotational_legendre_matrix(thresholds, l_max = Lmax + 1L)
  basis <- matrix(0, nrow = length(active), ncol = Lmax)
  for (ell in seq_len(Lmax)) {
    basis[, ell] <- (p_x[, ell] - p_x[, ell + 2L]) / (2 * (2 * ell + 1))
  }

  out[, active] <- matrix((1 - thresholds) / 2, nrow = n_omega, ncol = length(active), byrow = TRUE) +
    p_r[, -1L, drop = FALSE] %*% t(sweep(basis, 2L, coeffs[-1L], FUN = "*"))
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(n_omega)) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = t_grid,
      values = out[i, ],
      upper_bound = upper_bound
    )
  }
  out
}
rotational_distance_profile_integral <- function(omega,
                                                 t_values,
                                                 mu,
                                                 density_gz,
                                                 distance_type = c("geodesic", "chordal"),
                                                 quad_n = 1000L) {
  distance_type <- match.arg(distance_type)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = length(omega))
  if (!is.function(density_gz)) {
    stop("`density_gz` must be a function.")
  }

  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  quad <- rotational_gauss_legendre(as.integer(quad_n))
  density_z <- as.numeric(density_gz(quad$nodes))
  if (length(density_z) != length(quad$nodes) || any(!is.finite(density_z)) || any(density_z < 0)) {
    stop("`density_gz` must return finite nonnegative values of the same length as its input.")
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  kernel <- small_circle_projection_kernel(
    c_thresholds = thresholds,
    z_nodes = quad$nodes,
    r = sum(omega * mu)
  )
  out[active] <- as.numeric(kernel %*% (quad$weights * density_z))
  pmin(pmax(out, 0), 1)
}
rotational_sample_profile_matrix <- function(X,
                                             mu,
                                             profile_fun,
                                             ...) {
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  t(vapply(seq_len(nrow(X)), function(i) {
    as.numeric(profile_fun(
      omega = X[i, ],
      t_values = acos(dot_products[i, ]),
      mu = mu,
      ...
    ))
  }, numeric(nrow(X))))
}
