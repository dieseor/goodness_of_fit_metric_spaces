# Mathematical and numerical helpers: geometry.

euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}
hyperbolic_geodesic_distance_h2 <- function(x, y, tol = 1e-12) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (length(x) != 3L || length(y) != 3L) {
    stop("`x` and `y` must both have length 3.")
  }

  inner_xy <- minkowski_inner_product(x, y)
  acosh(pmax(-inner_xy, 1))
}
s1_chordal_half_width <- function(t) {
  t <- pmin(pmax(t, 0), 2)
  2 * asin(t / 2)
}
s1_event_segments_chordal <- function(omega, t, tol = 1e-12) {
  if (t <= tol) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }
  if (t >= 2 - tol) {
    return(matrix(c(0, 2 * pi), ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("start", "end"))))
  }
  theta <- circle_angle_from_point(omega)
  delta <- s1_chordal_half_width(t)
  s1_interval_to_segments(theta - delta, theta + delta, tol = tol)
}
distance_profile_jp <- function(omega,
                                t_values,
                                mu,
                                kappa,
                                psi,
                                distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  t_values <- as.numeric(t_values)

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = length(params$mu))
    if (ncol(omega) != length(params$mu)) {
      stop("`omega` and `mu` must have compatible ambient dimensions.")
    }
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }

    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_jp(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = params$mu,
        kappa = params$kappa,
        psi = params$psi,
        distance_type = distance_type
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = length(params$mu))
  if (length(omega) != length(params$mu)) {
    stop("`omega` and `mu` must have the same length.")
  }

  if (params$ambient_dim == 3L && params$psi != 0 && abs(params$alpha) <= 1e-15) {
    return(jp_uniform_s2_distance_profile(t_values, distance_type = distance_type))
  }

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # Theoretical profiles near the JP-vMF boundary are evaluated through the
    # exact vMF profile. This preserves numerical stability in S^2 at the cost
    # of replacing JP by its limit model when |kappa * psi| is small.
    if (params$ambient_dim == 3L) {
      return(theoretical_distance_profile_vmf_s2_fast(
        omega = omega,
        mu = params$mu,
        kappa = params$kappa,
        t_values = t_values,
        distance_type = distance_type
      ))
    }

    return(theoretical_distance_profile_vmf(
      omega = omega,
      mu = params$mu,
      kappa = params$kappa,
      t_values = t_values,
      distance_type = distance_type
    ))
  }

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

  rho <- pmin(pmax(sum(omega * params$mu), -1), 1)
  cdf_object <- build_jp_projected_cdf_matrix(
    rho = rho,
    q = params$q,
    alpha = params$alpha,
    beta = params$beta,
    n_u = 4097L,
    n_delta = 1025L
  )
  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  thresholds <- pmin(pmax(thresholds, -1), 1)
  cdf_interpolator <- jp_prepare_cdf_interpolator(
    x_grid = cdf_object$u,
    cdf_grid = cdf_object$cdf[, 1L]
  )
  output[active] <- jp_interpolate_upper_tail(
    threshold = thresholds,
    x_grid = cdf_object$u,
    cdf_grid = cdf_object$cdf[, 1L],
    interpolator = cdf_interpolator
  )

  output
}
distance_profile_beta_mixture2 <- function(omega,
                                                      t_values,
                                                      mu,
                                                      weight1,
                                                      alpha1,
                                                      beta1,
                                                      alpha2,
                                                      beta2,
                                                      distance_type = c("geodesic", "chordal"),
                                                      method = c("legendre", "integral"),
                                                      l_max = 150L,
                                                      quad_n = 1000L,
                                                      tol = 1e-10,
                                                      validate_against_integral = FALSE,
                                                      validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
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
      distance_profile_beta_mixture2(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
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
    out[active] <- 1 - beta_mixture2_cdf_y(
      y = (thresholds + 1) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- beta_mixture2_cdf_y(
      y = (1 - thresholds) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }

  if (identical(method, "integral")) {
    return(rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        beta_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  coeffs <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
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
        beta_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Beta-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}
distance_profile_beta_mixture2_grid <- function(omega_grid,
                                                           mu,
                                                           weight1,
                                                           alpha1,
                                                           beta1,
                                                           alpha2,
                                                           beta2,
                                                           t_grid,
                                                           distance_type = c("geodesic", "chordal"),
                                                           method = c("legendre", "integral"),
                                                           l_max = 150L,
                                                           quad_n = 1000L,
                                                           tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  if (identical(method, "integral")) {
    omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * beta_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    # At the two axial directions the projected law is exactly a beta
    # mixture CDF.  Gauss--Legendre quadrature converges slowly there because
    # the projection kernel becomes an indicator and beta shapes need not be
    # integer.  Use the closed form, consistently with the scalar evaluator.
    thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
    r_values <- as.numeric(omega_grid %*% theta$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - beta_mixture2_cdf_y(
        y = (thresholds + 1) / 2,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- beta_mixture2_cdf_y(
        y = (1 - thresholds) / 2,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2
      )
    }
    return(out)
  }

  coeffs <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients

  out <- rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = theta$mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )

  out
}
distance_profile_beta_mixture2_cvm_grid <- function(X,
                                                               mu,
                                                               weight1,
                                                               alpha1,
                                                               beta1,
                                                               alpha2,
                                                               beta2,
                                                               method = c("legendre", "integral"),
                                                               l_max = 150L,
                                                               quad_n = 1000L,
                                                               tol = 1e-10) {
  method <- match.arg(method)
  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)

  if (identical(method, "integral")) {
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * beta_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    out <- projection_sample_profile_matrix_integral(
      X = X,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    )
    r_values <- as.numeric(X %*% theta$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - beta_mixture2_cdf_y(
        y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- beta_mixture2_cdf_y(
        y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2
      )
    }
    return(pmin(pmax(out, 0), 1))
  }

  coeffs <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  out <- 1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )

  pos_idx <- which(abs(r_values - 1) <= 1e-12)
  if (length(pos_idx) > 0L) {
    out[pos_idx, ] <- 1 - beta_mixture2_cdf_y(
      y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
  }

  neg_idx <- which(abs(r_values + 1) <= 1e-12)
  if (length(neg_idx) > 0L) {
    out[neg_idx, ] <- beta_mixture2_cdf_y(
      y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
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
distance_profile_uniform_beta_mixture <- function(omega,
                                                  t_values,
                                                  mu,
                                                  weight_uniform,
                                                  alpha,
                                                  beta,
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
  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
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
      distance_profile_uniform_beta_mixture(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta,
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
    out[active] <- 1 - uniform_beta_mixture_cdf_y(
      y = (thresholds + 1) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- uniform_beta_mixture_cdf_y(
      y = (1 - thresholds) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }

  if (identical(method, "integral")) {
    return(rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        uniform_beta_mixture_density_gz(
          z = z,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  coeffs <- uniform_beta_mixture_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
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
        uniform_beta_mixture_density_gz(
          z = z,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Uniform-beta-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}
distance_profile_uniform_beta_mixture_grid <- function(omega_grid,
                                                       mu,
                                                       weight_uniform,
                                                       alpha,
                                                       beta,
                                                       t_grid,
                                                       distance_type = c("geodesic", "chordal"),
                                                       method = c("legendre", "integral"),
                                                       l_max = 150L,
                                                       quad_n = 1000L,
                                                       tol = 1e-10,
                                                       eps = 1e-12) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  if (identical(method, "integral")) {
    omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * uniform_beta_mixture_density_gz(
      z = quad$nodes,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
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
      out[pos_idx, ] <- 1 - uniform_beta_mixture_cdf_y(
        y = (thresholds + 1) / 2,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- uniform_beta_mixture_cdf_y(
        y = (1 - thresholds) / 2,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta
      )
    }
    return(out)
  }

  coeffs <- uniform_beta_mixture_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
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
distance_profile_uniform_beta_mixture_cvm_grid <- function(X,
                                                           mu,
                                                           weight_uniform,
                                                           alpha,
                                                           beta,
                                                           distance_matrix = NULL,
                                                           method = c("legendre", "integral"),
                                                           l_max = 150L,
                                                           quad_n = 1000L,
                                                           tol = 1e-10,
                                                           eps = 1e-12) {
  method <- match.arg(method)
  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  if (is.null(distance_matrix)) {
    dot_products <- pmin(pmax(X %*% t(X), -1), 1)
    distance_matrix <- acos(dot_products)
  } else {
    distance_matrix <- as.matrix(distance_matrix)
    if (!all(dim(distance_matrix) == c(nrow(X), nrow(X)))) {
      stop("`distance_matrix` must be an n x n matrix compatible with `X`.")
    }
    dot_products <- cos(pmin(pmax(distance_matrix, 0), pi))
  }

  if (identical(method, "integral")) {
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * uniform_beta_mixture_density_gz(
      z = quad$nodes,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      eps = eps
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- uniform_beta_mixture_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  out <- 1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )

  pos_idx <- which(abs(r_values - 1) <= 1e-12)
  if (length(pos_idx) > 0L) {
    out[pos_idx, ] <- 1 - uniform_beta_mixture_cdf_y(
      y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
  }

  neg_idx <- which(abs(r_values + 1) <= 1e-12)
  if (length(neg_idx) > 0L) {
    out[neg_idx, ] <- uniform_beta_mixture_cdf_y(
      y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
  }

  out <- pmin(pmax(out, 0), 1)
  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = distance_matrix[i, ],
      values = out[i, ],
      upper_bound = pi
    )
  }
  out
}
