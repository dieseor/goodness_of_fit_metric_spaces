# Deterministic profile derivatives for vMF and HvMF models.
#
# The derivatives in this file are with respect to the canonical Euclidean
# coordinates xi.  For HvMF, the model exponent is x' J xi with
# J = diag(-1, 1, ..., 1).

profile_derivative_cumtrapz <- function(grid, integrands) {
  grid <- as.numeric(grid)
  integrands <- as.matrix(integrands)
  if (length(grid) != nrow(integrands) || length(grid) < 2L) {
    stop("The quadrature grid and integrands have incompatible dimensions.")
  }
  increments <- 0.5 * diff(grid) *
    (integrands[-1L, , drop = FALSE] +
       integrands[-nrow(integrands), , drop = FALSE])
  rbind(rep.int(0, ncol(integrands)), apply(increments, 2L, cumsum))
}

profile_derivative_interpolate_columns <- function(grid, cumulative, xout) {
  grid <- as.numeric(grid)
  cumulative <- as.matrix(cumulative)
  xout <- as.numeric(xout)
  values <- vapply(seq_len(ncol(cumulative)), function(j) {
    stats::approx(
      x = grid,
      y = cumulative[, j],
      xout = xout,
      rule = 2L,
      ties = "ordered"
    )$y
  }, numeric(length(xout)))
  matrix(
    values,
    nrow = length(xout),
    ncol = ncol(cumulative),
    dimnames = list(NULL, colnames(cumulative))
  )
}

profile_derivative_nonnegative_square <- function(value,
                                                  scale,
                                                  label,
                                                  relative_tolerance = 1e-10) {
  value <- as.numeric(value)
  scale <- max(1, abs(as.numeric(scale)))
  if (!is.finite(value)) {
    stop(sprintf("The computed %s is not finite.", label))
  }
  if (value < -relative_tolerance * scale) {
    stop(sprintf(
      "The computed %s = %.17g is substantially negative; the supplied parameters are incompatible.",
      label,
      value
    ))
  }
  max(value, 0)
}

# Stable continuous extension of A_{q-1}(u) / u.
profile_derivative_bessel_ratio_over_argument <- function(u, q) {
  u <- as.numeric(u)
  q <- as.integer(q)
  if (length(q) != 1L || !is.finite(q) || q < 2L) {
    stop("`q` must be an integer of at least two.")
  }
  if (any(!is.finite(u)) || any(u < 0)) {
    stop("Bessel-ratio arguments must be finite and nonnegative.")
  }

  output <- numeric(length(u))
  small <- u <= 1e-3
  if (any(small)) {
    u2 <- u[small]^2
    output[small] <- 1 / q - u2 / (q^2 * (q + 2))
  }
  if (any(!small)) {
    u_regular <- u[!small]
    numerator <- besselI(u_regular, nu = q / 2, expon.scaled = TRUE)
    denominator <- besselI(
      u_regular,
      nu = q / 2 - 1,
      expon.scaled = TRUE
    )
    ratio <- numerator / denominator
    if (any(!is.finite(ratio)) || any(ratio <= 0)) {
      stop("Could not evaluate a stable scaled-Bessel ratio.")
    }
    output[!small] <- ratio / u_regular
  }
  output
}

vmf_log_normalizing_constant_intrinsic <- function(q, kappa) {
  q <- as.integer(q)
  kappa <- as.numeric(kappa)
  if (length(q) != 1L || !is.finite(q) || q < 1L ||
      length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("Invalid vMF dimension or concentration.")
  }
  if (kappa == 0) {
    return(lgamma((q + 1) / 2) - log(2) - ((q + 1) / 2) * log(pi))
  }
  nu <- (q - 1) / 2
  scaled_bessel <- besselI(kappa, nu = nu, expon.scaled = TRUE)
  if (!is.finite(scaled_bessel) || scaled_bessel <= 0) {
    stop("Could not evaluate the scaled Bessel function for the vMF normalising constant.")
  }
  nu * log(kappa) - ((q + 1) / 2) * log(2 * pi) -
    log(scaled_bessel) - kappa
}

vmf_projected_density_canonical <- function(s, xi, omega) {
  s <- as.numeric(s)
  xi <- as.numeric(xi)
  omega <- as.numeric(omega)
  q <- length(xi) - 1L
  if (q < 2L || length(omega) != length(xi) ||
      any(!is.finite(c(xi, omega, s)))) {
    stop("Deterministic vMF derivatives require finite vectors on S^q with q >= 2.")
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
  b <- sqrt(b_sq)
  one_minus_s2 <- pmax(0, 1 - s^2)
  u <- b * sqrt(one_minus_s2)

  log_density <- rep.int(-Inf, length(s))
  interior_or_q2 <- one_minus_s2 > 0 | q == 2L
  if (any(interior_or_q2)) {
    power_term <- if (q == 2L) {
      rep.int(0, sum(interior_or_q2))
    } else {
      ((q - 2) / 2) * log(one_minus_s2[interior_or_q2])
    }
    log_density[interior_or_q2] <-
      vmf_log_normalizing_constant_intrinsic(q, kappa) -
      vapply(
        u[interior_or_q2],
        vmf_log_normalizing_constant_intrinsic,
        numeric(1),
        q = q - 1L
      ) +
      a * s[interior_or_q2] +
      power_term
  }

  list(
    density = exp(log_density),
    a = a,
    b = b,
    u = u,
    one_minus_s2 = one_minus_s2,
    omega = omega,
    kappa = kappa,
    q = q
  )
}

vmf_profile_derivative_table_xi <- function(omega,
                                            xi,
                                            grid_size = 4097L) {
  grid_size <- as.integer(grid_size)
  if (!is.finite(grid_size) || grid_size < 17L || grid_size %% 2L != 1L) {
    stop("The vMF derivative grid size must be an odd integer of at least 17.")
  }
  s_grid <- seq(-1, 1, length.out = grid_size)
  projected <- vmf_projected_density_canonical(s_grid, xi = xi, omega = omega)
  ratio_over_u <- profile_derivative_bessel_ratio_over_argument(
    projected$u,
    q = projected$q
  )
  C <- projected$one_minus_s2 * ratio_over_u
  integrands <- cbind(
    F = projected$density,
    M1 = s_grid * projected$density,
    MC = C * projected$density
  )
  cumulative <- profile_derivative_cumtrapz(s_grid, integrands)
  colnames(cumulative) <- colnames(integrands)

  list(
    grid = s_grid,
    cumulative = cumulative,
    a = projected$a,
    b = projected$b,
    omega = projected$omega,
    xi = as.numeric(xi),
    kappa = projected$kappa,
    mu = if (projected$kappa > 0) {
      as.numeric(xi) / projected$kappa
    } else {
      c(1, rep.int(0, projected$q))
    },
    q = projected$q,
    density = projected$density,
    C = C
  )
}

vmf_profile_and_derivative_xi <- function(omega,
                                          xi,
                                          t_values,
                                          distance_type = c("geodesic", "chordal"),
                                          grid_size = 4097L,
                                          table = NULL) {
  distance_type <- match.arg(distance_type)
  t_values <- as.numeric(t_values)
  if (any(is.na(t_values))) {
    stop("`t_values` cannot contain missing values.")
  }
  if (is.null(table)) {
    table <- vmf_profile_derivative_table_xi(
      omega = omega,
      xi = xi,
      grid_size = grid_size
    )
  }
  threshold <- if (identical(distance_type, "geodesic")) {
    cos(t_values)
  } else {
    1 - t_values^2 / 2
  }
  threshold <- pmin(pmax(threshold, -1), 1)
  lower_moments <- profile_derivative_interpolate_columns(
    table$grid,
    table$cumulative,
    threshold
  )
  totals <- table$cumulative[nrow(table$cumulative), ]
  upper_moments <- matrix(
    totals,
    nrow = nrow(lower_moments),
    ncol = length(totals),
    byrow = TRUE
  ) - lower_moments
  F_value <- upper_moments[, 1L]
  M1 <- upper_moments[, 2L]
  MC <- upper_moments[, 3L]
  A_value <- if (table$kappa == 0) 0 else A_q(table$kappa, table$q)
  derivative <- outer(M1 - table$a * MC, table$omega) +
    outer(table$kappa * MC - A_value * F_value, table$mu)

  below <- !is.finite(t_values) & t_values < 0 | t_values <= 0
  full <- if (identical(distance_type, "geodesic")) {
    t_values >= pi
  } else {
    t_values >= 2
  }
  full[is.infinite(t_values) & t_values > 0] <- TRUE
  if (any(below)) {
    F_value[below] <- 0
    derivative[below, ] <- 0
  }
  if (any(full)) {
    F_value[full] <- 1
    derivative[full, ] <- 0
  }

  list(
    F = pmin(pmax(F_value, 0), 1),
    derivative = derivative,
    M1 = M1,
    MC = MC,
    table = table
  )
}

hvmf_canonical_score_matrix <- function(sample, xi) {
  sample <- normalize_hvmf_hq_data(sample)
  xi <- as.numeric(xi)
  q <- length(xi) - 1L
  if (ncol(sample) != length(xi)) {
    stop("The HvMF sample and canonical parameter have incompatible dimensions.")
  }
  kappa_sq <- -hvmf_minkowski_inner_product(xi, xi)
  if (!is.finite(kappa_sq) || kappa_sq <= 0 || xi[[1L]] <= 0) {
    stop("HvMF `xi` must lie in the future timelike cone.")
  }
  kappa <- sqrt(kappa_sq)
  mu <- xi / kappa
  B_value <- hvmf_mean_resultant_ratio(q, kappa)
  score <- sweep(sample, 2L, B_value * mu, "-")
  score[, 1L] <- -score[, 1L]
  score
}

hvmf_canonical_information <- function(xi) {
  xi <- as.numeric(xi)
  q <- length(xi) - 1L
  kappa_sq <- -hvmf_minkowski_inner_product(xi, xi)
  if (!is.finite(kappa_sq) || kappa_sq <= 0 || xi[[1L]] <= 0) {
    stop("HvMF `xi` must lie in the future timelike cone.")
  }
  kappa <- sqrt(kappa_sq)
  mu <- xi / kappa
  J <- diag(c(-1, rep.int(1, q)))
  J_mu <- drop(J %*% mu)
  B_value <- hvmf_mean_resultant_ratio(q, kappa)
  B_prime <- B_value^2 - 1 - q * B_value / kappa
  information <- (B_value / kappa) * J +
    (B_value / kappa - B_prime) * tcrossprod(J_mu)
  (information + t(information)) / 2
}

hvmf_profile_derivative_table_xi <- function(omega,
                                             xi,
                                             upper,
                                             grid_size = 4097L) {
  omega <- as.numeric(omega)
  xi <- as.numeric(xi)
  q <- length(xi) - 1L
  grid_size <- as.integer(grid_size)
  upper <- as.numeric(upper)
  if (q < 2L || length(omega) != length(xi)) {
    stop("Deterministic HvMF derivatives require vectors in R^(q+1), q >= 2.")
  }
  if (length(upper) != 1L || !is.finite(upper) || upper < 0) {
    stop("The HvMF derivative-table upper endpoint must be finite and nonnegative.")
  }
  if (!is.finite(grid_size) || grid_size < 3L) {
    stop("The HvMF derivative grid size must be an integer of at least three.")
  }
  normalize_hvmf_hq_data(omega, q = q)
  kappa_sq <- -hvmf_minkowski_inner_product(xi, xi)
  if (!is.finite(kappa_sq) || kappa_sq <= 0 || xi[[1L]] <= 0) {
    stop("HvMF `xi` must lie in the future timelike cone.")
  }
  kappa <- sqrt(kappa_sq)
  mu <- xi / kappa
  normalize_hvmf_hq_data(mu, q = q)
  a <- hvmf_minkowski_inner_product(xi, omega)
  if (a > -kappa + 1e-10 * max(1, kappa)) {
    stop("The HvMF center and canonical parameter violate the reverse Cauchy--Schwarz inequality.")
  }
  b_sq <- profile_derivative_nonnegative_square(
    a^2 - kappa^2,
    scale = max(a^2, kappa^2),
    label = "HvMF b^2"
  )
  b <- sqrt(b_sq)
  chi <- asinh(b / kappa)
  r_grid <- if (upper == 0) c(0, 0) else seq(0, upper, length.out = grid_size)
  density <- hvmf_radial_density(
    r_grid,
    q = q,
    kappa = kappa,
    chi = chi
  )
  sinh_r <- sinh(r_grid)
  u <- b * sinh_r
  D <- sinh_r^2 * profile_derivative_bessel_ratio_over_argument(u, q = q)
  integrands <- cbind(
    F = density,
    N1 = cosh(r_grid) * density,
    ND = D * density
  )
  cumulative <- profile_derivative_cumtrapz(r_grid, integrands)
  colnames(cumulative) <- colnames(integrands)

  list(
    grid = r_grid,
    cumulative = cumulative,
    a = a,
    b = b,
    omega = omega,
    xi = xi,
    kappa = kappa,
    mu = mu,
    q = q,
    density = density,
    D = D
  )
}

hvmf_profile_and_derivative_xi <- function(omega,
                                           xi,
                                           t_values,
                                           grid_size = 4097L,
                                           table = NULL) {
  t_values <- as.numeric(t_values)
  if (any(is.na(t_values)) || any(t_values < 0)) {
    stop("HvMF `t_values` must be nonnegative and cannot contain missing values.")
  }
  finite_values <- t_values[is.finite(t_values)]
  upper <- if (length(finite_values)) max(finite_values) else 0
  if (is.null(table)) {
    table <- hvmf_profile_derivative_table_xi(
      omega = omega,
      xi = xi,
      upper = upper,
      grid_size = grid_size
    )
  }
  moments <- profile_derivative_interpolate_columns(
    table$grid,
    table$cumulative,
    pmin(t_values, max(table$grid))
  )
  F_value <- moments[, 1L]
  N1 <- moments[, 2L]
  ND <- moments[, 3L]
  B_value <- hvmf_mean_resultant_ratio(table$q, table$kappa)
  J_omega <- table$omega
  J_omega[[1L]] <- -J_omega[[1L]]
  J_mu <- table$mu
  J_mu[[1L]] <- -J_mu[[1L]]

  if (table$b == 0) {
    derivative <- outer(N1 - B_value * F_value, J_mu)
  } else {
    derivative <- outer(N1 + table$a * ND, J_omega) +
      outer(table$kappa * ND - B_value * F_value, J_mu)
  }
  at_zero <- t_values == 0
  at_infinity <- is.infinite(t_values) & t_values > 0
  if (any(at_zero)) {
    F_value[at_zero] <- 0
    derivative[at_zero, ] <- 0
  }
  if (any(at_infinity)) {
    F_value[at_infinity] <- 1
    derivative[at_infinity, ] <- 0
  }

  list(
    F = pmin(pmax(F_value, 0), 1),
    derivative = derivative,
    N1 = N1,
    ND = ND,
    table = table
  )
}

profile_derivative_stack_centers <- function(centers,
                                             thresholds,
                                             evaluator) {
  centers <- as.matrix(centers)
  thresholds <- as.matrix(thresholds)
  if (nrow(centers) != nrow(thresholds)) {
    stop("Derivative centers and threshold rows have incompatible dimensions.")
  }
  do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
    evaluator(centers[i, ], thresholds[i, ])
  }))
}

fast_multiplier_deterministic_vhat_diagnostics <- function(S_obs,
                                                            Vhat,
                                                            par0) {
  S_obs <- as.matrix(S_obs)
  Vhat <- as.matrix(Vhat)
  eigenvalues <- eigen(
    (Vhat + t(Vhat)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  rcond_value <- rcond(Vhat)
  list(
    S_obs_dim = dim(S_obs),
    Psi_aux_dim = c(0L, ncol(S_obs)),
    Vhat_dim = dim(Vhat),
    score_mean_aux = rep.int(0, ncol(S_obs)),
    score_mean_aux_norm = 0,
    Vhat_eigenvalues = as.numeric(eigenvalues),
    Vhat_rcond = as.numeric(rcond_value),
    Vhat_condition_number = if (rcond_value > 0) 1 / rcond_value else Inf,
    par0 = as.numeric(par0)
  )
}

# Retained temporarily as a private validation reference.  Production uses the
# direct joint contrast inversion defined below; in particular it never forms
# derivatives by subtracting nearly equal CDFs.

gaussian_ruben_settings_legacy <- function(control = list()) {
  scalar <- function(name, default, lower = 0, integer = FALSE) {
    value <- suppressWarnings(as.numeric(control[[name]] %||% default)[1L])
    if (!is.finite(value) || value <= lower) {
      stop(sprintf("`control$%s` must be greater than %g.", name, lower))
    }
    if (integer) as.integer(value) else value
  }
  list(
    algorithm = "joint_ruben_gamma_mixture",
    abs_tol = scalar("gaussian_quadrature_abs_tol", 1e-9),
    max_terms = scalar(
      "gaussian_quadrature_max_terms", 50000L, lower = 2, integer = TRUE
    ),
    min_terms = scalar(
      "gaussian_quadrature_min_terms", 8L, lower = 0, integer = TRUE
    ),
    eigen_rel_tol = scalar("gaussian_quadrature_eigen_rel_tol", 1e-12),
    clip_tol = scalar("gaussian_quadrature_clip_tol", 1e-9)
  )
}

gaussian_profile_variant_specification <- function(q) {
  q <- as.integer(q)
  if (length(q) != 1L || !is.finite(q) || q < 1L) {
    stop("The Gaussian dimension must be a positive integer.")
  }
  h <- matrix(1, nrow = q, ncol = 1L)
  labels <- "base"
  plus2 <- plus4 <- integer(q)
  for (i in seq_len(q)) {
    candidate <- rep.int(1, q)
    candidate[[i]] <- 3
    h <- cbind(h, candidate)
    plus2[[i]] <- ncol(h)
    labels <- c(labels, sprintf("df_%d_plus_2", i))
  }
  for (i in seq_len(q)) {
    candidate <- rep.int(1, q)
    candidate[[i]] <- 5
    h <- cbind(h, candidate)
    plus4[[i]] <- ncol(h)
    labels <- c(labels, sprintf("df_%d_plus_4", i))
  }
  cross <- matrix(NA_integer_, q, q)
  if (q >= 2L) {
    for (i in seq_len(q - 1L)) {
      for (j in seq.int(i + 1L, q)) {
        candidate <- rep.int(1, q)
        candidate[c(i, j)] <- 3
        h <- cbind(h, candidate)
        cross[i, j] <- cross[j, i] <- ncol(h)
        labels <- c(labels, sprintf("df_%d_%d_plus_2", i, j))
      }
    }
  }
  colnames(h) <- labels
  list(
    h = h,
    total_df = colSums(h),
    base = 1L,
    plus2 = plus2,
    plus4 = plus4,
    cross = cross,
    labels = labels
  )
}

gaussian_joint_ruben_coefficients <- function(lambda,
                                               delta,
                                               variants,
                                               settings) {
  lambda <- as.numeric(lambda)
  delta <- as.numeric(delta)
  h <- as.matrix(variants$h)
  q <- length(lambda)
  if (length(delta) != q || nrow(h) != q ||
      any(!is.finite(lambda)) || any(lambda <= 0) ||
      any(!is.finite(delta)) || any(delta < 0)) {
    stop("Invalid weights or noncentralities for Gaussian quadrature.")
  }

  beta <- min(lambda)
  ratio <- beta / lambda
  r <- 1 - ratio
  log_a0 <- colSums(0.5 * h * log(ratio)) - 0.5 * sum(delta)
  if (any(log_a0 < log(.Machine$double.xmin))) {
    stop(paste(
      "Gaussian quadrature coefficient initialization underflowed.",
      "The covariance/noncentrality configuration is too extreme for the",
      "positive Ruben expansion at the requested tolerance."
    ))
  }

  n_variants <- ncol(h)
  count_mean <- sum(0.5 * r / ratio + 0.5 * delta / ratio)
  count_variance <- sum(
    0.5 * r / ratio^2 +
      0.5 * delta * (1 + r) / ratio^2
  )
  max_r <- max(r)
  geometric_terms <- if (max_r > 0) {
    ceiling(log(settings$abs_tol / (10 * n_variants)) / log(max_r))
  } else {
    settings$min_terms
  }
  estimated_terms <- ceiling(max(
    128,
    settings$min_terms,
    geometric_terms,
    count_mean + 10 * sqrt(max(count_variance, 0)) +
      if (max_r > 0) 2 * max_r / (1 - max_r) else 0
  ))
  calculation_terms <- if (isTRUE(settings$force_max_terms)) {
    settings$max_terms
  } else {
    min(settings$max_terms, estimated_terms)
  }
  max_length <- calculation_terms + 1L
  component_coefficients <- vector("list", q)
  for (j in seq_len(q)) {
    alpha <- 0.5
    d_value <- 0.5 * delta[[j]] * ratio[[j]]
    component <- numeric(max_length)
    component[[1L]] <- exp(
      alpha * log(ratio[[j]]) - 0.5 * delta[[j]]
    )
    if (calculation_terms >= 1L) {
      component[[2L]] <-
        (alpha * r[[j]] + d_value) * component[[1L]]
    }
    if (calculation_terms >= 2L) {
      for (n in seq_len(calculation_terms - 1L)) {
        component[[n + 2L]] <- (
          (r[[j]] * (2 * n + alpha) + d_value) * component[[n + 1L]] -
            r[[j]]^2 * (n + alpha - 1) * component[[n]]
        ) / (n + 1)
      }
    }
    if (any(!is.finite(component)) ||
        any(component < -1000 * .Machine$double.eps)) {
      stop("Gaussian quadrature produced invalid component coefficients.")
    }
    component_coefficients[[j]] <- pmax(component, 0)
  }

  base_coefficients <- component_coefficients[[1L]]
  if (q >= 2L) {
    for (j in 2:q) {
      convolution_length <- length(base_coefficients) +
        length(component_coefficients[[j]]) - 1L
      fft_length <- 2^ceiling(log2(convolution_length))
      left_fft <- fft(c(
        base_coefficients,
        numeric(fft_length - length(base_coefficients))
      ))
      right_fft <- fft(c(
        component_coefficients[[j]],
        numeric(fft_length - length(component_coefficients[[j]]))
      ))
      base_coefficients <- Re(fft(
        left_fft * right_fft,
        inverse = TRUE
      ) / fft_length)[seq_len(max_length)]
      small_negative <- base_coefficients < 0 &
        base_coefficients >= -1000 * .Machine$double.eps
      base_coefficients[small_negative] <- 0
      if (any(!is.finite(base_coefficients)) ||
          any(base_coefficients < 0)) {
        stop("Gaussian quadrature FFT convolution produced invalid coefficients.")
      }
    }
  }

  coefficients <- matrix(0, nrow = max_length, ncol = n_variants)
  coefficients[, variants$base] <- base_coefficients
  if (n_variants > 1L) {
    filtered_once <- matrix(0, nrow = max_length, ncol = q)
    for (i in seq_len(q)) {
      filtered_once[1L, i] <- base_coefficients[[1L]]
      for (k in seq_len(calculation_terms)) {
        filtered_once[k + 1L, i] <- base_coefficients[k + 1L] +
          r[[i]] * filtered_once[k, i]
      }
      coefficients[, variants$plus2[[i]]] <-
        ratio[[i]] * filtered_once[, i]

      filtered_twice <- numeric(max_length)
      filtered_twice[[1L]] <- filtered_once[1L, i]
      for (k in seq_len(calculation_terms)) {
        filtered_twice[k + 1L] <- filtered_once[k + 1L, i] +
          r[[i]] * filtered_twice[[k]]
      }
      coefficients[, variants$plus4[[i]]] <- ratio[[i]]^2 * filtered_twice
    }
    if (q >= 2L) {
      for (i in seq_len(q - 1L)) {
        for (j in seq.int(i + 1L, q)) {
          filtered_cross <- numeric(max_length)
          filtered_cross[[1L]] <- filtered_once[1L, i]
          for (k in seq_len(calculation_terms)) {
            filtered_cross[k + 1L] <- filtered_once[k + 1L, i] +
              r[[j]] * filtered_cross[[k]]
          }
          coefficients[, variants$cross[i, j]] <-
            ratio[[i]] * ratio[[j]] * filtered_cross
        }
      }
    }
  }
  if (any(!is.finite(coefficients)) ||
      any(coefficients < -1000 * .Machine$double.eps)) {
    stop("Gaussian quadrature produced invalid joint variant coefficients.")
  }
  coefficients <- pmax(coefficients, 0)
  cumulative_mass <- apply(coefficients, 2L, cumsum)
  if (is.null(dim(cumulative_mass))) {
    cumulative_mass <- matrix(cumulative_mass, ncol = n_variants)
  }
  converged_rows <- which(
    apply(cumulative_mass >= 1 - settings$abs_tol, 1L, all)
  )
  converged_rows <- converged_rows[converged_rows >= settings$min_terms + 1L]
  if (!length(converged_rows)) {
    residual <- pmax(0, 1 - cumulative_mass[nrow(cumulative_mass), ])
    if (calculation_terms < settings$max_terms) {
      expanded_settings <- settings
      expanded_settings$force_max_terms <- TRUE
      return(gaussian_joint_ruben_coefficients(
        lambda = lambda,
        delta = delta,
        variants = variants,
        settings = expanded_settings
      ))
    }
    stop(sprintf(
      paste(
        "Gaussian quadrature did not attain tolerance %.3g within %d terms;",
        "maximum omitted coefficient mass is %.3g."
      ),
      settings$abs_tol, calculation_terms, max(residual)
    ))
  }
  terms_used <- converged_rows[[1L]] - 1L
  keep <- seq_len(terms_used + 1L)
  coefficient_mass <- cumulative_mass[terms_used + 1L, ]
  residual <- pmax(0, 1 - coefficient_mass)
  log_coefficients <- matrix(NA_real_, nrow = length(keep),
                             ncol = n_variants)
  positive <- coefficients[keep, , drop = FALSE] > 0
  log_coefficients[positive] <-
    log(coefficients[keep, , drop = FALSE][positive])
  list(
    beta = beta,
    coefficients = coefficients[keep, , drop = FALSE],
    log_coefficients = log_coefficients,
    coefficient_mass = coefficient_mass,
    residual_bound = residual,
    terms_used = terms_used,
    max_r = max_r,
    log_a0 = log_a0
  )
}

gaussian_joint_ruben_probabilities <- function(x,
                                               variants,
                                               table,
                                               settings) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("Finite squared Gaussian-ball thresholds must be nonnegative.")
  }
  n_thresholds <- length(x)
  n_variants <- ncol(variants$h)
  lower <- matrix(0, nrow = n_thresholds, ncol = n_variants,
                  dimnames = list(NULL, variants$labels))
  upper <- matrix(0, nrow = n_thresholds, ncol = n_variants,
                  dimnames = list(NULL, variants$labels))
  term_index <- 0:table$terms_used
  for (total_df in sort(unique(variants$total_df))) {
    variant_index <- which(variants$total_df == total_df)
    dfs <- total_df + 2 * term_index
    lower_basis <- vapply(dfs, function(df_value) {
      stats::pchisq(x / table$beta, df = df_value)
    }, numeric(n_thresholds))
    upper_basis <- vapply(dfs, function(df_value) {
      stats::pchisq(x / table$beta, df = df_value, lower.tail = FALSE)
    }, numeric(n_thresholds))
    if (is.null(dim(lower_basis))) {
      lower_basis <- matrix(lower_basis, nrow = n_thresholds)
      upper_basis <- matrix(upper_basis, nrow = n_thresholds)
    }
    lower[, variant_index] <- lower_basis %*%
      table$coefficients[, variant_index, drop = FALSE]
    upper[, variant_index] <- upper_basis %*%
      table$coefficients[, variant_index, drop = FALSE]
  }
  excursion <- max(c(-lower, lower - 1, -upper, upper - 1))
  if (is.finite(excursion) && excursion > settings$clip_tol) {
    stop(sprintf(
      "Gaussian quadrature probability excursion %.3g exceeds clipping tolerance.",
      excursion
    ))
  }
  list(
    lower = pmin(pmax(lower, 0), 1),
    upper = pmin(pmax(upper, 0), 1)
  )
}

gaussian_ball_profile_ruben_legacy <- function(omega,
                                             mu,
                                             Sigma,
                                             t_values,
                                             control = list(),
                                             spectral = NULL,
                                             compute_derivative = TRUE) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  t_values <- as.numeric(t_values)
  q <- length(mu)
  if (length(omega) != q || !identical(dim(Sigma), c(q, q)) ||
      any(!is.finite(c(omega, mu, Sigma))) || any(is.na(t_values))) {
    stop("Gaussian quadrature received incompatible or non-finite inputs.")
  }
  settings <- gaussian_quadrature_settings(control)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  spectral <- spectral %||% eigen(Sigma, symmetric = TRUE)
  lambda <- as.numeric(spectral$values)
  U <- as.matrix(spectral$vectors)
  eigen_floor <- settings$eigen_rel_tol * max(lambda)
  if (any(!is.finite(lambda)) || max(lambda) <= 0 ||
      min(lambda) <= eigen_floor) {
    stop(sprintf(
      paste(
        "Gaussian quadrature requires a positive-definite covariance;",
        "minimum eigenvalue %.6g is not above the relative floor %.6g."
      ), min(lambda), eigen_floor
    ))
  }

  nu <- drop(crossprod(U, mu - omega))
  delta <- nu^2 / lambda
  compute_derivative <- isTRUE(compute_derivative)
  variants <- if (compute_derivative) {
    gaussian_profile_variant_specification(q)
  } else {
    list(
      h = matrix(1, nrow = q, ncol = 1L,
                 dimnames = list(NULL, "base")),
      total_df = q,
      base = 1L,
      labels = "base"
    )
  }
  table <- gaussian_joint_ruben_coefficients(
    lambda = lambda,
    delta = delta,
    variants = variants,
    settings = settings
  )

  n_t <- length(t_values)
  F_value <- numeric(n_t)
  gradient_mu <- matrix(0, nrow = n_t, ncol = q)
  gradient_sigma <- matrix(0, nrow = n_t,
                           ncol = q * (q + 1L) / 2L)
  hessian_nu <- array(0, dim = c(q, q, n_t))
  finite_positive <- is.finite(t_values) & t_values > 0
  infinite_positive <- is.infinite(t_values) & t_values > 0
  F_value[infinite_positive] <- 1
  probability_details <- NULL

  if (any(finite_positive)) {
    t_positive <- t_values[finite_positive]
    x <- t_positive^2
    if (any(!is.finite(x))) {
      stop("A finite Gaussian-ball threshold overflowed when squared.")
    }
    probability_details <- gaussian_joint_ruben_probabilities(
      x = x,
      variants = variants,
      table = table,
      settings = settings
    )
    lower <- probability_details$lower
    upper <- probability_details$upper
    use_upper <- lower[, variants$base] > 0.5
    F_positive <- lower[, variants$base]
    if (any(use_upper)) {
      F_positive[use_upper] <- 1 - upper[use_upper, variants$base]
    }
    F_value[finite_positive] <- F_positive

    if (compute_derivative) for (row in seq_along(x)) {
      values <- if (use_upper[[row]]) upper[row, ] else lower[row, ]
      tail_sign <- if (use_upper[[row]]) -1 else 1
      delta_plus2 <- numeric(q)
      second_diagonal <- numeric(q)
      for (i in seq_len(q)) {
        delta_plus2[[i]] <- tail_sign * (
          values[[variants$plus2[[i]]]] - values[[variants$base]]
        )
        second_diagonal[[i]] <- tail_sign * (
          values[[variants$plus4[[i]]]] -
            2 * values[[variants$plus2[[i]]]] +
            values[[variants$base]]
        )
      }
      g <- (nu / lambda) * delta_plus2
      H <- matrix(0, q, q)
      diag(H) <- delta_plus2 / lambda +
        (nu^2 / lambda^2) * second_diagonal
      if (q >= 2L) {
        for (i in seq_len(q - 1L)) {
          for (j in seq.int(i + 1L, q)) {
            cross_difference <- tail_sign * (
              values[[variants$cross[i, j]]] -
                values[[variants$plus2[[i]]]] -
                values[[variants$plus2[[j]]]] +
                values[[variants$base]]
            )
            H[i, j] <- H[j, i] <-
              (nu[[i]] * nu[[j]] / (lambda[[i]] * lambda[[j]])) *
              cross_difference
          }
        }
      }
      grad_mu_row <- drop(U %*% g)
      K <- 0.5 * U %*% H %*% t(U)
      grad_sigma_row <- fast_multiplier_sym_score_to_vech(K)
      result_row <- which(finite_positive)[[row]]
      gradient_mu[result_row, ] <- grad_mu_row
      gradient_sigma[result_row, ] <- grad_sigma_row
      hessian_nu[, , result_row] <- H
    }
  }

  derivative <- cbind(gradient_mu, gradient_sigma)
  list(
    F = F_value,
    derivative = derivative,
    gradient_mu = gradient_mu,
    gradient_vech_sigma = gradient_sigma,
    hessian_nu = hessian_nu,
    nu = nu,
    lambda = lambda,
    eigenvectors = U,
    delta = delta,
    diagnostics = list(
      algorithm = settings$algorithm,
      abs_tol = settings$abs_tol,
      max_terms = settings$max_terms,
      terms_used = table$terms_used,
      residual_bound = max(table$residual_bound),
      residual_bound_by_variant = table$residual_bound,
      beta = table$beta,
      max_r = table$max_r,
      condition_number = max(lambda) / min(lambda),
      eigen_rel_tol = settings$eigen_rel_tol,
      clip_tol = settings$clip_tol,
      variants = variants$labels
    ),
    coefficient_table = table
  )
}

# Direct joint Gil--Pelaez inversion for the Gaussian profile and the
# contrasts required by its derivatives.  These definitions deliberately
# supersede the Ruben implementation above: forming differences of separately
# computed CDFs is both slower and less accurate in the tails.

gaussian_quadrature_settings <- function(control = list()) {
  positive_scalar <- function(name, default) {
    value <- suppressWarnings(as.numeric(control[[name]] %||% default)[1L])
    if (!is.finite(value) || value <= 0) {
      stop(sprintf("`control$%s` must be positive and finite.", name))
    }
    value
  }
  positive_integer <- function(name, default) {
    as.integer(ceiling(positive_scalar(name, default)))
  }
  list(
    algorithm = "joint_gil_pelaez_contrasts_nested_infinite_rule",
    abs_tol = positive_scalar("gaussian_quadrature_abs_tol", 1e-5),
    initial_upper = positive_scalar("gaussian_quadrature_initial_upper", 8),
    max_upper = positive_scalar("gaussian_quadrature_max_upper", 8192),
    max_intervals = positive_integer(
      "gaussian_quadrature_max_intervals",
      control$gaussian_quadrature_max_terms %||% 1000000L
    ),
    max_terms = positive_integer(
      "gaussian_quadrature_max_intervals",
      control$gaussian_quadrature_max_terms %||% 1000000L
    ),
    tail_consecutive = positive_integer(
      "gaussian_quadrature_tail_consecutive", 3L
    ),
    eigen_rel_tol = positive_scalar("gaussian_quadrature_eigen_rel_tol", 1e-12),
    clip_tol = positive_scalar("gaussian_quadrature_clip_tol", 1e-9)
  )
}

gaussian_vector_gk15_interval <- function(fun, left, right) {
  xgk <- c(
    0.9914553711208126, 0.9491079123427585, 0.8648644233597691,
    0.7415311855993944, 0.5860872354676911, 0.4058451513773972,
    0.2077849550078985, 0
  )
  wgk <- c(
    0.02293532201052922, 0.06309209262997855, 0.1047900103222502,
    0.1406532597155259, 0.1690047266392679, 0.1903505780647854,
    0.2044329400752989, 0.2094821410847278
  )
  wg <- c(0.1294849661688697, 0.2797053914892767,
          0.3818300505051189, 0.4179591836734694)
  midpoint <- 0.5 * (left + right)
  half_width <- 0.5 * (right - left)
  centre <- fun(midpoint)
  kronrod <- wgk[[8L]] * centre
  gauss <- wg[[4L]] * centre
  for (k in seq_len(7L)) {
    pair_sum <- fun(midpoint - half_width * xgk[[k]]) +
      fun(midpoint + half_width * xgk[[k]])
    kronrod <- kronrod + wgk[[k]] * pair_sum
    if (k %in% c(2L, 4L, 6L)) {
      gauss <- gauss + wg[[k / 2L]] * pair_sum
    }
  }
  value <- half_width * kronrod
  list(value = value, error = abs(value - half_width * gauss),
       left = left, right = right, evaluations = 15L)
}

gaussian_vector_gk15 <- function(fun, left, right, target,
                                 max_intervals, max_panel_width) {
  breaks <- seq(left, right, by = max_panel_width)
  if (!length(breaks) || tail(breaks, 1L) < right) breaks <- c(breaks, right)
  if (length(breaks) == 1L) breaks <- c(left, right)
  intervals <- lapply(seq_len(length(breaks) - 1L), function(k) {
    gaussian_vector_gk15_interval(fun, breaks[[k]], breaks[[k + 1L]])
  })
  total_value <- Reduce(`+`, lapply(intervals, `[[`, "value"))
  total_error <- Reduce(`+`, lapply(intervals, `[[`, "error"))
  evaluations <- 15L * length(intervals)
  while (max(total_error / target) > 1) {
    if (length(intervals) >= max_intervals) {
      stop(sprintf(
        "Joint Gaussian quadrature exhausted %d intervals (scaled error %.3g).",
        max_intervals, max(total_error / target)
      ))
    }
    scores <- vapply(intervals, function(z) max(z$error / target), numeric(1L))
    index <- which.max(scores)
    old <- intervals[[index]]
    midpoint <- 0.5 * (old$left + old$right)
    first <- gaussian_vector_gk15_interval(fun, old$left, midpoint)
    second <- gaussian_vector_gk15_interval(fun, midpoint, old$right)
    total_value <- total_value - old$value + first$value + second$value
    total_error <- total_error - old$error + first$error + second$error
    intervals[[index]] <- first
    intervals[[length(intervals) + 1L]] <- second
    evaluations <- evaluations + 30L
  }
  list(value = total_value, error = total_error,
       intervals = length(intervals), evaluations = evaluations)
}

gaussian_joint_contrast_integrand <- function(u, y, a, delta, pairs,
                                               compute_derivative = TRUE,
                                               include_base = TRUE) {
  q <- length(a)
  first_moment <- sum(a * (1 + delta))
  if (u == 0) {
    base <- (y - first_moment) / pi
    if (!isTRUE(compute_derivative)) return(base)
    A <- matrix(rep(-2 * a / pi, each = length(y)), nrow = length(y))
    B <- matrix(0, length(y), q)
    C <- matrix(0, length(y), ncol(pairs))
    return(c(if (isTRUE(include_base)) base else numeric(), A, B, C))
  }
  denominator <- 1 - 2i * u * a
  log_phi <- sum(-0.5 * log(denominator) +
                   1i * u * a * delta / denominator)
  z <- exp(-1i * u * y + log_phi)
  d <- 1 / denominator - 1
  base <- -Im(z) / (pi * u)
  if (!isTRUE(compute_derivative)) return(base)
  A <- -Im(outer(z, d)) / (pi * u)
  B <- -Im(outer(z, d^2)) / (pi * u)
  C <- if (ncol(pairs)) {
    factors <- d[pairs[1L, ]] * d[pairs[2L, ]]
    -Im(outer(z, factors)) / (pi * u)
  } else matrix(0, length(y), 0L)
  c(if (isTRUE(include_base)) base else numeric(), A, B, C)
}

gaussian_joint_contrast_trapezoid <- function(upper, step, y, a, delta,
                                               pairs,
                                               compute_derivative = TRUE,
                                               include_base = TRUE) {
  number_steps <- as.integer(ceiling(upper / step))
  step <- upper / number_steps
  u <- seq.int(0, number_steps) * step
  weights <- rep(step, number_steps + 1L)
  weights[c(1L, number_steps + 1L)] <- step / 2
  positive <- seq.int(2L, number_steps + 1L)
  up <- u[positive]
  denominator <- outer(up, a, function(left, right) 1 - 2i * left * right)
  log_phi <- rowSums(-0.5 * log(denominator) +
                       1i * up * rep(a * delta, each = length(up)) /
                       denominator)
  z <- exp(-1i * outer(up, y) + log_phi)
  weighted_z <- z * (weights[positive] / (pi * up))
  first_moment <- sum(a * (1 + delta))
  base <- weights[[1L]] * (y - first_moment) / pi - colSums(Im(weighted_z))
  if (!isTRUE(compute_derivative)) {
    return(list(value = base, nodes = number_steps + 1L,
                step = step, upper = upper))
  }
  d <- 1 / denominator - 1
  q <- length(a); nt <- length(y)
  A <- B <- matrix(0, nt, q)
  for (i in seq_len(q)) {
    A[, i] <- -weights[[1L]] * 2 * a[[i]] / pi -
      colSums(Im(weighted_z * d[, i]))
    B[, i] <- -colSums(Im(weighted_z * d[, i]^2))
  }
  C <- matrix(0, nt, ncol(pairs))
  if (ncol(pairs)) for (k in seq_len(ncol(pairs))) {
    C[, k] <- -colSums(Im(
      weighted_z * d[, pairs[1L, k]] * d[, pairs[2L, k]]
    ))
  }
  list(value = c(if (isTRUE(include_base)) base else numeric(), A, B, C),
       nodes = number_steps + 1L,
       step = step, upper = upper)
}

gaussian_joint_contrast_quadinf <- function(fun, target) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Joint Gaussian infinite-interval quadrature requires `pracma`.")
  }
  rule <- getFromNamespace(".quadinf_pre", "pracma")()
  transformed <- function(node) {
    z <- (node + 1) / 2
    u <- z / (1 - z)
    2 / (1 - node)^2 * fun(u)
  }
  nodes <- rule$nodes[[1L]]; weights <- rule$weights[[1L]]
  value <- weights[[7L]] * transformed(nodes[[7L]])
  for (j in seq_len(6L)) {
    value <- value + weights[[j]] *
      (transformed(nodes[[j]]) + transformed(-nodes[[j]]))
  }
  h <- 0.5; value <- h * value
  evaluations <- 13L
  error <- rep(Inf, length(value))
  iteration <- 1L
  for (k in 2:7) {
    nodes <- rule$nodes[[k]]; weights <- rule$weights[[k]]
    increment <- numeric(length(value))
    for (j in seq_along(weights)) {
      increment <- increment + weights[[j]] *
        (transformed(nodes[[j]]) + transformed(-nodes[[j]]))
    }
    h <- h / 2
    refined <- h * increment + value / 2
    error <- abs(refined - value)
    value <- refined
    evaluations <- evaluations + 2L * length(weights)
    iteration <- k
    if (max(error / target) <= 1) break
  }
  if (max(error / target) > 1) {
    stop(sprintf(
      "Joint Gaussian infinite quadrature did not attain tolerance (scaled estimate %.3g).",
      max(error / target)
    ))
  }
  list(value = value, error = error, evaluations = evaluations,
       iterations = iteration)
}

gaussian_joint_contrast_trapezoid_refine <- function(y, a, delta, pairs,
                                                      target, settings,
                                                      compute_derivative) {
  frequency <- max(c(abs(y), sum(a * (1 + delta)), 1))
  step <- min(0.05, pi / (8 * frequency))
  upper <- settings$initial_upper
  evaluations <- 0L
  repeat {
    if (ceiling(upper / step) + 1L > settings$max_intervals) {
      stop(sprintf(
        "Joint Gaussian quadrature requires more than %d shared nodes.",
        settings$max_intervals
      ))
    }
    fine <- gaussian_joint_contrast_trapezoid(
      upper, step, y, a, delta, pairs, compute_derivative,
      include_base = FALSE
    )
    coarse <- gaussian_joint_contrast_trapezoid(
      upper, 2 * step, y, a, delta, pairs, compute_derivative,
      include_base = FALSE
    )
    half_domain <- gaussian_joint_contrast_trapezoid(
      upper / 2, step, y, a, delta, pairs, compute_derivative,
      include_base = FALSE
    )
    discretization <- abs(fine$value - coarse$value) / 3
    truncation <- 2 * abs(fine$value - half_domain$value)
    error <- discretization + truncation
    evaluations <- evaluations + fine$nodes + coarse$nodes + half_domain$nodes
    if (max(error / target) <= 1) break
    if (max(discretization / target) > 0.5) step <- step / 2
    if (max(truncation / target) > 0.5) upper <- 2 * upper
    if (upper > settings$max_upper) {
      stop(sprintf(
        "Joint Gaussian quadrature did not attain tolerance by u=%g (scaled estimate %.3g).",
        settings$max_upper, max(error / target)
      ))
    }
  }
  list(value = fine$value, error = error, evaluations = evaluations,
       iterations = NA_integer_, nodes = fine$nodes, upper = upper,
       step = fine$step,
       estimator = "Richardson_plus_successive_truncation_difference")
}

gaussian_univariate_ball_profile <- function(omega, mu, variance, t_values) {
  sd_value <- sqrt(variance)
  finite_positive <- is.finite(t_values) & t_values > 0
  F_value <- gradient_mu <- gradient_sigma <- numeric(length(t_values))
  F_value[is.infinite(t_values) & t_values > 0] <- 1
  if (any(finite_positive)) {
    upper <- (omega + t_values[finite_positive] - mu) / sd_value
    lower <- (omega - t_values[finite_positive] - mu) / sd_value
    phi_upper <- stats::dnorm(upper)
    phi_lower <- stats::dnorm(lower)
    F_value[finite_positive] <- stats::pnorm(upper) - stats::pnorm(lower)
    gradient_mu[finite_positive] <- (phi_lower - phi_upper) / sd_value
    gradient_sigma[finite_positive] <-
      (lower * phi_lower - upper * phi_upper) / (2 * variance)
  }
  list(F = F_value, derivative = cbind(gradient_mu, gradient_sigma),
       gradient_mu = matrix(gradient_mu, ncol = 1L),
       gradient_vech_sigma = matrix(gradient_sigma, ncol = 1L),
       hessian_nu = array(2 * gradient_sigma, c(1L, 1L, length(t_values))))
}

gaussian_ball_profile_quadrature <- function(omega, mu, Sigma, t_values,
                                             control = list(), spectral = NULL,
                                             compute_derivative = TRUE) {
  omega <- as.numeric(omega); mu <- as.numeric(mu)
  Sigma <- 0.5 * (as.matrix(Sigma) + t(as.matrix(Sigma)))
  t_values <- as.numeric(t_values); q <- length(mu)
  if (length(omega) != q || !identical(dim(Sigma), c(q, q)) ||
      any(!is.finite(c(omega, mu, Sigma))) || any(is.na(t_values))) {
    stop("Invalid input to Gaussian quadrature.")
  }
  settings <- gaussian_quadrature_settings(control)
  spectral <- spectral %||% eigen(Sigma, symmetric = TRUE)
  lambda <- as.numeric(spectral$values); U <- as.matrix(spectral$vectors)
  if (max(lambda) <= 0 || min(lambda) <= settings$eigen_rel_tol * max(lambda)) {
    stop("Gaussian quadrature requires a positive-definite covariance above the relative eigenvalue floor.")
  }
  nu <- drop(crossprod(U, mu - omega)); delta <- nu^2 / lambda
  if (q == 1L) {
    answer <- gaussian_univariate_ball_profile(omega, mu, Sigma[[1L]], t_values)
    answer$nu <- nu; answer$lambda <- lambda; answer$eigenvectors <- U
    answer$delta <- delta
    answer$diagnostics <- list(
      algorithm = "univariate_closed_form", abs_tol = settings$abs_tol,
      residual_error_estimate = 0, propagated_error_estimate = 0,
      terms_used = 0L,
      evaluations = 0L, upper_limit = Inf,
      condition_number = 1, eigen_rel_tol = settings$eigen_rel_tol,
      clip_tol = settings$clip_tol
    )
    answer$coefficient_table <- NULL
    return(answer)
  }

  finite_positive <- is.finite(t_values) & t_values > 0
  n_t <- length(t_values)
  F_value <- numeric(n_t); F_value[is.infinite(t_values) & t_values > 0] <- 1
  gradient_mu <- matrix(0, n_t, q)
  gradient_sigma <- matrix(0, n_t, q * (q + 1L) / 2L)
  hessian_nu <- array(0, c(q, q, n_t))
  diagnostics <- list(algorithm = settings$algorithm, abs_tol = settings$abs_tol,
                      residual_error_estimate = 0,
                      propagated_error_estimate = 0,
                      terms_used = 0L, evaluations = 0L, upper_limit = 0)
  if (any(finite_positive)) {
    radii <- t_values[finite_positive]
    density_bound <- (2 * pi)^(-q / 2) / sqrt(prod(lambda))
    probability_bound <- pi^(q / 2) / gamma(q / 2 + 1) *
      radii^q * density_bound
    inverse_norm <- 1 / min(lambda)
    displacement <- sqrt(sum((mu - omega)^2)) + radii
    score_mu_bound <- inverse_norm * displacement
    score_sigma_bound <- inverse_norm^2 *
      (displacement^2 + sqrt(sum(Sigma^2)))
    output_bound <- probability_bound *
      pmax(1, score_mu_bound, score_sigma_bound)
    negligible_local <- output_bound <= settings$abs_tol / 4
    negligible <- which(finite_positive)[negligible_local]
    finite_positive[negligible] <- FALSE
    diagnostics$small_ball_analytic_bound <- if (length(negligible)) {
      max(output_bound[negligible_local])
    } else 0
  }
  if (any(finite_positive)) {
    scale <- max(lambda); a <- lambda / scale
    y <- t_values[finite_positive]^2 / scale
    if (!exists("mvnormal_quadform_cdf", mode = "function")) {
      stop("Gaussian quadrature requires the shared MVN quadratic-form CDF engine.")
    }
    F_value[finite_positive] <- mvnormal_quadform_cdf(
      q = t_values[finite_positive]^2,
      lambda = lambda,
      h = rep.int(1, q),
      delta = delta,
      control = control
    )
    diagnostics$profile_cdf_algorithm <- "mvnormal_quadform_cdf"
    if (!isTRUE(compute_derivative)) {
      diagnostics$condition_number <- max(lambda) / min(lambda)
      diagnostics$eigen_rel_tol <- settings$eigen_rel_tol
      diagnostics$clip_tol <- settings$clip_tol
      return(list(
        F = F_value,
        derivative = cbind(gradient_mu, gradient_sigma),
        gradient_mu = gradient_mu,
        gradient_vech_sigma = gradient_sigma,
        hessian_nu = hessian_nu,
        nu = nu, lambda = lambda, eigenvectors = U, delta = delta,
        diagnostics = diagnostics, coefficient_table = NULL
      ))
    }
    pairs <- if (q >= 2L) utils::combn(q, 2L) else matrix(integer(), 2L, 0L)
    nt <- length(y); np <- ncol(pairs)
    if (isTRUE(compute_derivative)) {
      amplification_A <- pmax(abs(nu / lambda), 1 / lambda)
      amplification_B <- pmax(1, nu^2 / lambda^2)
      amplification_C <- if (np) abs(nu[pairs[1L, ]] * nu[pairs[2L, ]] /
                                        (lambda[pairs[1L, ]] * lambda[pairs[2L, ]])) else numeric()
      amplification <- c(rep(amplification_A, each = nt),
                         rep(amplification_B, each = nt),
                         rep(pmax(1, amplification_C), each = nt))
    } else {
      amplification <- rep(1, nt)
    }
    target <- settings$abs_tol / pmax(1, amplification)
    integrand <- function(u) gaussian_joint_contrast_integrand(
      u, y, a, delta, pairs, compute_derivative, include_base = FALSE
    )
    infinite_attempt <- tryCatch(
      gaussian_joint_contrast_quadinf(integrand, target),
      error = function(error) error
    )
    if (inherits(infinite_attempt, "error")) {
      integral <- gaussian_joint_contrast_trapezoid_refine(
        y, a, delta, pairs, target, settings, compute_derivative
      )
      diagnostics$algorithm <-
        "joint_gil_pelaez_contrasts_shared_trapezoid"
      diagnostics$primary_attempt_error <- conditionMessage(infinite_attempt)
    } else {
      integral <- infinite_attempt
      integral$nodes <- integral$evaluations
      integral$upper <- Inf
      integral$step <- NA_real_
      integral$estimator <- "successive_nested_infinite_rule_difference"
    }
    total <- integral$value
    total_error <- integral$error
    diagnostics$evaluations <- integral$evaluations
    diagnostics$terms_used <- integral$evaluations
    diagnostics$residual_error_estimate <- max(total_error)
    diagnostics$propagated_error_estimate <-
      max(total_error * amplification)
    diagnostics$upper_limit <- integral$upper
    diagnostics$step <- integral$step
    diagnostics$iterations <- integral$iterations
    diagnostics$error_estimator <- integral$estimator
    diagnostics$infinite_interval_transform <- if (is.infinite(integral$upper)) {
      "u=z/(1-z)"
    } else NA_character_
    offset <- 0L
    take <- function(number) {
      index <- offset + seq_len(number); offset <<- offset + number; total[index]
    }
    A <- matrix(take(nt * q), nt, q)
    B <- matrix(take(nt * q), nt, q)
    C <- if (np) matrix(take(nt * np), nt, np) else matrix(0, nt, 0L)
    rows <- which(finite_positive)
    if (isTRUE(compute_derivative)) for (r in seq_len(nt)) {
      g <- (nu / lambda) * A[r, ]
      H <- matrix(0, q, q)
      diag(H) <- A[r, ] / lambda + (nu^2 / lambda^2) * B[r, ]
      if (np) for (k in seq_len(np)) {
        i <- pairs[1L, k]; j <- pairs[2L, k]
        H[i, j] <- H[j, i] <-
          nu[[i]] * nu[[j]] / (lambda[[i]] * lambda[[j]]) * C[r, k]
      }
      gradient_mu[rows[[r]], ] <- drop(U %*% g)
      K <- 0.5 * U %*% H %*% t(U)
      gradient_sigma[rows[[r]], ] <- fast_multiplier_sym_score_to_vech(K)
      hessian_nu[, , rows[[r]]] <- H
    }
  }
  diagnostics$condition_number <- max(lambda) / min(lambda)
  diagnostics$eigen_rel_tol <- settings$eigen_rel_tol
  diagnostics$clip_tol <- settings$clip_tol
  list(F = F_value, derivative = cbind(gradient_mu, gradient_sigma),
       gradient_mu = gradient_mu, gradient_vech_sigma = gradient_sigma,
       hessian_nu = hessian_nu, nu = nu, lambda = lambda,
       eigenvectors = U, delta = delta, diagnostics = diagnostics,
       coefficient_table = NULL)
}

gaussian_fast_quadrature_tables <- function(data_centers,
                                            ks_centers,
                                            ks_prep,
                                            cvm_prep,
                                            observed_distance_matrix = NULL,
                                            evaluator) {
  data_centers <- as.matrix(data_centers)
  ks_centers <- if (is.null(ks_centers)) NULL else as.matrix(ks_centers)
  diagnostics <- new.env(parent = emptyenv())
  diagnostics$calls <- 0L
  diagnostics$max_terms <- 0L
  diagnostics$max_residual <- 0
  diagnostics$max_condition_number <- 0
  diagnostics$max_propagated_error <- 0
  diagnostics$max_upper_limit <- 0
  diagnostics$max_evaluations <- 0L
  diagnostics$algorithms <- character()

  evaluate_derivative <- function(center, thresholds) {
    result <- evaluator(center, thresholds)
    diagnostics$calls <- diagnostics$calls + 1L
    diagnostics$algorithms <- unique(c(
      diagnostics$algorithms, result$diagnostics$algorithm
    ))
    diagnostics$max_terms <- max(
      diagnostics$max_terms, result$diagnostics$terms_used
    )
    diagnostics$max_residual <- max(
      diagnostics$max_residual,
      result$diagnostics$residual_error_estimate
    )
    diagnostics$max_condition_number <- max(
      diagnostics$max_condition_number,
      result$diagnostics$condition_number
    )
    diagnostics$max_propagated_error <- max(
      diagnostics$max_propagated_error,
      result$diagnostics$propagated_error_estimate %||% 0
    )
    diagnostics$max_upper_limit <- max(
      diagnostics$max_upper_limit,
      result$diagnostics$upper_limit %||% 0
    )
    diagnostics$max_evaluations <- max(
      diagnostics$max_evaluations,
      result$diagnostics$evaluations %||% 0L
    )
    result$derivative
  }

  D_ks <- if (is.null(ks_prep)) {
    NULL
  } else if (identical(
    ks_prep$ks_grid_mode %||% "fixed",
    "sample_points_unique_distances"
  )) {
    list(
      mode = "sample_points_unique_distances",
      derivative_sorted = profile_derivative_stack_centers(
        centers = ks_centers,
        thresholds = ks_prep$sorted_distance_matrix,
        evaluator = evaluate_derivative
      )
    )
  } else {
    do.call(rbind, lapply(seq_len(nrow(ks_centers)), function(i) {
      evaluate_derivative(ks_centers[i, ], ks_prep$t_grid)
    }))
  }

  D_cvm <- if (is.null(cvm_prep)) {
    NULL
  } else if (isTRUE(cvm_prep$shared_with_ks) &&
             is.list(D_ks) && !is.null(D_ks$derivative_sorted)) {
    list(
      mode = "sample_points_unique_distances_sorted_rows",
      derivative_sorted = D_ks$derivative_sorted,
      shared_with_ks = TRUE
    )
  } else if (isTRUE(cvm_prep$light) &&
             !is.null(cvm_prep$sorted_distance_matrix)) {
    list(
      mode = "sample_points_unique_distances_sorted_rows",
      derivative_sorted = profile_derivative_stack_centers(
        centers = data_centers,
        thresholds = cvm_prep$sorted_distance_matrix,
        evaluator = evaluate_derivative
      )
    )
  } else {
    if (is.null(observed_distance_matrix)) {
      stop("Dense Gaussian CvM quadrature requires observed distances.")
    }
    profile_derivative_stack_centers(
      centers = data_centers,
      thresholds = observed_distance_matrix,
      evaluator = evaluate_derivative
    )
  }
  list(
    D_ks = D_ks,
    D_cvm = D_cvm,
    diagnostics = list(
      center_evaluations = diagnostics$calls,
      max_terms_used = diagnostics$max_terms,
      max_residual_error_estimate = diagnostics$max_residual,
      max_condition_number = diagnostics$max_condition_number,
      max_propagated_error_estimate = diagnostics$max_propagated_error,
      max_upper_limit = diagnostics$max_upper_limit,
      max_evaluations = diagnostics$max_evaluations,
      algorithms_effective = diagnostics$algorithms
    )
  )
}
