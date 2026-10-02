
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
