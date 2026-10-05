# Helpers used by research scripts but outside the public package.

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
joint_probability_vmf_s1_chordal_exact <- function(omega1,
                                                   t1,
                                                   omega2,
                                                   t2,
                                                   mu,
                                                   kappa,
                                                   cdf_object = NULL,
                                                   cdf_grid_size = 16385) {
  if (t1 <= 0 || t2 <= 0) return(0)
  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }
  seg1 <- s1_event_segments_chordal(omega1, t1)
  seg2 <- s1_event_segments_chordal(omega2, t2)
  overlap <- s1_intersect_segments(seg1, seg2)
  vmf_s1_segments_probability(overlap, cdf_object)
}
invert_distance_profile_vmf_s1_chordal <- function(omega,
                                                   mu,
                                                   kappa,
                                                   u_values,
                                                   cdf_object = NULL,
                                                   cdf_grid_size = 16385,
                                                   tol = 1e-8) {
  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }

  u_values <- as.numeric(u_values)
  vapply(u_values, function(u) {
    if (u <= 0) return(0)
    if (u >= 1) return(2)
    root_fun <- function(t) {
      theoretical_distance_profile_vmf_s1_chordal(
        omega = omega,
        mu = mu,
        kappa = kappa,
        t_values = t,
        cdf_object = cdf_object,
        cdf_grid_size = cdf_grid_size
      ) - u
    }
    uniroot(root_fun, interval = c(0, 2), tol = tol)$root
  }, numeric(1))
}
cov_vmf_s1_simple_exact <- function(omega_grid,
                                    t_grid,
                                    mu,
                                    kappa,
                                    cdf_object = NULL,
                                    cdf_grid_size = 16385) {
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  if (ncol(omega_grid) != 2) {
    stop("`omega_grid` must have exactly 2 columns for exact S^1 covariance.")
  }
  if (any(t_grid < 0) || any(t_grid > 2)) {
    stop("`t_grid` must lie in [0, 2] for chordal distance on S^1.")
  }
  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }

  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  F_matrix <- t(vapply(seq_len(n_omega), function(i) {
    theoretical_distance_profile_vmf_s1_chordal(
      omega = omega_grid[i, ],
      mu = mu,
      kappa = kappa,
      t_values = t_grid,
      cdf_object = cdf_object,
      cdf_grid_size = cdf_grid_size
    )
  }, numeric(n_t)))
  F_vec <- as.vector(F_matrix)

  segment_list <- vector("list", n_total)
  for (t_idx in seq_len(n_t)) {
    for (omega_idx in seq_len(n_omega)) {
      flat_idx <- omega_idx + (t_idx - 1) * n_omega
      segment_list[[flat_idx]] <- s1_event_segments_chordal(
        omega = omega_grid[omega_idx, ],
        t = t_grid[t_idx]
      )
    }
  }

  cov_matrix <- matrix(0, nrow = n_total, ncol = n_total)
  for (i in seq_len(n_total)) {
    for (j in i:n_total) {
      joint_prob <- vmf_s1_segments_probability(
        s1_intersect_segments(segment_list[[i]], segment_list[[j]]),
        cdf_object = cdf_object
      )
      cov_value <- joint_prob - F_vec[i] * F_vec[j]
      cov_matrix[i, j] <- cov_value
      cov_matrix[j, i] <- cov_value
    }
  }

  cov_matrix
}
joint_probability_vmf_s2_simple_integral <- function(omega1,
                                                     t1,
                                                     omega2,
                                                     t2,
                                                     mu,
                                                     kappa,
                                                     distance_type = "chordal",
                                                     rel.tol_outer = 1e-7,
                                                     abs.tol_outer = 1e-9,
                                                     rel.tol_inner = 1e-8,
                                                     abs.tol_inner = 1e-10,
                                                     subdivisions_outer = 200L,
                                                     subdivisions_inner = 200L,
                                                     tol = 1e-10) {
  omega1 <- as.numeric(omega1)
  omega2 <- as.numeric(omega2)
  mu <- as.numeric(mu)
  if (length(omega1) != 3 || length(omega2) != 3 || length(mu) != 3) {
    stop("`omega1`, `omega2`, and `mu` must all have length 3 for S^2 exact probabilities.")
  }
  omega1 <- omega1 / sqrt(sum(omega1^2))
  omega2 <- omega2 / sqrt(sum(omega2^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))

  a1 <- sphere_distance_to_dot_threshold(t1, distance_type)
  a2 <- sphere_distance_to_dot_threshold(t2, distance_type)
  if (a1 >= 1 || a2 >= 1) return(0)
  if (a1 <= -1 && a2 <= -1) return(1)

  rho <- sum(omega1 * omega2)
  rho <- pmin(pmax(rho, -1), 1)
  m1 <- sum(mu * omega1)
  m1 <- pmin(pmax(m1, -1), 1)

  density_u <- function(u) {
    vmf_s2_projected_density(u, omega = omega1, mu = mu, kappa = kappa)
  }

  if (rho >= 1 - tol) {
    lower <- max(a1, a2, -1)
    if (lower >= 1) return(0)
    return(integrate(
      f = density_u,
      lower = lower,
      upper = 1,
      rel.tol = rel.tol_outer,
      abs.tol = abs.tol_outer,
      subdivisions = subdivisions_outer
    )$value)
  }

  if (rho <= -1 + tol) {
    lower <- max(a1, -1)
    upper <- min(1, -a2)
    if (upper <= lower + tol) return(0)
    return(integrate(
      f = density_u,
      lower = lower,
      upper = upper,
      rel.tol = rel.tol_outer,
      abs.tol = abs.tol_outer,
      subdivisions = subdivisions_outer
    )$value)
  }

  denom_rho <- sqrt(max(0, 1 - rho^2))
  denom_mu <- sqrt(max(0, 1 - m1^2))
  cos_alpha <- if (denom_mu <= tol) {
    1
  } else {
    (sum(mu * omega2) - m1 * rho) / (denom_mu * denom_rho)
  }
  alpha <- acos(pmin(pmax(cos_alpha, -1), 1))

  lower <- max(a1, -1)
  if (lower >= 1) return(0)

  integrand <- function(u_vec) {
    vapply(u_vec, function(u) {
      radial <- sqrt(max(0, 1 - u^2))
      if (radial <= tol) {
        inner_prob <- as.numeric(rho * u >= a2 - tol)
      } else {
        b <- (a2 - rho * u) / (denom_rho * radial)
        lambda <- kappa * radial * denom_mu
        inner_prob <- vmf_s1_cap_probability_integral(
          lambda = lambda,
          alpha = alpha,
          b = b,
          rel.tol = rel.tol_inner,
          abs.tol = abs.tol_inner,
          subdivisions = subdivisions_inner,
          tol = tol
        )
      }
      density_u(u) * inner_prob
    }, numeric(1))
  }

  integrate(
    f = integrand,
    lower = lower,
    upper = 1,
    rel.tol = rel.tol_outer,
    abs.tol = abs.tol_outer,
    subdivisions = subdivisions_outer
  )$value
}
cov_vmf_s2_simple_integral <- function(omega_grid,
                                       t_grid,
                                       mu,
                                       kappa,
                                       distance_type = "chordal",
                                       rel.tol_outer = 1e-7,
                                       abs.tol_outer = 1e-9,
                                       rel.tol_inner = 1e-8,
                                       abs.tol_inner = 1e-10,
                                       subdivisions_outer = 200L,
                                       subdivisions_inner = 200L,
                                       tol = 1e-10) {
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  mu <- as.numeric(mu)
  if (ncol(omega_grid) != 3 || length(mu) != 3) {
    stop("`omega_grid` must have 3 columns and `mu` must have length 3 for S^2 exact covariance.")
  }

  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t

  F_matrix <- t(vapply(seq_len(n_omega), function(i) {
    distance_profile_vmf_s2_integral(
      omega = omega_grid[i, ],
      mu = mu,
      kappa = kappa,
      t_values = t_grid,
      distance_type = distance_type,
      rel.tol = rel.tol_outer,
      abs.tol = abs.tol_outer,
      subdivisions = subdivisions_outer
    )
  }, numeric(n_t)))
  F_vec <- as.vector(F_matrix)

  omega_idx_vec <- ((0:(n_total - 1)) %% n_omega) + 1
  t_idx_vec <- ((0:(n_total - 1)) %/% n_omega) + 1
  cov_matrix <- matrix(0, nrow = n_total, ncol = n_total)

  for (i in seq_len(n_total)) {
    omega1 <- omega_grid[omega_idx_vec[i], ]
    t1 <- t_grid[t_idx_vec[i]]
    for (j in i:n_total) {
      joint_prob <- joint_probability_vmf_s2_simple_integral(
        omega1 = omega1,
        t1 = t1,
        omega2 = omega_grid[omega_idx_vec[j], ],
        t2 = t_grid[t_idx_vec[j]],
        mu = mu,
        kappa = kappa,
        distance_type = distance_type,
        rel.tol_outer = rel.tol_outer,
        abs.tol_outer = abs.tol_outer,
        rel.tol_inner = rel.tol_inner,
        abs.tol_inner = abs.tol_inner,
        subdivisions_outer = subdivisions_outer,
        subdivisions_inner = subdivisions_inner,
        tol = tol
      )
      cov_value <- joint_prob - F_vec[i] * F_vec[j]
      cov_matrix[i, j] <- cov_value
      cov_matrix[j, i] <- cov_value
    }
  }

  cov_matrix
}
compute_conditional_expectation_vmf <- function(omega, t, mu, kappa,
                                               distance_type, mc_samples) {
  if (t <= 0) return(rep(0, length(mu)))
  dots <- mc_samples %*% omega
  if (distance_type == "chordal") {
    dists <- sqrt(2 * (1 - dots))
  } else {
    dots <- check_dot_products(dots)
    dists <- acos(dots)
  }
  in_ball <- dists <= t
  if (sum(in_ball) == 0) return(rep(0, length(mu)))
  filtered_samples <- mc_samples[in_ball, , drop = FALSE]
  return(colMeans(filtered_samples))
}
compute_precomp_vmf <- function(mc_samples, omega_grid, t_grid, mu, kappa, distance_type, A_q_kappa) {
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  q <- length(mu) - 1
  n_mc <- nrow(mc_samples)

  # Dot products and distances from all MC samples to all omegas
  dots_all <- mc_samples %*% t(omega_grid)
  if (distance_type != "chordal") dots_all <- check_dot_products(dots_all)
  dists_all <- if (distance_type == "chordal") sqrt(2 * (1 - dots_all)) else acos(dots_all)

  # F2_matrix: marginal probabilities for each omega2 and t
  F2_matrix <- t(apply(omega_grid, 1, function(omega2) {
    theoretical_distance_profile_vmf(omega2, mu, kappa, t_grid, distance_type)
  }))

  # Build E2_mat (n_total x q+1) using a simpler and explicit rbind-by-t approach
  in_ball_list <- vector("list", n_t)
  E2_rows_by_t <- lapply(seq_len(n_t), function(k) {
    in_ball_k <- dists_all <= t_grid[k]
    in_ball_list[[k]] <<- in_ball_k
    counts_k <- colSums(in_ball_k)
    if (all(counts_k == 0)) return(matrix(0, nrow = n_omega, ncol = q + 1))
    sums_mat <- t(mc_samples) %*% (in_ball_k * 1) # (q+1) x n_omega
    counts_k_safe <- counts_k
    counts_k_safe[counts_k_safe == 0] <- 1
    means_mat <- sweep(sums_mat, 2, counts_k_safe, "/")
    if (any(counts_k == 0)) means_mat[, counts_k == 0] <- 0
    return(t(means_mat))
  })
  E2_mat <- do.call(rbind, E2_rows_by_t)
  E2_array <- array(NA, dim = c(n_omega, n_t, q + 1))
  for (k in seq_len(n_t)) {
    E2_array[, k, ] <- E2_rows_by_t[[k]]
  }

  F2_vec <- as.vector(F2_matrix)
  m2_mat <- sweep(E2_mat, 2, A_q_kappa * mu, "-")

  return(list(
    dists_all = dists_all,
    F2_matrix = F2_matrix,
    E2_array = E2_array,
    E2_mat = E2_mat,
    F2_vec = F2_vec,
    m2_mat = m2_mat,
    in_ball_list = in_ball_list
  ))
}
small_circle_benchmark_nu_optimization <- function(n = 200L,
                                                   n_rep = 40L,
                                                   true_mu = c(0, 0, 1),
                                                   true_kappa = 8,
                                                   true_nu = 0.15,
                                                   seed = 1,
                                                   control_logistic = list(),
                                                   control_box = list()) {
  n <- as.integer(n)
  n_rep <- as.integer(n_rep)
  if (length(n) != 1L || !is.finite(n) || n < 5L) {
    stop("`n` must be an integer >= 5.")
  }
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1L) {
    stop("`n_rep` must be a strictly positive integer.")
  }

  true_params <- small_circle_validate_parameters(mu = true_mu, kappa = true_kappa, nu = true_nu)
  set.seed(as.integer(seed))

  evaluate_fit <- function(method_name,
                           x,
                           theta_start,
                           base_control) {
    control_method <- utils::modifyList(
      list(
        small_circle_nu_optim = method_name,
        small_circle_mle_start_theta = theta_start
      ),
      base_control
    )

    start_time <- proc.time()[[3L]]
    fit <- try(small_circle_mle_s2_weighted(x = x, control = control_method), silent = TRUE)
    elapsed <- proc.time()[[3L]] - start_time

    if (inherits(fit, "try-error")) {
      return(data.frame(
        method = method_name,
        elapsed_sec = elapsed,
        converged = FALSE,
        loglik = NA_real_,
        abs_kappa_error = NA_real_,
        abs_nu_error = NA_real_,
        mu_angle_error = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    dot_mu <- pmin(pmax(sum(fit$mu * true_params$mu), -1), 1)
    data.frame(
      method = method_name,
      elapsed_sec = elapsed,
      converged = is.list(fit$opt) && isTRUE(fit$opt$convergence == 0L),
      loglik = fit$loglik,
      abs_kappa_error = abs(fit$kappa - true_params$kappa),
      abs_nu_error = abs(fit$nu - true_params$nu),
      mu_angle_error = acos(dot_mu),
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(seq_len(n_rep), function(rep_id) {
    x <- r_sph_small_circle(
      n = n,
      mu = true_params$mu,
      kappa = true_params$kappa,
      nu = true_params$nu,
      check = FALSE
    )
    theta_start <- small_circle_start_theta_s2(x = x)

    row_logistic <- evaluate_fit(
      method_name = "logistic",
      x = x,
      theta_start = theta_start,
      base_control = control_logistic
    )
    row_box <- evaluate_fit(
      method_name = "box",
      x = x,
      theta_start = theta_start,
      base_control = control_box
    )

    cbind(rep = rep_id, rbind(row_logistic, row_box), stringsAsFactors = FALSE)
  })

  details <- do.call(rbind, rows)
  split_details <- split(details, details$method)

  summary <- do.call(rbind, lapply(names(split_details), function(method_name) {
    d <- split_details[[method_name]]
    converged_idx <- which(d$converged)
    d_conv <- if (length(converged_idx) == 0L) d[FALSE, , drop = FALSE] else d[converged_idx, , drop = FALSE]

    data.frame(
      method = method_name,
      n_rep = nrow(d),
      convergence_rate = mean(d$converged),
      mean_elapsed_sec = mean(d$elapsed_sec, na.rm = TRUE),
      median_elapsed_sec = stats::median(d$elapsed_sec, na.rm = TRUE),
      mean_loglik = if (nrow(d_conv) > 0L) mean(d_conv$loglik, na.rm = TRUE) else NA_real_,
      mean_abs_kappa_error = if (nrow(d_conv) > 0L) mean(d_conv$abs_kappa_error, na.rm = TRUE) else NA_real_,
      mean_abs_nu_error = if (nrow(d_conv) > 0L) mean(d_conv$abs_nu_error, na.rm = TRUE) else NA_real_,
      mean_mu_angle_error = if (nrow(d_conv) > 0L) mean(d_conv$mu_angle_error, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(
    setup = list(
      n = n,
      n_rep = n_rep,
      true_mu = true_params$mu,
      true_kappa = true_params$kappa,
      true_nu = true_params$nu,
      seed = seed
    ),
    summary = summary,
    details = details
  )
}
small_circle_compare_profile_methods <- function(mu,
                                                 kappa,
                                                 nu,
                                                 omega_list,
                                                 t_grid,
                                                 distance_type = c("geodesic", "chordal"),
                                                 l_max = 150L,
                                                 quad_n = 1000L,
                                                 tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  comparison_rows <- lapply(seq_along(omega_list), function(i) {
    omega <- omega_list[[i]]
    legendre <- distance_profile_small_circle(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      kappa = kappa,
      nu = nu,
      distance_type = distance_type,
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    )
    integral <- distance_profile_small_circle(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      kappa = kappa,
      nu = nu,
      distance_type = distance_type,
      method = "integral",
      quad_n = quad_n,
      tol = tol
    )

    data.frame(
      omega_id = i,
      max_abs_diff = max(abs(legendre - integral)),
      mean_abs_diff = mean(abs(legendre - integral)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, comparison_rows)
}
rotational_gauss_hermite <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  key <- as.character(n)
  if (exists(key, envir = rotational_gauss_hermite_cache, inherits = FALSE)) {
    return(get(key, envir = rotational_gauss_hermite_cache, inherits = FALSE))
  }

  off_diag <- sqrt(seq_len(n - 1L) / 2)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off_diag
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off_diag
  eig <- eigen(jacobi, symmetric = TRUE)
  order_idx <- order(eig$values)

  out <- list(
    nodes = eig$values[order_idx],
    weights = sqrt(pi) * (eig$vectors[1L, order_idx]^2)
  )
  assign(key, out, envir = rotational_gauss_hermite_cache)
  out
}
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
validate_covariance_matrix <- function(cov_matrix,
                                       symmetry_tol = 1e-8,
                                       psd_tol = 1e-8,
                                       stop_on_failure = TRUE) {
  cov_matrix <- as.matrix(cov_matrix)
  if (nrow(cov_matrix) != ncol(cov_matrix)) {
    stop("Covariance matrix must be square.")
  }
  if (!all(is.finite(cov_matrix))) {
    stop("Covariance matrix contains non-finite entries.")
  }

  symmetry_gap <- max(abs(cov_matrix - t(cov_matrix)))
  eigenvalues <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  min_eigenvalue <- min(eigenvalues)
  max_eigenvalue <- max(eigenvalues)
  negative_eigenvalues <- sum(eigenvalues < -psd_tol)
  is_valid <- symmetry_gap <= symmetry_tol && min_eigenvalue >= -psd_tol

  diagnostics <- list(
    valid = is_valid,
    symmetry_gap = symmetry_gap,
    min_eigenvalue = min_eigenvalue,
    max_eigenvalue = max_eigenvalue,
    negative_eigenvalues = negative_eigenvalues,
    symmetry_tol = symmetry_tol,
    psd_tol = psd_tol
  )

  if (!is_valid && isTRUE(stop_on_failure)) {
    reasons <- character(0)
    if (symmetry_gap > symmetry_tol) {
      reasons <- c(
        reasons,
        sprintf("symmetry gap %.3e exceeds tolerance %.3e", symmetry_gap, symmetry_tol)
      )
    }
    if (min_eigenvalue < -psd_tol) {
      reasons <- c(
        reasons,
        sprintf("minimum eigenvalue %.3e is below tolerance %.3e", min_eigenvalue, -psd_tol)
      )
    }
    stop(
      sprintf(
        "Covariance matrix rejected for simulation: %s.",
        paste(reasons, collapse = "; ")
      )
    )
  }

  diagnostics
}
safe_ks_test <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) {
    warning('safe_ks_test: inputs are not numeric; returning NA p.value')
    return(list(statistic = NA_real_, p.value = NA_real_, success = FALSE, message = 'non-numeric inputs'))
  }
  if (any(is.na(x)) || any(is.na(y))) {
    warning('safe_ks_test: inputs contain NA; returning NA p.value')
    return(list(statistic = NA_real_, p.value = NA_real_, success = FALSE, message = 'missing values'))
  }
  res <- tryCatch({
    ks <- ks.test(x, y)
    list(statistic = ks$statistic, p.value = ks$p.value, success = TRUE)
  }, error = function(e) {
    warning(sprintf('safe_ks_test: ks.test failed: %s', e$message))
    list(statistic = NA_real_, p.value = NA_real_, success = FALSE, message = e$message)
  })
  return(res)
}
analyze_single_trajectory <- function(mu_true, kappa_true, sample_sizes, trajectory_id) {
  q <- length(mu_true) - 1
  xi_true <- kappa_true * mu_true

  results <- data.frame()

  for (n in sample_sizes) {
    # Draw an independent sample for each n (no cumulative sample growth)
    current_sample <- r_vMF(n, mu_true, kappa_true)
    xi_hat <- compute_mle_xi(current_sample)
    quantity_1 <- sqrt(n) * (xi_hat - xi_true)

    score_sum <- rep(0, q + 1)
    for (i in 1:n) {
      score_sum <- score_sum + psi_xi(current_sample[i, ], xi_true, q)
    }

    dot_psi_matrix <- dot_psi_xi(xi_true, q)

    if (det(dot_psi_matrix) != 0) {
      dot_psi_inv <- solve(dot_psi_matrix)
      quantity_2 <- -dot_psi_inv %*% (score_sum / sqrt(n))
    } else {
      stop(sprintf("Singular matrix for trajectory %d at n = %d", trajectory_id, n))
    }

    difference_norm <- norm(quantity_1 - as.vector(quantity_2), type = "2")

    # Print progress messages at key sample sizes
    if (n %in% c(10000, 25000, 50000, 75000)) {
      cat("  >> Trajectory", trajectory_id, "reached n =", n, "- Difference norm:",
          round(difference_norm, 6), "\n")
    }

    results <- rbind(results, data.frame(
      trajectory = trajectory_id,
      n = n,
      difference_norm = difference_norm
    ))
  }

  return(results)
}

hvmf_mean_projection_cdf <- function(y, q, kappa) {
  q <- hvmf_validate_dimension(q)
  kappa <- hvmf_validate_kappa(kappa)
  y <- as.numeric(y)

  output <- numeric(length(y))
  output[is.infinite(y) & y > 0] <- 1
  finite_valid <- is.finite(y) & y > 1
  if (any(finite_valid)) {
    output[finite_valid] <- vapply(y[finite_valid], function(upper_bound) {
      integral <- stats::integrate(
        hvmf_mean_projection_density,
        lower = 1,
        upper = upper_bound,
        q = q,
        kappa = kappa,
        rel.tol = 1e-8,
        abs.tol = 1e-10,
        stop.on.error = FALSE
      )
      if (!identical(integral$message, "OK") || !is.finite(integral$value)) {
        stop(sprintf("Could not integrate the HvMF projection density up to y = %.8g.", upper_bound))
      }
      integral$value
    }, numeric(1))
  }
  pmin(pmax(output, 0), 1)
}

small_circle_monotone_clip_dot <- function(dot_values, values) {
  if (length(values) <= 1L) {
    out <- pmin(pmax(values, 0), 1)
    out[dot_values >= 1] <- 0
    out[dot_values <= -1] <- 1
    return(out)
  }

  order_idx <- order(dot_values, decreasing = TRUE)
  sorted_values <- pmin(pmax(values[order_idx], 0), 1)
  sorted_values <- cummax(sorted_values)
  sorted_dots <- dot_values[order_idx]
  sorted_values[sorted_dots >= 1] <- 0
  sorted_values[sorted_dots <= -1] <- 1
  out <- numeric(length(sorted_values))
  out[order_idx] <- sorted_values
  out
}

rotational_legendre_coefficients <- function(density_h,
                                             Lmax,
                                             quad_n = 1000L,
                                             tol = 1e-10) {
  if (!is.function(density_h)) {
    stop("`density_h` must be a function.")
  }

  Lmax <- as.integer(Lmax)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(Lmax) != 1L || !is.finite(Lmax) || Lmax < 0L) {
    stop("`Lmax` must be a nonnegative integer.")
  }

  quad <- rotational_gauss_legendre(quad_n)
  legendre_matrix <- rotational_legendre_matrix(quad$nodes, l_max = Lmax)
  h_values <- as.numeric(density_h(quad$nodes))
  if (length(h_values) != length(quad$nodes) || any(!is.finite(h_values)) || any(h_values < 0)) {
    stop("`density_h` must return finite nonnegative values of the same length as its input.")
  }

  raw_moments <- as.numeric(crossprod(legendre_matrix, quad$weights * h_values))
  ell <- 0:Lmax
  coeffs <- ((2 * ell + 1) / 2) * raw_moments
  a0_error <- abs(raw_moments[[1L]] / 2 - 1)
  coeffs[[1L]] <- 1

  if (a0_error > tol) {
    stop(sprintf("Rotational Legendre coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error
  )
}

vmf_s1_cap_probability_integral <- function(lambda,
                                            alpha,
                                            b,
                                            rel.tol = 1e-8,
                                            abs.tol = 1e-10,
                                            subdivisions = 200L,
                                            tol = 1e-10) {
  lambda <- as.numeric(lambda)
  alpha <- as.numeric(alpha)
  b <- as.numeric(b)

  if (b <= -1 + tol) return(1)
  if (b >= 1 - tol) return(0)

  delta <- acos(pmin(pmax(b, -1), 1))
  if (delta >= pi - tol) return(1)
  if (lambda <= tol) {
    return(delta / pi)
  }

  i0_scaled <- besselI(lambda, nu = 0, expon.scaled = TRUE)
  integrand <- function(phi) {
    exp(lambda * (cos(phi - alpha) - 1)) / (2 * pi * i0_scaled)
  }
  integrate(
    f = integrand,
    lower = -delta,
    upper = delta,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions
  )$value
}

s1_intersect_segments <- function(seg1, seg2, tol = 1e-12) {
  if (length(seg1) == 0 || length(seg2) == 0 || nrow(seg1) == 0 || nrow(seg2) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }

  intersections <- list()
  idx <- 0
  for (i in seq_len(nrow(seg1))) {
    for (j in seq_len(nrow(seg2))) {
      lower <- max(seg1[i, 1], seg2[j, 1])
      upper <- min(seg1[i, 2], seg2[j, 2])
      if (upper - lower > tol) {
        idx <- idx + 1
        intersections[[idx]] <- c(lower, upper)
      }
    }
  }

  if (length(intersections) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }

  do.call(rbind, intersections)
}
