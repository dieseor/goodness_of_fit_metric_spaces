# Mathematical and numerical helpers: numerics.

check_dot_products <- function(dot_products, tolerance = 1e-10) {
  if (any(dot_products < -1 - tolerance) || any(dot_products > 1 + tolerance)) {
    bad_values <- dot_products[dot_products < -1 - tolerance | dot_products > 1 + tolerance]
    stop(paste0("Invalid dot products detected: range [", 
                round(min(bad_values), 10), ", ", round(max(bad_values), 10), 
                "]. Check that inputs are unit vectors."))
  }
  return(dot_products)
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
minkowski_inner_product <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (length(x) != 3L || length(y) != 3L) {
    stop("`x` and `y` must both have length 3.")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("`x` and `y` must be finite.")
  }

  -x[[1]] * y[[1]] + sum(x[-1L] * y[-1L])
}
generate_canonical_lattice <- function(n, dim = 3) {
  if (dim == 3) {
    golden_ratio <- (1 + sqrt(5)) / 2
    i <- seq(0, n - 1)
    theta <- 2 * pi * i / golden_ratio
    phi <- acos(1 - 2 * (i + 0.5) / n)
    x <- cos(theta) * sin(phi)
    y <- sin(theta) * sin(phi)
    z <- cos(phi)
    mat <- cbind(x, y, z)
    return(mat)
  } else if (dim > 1) {
    mat <- matrix(rnorm(n * dim), nrow = n, ncol = dim)
    mat <- t(apply(mat, 1, function(v) v / sqrt(sum(v^2))))
    return(mat)
  } else {
    stop("Dimension must be >= 2")
  }
}
wrap_angle_2pi <- function(theta) {
  wrapped <- theta %% (2 * pi)
  wrapped[wrapped < 0] <- wrapped[wrapped < 0] + 2 * pi
  wrapped
}
generate_circle_grid <- function(n_angles, theta0 = 0) {
  if (!is.numeric(n_angles) || length(n_angles) != 1 || n_angles < 2) {
    stop("`n_angles` must be a single integer >= 2.")
  }
  theta <- wrap_angle_2pi(theta0 + 2 * pi * (0:(n_angles - 1)) / n_angles)
  theta <- sort(theta)
  data.frame(
    theta = theta,
    x = cos(theta),
    y = sin(theta),
    label = sprintf("%.2f", theta),
    stringsAsFactors = FALSE
  )
}
circle_angle_from_point <- function(omega, tol = 1e-8) {
  omega <- as.numeric(omega)
  if (length(omega) != 2) {
    stop("`omega` must have length 2 for S^1 computations.")
  }
  omega_norm <- sqrt(sum(omega^2))
  if (abs(omega_norm - 1) > tol) {
    stop("`omega` must have unit norm for S^1 computations.")
  }
  wrap_angle_2pi(atan2(omega[2], omega[1]))
}
circle_angles_from_matrix <- function(omega_grid) {
  omega_grid <- as.matrix(omega_grid)
  if (ncol(omega_grid) != 2) {
    stop("`omega_grid` must have exactly 2 columns for S^1 computations.")
  }
  apply(omega_grid, 1, circle_angle_from_point)
}
s1_interval_to_segments <- function(start, end, tol = 1e-12) {
  width <- end - start
  if (width <= tol) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }
  if (width >= 2 * pi - tol) {
    return(matrix(c(0, 2 * pi), ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("start", "end"))))
  }

  start_wrapped <- wrap_angle_2pi(start)
  end_wrapped <- wrap_angle_2pi(end)
  if (start_wrapped < end_wrapped) {
    return(matrix(c(start_wrapped, end_wrapped), ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("start", "end"))))
  }

  rbind(
    c(0, end_wrapped),
    c(start_wrapped, 2 * pi)
  )
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
compute_mle_xi <- function(data) {
  # Require movMF package for MLE
  if (!requireNamespace('movMF', quietly = TRUE)) {
    stop('movMF package is required to compute vMF MLE. Please install movMF.')
  }
  # movMF expects rows as observations; it returns theta matrix with estimated parameters
  fit <- movMF::movMF(data, k = 1)
  xi_hat <- as.vector(fit$theta[1, ])
  return(xi_hat)
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
d_proj_jp <- function(t, mu, kappa, psi, log = FALSE) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  t <- as.numeric(t)
  out <- rep(if (log) -Inf else 0, length(t))
  valid <- is.finite(t) & t >= -1 & t <= 1
  if (!any(valid)) {
    return(out)
  }

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # The axis density inherits the same S^2 regularization: near psi = 0 we
    # evaluate the exact vMF projected density instead of the JP one.
    return(jp_vmf_axis_density(
      t = t,
      q = params$q,
      kappa = params$kappa,
      log = log
    ))
  }

  log_density <- log(sphere_surface_area(params$q - 1L)) +
    c_jp_sphere(q = params$q, alpha = params$alpha, beta = params$beta, log = TRUE) +
    params$beta * log1p(params$alpha * t[valid]) +
    jp_log_one_minus_u2_term(t[valid], params$q / 2 - 1)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}
p_proj_jp <- function(s, omega, mu, kappa, psi, log = FALSE) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = length(params$mu))
  if (length(omega) != length(params$mu)) {
    stop("`omega` and `mu` must have the same length.")
  }

  s <- as.numeric(s)
  out <- rep(if (log) -Inf else 0, length(s))
  valid <- is.finite(s) & s >= -1 & s <= 1
  if (!any(valid)) {
    return(out)
  }

  rho <- pmin(pmax(sum(omega * params$mu), -1), 1)
  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # Same regularized vMF limit for general one-dimensional projections.
    return(jp_vmf_projected_density(
      s = s,
      rho = rho,
      ambient_dim = params$ambient_dim,
      kappa = params$kappa,
      log = log
    ))
  }

  reciprocal_table <- get_jp_delta_reciprocal_table(
    q = params$q - 1L,
    beta = params$beta,
    delta_max = jp_max_delta_from_rho(alpha = params$alpha, rho = rho)
  )
  num0 <- 1 + params$alpha * rho * s[valid]
  delta <- abs(params$alpha) * sqrt(pmax(0, 1 - rho^2)) * sqrt(pmax(0, 1 - s[valid]^2)) / num0
  reciprocal_values <- jp_eval_delta_reciprocal(delta, reciprocal_table)
  log_density <- c_jp_sphere(q = params$q, alpha = params$alpha, beta = params$beta, log = TRUE) +
    jp_log_one_minus_u2_term(s[valid], params$q / 2 - 1) +
    params$beta * log(num0) +
    log(reciprocal_values)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}
r_sph_jp <- function(n,
                     mu,
                     kappa,
                     psi,
                     grid_size = 8193L,
                     check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # The sampler follows the same regularization policy near psi = 0 on S^2,
    # drawing from the exact vMF limit instead of inverting the JP CDF.
    return(rotasym::r_vMF(n = n, mu = params$mu, kappa = params$kappa))
  }
  if (abs(params$alpha) <= 1e-15) {
    return(jp_uniform_sphere(n = n, ambient_dim = params$ambient_dim))
  }

  cdf_table <- build_jp_axis_cdf_table(
    mu = params$mu,
    kappa = params$kappa,
    psi = params$psi,
    grid_size = grid_size
  )
  inverse_table <- jp_prepare_inverse_cdf_table(u = cdf_table$u, cdf = cdf_table$cdf)
  inverse_interpolator <- jp_prepare_inverse_cdf_interpolator(inverse_table)
  u_samples <- stats::runif(n)
  t_samples <- inverse_interpolator$fn(
    pmin(pmax(u_samples, inverse_interpolator$x_min), inverse_interpolator$x_max)
  )
  t_samples <- pmin(pmax(t_samples, -1), 1)

  xi <- jp_uniform_sphere(n = n, ambient_dim = params$q)
  b_mu <- jp_orthonormal_complement(params$mu)
  tangent_part <- t(b_mu %*% t(xi))
  radial_scales <- sqrt(pmax(0, 1 - t_samples^2))
  x <- tcrossprod(t_samples, params$mu) + sweep(tangent_part, 1, radial_scales, "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Jones-Pewsey sampler returned non-unit vectors.")
    }
  }

  x
}
projection_profile_matrix_integral <- function(omega_grid,
                                               t_grid,
                                               mu,
                                               z_nodes,
                                               weighted_density,
                                               distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(omega_grid))
  t_grid <- as.numeric(t_grid)
  z_nodes <- as.numeric(z_nodes)
  weighted_density <- as.numeric(weighted_density)

  if (length(z_nodes) != length(weighted_density)) {
    stop("`z_nodes` and `weighted_density` must have the same length.")
  }

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
  r_values <- as.numeric(omega_grid %*% mu)

  for (i in seq_len(n_omega)) {
    kernel <- small_circle_projection_kernel(
      c_thresholds = thresholds,
      z_nodes = z_nodes,
      r = r_values[[i]]
    )
    out[i, active] <- as.numeric(kernel %*% weighted_density)
  }

  pmin(pmax(out, 0), 1)
}
projection_sample_profile_matrix_integral <- function(X,
                                                      mu,
                                                      z_nodes,
                                                      weighted_density) {
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(X))
  z_nodes <- as.numeric(z_nodes)
  weighted_density <- as.numeric(weighted_density)

  if (length(z_nodes) != length(weighted_density)) {
    stop("`z_nodes` and `weighted_density` must have the same length.")
  }

  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  r_values <- as.numeric(X %*% mu)
  out <- matrix(0, nrow = nrow(X), ncol = nrow(X))

  for (i in seq_len(nrow(X))) {
    kernel <- small_circle_projection_kernel(
      c_thresholds = dot_products[i, ],
      z_nodes = z_nodes,
      r = r_values[[i]]
    )
    out[i, ] <- as.numeric(kernel %*% weighted_density)
  }

  pmin(pmax(out, 0), 1)
}
debug_memory_log <- function(control = list(),
                             label,
                             objects = list(),
                             force = FALSE) {
  enabled <- isTRUE(force) || isTRUE(control$cvm_memory_debug %||% FALSE)
  if (!enabled) {
    return(invisible(NULL))
  }

  gc_info <- tryCatch(gc(), error = function(e) NULL)
  gc_summary <- if (is.null(gc_info)) {
    "gc=unavailable"
  } else {
    used_col <- grep("^used", colnames(gc_info), value = TRUE)[1L]
    trigger_col <- grep("^gc trigger", colnames(gc_info), value = TRUE)[1L]
    paste(
      sprintf(
        "%s used=%sMb gc_trigger=%sMb",
        rownames(gc_info),
        format(round(gc_info[, used_col], 1), trim = TRUE),
        format(round(gc_info[, trigger_col], 1), trim = TRUE)
      ),
      collapse = " | "
    )
  }

  message(sprintf("[CvM debug] %s | %s", label, gc_summary))
  if (length(objects) > 0L) {
    for (object_name in names(objects)) {
      obj <- objects[[object_name]]
      dims <- dim(obj)
      dim_text <- if (is.null(dims)) {
        sprintf("length=%d", length(obj))
      } else {
        sprintf("dim=%s", paste(dims, collapse = "x"))
      }
      size_mb <- as.numeric(utils::object.size(obj)) / 1024^2
      message(sprintf(
        "[CvM debug]   %s: %s, size=%.2f Mb, class=%s",
        object_name,
        dim_text,
        size_mb,
        paste(class(obj), collapse = "/")
      ))
    }
  }

  invisible(NULL)
}
compute_weighted_sample_profile_block <- function(order_matrix_block,
                                                  rank_linear_index_block,
                                                  normalized_weights,
                                                  n_total) {
  order_matrix_block <- as.matrix(order_matrix_block)
  rank_linear_index_block <- as.matrix(rank_linear_index_block)
  block_rows <- nrow(order_matrix_block)
  n <- ncol(order_matrix_block)

  if (length(normalized_weights) != n || n_total != n) {
    stop("Incompatible dimensions in `compute_weighted_sample_profile_block()`.")
  }
  if (any(!is.finite(normalized_weights)) || any(normalized_weights < 0)) {
    stop("`normalized_weights` must be finite and nonnegative.")
  }

  total_weight <- sum(normalized_weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("`normalized_weights` must have strictly positive finite sum.")
  }

  ordered_weights_matrix <- matrix(
    normalized_weights[order_matrix_block],
    nrow = block_rows,
    ncol = n
  )
  cumulative_weights_matrix <- ordered_weights_matrix

  if (n >= 2L) {
    for (j in 2:n) {
      cumulative_weights_matrix[, j] <- cumulative_weights_matrix[, j] +
        cumulative_weights_matrix[, j - 1L]
    }
  }

  global_linear_index <- as.integer(rank_linear_index_block)
  global_row_index <- ((global_linear_index - 1L) %% n_total) + 1L
  global_rank_index <- ((global_linear_index - global_row_index) %/% n_total) + 1L
  local_row_index <- matrix(rep.int(seq_len(block_rows), n), nrow = block_rows, ncol = n)
  local_linear_index <- local_row_index + (matrix(global_rank_index, nrow = block_rows, ncol = n) - 1L) * block_rows

  out <- matrix(cumulative_weights_matrix[local_linear_index] / total_weight, nrow = block_rows, ncol = n)
  if (any(!is.finite(out)) || any(out < -1e-12) || any(out > 1 + 1e-12)) {
    stop("`compute_weighted_sample_profile_block()` produced values outside [0, 1].")
  }
  out
}
beta_mixture2_validate_parameters <- function(mu,
                                                         weight1,
                                                         alpha1,
                                                         beta1,
                                                         alpha2,
                                                         beta2) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Rotational beta-mixture utilities currently support only S^2.")
  }

  weight1 <- as.numeric(weight1)
  alpha1 <- as.numeric(alpha1)
  beta1 <- as.numeric(beta1)
  alpha2 <- as.numeric(alpha2)
  beta2 <- as.numeric(beta2)

  if (length(weight1) != 1L || !is.finite(weight1) || weight1 <= 0 || weight1 >= 1) {
    stop("`weight1` must be a finite scalar in (0, 1).")
  }
  positive_shapes <- c(alpha1, beta1, alpha2, beta2)
  if (any(!is.finite(positive_shapes)) || any(positive_shapes <= 0)) {
    stop("Beta-mixture shape parameters must be finite and strictly positive.")
  }

  list(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2,
    ambient_dim = 3L
  )
}
beta_mixture2_canonicalize_theta <- function(theta) {
  params <- beta_mixture2_validate_parameters(
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )

  mean1 <- params$alpha1 / (params$alpha1 + params$beta1)
  mean2 <- params$alpha2 / (params$alpha2 + params$beta2)
  if (mean1 <= mean2) {
    return(params)
  }

  list(
    mu = params$mu,
    weight1 = 1 - params$weight1,
    alpha1 = params$alpha2,
    beta1 = params$beta2,
    alpha2 = params$alpha1,
    beta2 = params$beta1,
    ambient_dim = params$ambient_dim
  )
}
beta_mixture2_density_y <- function(y,
                                               weight1,
                                               alpha1,
                                               beta1,
                                               alpha2,
                                               beta2,
                                               log = FALSE,
                                               eps = 1e-12) {
  y <- rotational_clamp_unit_interval(as.numeric(y), eps = eps)
  density1 <- stats::dbeta(y, shape1 = alpha1, shape2 = beta1, log = TRUE)
  density2 <- stats::dbeta(y, shape1 = alpha2, shape2 = beta2, log = TRUE)
  log_density <- rotational_logsumexp2(
    log(weight1) + density1,
    log1p(-weight1) + density2
  )
  if (log) log_density else exp(log_density)
}
beta_mixture2_cdf_y <- function(y,
                                           weight1,
                                           alpha1,
                                           beta1,
                                           alpha2,
                                           beta2) {
  y <- as.numeric(y)
  out <- numeric(length(y))
  out[y <= 0] <- 0
  out[y >= 1] <- 1
  active <- which(y > 0 & y < 1)
  if (length(active) == 0L) {
    return(out)
  }

  out[active] <- weight1 * stats::pbeta(y[active], alpha1, beta1) +
    (1 - weight1) * stats::pbeta(y[active], alpha2, beta2)
  pmin(pmax(out, 0), 1)
}
beta_mixture2_density_h <- function(z,
                                               weight1,
                                               alpha1,
                                               beta1,
                                               alpha2,
                                               beta2,
                                               eps = 1e-12) {
  z <- as.numeric(z)
  y <- rotational_clamp_unit_interval((z + 1) / 2, eps = eps)
  out <- numeric(length(z))
  valid <- z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- beta_mixture2_density_y(
    y = y[valid],
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2,
    log = FALSE,
    eps = eps
  )
  out
}
beta_mixture2_density_gz <- function(z,
                                                weight1,
                                                alpha1,
                                                beta1,
                                                alpha2,
                                                beta2) {
  0.5 * beta_mixture2_density_h(
    z = z,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
}
d_sph_beta_mixture2_s2 <- function(x,
                                              mu,
                                              weight1,
                                              alpha1,
                                              beta1,
                                              alpha2,
                                              beta2,
                                              log = FALSE,
                                              eps = 1e-12) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  log_density <- -log(4 * pi) + beta_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    alpha1 = params$alpha1,
    beta1 = params$beta1,
    alpha2 = params$alpha2,
    beta2 = params$beta2,
    log = TRUE,
    eps = eps
  )
  if (log) log_density else exp(log_density)
}
beta_mixture2_weighted_loglik_s2 <- function(mu,
                                                        weight1,
                                                        alpha1,
                                                        beta1,
                                                        alpha2,
                                                        beta2,
                                                        x,
                                                        prob_weights = NULL,
                                                        eps = 1e-12) {
  params <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  sum(prob_weights * beta_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    alpha1 = params$alpha1,
    beta1 = params$beta1,
    alpha2 = params$alpha2,
    beta2 = params$beta2,
    log = TRUE,
    eps = eps
  ))
}
beta_mixture2_moment_match <- function(y,
                                                  weights,
                                                  shape_floor = 1e-3,
                                                  shape_ceiling = 1e4) {
  y <- rotational_clamp_unit_interval(y, eps = 1e-8)
  weights <- jp_normalize_probability_weights(weights, length(y))
  m <- min(max(rotational_weighted_mean(y, weights), 1e-4), 1 - 1e-4)
  v <- rotational_weighted_variance(y, weights, center = m)
  max_var <- m * (1 - m) * (1 - 1e-6)
  v <- min(max(v, 1e-6), max_var)
  precision <- min(max(m * (1 - m) / v - 1, 2 * shape_floor), shape_ceiling)

  list(
    alpha = min(max(m * precision, shape_floor), shape_ceiling),
    beta = min(max((1 - m) * precision, shape_floor), shape_ceiling)
  )
}
beta_mixture2_start_thetas_s2 <- function(x,
                                                      weights = NULL,
                                                      control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  warm_start <- control$beta_mixture2_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- beta_mixture2_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$beta_mixture2_warm_start_only %||% FALSE)) {
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
    y <- rotational_clamp_unit_interval((pmin(pmax(as.numeric(x %*% mu0), -1), 1) + 1) / 2, eps = 1e-8)
    for (split_prob in split_probs) {
      threshold <- rotational_weighted_quantile(y, prob_weights, prob = split_prob)
      group1 <- y <= threshold
      if (all(group1) || !any(group1)) {
        next
      }

      w1 <- prob_weights[group1]
      w2 <- prob_weights[!group1]
      p1 <- sum(w1)
      if (p1 <= 1e-6 || p1 >= 1 - 1e-6) {
        next
      }

      comp1 <- beta_mixture2_moment_match(y[group1], w1 / sum(w1))
      comp2 <- beta_mixture2_moment_match(y[!group1], w2 / sum(w2))
      out[[length(out) + 1L]] <- beta_mixture2_canonicalize_theta(list(
        mu = mu0,
        weight1 = p1,
        alpha1 = comp1$alpha,
        beta1 = comp1$beta,
        alpha2 = comp2$alpha,
        beta2 = comp2$beta
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- beta_mixture2_canonicalize_theta(list(
      mu = c(0, 0, 1),
      weight1 = 0.5,
      alpha1 = 2,
      beta1 = 5,
      alpha2 = 5,
      beta2 = 2
    ))
  }

  out
}
beta_mixture2_pack_par <- function(theta,
                                              control = list()) {
  theta <- beta_mixture2_normalize_theta(theta, ambient_dim = 3L)
  list(
    par = c(
      stats::qlogis(theta$weight1),
      log(theta$alpha1),
      log(theta$beta1),
      log(theta$alpha2),
      log(theta$beta2),
      theta$mu
    ),
    theta = theta
  )
}
beta_mixture2_unpack_par <- function(par,
                                                control = list()) {
  shape_lower <- as.numeric(control$beta_mixture2_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$beta_mixture2_shape_upper %||% 1e3)
  weight_eps <- as.numeric(control$beta_mixture2_weight_eps %||% 0.01)

  mu_raw <- par[6:8]
  mu_hat <- rotational_unit_vector_fallback(mu_raw)
  beta_mixture2_canonicalize_theta(list(
    mu = mu_hat,
    weight1 = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    alpha1 = rotational_positive_parameter(par[[2L]], lower = shape_lower, upper = shape_upper),
    beta1 = rotational_positive_parameter(par[[3L]], lower = shape_lower, upper = shape_upper),
    alpha2 = rotational_positive_parameter(par[[4L]], lower = shape_lower, upper = shape_upper),
    beta2 = rotational_positive_parameter(par[[5L]], lower = shape_lower, upper = shape_upper)
  ))
}
beta_mixture2_mle_s2_weighted <- function(x,
                                                      weights = NULL,
                                                      control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- beta_mixture2_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(length(candidate_thetas), as.integer(control$beta_mixture2_n_starts %||% 12L)))]

  objective <- function(par) {
    theta <- beta_mixture2_unpack_par(par, control = control)
    value <- -beta_mixture2_weighted_loglik_s2(
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$beta_mixture2_optim_method %||% "BFGS"
  optim_control <- control$beta_mixture2_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- beta_mixture2_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- beta_mixture2_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$beta_mixture2_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$beta_mixture2_warm_start_only <- FALSE
    return(beta_mixture2_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Beta-mixture weighted MLE failed for all starting values.")
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
beta_mixture2_legendre_expectations_one_beta <- function(alpha, beta, l_max, quad_n) {
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 ||
      length(beta) != 1L || !is.finite(beta) || beta <= 0) {
    stop("Beta shape parameters must be finite and strictly positive.")
  }
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  n_nodes <- max(quad_n, as.integer(ceiling((l_max + 1L) / 2)))
  jacobi <- beta_mixture2_gauss_jacobi(
    n = n_nodes,
    alpha = beta - 1,
    beta = alpha - 1
  )
  normalized_weights <- jacobi$weights / jacobi$total_mass
  legendre_matrix <- rotational_legendre_matrix(jacobi$nodes, l_max = l_max)
  expectations <- as.numeric(crossprod(legendre_matrix, normalized_weights))

  list(
    expectations = expectations,
    mass_error = abs(sum(normalized_weights) - 1),
    n_nodes = n_nodes
  )
}
beta_mixture2_legendre_coefficients <- function(theta,
                                                l_max = 150L,
                                                quad_n = 1000L,
                                                tol = 1e-10) {
  theta <- beta_mixture2_normalize_theta(theta, ambient_dim = 3L)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }
  if (length(quad_n) != 1L || !is.finite(quad_n) || quad_n < 1L) {
    stop("`quad_n` must be a strictly positive integer.")
  }

  comp1 <- beta_mixture2_legendre_expectations_one_beta(
    alpha = theta$alpha1,
    beta = theta$beta1,
    l_max = l_max,
    quad_n = quad_n
  )
  comp2 <- beta_mixture2_legendre_expectations_one_beta(
    alpha = theta$alpha2,
    beta = theta$beta2,
    l_max = l_max,
    quad_n = quad_n
  )

  expectations <- theta$weight1 * comp1$expectations +
    (1 - theta$weight1) * comp2$expectations
  ell <- 0:l_max
  coeffs <- (2 * ell + 1) * expectations
  coeffs[[1L]] <- 1

  if (any(!is.finite(coeffs))) {
    stop("Beta-mixture Gauss-Jacobi Legendre coefficient computation produced nonfinite coefficients.")
  }

  a0_error <- abs(expectations[[1L]] - 1)
  if (a0_error > tol) {
    stop(sprintf("Beta-mixture Gauss-Jacobi coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error,
    mass_error = max(comp1$mass_error, comp2$mass_error),
    method = "gauss_jacobi",
    quad_n = max(comp1$n_nodes, comp2$n_nodes)
  )
}
r_sph_beta_mixture2 <- function(n,
                                           mu,
                                           weight1,
                                           alpha1,
                                           beta1,
                                           alpha2,
                                           beta2,
                                           check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  component1 <- stats::runif(n) <= theta$weight1
  y <- numeric(n)
  n1 <- sum(component1)
  n2 <- n - n1
  if (n1 > 0L) {
    y[component1] <- stats::rbeta(n1, theta$alpha1, theta$beta1)
  }
  if (n2 > 0L) {
    y[!component1] <- stats::rbeta(n2, theta$alpha2, theta$beta2)
  }

  z <- pmin(pmax(2 * y - 1, -1), 1)
  phi <- stats::runif(n, min = 0, max = 2 * pi)
  basis <- jp_orthonormal_complement(theta$mu)
  tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
  radial <- sqrt(pmax(0, 1 - z^2))
  x <- tcrossprod(z, theta$mu) + sweep(tangent, 1L, radial, FUN = "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Beta-mixture sampler returned non-unit vectors.")
    }
  }
  x
}
uniform_beta_mixture_validate_parameters <- function(mu,
                                                     weight_uniform,
                                                     alpha,
                                                     beta) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Rotational uniform-beta-mixture utilities currently support only S^2.")
  }

  weight_uniform <- as.numeric(weight_uniform)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)

  if (length(weight_uniform) != 1L || !is.finite(weight_uniform) ||
      weight_uniform <= 0 || weight_uniform >= 1) {
    stop("`weight_uniform` must be a finite scalar in (0, 1).")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 ||
      length(beta) != 1L || !is.finite(beta) || beta <= 0) {
    stop("Uniform-beta-mixture shape parameters must be finite and strictly positive.")
  }

  list(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta,
    ambient_dim = 3L
  )
}
uniform_beta_mixture_density_y <- function(y,
                                           weight_uniform,
                                           alpha,
                                           beta,
                                           log = FALSE,
                                           eps = 1e-12) {
  y <- rotational_clamp_unit_interval(as.numeric(y), eps = eps)
  log_density <- rotational_logsumexp2(
    log(weight_uniform),
    log1p(-weight_uniform) + stats::dbeta(y, shape1 = alpha, shape2 = beta, log = TRUE)
  )
  if (log) log_density else exp(log_density)
}
uniform_beta_mixture_cdf_y <- function(y,
                                       weight_uniform,
                                       alpha,
                                       beta) {
  y <- as.numeric(y)
  out <- numeric(length(y))
  out[y <= 0] <- 0
  out[y >= 1] <- 1
  active <- which(y > 0 & y < 1)
  if (length(active) == 0L) {
    return(out)
  }

  out[active] <- weight_uniform * y[active] +
    (1 - weight_uniform) * stats::pbeta(y[active], shape1 = alpha, shape2 = beta)
  pmin(pmax(out, 0), 1)
}
uniform_beta_mixture_density_h <- function(z,
                                           weight_uniform,
                                           alpha,
                                           beta,
                                           eps = 1e-12) {
  z <- as.numeric(z)
  y <- rotational_clamp_unit_interval((z + 1) / 2, eps = eps)
  out <- numeric(length(z))
  valid <- z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- uniform_beta_mixture_density_y(
    y = y[valid],
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta,
    log = FALSE,
    eps = eps
  )
  out
}
uniform_beta_mixture_density_gz <- function(z,
                                            weight_uniform,
                                            alpha,
                                            beta,
                                            eps = 1e-12) {
  0.5 * uniform_beta_mixture_density_h(
    z = z,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta,
    eps = eps
  )
}
d_sph_uniform_beta_mixture_s2 <- function(x,
                                          mu,
                                          weight_uniform,
                                          alpha,
                                          beta,
                                          log = FALSE,
                                          eps = 1e-12) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )
  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  log_density <- -log(4 * pi) + uniform_beta_mixture_density_y(
    y = y,
    weight_uniform = params$weight_uniform,
    alpha = params$alpha,
    beta = params$beta,
    log = TRUE,
    eps = eps
  )
  if (log) log_density else exp(log_density)
}
uniform_beta_mixture_weighted_loglik_s2 <- function(mu,
                                                    weight_uniform,
                                                    alpha,
                                                    beta,
                                                    x,
                                                    prob_weights = NULL,
                                                    eps = 1e-12) {
  params <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  sum(prob_weights * uniform_beta_mixture_density_y(
    y = y,
    weight_uniform = params$weight_uniform,
    alpha = params$alpha,
    beta = params$beta,
    log = TRUE,
    eps = eps
  ))
}
uniform_beta_mixture_start_thetas_s2 <- function(x,
                                                 weights = NULL,
                                                 control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  warm_start <- control$uniform_beta_mixture_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- uniform_beta_mixture_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$uniform_beta_mixture_warm_start_only %||% FALSE)) {
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

  split_probs <- c(0.15, 0.25, 0.35, 0.5)
  for (mu0 in mu_candidates) {
    y <- rotational_clamp_unit_interval((pmin(pmax(as.numeric(x %*% mu0), -1), 1) + 1) / 2, eps = 1e-8)
    for (split_prob in split_probs) {
      threshold <- rotational_weighted_quantile(y, prob_weights, prob = split_prob)
      beta_group <- y > threshold
      if (all(beta_group) || !any(beta_group)) {
        next
      }

      w_beta <- prob_weights[beta_group]
      weight_uniform0 <- sum(prob_weights[!beta_group])
      if (weight_uniform0 <= 1e-6 || weight_uniform0 >= 1 - 1e-6) {
        next
      }

      beta_start <- beta_mixture2_moment_match(y[beta_group], w_beta / sum(w_beta))
      out[[length(out) + 1L]] <- uniform_beta_mixture_normalize_theta(list(
        mu = mu0,
        weight_uniform = weight_uniform0,
        alpha = beta_start$alpha,
        beta = beta_start$beta
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- uniform_beta_mixture_normalize_theta(list(
      mu = c(0, 0, 1),
      weight_uniform = 0.2,
      alpha = 8,
      beta = 2
    ))
  }

  out
}
uniform_beta_mixture_pack_par <- function(theta,
                                          control = list()) {
  theta <- uniform_beta_mixture_normalize_theta(theta, ambient_dim = 3L)
  list(
    par = c(
      stats::qlogis(theta$weight_uniform),
      log(theta$alpha),
      log(theta$beta),
      theta$mu
    ),
    theta = theta
  )
}
uniform_beta_mixture_unpack_par <- function(par,
                                            control = list()) {
  shape_lower <- as.numeric(control$uniform_beta_mixture_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$uniform_beta_mixture_shape_upper %||% 1e3)
  weight_eps <- as.numeric(control$uniform_beta_mixture_weight_eps %||% 0.01)

  mu_raw <- par[4:6]
  mu_hat <- rotational_unit_vector_fallback(mu_raw)
  uniform_beta_mixture_normalize_theta(list(
    mu = mu_hat,
    weight_uniform = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    alpha = rotational_positive_parameter(par[[2L]], lower = shape_lower, upper = shape_upper),
    beta = rotational_positive_parameter(par[[3L]], lower = shape_lower, upper = shape_upper)
  ))
}
uniform_beta_mixture_mle_s2_weighted <- function(x,
                                                 weights = NULL,
                                                 control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- uniform_beta_mixture_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(length(candidate_thetas), as.integer(control$uniform_beta_mixture_n_starts %||% 12L)))]

  objective <- function(par) {
    theta <- uniform_beta_mixture_unpack_par(par, control = control)
    value <- -uniform_beta_mixture_weighted_loglik_s2(
      mu = theta$mu,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      x = x,
      prob_weights = prob_weights,
      eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$uniform_beta_mixture_optim_method %||% "BFGS"
  optim_control <- control$uniform_beta_mixture_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- uniform_beta_mixture_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- uniform_beta_mixture_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$uniform_beta_mixture_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$uniform_beta_mixture_warm_start_only <- FALSE
    return(uniform_beta_mixture_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Uniform-beta-mixture weighted MLE failed for all starting values.")
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
uniform_beta_mixture_legendre_coefficients <- function(theta,
                                                       l_max = 150L,
                                                       quad_n = 1000L,
                                                       tol = 1e-10) {
  theta <- uniform_beta_mixture_normalize_theta(theta, ambient_dim = 3L)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }
  if (length(quad_n) != 1L || !is.finite(quad_n) || quad_n < 1L) {
    stop("`quad_n` must be a strictly positive integer.")
  }

  beta_part <- beta_mixture2_legendre_expectations_one_beta(
    alpha = theta$alpha,
    beta = theta$beta,
    l_max = l_max,
    quad_n = quad_n
  )
  ell <- 0:l_max
  coeffs <- numeric(l_max + 1L)
  coeffs[[1L]] <- 1
  if (l_max >= 1L) {
    coeffs[-1L] <- (1 - theta$weight_uniform) *
      (2 * ell[-1L] + 1) * beta_part$expectations[-1L]
  }

  if (any(!is.finite(coeffs))) {
    stop("Uniform-beta-mixture Legendre coefficient computation produced nonfinite coefficients.")
  }

  a0_error <- abs(beta_part$expectations[[1L]] - 1)
  if (a0_error > tol) {
    stop(sprintf("Uniform-beta-mixture coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error,
    mass_error = beta_part$mass_error,
    method = "gauss_jacobi",
    quad_n = beta_part$n_nodes
  )
}
r_sph_uniform_beta_mixture <- function(n,
                                       mu,
                                       weight_uniform,
                                       alpha,
                                       beta,
                                       check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  x <- matrix(0, nrow = n, ncol = 3L)
  component_uniform <- stats::runif(n) <= theta$weight_uniform
  n_uniform <- sum(component_uniform)
  n_beta <- n - n_uniform

  if (n_uniform > 0L) {
    x[component_uniform, ] <- jp_uniform_sphere(n = n_uniform, ambient_dim = 3L)
  }
  if (n_beta > 0L) {
    y <- stats::rbeta(n_beta, shape1 = theta$alpha, shape2 = theta$beta)
    z <- pmin(pmax(2 * y - 1, -1), 1)
    phi <- stats::runif(n_beta, min = 0, max = 2 * pi)
    basis <- jp_orthonormal_complement(theta$mu)
    tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
    radial <- sqrt(pmax(0, 1 - z^2))
    x[!component_uniform, ] <- tcrossprod(z, theta$mu) + sweep(tangent, 1L, radial, FUN = "*")
  }

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Uniform-beta-mixture sampler returned non-unit vectors.")
    }
  }

  x
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
A_q <- function(kappa, q) {
  if (kappa == 0) return(0)
  nu1 <- (q + 1) / 2
  nu2 <- (q - 1) / 2
  I_nu1 <- besselI(kappa, nu = nu1, expon.scaled = TRUE)
  I_nu2 <- besselI(kappa, nu = nu2, expon.scaled = TRUE)
  return(I_nu1 / I_nu2)
}
A_q_prime <- function(kappa, q) {
  if (kappa == 0) return(0)
  nu1 <- (q + 1) / 2
  nu2 <- (q - 1) / 2
  nu3 <- (q - 3) / 2
  I_nu1 <- besselI(kappa, nu = nu1, expon.scaled = TRUE)
  I_nu2 <- besselI(kappa, nu = nu2, expon.scaled = TRUE)
  I_nu3 <- besselI(kappa, nu = nu3, expon.scaled = TRUE)
  I_nu1_prime <- I_nu2 - (nu1 / kappa) * I_nu1
  I_nu2_prime <- I_nu3 - (nu2 / kappa) * I_nu2
  numerator <- I_nu1_prime * I_nu2 - I_nu1 * I_nu2_prime
  denominator <- I_nu2^2
  return(numerator / denominator)
}
psi_xi <- function(x, xi, q) {
  kappa <- norm(xi, type = "2")
  if (kappa == 0) return(x)
  mu <- xi / kappa
  A_q_val <- A_q(kappa, q)
  return(x - A_q_val * mu)
}
dot_psi_xi <- function(xi, q) {
  kappa <- norm(xi, type = "2")
  if (kappa == 0) {
    return(-diag(q + 1))
  }
  mu <- xi / kappa
  A_q_val <- A_q(kappa, q)
  A_q_prime_val <- A_q_prime(kappa, q)
  I_q <- diag(q + 1)
  mu_outer <- outer(mu, mu)
  return(-(A_q_prime_val * mu_outer + (A_q_val / kappa) * (I_q - mu_outer)))
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
