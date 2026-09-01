# Fixed-concentration von Mises--Fisher adapter for multiplier-bootstrap GOF.
#
# This is intentionally separate from vmf_model_spec.R.  It implements the
# composite model {vMF(mu, kappa): mu in S^q}, for a supplied kappa > 0, so
# that only the location direction is fitted.  The usual vMF adapter continues
# to fit xi = kappa * mu and is not modified by this file.

if (!exists("make_vmf_spec", mode = "function")) {
  vmf_candidates_fixed_kappa <- c(
    file.path("bootstrap", "vmf_model_spec.R"), "vmf_model_spec.R",
    file.path("..", "bootstrap", "vmf_model_spec.R"),
    file.path("..", "..", "bootstrap", "vmf_model_spec.R")
  )
  vmf_path_fixed_kappa <- vmf_candidates_fixed_kappa[file.exists(vmf_candidates_fixed_kappa)][1L]
  if (is.na(vmf_path_fixed_kappa)) stop("Could not locate bootstrap/vmf_model_spec.R.")
  source(vmf_path_fixed_kappa)
}

vmf_fixed_kappa_validate <- function(kappa) {
  kappa <- as.numeric(kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("The fixed vMF concentration `kappa` must be a strictly positive finite scalar.", call. = FALSE)
  }
  kappa
}

# Deterministic orthonormal basis of the tangent space at mu.  Dropping the
# coordinate most aligned with mu leaves a full-rank projected canonical basis.
vmf_fixed_kappa_tangent_basis <- function(mu) {
  mu <- as.numeric(mu)
  if (length(mu) < 2L || any(!is.finite(mu))) {
    stop("`mu` must be a finite vector of length at least two.", call. = FALSE)
  }
  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 0) stop("`mu` must have positive norm.", call. = FALSE)
  mu <- mu / mu_norm
  ambient_dim <- length(mu)
  anchor <- which.max(abs(mu))
  canonical <- diag(ambient_dim)[, -anchor, drop = FALSE]
  projected <- canonical - tcrossprod(mu, drop(crossprod(mu, canonical)))
  decomposition <- qr(projected, tol = 1e-12)
  if (decomposition$rank != ambient_dim - 1L) {
    stop("Could not construct a full-rank tangent basis for the fixed-kappa vMF model.", call. = FALSE)
  }
  qr.Q(decomposition, complete = FALSE)
}

normalize_vmf_fixed_kappa_theta <- function(theta, kappa, ambient_dim = NULL) {
  kappa <- vmf_fixed_kappa_validate(kappa)
  mu <- if (is.list(theta)) theta$mu else theta
  mu <- as.numeric(mu)
  if (length(mu) < 2L || any(!is.finite(mu))) {
    stop("Fixed-kappa vMF theta requires a finite location vector `mu`.", call. = FALSE)
  }
  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 1e-12) {
    stop("Fixed-kappa vMF MLE has a zero resultant direction.", call. = FALSE)
  }
  mu <- mu / mu_norm
  if (!is.null(ambient_dim) && length(mu) != as.integer(ambient_dim)) {
    stop("Fixed-kappa vMF theta has incompatible ambient dimension.", call. = FALSE)
  }
  list(mu = mu, kappa = kappa, q = length(mu) - 1L)
}

fit_vmf_fixed_kappa_theta <- function(data, weights = NULL, null, control = list()) {
  x <- normalize_vmf_data(data, control)
  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must specify `type` for the fixed-kappa vMF model.", call. = FALSE)
  }
  if (identical(null$type, "simple")) {
    if (is.null(null$theta)) stop("A simple fixed-kappa vMF null requires `null$theta`.", call. = FALSE)
    theta_kappa <- if (is.list(null$theta)) null$theta$kappa else NULL
    kappa <- null$kappa %||% theta_kappa
    return(normalize_vmf_fixed_kappa_theta(null$theta, kappa, ncol(x)))
  }
  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.", call. = FALSE)
  }
  kappa <- vmf_fixed_kappa_validate(null$kappa)
  probabilities <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }
  resultant <- colSums(x * probabilities)
  normalize_vmf_fixed_kappa_theta(resultant, kappa, ncol(x))
}

vmf_fixed_kappa_score_matrix <- function(data, theta) {
  x <- normalize_vmf_data(data)
  theta <- normalize_vmf_fixed_kappa_theta(theta, theta$kappa, ncol(x))
  tangent_basis <- vmf_fixed_kappa_tangent_basis(theta$mu)
  theta$kappa * x %*% tangent_basis
}

vmf_fixed_kappa_information <- function(theta) {
  theta <- normalize_vmf_fixed_kappa_theta(theta, theta$kappa)
  scalar <- theta$kappa * A_q(theta$kappa, theta$q)
  if (!is.finite(scalar) || scalar <= 0) {
    stop("Fixed-kappa vMF Fisher information is non-positive or non-finite.", call. = FALSE)
  }
  diag(scalar, theta$q)
}

vmf_fixed_kappa_profile_derivative <- function(omega, thresholds, theta,
                                               distance_type, grid_size) {
  theta <- normalize_vmf_fixed_kappa_theta(theta, theta$kappa)
  tangent_basis <- vmf_fixed_kappa_tangent_basis(theta$mu)
  derivative_xi <- vmf_profile_and_derivative_xi(
    omega = omega,
    xi = theta$kappa * theta$mu,
    t_values = thresholds,
    distance_type = distance_type,
    grid_size = grid_size
  )$derivative
  derivative_xi %*% (theta$kappa * tangent_basis)
}

prepare_vmf_fixed_kappa_fast_multiplier <- function(data, theta_hat, spec,
                                                      ks_prep = NULL, cvm_prep = NULL,
                                                      control = list(),
                                                      distance_type = "geodesic") {
  x <- normalize_vmf_data(data, control)
  theta_hat <- normalize_vmf_fixed_kappa_theta(theta_hat, theta_hat$kappa, ncol(x))
  grid_size <- as.integer(control$vmf_derivative_n_u %||% control$vmf_profile_n_u %||% 4097L)
  if (!is.finite(grid_size) || grid_size < 33L) {
    stop("`vmf_derivative_n_u` must be an integer of at least 33.", call. = FALSE)
  }

  S_obs <- vmf_fixed_kappa_score_matrix(x, theta_hat)
  Vhat <- vmf_fixed_kappa_information(theta_hat)
  evaluate_derivative <- function(omega, thresholds) {
    vmf_fixed_kappa_profile_derivative(
      omega = omega, thresholds = thresholds, theta = theta_hat,
      distance_type = distance_type, grid_size = grid_size
    )
  }
  stack_common_thresholds <- function(centers, thresholds) {
    centers <- normalize_vmf_data(centers, control)
    do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
      evaluate_derivative(centers[i, ], thresholds)
    }))
  }

  D_ks <- if (is.null(ks_prep)) {
    NULL
  } else if (identical(ks_prep$ks_grid_mode %||% "fixed", "sample_points_unique_distances")) {
    centers <- normalize_vmf_data(ks_prep$omega_grid, control)
    list(
      mode = "sample_points_unique_distances",
      derivative_sorted = profile_derivative_stack_centers(
        centers = centers,
        thresholds = ks_prep$sorted_distance_matrix,
        evaluator = evaluate_derivative
      )
    )
  } else {
    stack_common_thresholds(ks_prep$omega_grid, ks_prep$t_grid)
  }

  D_cvm <- if (is.null(cvm_prep)) {
    NULL
  } else if (isTRUE(cvm_prep$shared_with_ks) && is.list(D_ks) && !is.null(D_ks$derivative_sorted)) {
    list(
      mode = "sample_points_unique_distances_sorted_rows",
      derivative_sorted = D_ks$derivative_sorted,
      shared_with_ks = TRUE
    )
  } else if (isTRUE(cvm_prep$light) && !is.null(cvm_prep$sorted_distance_matrix)) {
    list(
      mode = "sample_points_unique_distances_sorted_rows",
      derivative_sorted = profile_derivative_stack_centers(
        centers = x,
        thresholds = cvm_prep$sorted_distance_matrix,
        evaluator = evaluate_derivative
      )
    )
  } else {
    profile_derivative_stack_centers(
      centers = x,
      thresholds = cvm_prep$distance_matrix %||% spec$distance_matrix(x, x, control),
      evaluator = evaluate_derivative
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = matrix(numeric(0), nrow = 0L, ncol = ncol(S_obs)),
    vhat_method = "analytic_fixed_kappa_fisher_information",
    vhat_diagnostics = fast_multiplier_deterministic_vhat_diagnostics(
      S_obs = S_obs, Vhat = Vhat, par0 = theta_hat$mu
    ),
    correction_representation = "score",
    derivative_method = "quadrature",
    derivative_mc_size = NA_integer_,
    derivative_mc_seed = NA_integer_,
    derivative_grid_size = grid_size,
    observed_cvm_distance_matrix = if (!is.null(cvm_prep) && !isTRUE(cvm_prep$light)) {
      cvm_prep$distance_matrix %||% spec$distance_matrix(x, x, control)
    } else {
      NULL
    },
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

make_vmf_fixed_kappa_spec <- function(kappa,
                                      distance_type = c("chordal", "geodesic")) {
  kappa <- vmf_fixed_kappa_validate(kappa)
  distance_type <- match.arg(distance_type)
  base_spec <- make_vmf_spec(distance_type = distance_type, unknown_param = "xi")
  new_model_spec(
    name = sprintf("vmf_fixed_kappa_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      if (identical(null$type, "composite")) null$kappa <- kappa
      fit_vmf_fixed_kappa_theta(data, weights = weights, null = null, control = control)
    },
    distance_matrix = base_spec$distance_matrix,
    profile_eval = base_spec$profile_eval,
    normalize_data = normalize_vmf_data,
    n_obs = base_spec$n_obs,
    observation_at = base_spec$observation_at,
    extras = list(
      profile_matrix_eval = base_spec$profile_matrix_eval,
      sample_profile_matrix_eval = base_spec$sample_profile_matrix_eval,
      fast_multiplier_prepare = function(data, theta_hat, ks_prep = NULL, cvm_prep = NULL,
                                         control = list()) {
        prepare_vmf_fixed_kappa_fast_multiplier(
          data = data, theta_hat = theta_hat,
          spec = make_vmf_fixed_kappa_spec(kappa, distance_type),
          ks_prep = ks_prep, cvm_prep = cvm_prep, control = control,
          distance_type = distance_type
        )
      },
      distance_type = distance_type,
      unknown_param = "mu_fixed_kappa",
      fixed_kappa = kappa
    )
  )
}

multiplier_bootstrap_vmf_fixed_kappa <- function(data,
                                                 kappa,
                                                 null = list(type = "composite"),
                                                 statistics = c("ks", "cvm"),
                                                 ks_grid = NULL,
                                                 B = 5000L,
                                                 alpha = 0.05,
                                                 multipliers = NULL,
                                                 n_cores = 1L,
                                                 seed = NULL,
                                                 bootstrap_method = c("reestimated", "fast_multiplier"),
                                                 keep = list(observed_process = FALSE,
                                                             bootstrap_statistics = FALSE,
                                                             bootstrap_thetas = FALSE),
                                                 control = list(),
                                                 distance_type = c("chordal", "geodesic"),
                                                 distance_profile_backend = c("r", "cpp"),
                                                 fast_multiplier_backend = c("cpp", "r"),
                                                 fast_multiplier_cpp_kernel = c("contiguous_double", "legacy"),
                                                 fuse_ks_cvm = TRUE,
                                                 cache_block_corrections = c("auto", "true", "false")) {
  kappa <- vmf_fixed_kappa_validate(kappa)
  distance_type <- match.arg(distance_type)
  distance_profile_backend <- match.arg(distance_profile_backend)
  fast_multiplier_backend <- normalize_fast_multiplier_backend(fast_multiplier_backend[[1L]])
  fast_multiplier_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(fast_multiplier_cpp_kernel[[1L]])
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(fuse_ks_cvm)
  cache_block_corrections <- normalize_fast_multiplier_cache(cache_block_corrections[[1L]])
  control$derivative_method <- "quadrature"
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_cpp_kernel <- fast_multiplier_cpp_kernel
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  null$type <- null$type %||% "composite"
  if (identical(null$type, "composite")) null$kappa <- kappa

  result <- multiplier_bootstrap_gof(
    data = data, spec = make_vmf_fixed_kappa_spec(kappa, distance_type), null = null,
    statistics = statistics, ks_grid = ks_grid, B = B, alpha = alpha,
    multipliers = multipliers, n_cores = n_cores, seed = seed,
    bootstrap_method = bootstrap_method, keep = keep, control = control,
    distance_profile_backend = distance_profile_backend
  )
  result$diagnostics$fixed_kappa <- kappa
  result$diagnostics$unknown_param <- "mu_fixed_kappa"
  result$diagnostics$derivative_method_requested <- "quadrature"
  result$diagnostics$fast_multiplier_backend_requested <- fast_multiplier_backend
  result$diagnostics$fast_multiplier_cpp_kernel_requested <- fast_multiplier_cpp_kernel
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <- fuse_ks_cvm
  result
}
