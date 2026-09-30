# Functions required by the joint sunspot model adapter.
# Extracted from the three real_data/sunspots model modules; unused runner code omitted.

fit_sunspots_cycle23_joint_time_space <- function(x, s,
                                                   hemisphere_regression = "asymmetric",
                                                   control = list(),
                                                   eta_hat = NULL,
                                                   theta_hat = NULL) {
  data <- sunspots_joint_validate_data(x, s)
  if (is.null(eta_hat)) eta_hat <- fit_sunspots_joint_time_beta_mixture2(data$s, control)
  if (is.null(theta_hat)) {
    theta_hat <- fit_sunspots_time_varying_asymmetric_mixture(
      data$x, data$s, hemisphere_regression = hemisphere_regression, control = control
    )
  }
  temporal_loglik <- sum(sunspots_joint_time_log_density(data$s, eta_hat, control))
  conditional_loglik <- sunspots_time_varying_loglik(data$x, data$s, theta_hat)
  n_parameters <- 5L + theta_hat$n_parameters
  list(
    eta_hat = eta_hat,
    theta_hat = theta_hat,
    temporal_loglik = temporal_loglik,
    conditional_loglik = conditional_loglik,
    loglik = temporal_loglik + conditional_loglik,
    n_parameters = n_parameters,
    aic = 2 * n_parameters - 2 * (temporal_loglik + conditional_loglik),
    bic = n_parameters * log(nrow(data$x)) - 2 * (temporal_loglik + conditional_loglik)
  )
}
sunspots_time_varying_normalize_hemisphere_regression <- function(hemisphere_regression = "asymmetric") {
  hemisphere_regression <- tolower(as.character(hemisphere_regression))
  if (length(hemisphere_regression) != 1L ||
      !hemisphere_regression %in% c("asymmetric", "shared")) {
    stop("`hemisphere_regression` must be either 'asymmetric' or 'shared'.")
  }
  hemisphere_regression
}
sunspots_joint_validate_data <- function(x, s) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  if (ncol(x) != 3L) stop("The joint sunspots model supports only S^2 data.")
  s <- as.numeric(s)
  if (length(s) != nrow(x) || any(!is.finite(s)) || any(s <= 0 | s >= 1)) {
    stop("`s` must contain one finite value in (0, 1) for each row of `x`.")
  }
  list(x = x, s = s)
}
sunspots_joint_pack_par <- function(fit, control = list()) {
  c(
    sunspots_joint_time_pack_eta(fit$eta_hat, control),
    sunspots_time_varying_pack_theta(
      fit$theta_hat,
      nu_eps = as.numeric(control$nu_eps %||% 1e-6),
      c_min = as.numeric(control$c_min %||% 1e-8),
      c_max = as.numeric(control$c_max %||% 1e6),
      hemisphere_regression = control$hemisphere_regression %||% "asymmetric"
    )
  )
}
sunspots_joint_score_matrix <- function(data, par, control = list()) {
  data <- sunspots_joint_validate_data(data$x, data$s)
  state <- sunspots_joint_state_from_par(par, control)
  cbind(
    sunspots_joint_time_score_matrix(data$s, state$eta_par, control),
    sunspots_time_gof_score_matrix(list(x = data$x, u = data$s), state$theta_par, control)
  )
}
sample_sunspots_joint_time_space <- function(n, par, control = list()) {
  n <- as.integer(n)
  if (!is.finite(n) || n <= 0L) stop("`n` must be a positive integer.")
  state <- sunspots_joint_state_from_par(par, control)
  s <- sample_sunspots_joint_time_beta_mixture2(n, state$eta, control)
  nu <- sunspots_time_varying_nu(s, state$theta)
  north <- stats::runif(n) <= 0.5
  x <- matrix(0, nrow = n, ncol = 3L)
  for (i in seq_len(n)) {
    x[i, ] <- r_sph_small_circle(
      n = 1L, mu = if (north[[i]]) c(0, 0, 1) else c(0, 0, -1),
      kappa = state$theta$c, nu = if (north[[i]]) nu$north[[i]] else nu$south[[i]], check = FALSE
    )
  }
  list(x = x, s = s)
}
sunspots_joint_time_quadrature <- function(eta, n_nodes = 64L, control = list()) {
  n_nodes <- as.integer(n_nodes)
  if (!is.finite(n_nodes) || n_nodes < 2L) stop("`n_nodes` must be an integer >= 2.")
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  component_rule <- function(alpha, beta, mixture_weight) {
    rule <- beta_mixture2_gauss_jacobi(n_nodes, alpha = beta - 1, beta = alpha - 1)
    list(
      nodes = (rule$nodes + 1) / 2,
      weights = mixture_weight * rule$weights / rule$total_mass
    )
  }
  first <- component_rule(eta$alpha1, eta$beta1, eta$weight1)
  second <- component_rule(eta$alpha2, eta$beta2, 1 - eta$weight1)
  order_index <- order(c(first$nodes, second$nodes))
  list(
    nodes = c(first$nodes, second$nodes)[order_index],
    weights = c(first$weights, second$weights)[order_index],
    n_nodes_per_component = n_nodes,
    mass_error = abs(sum(c(first$weights, second$weights)) - 1)
  )
}
sunspots_joint_conditional_legendre_coefficients <- function(theta, time_nodes,
                                                             l_max = 100L, quad_n = 400L) {
  theta <- sunspots_time_varying_validate_theta(theta)
  time_nodes <- as.numeric(time_nodes)
  if (length(time_nodes) == 0L || any(!is.finite(time_nodes)) || any(time_nodes <= 0 | time_nodes >= 1)) {
    stop("`time_nodes` must be finite and lie in (0, 1).")
  }
  quadrature <- small_circle_gauss_legendre(as.integer(quad_n))
  z_nodes <- quadrature$nodes
  legendre <- small_circle_legendre_matrix(z_nodes, l_max = as.integer(l_max))
  nu <- sunspots_time_varying_nu(time_nodes, theta)
  north_log_density <- sweep(
    -theta$c * outer(z_nodes, nu$north, FUN = "-")^2,
    2L, log(2) + small_circle_log_norm_constant(theta$c, nu$north), FUN = "-"
  )
  south_log_density <- sweep(
    -theta$c * outer(-z_nodes, nu$south, FUN = "-")^2,
    2L, log(2) + small_circle_log_norm_constant(theta$c, nu$south), FUN = "-"
  )
  rotational_density <- exp(north_log_density) + exp(south_log_density)
  ell <- 0:as.integer(l_max)
  moments <- crossprod(legendre, sweep(rotational_density, 1L, quadrature$weights, FUN = "*"))
  coefficients <- t(sweep(moments, 1L, (2 * ell + 1) / 2, FUN = "*"))
  coefficients[, 1L] <- 1
  coefficients
}
sunspots_joint_profile_block <- function(radii, rho, center_s, time_nodes,
                                         time_weights, coefficients, backend = "auto") {
  backend <- sunspots_joint_effective_backend(backend)
  if (identical(backend, "r")) {
    return(sunspots_joint_profile_block_r(radii, rho, center_s, time_nodes, time_weights, coefficients))
  }
  with_distance_profile_backend("cpp", distance_profile_cpp_call(
    "cpp_sunspots_joint_profile_block",
    as.matrix(radii), as.numeric(rho), as.numeric(center_s), as.numeric(time_nodes),
    as.numeric(time_weights), as.matrix(coefficients)
  ))
}
sunspots_joint_profile_block_sorted <- function(
    radii,
    rho,
    center_s,
    time_nodes,
    time_weights,
    coefficients,
    backend = "auto") {
  backend <- sunspots_joint_effective_backend(backend)
  if (identical(backend, "r")) {
    return(sunspots_joint_profile_block_r(
      radii,
      rho,
      center_s,
      time_nodes,
      time_weights,
      coefficients
    ))
  }

  with_distance_profile_backend(
    "cpp",
    distance_profile_cpp_call(
      "cpp_sunspots_joint_profile_block_sorted",
      as.matrix(radii),
      as.numeric(rho),
      as.numeric(center_s),
      as.numeric(time_nodes),
      as.numeric(time_weights),
      as.matrix(coefficients)
    )
  )
}
sunspots_joint_effective_backend <- function(backend = "auto") {
  backend <- sunspots_joint_normalize_backend(backend)
  if (!identical(backend, "auto")) return(backend)
  if (!requireNamespace("Rcpp", quietly = TRUE)) return("r")
  loaded <- try(ensure_distance_profile_cpp_loaded(), silent = TRUE)
  if (inherits(loaded, "try-error")) "r" else "cpp"
}
fit_sunspots_joint_time_beta_mixture2 <- function(s, control = list()) {
  s <- as.numeric(s)
  if (length(s) < 4L || any(!is.finite(s)) || any(s <= 0 | s >= 1)) {
    stop("`s` must contain at least four values in (0, 1).")
  }
  starts <- sunspots_joint_time_start_etas(s, control)
  n_starts <- as.integer(control$time_beta_n_starts %||% length(starts))
  if (!is.finite(n_starts) || n_starts < 1L) stop("`time_beta_n_starts` must be a positive integer.")
  starts <- starts[seq_len(min(length(starts), n_starts))]
  nelder_mead_control <- control$time_beta_nelder_mead_control %||%
    list(maxit = 4000L, reltol = 1e-10)
  lbfgsb_control <- control$time_beta_optim_control %||%
    list(maxit = 1500L, factr = 1e7, pgtol = 1e-8)
  bounds <- sunspots_joint_time_control(control)
  lower <- c(stats::qlogis(bounds$weight_eps), rep(log(bounds$shape_lower), 4L))
  upper <- c(stats::qlogis(1 - bounds$weight_eps), rep(log(bounds$shape_upper), 4L))
  objective <- function(par) {
    eta <- sunspots_joint_time_unpack_eta(par, control)
    value <- -sum(sunspots_joint_time_log_density(s, eta, control))
    if (is.finite(value)) value else .Machine$double.xmax / 100
  }
  fits <- unlist(lapply(starts, function(start) {
    par0 <- pmin(pmax(sunspots_joint_time_pack_eta(start, control), lower), upper)
    exploratory <- try(stats::optim(
      par0, objective, method = "Nelder-Mead", control = nelder_mead_control
    ), silent = TRUE)
    candidates <- list(par0)
    if (!inherits(exploratory, "try-error") && is.finite(exploratory$value)) {
      candidates[[length(candidates) + 1L]] <- pmin(pmax(exploratory$par, lower), upper)
    }
    lapply(candidates, function(candidate) {
      try(stats::optim(
        candidate, objective, method = "L-BFGS-B", lower = lower, upper = upper,
        control = lbfgsb_control
      ), silent = TRUE)
    })
  }), recursive = FALSE)
  selection <- sunspots_joint_time_select_fit(fits)
  best <- selection$fit
  eta_hat <- sunspots_joint_time_unpack_eta(best$par, control)
  boundary <- sunspots_joint_time_boundary_flags(eta_hat, control)
  boundary_diagnostics <- sunspots_joint_time_parameter_boundaries(eta_hat, control)
  if (any(unlist(boundary, use.names = FALSE))) {
    warning(
      "The temporal beta-mixture MLE reached an admissible boundary; inspect the fit before using regular asymptotics.",
      call. = FALSE
    )
  }
  c(eta_hat, list(
    loglik = -best$value,
    opt = best,
    n_starts = length(starts),
    n_successful_starts = selection$n_finite_fits,
    n_converged_starts = selection$n_converged_fits,
    selected_converged = selection$selected_converged,
    boundary_flags = boundary,
    boundary_diagnostics = boundary_diagnostics
  ))
}
fit_sunspots_time_varying_asymmetric_mixture <- function(x,
                                                          u,
                                                          hemisphere_regression = "asymmetric",
                                                          control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  u <- as.numeric(u)
  if (length(u) != nrow(x)) stop("`u` must have one entry per row of `x`.")
  if (any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) stop("`u` must lie in (0, 1).")

  nu_eps <- as.numeric(control$nu_eps %||% 1e-6)
  c_min <- as.numeric(control$c_min %||% 1e-8)
  c_max <- as.numeric(control$c_max %||% 1e6)
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  optim_control <- control$optim_control %||% list(maxit = 500L, reltol = 1e-10)
  base_start <- control$start_theta %||%
    sunspots_time_varying_initial_theta(x = x, u = u, nu_eps = nu_eps)
  if (identical(hemisphere_regression, "shared")) {
    base_start <- sunspots_time_varying_shared_theta(
      base_start, nu_eps = nu_eps, c_min = c_min, c_max = c_max
    )
  }
  base_start <- sunspots_time_varying_validate_theta(
    base_start, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )

  start_thetas <- list(
    base_start,
    modifyList(base_start, list(c = 10)),
    modifyList(base_start, list(c = 30)),
    modifyList(base_start, list(c = 60))
  )
  start_thetas <- lapply(start_thetas, sunspots_time_varying_validate_theta,
                         nu_eps = nu_eps, c_min = c_min, c_max = c_max,
                         hemisphere_regression = hemisphere_regression)

  objective <- function(par) {
    theta <- sunspots_time_varying_unpack_par(
      par, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
      hemisphere_regression = hemisphere_regression
    )
    value <- -sunspots_time_varying_loglik(x = x, u = u, theta = theta)
    if (is.finite(value)) value else .Machine$double.xmax / 100
  }

  fits <- lapply(start_thetas, function(start_theta) {
    par0 <- sunspots_time_varying_pack_theta(
      start_theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
      hemisphere_regression = hemisphere_regression
    )
    try(stats::optim(par = par0, fn = objective, method = "BFGS", control = optim_control), silent = TRUE)
  })
  fits <- Filter(function(fit) !inherits(fit, "try-error") && is.finite(fit$value), fits)
  if (length(fits) == 0L) stop("All time-varying mixture optimizations failed.")

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "value"))]]
  theta_hat <- sunspots_time_varying_unpack_par(
    best$par, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
  c(theta_hat, list(
    loglik = -best$value,
    opt = best,
    n_starts = length(start_thetas),
    n_successful_starts = length(fits),
    start_theta = base_start
  ))
}
sunspots_joint_time_log_density <- function(s, eta, control = list()) {
  s <- as.numeric(s)
  if (any(!is.finite(s)) || any(s <= 0 | s >= 1)) stop("`s` must lie in (0, 1).")
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  log1 <- stats::dbeta(s, eta$alpha1, eta$beta1, log = TRUE)
  log2 <- stats::dbeta(s, eta$alpha2, eta$beta2, log = TRUE)
  rotational_logsumexp2(log(eta$weight1) + log1, log1p(-eta$weight1) + log2)
}
sunspots_time_varying_loglik <- function(x, u, theta) {
  sum(sunspots_time_varying_log_density(x = x, u = u, theta = theta))
}
sunspots_joint_time_pack_eta <- function(eta, control = list()) {
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  c(stats::qlogis(eta$weight1), log(eta$alpha1), log(eta$beta1), log(eta$alpha2), log(eta$beta2))
}
sunspots_time_varying_pack_theta <- function(theta,
                                             nu_eps = 1e-6,
                                             c_min = 1e-8,
                                             c_max = 1e6,
                                             hemisphere_regression = "asymmetric") {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  theta <- sunspots_time_varying_validate_theta(
    theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
  encode_path <- function(start, end) {
    start_prob <- (start - nu_eps) / (1 - 2 * nu_eps)
    end_prob <- (end - nu_eps) / (start - nu_eps)
    c(stats::qlogis(start_prob), stats::qlogis(end_prob))
  }
  if (identical(hemisphere_regression, "shared")) {
    return(c(
      encode_path(theta$a_N, theta$nu_N_end),
      sunspots_time_varying_inverse_softplus(theta$c)
    ))
  }
  c(encode_path(theta$a_N, theta$nu_N_end), encode_path(theta$a_S, theta$nu_S_end),
    sunspots_time_varying_inverse_softplus(theta$c))
}
sunspots_joint_state_from_par <- function(par, control = list()) {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(
    control$hemisphere_regression %||% "asymmetric"
  )
  spatial_length <- if (identical(hemisphere_regression, "shared")) 3L else 5L
  par <- as.numeric(par)
  if (length(par) != 5L + spatial_length) stop("Joint parameter vector has an incompatible length.")
  list(
    eta = sunspots_joint_time_unpack_eta(par[seq_len(5L)], control),
    theta = sunspots_time_varying_unpack_par(
      par[5L + seq_len(spatial_length)],
      nu_eps = as.numeric(control$nu_eps %||% 1e-6),
      c_min = as.numeric(control$c_min %||% 1e-8),
      c_max = as.numeric(control$c_max %||% 1e6),
      hemisphere_regression = hemisphere_regression
    ),
    eta_par = par[seq_len(5L)],
    theta_par = par[5L + seq_len(spatial_length)]
  )
}
sunspots_joint_time_score_matrix <- function(s, par, control = list()) {
  s <- as.numeric(s)
  eta <- sunspots_joint_time_unpack_eta(par, control)
  log1 <- stats::dbeta(s, eta$alpha1, eta$beta1, log = TRUE)
  log2 <- stats::dbeta(s, eta$alpha2, eta$beta2, log = TRUE)
  log_mix <- rotational_logsumexp2(log(eta$weight1) + log1, log1p(-eta$weight1) + log2)
  responsibility1 <- exp(log(eta$weight1) + log1 - log_mix)
  responsibility2 <- 1 - responsibility1
  cbind(
    responsibility1 - eta$weight1,
    eta$alpha1 * responsibility1 * (log(s) - digamma(eta$alpha1) + digamma(eta$alpha1 + eta$beta1)),
    eta$beta1 * responsibility1 * (log1p(-s) - digamma(eta$beta1) + digamma(eta$alpha1 + eta$beta1)),
    eta$alpha2 * responsibility2 * (log(s) - digamma(eta$alpha2) + digamma(eta$alpha2 + eta$beta2)),
    eta$beta2 * responsibility2 * (log1p(-s) - digamma(eta$beta2) + digamma(eta$alpha2 + eta$beta2))
  )
}
sunspots_time_gof_score_matrix <- function(data, par, control = list()) {
  data <- normalize_sunspots_time_gof_data(data$x, data$u)
  state <- sunspots_time_gof_state_from_par(par, control)
  theta <- state$theta
  north_path <- sunspots_time_gof_path_state(par[[1L]], par[[2L]], data$u, state$nu_eps)
  south_path <- if (identical(state$hemisphere_regression, "shared")) {
    north_path
  } else {
    sunspots_time_gof_path_state(par[[3L]], par[[4L]], data$u, state$nu_eps)
  }
  z <- pmin(pmax(data$x[, 3L], -1), 1)

  log_north <- sunspots_time_varying_component_log_axis_density(z, theta$c, north_path$nu)
  log_south <- sunspots_time_varying_component_log_axis_density(-z, theta$c, south_path$nu)
  log_mix <- rotational_logsumexp2(log(0.5) + log_north, log(0.5) + log_south)
  responsibility_north <- exp(log(0.5) + log_north - log_mix)
  responsibility_south <- 1 - responsibility_north

  score_north <- fast_multiplier_small_circle_component_scores_natural(
    s = z, kappa = theta$c, nu = north_path$nu
  )
  score_south <- fast_multiplier_small_circle_component_scores_natural(
    s = -z, kappa = theta$c, nu = south_path$nu
  )
  c_derivative <- stats::plogis(par[[if (identical(state$hemisphere_regression, "shared")) 3L else 5L]])

  score_c <- (responsibility_north * score_north[, 1L] +
    responsibility_south * score_south[, 1L]) * c_derivative
  if (identical(state$hemisphere_regression, "shared")) {
    return(cbind(
      (responsibility_north * score_north[, 2L] + responsibility_south * score_south[, 2L]) *
        north_path$d_start_raw,
      (responsibility_north * score_north[, 2L] + responsibility_south * score_south[, 2L]) *
        north_path$d_end_raw,
      score_c
    ))
  }
  cbind(
    responsibility_north * score_north[, 2L] * north_path$d_start_raw,
    responsibility_north * score_north[, 2L] * north_path$d_end_raw,
    responsibility_south * score_south[, 2L] * south_path$d_start_raw,
    responsibility_south * score_south[, 2L] * south_path$d_end_raw,
    score_c
  )
}
sample_sunspots_joint_time_beta_mixture2 <- function(n, eta, control = list()) {
  n <- as.integer(n)
  if (!is.finite(n) || n <= 0L) stop("`n` must be a positive integer.")
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  component1 <- stats::runif(n) <= eta$weight1
  out <- numeric(n)
  out[component1] <- stats::rbeta(sum(component1), eta$alpha1, eta$beta1)
  out[!component1] <- stats::rbeta(sum(!component1), eta$alpha2, eta$beta2)
  rotational_clamp_unit_interval(out, eps = .Machine$double.eps)
}
sunspots_time_varying_nu <- function(u, theta) {
  u <- as.numeric(u)
  if (any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) {
    stop("`u` must contain finite values in (0, 1).")
  }
  theta <- sunspots_time_varying_validate_theta(theta)
  list(
    north = theta$a_N + theta$b_N * u,
    south = theta$a_S + theta$b_S * u
  )
}
sunspots_joint_time_canonicalize_eta <- function(eta, control = list()) {
  eta <- sunspots_joint_time_validate_eta(eta, control)
  mean1 <- eta$alpha1 / (eta$alpha1 + eta$beta1)
  mean2 <- eta$alpha2 / (eta$alpha2 + eta$beta2)
  if (mean1 <= mean2) {
    return(c(eta, list(mean1 = mean1, mean2 = mean2, component_swapped = FALSE)))
  }
  out <- list(
    weight1 = 1 - eta$weight1,
    alpha1 = eta$alpha2,
    beta1 = eta$beta2,
    alpha2 = eta$alpha1,
    beta2 = eta$beta1,
    n_parameters = 5L
  )
  c(out, list(
    mean1 = out$alpha1 / (out$alpha1 + out$beta1),
    mean2 = out$alpha2 / (out$alpha2 + out$beta2),
    component_swapped = TRUE
  ))
}
sunspots_time_varying_validate_theta <- function(theta,
                                                  nu_eps = 1e-6,
                                                  c_min = 1e-8,
                                                  c_max = 1e6,
                                                  hemisphere_regression = "asymmetric") {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  if (identical(hemisphere_regression, "shared") &&
      all(c("a", "b", "c") %in% names(theta)) &&
      !all(c("a_N", "b_N", "a_S", "b_S") %in% names(theta))) {
    theta <- c(list(a_N = theta$a, b_N = theta$b, a_S = theta$a, b_S = theta$b), theta["c"])
  }
  required <- c("a_N", "b_N", "a_S", "b_S", "c")
  if (!is.list(theta) || !all(required %in% names(theta))) {
    stop("`theta` must contain a_N, b_N, a_S, b_S, and c.")
  }

  values <- vapply(theta[required], as.numeric, numeric(1L))
  if (any(!is.finite(values))) {
    stop("All time-varying mixture parameters must be finite scalars.")
  }
  if (!is.finite(nu_eps) || nu_eps <= 0 || nu_eps >= 0.5) {
    stop("`nu_eps` must lie in (0, 0.5).")
  }
  if (!is.finite(c_min) || !is.finite(c_max) || c_min <= 0 || c_max <= c_min) {
    stop("`c_min` and `c_max` must satisfy 0 < c_min < c_max.")
  }

  nu_N_end <- values[["a_N"]] + values[["b_N"]]
  nu_S_end <- values[["a_S"]] + values[["b_S"]]
  if (values[["b_N"]] >= 0 || values[["b_S"]] >= 0) {
    stop("Both latitude regressions must be strictly decreasing: b_N < 0 and b_S < 0.")
  }
  if (identical(hemisphere_regression, "shared") &&
      (!isTRUE(all.equal(values[["a_N"]], values[["a_S"]], tolerance = 0)) ||
       !isTRUE(all.equal(values[["b_N"]], values[["b_S"]], tolerance = 0)))) {
    stop("The shared hemisphere regression requires a_N = a_S and b_N = b_S.")
  }
  if (any(c(values[["a_N"]], nu_N_end, values[["a_S"]], nu_S_end) <= nu_eps) ||
      any(c(values[["a_N"]], nu_N_end, values[["a_S"]], nu_S_end) >= 1 - nu_eps)) {
    stop("Both linear latitude paths must remain in (nu_eps, 1 - nu_eps) on [0, 1].")
  }
  if (values[["c"]] < c_min || values[["c"]] > c_max) {
    stop("`c` is outside the configured admissible interval.")
  }

  c(values, list(
    nu_N_end = nu_N_end,
    nu_S_end = nu_S_end,
    hemisphere_regression = hemisphere_regression,
    n_parameters = if (identical(hemisphere_regression, "shared")) 3L else 5L
  ))
}
sunspots_joint_profile_block_r <- function(radii, rho, center_s, time_nodes,
                                            time_weights, coefficients) {
  radii <- as.matrix(radii)
  rho <- pmin(pmax(as.numeric(rho), -1), 1)
  center_s <- as.numeric(center_s)
  time_nodes <- as.numeric(time_nodes)
  time_weights <- as.numeric(time_weights)
  coefficients <- as.matrix(coefficients)
  if (length(rho) != nrow(radii) || length(center_s) != nrow(radii) ||
      length(time_nodes) != length(time_weights) || nrow(coefficients) != length(time_nodes)) {
    stop("Joint profile block inputs have incompatible dimensions.")
  }
  out <- matrix(0, nrow = nrow(radii), ncol = ncol(radii))
  for (q in seq_along(time_nodes)) {
    sphere_radius <- pi * pmin(pmax(2 * radii - abs(center_s - time_nodes[[q]]), 0), 1)
    cdf <- small_circle_projection_cdf_legendre_matrix(
      x_matrix = cos(sphere_radius), r = rho, coefficients = coefficients[q, ], enforce_bounds = TRUE
    )
    out <- out + time_weights[[q]] * (1 - cdf)
  }
  pmin(pmax(out, 0), 1)
}
sunspots_joint_normalize_backend <- function(backend = c("auto", "r", "cpp")) {
  backend <- tolower(as.character(backend))
  if (length(backend) > 1L) backend <- backend[[1L]]
  if (length(backend) != 1L || is.na(backend) || !backend %in% c("auto", "r", "cpp")) {
    stop("`distance_profile_backend` must be one of 'auto', 'r', or 'cpp'.")
  }
  backend
}
sunspots_joint_time_boundary_flags <- function(eta, control = list(), tol = 1e-5) {
  diagnostics <- sunspots_joint_time_parameter_boundaries(eta, control, tol = tol)
  list(
    weight = diagnostics$near_any_bound[diagnostics$parameter == "weight1"],
    shape_lower = any(diagnostics$near_lower_bound[diagnostics$parameter != "weight1"]),
    shape_upper = any(diagnostics$near_upper_bound[diagnostics$parameter != "weight1"])
  )
}
sunspots_joint_time_control <- function(control = list()) {
  shape_lower <- as.numeric(control$time_beta_shape_lower %||% 1e-6)
  shape_upper <- as.numeric(control$time_beta_shape_upper %||% 1e3)
  weight_eps <- as.numeric(control$time_beta_weight_eps %||% 0.01)

  if (!is.finite(shape_lower) ||
      shape_lower <= 0 ||
      !is.finite(shape_upper) ||
      shape_upper <= shape_lower) {
    stop("Time-beta shape bounds must satisfy 0 < lower < upper.")
  }
  if (!is.finite(weight_eps) || weight_eps <= 0 || weight_eps >= 0.5) {
    stop("`time_beta_weight_eps` must lie in (0, 0.5).")
  }

  list(
    shape_lower = shape_lower,
    shape_upper = shape_upper,
    weight_eps = weight_eps
  )
}
sunspots_joint_time_parameter_boundaries <- function(eta, control = list(), tol = 1e-5) {
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  bounds <- sunspots_joint_time_control(control)
  tol <- as.numeric(tol)
  if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be one nonnegative finite value.")
  }
  parameter <- c("weight1", "alpha1", "beta1", "alpha2", "beta2")
  value <- unlist(eta[parameter], use.names = FALSE)
  lower <- c(bounds$weight_eps, rep(bounds$shape_lower, 4L))
  upper <- c(1 - bounds$weight_eps, rep(bounds$shape_upper, 4L))
  distance_lower <- value - lower
  distance_upper <- upper - value
  data.frame(
    parameter = parameter,
    value = value,
    lower_bound = lower,
    upper_bound = upper,
    distance_to_lower = distance_lower,
    distance_to_upper = distance_upper,
    near_lower_bound = distance_lower <= tol,
    near_upper_bound = distance_upper <= tol,
    near_any_bound = distance_lower <= tol | distance_upper <= tol,
    stringsAsFactors = FALSE
  )
}
sunspots_joint_time_select_fit <- function(fits) {
  finite_fits <- Filter(function(fit) {
    !inherits(fit, "try-error") &&
      is.list(fit) &&
      length(fit$value) == 1L &&
      is.finite(fit$value)
  }, fits)

  is_converged <- function(fit) {
    convergence <- as.integer(fit$convergence %||% NA_integer_)
    length(convergence) == 1L &&
      !is.na(convergence) &&
      convergence == 0L
  }
  converged_fits <- Filter(is_converged, finite_fits)

  if (length(converged_fits) == 0L) {
    stop(
      paste(
        "No temporal two-beta-mixture optimization converged",
        "with `optim$convergence == 0`."
      ),
      call. = FALSE
    )
  }

  best <- converged_fits[[
    which.min(vapply(converged_fits, `[[`, numeric(1L), "value"))
  ]]

  list(
    fit = best,
    n_finite_fits = length(finite_fits),
    n_converged_fits = length(converged_fits),
    selected_converged = TRUE
  )
}
sunspots_joint_time_start_etas <- function(s, control = list()) {
  s <- as.numeric(s)
  split_probs <- c(0.2, 0.35, 0.5, 0.65, 0.8)
  candidates <- list()
  for (split_prob in split_probs) {
    cutoff <- stats::quantile(s, probs = split_prob, names = FALSE)
    group1 <- s <= cutoff
    if (all(group1) || !any(group1)) next
    comp1 <- sunspots_joint_time_moment_match(s[group1], control)
    comp2 <- sunspots_joint_time_moment_match(s[!group1], control)
    candidates[[length(candidates) + 1L]] <- sunspots_joint_time_canonicalize_eta(list(
      weight1 = mean(group1), alpha1 = comp1$alpha, beta1 = comp1$beta,
      alpha2 = comp2$alpha, beta2 = comp2$beta
    ), control)
  }
  if (length(candidates) == 0L) {
    component <- sunspots_joint_time_moment_match(s, control)
    candidates[[1L]] <- sunspots_joint_time_canonicalize_eta(list(
      weight1 = 0.5, alpha1 = component$alpha, beta1 = component$beta,
      alpha2 = component$alpha, beta2 = component$beta
    ), control)
  }
  candidates
}
sunspots_joint_time_unpack_eta <- function(par, control = list()) {
  par <- as.numeric(par)
  if (length(par) != 5L || any(!is.finite(par))) stop("The time-beta parameter vector must have five finite entries.")
  bounds <- sunspots_joint_time_control(control)
  sunspots_joint_time_canonicalize_eta(list(
    weight1 = rotational_bounded_weight(par[[1L]], weight_eps = bounds$weight_eps),
    alpha1 = rotational_positive_parameter(par[[2L]], lower = bounds$shape_lower, upper = bounds$shape_upper),
    beta1 = rotational_positive_parameter(par[[3L]], lower = bounds$shape_lower, upper = bounds$shape_upper),
    alpha2 = rotational_positive_parameter(par[[4L]], lower = bounds$shape_lower, upper = bounds$shape_upper),
    beta2 = rotational_positive_parameter(par[[5L]], lower = bounds$shape_lower, upper = bounds$shape_upper)
  ), control)
}
sunspots_time_varying_initial_theta <- function(x, u, nu_eps = 1e-6) {
  z <- pmin(pmax(as.numeric(x[, 3L]), -1), 1)
  early <- u <= stats::quantile(u, probs = 1 / 3, names = FALSE)
  late <- u >= stats::quantile(u, probs = 2 / 3, names = FALSE)
  path_start <- function(values, index_early, index_late) {
    fallback <- stats::median(values, na.rm = TRUE)
    start <- stats::median(values[index_early], na.rm = TRUE)
    end <- stats::median(values[index_late], na.rm = TRUE)
    if (!is.finite(start)) start <- fallback
    if (!is.finite(end)) end <- fallback * 0.8
    start <- min(max(start, 0.08), 0.85)
    end <- min(max(end, nu_eps * 2), start - 0.01)
    if (end <= nu_eps) end <- max(nu_eps * 2, start * 0.8)
    c(a = start, b = end - start)
  }

  north_values <- z[z >= 0]
  south_values <- -z[z < 0]
  north <- path_start(z[z >= 0], early[z >= 0], late[z >= 0])
  south <- path_start(-z[z < 0], early[z < 0], late[z < 0])
  fitted_mean <- ifelse(z >= 0, north[["a"]] + north[["b"]] * u,
                        -(south[["a"]] + south[["b"]] * u))
  variance <- stats::var(z - fitted_mean)
  c_value <- min(max(1 / (2 * max(variance, 0.01)), 1), 200)

  list(a_N = north[["a"]], b_N = north[["b"]],
       a_S = south[["a"]], b_S = south[["b"]], c = c_value)
}
sunspots_time_varying_shared_theta <- function(theta,
                                                nu_eps = 1e-6,
                                                c_min = 1e-8,
                                                c_max = 1e6) {
  theta <- sunspots_time_varying_validate_theta(
    theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max
  )
  list(
    a_N = mean(c(theta$a_N, theta$a_S)),
    b_N = mean(c(theta$b_N, theta$b_S)),
    a_S = mean(c(theta$a_N, theta$a_S)),
    b_S = mean(c(theta$b_N, theta$b_S)),
    c = theta$c
  )
}
sunspots_time_varying_unpack_par <- function(par,
                                             nu_eps = 1e-6,
                                             c_min = 1e-8,
                                             c_max = 1e6,
                                             hemisphere_regression = "asymmetric") {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  par <- as.numeric(par)
  expected_length <- if (identical(hemisphere_regression, "shared")) 3L else 5L
  if (length(par) != expected_length || any(!is.finite(par))) {
    stop(sprintf("The unconstrained parameter vector must contain %d finite values.", expected_length))
  }

  nu_upper <- 1 - nu_eps
  decode_path <- function(start_raw, end_raw) {
    # BFGS can evaluate very large finite values, for which plogis rounds to 0 or 1.
    prob_eps <- sqrt(.Machine$double.eps)
    start_prob <- pmin(pmax(stats::plogis(start_raw), prob_eps), 1 - prob_eps)
    end_prob <- pmin(pmax(stats::plogis(end_raw), prob_eps), 1 - prob_eps)
    start <- nu_eps + (nu_upper - nu_eps) * start_prob
    end <- nu_eps + (start - nu_eps) * end_prob
    c(a = start, b = end - start)
  }

  north <- decode_path(par[[1L]], par[[2L]])
  south <- if (identical(hemisphere_regression, "shared")) north else decode_path(par[[3L]], par[[4L]])
  c_index <- if (identical(hemisphere_regression, "shared")) 3L else 5L
  theta <- list(
    a_N = north[["a"]],
    b_N = north[["b"]],
    a_S = south[["a"]],
    b_S = south[["b"]],
    c = min(max(sunspots_time_varying_softplus(par[[c_index]]), c_min), c_max)
  )
  sunspots_time_varying_validate_theta(
    theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
}
sunspots_time_varying_log_density <- function(x, u, theta) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  u <- as.numeric(u)
  if (length(u) != nrow(x)) {
    stop("`u` must have one entry per row of `x`.")
  }
  theta <- sunspots_time_varying_validate_theta(theta)
  nu <- sunspots_time_varying_nu(u, theta)
  z <- pmin(pmax(x[, 3L], -1), 1)

  log_north <- sunspots_time_varying_component_log_axis_density(z, theta$c, nu$north)
  log_south <- sunspots_time_varying_component_log_axis_density(-z, theta$c, nu$south)
  -log(2 * pi) + rotational_logsumexp2(log(0.5) + log_north, log(0.5) + log_south)
}
sunspots_time_varying_inverse_softplus <- function(x) {
  x <- as.numeric(x)
  ifelse(x > 30, x, log(expm1(x)))
}
normalize_sunspots_time_gof_data <- function(x, u) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  if (ncol(x) != 3L) stop("The conditional sunspots GOF supports only S^2 data.")
  u <- as.numeric(u)
  if (length(u) != nrow(x) || any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) {
    stop("`u` must have one finite value in (0, 1) for each row of `x`.")
  }
  list(x = x, u = u)
}
sunspots_time_gof_path_state <- function(start_raw, end_raw, u, nu_eps) {
  prob_eps <- sqrt(.Machine$double.eps)
  p_start <- pmin(pmax(stats::plogis(start_raw), prob_eps), 1 - prob_eps)
  p_end <- pmin(pmax(stats::plogis(end_raw), prob_eps), 1 - prob_eps)
  range_nu <- 1 - 2 * nu_eps
  start <- nu_eps + range_nu * p_start
  end <- nu_eps + (start - nu_eps) * p_end
  d_start <- range_nu * p_start * (1 - p_start)
  d_end_start <- p_end * d_start
  d_end_end <- (start - nu_eps) * p_end * (1 - p_end)

  list(
    nu = (1 - u) * start + u * end,
    d_start_raw = (1 - u) * d_start + u * d_end_start,
    d_end_raw = u * d_end_end
  )
}
sunspots_time_gof_state_from_par <- function(par, control = list()) {
  nu_eps <- as.numeric(control$nu_eps %||% 1e-6)
  c_min <- as.numeric(control$c_min %||% 1e-8)
  c_max <- as.numeric(control$c_max %||% 1e6)
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(
    control$hemisphere_regression %||% "asymmetric"
  )
  theta <- sunspots_time_varying_unpack_par(
    par, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
  list(
    theta = theta,
    nu_eps = nu_eps,
    c_min = c_min,
    c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
}
sunspots_time_varying_component_log_axis_density <- function(z, c_value, nu) {
  z <- as.numeric(z)
  nu <- as.numeric(nu)
  if (length(z) != length(nu)) {
    stop("`z` and `nu` must have the same length.")
  }
  log_norm <- small_circle_log_norm_constant(kappa = c_value, nu = nu)
  -log(2) - log_norm - c_value * (z - nu)^2
}
sunspots_joint_time_validate_eta <- function(eta, control = list()) {
  bounds <- sunspots_joint_time_control(control)
  required <- c("weight1", "alpha1", "beta1", "alpha2", "beta2")
  if (!is.list(eta) || !all(required %in% names(eta))) {
    stop("`eta` must contain weight1, alpha1, beta1, alpha2, and beta2.")
  }
  values <- vapply(eta[required], as.numeric, numeric(1L))
  if (any(!is.finite(values))) stop("All time-beta parameters must be finite scalars.")
  if (values[["weight1"]] < bounds$weight_eps || values[["weight1"]] > 1 - bounds$weight_eps) {
    stop("`weight1` is outside its admissible interval.")
  }
  shapes <- values[c("alpha1", "beta1", "alpha2", "beta2")]
  if (any(shapes < bounds$shape_lower) || any(shapes > bounds$shape_upper)) {
    stop("Time-beta shapes are outside their admissible interval.")
  }
  c(as.list(values), list(n_parameters = 5L))
}
sunspots_joint_time_moment_match <- function(s, control = list()) {
  bounds <- sunspots_joint_time_control(control)
  s <- rotational_clamp_unit_interval(s, eps = 1e-10)
  mean_s <- min(max(mean(s), 1e-5), 1 - 1e-5)
  variance_s <- stats::var(s)
  max_variance <- mean_s * (1 - mean_s) * (1 - 1e-6)
  variance_s <- min(max(variance_s, 1e-6), max_variance)
  precision <- min(max(mean_s * (1 - mean_s) / variance_s - 1, 2 * bounds$shape_lower),
                   bounds$shape_upper / max(mean_s, 1 - mean_s))
  list(
    alpha = min(max(mean_s * precision, bounds$shape_lower), bounds$shape_upper),
    beta = min(max((1 - mean_s) * precision, bounds$shape_lower), bounds$shape_upper)
  )
}
sunspots_time_varying_softplus <- function(x) {
  x <- as.numeric(x)
  pmax(x, 0) + log1p(exp(-abs(x)))
}
