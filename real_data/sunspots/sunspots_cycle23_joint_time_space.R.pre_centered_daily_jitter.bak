#!/usr/bin/env Rscript

# Joint time--space model for the first recorded observations of cycle-23
# sunspot groups.  The legacy rank-based conditional analysis is deliberately
# kept separate in run_sunspots_cycle23_time_varying_asymmetric_mixture_gof.R.

resolve_sunspots_joint_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

sunspots_joint_time_model_path <- resolve_sunspots_joint_path(
  "real_data", "sunspots", "run_sunspots_cycle23_time_varying_asymmetric_mixture.R"
)
sunspots_joint_gof_path <- resolve_sunspots_joint_path(
  "real_data", "sunspots", "run_sunspots_cycle23_time_varying_asymmetric_mixture_gof.R"
)
source(sunspots_joint_time_model_path)
source(sunspots_joint_gof_path)

sunspots_joint_with_seed <- function(seed, expr) {
  seed <- as.integer(seed)
  if (length(seed) != 1L || !is.finite(seed)) stop("`seed` must be one finite integer.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  eval.parent(substitute(expr))
}

sunspots_joint_normalize_backend <- function(backend = c("auto", "r", "cpp")) {
  backend <- tolower(as.character(backend))
  if (length(backend) > 1L) backend <- backend[[1L]]
  if (length(backend) != 1L || is.na(backend) || !backend %in% c("auto", "r", "cpp")) {
    stop("`distance_profile_backend` must be one of 'auto', 'r', or 'cpp'.")
  }
  backend
}

sunspots_joint_effective_backend <- function(backend = "auto") {
  backend <- sunspots_joint_normalize_backend(backend)
  if (!identical(backend, "auto")) return(backend)
  if (!requireNamespace("Rcpp", quietly = TRUE)) return("r")
  loaded <- try(ensure_distance_profile_cpp_loaded(), silent = TRUE)
  if (inherits(loaded, "try-error")) "r" else "cpp"
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

prepare_sunspots_cycle23_joint_time_space_data <- function(
    input_csv,
    start_date = "1997-06-01",
    end_date = "2006-01-01",
    dequantization_seed = 20260711L) {
  if (!file.exists(input_csv)) source(prep_path_sunspots_time_varying)
  if (!file.exists(input_csv)) stop(sprintf("Input CSV not found: %s", input_csv))

  df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  required <- c("cycle", "date", "x1", "x2", "x3")
  if (!all(required %in% names(df))) {
    stop(sprintf("Input CSV is missing: %s", paste(setdiff(required, names(df)), collapse = ", ")))
  }
  if (!all(df$cycle == 23L)) stop("The joint runner is restricted to cycle 23 data.")

  recorded_timestamp <- as.POSIXct(df$date, tz = "UTC")
  calendar_day <- as.Date(recorded_timestamp, tz = "UTC")
  start_time <- as.POSIXct(as.Date(start_date), tz = "UTC")
  end_time <- as.POSIXct(as.Date(end_date), tz = "UTC")
  if (is.na(start_time) || is.na(end_time) || start_time >= end_time) {
    stop("`start_date` and `end_date` must define a non-empty chronological interval.")
  }
  day_start <- as.POSIXct(calendar_day, tz = "UTC")
  keep <- !is.na(recorded_timestamp) & day_start >= start_time & day_start < end_time
  retained <- df[keep, , drop = FALSE]
  retained$recorded_timestamp <- recorded_timestamp[keep]
  retained$calendar_day <- calendar_day[keep]
  retained$calendar_day_start <- day_start[keep]
  retained <- retained[order(retained$calendar_day_start, retained$NOAA), , drop = FALSE]
  if (nrow(retained) == 0L) stop("The requested date interval retains no observations.")

  jitter <- sunspots_joint_with_seed(dequantization_seed, stats::runif(nrow(retained)))
  jitter <- pmin(pmax(jitter, .Machine$double.eps), 1 - .Machine$double.eps)
  dequantized_seconds <- as.numeric(retained$calendar_day_start) + jitter * 86400
  span_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  s <- (dequantized_seconds - as.numeric(start_time)) / span_seconds
  if (any(!is.finite(s)) || any(s <= 0 | s >= 1)) {
    stop("Dequantized times must lie strictly inside the analysis interval.")
  }

  retained$dequantization_jitter_day <- jitter
  retained$dequantization_jitter_centered_day <- jitter - 0.5
  retained$dequantized_timestamp <- as.POSIXct(dequantized_seconds, origin = "1970-01-01", tz = "UTC")
  retained$s <- s
  retained$start_date <- start_date
  retained$end_date_exclusive <- end_date
  retained
}

# Short alias for temporal-only diagnostics; both names use the same loader.
prepare_sunspots_joint_time_data <- prepare_sunspots_cycle23_joint_time_space_data

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

sunspots_joint_time_pack_eta <- function(eta, control = list()) {
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  c(stats::qlogis(eta$weight1), log(eta$alpha1), log(eta$beta1), log(eta$alpha2), log(eta$beta2))
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

sunspots_joint_time_log_density <- function(s, eta, control = list()) {
  s <- as.numeric(s)
  if (any(!is.finite(s)) || any(s <= 0 | s >= 1)) stop("`s` must lie in (0, 1).")
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  log1 <- stats::dbeta(s, eta$alpha1, eta$beta1, log = TRUE)
  log2 <- stats::dbeta(s, eta$alpha2, eta$beta2, log = TRUE)
  rotational_logsumexp2(log(eta$weight1) + log1, log1p(-eta$weight1) + log2)
}

sunspots_joint_time_density <- function(s, eta, control = list()) {
  exp(sunspots_joint_time_log_density(s, eta, control))
}

sunspots_joint_time_cdf <- function(s, eta, control = list()) {
  eta <- sunspots_joint_time_canonicalize_eta(eta, control)
  s <- as.numeric(s)
  out <- numeric(length(s))
  out[s >= 1] <- 1
  active <- s > 0 & s < 1
  out[active] <- eta$weight1 * stats::pbeta(s[active], eta$alpha1, eta$beta1) +
    (1 - eta$weight1) * stats::pbeta(s[active], eta$alpha2, eta$beta2)
  pmin(pmax(out, 0), 1)
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

sunspots_joint_time_boundary_flags <- function(eta, control = list(), tol = 1e-5) {
  diagnostics <- sunspots_joint_time_parameter_boundaries(eta, control, tol = tol)
  list(
    weight = diagnostics$near_any_bound[diagnostics$parameter == "weight1"],
    shape_lower = any(diagnostics$near_lower_bound[diagnostics$parameter != "weight1"]),
    shape_upper = any(diagnostics$near_upper_bound[diagnostics$parameter != "weight1"])
  )
}

sunspots_joint_time_information_criteria <- function(loglik, n, n_parameters = 5L) {
  loglik <- as.numeric(loglik)
  n <- as.integer(n)
  n_parameters <- as.integer(n_parameters)
  if (length(loglik) != 1L || !is.finite(loglik) ||
      length(n) != 1L || !is.finite(n) || n < 1L ||
      length(n_parameters) != 1L || !is.finite(n_parameters) || n_parameters < 1L) {
    stop("`loglik`, `n`, and `n_parameters` must be finite scalar values with positive integer counts.")
  }
  list(
    aic = 2 * n_parameters - 2 * loglik,
    bic = n_parameters * log(n) - 2 * loglik
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

sunspots_joint_distance <- function(x, s, omega, center_s) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  s <- as.numeric(s)
  center_s <- as.numeric(center_s)
  if (any(!is.finite(s)) || any(s < 0 | s > 1) ||
      length(center_s) != 1L || !is.finite(center_s) || center_s < 0 || center_s > 1) {
    stop("`s` and `center_s` must lie in [0, 1].")
  }
  0.5 * (acos(pmin(pmax(as.numeric(x %*% omega), -1), 1)) / pi + abs(s - center_s))
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

sunspots_joint_select_centers <- function(n, n_sample_centers = 100L, seed = 20260711L) {
  n <- as.integer(n)
  if (!is.finite(n) || n <= 0L) stop("`n` must be a positive integer.")
  # TODO: all centers retain O(n^2) distance/profile/correction state; keep this
  # option for small validation jobs until a streaming implementation is added.
  if (is.infinite(n_sample_centers) || n_sample_centers >= n) return(seq_len(n))
  sunspots_joint_with_seed(
    seed,
    sunspots_time_gof_select_centers(
      n,
      n_sample_centers = n_sample_centers,
      seed = seed
    )
  )
}

sunspots_joint_prepare_centers <- function(data, fit, center_indices,
                                            time_quad_n = 64L, l_max = 100L,
                                            spatial_quad_n = 400L, center_block_size = 8L,
                                            distance_profile_backend = "auto") {
  data <- sunspots_joint_validate_data(data$x, data$s)
  center_indices <- as.integer(center_indices)
  if (length(center_indices) == 0L || any(center_indices < 1L | center_indices > nrow(data$x))) {
    stop("`center_indices` must index at least one observation.")
  }
  quadrature <- sunspots_joint_time_quadrature(fit$eta_hat, n_nodes = time_quad_n)
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat, quadrature$nodes, l_max = l_max, quad_n = spatial_quad_n
  )
  center_block_size <- as.integer(center_block_size)
  if (!is.finite(center_block_size) || center_block_size <= 0L) stop("`center_block_size` must be positive.")
  center_blocks <- split(center_indices, ceiling(seq_along(center_indices) / center_block_size))
  centers <- vector("list", length(center_indices))
  center_position <- 0L
  for (block_indices in center_blocks) {
    omega <- data$x[block_indices, , drop = FALSE]
    center_s <- data$s[block_indices]
    spatial <- acos(pmin(pmax(omega %*% t(data$x), -1), 1)) / pi
    radii <- 0.5 * (spatial + abs(outer(center_s, data$s, FUN = "-")))
    order_matrix <- t(apply(radii, 1L, order))
    if (is.null(dim(order_matrix))) order_matrix <- matrix(order_matrix, nrow = 1L)
    sorted_radii <- matrix(0, nrow = nrow(radii), ncol = ncol(radii))
    for (row in seq_len(nrow(radii))) sorted_radii[row, ] <- radii[row, order_matrix[row, ]]
    theoretical <- sunspots_joint_profile_block(
      radii = sorted_radii, rho = omega[, 3L], center_s = center_s,
      time_nodes = quadrature$nodes, time_weights = quadrature$weights,
      coefficients = coefficients, backend = distance_profile_backend
    )
    for (row in seq_len(nrow(radii))) {
      sorted_distances <- sorted_radii[row, ]
      tie_end_positions <- sorted_tie_end_positions(sorted_distances)
      empirical <- tie_end_positions / nrow(data$x)
      profile <- theoretical[row, ]
      profile <- cummax(profile)
      observed_process <- sqrt(nrow(data$x)) * (empirical - profile)
      center_position <- center_position + 1L
      centers[[center_position]] <- list(
        center_index = block_indices[[row]], omega = omega[row, ], s = center_s[[row]],
        order_index = order_matrix[row, ], sorted_distances = sorted_distances,
        tie_end_positions = tie_end_positions, theoretical = profile,
        ks_observed_statistic = max(abs(observed_process)),
        cvm_observed_sum = sum(observed_process^2)
      )
    }
  }
  n_centers <- length(centers)
  list(
    quadrature = quadrature,
    coefficients = coefficients,
    centers = centers,
    ks_statistic = max(vapply(centers, `[[`, numeric(1L), "ks_observed_statistic")),
    cvm_statistic = sum(vapply(centers, `[[`, numeric(1L), "cvm_observed_sum")) / (n_centers * nrow(data$x))
  )
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

sunspots_joint_prepare_fast_corrections <- function(data, fit, centers,
                                                     derivative_mc_size, seed, control = list()) {
  data <- sunspots_joint_validate_data(data$x, data$s)
  boundary_flags <- fit$eta_hat$boundary_flags
  if (isTRUE(boundary_flags$weight) || isTRUE(boundary_flags$shape_lower) ||
      isTRUE(boundary_flags$shape_upper)) {
    message <- paste0(
      "The temporal beta-mixture MLE is at an admissible boundary; the regular fast multiplier correction is unavailable. ",
      "Set `control$allow_boundary_fast = TRUE` only for exploratory numerical output."
    )
    if (!isTRUE(control$allow_boundary_fast)) stop(message, call. = FALSE)
    warning(message, call. = FALSE)
  }
  par_hat <- sunspots_joint_pack_par(fit, control)
  auxiliary <- sunspots_joint_with_seed(seed, sample_sunspots_joint_time_space(
    derivative_mc_size, par_hat, control
  ))
  score_observed <- sunspots_joint_score_matrix(data, par_hat, control)
  score_auxiliary <- sunspots_joint_score_matrix(auxiliary, par_hat, control)
  vhat <- crossprod(score_auxiliary) / nrow(score_auxiliary)
  diagnostics <- fast_multiplier_vhat_diagnostics(
    S_obs = score_observed, Psi_aux = score_auxiliary, Vhat = vhat, par0 = par_hat
  )
  fast_multiplier_validate_vhat(
    Vhat = vhat, diagnostics = diagnostics, label = "joint sunspots fast multiplier preparation",
    rcond_tol = as.numeric(control$fast_multiplier_vhat_rcond_tol %||% 1e-12)
  )
  auxiliary_scores_solved <- score_auxiliary %*% t(fast_multiplier_solve_vhat(
    vhat, diag(ncol(vhat)), label = "the joint sunspots score correction"
  ))
  correction_rows <- vector("list", length(centers))
  for (i in seq_along(centers)) {
    distances <- sunspots_joint_distance(auxiliary$x, auxiliary$s, centers[[i]]$omega, centers[[i]]$s)
    order_index <- order(distances)
    cumulative_scores <- col_cumsums_base(auxiliary_scores_solved[order_index, , drop = FALSE]) /
      nrow(auxiliary_scores_solved)
    score_basis <- rbind(rep(0, ncol(cumulative_scores)), cumulative_scores)
    selected_counts <- findInterval(centers[[i]]$sorted_distances, distances[order_index])
    centers[[i]]$correction <- score_basis[selected_counts + 1L, , drop = FALSE]
    correction_rows[[i]] <- data.frame(
      center_rank = i,
      center_index = centers[[i]]$center_index,
      correction_all_finite = all(is.finite(centers[[i]]$correction)),
      correction_max_abs = suppressWarnings(max(abs(centers[[i]]$correction))),
      stringsAsFactors = FALSE
    )
  }
  correction_diagnostics <- do.call(rbind, correction_rows)
  positive_eigen <- diagnostics$Vhat_eigenvalues[is.finite(diagnostics$Vhat_eigenvalues) & diagnostics$Vhat_eigenvalues > 0]
  min_positive_eigen <- if (length(positive_eigen) > 0L) min(positive_eigen) else NA_real_
  max_positive_eigen <- if (length(positive_eigen) > 0L) max(positive_eigen) else NA_real_
  list(
    centers = centers, score_observed = score_observed, vhat = vhat,
    diagnostics = diagnostics, derivative_mc_size = nrow(auxiliary$x),
    inversion_method = "solve",
    regularization_added = 0,
    correction_diagnostics = correction_diagnostics,
    correction_all_finite = all(correction_diagnostics$correction_all_finite),
    correction_any_nonfinite = any(!correction_diagnostics$correction_all_finite),
    Vhat_min_positive_eigenvalue = min_positive_eigen,
    Vhat_max_positive_eigenvalue = max_positive_eigen
  )
}

# Reference parametric bootstrap used only for small validation runs.  Unlike
# the fast multiplier approximation, every replicate is refitted from scratch.
sunspots_joint_slow_reestimated_statistics <- function(
    n, par, center_indices, B, seed,
    hemisphere_regression = "asymmetric", statistics = c("ks", "cvm"),
    time_quad_n = 32L, l_max = 60L, spatial_quad_n = 200L,
    center_block_size = 8L, distance_profile_backend = "r", control = list()) {
  n <- as.integer(n)
  B <- as.integer(B)
  statistics <- unique(tolower(as.character(statistics)))
  if (!is.finite(n) || n < 4L || !is.finite(B) || B < 1L) {
    stop("`n` must be at least 4 and `B` must be a positive integer.")
  }
  if (length(statistics) == 0L || any(!statistics %in% c("ks", "cvm"))) {
    stop("`statistics` must contain one or both of 'ks' and 'cvm'.")
  }
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  control <- utils::modifyList(control, list(hemisphere_regression = hemisphere_regression))
  center_indices <- as.integer(center_indices)
  if (length(center_indices) == 0L || any(center_indices < 1L | center_indices > n)) {
    stop("`center_indices` must index the simulated sample.")
  }
  results <- sunspots_joint_with_seed(seed, {
    out <- matrix(NA_real_, nrow = B, ncol = length(statistics),
                  dimnames = list(NULL, statistics))
    temporal_boundaries <- matrix(FALSE, nrow = B, ncol = 3L,
                                  dimnames = list(NULL, c("weight", "shape_lower", "shape_upper")))
    for (replicate_index in seq_len(B)) {
      simulated <- sample_sunspots_joint_time_space(n, par, control)
      fit <- suppressWarnings(fit_sunspots_cycle23_joint_time_space(
        simulated$x, simulated$s, hemisphere_regression = hemisphere_regression, control = control
      ))
      temporal_boundaries[replicate_index, ] <- unlist(fit$eta_hat$boundary_flags, use.names = FALSE)
      prepared <- sunspots_joint_prepare_centers(
        data = simulated, fit = fit, center_indices = center_indices,
        time_quad_n = time_quad_n, l_max = l_max, spatial_quad_n = spatial_quad_n,
        center_block_size = center_block_size, distance_profile_backend = distance_profile_backend
      )
      out[replicate_index, ] <- c(ks = prepared$ks_statistic, cvm = prepared$cvm_statistic)[statistics]
    }
    list(statistics = out, temporal_boundaries = temporal_boundaries)
  })
  output <- as.data.frame(results$statistics)
  attr(output, "temporal_boundary_flags") <- results$temporal_boundaries
  output
}

plot_sunspots_cycle23_joint_time_space_diagnostics <- function(data, fit, output_dir) {
  data <- sunspots_joint_validate_data(data$x, data$s)
  s_grid <- seq(1e-5, 1 - 1e-5, length.out = 501L)
  temporal_path <- file.path(output_dir, "cycle23_joint_time_density_pit.png")
  grDevices::png(temporal_path, width = 1400, height = 700, res = 140)
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(data$s, breaks = 40, freq = FALSE, col = "#9ecae1", border = "white",
       xlab = "Dequantized first-record day", main = "Temporal density")
  lines(s_grid, sunspots_joint_time_density(s_grid, fit$eta_hat), col = "#8b0000", lwd = 2)
  time_pit <- sunspots_joint_time_cdf(data$s, fit$eta_hat)
  qqplot(stats::ppoints(length(time_pit)), sort(time_pit), xlim = c(0, 1), ylim = c(0, 1),
         pch = 16, cex = 0.45, col = "#2b6cb0", xlab = "Uniform quantiles",
         ylab = "Temporal PIT quantiles", main = "Temporal PIT Q-Q plot")
  abline(0, 1, col = "#8b0000", lwd = 2)
  par(old_par)
  grDevices::dev.off()

  latitude_path <- file.path(output_dir, "cycle23_joint_latitude_paths.png")
  latitude <- asin(pmin(pmax(data$x[, 3L], -1), 1)) * 180 / pi
  nu <- sunspots_time_varying_nu(s_grid, fit$theta_hat)
  grDevices::png(latitude_path, width = 1400, height = 900, res = 140)
  plot(data$s, latitude, pch = 16, cex = 0.35,
       col = grDevices::adjustcolor(ifelse(latitude >= 0, "#c43c39", "#2b6cb0"), alpha.f = 0.22),
       xlab = "Dequantized first-record day", ylab = "Latitude (degrees)",
       main = "Cycle 23: fitted north and south small-circle paths")
  lines(s_grid, asin(nu$north) * 180 / pi, col = "#8b0000", lwd = 3)
  lines(s_grid, -asin(nu$south) * 180 / pi, col = "#003f7f", lwd = 3)
  grid(col = "#d9d9d9")
  grDevices::dev.off()

  axial_pit <- sunspots_time_varying_conditional_pit(data$x[, 3L], data$s, fit$theta_hat)
  axial_path <- file.path(output_dir, "cycle23_joint_conditional_axial_pit.png")
  grDevices::png(axial_path, width = 1400, height = 700, res = 140)
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(axial_pit, breaks = 20, freq = FALSE, col = "#9ecae1", border = "white",
       xlab = "Conditional axial PIT", main = "Conditional axial PIT histogram", xlim = c(0, 1))
  abline(h = 1, col = "#8b0000", lwd = 2)
  qqplot(stats::ppoints(length(axial_pit)), sort(axial_pit), xlim = c(0, 1), ylim = c(0, 1),
         pch = 16, cex = 0.45, col = "#2b6cb0", xlab = "Uniform quantiles",
         ylab = "Conditional PIT quantiles", main = "Conditional axial PIT Q-Q plot")
  abline(0, 1, col = "#8b0000", lwd = 2)
  par(old_par)
  grDevices::dev.off()
  list(temporal_path = temporal_path, latitude_path = latitude_path, axial_path = axial_path,
       time_pit = time_pit, axial_pit = axial_pit)
}
