#!/usr/bin/env Rscript

# Parsimonious temporal laws for the joint cycle-23 sunspot model.
#
# This module deliberately leaves the validated two-beta implementation
# untouched. It adds two separate temporal models:
#   1. beta: a single Beta(alpha, beta) law;
#   2. uniform_beta: (1 - lambda) Unif(0, 1) + lambda Beta(alpha, beta).

resolve_sunspots_joint_parsimonious_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

if (!exists("prepare_sunspots_cycle23_joint_time_space_data", mode = "function")) {
  source(resolve_sunspots_joint_parsimonious_path(
    "real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"
  ))
}

normalize_sunspots_joint_parsimonious_time_model <- function(
    time_model = c("beta", "uniform_beta")) {
  match.arg(time_model)
}

sunspots_joint_parsimonious_time_control <- function(control = list()) {
  shape_lower <- as.numeric(
    control$parsimonious_time_shape_lower %||% 1e-3
  )
  shape_upper <- as.numeric(
    control$parsimonious_time_shape_upper %||% 1e3
  )
  weight_eps <- as.numeric(
    control$parsimonious_time_weight_eps %||% 1e-4
  )
  boundary_tol <- as.numeric(
    control$parsimonious_time_boundary_tol %||% 1e-5
  )
  uniform_identifiability_tol <- as.numeric(
    control$parsimonious_uniform_identifiability_tol %||% 1e-3
  )

  if (length(shape_lower) != 1L ||
      !is.finite(shape_lower) ||
      shape_lower <= 0 ||
      length(shape_upper) != 1L ||
      !is.finite(shape_upper) ||
      shape_upper <= shape_lower) {
    stop("Parsimonious time-shape bounds must satisfy 0 < lower < upper.")
  }
  if (length(weight_eps) != 1L ||
      !is.finite(weight_eps) ||
      weight_eps <= 0 ||
      weight_eps >= 0.5) {
    stop("`parsimonious_time_weight_eps` must lie in (0, 0.5).")
  }
  if (length(boundary_tol) != 1L ||
      !is.finite(boundary_tol) ||
      boundary_tol < 0) {
    stop("`parsimonious_time_boundary_tol` must be nonnegative.")
  }
  if (length(uniform_identifiability_tol) != 1L ||
      !is.finite(uniform_identifiability_tol) ||
      uniform_identifiability_tol < 0) {
    stop("`parsimonious_uniform_identifiability_tol` must be nonnegative.")
  }

  list(
    shape_lower = shape_lower,
    shape_upper = shape_upper,
    weight_eps = weight_eps,
    boundary_tol = boundary_tol,
    uniform_identifiability_tol = uniform_identifiability_tol
  )
}

sunspots_joint_parsimonious_time_n_parameters <- function(time_model) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  if (identical(time_model, "beta")) 2L else 3L
}

sunspots_joint_parsimonious_time_validate_eta <- function(
    eta,
    time_model = NULL,
    control = list()) {
  if (!is.list(eta)) stop("`eta` must be a list.")
  if (is.null(time_model)) time_model <- eta$time_model
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  bounds <- sunspots_joint_parsimonious_time_control(control)

  required <- if (identical(time_model, "beta")) {
    c("alpha", "beta")
  } else {
    c("beta_weight", "alpha", "beta")
  }
  if (!all(required %in% names(eta))) {
    stop(sprintf(
      "`eta` for time_model = '%s' must contain: %s.",
      time_model,
      paste(required, collapse = ", ")
    ))
  }

  values <- vapply(eta[required], as.numeric, numeric(1L))
  if (any(!is.finite(values))) {
    stop("All parsimonious temporal parameters must be finite scalars.")
  }
  alpha <- values[["alpha"]]
  beta <- values[["beta"]]
  if (alpha < bounds$shape_lower || alpha > bounds$shape_upper ||
      beta < bounds$shape_lower || beta > bounds$shape_upper) {
    stop("Parsimonious time-shape parameters are outside their bounds.")
  }

  beta_weight <- if (identical(time_model, "beta")) {
    1
  } else {
    values[["beta_weight"]]
  }
  if (identical(time_model, "uniform_beta") &&
      (beta_weight < bounds$weight_eps ||
       beta_weight > 1 - bounds$weight_eps)) {
    stop("`beta_weight` is outside its admissible interval.")
  }

  list(
    time_model = time_model,
    beta_weight = beta_weight,
    uniform_weight = 1 - beta_weight,
    alpha = alpha,
    beta = beta,
    mean = alpha / (alpha + beta),
    n_parameters = sunspots_joint_parsimonious_time_n_parameters(time_model)
  )
}

sunspots_joint_parsimonious_time_pack_eta <- function(
    eta,
    time_model = NULL,
    control = list()) {
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  if (identical(eta$time_model, "beta")) {
    return(c(log(eta$alpha), log(eta$beta)))
  }
  c(
    stats::qlogis(eta$beta_weight),
    log(eta$alpha),
    log(eta$beta)
  )
}

sunspots_joint_parsimonious_time_unpack_eta <- function(
    par,
    time_model = c("beta", "uniform_beta"),
    control = list()) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  par <- as.numeric(par)
  expected <- sunspots_joint_parsimonious_time_n_parameters(time_model)
  if (length(par) != expected || any(!is.finite(par))) {
    stop(sprintf(
      "The '%s' temporal parameter vector must have %d finite entries.",
      time_model, expected
    ))
  }

  eta <- if (identical(time_model, "beta")) {
    list(
      time_model = time_model,
      alpha = exp(par[[1L]]),
      beta = exp(par[[2L]])
    )
  } else {
    list(
      time_model = time_model,
      beta_weight = stats::plogis(par[[1L]]),
      alpha = exp(par[[2L]]),
      beta = exp(par[[3L]])
    )
  }
  sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
}

sunspots_joint_parsimonious_time_log_density <- function(
    s,
    eta,
    time_model = NULL,
    control = list()) {
  s <- as.numeric(s)
  if (any(!is.finite(s)) || any(s <= 0 | s >= 1)) {
    stop("`s` must contain finite values in (0, 1).")
  }
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  log_beta <- stats::dbeta(s, eta$alpha, eta$beta, log = TRUE)
  if (identical(eta$time_model, "beta")) return(log_beta)

  rotational_logsumexp2(
    rep(log1p(-eta$beta_weight), length(s)),
    log(eta$beta_weight) + log_beta
  )
}

sunspots_joint_parsimonious_time_density <- function(
    s,
    eta,
    time_model = NULL,
    control = list()) {
  exp(sunspots_joint_parsimonious_time_log_density(
    s, eta, time_model = time_model, control = control
  ))
}

sunspots_joint_parsimonious_time_cdf <- function(
    s,
    eta,
    time_model = NULL,
    control = list()) {
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  s <- as.numeric(s)
  out <- numeric(length(s))
  out[s >= 1] <- 1
  active <- is.finite(s) & s > 0 & s < 1
  beta_cdf <- stats::pbeta(s[active], eta$alpha, eta$beta)
  out[active] <- if (identical(eta$time_model, "beta")) {
    beta_cdf
  } else {
    eta$uniform_weight * s[active] + eta$beta_weight * beta_cdf
  }
  out[!is.finite(s) & s > 0] <- 1
  pmin(pmax(out, 0), 1)
}

sunspots_joint_parsimonious_time_moment_match <- function(
    s,
    control = list()) {
  bounds <- sunspots_joint_parsimonious_time_control(control)
  s <- rotational_clamp_unit_interval(as.numeric(s), eps = 1e-10)
  mean_s <- min(max(mean(s), 1e-6), 1 - 1e-6)
  variance_s <- stats::var(s)
  max_variance <- mean_s * (1 - mean_s) * (1 - 1e-8)
  variance_s <- min(max(variance_s, 1e-8), max_variance)

  raw_precision <- mean_s * (1 - mean_s) / variance_s - 1
  precision_lower <- max(
    bounds$shape_lower / mean_s,
    bounds$shape_lower / (1 - mean_s)
  )
  precision_upper <- min(
    bounds$shape_upper / mean_s,
    bounds$shape_upper / (1 - mean_s)
  )
  precision <- min(max(raw_precision, precision_lower), precision_upper)

  list(
    alpha = mean_s * precision,
    beta = (1 - mean_s) * precision
  )
}

sunspots_joint_parsimonious_time_start_etas <- function(
    s,
    time_model = c("beta", "uniform_beta"),
    control = list()) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  base <- sunspots_joint_parsimonious_time_moment_match(s, control)
  bounds <- sunspots_joint_parsimonious_time_control(control)
  mean_s <- base$alpha / (base$alpha + base$beta)
  precision <- base$alpha + base$beta

  scaled_beta <- function(scale) {
    candidate_precision <- precision * scale
    candidate_precision <- min(
      max(
        candidate_precision,
        bounds$shape_lower / mean_s,
        bounds$shape_lower / (1 - mean_s)
      ),
      bounds$shape_upper / mean_s,
      bounds$shape_upper / (1 - mean_s)
    )
    list(
      alpha = mean_s * candidate_precision,
      beta = (1 - mean_s) * candidate_precision
    )
  }

  beta_candidates <- list(
    base,
    scaled_beta(0.5),
    scaled_beta(2),
    list(alpha = 1, beta = 1)
  )
  beta_candidates <- lapply(beta_candidates, function(candidate) {
    sunspots_joint_parsimonious_time_validate_eta(
      c(list(time_model = "beta"), candidate),
      time_model = "beta",
      control = control
    )
  })

  if (identical(time_model, "beta")) return(beta_candidates)

  weights <- c(0.2, 0.5, 0.8)
  out <- list()
  for (candidate in beta_candidates[seq_len(min(3L, length(beta_candidates)))]) {
    for (weight in weights) {
      out[[length(out) + 1L]] <-
        sunspots_joint_parsimonious_time_validate_eta(
          list(
            time_model = "uniform_beta",
            beta_weight = weight,
            alpha = candidate$alpha,
            beta = candidate$beta
          ),
          time_model = "uniform_beta",
          control = control
        )
    }
  }
  out
}

sunspots_joint_parsimonious_time_select_fit <- function(
    fits,
    label = "parsimonious temporal") {
  fits <- Filter(function(fit) {
    !inherits(fit, "try-error") &&
      is.list(fit) &&
      length(fit$value) == 1L &&
      is.finite(fit$value)
  }, fits)
  if (length(fits) == 0L) {
    stop(sprintf("All %s optimizations failed.", label))
  }

  is_converged <- function(fit) {
    convergence <- as.integer(fit$convergence %||% NA_integer_)
    length(convergence) == 1L &&
      !is.na(convergence) &&
      convergence == 0L
  }
  converged <- Filter(is_converged, fits)
  candidate_pool <- if (length(converged) > 0L) converged else fits
  best <- candidate_pool[[
    which.min(vapply(candidate_pool, `[[`, numeric(1L), "value"))
  ]]
  selected_converged <- is_converged(best)
  if (!selected_converged) {
    warning(
      sprintf(
        "The selected %s optimization has convergence != 0; inspect the optimizer result.",
        label
      ),
      call. = FALSE
    )
  }

  list(
    fit = best,
    n_finite_fits = length(fits),
    n_converged_fits = length(converged),
    selected_converged = selected_converged
  )
}

sunspots_joint_parsimonious_time_parameter_boundaries <- function(
    eta,
    time_model = NULL,
    control = list()) {
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  bounds <- sunspots_joint_parsimonious_time_control(control)
  tol <- bounds$boundary_tol

  if (identical(eta$time_model, "beta")) {
    parameter <- c("alpha", "beta")
    value <- c(eta$alpha, eta$beta)
    lower <- rep(bounds$shape_lower, 2L)
    upper <- rep(bounds$shape_upper, 2L)
  } else {
    parameter <- c("beta_weight", "alpha", "beta")
    value <- c(eta$beta_weight, eta$alpha, eta$beta)
    lower <- c(bounds$weight_eps, rep(bounds$shape_lower, 2L))
    upper <- c(1 - bounds$weight_eps, rep(bounds$shape_upper, 2L))
  }

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

sunspots_joint_parsimonious_time_boundary_flags <- function(
    eta,
    time_model = NULL,
    control = list()) {
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  diagnostics <- sunspots_joint_parsimonious_time_parameter_boundaries(
    eta, time_model = eta$time_model, control = control
  )
  bounds <- sunspots_joint_parsimonious_time_control(control)
  identification_distance <- if (identical(eta$time_model, "uniform_beta")) {
    sqrt((eta$alpha - 1)^2 + (eta$beta - 1)^2)
  } else {
    NA_real_
  }

  list(
    weight = if (identical(eta$time_model, "uniform_beta")) {
      diagnostics$near_any_bound[
        diagnostics$parameter == "beta_weight"
      ]
    } else {
      FALSE
    },
    shape_lower = any(
      diagnostics$near_lower_bound[
        diagnostics$parameter %in% c("alpha", "beta")
      ]
    ),
    shape_upper = any(
      diagnostics$near_upper_bound[
        diagnostics$parameter %in% c("alpha", "beta")
      ]
    ),
    identification = identical(eta$time_model, "uniform_beta") &&
      identification_distance <= bounds$uniform_identifiability_tol,
    identification_distance = identification_distance
  )
}

sunspots_joint_parsimonious_time_fast_regular <- function(
    eta,
    time_model = NULL,
    control = list()) {
  flags <- sunspots_joint_parsimonious_time_boundary_flags(
    eta, time_model = time_model, control = control
  )
  !isTRUE(flags$weight) &&
    !isTRUE(flags$shape_lower) &&
    !isTRUE(flags$shape_upper) &&
    !isTRUE(flags$identification)
}

fit_sunspots_joint_parsimonious_time <- function(
    s,
    time_model = c("beta", "uniform_beta"),
    control = list()) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  s <- as.numeric(s)
  if (length(s) < 4L || any(!is.finite(s)) || any(s <= 0 | s >= 1)) {
    stop("`s` must contain at least four finite values in (0, 1).")
  }

  starts <- sunspots_joint_parsimonious_time_start_etas(
    s, time_model = time_model, control = control
  )
  n_starts <- as.integer(
    control$parsimonious_time_n_starts %||% length(starts)
  )
  if (length(n_starts) != 1L ||
      !is.finite(n_starts) ||
      n_starts < 1L) {
    stop("`parsimonious_time_n_starts` must be a positive integer.")
  }
  starts <- starts[seq_len(min(length(starts), n_starts))]

  bounds <- sunspots_joint_parsimonious_time_control(control)
  lower <- if (identical(time_model, "beta")) {
    rep(log(bounds$shape_lower), 2L)
  } else {
    c(
      stats::qlogis(bounds$weight_eps),
      rep(log(bounds$shape_lower), 2L)
    )
  }
  upper <- if (identical(time_model, "beta")) {
    rep(log(bounds$shape_upper), 2L)
  } else {
    c(
      stats::qlogis(1 - bounds$weight_eps),
      rep(log(bounds$shape_upper), 2L)
    )
  }

  objective <- function(par) {
    value <- tryCatch({
      eta <- sunspots_joint_parsimonious_time_unpack_eta(
        par, time_model = time_model, control = control
      )
      -sum(sunspots_joint_parsimonious_time_log_density(
        s, eta, time_model = time_model, control = control
      ))
    }, error = function(error) Inf)
    if (is.finite(value)) value else .Machine$double.xmax / 100
  }

  nelder_mead_control <-
    control$parsimonious_time_nelder_mead_control %||%
    list(maxit = 3000L, reltol = 1e-10)
  lbfgsb_control <-
    control$parsimonious_time_optim_control %||%
    list(maxit = 1500L, factr = 1e7, pgtol = 1e-8)

  fits <- unlist(lapply(starts, function(start) {
    par0 <- sunspots_joint_parsimonious_time_pack_eta(
      start, time_model = time_model, control = control
    )
    par0 <- pmin(pmax(par0, lower), upper)
    exploratory <- try(
      stats::optim(
        par0,
        objective,
        method = "Nelder-Mead",
        control = nelder_mead_control
      ),
      silent = TRUE
    )
    candidates <- list(par0)
    if (!inherits(exploratory, "try-error") &&
        is.finite(exploratory$value)) {
      candidates[[length(candidates) + 1L]] <-
        pmin(pmax(exploratory$par, lower), upper)
    }
    refined <- lapply(candidates, function(candidate) {
      try(
        stats::optim(
          candidate,
          objective,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = lbfgsb_control
        ),
        silent = TRUE
      )
    })

    if (!inherits(exploratory, "try-error") &&
        is.list(exploratory) &&
        length(exploratory$value) == 1L &&
        is.finite(exploratory$value)) {
      c(list(exploratory), refined)
    } else {
      refined
    }
  }), recursive = FALSE)

  selection <- sunspots_joint_parsimonious_time_select_fit(
    fits,
    label = sprintf("'%s' temporal", time_model)
  )
  best <- selection$fit
  eta_hat <- sunspots_joint_parsimonious_time_unpack_eta(
    best$par, time_model = time_model, control = control
  )
  boundary_diagnostics <-
    sunspots_joint_parsimonious_time_parameter_boundaries(
      eta_hat, time_model = time_model, control = control
    )
  boundary_flags <- sunspots_joint_parsimonious_time_boundary_flags(
    eta_hat, time_model = time_model, control = control
  )

  if (!sunspots_joint_parsimonious_time_fast_regular(
    eta_hat, time_model = time_model, control = control
  )) {
    warning(
      sprintf(
        paste(
          "The '%s' temporal MLE is on a numerical boundary or is",
          "not locally identified; the regular fast multiplier correction",
          "must not be used."
        ),
        time_model
      ),
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
    boundary_flags = boundary_flags,
    boundary_diagnostics = boundary_diagnostics,
    fast_regular = sunspots_joint_parsimonious_time_fast_regular(
      eta_hat, time_model = time_model, control = control
    )
  ))
}

sunspots_joint_parsimonious_time_score_matrix <- function(
    s,
    par,
    time_model = c("beta", "uniform_beta"),
    control = list()) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  s <- as.numeric(s)
  if (any(!is.finite(s)) || any(s <= 0 | s >= 1)) {
    stop("`s` must contain finite values in (0, 1).")
  }
  eta <- sunspots_joint_parsimonious_time_unpack_eta(
    par, time_model = time_model, control = control
  )
  score_alpha <- eta$alpha * (
    log(s) -
      digamma(eta$alpha) +
      digamma(eta$alpha + eta$beta)
  )
  score_beta <- eta$beta * (
    log1p(-s) -
      digamma(eta$beta) +
      digamma(eta$alpha + eta$beta)
  )

  if (identical(time_model, "beta")) {
    return(cbind(score_alpha, score_beta))
  }

  log_beta <- stats::dbeta(s, eta$alpha, eta$beta, log = TRUE)
  log_mix <- rotational_logsumexp2(
    rep(log1p(-eta$beta_weight), length(s)),
    log(eta$beta_weight) + log_beta
  )
  responsibility_beta <- exp(
    log(eta$beta_weight) + log_beta - log_mix
  )
  cbind(
    responsibility_beta - eta$beta_weight,
    responsibility_beta * score_alpha,
    responsibility_beta * score_beta
  )
}

sample_sunspots_joint_parsimonious_time <- function(
    n,
    eta,
    time_model = NULL,
    control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a positive integer.")
  }
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  if (identical(eta$time_model, "beta")) {
    out <- stats::rbeta(n, eta$alpha, eta$beta)
  } else {
    beta_component <- stats::runif(n) <= eta$beta_weight
    out <- stats::runif(n)
    out[beta_component] <- stats::rbeta(
      sum(beta_component), eta$alpha, eta$beta
    )
  }
  rotational_clamp_unit_interval(out, eps = .Machine$double.eps)
}

sunspots_joint_parsimonious_time_quadrature <- function(
    eta,
    n_nodes = 64L,
    time_model = NULL,
    control = list()) {
  n_nodes <- as.integer(n_nodes)
  if (length(n_nodes) != 1L ||
      !is.finite(n_nodes) ||
      n_nodes < 2L) {
    stop("`n_nodes` must be an integer >= 2.")
  }
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )

  beta_rule <- beta_mixture2_gauss_jacobi(
    n_nodes,
    alpha = eta$beta - 1,
    beta = eta$alpha - 1
  )
  beta_nodes <- (beta_rule$nodes + 1) / 2
  beta_weights <- beta_rule$weights / beta_rule$total_mass

  if (identical(eta$time_model, "beta")) {
    return(list(
      nodes = beta_nodes,
      weights = beta_weights,
      n_nodes_per_component = n_nodes,
      n_total_nodes = n_nodes,
      mass_error = abs(sum(beta_weights) - 1)
    ))
  }

  uniform_rule <- small_circle_gauss_legendre(n_nodes)
  uniform_nodes <- (uniform_rule$nodes + 1) / 2
  uniform_weights <- uniform_rule$weights / 2

  nodes <- c(uniform_nodes, beta_nodes)
  weights <- c(
    eta$uniform_weight * uniform_weights,
    eta$beta_weight * beta_weights
  )
  order_index <- order(nodes)
  nodes <- nodes[order_index]
  weights <- weights[order_index]

  list(
    nodes = nodes,
    weights = weights,
    n_nodes_per_component = n_nodes,
    n_total_nodes = length(nodes),
    mass_error = abs(sum(weights) - 1)
  )
}

sunspots_joint_parsimonious_time_summary <- function(
    eta,
    time_model = NULL,
    control = list()) {
  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta, time_model = time_model, control = control
  )
  flags <- sunspots_joint_parsimonious_time_boundary_flags(
    eta, time_model = eta$time_model, control = control
  )
  list(
    temporal_model = eta$time_model,
    temporal_n_parameters = eta$n_parameters,
    temporal_uniform_weight = eta$uniform_weight,
    temporal_beta_weight = eta$beta_weight,
    temporal_alpha = eta$alpha,
    temporal_beta = eta$beta,
    temporal_mean = eta$mean,
    temporal_identification_distance =
      flags$identification_distance,
    temporal_boundary_weight = flags$weight,
    temporal_boundary_shape_lower = flags$shape_lower,
    temporal_boundary_shape_upper = flags$shape_upper,
    temporal_identification_failure = flags$identification,
    temporal_fast_regular =
      sunspots_joint_parsimonious_time_fast_regular(
        eta, time_model = eta$time_model, control = control
      )
  )
}

fit_sunspots_cycle23_joint_time_space_parsimonious <- function(
    x,
    s,
    time_model = c("beta", "uniform_beta"),
    hemisphere_regression = "shared",
    control = list(),
    eta_hat = NULL,
    theta_hat = NULL) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(time_model)
  data <- sunspots_joint_validate_data(x, s)

  if (is.null(eta_hat)) {
    eta_hat <- fit_sunspots_joint_parsimonious_time(
      data$s, time_model = time_model, control = control
    )
  } else {
    validated_eta <- sunspots_joint_parsimonious_time_validate_eta(
      eta_hat, time_model = time_model, control = control
    )
    eta_hat <- c(
      validated_eta,
      eta_hat[setdiff(names(eta_hat), names(validated_eta))]
    )
  }
  if (is.null(theta_hat)) {
    theta_hat <- fit_sunspots_time_varying_asymmetric_mixture(
      data$x,
      data$s,
      hemisphere_regression = hemisphere_regression,
      control = control
    )
  }

  temporal_loglik <- sum(
    sunspots_joint_parsimonious_time_log_density(
      data$s, eta_hat, time_model = time_model, control = control
    )
  )
  conditional_loglik <- sunspots_time_varying_loglik(
    data$x, data$s, theta_hat
  )
  n_parameters <- eta_hat$n_parameters + theta_hat$n_parameters
  joint_loglik <- temporal_loglik + conditional_loglik

  list(
    time_model = time_model,
    eta_hat = eta_hat,
    theta_hat = theta_hat,
    temporal_loglik = temporal_loglik,
    conditional_loglik = conditional_loglik,
    loglik = joint_loglik,
    n_parameters = n_parameters,
    aic = 2 * n_parameters - 2 * joint_loglik,
    bic = n_parameters * log(nrow(data$x)) - 2 * joint_loglik
  )
}

sunspots_joint_parsimonious_pack_par <- function(
    fit,
    control = list()) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(
    fit$time_model %||% fit$eta_hat$time_model
  )
  c(
    sunspots_joint_parsimonious_time_pack_eta(
      fit$eta_hat, time_model = time_model, control = control
    ),
    sunspots_time_varying_pack_theta(
      fit$theta_hat,
      nu_eps = as.numeric(control$nu_eps %||% 1e-6),
      c_min = as.numeric(control$c_min %||% 1e-8),
      c_max = as.numeric(control$c_max %||% 1e6),
      hemisphere_regression =
        control$hemisphere_regression %||% "shared"
    )
  )
}

sunspots_joint_parsimonious_state_from_par <- function(
    par,
    control = list()) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(
    control$time_model %||% "beta"
  )
  hemisphere_regression <-
    sunspots_time_varying_normalize_hemisphere_regression(
      control$hemisphere_regression %||% "shared"
    )
  temporal_length <-
    sunspots_joint_parsimonious_time_n_parameters(time_model)
  spatial_length <- if (identical(hemisphere_regression, "shared")) {
    3L
  } else {
    5L
  }
  par <- as.numeric(par)
  if (length(par) != temporal_length + spatial_length) {
    stop("Joint parsimonious parameter vector has an incompatible length.")
  }

  temporal_indices <- seq_len(temporal_length)
  spatial_indices <- temporal_length + seq_len(spatial_length)
  list(
    time_model = time_model,
    eta = sunspots_joint_parsimonious_time_unpack_eta(
      par[temporal_indices],
      time_model = time_model,
      control = control
    ),
    theta = sunspots_time_varying_unpack_par(
      par[spatial_indices],
      nu_eps = as.numeric(control$nu_eps %||% 1e-6),
      c_min = as.numeric(control$c_min %||% 1e-8),
      c_max = as.numeric(control$c_max %||% 1e6),
      hemisphere_regression = hemisphere_regression
    ),
    eta_par = par[temporal_indices],
    theta_par = par[spatial_indices]
  )
}

sunspots_joint_parsimonious_score_matrix <- function(
    data,
    par,
    control = list()) {
  data <- sunspots_joint_validate_data(data$x, data$s)
  state <- sunspots_joint_parsimonious_state_from_par(par, control)
  cbind(
    sunspots_joint_parsimonious_time_score_matrix(
      data$s,
      state$eta_par,
      time_model = state$time_model,
      control = control
    ),
    sunspots_time_gof_score_matrix(
      list(x = data$x, u = data$s),
      state$theta_par,
      control
    )
  )
}

sample_sunspots_joint_time_space_parsimonious <- function(
    n,
    par,
    control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a positive integer.")
  }
  state <- sunspots_joint_parsimonious_state_from_par(par, control)
  s <- sample_sunspots_joint_parsimonious_time(
    n,
    state$eta,
    time_model = state$time_model,
    control = control
  )
  nu <- sunspots_time_varying_nu(s, state$theta)
  north <- stats::runif(n) <= 0.5
  x <- matrix(0, nrow = n, ncol = 3L)
  for (i in seq_len(n)) {
    x[i, ] <- r_sph_small_circle(
      n = 1L,
      mu = if (north[[i]]) c(0, 0, 1) else c(0, 0, -1),
      kappa = state$theta$c,
      nu = if (north[[i]]) nu$north[[i]] else nu$south[[i]],
      check = FALSE
    )
  }
  list(x = x, s = s)
}

plot_sunspots_cycle23_joint_time_space_parsimonious_diagnostics <- function(
    data,
    fit,
    output_dir,
    control = list()) {
  data <- sunspots_joint_validate_data(data$x, data$s)
  time_model <- normalize_sunspots_joint_parsimonious_time_model(
    fit$time_model
  )
  s_grid <- seq(1e-5, 1 - 1e-5, length.out = 501L)

  temporal_path <- file.path(
    output_dir, "cycle23_joint_time_density_pit.png"
  )
  grDevices::png(
    temporal_path, width = 1400, height = 700, res = 140
  )
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(
    data$s,
    breaks = 40,
    freq = FALSE,
    col = "#9ecae1",
    border = "white",
    xlab = "Dequantized first-record day",
    main = sprintf("Temporal density: %s", time_model)
  )
  lines(
    s_grid,
    sunspots_joint_parsimonious_time_density(
      s_grid,
      fit$eta_hat,
      time_model = time_model,
      control = control
    ),
    col = "#8b0000",
    lwd = 2
  )
  time_pit <- sunspots_joint_parsimonious_time_cdf(
    data$s,
    fit$eta_hat,
    time_model = time_model,
    control = control
  )
  qqplot(
    stats::ppoints(length(time_pit)),
    sort(time_pit),
    xlim = c(0, 1),
    ylim = c(0, 1),
    pch = 16,
    cex = 0.45,
    col = "#2b6cb0",
    xlab = "Uniform quantiles",
    ylab = "Temporal PIT quantiles",
    main = "Temporal PIT Q-Q plot"
  )
  abline(0, 1, col = "#8b0000", lwd = 2)
  par(old_par)
  grDevices::dev.off()

  latitude_path <- file.path(
    output_dir, "cycle23_joint_latitude_paths.png"
  )
  latitude <- asin(pmin(pmax(data$x[, 3L], -1), 1)) * 180 / pi
  nu <- sunspots_time_varying_nu(s_grid, fit$theta_hat)
  grDevices::png(
    latitude_path, width = 1400, height = 900, res = 140
  )
  plot(
    data$s,
    latitude,
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor(
      ifelse(latitude >= 0, "#c43c39", "#2b6cb0"),
      alpha.f = 0.22
    ),
    xlab = "Dequantized first-record day",
    ylab = "Latitude (degrees)",
    main = "Cycle 23: fitted north and south small-circle paths"
  )
  lines(
    s_grid,
    asin(nu$north) * 180 / pi,
    col = "#8b0000",
    lwd = 3
  )
  lines(
    s_grid,
    -asin(nu$south) * 180 / pi,
    col = "#003f7f",
    lwd = 3
  )
  grid(col = "#d9d9d9")
  grDevices::dev.off()

  axial_pit <- sunspots_time_varying_conditional_pit(
    data$x[, 3L], data$s, fit$theta_hat
  )
  axial_path <- file.path(
    output_dir, "cycle23_joint_conditional_axial_pit.png"
  )
  grDevices::png(
    axial_path, width = 1400, height = 700, res = 140
  )
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(
    axial_pit,
    breaks = 20,
    freq = FALSE,
    col = "#9ecae1",
    border = "white",
    xlab = "Conditional axial PIT",
    main = "Conditional axial PIT histogram",
    xlim = c(0, 1)
  )
  abline(h = 1, col = "#8b0000", lwd = 2)
  qqplot(
    stats::ppoints(length(axial_pit)),
    sort(axial_pit),
    xlim = c(0, 1),
    ylim = c(0, 1),
    pch = 16,
    cex = 0.45,
    col = "#2b6cb0",
    xlab = "Uniform quantiles",
    ylab = "Conditional PIT quantiles",
    main = "Conditional axial PIT Q-Q plot"
  )
  abline(0, 1, col = "#8b0000", lwd = 2)
  par(old_par)
  grDevices::dev.off()

  list(
    temporal_path = temporal_path,
    latitude_path = latitude_path,
    axial_path = axial_path,
    time_pit = time_pit,
    axial_pit = axial_pit
  )
}
