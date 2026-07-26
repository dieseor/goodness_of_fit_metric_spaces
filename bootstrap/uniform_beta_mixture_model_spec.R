# Uniform-plus-beta-mixture model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_uniform_beta_mixture <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_uniform_beta_mixture <- model_specs_candidates_uniform_beta_mixture[file.exists(model_specs_candidates_uniform_beta_mixture)][1L]
  if (is.na(model_specs_path_uniform_beta_mixture)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_uniform_beta_mixture)
}

prepare_uniform_beta_mixture_fast_multiplier <- function(spec,
                                                         data,
                                                         theta_hat,
                                                         ks_prep = NULL,
                                                         cvm_prep = NULL,
                                                         control = list()) {
  x <- normalize_uniform_beta_mixture_data(data, control)
  theta_hat <- normalize_uniform_beta_mixture_theta(theta_hat, ambient_dim = ncol(x))
  fast_shape_regular_eps <- as.numeric(control$uniform_beta_mixture_fast_shape_regular_eps %||% 0)
  shape_values <- c(theta_hat$alpha, theta_hat$beta)
  if (any(!is.finite(shape_values))) {
    stop("The fitted uniform_beta_mixture shape parameters are not finite.")
  }
  if (min(shape_values) <= 1 + fast_shape_regular_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "uniform_beta_mixture_shape_nonregular"
    ))
  }

  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(
    0,
    0,
    stats::qlogis(theta_hat$weight_uniform),
    log(theta_hat$alpha),
    log(theta_hat$beta)
  )
  weight_eps <- as.numeric(control$uniform_beta_mixture_weight_eps %||% 0.01)
  shape_lower <- as.numeric(control$uniform_beta_mixture_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$uniform_beta_mixture_shape_upper %||% 1e3)
  y_eps <- as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_uniform_beta_mixture_theta(
      list(
        mu = mapped$mu,
        weight_uniform = rotational_bounded_weight(par[[3L]], weight_eps = weight_eps),
        alpha = rotational_positive_parameter(par[[4L]], lower = shape_lower, upper = shape_upper),
        beta = rotational_positive_parameter(par[[5L]], lower = shape_lower, upper = shape_upper)
      ),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_uniform_beta_mixture_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    y <- rotational_clamp_unit_interval((z + 1) / 2, eps = y_eps)
    log_b <- stats::dbeta(y, theta_state$alpha, theta_state$beta, log = TRUE)
    log_m <- rotational_logsumexp2(
      log(theta_state$weight_uniform),
      log1p(-theta_state$weight_uniform) + log_b
    )
    r_beta <- exp(log1p(-theta_state$weight_uniform) + log_b - log_m)
    r_uniform <- 1 - r_beta

    g_beta <- (theta_state$alpha - 1) / y - (theta_state$beta - 1) / (1 - y)
    coeff_mu <- 0.5 * r_beta * g_beta
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))

    cbind(
      score_mu,
      r_uniform - theta_state$weight_uniform,
      theta_state$alpha * r_beta * (
        log(y) - digamma(theta_state$alpha) + digamma(theta_state$alpha + theta_state$beta)
      ),
      theta_state$beta * r_beta * (
        log1p(-y) - digamma(theta_state$beta) + digamma(theta_state$alpha + theta_state$beta)
      )
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_uniform_beta_mixture(
      n = n_aux,
      mu = theta_state$mu,
      weight_uniform = theta_state$weight_uniform,
      alpha = theta_state$alpha,
      beta = theta_state$beta
    )
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}
normalize_uniform_beta_mixture_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Rotational uniform-beta-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_uniform_beta_mixture_theta <- function(theta, ambient_dim = 3L) {
  uniform_beta_mixture_normalize_theta(theta, ambient_dim = ambient_dim)
}

fit_uniform_beta_mixture_theta <- function(data,
                                           weights = NULL,
                                           null,
                                           control = list()) {
  x <- normalize_uniform_beta_mixture_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_uniform_beta_mixture_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- uniform_beta_mixture_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_uniform_beta_mixture_theta(fit, ambient_dim = ncol(x))
}

make_uniform_beta_mixture_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("uniform_beta_mixture_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_uniform_beta_mixture_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_uniform_beta_mixture_data(data, control)
      omega_matrix <- normalize_uniform_beta_mixture_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- pmin(pmax(x %*% t(omega_matrix), -1), 1)
      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_uniform_beta_mixture_theta(theta)
      distance_profile_uniform_beta_mixture(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta,
        distance_type = distance_type,
        method = control$uniform_beta_mixture_profile_method %||% "legendre",
        l_max = as.integer(control$uniform_beta_mixture_L_max %||% 150L),
        quad_n = as.integer(control$uniform_beta_mixture_quad_n %||% 1000L),
        tol = as.numeric(control$uniform_beta_mixture_tol %||% 1e-6),
        eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12),
        validate_against_integral = isTRUE(control$uniform_beta_mixture_validate_against_integral),
        validation_tol = as.numeric(control$uniform_beta_mixture_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_uniform_beta_mixture_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_uniform_beta_mixture_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_uniform_beta_mixture_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_uniform_beta_mixture_theta(theta)
        distance_profile_uniform_beta_mixture_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$uniform_beta_mixture_profile_method %||% "legendre",
          l_max = as.integer(control$uniform_beta_mixture_L_max %||% 150L),
          quad_n = as.integer(control$uniform_beta_mixture_quad_n %||% 1000L),
          tol = as.numeric(control$uniform_beta_mixture_tol %||% 1e-6),
          eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_uniform_beta_mixture_theta(theta)
        distance_profile_uniform_beta_mixture_cvm_grid(
          X = data,
          mu = theta$mu,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          distance_matrix = distance_matrix,
          method = control$uniform_beta_mixture_profile_method %||% "legendre",
          l_max = as.integer(control$uniform_beta_mixture_L_max %||% 150L),
          quad_n = as.integer(control$uniform_beta_mixture_quad_n %||% 1000L),
          tol = as.numeric(control$uniform_beta_mixture_tol %||% 1e-6),
          eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_uniform_beta_mixture_fast_multiplier(
          spec = make_uniform_beta_mixture_spec(distance_type = distance_type),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control
        )
      }
    )
  )
}
