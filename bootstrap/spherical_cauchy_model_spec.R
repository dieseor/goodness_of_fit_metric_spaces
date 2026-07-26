# Spherical-Cauchy model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_spherical_cauchy <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_spherical_cauchy <- model_specs_candidates_spherical_cauchy[file.exists(model_specs_candidates_spherical_cauchy)][1L]
  if (is.na(model_specs_path_spherical_cauchy)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_spherical_cauchy)
}

normalize_spherical_cauchy_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Spherical Cauchy data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_spherical_cauchy_theta <- function(theta, ambient_dim = 3L) {
  theta <- spherical_cauchy_normalize_theta(theta, ambient_dim = ambient_dim)
  if (theta$ambient_dim != 3L) {
    stop("Spherical Cauchy theta currently supports only S^2.")
  }
  theta
}

fit_spherical_cauchy_theta <- function(data,
                                       weights = NULL,
                                       null,
                                       control = list()) {
  x <- normalize_spherical_cauchy_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_spherical_cauchy_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- spherical_cauchy_mle_s2_weighted(
    data = x,
    weights = weights,
    control = control
  )

  normalize_spherical_cauchy_theta(fit, ambient_dim = ncol(x))
}

prepare_spherical_cauchy_fast_multiplier <- function(spec,
                                                     data,
                                                     theta_hat,
                                                     ks_prep = NULL,
                                                     cvm_prep = NULL,
                                                     control = list()) {
  x <- normalize_spherical_cauchy_data(data, control)
  theta_hat <- normalize_spherical_cauchy_theta(theta_hat, ambient_dim = ncol(x))
  rho_zero_eps <- as.numeric(control$spherical_cauchy_fast_boundary_eps %||%
    control$spherical_cauchy_fast_zero_eps %||% 1e-8)
  if (theta_hat$rho <= rho_zero_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "spherical_cauchy_rho_zero_nonidentification",
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_
    ))
  }
  par0 <- as.numeric(theta_hat$phi)

  state_from_par <- function(par) {
    normalize_spherical_cauchy_theta(list(phi = par), ambient_dim = ncol(x))
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_spherical_cauchy_data(sample, control)
    phi <- as.numeric(theta_state$phi)
    phi_norm_sq <- sum(phi^2)
    denom <- as.numeric(1 - 2 * (sample %*% phi) + phi_norm_sq)
    const_term <- matrix(
      rep(-2 * phi / (1 - phi_norm_sq), each = nrow(sample)),
      nrow = nrow(sample),
      ncol = length(phi)
    )
    phi_matrix <- matrix(
      rep(phi, each = nrow(sample)),
      nrow = nrow(sample),
      ncol = length(phi)
    )
    const_term + 3 * sample / denom - sweep(phi_matrix, 1L, 3 / denom, FUN = "*")
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_spherical_cauchy(
      n = n_aux,
      mu = theta_state$mu,
      rho = theta_state$rho
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

make_spherical_cauchy_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("spherical_cauchy_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_spherical_cauchy_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_spherical_cauchy_data(data, control)
      omega_matrix <- normalize_spherical_cauchy_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- x %*% t(omega_matrix)
      dot_products <- pmin(pmax(dot_products, -1), 1)

      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_spherical_cauchy_theta(theta)
      l_max <- control$spherical_cauchy_profile_l_max %||% NULL
      max_l_max <- control$spherical_cauchy_profile_max_l_max %||% NULL
      if (!is.null(l_max)) {
        l_max <- as.integer(l_max)
      }
      if (!is.null(max_l_max)) {
        max_l_max <- as.integer(max_l_max)
      }
      distance_profile_spherical_cauchy(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        rho = theta$rho,
        distance_type = distance_type,
        tail_tol = as.numeric(control$spherical_cauchy_profile_tol %||% 1e-10),
        l_max = l_max,
        max_l_max = max_l_max,
        warn = isTRUE(control$spherical_cauchy_profile_warn %||% FALSE)
      )
    },
    normalize_data = normalize_spherical_cauchy_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_spherical_cauchy_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_spherical_cauchy_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_spherical_cauchy_theta(theta)

        l_max <- control$spherical_cauchy_profile_l_max %||% NULL
        max_l_max <- control$spherical_cauchy_profile_max_l_max %||% NULL
        if (!is.null(l_max)) {
          l_max <- as.integer(l_max)
        }
        if (!is.null(max_l_max)) {
          max_l_max <- as.integer(max_l_max)
        }

        distance_profile_spherical_cauchy_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          rho = theta$rho,
          t_grid = t_grid,
          distance_type = distance_type,
          tail_tol = as.numeric(control$spherical_cauchy_profile_tol %||% 1e-10),
          l_max = l_max,
          max_l_max = max_l_max,
          warn = isTRUE(control$spherical_cauchy_profile_warn %||% FALSE)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_spherical_cauchy_theta(theta)

        l_max <- control$spherical_cauchy_profile_l_max %||% NULL
        max_l_max <- control$spherical_cauchy_profile_max_l_max %||% NULL
        if (!is.null(l_max)) {
          l_max <- as.integer(l_max)
        }
        if (!is.null(max_l_max)) {
          max_l_max <- as.integer(max_l_max)
        }

        distance_profile_spherical_cauchy_cvm_grid(
          data = data,
          mu = theta$mu,
          rho = theta$rho,
          distance_matrix = distance_matrix,
          distance_type = distance_type,
          tail_tol = as.numeric(control$spherical_cauchy_profile_tol %||% 1e-10),
          l_max = l_max,
          max_l_max = max_l_max,
          warn = isTRUE(control$spherical_cauchy_profile_warn %||% FALSE)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_spherical_cauchy_fast_multiplier(
          spec = make_spherical_cauchy_spec(distance_type = distance_type),
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
