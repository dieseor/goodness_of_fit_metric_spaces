# Dimroth--Watson model adapter for multiplier bootstrap GOF tests.
# This is the exact nu = 0 submodel of Bingham--Mardia Small Circle.

resolve_watson_model_spec_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

model_specs_path_watson <- resolve_watson_model_spec_path("bootstrap", "model_specs.R")
if (!exists("new_model_spec", mode = "function")) source(model_specs_path_watson)

normalize_watson_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) stop("Watson data must be a non-empty S^2 sample with three columns.")
  x
}

normalize_watson_theta <- function(theta, ambient_dim = 3L) {
  if (!is.list(theta)) stop("Watson theta must be a list containing `mu` and `kappa`.")
  params <- watson_validate_parameters(theta$mu, theta$kappa)
  if (length(params$mu) != ambient_dim) stop("Watson theta has incompatible ambient dimension.")
  list(mu = watson_canonicalize_axis(params$mu), kappa = params$kappa, ambient_dim = length(params$mu))
}

fit_watson_theta <- function(data, weights = NULL, null, control = list()) {
  x <- normalize_watson_data(data, control)
  if (!is.list(null) || is.null(null$type)) stop("`null` must be a list containing at least the field `type`.")
  if (identical(null$type, "simple")) return(normalize_watson_theta(null$theta, ncol(x)))
  if (!identical(null$type, "composite")) stop("`null$type` must be either `simple` or `composite`.")
  normalize_watson_theta(watson_mle_s2_weighted(x, weights = weights, control = control), ncol(x))
}

prepare_watson_fast_multiplier <- function(spec,
                                            data,
                                            theta_hat,
                                            ks_prep = NULL,
                                            cvm_prep = NULL,
                                            control = list()) {
  x <- normalize_watson_data(data, control)
  theta_hat <- normalize_watson_theta(theta_hat, ncol(x))
  regular_tol <- as.numeric(control$watson_fast_regular_kappa_tol %||% 1e-6)
  if (!is.finite(regular_tol) || regular_tol <= 0) stop("`watson_fast_regular_kappa_tol` must be positive and finite.")
  if (theta_hat$kappa <= regular_tol) {
    stop("The fast Watson multiplier bootstrap is not regular at kappa = 0; use a simple uniform null or a regular Watson fit.")
  }

  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(0, 0, theta_hat$kappa)
  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    # Do not canonicalize inside this local chart: mu and -mu are the same
    # Watson axis, but a sign convention would introduce an artificial chart
    # discontinuity when finite-difference V-hat is requested.
    watson_validate_parameters(mapped$mu, par[[3L]])
  }
  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_watson_data(sample, control)
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    score_mu <- sweep(sample, 1L, -2 * theta_state$kappa * z, "*") %*% jac_mu
    cbind(score_mu, watson_axis_second_moment(theta_state$kappa) - z^2)
  }
  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_watson(n_aux, theta_state$mu, theta_state$kappa)
  }
  prepare_fast_multiplier_score_model(
    spec = spec, data = x, theta_hat = theta_hat, ks_prep = ks_prep,
    cvm_prep = cvm_prep, control = control, par0 = par0,
    score_matrix_fn = score_matrix_fn, sample_fn = sample_fn
  )
}

make_watson_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)
  new_model_spec(
    name = sprintf("watson_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_watson_theta(data, weights, null, control)
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_watson_data(data, control)
      omega <- normalize_watson_data(omega, control)
      dots <- pmin(pmax(x %*% t(omega), -1), 1)
      if (identical(distance_type, "geodesic")) acos(dots) else sqrt(pmax(0, 2 * (1 - dots)))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_watson_theta(theta)
      distance_profile_watson(
        omega, as.numeric(t), theta$mu, theta$kappa, distance_type,
        method = control$watson_profile_method %||% "legendre",
        l_max = as.integer(control$watson_L_max %||% 200L),
        quad_n = as.integer(control$watson_quad_n %||% 400L),
        tol = as.numeric(control$watson_tol %||% 1e-10)
      )
    },
    normalize_data = normalize_watson_data,
    n_obs = function(data, control = list()) nrow(normalize_watson_data(data, control)),
    observation_at = function(data, idx, control = list()) normalize_watson_data(data, control)[idx, , drop = TRUE],
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_watson_theta(theta)
        distance_profile_watson_grid(
          omega_grid, theta$mu, theta$kappa, t_grid, distance_type,
          l_max = as.integer(control$watson_L_max %||% 200L),
          quad_n = as.integer(control$watson_quad_n %||% 400L),
          tol = as.numeric(control$watson_tol %||% 1e-10)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) return(NULL)
        theta <- normalize_watson_theta(theta)
        distance_profile_watson_cvm_grid(
          data, theta$mu, theta$kappa,
          l_max = as.integer(control$watson_L_max %||% 200L),
          quad_n = as.integer(control$watson_quad_n %||% 400L),
          tol = as.numeric(control$watson_tol %||% 1e-10)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data, theta_hat, ks_prep = NULL, cvm_prep = NULL, control = list()) {
        prepare_watson_fast_multiplier(make_watson_spec(distance_type), data, theta_hat, ks_prep, cvm_prep, control)
      }
    )
  )
}
