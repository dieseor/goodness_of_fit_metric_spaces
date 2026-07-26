# Hyperbolic von Mises--Fisher model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_hvmf <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_hvmf <- model_specs_candidates_hvmf[file.exists(model_specs_candidates_hvmf)][1L]
  if (is.na(model_specs_path_hvmf)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_hvmf)
}

normalize_hvmf_data <- function(data, control = list()) {
  tol <- as.numeric(control$hvmf_tol %||% 1e-10)
  normalize_hvmf_h2_data(data, tol = tol)
}

normalize_hvmf_theta <- function(theta, control = list()) {
  tol <- as.numeric(control$hvmf_tol %||% 1e-10)

  if (!is.list(theta)) {
    stop("HvMF theta must be a list containing `mu` and `kappa`.")
  }

  mu <- theta$mu
  kappa <- as.numeric(theta$kappa)

  mu_matrix <- normalize_hvmf_h2_data(mu, tol = tol)
  mu <- as.numeric(mu_matrix[1L, , drop = TRUE])

  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("HvMF theta requires a strictly positive finite scalar `kappa`.")
  }

  sinh_chi <- sqrt(sum(mu[-1L]^2))
  chi <- asinh(sinh_chi)
  theta_angle <- atan2(mu[[3L]], mu[[2L]])
  theta_deg <- (theta_angle * 180 / pi) %% 360

  list(
    mu = mu,
    kappa = kappa,
    chi = chi,
    sinh_chi = sinh_chi,
    theta = theta_angle,
    theta_deg = theta_deg,
    q = 2L
  )
}

fit_hvmf_theta <- function(data,
                           weights = NULL,
                           null,
                           unknown_param = "both",
                           control = list()) {
  x <- normalize_hvmf_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_hvmf_theta(null$theta, control = control))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (!is.null(unknown_param) && !identical(unknown_param, "both")) {
    stop("The HvMF adapter currently supports only `unknown_param = 'both'`.")
  }

  tol <- as.numeric(control$hvmf_tol %||% 1e-10)
  fit <- hvmf_mle_h2(x, weights = weights, tol = tol)

  normalize_hvmf_theta(
    list(mu = fit$mu, kappa = fit$kappa),
    control = control
  )
}

prepare_hvmf_fast_multiplier <- function(spec,
                                         data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list(),
                                         unknown_param = "both") {
  if (!identical(unknown_param, "both")) {
    stop("The fast HvMF multiplier bootstrap currently supports only `unknown_param = 'both'`.")
  }

  x <- normalize_hvmf_data(data, control)
  theta_hat <- normalize_hvmf_theta(theta_hat, control = control)
  par0 <- c(theta_hat$chi, theta_hat$theta, theta_hat$kappa)

  state_from_par <- function(par) {
    chi <- as.numeric(par[[1L]])
    theta_angle <- as.numeric(par[[2L]])
    kappa <- as.numeric(par[[3L]])
    mu <- c(
      cosh(chi),
      sinh(chi) * cos(theta_angle),
      sinh(chi) * sin(theta_angle)
    )
    normalize_hvmf_theta(list(mu = mu, kappa = kappa), control = control)
  }

  dmu_dchi <- function(chi, theta_angle) {
    c(
      sinh(chi),
      cosh(chi) * cos(theta_angle),
      cosh(chi) * sin(theta_angle)
    )
  }

  dmu_dtheta <- function(chi, theta_angle) {
    c(
      0,
      -sinh(chi) * sin(theta_angle),
      sinh(chi) * cos(theta_angle)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_hvmf_data(sample, control)
    chi <- theta_state$chi
    theta_angle <- theta_state$theta
    mu <- theta_state$mu
    kappa <- theta_state$kappa
    dchi_vec <- dmu_dchi(chi, theta_angle)
    dtheta_vec <- dmu_dtheta(chi, theta_angle)
    mink_mu <- apply(sample, 1L, function(row) minkowski_inner_product(row, mu))
    score_loc_chi <- kappa * apply(sample, 1L, function(row) minkowski_inner_product(row, dchi_vec))
    score_loc_theta <- kappa * apply(sample, 1L, function(row) minkowski_inner_product(row, dtheta_vec))
    cbind(
      score_loc_chi,
      score_loc_theta,
      1 / kappa + 1 + mink_mu
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    rhvmf_h2_polar(n_aux, mu = theta_state$mu, kappa = theta_state$kappa)
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

make_hvmf_spec <- function(unknown_param = "both") {
  new_model_spec(
    name = "hvmf_geodesic",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_hvmf_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_hvmf_data(data, control)
      omega_matrix <- normalize_hvmf_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      t(hvmf_distance_matrix(omega_matrix, x))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_hvmf_theta(theta, control = control)
      theoretical_distance_profile_hvmf(
        omega = omega,
        mu = theta$mu,
        kappa = theta$kappa,
        t_values = as.numeric(t)
      )
    },
    normalize_data = normalize_hvmf_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_hvmf_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_hvmf_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_hvmf_theta(theta, control = control)
        profile_method <- tolower(as.character(control$hvmf_profile_method %||% "tabulated"))
        n_y <- as.integer(control$hvmf_profile_n_y %||% control$hvmf_profile_grid_size %||% 4097L)

        if (!profile_method %in% c("exact", "tabulated")) {
          stop("`control$hvmf_profile_method` must be either 'exact' or 'tabulated'.")
        }

        if (identical(profile_method, "tabulated")) {
          return(hvmf_cvm_profile_matrix_tabulated(
            data = data,
            theta = theta,
            grid_size = n_y,
            distance_matrix = distance_matrix,
            tol = as.numeric(control$hvmf_tol %||% 1e-10)
          ))
        }

        NULL
      },
      distance_type = "geodesic",
      unknown_param = unknown_param,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_hvmf_fast_multiplier(
          spec = make_hvmf_spec(unknown_param = unknown_param),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          unknown_param = unknown_param
        )
      }
    )
  )
}
