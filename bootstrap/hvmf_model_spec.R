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
  normalize_hvmf_hq_data(data, tol = tol)
}

normalize_hvmf_theta <- function(theta, control = list()) {
  tol <- as.numeric(control$hvmf_tol %||% 1e-10)

  if (!is.list(theta)) {
    stop("HvMF theta must be a list containing `mu` and `kappa`.")
  }

  mu <- theta$mu
  kappa <- as.numeric(theta$kappa)

  mu_matrix <- normalize_hvmf_hq_data(mu, tol = tol)
  mu <- as.numeric(mu_matrix[1L, , drop = TRUE])
  q <- length(mu) - 1L

  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("HvMF theta requires a strictly positive finite scalar `kappa`.")
  }

  sinh_chi <- sqrt(sum(mu[-1L]^2))
  chi <- asinh(sinh_chi)
  output <- list(
    xi = kappa * mu,
    mu = mu,
    kappa = kappa,
    chi = chi,
    sinh_chi = sinh_chi,
    spatial_direction = if (sinh_chi > sqrt(.Machine$double.eps)) {
      mu[-1L] / sinh_chi
    } else {
      c(1, rep.int(0, q - 1L))
    },
    q = q
  )
  if (q == 2L) {
    theta_angle <- atan2(mu[[3L]], mu[[2L]])
    output$theta <- theta_angle
    output$theta_deg <- (theta_angle * 180 / pi) %% 360
  }
  output
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
  fit <- hvmf_mle_hq(x, weights = weights, tol = tol)

  normalize_hvmf_theta(
    list(mu = fit$mu, kappa = fit$kappa),
    control = control
  )
}

hvmf_hq_theta_from_coordinates <- function(par, q, control = list()) {
  par <- as.numeric(par)
  if (length(par) != q + 1L || any(!is.finite(par))) {
    stop("General HvMF coordinates must have length q + 1.")
  }
  eta <- par[seq_len(q)]
  kappa <- exp(par[[q + 1L]])
  mu <- c(sqrt(1 + sum(eta^2)), eta)
  normalize_hvmf_theta(list(mu = mu, kappa = kappa), control = control)
}

prepare_hvmf_hq_fast_multiplier <- function(spec,
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
  q <- theta_hat$q
  eta_hat <- theta_hat$mu[-1L]
  par0 <- c(eta_hat, log(theta_hat$kappa))

  state_from_par <- function(par) {
    hvmf_hq_theta_from_coordinates(par, q = q, control = control)
  }
  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_hvmf_data(sample, control)
    eta <- theta_state$mu[-1L]
    mu0 <- theta_state$mu[[1L]]
    kappa <- theta_state$kappa
    minkowski_mu <- -sample[, 1L] * mu0 +
      rowSums(sample[, -1L, drop = FALSE] * rep(eta, each = nrow(sample)))
    location_score <- kappa * (
      sample[, -1L, drop = FALSE] - tcrossprod(sample[, 1L], eta / mu0)
    )
    kappa_score <- kappa * (
      hvmf_mean_resultant_ratio(q, kappa) + minkowski_mu
    )
    cbind(location_score, kappa_score)
  }
  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    rhvmf_polar(n_aux, mu = theta_state$mu, kappa = theta_state$kappa)
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

prepare_hvmf_quadrature_fast_multiplier <- function(spec,
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
  q <- theta_hat$q
  p <- q + 1L
  xi_hat <- theta_hat$xi
  S_obs <- hvmf_canonical_score_matrix(x, xi = xi_hat)

  evaluate_derivative <- function(omega, thresholds) {
    hvmf_profile_and_derivative_xi(
      omega = omega,
      xi = xi_hat,
      t_values = thresholds,
      grid_size = as.integer(
        control$hvmf_derivative_n_y %||%
          control$hvmf_profile_n_y %||%
          control$hvmf_profile_grid_size %||%
          4097L
      )
    )$derivative
  }

  grid_size <- as.integer(
    control$hvmf_derivative_n_y %||%
      control$hvmf_profile_n_y %||%
      control$hvmf_profile_grid_size %||%
      4097L
  )
  D_ks <- if (is.null(ks_prep)) {
    NULL
  } else if (identical(
    ks_prep$ks_grid_mode %||% "fixed",
    "sample_points_unique_distances"
  )) {
    centers <- normalize_hvmf_data(ks_prep$omega_grid, control)
    list(
      mode = "sample_points_unique_distances",
      derivative_sorted = profile_derivative_stack_centers(
        centers = centers,
        thresholds = ks_prep$sorted_distance_matrix,
        evaluator = evaluate_derivative
      )
    )
  } else {
    centers <- normalize_hvmf_data(ks_prep$omega_grid, control)
    do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
      evaluate_derivative(centers[i, ], ks_prep$t_grid)
    }))
  }

  D_cvm <- if (is.null(cvm_prep)) {
    NULL
  } else if (isTRUE(cvm_prep$shared_with_ks) &&
             is.list(D_ks) &&
             !is.null(D_ks$derivative_sorted)) {
    list(
      mode = "sample_points_unique_distances_sorted_rows",
      derivative_sorted = D_ks$derivative_sorted,
      shared_with_ks = TRUE
    )
  } else if (isTRUE(cvm_prep$light) &&
             !is.null(cvm_prep$sorted_distance_matrix)) {
    list(
      mode = "sample_points_unique_distances_sorted_rows",
      derivative_sorted = profile_derivative_stack_centers(
        centers = x,
        thresholds = cvm_prep$sorted_distance_matrix,
        evaluator = evaluate_derivative
      )
    )
  } else {
    observed_distances <- cvm_prep$distance_matrix
    if (is.null(observed_distances)) {
      observed_distances <- spec$distance_matrix(x, x, control)
    }
    profile_derivative_stack_centers(
      centers = x,
      thresholds = observed_distances,
      evaluator = evaluate_derivative
    )
  }
  Vhat <- hvmf_canonical_information(xi_hat)
  vhat_diagnostics <- fast_multiplier_deterministic_vhat_diagnostics(
    S_obs = S_obs,
    Vhat = Vhat,
    par0 = xi_hat
  )
  fast_multiplier_validate_vhat(
    Vhat = Vhat,
    diagnostics = vhat_diagnostics,
    label = "Monte Carlo HvMF fast multiplier preparation",
    rcond_tol = as.numeric(control$fast_multiplier_vhat_rcond_tol %||% 1e-12)
  )
  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = matrix(numeric(0), nrow = 0L, ncol = p),
    vhat_method = "analytic_fisher_information",
    vhat_diagnostics = vhat_diagnostics,
    observed_cvm_distance_matrix = if (!is.null(cvm_prep) && !isTRUE(cvm_prep$light)) {
      cvm_prep$distance_matrix %||% spec$distance_matrix(x, x, control)
    } else {
      NULL
    },
    derivative_method = "quadrature",
    derivative_mc_size = NA_integer_,
    derivative_mc_seed = NA_integer_,
    derivative_grid_size = grid_size,
    canonical_parameter = "xi",
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

prepare_hvmf_fast_multiplier <- function(spec,
                                         data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list(),
                                         unknown_param = "both") {
  legacy_mc_control <- !is.null(control$derivative_mc_size) ||
    !is.null(control$derivative_mc_seed)
  if (is.null(control$derivative_method) && legacy_mc_control) {
    warning(
      paste(
        "HvMF fast multiplier: `derivative_method` was not supplied, but",
        "legacy `derivative_mc_size`/`derivative_mc_seed` controls were found.",
        "They no longer change the default: selecting `quadrature`.",
        "Set `derivative_method = 'score_mc'` explicitly to use Monte Carlo."
      ),
      call. = FALSE
    )
  }
  derivative_control <- fast_multiplier_parse_derivative_control(
    control,
    default_method = "quadrature"
  )
  if (identical(derivative_control$derivative_method, "quadrature")) {
    return(prepare_hvmf_quadrature_fast_multiplier(
      spec = spec,
      data = data,
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control,
      unknown_param = unknown_param
    ))
  }

  if (!identical(unknown_param, "both")) {
    stop("The fast HvMF multiplier bootstrap currently supports only `unknown_param = 'both'`.")
  }
  x <- normalize_hvmf_data(data, control)
  theta_hat <- normalize_hvmf_theta(theta_hat, control = control)
  if (theta_hat$q != 2L) {
    return(prepare_hvmf_hq_fast_multiplier(
      spec = spec,
      data = x,
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control,
      unknown_param = unknown_param
    ))
  }
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
    mink_mu <- apply(sample, 1L, function(row) {
      minkowski_inner_product(row, mu)
    })
    score_loc_chi <- kappa * apply(sample, 1L, function(row) {
      minkowski_inner_product(row, dchi_vec)
    })
    score_loc_theta <- kappa * apply(sample, 1L, function(row) {
      minkowski_inner_product(row, dtheta_vec)
    })
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

      t(hvmf_distance_matrix_hq(omega_matrix, x))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_hvmf_theta(theta, control = control)
      grid_size <- as.integer(
        control$hvmf_profile_n_y %||% control$hvmf_profile_grid_size %||% 4097L
      )
      if (!is.finite(grid_size) || grid_size < 3L) {
        stop("`control$hvmf_profile_n_y` must be an integer of at least three.")
      }
      if (theta$q == 5L) {
        return(hvmf_distance_profile_hq_integral(
          omega = omega, mu = theta$mu, kappa = theta$kappa,
          t_values = as.numeric(t)
        ))
      }
      if (theta$q != 2L) {
        return(hvmf_distance_profile_hq(
          omega = omega, mu = theta$mu, kappa = theta$kappa,
          t_values = as.numeric(t), grid_size = grid_size
        ))
      }
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
        if (theta$q == 5L) {
          # Returning NULL delegates to `profile_eval()` above, which uses
          # adaptive integration for the Section 6 H^5 experiment.  In
          # particular, its observed KS/CvM statistics do not use the
          # tabulated radial profile.
          return(NULL)
        }
        if (theta$q != 2L) {
          x <- normalize_hvmf_data(data, control)
          output <- matrix(0, nrow = nrow(x), ncol = ncol(distance_matrix))
          for (i in seq_len(nrow(x))) {
            output[i, ] <- hvmf_distance_profile_hq(
              omega = x[i, ], mu = theta$mu, kappa = theta$kappa,
              t_values = distance_matrix[i, ],
              grid_size = as.integer(
                control$hvmf_profile_n_y %||% control$hvmf_profile_grid_size %||% 4097L
              )
            )
          }
          return(output)
        }
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
