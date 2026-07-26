# von Mises--Fisher model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_vmf <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_vmf <- model_specs_candidates_vmf[file.exists(model_specs_candidates_vmf)][1L]
  if (is.na(model_specs_path_vmf)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_vmf)
}

normalize_vmf_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0 || ncol(data) < 2L) {
    stop("vMF data must be a non-empty matrix with at least two columns.")
  }
  if (any(!is.finite(data))) {
    stop("vMF data must be finite.")
  }

  norms <- sqrt(rowSums(data^2))
  if (any(norms <= 0)) {
    stop("vMF data rows must have strictly positive norm.")
  }

  data / norms
}

default_unit_vector <- function(dim) {
  output <- rep(0, dim)
  output[[1]] <- 1
  output
}

normalize_vmf_theta <- function(theta, ambient_dim = NULL) {
  if (is.numeric(theta) && !is.list(theta)) {
    xi <- as.numeric(theta)
    if (length(xi) < 2L) {
      stop("vMF xi must have length at least 2.")
    }
    if (any(!is.finite(xi))) {
      stop("vMF xi must be finite.")
    }

    kappa <- sqrt(sum(xi^2))
    mu <- if (kappa > 0) {
      xi / kappa
    } else {
      default_unit_vector(length(xi))
    }

    return(list(
      xi = xi,
      mu = mu,
      kappa = kappa,
      q = length(xi) - 1L
    ))
  }

  if (!is.list(theta)) {
    stop("vMF theta must be either a xi vector or a list containing `mu` and `kappa`.")
  }

  if (!is.null(theta$xi)) {
    return(normalize_vmf_theta(theta$xi, ambient_dim = ambient_dim))
  }

  mu <- as.numeric(theta$mu)
  kappa <- as.numeric(theta$kappa)

  if (length(mu) < 2L || any(!is.finite(mu))) {
    stop("vMF theta requires a finite vector `mu` of length at least 2.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("vMF theta requires a nonnegative finite scalar `kappa`.")
  }

  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 0) {
    stop("vMF theta requires `mu` with strictly positive norm.")
  }

  mu <- mu / mu_norm
  xi <- if (kappa > 0) kappa * mu else rep(0, length(mu))
  if (!is.null(ambient_dim) && length(mu) != ambient_dim) {
    stop("vMF theta has incompatible dimension.")
  }

  list(
    xi = xi,
    mu = mu,
    kappa = kappa,
    q = length(mu) - 1L
  )
}

solve_vmf_kappa_from_rbar <- function(r_bar,
                                      q,
                                      tol = 1e-10,
                                      max_kappa = 1e6) {
  r_bar <- as.numeric(r_bar)
  if (length(r_bar) != 1L || !is.finite(r_bar) || r_bar < 0 || r_bar > 1) {
    stop("`r_bar` must be a scalar in [0, 1].")
  }
  if (length(q) != 1L || !is.finite(q) || q < 1) {
    stop("`q` must be a positive scalar.")
  }

  if (r_bar <= tol) {
    return(0)
  }
  if (r_bar >= 1 - tol) {
    return(max_kappa)
  }

  objective <- function(kappa) A_q(kappa, q) - r_bar
  upper <- 1

  while (objective(upper) < 0 && upper < max_kappa) {
    upper <- upper * 2
  }

  if (upper >= max_kappa && objective(upper) < 0) {
    return(max_kappa)
  }

  stats::uniroot(objective, interval = c(0, upper), tol = tol)$root
}

fit_vmf_theta <- function(data,
                          weights = NULL,
                          null,
                          unknown_param = "xi",
                          control = list()) {
  x <- normalize_vmf_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_vmf_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (!identical(unknown_param, "xi")) {
    stop("The vMF adapter currently supports only `unknown_param = 'xi'`.")
  }

  if (is.null(weights)) {
    xi_hat <- compute_mle_xi(x)
    return(normalize_vmf_theta(xi_hat, ambient_dim = ncol(x)))
  }

  prob_weights <- normalize_probability_weights(weights, nrow(x))
  resultant <- colSums(x * prob_weights)
  r_bar <- sqrt(sum(resultant^2))

  if (r_bar <= 1e-12) {
    mu_hat <- default_unit_vector(ncol(x))
    kappa_hat <- 0
  } else {
    mu_hat <- resultant / r_bar
    kappa_hat <- solve_vmf_kappa_from_rbar(r_bar, q = ncol(x) - 1L)
  }

  normalize_vmf_theta(
    list(mu = mu_hat, kappa = kappa_hat),
    ambient_dim = ncol(x)
  )
}

prepare_vmf_fast_multiplier <- function(data,
                                        theta_hat,
                                        spec,
                                        ks_prep = NULL,
                                        cvm_prep = NULL,
                                        control = list(),
                                        distance_type = "geodesic") {
  x <- normalize_vmf_data(data, control)
  theta_hat <- normalize_vmf_theta(theta_hat, ambient_dim = ncol(x))
  p <- length(theta_hat$xi)
  q <- theta_hat$q
  derivative_control <- fast_multiplier_parse_derivative_control(control)

  S_obs <- t(vapply(seq_len(nrow(x)), function(i) {
    psi_xi(x[i, ], theta_hat$xi, q)
  }, numeric(p)))
  Vhat <- -dot_psi_xi(theta_hat$xi, q)

  if (isTRUE(any(!is.finite(Vhat)))) {
    stop("The vMF fast multiplier preparation produced a non-finite `Vhat`.")
  }

  if (!is.null(derivative_control$derivative_mc_seed)) {
    set.seed(derivative_control$derivative_mc_seed)
  }
  aux_sample <- rotasym::r_vMF(
    derivative_control$derivative_mc_size,
    mu = theta_hat$mu,
    kappa = theta_hat$kappa
  )
  aux_sample <- normalize_vmf_data(aux_sample, control)
  Psi_aux <- t(vapply(seq_len(nrow(aux_sample)), function(i) {
    psi_xi(aux_sample[i, ], theta_hat$xi, q)
  }, numeric(p)))

  D_ks <- fast_multiplier_compute_D_ks(
    spec = spec,
    aux_sample = aux_sample,
    Psi_aux = Psi_aux,
    ks_prep = ks_prep,
    control = control
  )
  D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
  if (is.null(D_cvm)) {
    D_cvm <- fast_multiplier_compute_D_cvm(
      spec = spec,
      aux_sample = aux_sample,
      Psi_aux = Psi_aux,
      data = x,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = Psi_aux,
    derivative_method = derivative_control$derivative_method,
    derivative_mc_size = derivative_control$derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_control$derivative_mc_seed)) NA_integer_ else derivative_control$derivative_mc_seed,
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

make_vmf_spec <- function(distance_type = c("chordal", "geodesic"),
                          unknown_param = "xi") {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("vmf_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_vmf_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_vmf_data(data, control)
      omega_matrix <- normalize_vmf_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- x %*% t(omega_matrix)
      dot_products <- pmax(pmin(dot_products, 1), -1)

      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_vmf_theta(theta)
      ambient_dim <- length(theta$mu)
      profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
      l_max <- control$vmf_profile_l_max %||% NULL
      tail_tol <- as.numeric(control$vmf_profile_legendre_tail_tol %||% 1e-10)

      if (ambient_dim == 3L) {
        if (identical(profile_method, "legendre")) {
          return(distance_profile_vmf_s2_legendre(
            omega = omega,
            mu = theta$mu,
            kappa = theta$kappa,
            t_values = as.numeric(t),
            distance_type = distance_type,
            l_max = l_max,
            tail_tol = tail_tol
          ))
        }

        return(theoretical_distance_profile_vmf_s2_fast(
          omega = omega,
          mu = theta$mu,
          kappa = theta$kappa,
          t_values = as.numeric(t),
          distance_type = distance_type
        ))
      }

      theoretical_distance_profile_vmf(
        omega = omega,
        mu = theta$mu,
        kappa = theta$kappa,
        t_values = as.numeric(t),
        distance_type = distance_type
      )
    },
    normalize_data = normalize_vmf_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_vmf_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_vmf_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_vmf_theta(theta)
        ambient_dim <- length(theta$mu)
        profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
        n_u <- as.integer(control$vmf_profile_n_u %||% 4097L)
        l_max <- control$vmf_profile_l_max %||% NULL
        tail_tol <- as.numeric(control$vmf_profile_legendre_tail_tol %||% 1e-10)

        if (ambient_dim == 3L && identical(profile_method, "tabulated")) {
          return(distance_profile_vmf_s2_grid(
            omega_grid = omega_grid,
            mu = theta$mu,
            kappa = theta$kappa,
            t_grid = t_grid,
            distance_type = distance_type,
            n_u = n_u
          ))
        }
        if (ambient_dim == 3L && identical(profile_method, "legendre")) {
          return(distance_profile_vmf_s2_legendre_grid(
            omega_grid = omega_grid,
            mu = theta$mu,
            kappa = theta$kappa,
            t_grid = t_grid,
            distance_type = distance_type,
            l_max = l_max,
            tail_tol = tail_tol
          ))
        }

        NULL
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_vmf_theta(theta)
        ambient_dim <- length(theta$mu)
        profile_method <- tolower(as.character(control$vmf_profile_method %||% "tabulated"))
        n_u <- as.integer(control$vmf_profile_n_u %||% 4097L)
        l_max <- control$vmf_profile_l_max %||% NULL
        tail_tol <- as.numeric(control$vmf_profile_legendre_tail_tol %||% 1e-10)

        if (ambient_dim == 3L && identical(profile_method, "tabulated")) {
          return(distance_profile_vmf_s2_cvm_grid(
            X = data,
            mu = theta$mu,
            kappa = theta$kappa,
            n_u = n_u
          ))
        }
        if (ambient_dim == 3L && identical(profile_method, "legendre")) {
          return(distance_profile_vmf_s2_legendre_cvm_grid(
            X = data,
            mu = theta$mu,
            kappa = theta$kappa,
            l_max = l_max,
            tail_tol = tail_tol
          ))
        }

        NULL
      },
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        if (!identical(unknown_param, "xi")) {
          stop("The fast vMF multiplier bootstrap currently supports only `unknown_param = 'xi'`.")
        }

        prepare_vmf_fast_multiplier(
          data = data,
          theta_hat = theta_hat,
          spec = make_vmf_spec(distance_type = distance_type, unknown_param = unknown_param),
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          distance_type = distance_type
        )
      },
      distance_type = distance_type,
      unknown_param = unknown_param
    )
  )
}
