# Jones--Pewsey model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_jp <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_jp <- model_specs_candidates_jp[file.exists(model_specs_candidates_jp)][1L]
  if (is.na(model_specs_path_jp)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_jp)
}

normalize_jp_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) < 3L) {
    stop("JP data must be a non-empty matrix with at least three columns.")
  }
  if (any(!is.finite(data))) {
    stop("JP data must be finite.")
  }

  norms <- sqrt(rowSums(data^2))
  if (any(norms <= 0)) {
    stop("JP data rows must have strictly positive norm.")
  }

  data / norms
}

normalize_jp_theta <- function(theta, ambient_dim = NULL) {
  if (!is.list(theta)) {
    stop("JP theta must be a list containing `mu`, `kappa`, and `psi`.")
  }

  mu <- as.numeric(theta$mu)
  kappa <- as.numeric(theta$kappa)
  psi <- as.numeric(theta$psi)

  if (length(mu) < 3L || any(!is.finite(mu))) {
    stop("JP theta requires a finite vector `mu` of length at least 3.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("JP theta requires a nonnegative finite scalar `kappa`.")
  }
  if (length(psi) != 1L || !is.finite(psi)) {
    stop("JP theta requires a finite scalar `psi`.")
  }

  mu_norm <- sqrt(sum(mu^2))
  if (mu_norm <= 0) {
    stop("JP theta requires `mu` with strictly positive norm.")
  }
  mu <- mu / mu_norm

  if (!is.null(ambient_dim) && length(mu) != ambient_dim) {
    stop("JP theta has incompatible ambient dimension.")
  }

  q <- length(mu) - 1L
  if (q < 2L) {
    stop("JP theta currently supports only q >= 2.")
  }

  list(
    mu = mu,
    kappa = kappa,
    psi = psi,
    alpha = if (psi == 0) 0 else tanh(kappa * psi),
    beta = if (psi == 0) NA_real_ else 1 / psi,
    q = q
  )
}
fast_multiplier_jp_axial_expectations <- function(alpha,
                                                  beta,
                                                  quad_n = 1024L) {
  quad <- rotational_gauss_legendre(as.integer(quad_n))
  term <- 1 + alpha * quad$nodes
  if (any(!is.finite(term)) || any(term <= 0)) {
    stop("JP axial expectation quadrature encountered nonpositive support.")
  }

  weights <- quad$weights * term^beta
  denom <- sum(weights)
  if (!is.finite(denom) || denom <= 0) {
    stop("JP axial expectation quadrature produced a nonpositive normalizing constant.")
  }

  list(
    e_t_over_one_plus_alpha_t = sum(weights * (quad$nodes / term)) / denom,
    e_log_one_plus_alpha_t = sum(weights * log(term)) / denom
  )
}

prepare_jp_fast_multiplier <- function(spec,
                                       data,
                                       theta_hat,
                                       ks_prep = NULL,
                                       cvm_prep = NULL,
                                       control = list(),
                                       distance_type = "geodesic") {
  x <- normalize_jp_data(data, control)
  theta_hat <- normalize_jp_theta(theta_hat, ambient_dim = ncol(x))

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = length(theta_hat$mu),
    kappa = theta_hat$kappa,
    psi = theta_hat$psi,
    abs_kappa_psi_tol = as.numeric(control$jp_vmf_switch_abs_kappa_psi %||% jp_vmf_near_zero_abs_kappa_psi_default)
  )) {
    theta_vmf <- normalize_vmf_theta(list(mu = theta_hat$mu, kappa = theta_hat$kappa))
    return(prepare_vmf_fast_multiplier(
      data = x,
      theta_hat = theta_vmf,
      spec = make_vmf_spec(distance_type = distance_type, unknown_param = "xi"),
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control,
      distance_type = distance_type
    ))
  }

  jp_alpha_boundary_eps <- as.numeric(control$jp_fast_alpha_boundary_eps %||% 1e-8)
  if (abs(theta_hat$alpha) >= 1 - jp_alpha_boundary_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "jp_alpha_boundary",
      derivative_method = NA_character_,
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_
    ))
  }

  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(0, 0, theta_hat$kappa, theta_hat$psi)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_jp_theta(
      list(mu = mapped$mu, kappa = par[[3L]], psi = par[[4L]]),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_jp_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    alpha <- theta_state$alpha
    beta <- theta_state$beta
    term <- 1 + alpha * z
    axial <- fast_multiplier_jp_axial_expectations(
      alpha = alpha,
      beta = beta,
      quad_n = as.integer(control$jp_fast_axial_quad_n %||% 1024L)
    )

    coeff_mu <- alpha * beta / term
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))
    psi_alpha <- beta * (z / term - axial$e_t_over_one_plus_alpha_t)
    psi_beta <- log(term) - axial$e_log_one_plus_alpha_t

    cbind(
      score_mu,
      theta_state$psi * (1 - alpha^2) * psi_alpha,
      theta_state$kappa * (1 - alpha^2) * psi_alpha - (theta_state$psi^-2) * psi_beta
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_jp(
      n = n_aux,
      mu = theta_state$mu,
      kappa = theta_state$kappa,
      psi = theta_state$psi
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
fit_jp_theta <- function(data,
                         weights = NULL,
                         null,
                         control = list()) {
  x <- normalize_jp_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_jp_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (ncol(x) != 3L) {
    stop("The JP composite adapter currently supports only S^2, i.e. data with three columns.")
  }

  jp_control <- control
  jp_control$jp_data_already_normalized <- TRUE

  normalize_jp_theta(
    jp_mle_s2_weighted(
      data = x,
      weights = weights,
      control = jp_control
    ),
    ambient_dim = ncol(x)
  )
}

make_jp_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("jp_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_jp_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_jp_data(data, control)
      omega_matrix <- normalize_jp_data(omega, control)

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
      theta <- normalize_jp_theta(theta)
      distance_profile_jp(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        kappa = theta$kappa,
        psi = theta$psi,
        distance_type = distance_type
      )
    },
    normalize_data = normalize_jp_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_jp_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_jp_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_jp_theta(theta)
        profile_method <- tolower(as.character(control$jp_profile_method %||% "tabulated"))
        n_u <- as.integer(control$jp_profile_n_u %||% 1025L)
        n_delta <- as.integer(control$jp_profile_n_delta %||% 257L)

        if (!identical(profile_method, "tabulated")) {
          return(NULL)
        }
        if (theta$psi == 0 && length(theta$mu) != 3L) {
          return(NULL)
        }

        distance_profile_jp_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          kappa = theta$kappa,
          psi = theta$psi,
          t_grid = t_grid,
          distance_type = distance_type,
          n_u = n_u,
          n_delta = n_delta
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_jp_theta(theta)
        profile_method <- tolower(as.character(control$jp_profile_method %||% "tabulated"))
        n_u <- as.integer(control$jp_profile_n_u %||% 1025L)
        n_delta <- as.integer(control$jp_profile_n_delta %||% 257L)

        if (!identical(profile_method, "tabulated")) {
          return(NULL)
        }
        if (theta$psi == 0 && length(theta$mu) != 3L) {
          return(NULL)
        }

        distance_profile_jp_cvm_grid(
          X = data,
          mu = theta$mu,
          kappa = theta$kappa,
          psi = theta$psi,
          n_u = n_u,
          n_delta = n_delta
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_jp_fast_multiplier(
          spec = make_jp_spec(distance_type = distance_type),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          distance_type = distance_type
        )
      }
    )
  )
}
