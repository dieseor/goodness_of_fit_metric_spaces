# Multivariate-normal model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_mvnormal <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_mvnormal <- model_specs_candidates_mvnormal[file.exists(model_specs_candidates_mvnormal)][1L]
  if (is.na(model_specs_path_mvnormal)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_mvnormal)
}

normalize_mvnormal_data <- function(data, control = list()) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (!is.matrix(data) || nrow(data) == 0L || ncol(data) == 0L) {
    stop("Multivariate normal data must be a non-empty matrix.")
  }

  data <- matrix(
    as.numeric(data),
    nrow = nrow(data),
    ncol = ncol(data)
  )

  if (any(!is.finite(data))) {
    stop("Multivariate normal data must be finite.")
  }

  data
}

normalize_mvnormal_theta <- function(theta,
                                     ambient_dim = NULL,
                                     control = list()) {
  if (!is.list(theta)) {
    stop("Multivariate normal theta must be a list containing `mu` and `Sigma`.")
  }

  mu <- as.numeric(theta$mu)
  Sigma <- theta$Sigma

  if (length(mu) == 0L || any(!is.finite(mu))) {
    stop("Multivariate normal theta requires a finite vector `mu`.")
  }

  ambient_dim <- ambient_dim %||% theta$ambient_dim %||% length(mu)
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 1L) {
    stop("`ambient_dim` must be a strictly positive integer.")
  }
  if (length(mu) != ambient_dim) {
    stop("Multivariate normal theta has incompatible dimension.")
  }

  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(ambient_dim, ambient_dim))) {
    stop("Multivariate normal theta requires a square covariance matrix `Sigma`.")
  }
  if (any(!is.finite(Sigma))) {
    stop("Multivariate normal covariance must be finite.")
  }

  Sigma <- 0.5 * (Sigma + t(Sigma))
  tol <- as.numeric(control$mvnormal_cov_tol %||% control$logistic_gaussian_cov_tol %||% 1e-10)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`control$mvnormal_cov_tol` must be a strictly positive finite scalar.")
  }

  eigen_sigma <- eigen(Sigma, symmetric = TRUE)
  if (any(eigen_sigma$values < -tol)) {
    stop("Multivariate normal covariance must be positive semidefinite.")
  }

  eigenvalues_full <- pmax(eigen_sigma$values, 0)
  positive_idx <- eigenvalues_full > tol

  list(
    mu = mu,
    Sigma = Sigma,
    ambient_dim = ambient_dim,
    eigenvalues_full = eigenvalues_full,
    eigenvectors_full = eigen_sigma$vectors,
    positive_idx = positive_idx,
    rank = sum(positive_idx)
  )
}

mvnormal_weighted_covariance <- function(x, prob_weights, center) {
  centered <- sweep(x, 2L, center, FUN = "-")
  crossprod(centered * sqrt(prob_weights))
}

rmvnormal_euclidean <- function(n, mu, Sigma, control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- normalize_mvnormal_theta(
    list(mu = mu, Sigma = Sigma),
    control = control
  )

  mvtnorm::rmvnorm(n = n, mean = theta$mu, sigma = theta$Sigma)
}

fit_mvnormal_theta <- function(data,
                               weights = NULL,
                               null,
                               unknown_param = "both",
                               control = list()) {
  x <- normalize_mvnormal_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_mvnormal_theta(
      null$theta,
      ambient_dim = ncol(x),
      control = control
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the multivariate normal model.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }

  fixed <- null$fixed %||% list()
  fixed_mu_raw <- fixed$mu
  fixed_sigma_raw <- fixed$Sigma
  weighted_mean <- colSums(x * prob_weights)

  if (identical(unknown_param, "mu")) {
    if (is.null(fixed_sigma_raw)) {
      stop("Composite multivariate normal null with unknown `mu` requires `null$fixed$Sigma`.")
    }
    theta_out <- list(
      mu = weighted_mean,
      Sigma = normalize_mvnormal_theta(
        list(mu = rep.int(0, ncol(x)), Sigma = fixed_sigma_raw),
        ambient_dim = ncol(x),
        control = control
      )$Sigma
    )
    return(normalize_mvnormal_theta(theta_out, ambient_dim = ncol(x), control = control))
  }

  if (identical(unknown_param, "sigma")) {
    if (is.null(fixed_mu_raw)) {
      stop("Composite multivariate normal null with unknown `Sigma` requires `null$fixed$mu`.")
    }
    fixed_mu <- normalize_mvnormal_theta(
      list(mu = fixed_mu_raw, Sigma = diag(ncol(x))),
      ambient_dim = ncol(x),
      control = control
    )$mu
    theta_out <- list(
      mu = fixed_mu,
      Sigma = mvnormal_weighted_covariance(x, prob_weights, fixed_mu)
    )
    return(normalize_mvnormal_theta(theta_out, ambient_dim = ncol(x), control = control))
  }

  theta_out <- list(
    mu = weighted_mean,
    Sigma = mvnormal_weighted_covariance(x, prob_weights, weighted_mean)
  )
  normalize_mvnormal_theta(theta_out, ambient_dim = ncol(x), control = control)
}

prepare_mvnormal_fast_multiplier <- function(spec,
                                             data,
                                             theta_hat,
                                             ks_prep = NULL,
                                             cvm_prep = NULL,
                                             control = list(),
                                             unknown_param = "both") {
  normalized <- normalize_mvnormal_data(data, control)
  theta_hat <- normalize_mvnormal_theta(theta_hat, ambient_dim = ncol(normalized), control = control)
  d <- theta_hat$ambient_dim
  sigma_vech_hat <- fast_multiplier_vech(theta_hat$Sigma)
  unknown_param <- tolower(as.character(unknown_param %||% "both"))

  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the fast multivariate normal multiplier bootstrap.")
  }

  par0 <- switch(
    unknown_param,
    mu = theta_hat$mu,
    sigma = sigma_vech_hat,
    both = c(theta_hat$mu, sigma_vech_hat)
  )

  gaussian_score_matrix_from_matrix <- function(x, theta) {
    gaussian_score_matrix_vech(
      x, theta$mu, theta$Sigma, unknown_param = unknown_param
    )
  }

  gaussian_influence_matrix_from_matrix <- function(x, theta) {
    gaussian_mle_influence_matrix_vech(
      x, theta$mu, theta$Sigma, unknown_param = unknown_param
    )
  }

  # Score-MC is the production default for the composite multivariate-normal
  # model. The quadrature branch below is retained only for explicit
  # reproducibility/diagnostic calls.
  derivative_control <- fast_multiplier_parse_derivative_control(
    control,
    default_method = "score_mc",
    default_mc_size = 10000L
  )
  store_paper_quantities <- isTRUE(control$fast_multiplier_store_paper_quantities)
  paper_score_obs <- if (store_paper_quantities) {
    gaussian_score_matrix_from_matrix(normalized, theta_hat)
  } else {
    NULL
  }
  paper_Vhat <- if (store_paper_quantities) {
    fast_multiplier_gaussian_paper_vhat(
      theta_hat$Sigma,
      unknown_param = unknown_param
    )
  } else {
    NULL
  }
  paper_Vhat_diagnostics <- if (store_paper_quantities) {
    fast_multiplier_matrix_condition_diagnostics(paper_Vhat)
  } else {
    NULL
  }
  paper_influence_obs <- gaussian_influence_matrix_from_matrix(normalized, theta_hat)

  observed_cvm_distance_matrix <- NULL
  if (!is.null(cvm_prep) && !isTRUE(cvm_prep$light)) {
    observed_cvm_distance_matrix <- cvm_prep$distance_matrix
    if (is.null(observed_cvm_distance_matrix)) {
      observed_cvm_distance_matrix <- spec$distance_matrix(normalized, normalized, control)
    }
    observed_cvm_distance_matrix <- as.matrix(observed_cvm_distance_matrix)
  }

  if (identical(derivative_control$derivative_method, "quadrature")) {
    parameter_index <- switch(
      unknown_param,
      mu = seq_len(d),
      sigma = d + seq_len(d * (d + 1L) / 2L),
      both = seq_len(d + d * (d + 1L) / 2L)
    )
    spectral <- list(
      values = theta_hat$eigenvalues_full,
      vectors = theta_hat$eigenvectors_full
    )
    evaluator <- function(omega, thresholds) {
      result <- gaussian_ball_profile_quadrature(
        omega = omega,
        mu = theta_hat$mu,
        Sigma = theta_hat$Sigma,
        t_values = thresholds,
        control = control,
        spectral = spectral
      )
      result$derivative <- result$derivative[, parameter_index, drop = FALSE]
      result
    }
    table_result <- gaussian_fast_quadrature_tables(
      data_centers = normalized,
      ks_centers = if (is.null(ks_prep)) NULL else
        normalize_mvnormal_data(ks_prep$omega_grid, control),
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      observed_distance_matrix = observed_cvm_distance_matrix,
      evaluator = evaluator
    )
    S_obs <- paper_influence_obs
    Vhat <- diag(ncol(S_obs))
    vhat_diagnostics <- fast_multiplier_deterministic_vhat_diagnostics(
      S_obs = S_obs,
      Vhat = Vhat,
      par0 = par0
    )
    quadrature_settings <- gaussian_quadrature_settings(control)
    return(list(
      S_obs = S_obs,
      Vhat = Vhat,
      Psi_aux = matrix(numeric(0), nrow = 0L, ncol = ncol(S_obs)),
      correction_representation = "fitted_mle_influence",
      paper_score_obs = gaussian_score_matrix_vech(
        normalized, theta_hat$mu, theta_hat$Sigma, unknown_param
      ),
      paper_Vhat = -gaussian_fisher_information_vech(
        theta_hat$Sigma, unknown_param
      ),
      paper_Vhat_inverse = -solve(gaussian_fisher_information_vech(
        theta_hat$Sigma, unknown_param
      )),
      paper_Vhat_method = "analytic_expected_score_jacobian",
      paper_Vhat_diagnostics = fast_multiplier_matrix_condition_diagnostics(
        -gaussian_fisher_information_vech(theta_hat$Sigma, unknown_param)
      ),
      paper_influence_obs = paper_influence_obs,
      vhat_method = "fitted_gaussian_influence_reparameterization",
      vhat_diagnostics = vhat_diagnostics,
      observed_cvm_distance_matrix = observed_cvm_distance_matrix,
      derivative_method_requested = derivative_control$derivative_method_requested,
      derivative_method_effective = "quadrature",
      derivative_method_selection_source =
        derivative_control$derivative_method_selection_source,
      derivative_method = "quadrature",
      derivative_mc_size = NA_integer_,
      derivative_mc_seed = NA_integer_,
      quadrature_settings = quadrature_settings,
      quadrature_diagnostics = table_result$diagnostics,
      D_ks = table_result$D_ks,
      D_cvm = table_result$D_cvm
    ))
  }

  if (!is.null(derivative_control$derivative_mc_seed)) {
    set.seed(derivative_control$derivative_mc_seed)
  }
  aux_sample <- rmvnormal_euclidean(
    n = derivative_control$derivative_mc_size,
    mu = theta_hat$mu,
    Sigma = theta_hat$Sigma,
    control = control
  )
  Psi_aux <- gaussian_score_matrix_from_matrix(aux_sample, theta_hat)
  # The generic fast engine predates the paper notation.  Its `S_obs` and
  # `Vhat` fields encode the correction in influence coordinates.  They are
  # not respectively psi_{theta_hat} and Vhat from the paper: explicitly,
  # S_obs = - paper_Vhat^{-1} paper_score_obs and the internal metric is I.
  S_obs <- paper_influence_obs
  Vhat <- diag(ncol(S_obs))
  vhat_diagnostics <- fast_multiplier_vhat_diagnostics(
    S_obs = S_obs,
    Psi_aux = Psi_aux,
    Vhat = Vhat,
    par0 = par0
  )

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
      data = normalized,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  list(
    S_obs = S_obs,
    Vhat = Vhat,
    Psi_aux = Psi_aux,
    correction_representation = "fitted_mle_influence",
    paper_score_obs = paper_score_obs,
    paper_Vhat = paper_Vhat,
    paper_Vhat_inverse = if (store_paper_quantities) solve(paper_Vhat) else NULL,
    paper_Vhat_method = "analytic_expected_score_jacobian",
    paper_Vhat_diagnostics = paper_Vhat_diagnostics,
    paper_influence_obs = paper_influence_obs,
    vhat_method = "fitted_gaussian_influence_reparameterization",
    vhat_diagnostics = vhat_diagnostics,
    observed_cvm_distance_matrix = observed_cvm_distance_matrix,
    derivative_method_requested = derivative_control$derivative_method_requested,
    derivative_method_effective = derivative_control$derivative_method_effective,
    derivative_method_selection_source =
      derivative_control$derivative_method_selection_source,
    derivative_method = derivative_control$derivative_method,
    derivative_mc_size = derivative_control$derivative_mc_size,
    derivative_mc_seed = if (is.null(derivative_control$derivative_mc_seed)) NA_integer_ else derivative_control$derivative_mc_seed,
    D_ks = D_ks,
    D_cvm = D_cvm
  )
}

make_mvnormal_spec <- function(unknown_param = "both") {
  new_model_spec(
    name = "mvnormal_euclidean",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_mvnormal_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_mvnormal_data(data, control)
      omega_matrix <- normalize_mvnormal_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      x_sq <- rowSums(x^2)
      omega_sq <- rowSums(omega_matrix^2)
      sq_distances <- outer(x_sq, omega_sq, FUN = "+") - 2 * (x %*% t(omega_matrix))
      sqrt(pmax(sq_distances, 0))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_mvnormal_theta(theta, control = control)
      omega_vec <- as.numeric(normalize_mvnormal_data(omega, control)[1L, , drop = TRUE])
      derivative_method <- tolower(as.character(
        control$derivative_method %||% ""
      ))
      if (derivative_method %in% c("quadrature", "auto", "deterministic")) {
        return(gaussian_ball_profile_quadrature(
          omega = omega_vec,
          mu = theta$mu,
          Sigma = theta$Sigma,
          t_values = as.numeric(t),
          control = control,
          spectral = list(
            values = theta$eigenvalues_full,
            vectors = theta$eigenvectors_full
          ),
          compute_derivative = FALSE
        )$F)
      }
      evaluate_mvnorm_distance_profile(
        shift = theta$mu - omega_vec,
        t_values = as.numeric(t),
        eigenvalues_full = theta$eigenvalues_full,
        eigenvectors_full = theta$eigenvectors_full,
        positive_idx = theta$positive_idx,
        control = mvnormal_quadform_with_label(control, "multivariate-normal distance profile")
      )
    },
    normalize_data = normalize_mvnormal_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_mvnormal_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_mvnormal_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_mvnormal_theta(theta, control = control)
        omega_matrix <- normalize_mvnormal_data(omega_grid, control)
        if (ncol(omega_matrix) != theta$ambient_dim) {
          stop("`omega_grid` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu, each = nrow(omega_matrix)),
          nrow = nrow(omega_matrix),
          ncol = length(theta$mu)
        ) - omega_matrix
        t_matrix <- matrix(
          rep(as.numeric(t_grid), each = nrow(shift_matrix)),
          nrow = nrow(shift_matrix),
          ncol = length(t_grid)
        )
        derivative_method <- tolower(as.character(
          control$derivative_method %||% ""
        ))
        if (derivative_method %in% c("quadrature", "auto", "deterministic")) {
          spectral <- list(
            values = theta$eigenvalues_full,
            vectors = theta$eigenvectors_full
          )
          return(t(vapply(seq_len(nrow(omega_matrix)), function(i) {
            gaussian_ball_profile_quadrature(
              omega = omega_matrix[i, ],
              mu = theta$mu,
              Sigma = theta$Sigma,
              t_values = t_matrix[i, ],
              control = control,
              spectral = spectral,
              compute_derivative = FALSE
            )$F
          }, numeric(ncol(t_matrix)))))
        }
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = t_matrix,
          eigenvalues_full = theta$eigenvalues_full,
          eigenvectors_full = theta$eigenvectors_full,
          positive_idx = theta$positive_idx,
          control = mvnormal_quadform_with_label(control, "multivariate-normal distance-profile grid")
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_mvnormal_theta(theta, control = control)
        x <- normalize_mvnormal_data(data, control)
        if (ncol(x) != theta$ambient_dim) {
          stop("`data` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu, each = nrow(x)),
          nrow = nrow(x),
          ncol = length(theta$mu)
        ) - x
        derivative_method <- tolower(as.character(
          control$derivative_method %||% ""
        ))
        if (derivative_method %in% c("quadrature", "auto", "deterministic")) {
          spectral <- list(
            values = theta$eigenvalues_full,
            vectors = theta$eigenvectors_full
          )
          return(t(vapply(seq_len(nrow(x)), function(i) {
            gaussian_ball_profile_quadrature(
              omega = x[i, ],
              mu = theta$mu,
              Sigma = theta$Sigma,
              t_values = distance_matrix[i, ],
              control = control,
              spectral = spectral,
              compute_derivative = FALSE
            )$F
          }, numeric(ncol(distance_matrix)))))
        }
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = distance_matrix,
          eigenvalues_full = theta$eigenvalues_full,
          eigenvectors_full = theta$eigenvectors_full,
          positive_idx = theta$positive_idx,
          control = mvnormal_quadform_with_label(control, "multivariate-normal sample distance-profile grid")
        )
      },
      score_matrix = function(data, theta, control = list()) {
        theta <- normalize_mvnormal_theta(theta, control = control)
        gaussian_score_matrix_vech(
          normalize_mvnormal_data(data, control),
          theta$mu,
          theta$Sigma,
          unknown_param = unknown_param
        )
      },
      fisher_information = function(theta, control = list()) {
        theta <- normalize_mvnormal_theta(theta, control = control)
        gaussian_fisher_information_vech(
          theta$Sigma, unknown_param = unknown_param
        )
      },
      unknown_param = unknown_param,
      distance_type = "euclidean",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_mvnormal_fast_multiplier(
          spec = make_mvnormal_spec(unknown_param = unknown_param),
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
