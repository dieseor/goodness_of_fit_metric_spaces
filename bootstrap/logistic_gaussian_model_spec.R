# Logistic-Gaussian simplex model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_logistic_gaussian <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_logistic_gaussian <- model_specs_candidates_logistic_gaussian[file.exists(model_specs_candidates_logistic_gaussian)][1L]
  if (is.na(model_specs_path_logistic_gaussian)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_logistic_gaussian)
}

logistic_gaussian_ilr_basis <- function(ambient_dim) {
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 2L) {
    stop("`ambient_dim` must be an integer greater than or equal to 2.")
  }

  basis <- matrix(0, nrow = ambient_dim, ncol = ambient_dim - 1L)
  for (j in seq_len(ambient_dim - 1L)) {
    basis[seq_len(j), j] <- 1 / sqrt(j * (j + 1))
    basis[j + 1L, j] <- -sqrt(j / (j + 1))
  }

  basis
}

row_softmax_matrix <- function(values) {
  values <- as.matrix(values)
  row_max <- apply(values, 1L, max)
  shifted <- sweep(values, 1L, row_max, FUN = "-")
  exp_shifted <- exp(shifted)
  exp_shifted / rowSums(exp_shifted)
}

coerce_logistic_gaussian_simplex_matrix <- function(data) {
  if (inherits(data, "logistic_gaussian_simplex_data")) {
    return(data$simplex)
  }

  if (is.list(data) && is.null(dim(data))) {
    if (length(data) == 0L) {
      stop("`data` cannot be empty.")
    }
    data <- do.call(rbind, lapply(data, as.numeric))
  } else if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (!is.matrix(data) || nrow(data) == 0L || ncol(data) < 2L) {
    stop("Logistic Gaussian data must be a non-empty matrix with at least two columns.")
  }

  data
}

normalize_logistic_gaussian_data <- function(data, control = list()) {
  if (inherits(data, "logistic_gaussian_simplex_data")) {
    return(data)
  }

  simplex <- coerce_logistic_gaussian_simplex_matrix(data)
  simplex <- matrix(
    as.numeric(simplex),
    nrow = nrow(simplex),
    ncol = ncol(simplex)
  )

  if (any(!is.finite(simplex))) {
    stop("Logistic Gaussian data must be finite.")
  }
  if (any(simplex <= 0)) {
    stop("Logistic Gaussian data must lie in the interior of the simplex.")
  }

  row_sums <- rowSums(simplex)
  tol <- as.numeric(control$logistic_gaussian_simplex_tol %||% 1e-8)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`control$logistic_gaussian_simplex_tol` must be a strictly positive finite scalar.")
  }
  if (any(abs(row_sums - 1) > tol)) {
    stop("Logistic Gaussian data rows must sum to 1 up to tolerance.")
  }

  simplex <- simplex / row_sums
  log_simplex <- log(simplex)
  clr_matrix <- log_simplex - rowMeans(log_simplex)
  basis <- logistic_gaussian_ilr_basis(ncol(simplex))
  ilr_matrix <- clr_matrix %*% basis

  structure(
    list(
      simplex = simplex,
      ilr = ilr_matrix,
      basis = basis,
      ambient_dim = ncol(simplex),
      ilr_dim = ncol(simplex) - 1L
    ),
    class = c("logistic_gaussian_simplex_data", "list")
  )
}

logistic_gaussian_simplex_matrix <- function(data, control = list()) {
  normalize_logistic_gaussian_data(data, control)$simplex
}

logistic_gaussian_ilr_matrix <- function(data, control = list()) {
  normalize_logistic_gaussian_data(data, control)$ilr
}

logistic_gaussian_point_to_ilr <- function(point, ambient_dim = NULL, control = list()) {
  normalized <- normalize_logistic_gaussian_data(point, control)
  if (!is.null(ambient_dim) && normalized$ambient_dim != as.integer(ambient_dim)) {
    stop("Logistic Gaussian point has incompatible ambient dimension.")
  }
  as.numeric(normalized$ilr[1L, , drop = TRUE])
}

logistic_gaussian_ilr_to_simplex <- function(z, ambient_dim = NULL) {
  if (is.vector(z)) {
    z <- matrix(as.numeric(z), nrow = 1L)
  } else {
    z <- as.matrix(z)
  }

  if (nrow(z) == 0L || ncol(z) == 0L || any(!is.finite(z))) {
    stop("`z` must be a non-empty finite matrix.")
  }

  ilr_dim <- ncol(z)
  ambient_dim <- ambient_dim %||% (ilr_dim + 1L)
  ambient_dim <- as.integer(ambient_dim)
  if (ambient_dim != ilr_dim + 1L) {
    stop("`ambient_dim` is incompatible with the number of ilr coordinates.")
  }

  basis <- logistic_gaussian_ilr_basis(ambient_dim)
  clr_matrix <- z %*% t(basis)
  row_softmax_matrix(clr_matrix)
}

rlogistic_gaussian_simplex <- function(n, mu_ilr, Sigma_ilr, control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- normalize_logistic_gaussian_theta(
    list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    control = control
  )
  z <- mvtnorm::rmvnorm(n = n, mean = theta$mu_ilr, sigma = theta$Sigma_ilr)
  logistic_gaussian_ilr_to_simplex(z, ambient_dim = theta$ambient_dim)
}

normalize_logistic_gaussian_theta <- function(theta,
                                              ambient_dim = NULL,
                                              control = list()) {
  if (!is.list(theta)) {
    stop("Logistic Gaussian theta must be a list.")
  }

  mu_ilr <- as.numeric(theta$mu_ilr %||% theta$mu)
  Sigma_ilr <- theta$Sigma_ilr %||% theta$Sigma

  if (length(mu_ilr) == 0L || any(!is.finite(mu_ilr))) {
    stop("Logistic Gaussian theta requires a finite vector `mu_ilr`.")
  }

  ambient_dim <- ambient_dim %||% theta$ambient_dim %||% (length(mu_ilr) + 1L)
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 2L) {
    stop("`ambient_dim` must be an integer greater than or equal to 2.")
  }
  if (length(mu_ilr) != ambient_dim - 1L) {
    stop("Logistic Gaussian theta has incompatible ilr dimension.")
  }

  Sigma_ilr <- as.matrix(Sigma_ilr)
  if (!all(dim(Sigma_ilr) == c(length(mu_ilr), length(mu_ilr)))) {
    stop("Logistic Gaussian theta requires a square covariance matrix `Sigma_ilr`.")
  }
  if (any(!is.finite(Sigma_ilr))) {
    stop("Logistic Gaussian covariance must be finite.")
  }

  Sigma_ilr <- 0.5 * (Sigma_ilr + t(Sigma_ilr))
  tol <- as.numeric(control$logistic_gaussian_cov_tol %||% 1e-10)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`control$logistic_gaussian_cov_tol` must be a strictly positive finite scalar.")
  }

  eigen_sigma <- eigen(Sigma_ilr, symmetric = TRUE)
  if (any(eigen_sigma$values < -tol)) {
    stop("Logistic Gaussian covariance must be positive semidefinite.")
  }

  eigenvalues_full <- pmax(eigen_sigma$values, 0)
  positive_idx <- eigenvalues_full > tol
  basis <- theta$basis %||% logistic_gaussian_ilr_basis(ambient_dim)

  list(
    mu_ilr = mu_ilr,
    Sigma_ilr = Sigma_ilr,
    basis = basis,
    ambient_dim = ambient_dim,
    ilr_dim = ambient_dim - 1L,
    eigenvalues_full = eigenvalues_full,
    eigenvectors_full = eigen_sigma$vectors,
    positive_idx = positive_idx,
    rank = sum(positive_idx)
  )
}

logistic_gaussian_weighted_covariance <- function(z, prob_weights, center) {
  centered <- sweep(z, 2L, center, FUN = "-")
  crossprod(centered * sqrt(prob_weights))
}

fit_logistic_gaussian_theta <- function(data,
                                        weights = NULL,
                                        null,
                                        unknown_param = "both",
                                        control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_logistic_gaussian_theta(
      null$theta,
      ambient_dim = normalized$ambient_dim,
      control = control
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  unknown_param <- tolower(as.character(unknown_param %||% "both"))
  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the Logistic Gaussian model.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(z), nrow(z))
  } else {
    normalize_probability_weights(weights, nrow(z))
  }

  fixed <- null$fixed %||% list()
  fixed_mu_raw <- fixed$mu_ilr %||% fixed$mu
  fixed_sigma_raw <- fixed$Sigma_ilr %||% fixed$Sigma
  weighted_mean <- colSums(z * prob_weights)

  if (identical(unknown_param, "mu")) {
    if (is.null(fixed_sigma_raw)) {
      stop("Composite Logistic Gaussian null with unknown `mu` requires `null$fixed$Sigma_ilr`.")
    }
    theta_out <- list(
      mu_ilr = weighted_mean,
      Sigma_ilr = normalize_logistic_gaussian_theta(
        list(mu_ilr = rep.int(0, normalized$ilr_dim), Sigma_ilr = fixed_sigma_raw),
        ambient_dim = normalized$ambient_dim,
        control = control
      )$Sigma_ilr
    )
    return(normalize_logistic_gaussian_theta(theta_out, ambient_dim = normalized$ambient_dim, control = control))
  }

  if (identical(unknown_param, "sigma")) {
    if (is.null(fixed_mu_raw)) {
      stop("Composite Logistic Gaussian null with unknown `Sigma` requires `null$fixed$mu_ilr`.")
    }
    fixed_mu <- normalize_logistic_gaussian_theta(
      list(mu_ilr = fixed_mu_raw, Sigma_ilr = diag(normalized$ilr_dim)),
      ambient_dim = normalized$ambient_dim,
      control = control
    )$mu_ilr
    theta_out <- list(
      mu_ilr = fixed_mu,
      Sigma_ilr = logistic_gaussian_weighted_covariance(z, prob_weights, fixed_mu)
    )
    return(normalize_logistic_gaussian_theta(theta_out, ambient_dim = normalized$ambient_dim, control = control))
  }

  theta_out <- list(
    mu_ilr = weighted_mean,
    Sigma_ilr = logistic_gaussian_weighted_covariance(z, prob_weights, weighted_mean)
  )
  normalize_logistic_gaussian_theta(theta_out, ambient_dim = normalized$ambient_dim, control = control)
}

prepare_logistic_gaussian_fast_multiplier <- function(spec,
                                                      data,
                                                      theta_hat,
                                                      ks_prep = NULL,
                                                      cvm_prep = NULL,
                                                      control = list(),
                                                      unknown_param = "both") {
  normalized <- normalize_logistic_gaussian_data(data, control)
  theta_hat <- normalize_logistic_gaussian_theta(theta_hat, ambient_dim = normalized$ambient_dim, control = control)
  d <- theta_hat$ilr_dim
  sigma_vech_hat <- fast_multiplier_vech(theta_hat$Sigma_ilr)
  unknown_param <- tolower(as.character(unknown_param %||% "both"))

  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the fast Logistic Gaussian multiplier bootstrap.")
  }

  par0 <- switch(
    unknown_param,
    mu = theta_hat$mu_ilr,
    sigma = sigma_vech_hat,
    both = c(theta_hat$mu_ilr, sigma_vech_hat)
  )

  gaussian_score_matrix_from_ilr <- function(z, theta) {
    gaussian_score_matrix_vech(
      z, theta$mu_ilr, theta$Sigma_ilr, unknown_param = unknown_param
    )
  }

  gaussian_influence_matrix_from_ilr <- function(z, theta) {
    gaussian_mle_influence_matrix_vech(
      z, theta$mu_ilr, theta$Sigma_ilr, unknown_param = unknown_param
    )
  }

  legacy_mc_control <- !is.null(control$derivative_mc_size) ||
    !is.null(control$derivative_mc_seed)
  if (is.null(control$derivative_method) && legacy_mc_control) {
    warning(
      paste(
        "Logistic-Gaussian fast multiplier: `derivative_method` was not",
        "supplied, but legacy `derivative_mc_size`/`derivative_mc_seed`",
        "controls were found. Selecting `score_mc`; set",
        "`derivative_method = 'score_mc'` or 'quadrature' explicitly."
      ),
      call. = FALSE
    )
  }
  derivative_control <- fast_multiplier_parse_derivative_control(
    control,
    default_method = if (legacy_mc_control) "score_mc" else "quadrature"
  )
  z_obs <- normalized$ilr
  store_paper_quantities <- isTRUE(control$fast_multiplier_store_paper_quantities)
  paper_score_obs <- if (store_paper_quantities) {
    gaussian_score_matrix_from_ilr(z_obs, theta_hat)
  } else {
    NULL
  }
  paper_Vhat <- if (store_paper_quantities) {
    fast_multiplier_gaussian_paper_vhat(
      theta_hat$Sigma_ilr,
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
  paper_influence_obs <- gaussian_influence_matrix_from_ilr(z_obs, theta_hat)

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
    evaluator <- function(omega_ilr, thresholds) {
      result <- gaussian_ball_profile_quadrature(
        omega = omega_ilr,
        mu = theta_hat$mu_ilr,
        Sigma = theta_hat$Sigma_ilr,
        t_values = thresholds,
        control = control,
        spectral = spectral
      )
      result$derivative <- result$derivative[, parameter_index, drop = FALSE]
      result
    }
    table_result <- gaussian_fast_quadrature_tables(
      data_centers = z_obs,
      ks_centers = if (is.null(ks_prep)) NULL else
        normalize_logistic_gaussian_data(ks_prep$omega_grid, control)$ilr,
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
    information <- gaussian_fisher_information_vech(
      theta_hat$Sigma_ilr, unknown_param
    )
    quadrature_settings <- gaussian_quadrature_settings(control)
    return(list(
      S_obs = S_obs,
      Vhat = Vhat,
      Psi_aux = matrix(numeric(0), nrow = 0L, ncol = ncol(S_obs)),
      correction_representation = "fitted_mle_influence",
      paper_score_obs = gaussian_score_matrix_vech(
        z_obs, theta_hat$mu_ilr, theta_hat$Sigma_ilr, unknown_param
      ),
      paper_Vhat = -information,
      paper_Vhat_inverse = -solve(information),
      paper_Vhat_method = "analytic_expected_score_jacobian",
      paper_Vhat_diagnostics = fast_multiplier_matrix_condition_diagnostics(
        -information
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
  aux_sample <- rlogistic_gaussian_simplex(
    n = derivative_control$derivative_mc_size,
    mu_ilr = theta_hat$mu_ilr,
    Sigma_ilr = theta_hat$Sigma_ilr,
    control = control
  )
  z_aux <- logistic_gaussian_ilr_matrix(aux_sample, control)
  Psi_aux <- gaussian_score_matrix_from_ilr(z_aux, theta_hat)
  # `S_obs` and `Vhat` are legacy names required by the shared fast engine.
  # They represent the correction in fitted MLE-influence coordinates, rather
  # than psi_{theta_hat} and Vhat in the paper's notation.
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

make_logistic_gaussian_spec <- function(unknown_param = "both") {
  new_model_spec(
    name = "logistic_gaussian_aitchison",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_logistic_gaussian_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_logistic_gaussian_data(data, control)
      omega_normalized <- normalize_logistic_gaussian_data(omega, control)

      if (x$ambient_dim != omega_normalized$ambient_dim) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      x_sq <- rowSums(x$ilr^2)
      omega_sq <- rowSums(omega_normalized$ilr^2)
      sq_distances <- outer(x_sq, omega_sq, FUN = "+") - 2 * (x$ilr %*% t(omega_normalized$ilr))
      sqrt(pmax(sq_distances, 0))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_logistic_gaussian_theta(theta, control = control)
      omega_ilr <- logistic_gaussian_point_to_ilr(
        omega,
        ambient_dim = theta$ambient_dim,
        control = control
      )
      derivative_method <- tolower(as.character(
        control$derivative_method %||% ""
      ))
      if (derivative_method %in% c("quadrature", "auto", "deterministic")) {
        return(gaussian_ball_profile_quadrature(
          omega = omega_ilr,
          mu = theta$mu_ilr,
          Sigma = theta$Sigma_ilr,
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
        shift = theta$mu_ilr - omega_ilr,
        t_values = as.numeric(t),
        eigenvalues_full = theta$eigenvalues_full,
        eigenvectors_full = theta$eigenvectors_full,
        positive_idx = theta$positive_idx,
        control = mvnormal_quadform_with_label(control, "logistic-Gaussian ilr distance profile")
      )
    },
    normalize_data = normalize_logistic_gaussian_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_logistic_gaussian_data(data, control)$simplex)
    },
    observation_at = function(data, idx, control = list()) {
      normalize_logistic_gaussian_data(data, control)$simplex[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_logistic_gaussian_theta(theta, control = control)
        omega_normalized <- normalize_logistic_gaussian_data(omega_grid, control)
        if (omega_normalized$ambient_dim != theta$ambient_dim) {
          stop("`omega_grid` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu_ilr, each = nrow(omega_normalized$ilr)),
          nrow = nrow(omega_normalized$ilr),
          ncol = length(theta$mu_ilr)
        ) - omega_normalized$ilr
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
          return(t(vapply(seq_len(nrow(omega_normalized$ilr)), function(i) {
            gaussian_ball_profile_quadrature(
              omega = omega_normalized$ilr[i, ],
              mu = theta$mu_ilr,
              Sigma = theta$Sigma_ilr,
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
          control = mvnormal_quadform_with_label(control, "logistic-Gaussian ilr distance-profile grid")
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        theta <- normalize_logistic_gaussian_theta(theta, control = control)
        data_normalized <- normalize_logistic_gaussian_data(data, control)
        if (data_normalized$ambient_dim != theta$ambient_dim) {
          stop("`data` has incompatible ambient dimension.")
        }

        shift_matrix <- matrix(
          rep(theta$mu_ilr, each = nrow(data_normalized$ilr)),
          nrow = nrow(data_normalized$ilr),
          ncol = length(theta$mu_ilr)
        ) - data_normalized$ilr
        derivative_method <- tolower(as.character(
          control$derivative_method %||% ""
        ))
        if (derivative_method %in% c("quadrature", "auto", "deterministic")) {
          spectral <- list(
            values = theta$eigenvalues_full,
            vectors = theta$eigenvectors_full
          )
          return(t(vapply(seq_len(nrow(data_normalized$ilr)), function(i) {
            gaussian_ball_profile_quadrature(
              omega = data_normalized$ilr[i, ],
              mu = theta$mu_ilr,
              Sigma = theta$Sigma_ilr,
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
          control = mvnormal_quadform_with_label(control, "logistic-Gaussian ilr sample distance-profile grid")
        )
      },
      score_matrix = function(data, theta, control = list()) {
        theta <- normalize_logistic_gaussian_theta(theta, control = control)
        gaussian_score_matrix_vech(
          logistic_gaussian_ilr_matrix(data, control),
          theta$mu_ilr,
          theta$Sigma_ilr,
          unknown_param = unknown_param
        )
      },
      fisher_information = function(theta, control = list()) {
        theta <- normalize_logistic_gaussian_theta(theta, control = control)
        gaussian_fisher_information_vech(
          theta$Sigma_ilr, unknown_param = unknown_param
        )
      },
      unknown_param = unknown_param,
      distance_type = "aitchison",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_logistic_gaussian_fast_multiplier(
          spec = make_logistic_gaussian_spec(unknown_param = unknown_param),
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
