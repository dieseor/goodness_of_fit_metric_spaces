#!/usr/bin/env Rscript

# Exploratory known-parameter population-signal scan for dimension-scaled
# Section 6 alternatives.  This script is intentionally separate from the
# production runner: it writes only diagnostic Monte Carlo summaries.
#
# Every candidate is governed by one scalar that is held fixed across d.
# Dimension enters through a stated normalization, rather than through a
# separately tuned parameter at each d.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source("scripts/estimate_section6_population_signal.R")

candidate_bulk_vector <- function(d) c(0, rep(1, d - 1L))
candidate_unit_e <- function(d, index = 1L) section6_e(d, index)
candidate_transverse_indices <- function(d, support = c("all", "half")) {
  support <- match.arg(support)
  available <- seq.int(2L, d)
  if (support == "all") return(available)
  available[seq_len(ceiling(length(available) / 2))]
}

candidate_catalog <- function() {
  rows <- list()
  add <- function(id, family, rule, value, description) {
    rows[[length(rows) + 1L]] <<- data.frame(
      candidate = id, family = family, rule = rule, value = value,
      description = description, stringsAsFactors = FALSE
    )
  }

  # These four candidates implement only intrinsic, fixed-size changes.  No
  # exponent is selected from the numerical output: the mean vector has fixed
  # Euclidean norm and the covariance perturbation fixed Frobenius norm.  The
  # all/half distinction is precisely the number of transverse coordinates
  # affected, as requested for the initial screening.
  for (family in c("normal", "lg")) {
    prefix <- if (family == "normal") "" else "ILR "
    for (support in c("all", "half")) {
      add(
        sprintf("%s_known_mean_%s_l2", family, support), family,
        paste0("known_mean_", support, "_l2"), 0.5,
        sprintf("%sopposite means of Euclidean norm 0.5 on %s transverse coordinates", prefix, support)
      )
      add(
        sprintf("%s_known_scale_%s_frobenius", family, support), family,
        paste0("known_scale_", support, "_frobenius"), 0.75,
        sprintf("%ssymmetric diagonal covariance perturbation of Frobenius norm 0.75 on %s transverse coordinates", prefix, support)
      )
    }
  }
  add(
    "normal_known_joint_all", "normal", "known_joint_all", 0.5,
    "opposite means of norm 0.5 and all-transverse covariance perturbation of Frobenius norm 0.75"
  )
  add(
    "lg_known_joint_all", "lg", "known_joint_all", 0.5,
    "ILR opposite means of norm 0.5 and all-transverse covariance perturbation of Frobenius norm 0.75"
  )

  # A per-coordinate mean displacement on coordinates 2,...,d.  The factor
  # (d-1)^(-1/4) makes the squared separation O(sqrt(d)), the natural scale of
  # a Euclidean squared-distance fluctuation.
  for (a in c(0.3, 0.4, 0.5)) {
    add(sprintf("normal_bulk_mean_a%.1f", a), "normal", "bulk_mean", a,
      sprintf("means +/- a(d-1)^(-1/4)(0,1,...,1), a=%.1f", a))
    add(sprintf("lg_bulk_mean_a%.1f", a), "lg", "bulk_mean", a,
      sprintf("ILR means +/- a(d-1)^(-1/4)(0,1,...,1), a=%.1f", a))
  }
  for (a in c(0.55, 0.6)) {
    add(sprintf("normal_bulk_mean_a%.2f", a), "normal", "bulk_mean", a,
      sprintf("means +/- a(d-1)^(-1/4)(0,1,...,1), a=%.2f", a))
    add(sprintf("lg_bulk_mean_a%.2f", a), "lg", "bulk_mean", a,
      sprintf("ILR means +/- a(d-1)^(-1/4)(0,1,...,1), a=%.2f", a))
  }

  # The Frobenius norm of Sigma_+ - Sigma_- equals 2c for every d.  This
  # changes all d-1 variances transverse to e1 while retaining a fixed total
  # covariance perturbation.
  for (c_value in c(0.3, 0.5, 0.7)) {
    add(sprintf("normal_bulk_scale_c%.1f", c_value), "normal", "bulk_scale", c_value,
      sprintf("Sigma_+/- = I +/- c(d-1)^(-1/2) diag(0,1,...,1), c=%.1f", c_value))
    add(sprintf("lg_bulk_scale_c%.1f", c_value), "lg", "bulk_scale", c_value,
      sprintf("ILR Sigma_+/- = I +/- c(d-1)^(-1/2) diag(0,1,...,1), c=%.1f", c_value))
  }
  for (c_value in c(0.8, 0.9)) {
    add(sprintf("normal_bulk_scale_c%.1f", c_value), "normal", "bulk_scale", c_value,
      sprintf("Sigma_+/- = I +/- c(d-1)^(-1/2) diag(0,1,...,1), c=%.1f", c_value))
    add(sprintf("lg_bulk_scale_c%.1f", c_value), "lg", "bulk_scale", c_value,
      sprintf("ILR Sigma_+/- = I +/- c(d-1)^(-1/2) diag(0,1,...,1), c=%.1f", c_value))
  }
  # The same all-transverse-coordinate covariance mixture with a slower
  # d^{-1/4} per-coordinate scale.  This is included because the composite
  # normal fit absorbs the second-order covariance change under d^{-1/2}.
  for (c_value in c(0.5, 0.6, 0.7, 0.8)) {
    add(sprintf("normal_bulk_scale_dquarter_c%.1f", c_value), "normal", "bulk_scale_dquarter", c_value,
      sprintf("Sigma_+/- = I +/- c(d-1)^(-1/4) diag(0,1,...,1), c=%.1f", c_value))
    add(sprintf("lg_bulk_scale_dquarter_c%.1f", c_value), "lg", "bulk_scale_dquarter", c_value,
      sprintf("ILR Sigma_+/- = I +/- c(d-1)^(-1/4) diag(0,1,...,1), c=%.1f", c_value))
  }
  # Intermediate exponent selected for the composite diagnostic: it lies
  # between constant-Frobenius (d^{-1/2}) and the d^{-1/4} screen.
  for (c_value in c(0.75, 0.9)) {
    add(sprintf("normal_bulk_scale_dthird_c%.2f", c_value), "normal", "bulk_scale_dthird", c_value,
      sprintf("Sigma_+/- = I +/- c(d-1)^(-1/3) diag(0,1,...,1), c=%.2f", c_value))
    add(sprintf("lg_bulk_scale_dthird_c%.2f", c_value), "lg", "bulk_scale_dthird", c_value,
      sprintf("ILR Sigma_+/- = I +/- c(d-1)^(-1/3) diag(0,1,...,1), c=%.2f", c_value))
  }

  # nu_d = 2 + a d is the simple scale at which the common random scale of a
  # multivariate t becomes less variable as the radial dimension grows.
  for (a in c(0.5, 0.75, 1)) {
    add(sprintf("normal_t_nu_2_plus_%.2fd", a), "normal", "t_df_linear", a,
      sprintf("standardized multivariate t with nu_d = 2 + %.2f d", a))
    add(sprintf("lg_t_nu_2_plus_%.2fd", a), "lg", "t_df_linear", a,
      sprintf("ILR standardized multivariate t with nu_d = 2 + %.2f d", a))
  }
  for (a in c(0.15, 0.25, 0.35)) {
    add(sprintf("normal_t_from_t3_slope%.2f", a), "normal", "t_df_from_t3", a,
      sprintf("standardized multivariate t with nu_d = 3 + %.2f(d-2)", a))
    add(sprintf("lg_t_from_t3_slope%.2f", a), "lg", "t_df_from_t3", a,
      sprintf("ILR standardized multivariate t with nu_d = 3 + %.2f(d-2)", a))
  }

  # This keeps the t alternative but gives its null a genuinely different
  # (non-translational) covariance: diag(1,tau^2,...,tau^2).  Both laws share
  # this covariance, so any discrepancy is due to the t radial scale rather
  # than a covariance mismatch.
  for (tau2 in c(0.5, 1.5, 2.0)) {
    add(sprintf("normal_t_from_t3_tailvar%.1f", tau2), "normal", "t_df_from_t3_tailvar", tau2,
      sprintf("N(0,diag(1,%.1f,...,%.1f)) versus covariance-matched t, nu_d=3+.25(d-2)", tau2, tau2))
    add(sprintf("lg_t_from_t3_tailvar%.1f", tau2), "lg", "t_df_from_t3_tailvar", tau2,
      sprintf("ILR N(0,diag(1,%.1f,...,%.1f)) versus covariance-matched t, nu_d=3+.25(d-2)", tau2, tau2))
  }

  # The vMF concentration mixture changes the dispersion in every tangent
  # direction.  The null concentration is 1.5d and Q is an equal mixture at
  # relative concentrations 1-rho and 1+rho.
  for (rho in c(0.2, 0.4, 0.6)) {
    add(sprintf("vmf_concentration_mixture_rho%.1f", rho), "vmf", "concentration_mixture", rho,
      sprintf("equal mixture vMF(e1,1.5d(1+/-rho)), rho=%.1f", rho))
  }
  for (c_value in c(0.9, 1.0)) {
    add(sprintf("vmf_concentration_dquarter_c%.1f", c_value), "vmf", "concentration_dquarter", c_value,
      sprintf("equal mixture vMF(e1,1.5d(1+/-c d^(-1/4))), c=%.1f", c_value))
  }

  # A local location mixture is included as a benchmark.  With kappa=1.5d,
  # delta_d = c/sqrt(d) keeps sqrt(kappa) delta_d constant.
  for (c_value in c(0.5, 1, 1.5)) {
    add(sprintf("vmf_local_location_c%.1f", c_value), "vmf", "local_location", c_value,
      sprintf("vMF location mixture with angular displacement c/sqrt(d), c=%.1f", c_value))
  }
  for (c_value in c(0.9, 1.1)) {
    add(sprintf("vmf_local_location_dquarter_c%.1f", c_value), "vmf", "local_location_dquarter", c_value,
      sprintf("vMF location mixture with angular displacement c d^(-1/4), c=%.1f", c_value))
  }
  add(
    "vmf_joint_location_concentration_fisher", "vmf", "joint_location_concentration_fisher", 1,
    "vMF(e1,1.5d) mixed with vMF(mu_d,d), with delta_d=(pi/5)sqrt(2/d)"
  )
  add(
    "vmf_projected_normal", "vmf", "projected_normal", 1,
    "vMF(e1,1.5d) versus projected N(sqrt(d+1)e1,I_{d+1})"
  )
  for (c_value in c(0.75, 1.0, 1.25, 1.5)) {
    add(sprintf("vmf_symmetric_location_sqrt_c%.2f", c_value), "vmf", "symmetric_location_sqrt", c_value,
      sprintf("symmetric vMF location mixture with displacement c/sqrt(d), c=%.2f", c_value))
  }
  for (c_value in c(0.9, 1.1)) {
    add(sprintf("vmf_symmetric_location_dquarter_c%.1f", c_value), "vmf", "symmetric_location_dquarter", c_value,
      sprintf("symmetric vMF location mixture with displacement c d^(-1/4), c=%.1f", c_value))
  }
  for (c_value in c(0.75, 1.0, 1.25, 1.5)) {
    add(sprintf("hvmf_symmetric_location_sqrt_c%.2f", c_value), "hvmf", "symmetric_location_sqrt", c_value,
      sprintf("symmetric HvMF location mixture with displacement c/sqrt(d), c=%.2f", c_value))
  }

  # HvMF analogues are diagnostic only until the full H^d campaign is chosen.
  for (c_value in c(0.5, 1, 1.5)) {
    add(sprintf("hvmf_local_location_c%.1f", c_value), "hvmf", "local_location", c_value,
      sprintf("HvMF location mixture at hyperbolic displacement c/sqrt(d), c=%.1f", c_value))
  }
  for (c_value in c(0.9, 1.1)) {
    add(sprintf("hvmf_local_location_dquarter_c%.1f", c_value), "hvmf", "local_location_dquarter", c_value,
      sprintf("HvMF location mixture at hyperbolic displacement c d^(-1/4), c=%.1f", c_value))
  }
  for (c_value in c(0.9, 1.1)) {
    add(sprintf("hvmf_symmetric_location_dquarter_c%.1f", c_value), "hvmf", "symmetric_location_dquarter", c_value,
      sprintf("symmetric HvMF location mixture with displacement c d^(-1/4), c=%.1f", c_value))
  }
  for (c_value in c(0.5, 0.75, 1)) {
    add(sprintf("hvmf_local_angular_c%.2f", c_value), "hvmf", "local_angular", c_value,
      sprintf("HvMF angular mixture with displacement c/sqrt(d), c=%.2f", c_value))
  }
  for (c_value in c(0.9, 1.1)) {
    add(sprintf("hvmf_local_angular_dquarter_c%.1f", c_value), "hvmf", "local_angular_dquarter", c_value,
      sprintf("HvMF angular mixture with displacement c d^(-1/4), c=%.1f", c_value))
  }

  do.call(rbind, rows)
}

candidate_draw <- function(candidate, d, n, law = c("p0", "p1")) {
  law <- match.arg(law)
  rule <- as.character(candidate$rule)
  value <- as.numeric(candidate$value)
  d <- as.integer(d)
  n <- as.integer(n)

  if (candidate$family %in% c("normal", "lg")) {
    to_observation_space <- function(z) {
      if (candidate$family == "normal") return(z)
      logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L)
    }
    if (grepl("^known_mean_(all|half)_l2$", rule)) {
      support <- sub("^known_mean_(all|half)_l2$", "\\1", rule)
      active <- candidate_transverse_indices(d, support)
      direction <- numeric(d)
      direction[active] <- 1 / sqrt(length(active))
      mean_0 <- value * direction
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = mean_0, sigma = diag(d))))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = mean_0, sigma = diag(d))
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = -mean_0, sigma = diag(d))
      return(to_observation_space(x))
    }
    if (grepl("^known_scale_(all|half)_frobenius$", rule)) {
      support <- sub("^known_scale_(all|half)_frobenius$", "\\1", rule)
      active <- candidate_transverse_indices(d, support)
      diagonal_shift <- numeric(d)
      diagonal_shift[active] <- value / sqrt(length(active))
      sigma_plus <- diag(1 + diagonal_shift, d)
      sigma_minus <- diag(1 - diagonal_shift, d)
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = rep(0, d), sigma = sigma_minus)
      return(to_observation_space(x))
    }
    if (rule == "known_joint_all") {
      active <- candidate_transverse_indices(d, "all")
      direction <- numeric(d)
      direction[active] <- 1 / sqrt(length(active))
      mean_0 <- 0.5 * direction
      diagonal_shift <- numeric(d)
      diagonal_shift[active] <- 0.75 / sqrt(length(active))
      sigma_plus <- diag(1 + diagonal_shift, d)
      sigma_minus <- diag(1 - diagonal_shift, d)
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = mean_0, sigma = sigma_plus)))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = mean_0, sigma = sigma_plus)
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = -mean_0, sigma = sigma_minus)
      return(to_observation_space(x))
    }
    if (rule == "bulk_mean") {
      a_d <- value / (d - 1)^(1 / 4)
      mean_0 <- a_d * candidate_bulk_vector(d)
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = mean_0, sigma = diag(d))))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = mean_0, sigma = diag(d))
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = -mean_0, sigma = diag(d))
      return(to_observation_space(x))
    }
    if (rule == "bulk_scale") {
      rho_d <- value / sqrt(d - 1)
      diagonal_shift <- c(0, rep(rho_d, d - 1L))
      sigma_plus <- diag(1 + diagonal_shift, d)
      sigma_minus <- diag(1 - diagonal_shift, d)
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = rep(0, d), sigma = sigma_minus)
      return(to_observation_space(x))
    }
    if (rule == "bulk_scale_dquarter") {
      rho_d <- value / (d - 1)^(1 / 4)
      diagonal_shift <- c(0, rep(rho_d, d - 1L))
      sigma_plus <- diag(1 + diagonal_shift, d)
      sigma_minus <- diag(1 - diagonal_shift, d)
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = rep(0, d), sigma = sigma_minus)
      return(to_observation_space(x))
    }
    if (rule == "bulk_scale_dthird") {
      rho_d <- value / (d - 1)^(1 / 3)
      diagonal_shift <- c(0, rep(rho_d, d - 1L))
      sigma_plus <- diag(1 + diagonal_shift, d)
      sigma_minus <- diag(1 - diagonal_shift, d)
      if (law == "p0") return(to_observation_space(mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)))
      take_q <- stats::runif(n) < 0.5
      x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = sigma_plus)
      if (any(take_q)) x[take_q, ] <- mvtnorm::rmvnorm(sum(take_q), mean = rep(0, d), sigma = sigma_minus)
      return(to_observation_space(x))
    }
    if (rule == "t_df_linear") {
      if (law == "p0") return(to_observation_space(matrix(stats::rnorm(n * d), nrow = n, ncol = d)))
      return(to_observation_space(r_standardized_multivariate_t(n, d, nu = 2 + value * d)))
    }
    if (rule == "t_df_from_t3") {
      if (law == "p0") return(to_observation_space(matrix(stats::rnorm(n * d), nrow = n, ncol = d)))
      return(to_observation_space(r_standardized_multivariate_t(n, d, nu = 3 + value * (d - 2))))
    }
    if (rule == "t_df_from_t3_tailvar") {
      sigma_sqrt <- diag(c(1, rep(sqrt(value), d - 1L)), d)
      if (law == "p0") {
        return(to_observation_space(matrix(stats::rnorm(n * d), nrow = n, ncol = d) %*% sigma_sqrt))
      }
      z <- r_standardized_multivariate_t(n, d, nu = 3 + 0.25 * (d - 2)) %*% sigma_sqrt
      return(to_observation_space(z))
    }
  }

  if (candidate$family == "vmf") {
    mu0 <- candidate_unit_e(d + 1L)
    kappa <- 1.5 * d
    if (rule == "concentration_mixture") {
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa))
      take_high <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu0, kappa = kappa * (1 - value))
      if (any(take_high)) x[take_high, ] <- rotasym::r_vMF(sum(take_high), mu = mu0, kappa = kappa * (1 + value))
      return(x)
    }
    if (rule == "joint_location_concentration_fisher") {
      kappa_0 <- 1.5 * d
      delta <- (pi / 5) * sqrt(2 / d)
      mu_1 <- cos(delta) * mu0 + sin(delta) * candidate_unit_e(d + 1L, 2L)
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa_0))
      take_q <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu0, kappa = kappa_0)
      if (any(take_q)) x[take_q, ] <- rotasym::r_vMF(sum(take_q), mu = mu_1, kappa = d)
      return(x)
    }
    if (rule == "projected_normal") {
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = 1.5 * d))
      return(projected_normal_on_sphere(n, d))
    }
    if (rule == "concentration_dquarter") {
      rho_d <- value / d^(1 / 4)
      if (rho_d >= 1) stop("The selected concentration mixture has nonpositive low concentration.")
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa))
      take_high <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu0, kappa = kappa * (1 - rho_d))
      if (any(take_high)) x[take_high, ] <- rotasym::r_vMF(sum(take_high), mu = mu0, kappa = kappa * (1 + rho_d))
      return(x)
    }
    if (rule == "local_location") {
      delta <- value / sqrt(d)
      mu1 <- cos(delta) * mu0 + sin(delta) * candidate_unit_e(d + 1L, 2L)
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa))
      take_q <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu0, kappa = kappa)
      if (any(take_q)) x[take_q, ] <- rotasym::r_vMF(sum(take_q), mu = mu1, kappa = kappa)
      return(x)
    }
    if (rule == "local_location_dquarter") {
      delta <- value / d^(1 / 4)
      mu1 <- cos(delta) * mu0 + sin(delta) * candidate_unit_e(d + 1L, 2L)
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa))
      take_q <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu0, kappa = kappa)
      if (any(take_q)) x[take_q, ] <- rotasym::r_vMF(sum(take_q), mu = mu1, kappa = kappa)
      return(x)
    }
    if (rule == "symmetric_location_dquarter") {
      delta <- value / d^(1 / 4)
      mu_plus <- cos(delta) * mu0 + sin(delta) * candidate_unit_e(d + 1L, 2L)
      mu_minus <- cos(delta) * mu0 - sin(delta) * candidate_unit_e(d + 1L, 2L)
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa))
      take_plus <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu_minus, kappa = kappa)
      if (any(take_plus)) x[take_plus, ] <- rotasym::r_vMF(sum(take_plus), mu = mu_plus, kappa = kappa)
      return(x)
    }
    if (rule == "symmetric_location_sqrt") {
      delta <- value / sqrt(d)
      mu_plus <- cos(delta) * mu0 + sin(delta) * candidate_unit_e(d + 1L, 2L)
      mu_minus <- cos(delta) * mu0 - sin(delta) * candidate_unit_e(d + 1L, 2L)
      if (law == "p0") return(rotasym::r_vMF(n, mu = mu0, kappa = kappa))
      take_plus <- stats::runif(n) < 0.5
      x <- rotasym::r_vMF(n, mu = mu_minus, kappa = kappa)
      if (any(take_plus)) x[take_plus, ] <- rotasym::r_vMF(sum(take_plus), mu = mu_plus, kappa = kappa)
      return(x)
    }
  }

  if (candidate$family == "hvmf") {
    mu0 <- c(sqrt(2), candidate_unit_e(d))
    if (rule == "local_location") {
      r <- value / sqrt(d)
      tangent <- c(0, candidate_unit_e(d, 2L))
      mu1 <- cosh(r) * mu0 + sinh(r) * tangent
      if (law == "p0") return(rhvmf_polar(n, mu = mu0, kappa = d))
      take_q <- stats::runif(n) < 0.5
      x <- rhvmf_polar(n, mu = mu0, kappa = d)
      if (any(take_q)) x[take_q, ] <- rhvmf_polar(sum(take_q), mu = mu1, kappa = d)
      return(x)
    }
    if (rule == "local_location_dquarter") {
      r <- value / d^(1 / 4)
      tangent <- c(0, candidate_unit_e(d, 2L))
      mu1 <- cosh(r) * mu0 + sinh(r) * tangent
      if (law == "p0") return(rhvmf_polar(n, mu = mu0, kappa = d))
      take_q <- stats::runif(n) < 0.5
      x <- rhvmf_polar(n, mu = mu0, kappa = d)
      if (any(take_q)) x[take_q, ] <- rhvmf_polar(sum(take_q), mu = mu1, kappa = d)
      return(x)
    }
    if (rule == "local_angular") {
      kappa <- 1.5 * d
      if (law == "p0") return(rhvmf_polar(n, mu = mu0, kappa = kappa))
      return(rhvmf_angular_mixture(
        n, mu = mu0, kappa = kappa, delta = value / sqrt(d),
        tangent = candidate_unit_e(d, 2L)
      ))
    }
    if (rule == "local_angular_dquarter") {
      kappa <- 1.5 * d
      if (law == "p0") return(rhvmf_polar(n, mu = mu0, kappa = kappa))
      return(rhvmf_angular_mixture(
        n, mu = mu0, kappa = kappa, delta = value / d^(1 / 4),
        tangent = candidate_unit_e(d, 2L)
      ))
    }
    if (rule == "symmetric_location_dquarter") {
      kappa <- 1.5 * d
      r <- value / d^(1 / 4)
      tangent <- c(0, candidate_unit_e(d, 2L))
      mu_plus <- cosh(r) * mu0 + sinh(r) * tangent
      mu_minus <- cosh(r) * mu0 - sinh(r) * tangent
      if (law == "p0") return(rhvmf_polar(n, mu = mu0, kappa = kappa))
      take_plus <- stats::runif(n) < 0.5
      x <- rhvmf_polar(n, mu = mu_minus, kappa = kappa)
      if (any(take_plus)) x[take_plus, ] <- rhvmf_polar(sum(take_plus), mu = mu_plus, kappa = kappa)
      return(x)
    }
    if (rule == "symmetric_location_sqrt") {
      kappa <- 1.5 * d
      r <- value / sqrt(d)
      tangent <- c(0, candidate_unit_e(d, 2L))
      mu_plus <- cosh(r) * mu0 + sinh(r) * tangent
      mu_minus <- cosh(r) * mu0 - sinh(r) * tangent
      if (law == "p0") return(rhvmf_polar(n, mu = mu0, kappa = kappa))
      take_plus <- stats::runif(n) < 0.5
      x <- rhvmf_polar(n, mu = mu_minus, kappa = kappa)
      if (any(take_plus)) x[take_plus, ] <- rhvmf_polar(sum(take_plus), mu = mu_plus, kappa = kappa)
      return(x)
    }
  }
  stop(sprintf("Unsupported candidate '%s'.", candidate$candidate))
}

candidate_one_batch <- function(candidate, d, n_centres, n_profile, n_eval, seed) {
  set.seed(seed)
  centres <- candidate_draw(candidate, d, n_centres, "p0")
  x0_a <- candidate_draw(candidate, d, n_profile, "p0")
  x1_a <- candidate_draw(candidate, d, n_profile, "p1")
  x0_b <- candidate_draw(candidate, d, n_profile, "p0")
  x1_b <- candidate_draw(candidate, d, n_profile, "p1")
  x0_eval <- candidate_draw(candidate, d, n_eval, "p0")
  spec <- population_signal_spec(candidate$family)
  control <- section6_control(derivative_mc_size = 1000L, derivative_seed = seed, cvm_block_size = 50L)

  d0_a <- spec$distance_matrix(x0_a, centres, control)
  d1_a <- spec$distance_matrix(x1_a, centres, control)
  d0_b <- spec$distance_matrix(x0_b, centres, control)
  d1_b <- spec$distance_matrix(x1_b, centres, control)
  d0_eval <- spec$distance_matrix(x0_eval, centres, control)
  values <- t(vapply(seq_len(ncol(d0_a)), function(j) {
    centre_signal(d0_a[, j], d1_a[, j], d0_b[, j], d1_b[, j], d0_eval[, j])
  }, numeric(2)))
  data.frame(
    S_mean_center = mean(values[, "S"]), S_median_center = stats::median(values[, "S"]),
    S_max_center = max(values[, "S"]), J = mean(values[, "J"]), stringsAsFactors = FALSE
  )
}

candidate_summarize <- function(x) {
  group <- interaction(x$candidate, x$d, drop = TRUE)
  pieces <- lapply(split(x, group), function(cell) {
    first <- cell[1L, c("candidate", "family", "rule", "value", "description", "d"), drop = FALSE]
    for (metric in c("S_mean_center", "S_median_center", "S_max_center", "J")) {
      first[[metric]] <- mean(cell[[metric]])
      first[[paste0(metric, "_mcse")]] <- if (nrow(cell) > 1L) stats::sd(cell[[metric]]) / sqrt(nrow(cell)) else NA_real_
    }
    first
  })
  out <- do.call(rbind, pieces)
  out[order(out$family, out$candidate, out$d), , drop = FALSE]
}

run_candidate_scan <- function(dimensions = c(2L, 5L, 10L),
                               families = c("normal", "lg", "vmf", "hvmf"),
                               candidates = NULL,
                               n_centres = 96L,
                               n_profile = 6000L,
                               n_eval = 3000L,
                               batches = 3L,
                               cores = 1L,
                               seed = 20260729L,
                               output_dir,
                               show_progress = TRUE) {
  catalog <- candidate_catalog()
  catalog <- catalog[catalog$family %in% families, , drop = FALSE]
  if (!is.null(candidates)) {
    candidates <- unique(as.character(candidates))
    missing <- setdiff(candidates, catalog$candidate)
    if (length(missing)) stop(sprintf("Unknown candidate(s): %s", paste(missing, collapse = ", ")))
    catalog <- catalog[catalog$candidate %in% candidates, , drop = FALSE]
  }
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) stop("`cores` must be a positive integer.")
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Parallel candidate scans require a Unix platform.")
  }
  dimensions <- sort(unique(as.integer(dimensions)))
  jobs <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    candidate <- catalog[i, , drop = FALSE]
    do.call(rbind, lapply(dimensions, function(d) {
      cbind(candidate[rep(1L, batches), , drop = FALSE], d = d, batch = seq_len(batches))
    }))
  }))
  rownames(jobs) <- NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "batch_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  manifest_path <- file.path(output_dir, "manifest.txt")
  if (!file.exists(manifest_path)) writeLines(c(
    "Dimension-scaled Section 6 candidate scan; known-parameter population signal only.",
    "S is the centre-averaged profile KS discrepancy; J is split-sample unbiased for the CvM integrand.",
    sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
    sprintf("families: %s", paste(families, collapse = ",")),
    sprintf("n_centres: %d; n_profile: %d; n_eval: %d; batches: %d", n_centres, n_profile, n_eval, batches),
    sprintf("seed: %d", seed)
  ), manifest_path)
  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else data.frame()
  key <- function(x) paste(x$candidate, x$d, x$batch, sep = "|")
  pending <- jobs[!key(jobs) %in% if (nrow(existing)) key(existing) else character(), , drop = FALSE]
  if (isTRUE(show_progress)) message(sprintf("Candidate signal scan: %d pending batches.", nrow(pending)))
  completed_now <- 0L
  for (first in seq.int(1L, nrow(pending), by = cores)) {
    indices <- seq.int(first, min(first + cores - 1L, nrow(pending)))
    evaluate <- function(index) {
      job <- pending[index, , drop = FALSE]
      job_seed <- section6_seed(
        seed, match(as.character(job$candidate), catalog$candidate),
        job$batch + 100L * job$d, stream = 92L
      )
      metrics <- candidate_one_batch(job, job$d, n_centres, n_profile, n_eval, job_seed)
      cbind(job, seed = job_seed, metrics)
    }
    rows <- if (length(indices) == 1L) list(evaluate(indices)) else {
      parallel::mclapply(indices, evaluate,
        mc.cores = min(cores, length(indices)), mc.preschedule = FALSE
      )
    }
    existing <- rbind(existing, do.call(rbind, rows))
    utils::write.csv(existing, result_path, row.names = FALSE)
    utils::write.csv(candidate_summarize(existing), summary_path, row.names = FALSE)
    completed_now <- completed_now + length(indices)
    if (isTRUE(show_progress)) {
      message(sprintf("completed batches: %d/%d", nrow(jobs) - nrow(pending) + completed_now, nrow(jobs)))
    }
  }
  invisible(list(batch_results = existing, summary = candidate_summarize(existing)))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  families <- parse_section6_csv(args$families, c("normal", "lg", "vmf", "hvmf"), "character")
  result <- run_candidate_scan(
    dimensions = parse_section6_csv(args$dimensions, c(2L, 5L, 10L), "integer"),
    families = families,
    candidates = parse_section6_csv(args$candidates, NULL, "character"),
    n_centres = as.integer(args$n_centres %||% 96L),
    n_profile = as.integer(args$n_profile %||% 6000L),
    n_eval = as.integer(args$n_eval %||% 3000L),
    batches = as.integer(args$batches %||% 3L),
    cores = as.integer(args$cores %||% 1L),
    seed = as.integer(args$seed %||% 20260729L),
    output_dir = args$output_dir %||% file.path("simulation_results", "section6_new_scenarios", "candidate_signal_scan"),
    show_progress = parse_bool(args$show_progress, TRUE)
  )
  print(result$summary)
}
