#!/usr/bin/env Rscript

# Section 6 power experiments for the revised, dimension-indexed scenarios.
#
# This runner covers the dimension-indexed Normal, logistic Gaussian, vMF,
# and HvMF scenarios of Section 6. Each Monte Carlo replication requests the
# sample-point KS statistic and the fast multiplier bootstrap. The common
# sample-point KS/CvM bootstrap loop is fused and evaluated by the C++ kernel.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source("bootstrap/multiplier_bootstrap.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

section6_families <- c("vmf", "hvmf", "normal", "lg")

parse_section6_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    fields <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[fields[[1L]]]] <- if (length(fields) == 1L) "TRUE" else paste(fields[-1L], collapse = "=")
  }
  out
}

parse_section6_csv <- function(value, default, storage.mode = "numeric") {
  if (is.null(value) || !nzchar(value)) return(default)
  values <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  switch(
    storage.mode,
    integer = as.integer(values),
    numeric = as.numeric(values),
    character = as.character(values),
    stop("Unsupported storage mode.")
  )
}

section6_seed <- function(base_seed, design_id, rep, stream = 0L) {
  modulus <- 2147483647
  value <- (as.numeric(base_seed) + 1000003 * as.numeric(design_id) +
    1009 * as.numeric(rep) + 10000019 * as.numeric(stream)) %% modulus
  as.integer(value) + 1L
}

section6_e <- function(d, index = 1L) {
  out <- numeric(d)
  out[[as.integer(index)]] <- 1
  out
}

section6_hvmf_radial_mu1 <- function(d, c_shift) {
  rho0 <- asinh(1)
  rho1 <- rho0 + c_shift / sqrt(d)
  c(cosh(rho1), sinh(rho1) * section6_e(d))
}

section6_sigma <- function(d, sign = c("plus", "minus")) {
  sign <- match.arg(sign)
  if (d < 2L) stop("The Section 6 covariance scenarios require d >= 2.")
  rho <- if (identical(sign, "plus")) 0.75 else -0.75
  sigma <- diag(d)
  sigma[1L, 2L] <- rho
  sigma[2L, 1L] <- rho
  sigma
}

# Dimension-scaled version of the Normal/LG mixture used only in the
# Section 6 prepilot.  It keeps the Euclidean norm of the location change
# equal to 0.5 and the Frobenius norm of each covariance perturbation equal
# to 0.75, while distributing the covariance change over all transverse
# directions.
section6_transverse_direction <- function(d) {
  if (d < 2L) stop("The transverse Section 6 scenarios require d >= 2.")
  out <- numeric(d)
  out[2:d] <- 1 / sqrt(d - 1L)
  out
}

section6_sigma_transverse <- function(d, sign = c("plus", "minus")) {
  sign <- match.arg(sign)
  shift <- if (identical(sign, "plus")) 0.75 else -0.75
  diagonal_shift <- numeric(d)
  diagonal_shift[2:d] <- shift / sqrt(d - 1L)
  diag(1 + diagonal_shift, d)
}

r_standardized_multivariate_t <- function(n, d, nu = 3) {
  n <- as.integer(n)
  d <- as.integer(d)
  if (n == 0L) return(matrix(numeric(), nrow = 0L, ncol = d))
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  w <- stats::rchisq(n, df = nu)
  z * sqrt((nu - 2) / w)
}

projected_normal_on_sphere <- function(n, d) {
  ambient_dim <- as.integer(d) + 1L
  z <- matrix(stats::rnorm(n * ambient_dim), nrow = n, ncol = ambient_dim)
  z[, 1L] <- z[, 1L] + sqrt(ambient_dim)
  z / sqrt(rowSums(z^2))
}

section6_scenario_catalog <- function() {
  list(
    normal_1_mixture = list(
      family = "normal",
      alternative = "opposite_location_correlation_mixture",
      description = "(1-beta/2) N_d(0.5 e1,Sigma_+) + (beta/2) N_d(-0.5 e1,Sigma_-)",
      generator = "normal_mixture"
    ),
    normal_2_t3 = list(
      family = "normal",
      alternative = "standardized_multivariate_t3",
      description = "(1-beta) N_d(0,I_d) + beta t_{3,d}^{std}",
      generator = "normal_t3"
    ),
    normal_1_mixture_transverse = list(
      family = "normal",
      alternative = "opposite_location_transverse_scale_mixture",
      description = "(1-beta/2) N_d(0.5 v_d,Sigma_+) + (beta/2) N_d(-0.5 v_d,Sigma_-); transverse fixed-norm perturbations",
      generator = "normal_mixture_transverse",
      experimental = TRUE
    ),
    lg_1_mixture = list(
      family = "lg",
      alternative = "ilr_opposite_location_correlation_mixture",
      description = "ilr analogue of normal_1_mixture",
      generator = "lg_mixture"
    ),
    lg_2_t3 = list(
      family = "lg",
      alternative = "inverse_ilr_standardized_multivariate_t3",
      description = "inverse-ilr image of (1-beta) N_d(0,I_d) + beta t_{3,d}^{std}",
      generator = "lg_t3"
    ),
    lg_1_mixture_transverse = list(
      family = "lg",
      alternative = "ilr_opposite_location_transverse_scale_mixture",
      description = "ILR analogue of normal_1_mixture_transverse",
      generator = "lg_mixture_transverse",
      experimental = TRUE
    ),
    vmf_1_antipodal = list(
      family = "vmf",
      alternative = "antipodal_vmf_mixture",
      description = "(1-beta/2) vMF(e1,d) + (beta/2) vMF(-e1,d)",
      generator = "vmf_antipodal"
    ),
    vmf_2_projected_normal = list(
      family = "vmf",
      alternative = "projected_normal_mixture",
      description = "(1-beta) vMF(e1,1.5d) + beta Law(Z/||Z||), Z~N(sqrt(d+1)e1,I)",
      generator = "vmf_projected_normal"
    ),
    hvmf_1_mixture = list(
      family = "hvmf",
      alternative = "halfway_hvmf_location_mixture",
      description = "(1-beta/2) HvMF((sqrt(2),e1),d) + (beta/2) HvMF((sqrt(2),e2),d)",
      generator = "hvmf_location_mixture"
    ),
    hvmf_2_angular = list(
      family = "hvmf",
      alternative = "conditional_angular_hvmf_mixture",
      description = "(1-beta) HvMF((sqrt(2),e1),1.5d) + beta Q_delta, delta=pi/5",
      generator = "hvmf_angular_mixture"
    ),
    hvmf_1_dimension_scaled_location_mixture = list(
      family = "hvmf",
      alternative = "dimension_scaled_hvmf_location_mixture",
      description = "(1-beta/2) HvMF((sqrt(2),e1),d) + (beta/2) HvMF((sqrt(1+2/d),sqrt(2/d)e2),d)",
      generator = "hvmf_dimension_scaled_location_mixture",
      experimental = TRUE
    ),
    hvmf_1_radial_c_inv_sqrt2 = list(
      family = "hvmf",
      alternative = "radial_hvmf_location_mixture_c_inv_sqrt2",
      description = "(1-beta/2) HvMF((sqrt(2),e1),d) + (beta/2) HvMF(mu1,d), d_H(mu0,mu1)=1/sqrt(2d)",
      generator = "hvmf_radial_location_mixture",
      experimental = TRUE
    ),
    hvmf_1_radial_c_half = list(
      family = "hvmf",
      alternative = "radial_hvmf_location_mixture_c_half",
      description = "(1-beta/2) HvMF((sqrt(2),e1),d) + (beta/2) HvMF(mu1,d), d_H(mu0,mu1)=1/(2sqrt(d))",
      generator = "hvmf_radial_location_mixture",
      experimental = TRUE
    ),
    hvmf_2_angular_sqrt_d_concentration = list(
      family = "hvmf",
      alternative = "conditional_angular_hvmf_mixture_sqrt_d_concentration",
      description = "(1-beta) HvMF((sqrt(2),e1),2sqrt(d)) + beta Q_delta, delta=pi/5",
      generator = "hvmf_angular_sqrt_d_concentration",
      experimental = TRUE
    )
  )
}

section6_family_scenarios <- function(family, include_experimental = FALSE) {
  catalog <- section6_scenario_catalog()
  names(catalog)[vapply(catalog, function(x) {
    identical(x$family, family) && (isTRUE(include_experimental) || !isTRUE(x$experimental))
  }, logical(1))]
}

make_section6_design <- function(family,
                                 dimensions = c(2L, 10L),
                                 n_values = c(50L, 100L, 200L, 400L),
                                 beta_values = c(0, 0.5, 1),
                                 scenarios = NULL) {
  if (!family %in% section6_families) stop("Unknown Section 6 family.")
  catalog <- section6_scenario_catalog()
  available <- names(catalog)[vapply(catalog, function(x) identical(x$family, family), logical(1))]
  scenario_names <- if (is.null(scenarios)) {
    section6_family_scenarios(family)
  } else {
    scenarios <- unique(as.character(scenarios))
    invalid <- setdiff(scenarios, available)
    if (length(invalid)) {
      stop(sprintf("Invalid scenario(s) for family '%s': %s", family, paste(invalid, collapse = ", ")))
    }
    scenarios
  }
  rows <- list()
  index <- 1L
  for (scenario in scenario_names) {
    for (d in sort(unique(as.integer(dimensions)))) {
      for (n in sort(unique(as.integer(n_values)))) {
        for (beta in sort(unique(as.numeric(beta_values)))) {
          rows[[index]] <- data.frame(
            scenario = scenario,
            family = family,
            alternative = catalog[[scenario]]$alternative,
            description = catalog[[scenario]]$description,
            d = d,
            n = n,
            beta = beta,
            stringsAsFactors = FALSE
          )
          index <- index + 1L
        }
      }
    }
  }
  design <- do.call(rbind, rows)
  design$design_id <- seq_len(nrow(design))
  design
}

generate_section6_sample <- function(design_row) {
  scenario <- as.character(design_row$scenario)
  d <- as.integer(design_row$d)
  n <- as.integer(design_row$n)
  beta <- as.numeric(design_row$beta)
  e1 <- section6_e(d)

  if (identical(scenario, "normal_1_mixture")) {
    choose_alt <- stats::runif(n) < beta / 2
    x <- mvtnorm::rmvnorm(n, mean = 0.5 * e1, sigma = section6_sigma(d, "plus"))
    if (any(choose_alt)) {
      x[choose_alt, ] <- mvtnorm::rmvnorm(
        sum(choose_alt), mean = -0.5 * e1, sigma = section6_sigma(d, "minus")
      )
    }
    return(x)
  }

  if (identical(scenario, "normal_2_t3")) {
    choose_alt <- stats::runif(n) < beta
    x <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    if (any(choose_alt)) x[choose_alt, ] <- r_standardized_multivariate_t(sum(choose_alt), d)
    return(x)
  }

  if (identical(scenario, "normal_1_mixture_transverse")) {
    direction <- section6_transverse_direction(d)
    choose_alt <- stats::runif(n) < beta / 2
    x <- mvtnorm::rmvnorm(n, mean = 0.5 * direction, sigma = section6_sigma_transverse(d, "plus"))
    if (any(choose_alt)) {
      x[choose_alt, ] <- mvtnorm::rmvnorm(
        sum(choose_alt), mean = -0.5 * direction, sigma = section6_sigma_transverse(d, "minus")
      )
    }
    return(x)
  }

  if (identical(scenario, "lg_1_mixture")) {
    choose_alt <- stats::runif(n) < beta / 2
    z <- mvtnorm::rmvnorm(n, mean = 0.5 * e1, sigma = section6_sigma(d, "plus"))
    if (any(choose_alt)) {
      z[choose_alt, ] <- mvtnorm::rmvnorm(
        sum(choose_alt), mean = -0.5 * e1, sigma = section6_sigma(d, "minus")
      )
    }
    return(logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L))
  }

  if (identical(scenario, "lg_2_t3")) {
    choose_alt <- stats::runif(n) < beta
    z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    if (any(choose_alt)) z[choose_alt, ] <- r_standardized_multivariate_t(sum(choose_alt), d)
    return(logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L))
  }

  if (identical(scenario, "lg_1_mixture_transverse")) {
    direction <- section6_transverse_direction(d)
    choose_alt <- stats::runif(n) < beta / 2
    z <- mvtnorm::rmvnorm(n, mean = 0.5 * direction, sigma = section6_sigma_transverse(d, "plus"))
    if (any(choose_alt)) {
      z[choose_alt, ] <- mvtnorm::rmvnorm(
        sum(choose_alt), mean = -0.5 * direction, sigma = section6_sigma_transverse(d, "minus")
      )
    }
    return(logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L))
  }

  if (identical(scenario, "vmf_1_antipodal")) {
    mu <- section6_e(d + 1L)
    choose_alt <- stats::runif(n) < beta / 2
    x <- rotasym::r_vMF(n, mu = mu, kappa = d)
    if (any(choose_alt)) x[choose_alt, ] <- rotasym::r_vMF(sum(choose_alt), mu = -mu, kappa = d)
    return(x)
  }

  if (identical(scenario, "vmf_2_projected_normal")) {
    mu <- section6_e(d + 1L)
    choose_alt <- stats::runif(n) < beta
    x <- rotasym::r_vMF(n, mu = mu, kappa = 1.5 * d)
    if (any(choose_alt)) x[choose_alt, ] <- projected_normal_on_sphere(sum(choose_alt), d)
    return(x)
  }

  if (identical(scenario, "hvmf_1_mixture")) {
    mu0 <- c(sqrt(2), section6_e(d))
    mu1 <- c(sqrt(2), section6_e(d, index = 2L))
    choose_alt <- stats::runif(n) < beta / 2
    x <- rhvmf_polar(n, mu = mu0, kappa = d)
    if (any(choose_alt)) {
      x[choose_alt, ] <- rhvmf_polar(sum(choose_alt), mu = mu1, kappa = d)
    }
    return(x)
  }

  if (identical(scenario, "hvmf_2_angular")) {
    mu0 <- c(sqrt(2), section6_e(d))
    tangent <- section6_e(d, index = 2L)
    choose_alt <- stats::runif(n) < beta
    x <- rhvmf_polar(n, mu = mu0, kappa = 1.5 * d)
    if (any(choose_alt)) {
      x[choose_alt, ] <- rhvmf_angular_mixture(
        sum(choose_alt), mu = mu0, kappa = 1.5 * d,
        delta = pi / 5, tangent = tangent
      )
    }
    return(x)
  }

  if (identical(scenario, "hvmf_1_dimension_scaled_location_mixture")) {
    mu0 <- c(sqrt(2), section6_e(d))
    spatial_scale <- sqrt(2 / d)
    mu1 <- c(
      sqrt(1 + spatial_scale^2),
      spatial_scale * section6_e(d, index = 2L)
    )
    choose_alt <- stats::runif(n) < beta / 2
    x <- rhvmf_polar(n, mu = mu0, kappa = d)
    if (any(choose_alt)) {
      x[choose_alt, ] <- rhvmf_polar(sum(choose_alt), mu = mu1, kappa = d)
    }
    return(x)
  }

  if (scenario %in% c(
      "hvmf_1_radial_c_inv_sqrt2",
      "hvmf_1_radial_c_half"
  )) {
    mu0 <- c(sqrt(2), section6_e(d))
    c_shift <- switch(
      scenario,
      hvmf_1_radial_c_inv_sqrt2 = 1 / sqrt(2),
      hvmf_1_radial_c_half = 1 / 2
    )
    mu1 <- section6_hvmf_radial_mu1(d, c_shift)
    choose_alt <- stats::runif(n) < beta / 2
    x <- rhvmf_polar(n, mu = mu0, kappa = d)
    if (any(choose_alt)) {
      x[choose_alt, ] <- rhvmf_polar(sum(choose_alt), mu = mu1, kappa = d)
    }
    return(x)
  }

  if (identical(scenario, "hvmf_2_angular_sqrt_d_concentration")) {
    mu0 <- c(sqrt(2), section6_e(d))
    tangent <- section6_e(d, index = 2L)
    kappa <- 2 * sqrt(d)
    choose_alt <- stats::runif(n) < beta
    x <- rhvmf_polar(n, mu = mu0, kappa = kappa)
    if (any(choose_alt)) {
      x[choose_alt, ] <- rhvmf_angular_mixture(
        sum(choose_alt), mu = mu0, kappa = kappa,
        delta = pi / 5, tangent = tangent
      )
    }
    return(x)
  }

  stop(sprintf("Unsupported scenario '%s'.", scenario))
}

section6_control <- function(derivative_mc_size, derivative_seed, cvm_block_size,
                             derivative_method = "score_mc") {
  derivative_method <- tolower(as.character(derivative_method))
  if (length(derivative_method) != 1L ||
      !derivative_method %in% c("score_mc", "quadrature")) {
    stop("`derivative_method` must be either 'score_mc' or 'quadrature'.")
  }
  list(
    derivative_method = derivative_method,
    derivative_mc_size = as.integer(derivative_mc_size),
    derivative_mc_seed = as.integer(derivative_seed),
    fast_multiplier_cvm_block_size = as.integer(cvm_block_size),
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fast_multiplier_fuse_ks_cvm = TRUE,
    fast_multiplier_cache_corrections = "auto",
    fast_multiplier_stream_chunk_size = 100L,
    vmf_profile_method = "tabulated",
    vmf_profile_n_u = 4097L,
    hvmf_profile_method = "tabulated",
    hvmf_profile_n_y = 4097L,
    # EN DUDA (2026-07-26): generic MVN/LG dispatcher control.
    mvnormal_quadform_method = "auto"
  )
}

run_section6_bootstrap <- function(design_row, x, B, seed, derivative_seed,
                                   derivative_mc_size, cvm_block_size,
                                   derivative_method = "score_mc",
                                   bootstrap_method = "fast_multiplier",
                                   n_cores = 1L) {
  family <- as.character(design_row$family)
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    bootstrap_method = bootstrap_method,
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
    control = section6_control(
      derivative_mc_size, derivative_seed, cvm_block_size,
      derivative_method = derivative_method
    ),
    distance_profile_backend = "r"
  )
  if (identical(family, "normal")) {
    return(do.call(multiplier_bootstrap_mvnormal, c(common, list(
      unknown_param = "both", fast_multiplier_backend = "cpp", fuse_ks_cvm = TRUE,
      cache_block_corrections = "auto"
    ))))
  }
  if (identical(family, "lg")) {
    return(do.call(multiplier_bootstrap_logistic_gaussian, c(common, list(unknown_param = "both"))))
  }
  if (identical(family, "vmf")) {
    return(do.call(multiplier_bootstrap_vmf, c(common, list(
      unknown_param = "xi", distance_type = "geodesic"
    ))))
  }
  if (identical(family, "hvmf")) {
    return(do.call(multiplier_bootstrap_hvmf, c(common, list(
      unknown_param = "both", fast_multiplier_backend = "cpp", fuse_ks_cvm = TRUE,
      cache_block_corrections = "auto"
    ))))
  }
  stop(sprintf("Unsupported family '%s'.", family))
}

empty_section6_results <- function() {
  data.frame(
    scenario = character(), family = character(), alternative = character(), d = integer(), n = integer(), beta = numeric(),
    design_id = integer(), rep = integer(), status = character(), error_message = character(),
    ks_pvalue = numeric(), cvm_pvalue = numeric(), ks_reject = logical(), cvm_reject = logical(),
    bootstrap_method_requested = character(), bootstrap_method_effective = character(), fallback_to_reestimated = logical(),
    derivative_method_requested = character(), derivative_method_effective = character(),
    derivative_method_selection_source = character(),
    quadrature_algorithm = character(), quadrature_abs_tol = numeric(),
    quadrature_max_terms = integer(), quadrature_max_terms_used = integer(),
    quadrature_max_residual_error_estimate = numeric(),
    ks_grid = character(), fast_multiplier_backend_requested = character(), fast_multiplier_backend_effective = character(),
    fast_multiplier_cpp_kernel_requested = character(), fast_multiplier_cpp_kernel_effective = character(),
    fast_multiplier_fuse_ks_cvm_requested = logical(), fast_multiplier_fuse_ks_cvm_effective = logical(),
    seed_data = integer(), seed_bootstrap = integer(), seed_derivative = integer(), elapsed_seconds = numeric(),
    stringsAsFactors = FALSE
  )
}

run_section6_job <- function(job, B, base_seed, derivative_mc_size, cvm_block_size,
                             derivative_method = "score_mc") {
  started <- proc.time()[["elapsed"]]
  data_seed <- section6_seed(base_seed, job$design_id, job$rep, 0L)
  bootstrap_seed <- section6_seed(base_seed, job$design_id, job$rep, 1L)
  derivative_seed <- section6_seed(base_seed, job$design_id, job$rep, 2L)
  out <- data.frame(
    scenario = as.character(job$scenario), family = as.character(job$family), alternative = as.character(job$alternative),
    d = as.integer(job$d), n = as.integer(job$n), beta = as.numeric(job$beta), design_id = as.integer(job$design_id),
    rep = as.integer(job$rep), status = "ok", error_message = NA_character_,
    ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
    bootstrap_method_requested = "fast_multiplier", bootstrap_method_effective = NA_character_, fallback_to_reestimated = NA,
    derivative_method_requested = derivative_method, derivative_method_effective = NA_character_,
    derivative_method_selection_source = "explicit",
    quadrature_algorithm = NA_character_, quadrature_abs_tol = NA_real_,
    quadrature_max_terms = NA_integer_, quadrature_max_terms_used = NA_integer_,
    quadrature_max_residual_error_estimate = NA_real_,
    ks_grid = "sample_unique_distances", fast_multiplier_backend_requested = "cpp", fast_multiplier_backend_effective = NA_character_,
    fast_multiplier_cpp_kernel_requested = "contiguous_double", fast_multiplier_cpp_kernel_effective = NA_character_,
    fast_multiplier_fuse_ks_cvm_requested = TRUE, fast_multiplier_fuse_ks_cvm_effective = NA,
    seed_data = data_seed, seed_bootstrap = bootstrap_seed, seed_derivative = derivative_seed,
    elapsed_seconds = NA_real_, stringsAsFactors = FALSE
  )
  out <- tryCatch({
    set.seed(data_seed)
    x <- generate_section6_sample(job)
    fit <- run_section6_bootstrap(
      job, x, B, bootstrap_seed, derivative_seed, derivative_mc_size,
      cvm_block_size, derivative_method = derivative_method
    )
    diagnostics <- fit$diagnostics
    out$ks_pvalue <- fit$inference$ks$p_value
    out$cvm_pvalue <- fit$inference$cvm$p_value
    out$ks_reject <- fit$inference$ks$reject
    out$cvm_reject <- fit$inference$cvm$reject
    out$bootstrap_method_effective <- diagnostics$effective_bootstrap_method %||% NA_character_
    out$derivative_method_effective <- diagnostics$derivative_method_effective %||%
      diagnostics$derivative_method %||% NA_character_
    out$derivative_method_selection_source <-
      diagnostics$derivative_method_selection_source %||% NA_character_
    out$quadrature_algorithm <- diagnostics$quadrature_algorithm %||% NA_character_
    out$quadrature_abs_tol <- diagnostics$quadrature_abs_tol %||% NA_real_
    out$quadrature_max_terms <- diagnostics$quadrature_max_terms %||% NA_integer_
    out$quadrature_max_terms_used <-
      diagnostics$quadrature_max_terms_used %||% NA_integer_
    out$quadrature_max_residual_error_estimate <-
      diagnostics$quadrature_max_residual_error_estimate %||% NA_real_
    out$fallback_to_reestimated <- isTRUE(diagnostics$fallback_to_reestimated)
    out$fast_multiplier_backend_effective <- diagnostics$fast_multiplier_backend_effective %||% NA_character_
    out$fast_multiplier_cpp_kernel_effective <-
      diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_
    out$fast_multiplier_fuse_ks_cvm_effective <- isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective)
    if (!identical(out$bootstrap_method_effective, "fast_multiplier") ||
        !identical(out$fast_multiplier_backend_effective, "cpp") ||
        !identical(out$fast_multiplier_cpp_kernel_effective, "contiguous_double") ||
        !isTRUE(out$fast_multiplier_fuse_ks_cvm_effective) ||
        !identical(out$derivative_method_effective, derivative_method)) {
      out$status <- "nonconforming"
      out$error_message <- paste(
        "Requested fast, fused C++ KS/CvM bootstrap or derivative method was not effective.",
        sprintf("requested_derivative_method=%s", derivative_method),
        sprintf("effective_derivative_method=%s", out$derivative_method_effective),
        sprintf("effective_cpp_kernel=%s", out$fast_multiplier_cpp_kernel_effective)
      )
    }
    out
  }, error = function(e) {
    out$status <- "error"
    out$error_message <- conditionMessage(e)
    out
  })
  out$elapsed_seconds <- proc.time()[["elapsed"]] - started
  out
}

section6_write_atomic_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop(sprintf("Could not update '%s'.", path))
}

section6_output_lock_path <- function(output_dir) {
  file.path(output_dir, ".section6.lock")
}

section6_acquire_output_lock <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock_path <- section6_output_lock_path(output_dir)
  acquired <- dir.create(lock_path, showWarnings = FALSE)
  if (!isTRUE(acquired)) {
    owner_path <- file.path(lock_path, "owner.txt")
    owner <- if (file.exists(owner_path)) {
      paste(readLines(owner_path, warn = FALSE), collapse = "; ")
    } else {
      "owner metadata unavailable"
    }
    stop(sprintf(
      "Section 6 output directory '%s' is already locked by another process or job (%s).",
      output_dir, owner
    ), call. = FALSE)
  }

  initialized <- FALSE
  on.exit({
    if (!initialized && dir.exists(lock_path)) {
      unlink(lock_path, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  hostname <- unname(Sys.info()[["nodename"]])
  if (is.null(hostname) || is.na(hostname) || !nzchar(hostname)) {
    hostname <- "unknown"
  }
  slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
  if (!nzchar(slurm_job_id)) slurm_job_id <- "unavailable"
  acquired_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
  owner_token <- paste(hostname, Sys.getpid(), acquired_at, sep = "|")
  owner <- c(
    sprintf("slurm_job_id: %s", slurm_job_id),
    sprintf("hostname: %s", hostname),
    sprintf("pid: %d", Sys.getpid()),
    sprintf("acquired_at_utc: %s", acquired_at),
    sprintf("owner_token: %s", owner_token)
  )
  writeLines(owner, file.path(lock_path, "owner.txt"))
  initialized <- TRUE

  structure(
    list(path = lock_path, owner_token = owner_token),
    class = "section6_output_lock"
  )
}

section6_release_output_lock <- function(lock) {
  if (is.null(lock) || is.null(lock$path) || !dir.exists(lock$path)) {
    return(invisible(FALSE))
  }
  owner_path <- file.path(lock$path, "owner.txt")
  owner <- if (file.exists(owner_path)) readLines(owner_path, warn = FALSE) else character()
  expected <- sprintf("owner_token: %s", lock$owner_token)
  if (!expected %in% owner) {
    warning(sprintf(
      "Refusing to release Section 6 lock '%s' because its owner token changed.",
      lock$path
    ), call. = FALSE)
    return(invisible(FALSE))
  }
  unlink(lock$path, recursive = TRUE, force = TRUE)
  invisible(!dir.exists(lock$path))
}

section6_design_key <- function(x, include_rep = FALSE) {
  required <- c("scenario", "d", "n", "beta")
  if (!all(required %in% names(x))) {
    stop("Section 6 design rows must contain scenario, d, n, and beta.")
  }
  pieces <- list(
    as.character(x$scenario),
    as.integer(x$d),
    as.integer(x$n),
    formatC(as.numeric(x$beta), digits = 16L, format = "fg", flag = "#")
  )
  if (isTRUE(include_rep)) {
    if (!"rep" %in% names(x)) stop("Replication rows must contain `rep`.")
    pieces <- c(pieces, list(as.integer(x$rep)))
  }
  do.call(paste, c(pieces, sep = "|"))
}

section6_validate_manifest_design <- function(manifest_path, design, M, B,
                                              base_seed,
                                              derivative_mc_size,
                                              cvm_block_size,
                                              derivative_method = "score_mc") {
  if (!file.exists(manifest_path)) return(invisible(TRUE))
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  required <- c(
    "scenario", "d", "n", "beta", "M", "B", "base_seed",
    "derivative_mc_size", "cvm_block_size", "ks_grid",
    "fast_multiplier_backend", "fast_multiplier_cpp_kernel", "fused_ks_cvm_kernel"
  )
  if (!all(required %in% names(manifest))) {
    stop("Existing Section 6 manifest is incomplete; use a new output directory.")
  }

  checks <- c(
    design = setequal(section6_design_key(manifest), section6_design_key(design)),
    M = all(as.integer(manifest$M) == as.integer(M)),
    B = all(as.integer(manifest$B) == as.integer(B)),
    base_seed = all(as.integer(manifest$base_seed) == as.integer(base_seed)),
    derivative_mc_size = all(
      as.integer(manifest$derivative_mc_size) == as.integer(derivative_mc_size)
    ),
    derivative_method = all(
      as.character(
        if ("derivative_method" %in% names(manifest)) {
          manifest$derivative_method
        } else {
          rep.int("score_mc", nrow(manifest))
        }
      ) == as.character(derivative_method)
    ),
    cvm_block_size = all(
      as.integer(manifest$cvm_block_size) == as.integer(cvm_block_size)
    ),
    ks_grid = all(as.character(manifest$ks_grid) == "sample_unique_distances"),
    fast_multiplier_backend = all(
      as.character(manifest$fast_multiplier_backend) == "cpp"
    ),
    fast_multiplier_cpp_kernel = all(
      as.character(manifest$fast_multiplier_cpp_kernel) == "contiguous_double"
    ),
    fused_ks_cvm_kernel = all(as.logical(manifest$fused_ks_cvm_kernel))
  )
  checks[is.na(checks)] <- FALSE
  if (!all(checks)) {
    stop(sprintf(
      paste(
        "The existing output directory has an incompatible Section 6 manifest",
        "(mismatched fields: %s). Do not resume into it; use a new output",
        "directory or recover matching rows first."
      ),
      paste(names(checks)[!checks], collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

section6_make_manifest <- function(design, M, B, cores, base_seed,
                                   derivative_mc_size, cvm_block_size,
                                   derivative_method = "score_mc") {
  transform(design,
    M = as.integer(M), B = as.integer(B), cores = as.integer(cores),
    base_seed = as.integer(base_seed), derivative_mc_size = as.integer(derivative_mc_size),
    cvm_block_size = as.integer(cvm_block_size),
    derivative_method = as.character(derivative_method),
    ks_grid = "sample_unique_distances", bootstrap_method = "fast_multiplier",
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fused_ks_cvm_kernel = TRUE
  )
}

# Copy only completed, verifiably conforming rows into a fresh directory whose
# manifest matches the requested design.  This is deliberately non-destructive:
# the source directory is never modified, and ambiguous duplicate replications
# cause an error instead of being silently chosen.
recover_section6_results <- function(source_dir, target_dir, family,
                                     M = 1000L, B = 5000L,
                                     dimensions = c(2L, 10L),
                                     n_values = c(50L, 100L, 200L, 400L),
                                     beta_values = c(0, 0.5, 1),
                                     cores = 8L,
                                     base_seed = 20260727L,
                                     derivative_mc_size = 10000L,
                                     cvm_block_size = 50L,
                                     derivative_method = "score_mc") {
  source_path <- file.path(source_dir, "raw_results.csv")
  target_result_path <- file.path(target_dir, "raw_results.csv")
  target_manifest_path <- file.path(target_dir, "manifest.csv")
  if (!file.exists(source_path)) stop("The source directory has no raw_results.csv.", call. = FALSE)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  output_lock <- section6_acquire_output_lock(target_dir)
  on.exit(section6_release_output_lock(output_lock), add = TRUE)
  if (file.exists(target_result_path) || file.exists(target_manifest_path)) {
    stop("Recovery target already contains results or a manifest; choose a new directory.", call. = FALSE)
  }

  design <- make_section6_design(family, dimensions, n_values, beta_values)
  jobs <- merge(design, data.frame(rep = seq_len(as.integer(M))), by = NULL)
  target_keys <- section6_design_key(jobs, include_rep = TRUE)
  source <- utils::read.csv(source_path, stringsAsFactors = FALSE)
  required <- names(empty_section6_results())
  missing <- setdiff(required, names(source))
  if (length(missing)) {
    stop(sprintf("Source results lack required columns: %s.", paste(missing, collapse = ", ")), call. = FALSE)
  }
  source_key <- section6_design_key(source, include_rep = TRUE)
  keep <- source$status == "ok" & source_key %in% target_keys
  recovered <- source[keep, required, drop = FALSE]
  recovered_key <- section6_design_key(recovered, include_rep = TRUE)
  if (anyDuplicated(recovered_key)) {
    stop("Recovery found duplicate completed replications; refusing to select among them.", call. = FALSE)
  }
  conforming <- recovered$bootstrap_method_effective == "fast_multiplier" &
    recovered$fast_multiplier_backend_effective == "cpp" &
    recovered$fast_multiplier_cpp_kernel_effective == "contiguous_double" &
    recovered$fast_multiplier_fuse_ks_cvm_effective &
    recovered$ks_grid == "sample_unique_distances" &
    !recovered$fallback_to_reestimated
  if (!all(conforming)) {
    stop("Recovery found completed target rows that do not use the required fused C++ fast bootstrap.", call. = FALSE)
  }
  if (!all(recovered$family == family)) {
    stop("Recovery source contains a family inconsistent with the requested target design.", call. = FALSE)
  }

  section6_write_atomic_csv(section6_make_manifest(
    design, M, B, cores, base_seed, derivative_mc_size, cvm_block_size,
    derivative_method = derivative_method
  ), target_manifest_path)
  section6_write_atomic_csv(recovered, target_result_path)
  section6_write_atomic_csv(summarize_section6_results(recovered), file.path(target_dir, "summary.csv"))
  started <- Sys.time()
  section6_write_status(
    file.path(target_dir, "progress_status.txt"), family, nrow(jobs), nrow(recovered),
    recovered, started, as.integer(cores)
  )
  writeLines(c(
    sprintf("recovered_at: %s", format(started, tz = "Europe/Madrid")),
    sprintf("source_dir: %s", normalizePath(source_dir)),
    sprintf("target_design_rows: %d", nrow(design)),
    sprintf("target_replications: %d", nrow(jobs)),
    sprintf("recovered_conforming_replications: %d", nrow(recovered)),
    "selection: status=ok; matching scenario,d,n,beta,rep; fused C++ fast bootstrap; sample KS",
    "source_directory_was_not_modified: TRUE"
  ), file.path(target_dir, "recovery_report.txt"))
  cat(sprintf("Recovered %d conforming replications into %s\n", nrow(recovered), target_dir))
  invisible(list(design = design, results = recovered))
}

summarize_section6_results <- function(results) {
  ok <- results[results$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) return(data.frame())
  keys <- interaction(ok$scenario, ok$d, ok$n, ok$beta, drop = TRUE)
  out <- lapply(split(ok, keys), function(x) data.frame(
    scenario = x$scenario[[1L]], family = x$family[[1L]], alternative = x$alternative[[1L]],
    d = x$d[[1L]], n = x$n[[1L]], beta = x$beta[[1L]], M = nrow(x),
    rejection_ks = mean(x$ks_reject), rejection_cvm = mean(x$cvm_reject), stringsAsFactors = FALSE
  ))
  do.call(rbind, out)
}

section6_write_status <- function(path, family, total, completed, results, started, cores) {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  successful <- results$elapsed_seconds[results$status == "ok"]
  per_job <- if (length(successful)) mean(successful) else NA_real_
  eta <- if (is.finite(per_job)) per_job * (total - completed) / cores else NA_real_
  writeLines(c(
    sprintf("family: %s", family),
    sprintf("updated: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("completed: %d/%d", completed, total),
    sprintf("remaining: %d", total - completed),
    sprintf("elapsed_seconds: %.1f", elapsed),
    sprintf("mean_seconds_per_successful_job: %s", if (is.finite(per_job)) format(round(per_job, 3), nsmall = 3) else "NA"),
    sprintf("eta_seconds: %s", if (is.finite(eta)) format(round(eta), scientific = FALSE) else "NA"),
    sprintf("ok: %d", sum(results$status == "ok")),
    sprintf("nonconforming: %d", sum(results$status == "nonconforming")),
    sprintf("errors: %d", sum(results$status == "error"))
  ), path)
}

section6_progress <- function(completed, total, started, results, cores, width = 34L) {
  proportion <- completed / total
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  successful <- results$elapsed_seconds[results$status == "ok"]
  per_job <- if (length(successful)) mean(successful) else NA_real_
  eta <- if (is.finite(per_job)) per_job * (total - completed) / cores else NA_real_
  bar <- paste0(strrep("=", floor(width * proportion)), strrep("-", width - floor(width * proportion)))
  cat(sprintf("\r[%s] %6.2f%%  %d/%d  elapsed %s  ETA %s", bar, 100 * proportion, completed, total,
    format(round(elapsed), units = "secs"), if (is.finite(eta)) format(round(eta), units = "secs") else "--"))
  if (completed == total) cat("\n")
  flush.console()
}

# A forked worker queue provides one shared job queue.  Whenever a worker
# returns a replication, a fresh worker is immediately given the next one.
# Hence heterogeneous n values cannot strand workers at the end of a static
# batch. Results are checkpointed in small groups. `mcparallel` is used rather
# than a socket cluster so this also works on macOS installations that disallow
# opening a local server socket.
section6_run_dynamic_queue <- function(pending, B, base_seed,
                                       derivative_mc_size, cvm_block_size,
                                       derivative_method,
                                       cores, checkpoint_results,
                                       on_checkpoint, on_progress) {
  if (nrow(pending) == 0L) return(invisible(NULL))
  workers <- min(as.integer(cores), nrow(pending))
  next_job <- 1L
  finished <- 0L
  buffer <- empty_section6_results()
  active <- list()
  on.exit({
    if (length(active)) {
      try(parallel:::mckill(unname(active), signal = 2L), silent = TRUE)
    }
  }, add = TRUE)
  submit <- function(index) {
    process <- parallel::mcparallel(
      run_section6_job(pending[index, , drop = FALSE], B, base_seed,
                       derivative_mc_size, cvm_block_size, derivative_method),
      mc.set.seed = FALSE, silent = TRUE
    )
    active[[as.character(process$pid)]] <<- process
  }
  for (worker in seq_len(workers)) {
    submit(next_job)
    next_job <- next_job + 1L
  }
  while (finished < nrow(pending)) {
    received <- parallel::mccollect(active, wait = FALSE, timeout = -1)
    if (is.null(received)) next
    for (pid in names(received)) {
      active[[pid]] <- NULL
      result <- received[[pid]]
      if (!is.data.frame(result) || nrow(result) != 1L) {
        stop("A Section 6 worker returned an invalid replication result.", call. = FALSE)
      }
      buffer <- rbind(buffer, result)
      finished <- finished + 1L
      if (next_job <= nrow(pending)) {
        submit(next_job)
        next_job <- next_job + 1L
      }
      if (nrow(buffer) >= checkpoint_results || finished == nrow(pending)) {
        on_checkpoint(buffer, finished)
        buffer <- empty_section6_results()
      }
      on_progress(finished)
    }
  }
  invisible(NULL)
}

run_section6_family <- function(family,
                                output_dir,
                                M = 1000L,
                                B = 5000L,
                                dimensions = c(2L, 10L),
                                n_values = c(50L, 100L, 200L, 400L),
                                beta_values = c(0, 0.5, 1),
                                cores = 8L,
                                base_seed = 20260727L,
                                derivative_mc_size = 10000L,
                                cvm_block_size = 50L,
                                derivative_method = "score_mc",
                                checkpoint_results = 64L,
                                show_progress = TRUE,
                                scenarios = NULL) {
  family <- tolower(family)
  if (!family %in% section6_families) stop("`family` must be one of vmf, hvmf, normal, or lg.")
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) stop("`cores` must be a positive integer.")
  if (.Platform$OS.type != "unix" && cores > 1L) stop("This runner requires a Unix platform for outer parallelism.")
  derivative_method <- tolower(as.character(derivative_method))
  if (identical(derivative_method, "auto")) {
    derivative_method <- "score_mc"
  }
  if (length(derivative_method) != 1L ||
      !derivative_method %in% c("score_mc", "quadrature")) {
    stop("`derivative_method` must be either 'score_mc' or 'quadrature'.")
  }
  checkpoint_results <- as.integer(checkpoint_results)
  if (!is.finite(checkpoint_results) || checkpoint_results < 1L) {
    stop("`checkpoint_results` must be a positive integer.")
  }
  ensure_distance_profile_cpp_loaded()

  design <- make_section6_design(family, dimensions, n_values, beta_values, scenarios = scenarios)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(output_lock), add = TRUE)
  result_path <- file.path(output_dir, "raw_results.csv")
  manifest_path <- file.path(output_dir, "manifest.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  summary_path <- file.path(output_dir, "summary.csv")
  log_path <- file.path(output_dir, "run.log")
  section6_validate_manifest_design(
    manifest_path, design, M = M, B = B,
    base_seed = base_seed,
    derivative_mc_size = derivative_mc_size,
    cvm_block_size = cvm_block_size,
    derivative_method = derivative_method
  )
  if (!file.exists(manifest_path)) {
    section6_write_atomic_csv(section6_make_manifest(
      design, M, B, cores, base_seed, derivative_mc_size, cvm_block_size,
      derivative_method = derivative_method
    ), manifest_path)
  }
  existing <- if (file.exists(result_path)) {
    prior <- utils::read.csv(result_path, stringsAsFactors = FALSE)
    if (!"derivative_method_requested" %in% names(prior)) {
      prior$derivative_method_requested <- derivative_method
    }
    if (!"derivative_method_effective" %in% names(prior)) {
      prior$derivative_method_effective <- derivative_method
    }
    prior
  } else {
    empty_section6_results()
  }
  jobs <- merge(design, data.frame(rep = seq_len(as.integer(M))), by = NULL)
  key <- function(x) paste(x$design_id, x$rep, sep = "|")
  done <- if (nrow(existing)) {
    section6_design_key(existing[existing$status == "ok", , drop = FALSE], include_rep = TRUE)
  } else {
    character()
  }
  pending <- jobs[!section6_design_key(jobs, include_rep = TRUE) %in% done, , drop = FALSE]
  started <- Sys.time()
  completed <- nrow(jobs) - nrow(pending)
  cat(sprintf("%s family=%s M=%d B=%d cores=%d derivative_method=%s pending=%d\n", format(started, tz = "Europe/Madrid"), family, M, B, cores, derivative_method, nrow(pending)), file = log_path, append = TRUE)
  section6_write_status(status_path, family, nrow(jobs), completed, existing, started, cores)
  if (isTRUE(show_progress)) section6_progress(completed, nrow(jobs), started, existing, cores)
  if (nrow(pending) == 0L) return(invisible(list(results = existing, summary = summarize_section6_results(existing))))

  section6_run_dynamic_queue(
    pending = pending, B = B, base_seed = base_seed,
    derivative_mc_size = derivative_mc_size, cvm_block_size = cvm_block_size,
    derivative_method = derivative_method,
    cores = cores, checkpoint_results = checkpoint_results,
    on_checkpoint = function(rows, finished) {
    existing <<- rbind(existing, rows)
    existing <- existing[order(existing$design_id, existing$rep), , drop = FALSE]
    section6_write_atomic_csv(existing, result_path)
    section6_write_atomic_csv(summarize_section6_results(existing), summary_path)
    completed <<- nrow(jobs) - nrow(pending) + finished
    section6_write_status(status_path, family, nrow(jobs), completed, existing, started, cores)
    cat(sprintf("%s completed=%d/%d\n", format(Sys.time(), tz = "Europe/Madrid"), completed, nrow(jobs)), file = log_path, append = TRUE)
    },
    on_progress = function(finished) {
      if (isTRUE(show_progress)) {
        section6_progress(nrow(jobs) - nrow(pending) + finished, nrow(jobs), started, existing, cores)
      }
    }
  )
  invisible(list(results = existing, summary = summarize_section6_results(existing)))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  family <- tolower(as.character(args$family %||% "vmf"))
  M <- as.integer(args$M %||% 1000L)
  B <- as.integer(args$B %||% 5000L)
  output_dir <- as.character(args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios",
    sprintf("final_%s_d2_d10_M%d_B%d", family, M, B)
  ))
  common <- list(
    family = family,
    output_dir = output_dir,
    M = M,
    B = B,
    dimensions = parse_section6_csv(args$dimensions, c(2L, 10L), "integer"),
    n_values = parse_section6_csv(args$n_values, c(50L, 100L, 200L, 400L), "integer"),
    beta_values = parse_section6_csv(args$beta_values, c(0, 0.5, 1), "numeric"),
    cores = as.integer(args$cores %||% 8L),
    base_seed = as.integer(args$seed %||% 20260727L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 10000L),
    cvm_block_size = as.integer(args$cvm_block_size %||% 50L),
    derivative_method = tolower(as.character(args$derivative_method %||% "score_mc")),
    scenarios = if (is.null(args$scenarios)) NULL else parse_section6_csv(args$scenarios, NULL, "character")
  )
  if (!is.null(args$recover_from)) {
    recovery_common <- common
    recovery_common$output_dir <- NULL
    do.call(recover_section6_results, c(
      list(source_dir = as.character(args$recover_from), target_dir = output_dir), recovery_common
    ))
  } else {
    do.call(run_section6_family, c(common, list(
      checkpoint_results = as.integer(args$checkpoint_results %||% 64L),
      show_progress = !identical(tolower(as.character(args$show_progress %||% "true")), "false")
    )))
  }
}
