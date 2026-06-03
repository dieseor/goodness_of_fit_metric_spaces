

#!/usr/bin/env Rscript

# Validation script for beta_mixture2 Legendre coefficients and profiles.
# Goal: check the new Gauss--Jacobi coefficient computation against
# exact identities, integral fallback profiles, and special cases omega = +/- mu.


set.seed(123)

n_cores <- as.integer(Sys.getenv("N_CORES", "12"))
if (!is.finite(n_cores) || n_cores < 1L) {
  n_cores <- 1L
}

source("utils.R")
source("bootstrap/calibration_study.R")

norm_vec <- function(x) x / sqrt(sum(x^2))
max_abs <- function(x) max(abs(x), na.rm = TRUE)

make_theta <- function(mu = c(0.2, -0.3, 0.93),
                       weight1 = 0.4,
                       alpha1 = 2,
                       beta1 = 8,
                       alpha2 = 8,
                       beta2 = 2) {
  beta_mixture2_normalize_theta(list(
    mu = norm_vec(mu),
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  ), ambient_dim = 3L)
}

random_unit <- function() norm_vec(rnorm(3))

orthogonal_unit <- function(mu) {
  v <- c(1, 0, 0)
  if (abs(sum(v * mu)) > 0.9) v <- c(0, 1, 0)
  v <- v - sum(v * mu) * mu
  norm_vec(v)
}

# Old generic Gauss--Legendre coefficient computation, kept only as a reference
# in easy non-singular cases. This is not used as the production method.
coefficients_legendre_generic_reference <- function(theta, l_max = 60L, quad_n = 5000L, tol = 1e-3) {
  rotational_legendre_coefficients(
    density_h = function(z) {
      beta_mixture2_density_h(
        z = z,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2
      )
    },
    Lmax = l_max,
    quad_n = quad_n,
    tol = tol
  )
}

validate_coefficients <- function(theta, l_max = 120L, quad_n = 100L, tol = 1e-10) {
  cat("\n--- Coefficient validation ---\n")
  print(theta[c("weight1", "alpha1", "beta1", "alpha2", "beta2")])

  gj <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )

  cat("method:", gj$method, "\n")
  cat("quad_n used:", gj$quad_n, "\n")
  cat("a0_error:", format(gj$a0_error, scientific = TRUE), "\n")
  cat("mass_error:", format(gj$mass_error, scientific = TRUE), "\n")
  stopifnot(is.finite(gj$a0_error), gj$a0_error < tol)
  stopifnot(abs(gj$coefficients[[1L]] - 1) < tol)

  invisible(gj)
}

validate_against_generic_reference <- function(theta, l_max = 40L) {
  cat("\n--- Gauss--Jacobi vs generic Gauss--Legendre reference ---\n")
  cat("This check is meaningful mainly for non-extreme beta shapes.\n")

  time_gj <- system.time({
    gj <- beta_mixture2_legendre_coefficients(theta, l_max = l_max, quad_n = 100L)
  })
  time_gl <- system.time({
    ref <- coefficients_legendre_generic_reference(theta, l_max = l_max, quad_n = 3000L, tol = 1e-2)
  })

  diff <- max_abs(gj$coefficients - ref$coefficients)

  cat("max abs coefficient difference:", format(diff, scientific = TRUE), "\n")
  cat("Gauss--Jacobi elapsed:", sprintf("%.4f", unname(time_gj[["elapsed"]])), "sec\n")
  cat("Gauss--Legendre reference elapsed:", sprintf("%.4f", unname(time_gl[["elapsed"]])), "sec\n")
  cat("time ratio GL/GJ:", sprintf("%.2f", unname(time_gl[["elapsed"]]) / max(unname(time_gj[["elapsed"]]), 1e-12)), "\n")

  invisible(list(
    max_abs_diff = diff,
    time_gauss_jacobi = unname(time_gj[["elapsed"]]),
    time_gauss_legendre = unname(time_gl[["elapsed"]])
  ))
}

benchmark_coefficients_time <- function(theta, l_max = 150L, gj_quad_n = 100L, gl_quad_n = 3000L) {
  cat("\n--- Coefficient timing benchmark ---\n")

  time_gj <- system.time({
    gj <- beta_mixture2_legendre_coefficients(theta, l_max = l_max, quad_n = gj_quad_n)
  })
  time_gl <- system.time({
    gl <- coefficients_legendre_generic_reference(theta, l_max = l_max, quad_n = gl_quad_n, tol = 1e-2)
  })

  diff <- max_abs(gj$coefficients - gl$coefficients)
  cat("l_max:", l_max, "\n")
  cat("Gauss--Jacobi quad_n:", gj_quad_n, "elapsed:", sprintf("%.4f", unname(time_gj[["elapsed"]])), "sec\n")
  cat("Gauss--Legendre quad_n:", gl_quad_n, "elapsed:", sprintf("%.4f", unname(time_gl[["elapsed"]])), "sec\n")
  cat("time ratio GL/GJ:", sprintf("%.2f", unname(time_gl[["elapsed"]]) / max(unname(time_gj[["elapsed"]]), 1e-12)), "\n")
  cat("max abs coefficient difference:", format(diff, scientific = TRUE), "\n")

  invisible(list(
    time_gauss_jacobi = unname(time_gj[["elapsed"]]),
    time_gauss_legendre = unname(time_gl[["elapsed"]]),
    max_abs_diff = diff
  ))
}

validate_special_cases <- function(theta, l_max = 120L, quad_n = 100L) {
  cat("\n--- Special cases omega = +/- mu ---\n")

  t_grid <- seq(0, pi, length.out = 101)
  c_grid <- cos(t_grid)

  prof_plus <- distance_profile_beta_mixture2(
    omega = theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    distance_type = "geodesic",
    method = "legendre",
    l_max = l_max,
    quad_n = quad_n
  )
  closed_plus <- 1 - beta_mixture2_cdf_y(
    y = (c_grid + 1) / 2,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )

  prof_minus <- distance_profile_beta_mixture2(
    omega = -theta$mu,
    t_values = t_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    distance_type = "geodesic",
    method = "legendre",
    l_max = l_max,
    quad_n = quad_n
  )
  closed_minus <- beta_mixture2_cdf_y(
    y = (1 - c_grid) / 2,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )

  err_plus <- max_abs(prof_plus - closed_plus)
  err_minus <- max_abs(prof_minus - closed_minus)

  cat("max abs error omega=mu:", format(err_plus, scientific = TRUE), "\n")
  cat("max abs error omega=-mu:", format(err_minus, scientific = TRUE), "\n")

  stopifnot(err_plus < 1e-10)
  stopifnot(err_minus < 1e-10)

  invisible(c(plus = err_plus, minus = err_minus))
}

validate_profile_vs_integral <- function(theta,
                                         l_max = 120L,
                                         quad_n_coeff = 100L,
                                         quad_n_integral = 1200L,
                                         tolerance = 5e-4) {
  cat("\n--- Legendre profile vs integral fallback ---\n")

  omega_list <- list(
    mu = theta$mu,
    minus_mu = -theta$mu,
    orthogonal = orthogonal_unit(theta$mu),
    random1 = random_unit(),
    random2 = random_unit()
  )
  t_grid <- seq(0, pi, length.out = 81)

  rows <- lapply(names(omega_list), function(nm) {
    omega <- omega_list[[nm]]
    leg <- distance_profile_beta_mixture2(
      omega = omega,
      t_values = t_grid,
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      distance_type = "geodesic",
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n_coeff
    )
    integ <- distance_profile_beta_mixture2(
      omega = omega,
      t_values = t_grid,
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      distance_type = "geodesic",
      method = "integral",
      quad_n = quad_n_integral
    )

    data.frame(
      omega = nm,
      max_abs_diff = max_abs(leg - integ),
      mean_abs_diff = mean(abs(leg - integ)),
      monotone_legendre = all(diff(leg) >= -1e-10),
      range_ok = min(leg) >= -1e-10 && max(leg) <= 1 + 1e-10,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  print(out)

  if (max(out$max_abs_diff) > tolerance) {
    stop(sprintf(
      "Legendre vs integral validation failed: max discrepancy %.3e exceeds %.3e.",
      max(out$max_abs_diff), tolerance
    ))
  }
  stopifnot(all(out$monotone_legendre), all(out$range_ok))

  invisible(out)
}

validate_grid_vs_scalar <- function(theta, l_max = 120L, quad_n = 100L) {
  cat("\n--- Grid evaluator vs scalar evaluator ---\n")

    omega_grid <- rbind(
    orthogonal_unit(theta$mu),
    random_unit(),
    random_unit(),
    random_unit()
    )

  t_grid <- seq(0, pi, length.out = 41)

  grid <- distance_profile_beta_mixture2_grid(
    omega_grid = omega_grid,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    t_grid = t_grid,
    distance_type = "geodesic",
    method = "legendre",
    l_max = l_max,
    quad_n = quad_n
  )

  scalar <- t(vapply(seq_len(nrow(omega_grid)), function(i) {
    distance_profile_beta_mixture2(
      omega = omega_grid[i, ],
      t_values = t_grid,
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      distance_type = "geodesic",
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n
    )
  }, numeric(length(t_grid))))

  diff <- max_abs(grid - scalar)
  cat("max abs grid-scalar diff:", format(diff, scientific = TRUE), "\n")
  stopifnot(diff < 1e-8)

  invisible(diff)
}

validate_cvm_grid_vs_naive <- function(theta, n = 25L, l_max = 120L, quad_n = 100L) {
  cat("\n--- CvM matrix evaluator vs naive scalar loop ---\n")

  X <- r_sph_beta_mixture2(
    n = n,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  t_matrix <- acos(dot_products)

  fast <- distance_profile_beta_mixture2_cvm_grid(
    X = X,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    method = "legendre",
    l_max = l_max,
    quad_n = quad_n
  )

  naive <- matrix(NA_real_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    naive[i, ] <- distance_profile_beta_mixture2(
      omega = X[i, ],
      t_values = t_matrix[i, ],
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      distance_type = "geodesic",
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n
    )
  }

  diff <- max_abs(fast - naive)
  cat("max abs CvM fast-naive diff:", format(diff, scientific = TRUE), "\n")
  stopifnot(diff < 1e-12)

  invisible(diff)
}

# Test scenarios. Include the calibration scenario and more edge-heavy shapes.
theta_list <- list(
  calibration = make_theta(weight1 = 0.4, alpha1 = 2, beta1 = 8, alpha2 = 8, beta2 = 2),
  mild = make_theta(weight1 = 0.55, alpha1 = 3, beta1 = 5, alpha2 = 6, beta2 = 2.5),
  edge_left_right = make_theta(weight1 = 0.5, alpha1 = 0.6, beta1 = 4, alpha2 = 5, beta2 = 0.7),
  concentrated = make_theta(weight1 = 0.35, alpha1 = 20, beta1 = 3, alpha2 = 3, beta2 = 18)
)

run_validation_scenario <- function(nm) {
  capture.output({
    cat("\n============================================================\n")
    cat("Scenario:", nm, "\n")
    cat("============================================================\n")

    theta <- theta_list[[nm]]
    validate_coefficients(theta, l_max = 150L, quad_n = 100L, tol = 1e-10)
    validate_special_cases(theta, l_max = 150L, quad_n = 100L)
    validate_profile_vs_integral(theta, l_max = 150L, quad_n_coeff = 100L, quad_n_integral = 1500L)
    validate_grid_vs_scalar(theta, l_max = 150L, quad_n = 100L)
    validate_cvm_grid_vs_naive(theta, n = 20L, l_max = 150L, quad_n = 100L)

    if (nm %in% c("calibration", "mild")) {
      validate_against_generic_reference(theta, l_max = 30L)
      benchmark_coefficients_time(theta, l_max = 150L, gj_quad_n = 100L, gl_quad_n = 3000L)
    }
  })
}

scenario_names <- names(theta_list)
n_cores_used <- n_cores
cat("Running", length(scenario_names), "validation scenarios with", n_cores_used, "cores.\n")

if (.Platform$OS.type == "unix" && n_cores_used > 1L) {
  outputs <- parallel::mclapply(
    scenario_names,
    run_validation_scenario,
    mc.cores = n_cores_used,
    mc.preschedule = FALSE
  )
} else if (n_cores_used > 1L) {
  cl <- parallel::makeCluster(n_cores_used)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(
    cl,
    varlist = setdiff(ls(envir = environment()), c("cl")),
    envir = environment()
  )
  outputs <- parallel::parLapply(cl, scenario_names, run_validation_scenario)
} else {
  outputs <- lapply(scenario_names, run_validation_scenario)
}

has_error <- vapply(outputs, inherits, logical(1), what = "try-error")

for (i in seq_along(outputs)) {
  out <- outputs[[i]]
  if (inherits(out, "try-error")) {
    cat("\n============================================================\n")
    cat("Scenario failed:", scenario_names[[i]], "\n")
    cat("============================================================\n")
    cat(as.character(out), "\n", sep = "")
  } else {
    cat(paste(out, collapse = "\n"), "\n", sep = "")
  }
}

if (any(has_error)) {
  stop(sprintf("%d beta-mixture Gauss--Jacobi validation scenario(s) failed.", sum(has_error)))
}

cat("\nAll beta-mixture Gauss--Jacobi coefficient/profile validations passed.\n")