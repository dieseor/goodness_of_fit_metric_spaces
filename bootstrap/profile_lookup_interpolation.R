# Experimental invariant lookup tables for deterministic vMF/HvMF profiles.
# This file is intentionally not wired into the production method selector yet.

# This module uses the C++ interpolation backend directly.  Load that narrow
# dependency explicitly rather than relying on another sourced module (or on
# the order in which test files happen to run) having loaded it already.
# Loading the backend only defines helpers; compilation remains lazy inside
# `distance_profile_cpp_call()`.
if (!exists("distance_profile_cpp_call", mode = "function") ||
    !exists("with_distance_profile_cpp_cache_lock", mode = "function")) {
  distance_profile_backend_candidates <- c(
    "distance_profile_backend.R",
    file.path("..", "distance_profile_backend.R"),
    file.path("..", "..", "distance_profile_backend.R")
  )
  distance_profile_backend_path <- distance_profile_backend_candidates[
    file.exists(distance_profile_backend_candidates)
  ][1L]
  if (is.na(distance_profile_backend_path)) {
    stop("Could not locate `distance_profile_backend.R` for profile lookup.")
  }
  source(distance_profile_backend_path)
  rm(distance_profile_backend_candidates, distance_profile_backend_path)
}

profile_lookup_build <- function(model,
                                 q,
                                 kappa_grid,
                                 geometry_grid,
                                 t_grid,
                                 integration_grid_size = 16385L,
                                 cores = 1L) {
  model <- match.arg(tolower(model), c("vmf", "hvmf"))
  q <- as.integer(q)
  integration_grid_size <- as.integer(integration_grid_size)
  cores <- as.integer(cores)
  validate_grid <- function(x, label) {
    x <- as.numeric(x)
    if (length(x) < 6L || any(!is.finite(x)) || any(diff(x) <= 0)) {
      stop(sprintf("%s must contain at least six finite, strictly increasing nodes.", label))
    }
    x
  }
  kappa_grid <- validate_grid(kappa_grid, "`kappa_grid`")
  geometry_grid <- validate_grid(geometry_grid, "`geometry_grid`")
  t_grid <- validate_grid(t_grid, "`t_grid`")
  if (q < 2L || kappa_grid[[1L]] <= 0 || t_grid[[1L]] < 0 ||
      integration_grid_size < 17L || cores < 1L) {
    stop("Invalid profile lookup-table arguments.")
  }
  if (identical(model, "vmf") &&
      (geometry_grid[[1L]] < -1 || tail(geometry_grid, 1L) > 1)) {
    stop("The vMF geometry grid must be contained in [-1, 1].")
  }
  if (identical(model, "hvmf") && geometry_grid[[1L]] < 0) {
    stop("The HvMF geometry grid must be nonnegative.")
  }

  evaluate_node <- function(kappa, geometry) {
    mu <- c(1, rep.int(0, q))
    if (identical(model, "vmf")) {
      omega <- c(
        geometry,
        sqrt(pmax(0, 1 - geometry^2)),
        rep.int(0, q - 1L)
      )
      value <- vmf_profile_and_derivative_xi(
        omega = omega,
        xi = kappa * mu,
        t_values = t_grid,
        distance_type = "geodesic",
        grid_size = integration_grid_size
      )
      return(cbind(
        F = value$F,
        coefficient_omega = value$M1 - kappa * geometry * value$MC,
        coefficient_mu = kappa * value$MC - A_q(kappa, q) * value$F
      ))
    }

    omega <- c(cosh(geometry), sinh(geometry), rep.int(0, q - 1L))
    value <- hvmf_profile_and_derivative_xi(
      omega = omega,
      xi = kappa * mu,
      t_values = t_grid,
      grid_size = integration_grid_size
    )
    cbind(
      F = value$F,
      coefficient_omega = value$N1 -
        kappa * cosh(geometry) * value$ND,
      coefficient_mu = kappa * value$ND -
        hvmf_mean_resultant_ratio(q, kappa) * value$F
    )
  }

  started <- proc.time()[["elapsed"]]
  build_slice <- function(kappa) {
    do.call(rbind, lapply(geometry_grid, function(geometry) {
      evaluate_node(kappa, geometry)
    }))
  }
  slices <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(
      kappa_grid,
      build_slice,
      mc.cores = min(cores, length(kappa_grid)),
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(kappa_grid, build_slice)
  }
  values <- do.call(rbind, slices)
  expected_rows <- length(t_grid) * length(geometry_grid) * length(kappa_grid)
  if (nrow(values) != expected_rows || ncol(values) != 3L ||
      any(!is.finite(values))) {
    stop("Profile lookup-table construction returned incompatible values.")
  }

  structure(list(
    model = model,
    q = q,
    kappa_grid = kappa_grid,
    geometry_grid = geometry_grid,
    t_grid = t_grid,
    values = unname(values),
    stencil_size = 6L,
    integration_grid_size = integration_grid_size,
    build_seconds = proc.time()[["elapsed"]] - started,
    bytes = as.numeric(object.size(values))
  ), class = "profile_lookup_table")
}

profile_lookup_evaluate <- function(table, xi, centers, thresholds) {
  if (!inherits(table, "profile_lookup_table")) {
    stop("`table` must be a profile lookup table.")
  }
  centers <- as.matrix(centers)
  xi <- as.numeric(xi)
  if (ncol(centers) != table$q + 1L || length(xi) != table$q + 1L) {
    stop("The centers, parameter, and lookup table have incompatible dimensions.")
  }
  n_centers <- nrow(centers)
  if (is.null(dim(thresholds))) {
    thresholds <- matrix(
      rep.int(as.numeric(thresholds), n_centers),
      nrow = n_centers,
      byrow = TRUE
    )
  } else {
    thresholds <- as.matrix(thresholds)
    if (nrow(thresholds) != n_centers) {
      stop("The threshold matrix must have one row per center.")
    }
  }
  n_thresholds <- ncol(thresholds)

  if (identical(table$model, "vmf")) {
    kappa <- sqrt(sum(xi^2))
    mu <- xi / kappa
    center_norms <- sqrt(rowSums(centers^2))
    if (any(abs(center_norms - 1) > 1e-8)) {
      stop("vMF lookup centers must lie on the unit sphere.")
    }
    geometry <- drop(centers %*% mu)
    geometry <- pmax(-1, pmin(1, geometry))
  } else {
    kappa_squared <- xi[[1L]]^2 - sum(xi[-1L]^2)
    if (!is.finite(kappa_squared) || kappa_squared <= 0 || xi[[1L]] <= 0) {
      stop("HvMF lookup requires a future-directed timelike `xi`.")
    }
    kappa <- sqrt(kappa_squared)
    mu <- xi / kappa
    minkowski_center_mu <- -centers[, 1L] * mu[[1L]] +
      drop(centers[, -1L, drop = FALSE] %*% mu[-1L])
    geometry <- acosh(pmax(1, -minkowski_center_mu))
  }
  if (!is.finite(kappa) || kappa <= 0) {
    stop("Profile lookup requires a strictly positive concentration.")
  }

  t_query <- as.vector(t(thresholds))
  geometry_query <- rep(geometry, each = n_thresholds)
  coefficients <- distance_profile_cpp_call(
    "cpp_profile_lookup_tensor_local_polynomial",
    table$values,
    table$t_grid,
    table$geometry_grid,
    table$kappa_grid,
    t_query,
    geometry_query,
    rep.int(kappa, length(t_query)),
    table$stencil_size
  )
  center_rows <- rep(seq_len(n_centers), each = n_thresholds)
  repeated_mu <- matrix(mu, nrow = length(t_query), ncol = length(mu), byrow = TRUE)
  repeated_centers <- centers[center_rows, , drop = FALSE]
  if (identical(table$model, "vmf")) {
    derivative <- coefficients[, 2L] * repeated_centers +
      coefficients[, 3L] * repeated_mu
  } else {
    repeated_centers[, 1L] <- -repeated_centers[, 1L]
    repeated_mu[, 1L] <- -repeated_mu[, 1L]
    derivative <- coefficients[, 2L] * repeated_centers +
      coefficients[, 3L] * repeated_mu
  }

  list(
    F = matrix(coefficients[, 1L], nrow = n_centers, byrow = TRUE),
    derivative_sorted = derivative,
    coefficients = coefficients,
    kappa = kappa,
    geometry = geometry
  )
}
