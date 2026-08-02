#!/usr/bin/env Rscript

# Accuracy and memory prototype for invariant vMF/HvMF profile lookup tables.
# The production profile code is not changed by this script.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("bootstrap/multiplier_bootstrap.R")
source("bootstrap/profile_lookup_interpolation.R")

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  args <- commandArgs(trailingOnly = TRUE)
  matches <- args[startsWith(args, prefix)]
  if (!length(matches)) return(default)
  substring(matches[[length(matches)]], nchar(prefix) + 1L)
}

model <- tolower(parse_option("model", "vmf"))
q <- as.integer(parse_option("q", "10"))
kappa_min <- as.numeric(parse_option("kappa-min", if (model == "vmf") "7" else "7"))
kappa_max <- as.numeric(parse_option("kappa-max", if (model == "vmf") "14" else "14"))
geometry_max <- as.numeric(parse_option("geometry-max", if (model == "vmf") "1" else "3"))
t_max <- as.numeric(parse_option("t-max", if (model == "vmf") format(pi, digits = 17) else "4"))
n_kappa <- as.integer(parse_option("n-kappa", "17"))
n_geometry <- as.integer(parse_option("n-geometry", "65"))
n_t <- as.integer(parse_option("n-t", "513"))
n_queries <- as.integer(parse_option("queries", "200"))
cores <- as.integer(parse_option("cores", "1"))
stencil_size <- as.integer(parse_option("stencil-size", "4"))
integration_grid_size <- as.integer(parse_option("integration-grid-size", "4097"))
reference_grid_size <- as.integer(parse_option("reference-grid-size", as.character(integration_grid_size)))
seed <- as.integer(parse_option("seed", "20260801"))
geometry_grid_type <- tolower(parse_option("geometry-grid", "uniform"))
kappa_grid_type <- tolower(parse_option("kappa-grid", "uniform"))
output_path <- parse_option(
  "output",
  file.path(
    "benchmarks",
    sprintf("%s_q%d_profile_lookup_interpolation.csv", model, q)
  )
)

if (!model %in% c("vmf", "hvmf")) stop("`model` must be vmf or hvmf.")
if (!geometry_grid_type %in% c("uniform", "chebyshev") ||
    !kappa_grid_type %in% c("uniform", "chebyshev")) {
  stop("Grid types must be `uniform` or `chebyshev`.")
}
if (q < 2L || kappa_min <= 0 || kappa_max <= kappa_min ||
    geometry_max <= 0 || t_max <= 0 ||
    any(c(n_kappa, n_geometry, n_t) < stencil_size + 1L) ||
    n_queries < 1L || cores < 1L || !stencil_size %in% c(4L, 6L) ||
    integration_grid_size < 17L || reference_grid_size < 17L) {
  stop("Invalid interpolation benchmark arguments.")
}

canonical_state <- function(model, q, kappa, geometry) {
  if (model == "vmf") {
    rho <- geometry
    mu <- c(1, rep.int(0, q))
    omega <- c(rho, sqrt(pmax(0, 1 - rho^2)), rep.int(0, q - 1L))
  } else {
    chi <- geometry
    mu <- c(1, rep.int(0, q))
    omega <- c(cosh(chi), sinh(chi), rep.int(0, q - 1L))
  }
  list(mu = mu, omega = omega, xi = kappa * mu)
}

evaluate_direct <- function(model, q, kappa, geometry, t_values,
                            grid_size = integration_grid_size) {
  state <- canonical_state(model, q, kappa, geometry)
  if (model == "vmf") {
    value <- vmf_profile_and_derivative_xi(
      omega = state$omega,
      xi = state$xi,
      t_values = t_values,
      distance_type = "geodesic",
      grid_size = grid_size
    )
    moments <- cbind(F = value$F, second = value$M1, third = value$MC)
    coefficients <- cbind(
      F = value$F,
      omega = value$M1 - kappa * geometry * value$MC,
      mu = kappa * value$MC - A_q(kappa, q) * value$F
    )
  } else {
    value <- hvmf_profile_and_derivative_xi(
      omega = state$omega,
      xi = state$xi,
      t_values = t_values,
      grid_size = grid_size
    )
    moments <- cbind(F = value$F, second = value$N1, third = value$ND)
    coefficients <- cbind(
      F = value$F,
      omega = value$N1 - kappa * cosh(geometry) * value$ND,
      mu = kappa * value$ND -
        hvmf_mean_resultant_ratio(q, kappa) * value$F
    )
  }
  list(moments = moments, coefficients = coefficients, derivative = value$derivative)
}

build_surface <- function(model, q, kappa_grid, geometry_grid, t_grid) {
  dimensions <- c(length(t_grid), length(geometry_grid), length(kappa_grid))
  table <- profile_lookup_build(
    model = model,
    q = q,
    kappa_grid = kappa_grid,
    geometry_grid = geometry_grid,
    t_grid = t_grid,
    integration_grid_size = integration_grid_size,
    cores = cores
  )
  moments <- lapply(seq_len(3L), function(index) {
    array(table$values[, index], dim = dimensions)
  })
  list(
    moments = moments,
    build_seconds = table$build_seconds,
    bytes = table$bytes
  )
}

subsample_surface <- function(surface, t_grid, geometry_grid, kappa_grid,
                              t_indices, geometry_indices, kappa_indices) {
  list(
    moments = lapply(surface$moments, function(values) {
      values[t_indices, geometry_indices, kappa_indices, drop = FALSE]
    }),
    t_grid = t_grid[t_indices],
    geometry_grid = geometry_grid[geometry_indices],
    kappa_grid = kappa_grid[kappa_indices]
  )
}

evaluate_variant <- function(name, variant, queries, references) {
  n <- nrow(queries)
  lookup_values <- do.call(cbind, lapply(variant$moments, as.vector))
  invisible(distance_profile_cpp_call(
    "cpp_profile_lookup_tensor_local_polynomial",
    lookup_values,
    variant$t_grid,
    variant$geometry_grid,
    variant$kappa_grid,
    queries$t[[1L]],
    queries$geometry[[1L]],
    queries$kappa[[1L]],
    stencil_size
  ))
  started <- proc.time()[["elapsed"]]
  interpolated_coefficients <- distance_profile_cpp_call(
    "cpp_profile_lookup_tensor_local_polynomial",
    lookup_values,
    variant$t_grid,
    variant$geometry_grid,
    variant$kappa_grid,
    queries$t,
    queries$geometry,
    queries$kappa,
    stencil_size
  )
  derivative <- matrix(0, nrow = n, ncol = q + 1L)
  if (model == "vmf") {
    derivative[, 1L] <- interpolated_coefficients[, 2L] * queries$geometry +
      interpolated_coefficients[, 3L]
    derivative[, 2L] <- interpolated_coefficients[, 2L] *
      sqrt(pmax(0, 1 - queries$geometry^2))
  } else {
    derivative[, 1L] <- -interpolated_coefficients[, 2L] *
      cosh(queries$geometry) - interpolated_coefficients[, 3L]
    derivative[, 2L] <- interpolated_coefficients[, 2L] *
      sinh(queries$geometry)
  }
  evaluation_seconds <- proc.time()[["elapsed"]] - started
  F_error <- interpolated_coefficients[, 1L] - references$F
  derivative_difference <- derivative - references$derivative
  derivative_max_error <- apply(abs(derivative_difference), 1L, max)
  derivative_norm_error <- sqrt(rowSums(derivative_difference^2))
  data.frame(
    model = model,
    q = q,
    variant = name,
    n_kappa = length(variant$kappa_grid),
    n_geometry = length(variant$geometry_grid),
    n_t = length(variant$t_grid),
    variant_MiB = sum(vapply(variant$moments, object.size, numeric(1))) / 1024^2,
    max_abs_F_error = max(abs(F_error)),
    rmse_F = sqrt(mean(F_error^2)),
    max_component_derivative_error = max(derivative_max_error),
    max_norm_derivative_error = max(derivative_norm_error),
    rmse_derivative_norm = sqrt(mean(derivative_norm_error^2)),
    evaluation_microseconds_per_query = 1e6 * evaluation_seconds / n,
    evaluation_backend = "cpp_local_polynomial"
  )
}

make_grid <- function(lower, upper, n, type) {
  if (type == "uniform") return(seq(lower, upper, length.out = n))
  midpoint <- (lower + upper) / 2
  half_width <- (upper - lower) / 2
  midpoint + half_width * sort(cos(pi * (0:(n - 1L)) / (n - 1L)))
}

kappa_grid <- make_grid(kappa_min, kappa_max, n_kappa, kappa_grid_type)
geometry_grid <- if (model == "vmf") {
  make_grid(-geometry_max, geometry_max, n_geometry, geometry_grid_type)
} else {
  make_grid(0, geometry_max, n_geometry, geometry_grid_type)
}
t_grid <- seq(0, t_max, length.out = n_t)

surface <- build_surface(model, q, kappa_grid, geometry_grid, t_grid)

set.seed(seed)
queries <- data.frame(
  kappa = stats::runif(n_queries, kappa_min, kappa_max),
  geometry = if (model == "vmf") {
    stats::runif(n_queries, -0.995 * geometry_max, 0.995 * geometry_max)
  } else {
    stats::runif(n_queries, 0.002 * geometry_max, 0.995 * geometry_max)
  },
  t = stats::runif(n_queries, 0.002 * t_max, 0.998 * t_max)
)
references <- list(F = numeric(n_queries), derivative = matrix(0, n_queries, q + 1L))
reference_started <- proc.time()[["elapsed"]]
for (i in seq_len(n_queries)) {
  value <- evaluate_direct(
    model, q, queries$kappa[[i]], queries$geometry[[i]], queries$t[[i]],
    grid_size = reference_grid_size
  )
  references$F[[i]] <- value$moments[1L, 1L]
  references$derivative[i, ] <- value$derivative[1L, ]
}
reference_evaluation_seconds <- proc.time()[["elapsed"]] - reference_started

all_indices <- list(
  t = seq_along(t_grid),
  geometry = seq_along(geometry_grid),
  kappa = seq_along(kappa_grid)
)
coarse_indices <- lapply(all_indices, function(indices) {
  unique(c(seq.int(1L, length(indices), by = 2L), length(indices)))
})
variants <- list(
  full = subsample_surface(
    surface, t_grid, geometry_grid, kappa_grid,
    all_indices$t, all_indices$geometry, all_indices$kappa
  ),
  coarse_kappa = subsample_surface(
    surface, t_grid, geometry_grid, kappa_grid,
    all_indices$t, all_indices$geometry, coarse_indices$kappa
  ),
  coarse_geometry = subsample_surface(
    surface, t_grid, geometry_grid, kappa_grid,
    all_indices$t, coarse_indices$geometry, all_indices$kappa
  ),
  coarse_t = subsample_surface(
    surface, t_grid, geometry_grid, kappa_grid,
    coarse_indices$t, all_indices$geometry, all_indices$kappa
  ),
  coarse_all = subsample_surface(
    surface, t_grid, geometry_grid, kappa_grid,
    coarse_indices$t, coarse_indices$geometry, coarse_indices$kappa
  )
)

results <- do.call(rbind, lapply(names(variants), function(name) {
  evaluate_variant(name, variants[[name]], queries, references)
}))
results$surface_build_seconds <- surface$build_seconds
results$surface_MiB <- as.numeric(surface$bytes) / 1024^2
results$kappa_min <- kappa_min
results$kappa_max <- kappa_max
results$geometry_min <- min(geometry_grid)
results$geometry_max <- max(geometry_grid)
results$t_max <- t_max
results$n_queries <- n_queries
results$cores <- cores
results$kappa_grid_type <- kappa_grid_type
results$geometry_grid_type <- geometry_grid_type
results$stencil_size <- stencil_size
results$integration_grid_size <- integration_grid_size
results$reference_grid_size <- reference_grid_size
results$reference_evaluation_seconds <- reference_evaluation_seconds
results$reference_microseconds_per_query <-
  1e6 * reference_evaluation_seconds / n_queries

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output_path, row.names = FALSE)
print(results, row.names = FALSE)
