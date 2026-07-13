#!/usr/bin/env Rscript

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "utils.R")
}
source(utils_path)

avg_time <- function(fun, n_rep) {
  gc(FALSE)
  system.time(for (i in seq_len(n_rep)) fun())[["elapsed"]] / n_rep
}

run_vmf_legendre_benchmark <- function() {
  set.seed(7)
  mu <- rnorm(3)
  mu <- mu / sqrt(sum(mu^2))
  omega <- rnorm(3)
  omega <- omega / sqrt(sum(omega^2))
  omega_grid <- matrix(rnorm(2000 * 3), ncol = 3)
  omega_grid <- omega_grid / sqrt(rowSums(omega_grid^2))
  x <- matrix(rnorm(800 * 3), ncol = 3)
  x <- x / sqrt(rowSums(x^2))
  t_grid <- seq(0, pi, length.out = 1001)

  scenarios <- list(
    list(kappa = 2, l_max = 10L),
    list(kappa = 10, l_max = 20L)
  )

  rows <- lapply(scenarios, function(scenario) {
    kappa <- scenario$kappa
    l_max <- scenario$l_max

    grid_error <- max(abs(
      distance_profile_vmf_s2_grid(
        omega_grid = omega_grid,
        mu = mu,
        kappa = kappa,
        t_grid = t_grid,
        distance_type = "geodesic",
        n_u = 4097L
      ) -
        distance_profile_vmf_s2_legendre_grid(
          omega_grid = omega_grid,
          mu = mu,
          kappa = kappa,
          t_grid = t_grid,
          distance_type = "geodesic",
          l_max = l_max
        )
    ))

    cvm_error <- max(abs(
      distance_profile_vmf_s2_cvm_grid(
        X = x,
        mu = mu,
        kappa = kappa,
        n_u = 4097L
      ) -
        distance_profile_vmf_s2_legendre_cvm_grid(
          X = x,
          mu = mu,
          kappa = kappa,
          l_max = l_max
        )
    ))

    scalar_exact <- avg_time(function() {
      theoretical_distance_profile_vmf_s2_fast(
        omega = omega,
        mu = mu,
        kappa = kappa,
        t_values = t_grid,
        distance_type = "geodesic"
      )
    }, n_rep = 3L)
    scalar_legendre <- avg_time(function() {
      distance_profile_vmf_s2_legendre(
        omega = omega,
        mu = mu,
        kappa = kappa,
        t_values = t_grid,
        distance_type = "geodesic",
        l_max = l_max
      )
    }, n_rep = 20L)
    grid_tabulated <- avg_time(function() {
      distance_profile_vmf_s2_grid(
        omega_grid = omega_grid,
        mu = mu,
        kappa = kappa,
        t_grid = t_grid,
        distance_type = "geodesic",
        n_u = 4097L
      )
    }, n_rep = 1L)
    grid_legendre <- avg_time(function() {
      distance_profile_vmf_s2_legendre_grid(
        omega_grid = omega_grid,
        mu = mu,
        kappa = kappa,
        t_grid = t_grid,
        distance_type = "geodesic",
        l_max = l_max
      )
    }, n_rep = 5L)
    cvm_tabulated <- avg_time(function() {
      distance_profile_vmf_s2_cvm_grid(
        X = x,
        mu = mu,
        kappa = kappa,
        n_u = 4097L
      )
    }, n_rep = 1L)
    cvm_legendre <- avg_time(function() {
      distance_profile_vmf_s2_legendre_cvm_grid(
        X = x,
        mu = mu,
        kappa = kappa,
        l_max = l_max
      )
    }, n_rep = 5L)

    data.frame(
      kappa = kappa,
      l_max = l_max,
      scalar_exact_sec = scalar_exact,
      scalar_legendre_sec = scalar_legendre,
      scalar_speedup = scalar_exact / scalar_legendre,
      grid_tabulated_sec = grid_tabulated,
      grid_legendre_sec = grid_legendre,
      grid_speedup = grid_tabulated / grid_legendre,
      grid_max_abs_error = grid_error,
      cvm_tabulated_sec = cvm_tabulated,
      cvm_legendre_sec = cvm_legendre,
      cvm_speedup = cvm_tabulated / cvm_legendre,
      cvm_max_abs_error = cvm_error,
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}

results <- run_vmf_legendre_benchmark()
print(results, digits = 6, row.names = FALSE)
