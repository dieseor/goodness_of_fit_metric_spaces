resolve_small_circle_distance_profile_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

utils_path_small_circle_distance_profile <- resolve_small_circle_distance_profile_path("utils.R")
if (!exists("distance_profile_small_circle", mode = "function")) {
  source(utils_path_small_circle_distance_profile)
}

run_small_circle_distance_profile_validation <- function(output_root = file.path("output", "small_circle_validation"),
                                                         l_max = 200L,
                                                         quad_n = 500L,
                                                         seed = 20260531L) {
  set.seed(as.integer(seed))
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  mu <- c(0, 0, 1)
  complement <- jp_orthonormal_complement(mu)
  omega_list <- list(
    mu,
    -mu,
    complement[, 1L],
    jp_normalize_unit_vector(c(0.3, -0.4, 0.85), arg_name = "omega", min_length = 3L)
  )
  names(omega_list) <- c("mu", "minus_mu", "orthogonal", "random")
  t_grid <- seq(0.05, pi - 0.05, length.out = 101)
  parameter_grid <- expand.grid(
    kappa = c(0, 1, 5, 20, 50),
    nu = c(0, 0.3, 0.7, 0.9),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  summary_rows <- vector("list", nrow(parameter_grid))
  for (i in seq_len(nrow(parameter_grid))) {
    comparison <- small_circle_compare_profile_methods(
      mu = mu,
      kappa = parameter_grid$kappa[[i]],
      nu = parameter_grid$nu[[i]],
      omega_list = omega_list,
      t_grid = t_grid,
      distance_type = "geodesic",
      l_max = l_max,
      quad_n = quad_n,
      tol = 1e-10
    )
    comparison$kappa <- parameter_grid$kappa[[i]]
    comparison$nu <- parameter_grid$nu[[i]]
    comparison$omega_label <- names(omega_list)[comparison$omega_id]
    summary_rows[[i]] <- comparison[, c("kappa", "nu", "omega_id", "omega_label", "max_abs_diff", "mean_abs_diff")]
  }

  summary_df <- do.call(rbind, summary_rows)
  utils::write.csv(summary_df, file = file.path(output_root, "series_vs_integral_summary.csv"), row.names = FALSE)
  saveRDS(summary_df, file = file.path(output_root, "series_vs_integral_summary.rds"))

  checks_df <- data.frame(
    check = c("a0_uniform", "profile_uniform", "omega_plus_mu", "omega_minus_mu"),
    passed = c(
      abs(small_circle_legendre_coefficients(0, 0.3, l_max = 8L)$coefficients[[1L]] - 1) <= 1e-14,
      max(abs(distance_profile_small_circle(mu, t_grid, mu, 0, 0.4) - (1 - cos(t_grid)) / 2)) <= 1e-12,
      max(abs(distance_profile_small_circle(mu, t_grid, mu, 5, 0.3) - (1 - small_circle_axis_cdf(cos(t_grid), 5, 0.3)))) <= 5e-8,
      max(abs(distance_profile_small_circle(-mu, t_grid, mu, 5, 0.3) - small_circle_axis_cdf(-cos(t_grid), 5, 0.3))) <= 5e-8
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(checks_df, file = file.path(output_root, "theoretical_checks.csv"), row.names = FALSE)

  list(summary = summary_df, checks = checks_df, output_root = output_root)
}

if (sys.nframe() == 0L) {
  run_small_circle_distance_profile_validation()
}