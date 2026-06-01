suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(gridExtra)
})

source(file.path("distance_profiles", "distance_profile_analysis.R"))
source(file.path("utils.R"))

compute_geodesic_distances_spherical_cauchy <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1L)
  }
  apply(data, 1, function(x) sphere_distance(omega, x, distance_type = "geodesic"))
}

compute_chordal_distances_spherical_cauchy <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1L)
  }
  apply(data, 1, function(x) sphere_distance(omega, x, distance_type = "chordal"))
}

spherical_cauchy_reference_omegas <- function(mu, seed = 123L) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  complement <- jp_orthonormal_complement(mu)
  orthogonal <- as.numeric(complement[, 1L])

  set.seed(as.integer(seed))
  generic <- stats::rnorm(length(mu))
  generic <- generic / sqrt(sum(generic^2))

  list(
    mu = mu,
    minus_mu = -mu,
    orthogonal = orthogonal,
    generic = generic
  )
}

spherical_cauchy_distance_profile_validation_summary <- function(sample,
                                                                 omega_values,
                                                                 mu,
                                                                 rho,
                                                                 t_grid,
                                                                 distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_names <- names(omega_values)
  distance_fun <- if (identical(distance_type, "geodesic")) {
    compute_geodesic_distances_spherical_cauchy
  } else {
    compute_chordal_distances_spherical_cauchy
  }

  summary_rows <- lapply(seq_along(omega_values), function(i) {
    omega <- omega_values[[i]]
    empirical <- empirical_distance_profile(
      distances = distance_fun(omega, sample),
      t_values = t_grid
    )
    theoretical <- distance_profile_spherical_cauchy(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      rho = rho,
      distance_type = distance_type,
      warn = FALSE
    )

    data.frame(
      rho = as.numeric(rho),
      distance_type = distance_type,
      omega_id = i,
      omega_name = omega_names[[i]],
      n = nrow(sample),
      grid_size = length(t_grid),
      max_abs_diff = max(abs(empirical - theoretical)),
      mean_abs_diff = mean(abs(empirical - theoretical)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_rows)
}

run_spherical_cauchy_distance_profile_analysis <- function(rho,
                                                           distance_type = c("geodesic", "chordal"),
                                                           output_dir = NULL,
                                                           mu = c(0, 0, 1),
                                                           sample_sizes = c(50, 200),
                                                           n_simulations = 10,
                                                           validation_n = 5000L,
                                                           n_points = 200L,
                                                           seed = 123L,
                                                           save_plots = TRUE) {
  distance_type <- match.arg(distance_type)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)

  rho_tag <- gsub("\\.", "p", format(as.numeric(rho), scientific = FALSE, trim = TRUE))
  if (is.null(output_dir)) {
    output_dir <- file.path("output", paste0("spherical_cauchy_", distance_type), paste0("rho_", rho_tag))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  omega_values <- spherical_cauchy_reference_omegas(mu = mu, seed = seed + 1L)
  omega_list <- unname(omega_values)
  legend_positions <- c("bottom_right", "top_left", "top_left", "top_left")

  data_generator <- function(n) {
    r_sph_spherical_cauchy(n = n, mu = mu, rho = rho)
  }

  theoretical_profile <- function(omega, t_values) {
    distance_profile_spherical_cauchy(
      omega = omega,
      t_values = t_values,
      mu = mu,
      rho = rho,
      distance_type = distance_type,
      warn = FALSE
    )
  }

  distance_fun <- if (identical(distance_type, "geodesic")) {
    compute_geodesic_distances_spherical_cauchy
  } else {
    compute_chordal_distances_spherical_cauchy
  }

  t_max <- if (identical(distance_type, "geodesic")) pi - 1e-8 else 2 - 1e-8

  plots <- create_distance_profile_analysis(
    omega_values = omega_list,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = distance_fun,
    sample_sizes = sample_sizes,
    n_simulations = n_simulations,
    t_max = t_max,
    save_plots = save_plots,
    output_dir = output_dir,
    file_prefix = paste0("spherical_cauchy_", distance_type, "_dp"),
    legend_positions = legend_positions
  )

  set.seed(seed + 2L)
  validation_sample <- r_sph_spherical_cauchy(
    n = validation_n,
    mu = mu,
    rho = rho
  )

  t_grid <- seq(0, t_max, length.out = n_points)
  summary_df <- spherical_cauchy_distance_profile_validation_summary(
    sample = validation_sample,
    omega_values = omega_values,
    mu = mu,
    rho = rho,
    t_grid = t_grid,
    distance_type = distance_type
  )

  utils::write.csv(
    summary_df,
    file = file.path(output_dir, "distance_profile_validation_summary.csv"),
    row.names = FALSE
  )

  invisible(list(
    output_dir = output_dir,
    plots = plots,
    validation_summary = summary_df
  ))
}

run_spherical_cauchy_distance_profile_suite <- function(rho_values = c(0.3, 0.7),
                                                        distance_types = c("geodesic", "chordal"),
                                                        seed = 123L,
                                                        sample_sizes = c(50, 200),
                                                        n_simulations = 10,
                                                        validation_n = 5000L,
                                                        n_points = 200L,
                                                        save_plots = TRUE) {
  results <- list()
  counter <- 1L

  for (distance_type in distance_types) {
    for (rho in rho_values) {
      results[[counter]] <- run_spherical_cauchy_distance_profile_analysis(
        rho = rho,
        distance_type = distance_type,
        sample_sizes = sample_sizes,
        n_simulations = n_simulations,
        validation_n = validation_n,
        n_points = n_points,
        seed = seed + counter,
        save_plots = save_plots
      )
      counter <- counter + 1L
    }
  }

  invisible(results)
}

if (sys.nframe() == 0L) {
  run_spherical_cauchy_distance_profile_suite()
}
