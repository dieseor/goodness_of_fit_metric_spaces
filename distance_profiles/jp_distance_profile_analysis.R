suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(gridExtra)
})

source(file.path("distance_profiles", "distance_profile_analysis.R"))
source(file.path("utils.R"))

format_jp_value_tag <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) {
    stop("`x` must be a finite scalar.")
  }

  tag <- format(x, scientific = FALSE, trim = TRUE)
  tag <- gsub("-", "neg", tag, fixed = TRUE)
  gsub("\\.", "p", tag)
}

jp_scenario_tag <- function(kappa, psi) {
  psi_tag <- format_jp_value_tag(psi)
  psi_prefix <- if (startsWith(psi_tag, "neg")) "psi_" else "psi"
  paste0(
    "kappa", format_jp_value_tag(kappa),
    "_", psi_prefix, psi_tag
  )
}

compute_geodesic_distances_jp <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1L)
  }
  apply(data, 1, function(x) sphere_distance(omega, x, distance_type = "geodesic"))
}

jp_reference_omegas <- function(mu, seed = 123L) {
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

default_jp_validation_scenarios <- function() {
  data.frame(
    kappa = c(0.5, 1, 2, 5),
    psi = c(-1, 0.5, -0.5, 1),
    stringsAsFactors = FALSE
  )
}

jp_projection_validation_plot_cdf <- function(cdf_df,
                                              kappa,
                                              psi) {
  ggplot(cdf_df, aes(x = t)) +
    geom_step(aes(y = empirical, color = "Empirical"), linewidth = 0.8, direction = "hv") +
    geom_line(aes(y = theoretical, color = "Theoretical"), linewidth = 0.9) +
    scale_color_manual(values = c(Empirical = "#252323", Theoretical = "red")) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1)) +
    labs(
      title = sprintf("JP projected CDF validation (kappa = %s, psi = %s)", kappa, psi),
      x = "t",
      y = "F_T(t)",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "bottom",
      plot.title = element_text(size = 12)
    )
}

jp_projection_validation_plot_density <- function(sample_projections,
                                                  density_df,
                                                  kappa,
                                                  psi) {
  sample_df <- data.frame(t = sample_projections)

  ggplot(sample_df, aes(x = t)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#6b6b9e", alpha = 0.35, color = "white") +
    geom_line(data = density_df, aes(x = t, y = theoretical_density), color = "red", linewidth = 0.9) +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(
      title = sprintf("JP projected density validation (kappa = %s, psi = %s)", kappa, psi),
      x = "t",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(size = 12)
    )
}

jp_projection_validation_summary <- function(sample,
                                             mu,
                                             kappa,
                                             psi,
                                             n_grid = 200L) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  projections <- as.numeric(sample %*% mu)
  axis_table <- build_jp_axis_cdf_table(mu = mu, kappa = kappa, psi = psi, grid_size = 8193L)
  t_grid <- seq(-1, 1, length.out = n_grid)
  empirical_cdf <- vapply(t_grid, function(t0) mean(projections <= t0), numeric(1))
  theoretical_cdf <- jp_interpolate_cdf(t_grid, x_grid = axis_table$u, cdf_grid = axis_table$cdf)
  theoretical_density <- d_proj_jp(t_grid, mu = mu, kappa = kappa, psi = psi)

  summary_df <- data.frame(
    kappa = as.numeric(kappa),
    psi = as.numeric(psi),
    n = nrow(sample),
    grid_size = as.integer(n_grid),
    max_abs_diff = max(abs(empirical_cdf - theoretical_cdf)),
    mean_abs_diff = mean(abs(empirical_cdf - theoretical_cdf)),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_df,
    cdf_df = data.frame(
      t = t_grid,
      empirical = empirical_cdf,
      theoretical = theoretical_cdf,
      stringsAsFactors = FALSE
    ),
    density_df = data.frame(
      t = t_grid,
      theoretical_density = theoretical_density,
      stringsAsFactors = FALSE
    ),
    projections = projections
  )
}

jp_distance_profile_validation_summary <- function(sample,
                                                   omega_values,
                                                   mu,
                                                   kappa,
                                                   psi,
                                                   t_grid) {
  omega_names <- names(omega_values)

  summary_rows <- lapply(seq_along(omega_values), function(i) {
    omega <- omega_values[[i]]
    empirical <- empirical_distance_profile(
      distances = compute_geodesic_distances_jp(omega, sample),
      t_values = t_grid
    )
    theoretical <- distance_profile_jp(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      kappa = kappa,
      psi = psi,
      distance_type = "geodesic"
    )

    data.frame(
      kappa = as.numeric(kappa),
      psi = as.numeric(psi),
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

run_jp_distance_profile_analysis <- function(kappa,
                                             psi,
                                             output_dir = NULL,
                                             mu = c(0, 0, 1),
                                             sample_sizes = c(50, 200),
                                             n_simulations = 10,
                                             validation_n = 5000L,
                                             n_points = 200L,
                                             seed = 123L,
                                             save_plots = TRUE) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  scenario_tag <- jp_scenario_tag(kappa = kappa, psi = psi)
  if (is.null(output_dir)) {
    output_dir <- file.path("output", "jp_geodesic", scenario_tag)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  omega_values <- jp_reference_omegas(mu = mu, seed = seed + 1L)
  omega_list <- unname(omega_values)
  legend_positions <- c("bottom_right", "top_left", "top_left", "top_left")

  data_generator <- function(n) {
    r_sph_jp(n = n, mu = mu, kappa = kappa, psi = psi)
  }
  theoretical_profile <- function(omega, t_values) {
    distance_profile_jp(
      omega = omega,
      t_values = t_values,
      mu = mu,
      kappa = kappa,
      psi = psi,
      distance_type = "geodesic"
    )
  }

  plots <- create_distance_profile_analysis(
    omega_values = omega_list,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = compute_geodesic_distances_jp,
    sample_sizes = sample_sizes,
    n_simulations = n_simulations,
    t_max = pi - 1e-8,
    save_plots = save_plots,
    output_dir = output_dir,
    file_prefix = "jp_geodesic_dp",
    legend_positions = legend_positions
  )

  set.seed(seed + 2L)
  validation_sample <- r_sph_jp(
    n = validation_n,
    mu = mu,
    kappa = kappa,
    psi = psi
  )

  t_grid <- seq(0, pi, length.out = n_points)
  distance_profile_summary <- jp_distance_profile_validation_summary(
    sample = validation_sample,
    omega_values = omega_values,
    mu = mu,
    kappa = kappa,
    psi = psi,
    t_grid = t_grid
  )
  utils::write.csv(
    distance_profile_summary,
    file = file.path(output_dir, "distance_profile_validation_summary.csv"),
    row.names = FALSE
  )

  projection_validation <- jp_projection_validation_summary(
    sample = validation_sample,
    mu = mu,
    kappa = kappa,
    psi = psi,
    n_grid = n_points
  )
  utils::write.csv(
    projection_validation$summary,
    file = file.path(output_dir, "projection_validation_summary.csv"),
    row.names = FALSE
  )

  cdf_plot <- jp_projection_validation_plot_cdf(
    cdf_df = projection_validation$cdf_df,
    kappa = kappa,
    psi = psi
  )
  density_plot <- jp_projection_validation_plot_density(
    sample_projections = projection_validation$projections,
    density_df = projection_validation$density_df,
    kappa = kappa,
    psi = psi
  )

  if (isTRUE(save_plots)) {
    ggsave(
      filename = file.path(output_dir, "projection_cdf_validation.png"),
      plot = cdf_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
    ggsave(
      filename = file.path(output_dir, "projection_density_validation.png"),
      plot = density_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  metadata_lines <- c(
    sprintf("kappa: %s", kappa),
    sprintf("psi: %s", psi),
    sprintf("output_dir: %s", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
    sprintf("sample_sizes: %s", paste(sample_sizes, collapse = ", ")),
    sprintf("n_simulations: %d", as.integer(n_simulations)),
    sprintf("validation_n: %d", as.integer(validation_n)),
    sprintf("n_points: %d", as.integer(n_points)),
    sprintf("seed: %d", as.integer(seed))
  )
  writeLines(metadata_lines, con = file.path(output_dir, "metadata.txt"))

  list(
    scenario_tag = scenario_tag,
    output_dir = output_dir,
    plots = plots,
    distance_profile_summary = distance_profile_summary,
    projection_validation_summary = projection_validation$summary,
    projection_cdf_plot = cdf_plot,
    projection_density_plot = density_plot
  )
}

run_jp_distance_profile_suite <- function(output_root = file.path("output", "jp_geodesic"),
                                          scenarios = default_jp_validation_scenarios(),
                                          mu = c(0, 0, 1),
                                          sample_sizes = c(50, 200),
                                          n_simulations = 10,
                                          validation_n = 5000L,
                                          n_points = 200L,
                                          seed = 123L,
                                          save_plots = TRUE) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  scenario_results <- vector("list", nrow(scenarios))

  for (i in seq_len(nrow(scenarios))) {
    scenario_kappa <- scenarios$kappa[[i]]
    scenario_psi <- scenarios$psi[[i]]
    scenario_dir <- file.path(output_root, jp_scenario_tag(scenario_kappa, scenario_psi))
    scenario_seed <- as.integer(seed + 100L * i)

    scenario_results[[i]] <- run_jp_distance_profile_analysis(
      kappa = scenario_kappa,
      psi = scenario_psi,
      output_dir = scenario_dir,
      mu = mu,
      sample_sizes = sample_sizes,
      n_simulations = n_simulations,
      validation_n = validation_n,
      n_points = n_points,
      seed = scenario_seed,
      save_plots = save_plots
    )
  }

  distance_profile_summary <- do.call(
    rbind,
    lapply(scenario_results, `[[`, "distance_profile_summary")
  )
  projection_summary <- do.call(
    rbind,
    lapply(scenario_results, `[[`, "projection_validation_summary")
  )

  utils::write.csv(
    distance_profile_summary,
    file = file.path(output_root, "suite_distance_profile_validation_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    projection_summary,
    file = file.path(output_root, "suite_projection_validation_summary.csv"),
    row.names = FALSE
  )

  list(
    output_root = output_root,
    scenarios = scenarios,
    results = scenario_results,
    distance_profile_summary = distance_profile_summary,
    projection_summary = projection_summary
  )
}
