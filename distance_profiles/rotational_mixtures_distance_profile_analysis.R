suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(gridExtra)
})

source(file.path("distance_profiles", "distance_profile_analysis.R"))
source(file.path("utils.R"))

format_rotmix_value_tag <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) {
    stop("`x` must be a finite scalar.")
  }

  tag <- format(x, scientific = FALSE, trim = TRUE)
  tag <- gsub("-", "neg", tag, fixed = TRUE)
  gsub("\\.", "p", tag)
}

rotmix_reference_omegas <- function(mu, seed = 123L) {
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

compute_rotmix_distances <- function(omega, data, distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1L)
  }
  apply(data, 1, function(x) sphere_distance(omega, x, distance_type = distance_type))
}

rotmix_default_theta <- function(model_name = c("rotational_beta_mixture2", "rotational_logitnormal_mixture2"),
                                 mu = c(0, 0, 1)) {
  model_name <- match.arg(model_name)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)

  if (identical(model_name, "rotational_beta_mixture2")) {
    return(list(
      mu = mu,
      weight1 = 0.4,
      alpha1 = 2,
      beta1 = 8,
      alpha2 = 8,
      beta2 = 2
    ))
  }

  list(
    mu = mu,
    weight1 = 0.45,
    mean1 = -1.2,
    sd1 = 0.5,
    mean2 = 1.0,
    sd2 = 0.6
  )
}

rotmix_model_tag <- function(model_name, theta) {
  if (identical(model_name, "rotational_beta_mixture2")) {
    return(sprintf(
      "w%s_a1_%s_b1_%s_a2_%s_b2_%s",
      format_rotmix_value_tag(theta$weight1),
      format_rotmix_value_tag(theta$alpha1),
      format_rotmix_value_tag(theta$beta1),
      format_rotmix_value_tag(theta$alpha2),
      format_rotmix_value_tag(theta$beta2)
    ))
  }

  sprintf(
    "w%s_m1_%s_sd1_%s_m2_%s_sd2_%s",
    format_rotmix_value_tag(theta$weight1),
    format_rotmix_value_tag(theta$mean1),
    format_rotmix_value_tag(theta$sd1),
    format_rotmix_value_tag(theta$mean2),
    format_rotmix_value_tag(theta$sd2)
  )
}

rotmix_sample_generator <- function(model_name, theta) {
  if (identical(model_name, "rotational_beta_mixture2")) {
    return(function(n) {
      r_sph_rotational_beta_mixture2(
        n = n,
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2
      )
    })
  }

  function(n) {
    r_sph_rotational_logitnormal_mixture2(
      n = n,
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2
    )
  }
}

rotmix_theoretical_profile <- function(model_name, theta, distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)

  if (identical(model_name, "rotational_beta_mixture2")) {
    return(function(omega, t_values) {
      distance_profile_rotational_beta_mixture2(
        omega = omega,
        t_values = t_values,
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
        distance_type = distance_type,
        method = "legendre"
      )
    })
  }

  function(omega, t_values) {
    distance_profile_rotational_logitnormal_mixture2(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      distance_type = distance_type,
      method = "legendre"
    )
  }
}

rotmix_distance_profile_validation_summary <- function(sample,
                                                       omega_values,
                                                       model_name,
                                                       theta,
                                                       t_grid,
                                                       distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_names <- names(omega_values)
  theoretical_profile <- rotmix_theoretical_profile(
    model_name = model_name,
    theta = theta,
    distance_type = distance_type
  )

  summary_rows <- lapply(seq_along(omega_values), function(i) {
    omega <- omega_values[[i]]
    empirical <- empirical_distance_profile(
      distances = compute_rotmix_distances(omega, sample, distance_type = distance_type),
      t_values = t_grid
    )
    theoretical <- theoretical_profile(omega, t_grid)

    data.frame(
      model = model_name,
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

run_rotational_mixture_distance_profile_analysis <- function(model_name = c("rotational_beta_mixture2", "rotational_logitnormal_mixture2"),
                                                             theta = NULL,
                                                             distance_type = c("geodesic", "chordal"),
                                                             output_dir = NULL,
                                                             sample_sizes = c(50, 200),
                                                             n_simulations = 10,
                                                             validation_n = 5000L,
                                                             n_points = 200L,
                                                             seed = 123L,
                                                             save_plots = TRUE) {
  model_name <- match.arg(model_name)
  distance_type <- match.arg(distance_type)
  theta <- theta %||% rotmix_default_theta(model_name = model_name)
  theta$mu <- jp_normalize_unit_vector(theta$mu, arg_name = "`theta$mu`", min_length = 3L)

  scenario_tag <- rotmix_model_tag(model_name = model_name, theta = theta)
  if (is.null(output_dir)) {
    output_dir <- file.path("output", paste0(model_name, "_", distance_type), scenario_tag)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  omega_values <- rotmix_reference_omegas(mu = theta$mu, seed = seed + 1L)
  omega_list <- unname(omega_values)
  legend_positions <- c("bottom_right", "top_left", "top_left", "top_left")
  data_generator <- rotmix_sample_generator(model_name = model_name, theta = theta)
  theoretical_profile <- rotmix_theoretical_profile(
    model_name = model_name,
    theta = theta,
    distance_type = distance_type
  )

  t_max <- if (identical(distance_type, "geodesic")) pi - 1e-8 else 2 - 1e-8

  plots <- create_distance_profile_analysis(
    omega_values = omega_list,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = function(omega, data) compute_rotmix_distances(omega, data, distance_type = distance_type),
    sample_sizes = sample_sizes,
    n_simulations = n_simulations,
    t_max = t_max,
    save_plots = save_plots,
    output_dir = output_dir,
    file_prefix = paste0(model_name, "_", distance_type, "_dp"),
    legend_positions = legend_positions
  )

  set.seed(seed + 2L)
  validation_sample <- data_generator(validation_n)
  t_grid <- seq(0, t_max, length.out = n_points)
  summary_df <- rotmix_distance_profile_validation_summary(
    sample = validation_sample,
    omega_values = omega_values,
    model_name = model_name,
    theta = theta,
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

run_rotational_mixture_distance_profile_suite <- function(distance_types = c("geodesic"),
                                                          seed = 123L,
                                                          sample_sizes = c(50, 200),
                                                          n_simulations = 10,
                                                          validation_n = 5000L,
                                                          n_points = 200L,
                                                          save_plots = TRUE) {
  models <- c("rotational_beta_mixture2", "rotational_logitnormal_mixture2")
  results <- list()
  counter <- 1L

  for (model_name in models) {
    for (distance_type in distance_types) {
      results[[counter]] <- run_rotational_mixture_distance_profile_analysis(
        model_name = model_name,
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
  run_rotational_mixture_distance_profile_suite()
}
