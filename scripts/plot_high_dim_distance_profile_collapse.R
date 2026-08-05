#!/usr/bin/env Rscript

# Exact distance-profile plots along increasing ambient dimension.
#
# The normal block uses N_d(3 e_1, I_d) and centres omega of norm one.  Thus,
# their direction is encoded by rho = <omega, e_1>.  The vMF block uses
# vMF_d(e_1, 5) on S^{d-1}, with the geodesic distance and the same values of
# rho.  The vMF profiles are evaluated by the one-dimensional integral in
# theoretical_distance_profile_vmf(), not by Monte Carlo.

resolve_high_dim_dp_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  existing <- candidates[file.exists(candidates) | dir.exists(candidates)]
  if (!length(existing)) {
    stop(sprintf("Could not resolve path: %s", file.path(...)))
  }
  existing[[1L]]
}

parse_arguments <- function(args) {
  values <- list(
    n_cores = 1L,
    output_dir = file.path("output", "high_dim_distance_profiles"),
    dimensions = c(3L, 5L, 10L, 20L, 50L, 100L, 200L),
    n_t = 501L
  )

  for (argument in args) {
    if (startsWith(argument, "--n-cores=")) {
      values$n_cores <- as.integer(sub("^--n-cores=", "", argument))
    } else if (startsWith(argument, "--output-dir=")) {
      values$output_dir <- sub("^--output-dir=", "", argument)
    } else if (startsWith(argument, "--dimensions=")) {
      values$dimensions <- as.integer(strsplit(sub("^--dimensions=", "", argument), ",", fixed = TRUE)[[1L]])
    } else if (startsWith(argument, "--n-t=")) {
      values$n_t <- as.integer(sub("^--n-t=", "", argument))
    } else {
      stop(sprintf("Unknown argument: %s", argument))
    }
  }

  values$dimensions <- sort(unique(values$dimensions))
  if (any(!is.finite(values$dimensions)) || any(values$dimensions < 3L)) {
    stop("`--dimensions` must contain integers at least 3.")
  }
  if (length(values$n_cores) != 1L || !is.finite(values$n_cores) || values$n_cores < 1L) {
    stop("`--n-cores` must be a positive integer.")
  }
  if (length(values$n_t) != 1L || !is.finite(values$n_t) || values$n_t < 11L) {
    stop("`--n-t` must be an integer at least 11.")
  }

  values
}

unit_direction <- function(dimension, rho) {
  c(rho, sqrt(pmax(0, 1 - rho^2)), rep.int(0, dimension - 2L))
}

normal_profile <- function(t_values, dimension, rho, mu_norm = 3, omega_norm = 1) {
  noncentrality <- mu_norm^2 + omega_norm^2 - 2 * mu_norm * omega_norm * rho
  stats::pchisq(t_values^2, df = dimension, ncp = noncentrality)
}

vmf_profile <- function(t_values, dimension, rho, kappa = 5) {
  mu <- c(1, rep.int(0, dimension - 1L))
  omega <- unit_direction(dimension, rho)
  upper_bound <- pi
  profile <- numeric(length(t_values))
  below <- t_values <= 0
  above <- t_values >= upper_bound
  active <- !(below | above)

  profile[below] <- 0
  profile[above] <- 1
  if (any(active)) {
    profile[active] <- theoretical_distance_profile_vmf(
      omega = omega,
      mu = mu,
      kappa = kappa,
      t_values = t_values[active],
      distance_type = "geodesic"
    )
  }
  pmin(1, pmax(0, profile))
}

rho_label <- function(rho) {
  sprintf("rho_%s", gsub("-", "minus_", formatC(rho, format = "f", digits = 1L)))
}

rho_plot_label <- function(rho) {
  sprintf("rho = %s", formatC(rho, format = "f", digits = 1L))
}

make_profile_data <- function(model, dimensions, rhos, t_values, n_cores) {
  jobs <- expand.grid(
    dimension = dimensions,
    rho = rhos,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  evaluate_job <- function(index) {
    dimension <- jobs$dimension[[index]]
    rho <- jobs$rho[[index]]
    values <- if (identical(model, "normal")) {
      normal_profile(t_values, dimension = dimension, rho = rho)
    } else {
      vmf_profile(t_values, dimension = dimension, rho = rho)
    }
    data.frame(
      model = model,
      dimension = dimension,
      rho = rho,
      t = t_values,
      profile = values
    )
  }

  indices <- seq_len(nrow(jobs))
  profiles <- if (.Platform$OS.type == "windows" || n_cores == 1L) {
    lapply(indices, evaluate_job)
  } else {
    parallel::mclapply(indices, evaluate_job, mc.cores = n_cores, mc.preschedule = FALSE)
  }
  do.call(rbind, profiles)
}

plot_profiles_by_rho <- function(profile_data,
                                 model,
                                 output_dir,
                                 palette,
                                 scale_normal_radius = FALSE) {
  if (scale_normal_radius && !identical(model, "normal")) {
    stop("Only the Euclidean normal profiles can be plotted on the t/sqrt(d) scale.")
  }
  model_data <- profile_data[profile_data$model == model, , drop = FALSE]
  model_data$dimension_label <- factor(model_data$dimension, levels = sort(unique(model_data$dimension)))
  model_data$x_value <- if (scale_normal_radius) model_data$t / sqrt(model_data$dimension) else model_data$t
  rhos <- sort(unique(model_data$rho), decreasing = TRUE)
  subtitle <- if (identical(model, "normal")) {
    bquote(X ~ "\u223c" ~ N[d](3 * e[1], I[d]) ~ "," ~ norm(omega) == 1)
  } else {
    bquote(X ~ "\u223c" ~ plain(vMF)[d](e[1], kappa == 5) ~ ", geodesic distance")
  }
  x_label <- if (scale_normal_radius) expression(t / sqrt(d)) else "t"
  scale_suffix <- if (scale_normal_radius) " on the t/sqrt(d) scale" else ""
  file_prefix <- if (scale_normal_radius) {
    sprintf("%s_distance_profiles_scaled", model)
  } else {
    sprintf("%s_distance_profiles", model)
  }

  for (rho in rhos) {
    plot_data <- model_data[abs(model_data$rho - rho) < 1e-12, , drop = FALSE]
    figure <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = x_value, y = profile, colour = dimension_label, group = dimension_label)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_colour_manual(values = palette, name = "ambient dimension d") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        x = x_label,
        y = expression(F[omega](t)),
        title = sprintf("%s distance profiles%s (%s)", model, scale_suffix, rho_plot_label(rho)),
        subtitle = subtitle
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(legend.position = "right")

    file_name <- sprintf("%s_%s.png", file_prefix, rho_label(rho))
    ggplot2::ggsave(
      filename = file.path(output_dir, file_name),
      plot = figure,
      width = 8.5,
      height = 5.5,
      dpi = 300,
      bg = "white"
    )
  }
}

main <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required.")
  }
  if (!requireNamespace("rotasym", quietly = TRUE)) {
    stop("Package `rotasym` is required.")
  }

  config <- parse_arguments(commandArgs(trailingOnly = TRUE))
  source(resolve_high_dim_dp_path("utils.R"))

  rhos <- c(-1, -0.5, 0, 0.5, 1)
  normal_mu_norm <- 3
  normal_omega_norm <- 1
  normal_max_ncp <- (normal_mu_norm + normal_omega_norm)^2
  normal_t <- seq(
    0,
    sqrt(stats::qchisq(0.999999, df = max(config$dimensions), ncp = normal_max_ncp)),
    length.out = config$n_t
  )
  vmf_t <- seq(0, pi, length.out = config$n_t)

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  normal_data <- make_profile_data(
    model = "normal",
    dimensions = config$dimensions,
    rhos = rhos,
    t_values = normal_t,
    n_cores = config$n_cores
  )
  vmf_data <- make_profile_data(
    model = "vMF",
    dimensions = config$dimensions,
    rhos = rhos,
    t_values = vmf_t,
    n_cores = config$n_cores
  )
  all_data <- rbind(normal_data, vmf_data)

  palette <- grDevices::hcl.colors(length(config$dimensions), palette = "Dynamic")
  names(palette) <- as.character(config$dimensions)
  plot_profiles_by_rho(
    normal_data,
    model = "normal",
    output_dir = config$output_dir,
    palette = palette
  )
  plot_profiles_by_rho(
    normal_data,
    model = "normal",
    output_dir = config$output_dir,
    palette = palette,
    scale_normal_radius = TRUE
  )
  plot_profiles_by_rho(
    vmf_data,
    model = "vMF",
    output_dir = config$output_dir,
    palette = palette
  )

  utils::write.csv(all_data, file.path(config$output_dir, "distance_profile_values.csv"), row.names = FALSE)
  writeLines(
    c(
      sprintf("n_cores: %d", config$n_cores),
      sprintf("dimensions: %s", paste(config$dimensions, collapse = ", ")),
      sprintf("n_t: %d", config$n_t),
      "normal: X ~ N_d(3 e_1, I_d); ||omega|| = 1; rho = <omega, e_1>.",
      "vMF: X ~ vMF_d(e_1, kappa = 5); geodesic distance; rho = <omega, e_1>.",
      "The additional normal plots on the t/sqrt(d) scale show the relative-distance collapse.",
      "The profiles are theoretical. The normal formula is a noncentral chi-square CDF; the vMF formula is evaluated by one-dimensional numerical integration."
    ),
    con = file.path(config$output_dir, "metadata.txt")
  )

  message(sprintf("Saved 15 plots and the profile values to %s", normalizePath(config$output_dir)))
}

main()
