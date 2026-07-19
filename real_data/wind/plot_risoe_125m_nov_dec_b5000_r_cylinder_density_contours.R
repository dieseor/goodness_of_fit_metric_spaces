#!/usr/bin/env Rscript

# Cylinder density plots in the original r = speed / mean(speed) coordinate,
# matched exactly to the strict-step-four November-December B=5000 analysis.
# Both densities are with respect to dtheta dr.  Hence the fitted HvMF density
# is multiplied by sinh(r), while the reflected directional-linear KDE already
# has dtheta dr as its dominating product measure.

if (!exists("hvmf_log_cylinder_density", mode = "function")) {
  source(file.path(
    "real_data", "wind", "plot_risoe_125m_nov_dec_b5000_log_cylinder_density_contours.R"
  ))
}

log_reflected_gaussian_kernel <- function(r, center, bandwidth) {
  a <- stats::dnorm((r - center) / bandwidth, log = TRUE)
  b <- stats::dnorm((r + center) / bandwidth, log = TRUE)
  pmax(a, b) + log1p(exp(-abs(a - b))) - log(bandwidth)
}

select_r_bandwidths <- function(r, theta) {
  n <- length(r)
  radial_scale <- min(stats::sd(r), stats::IQR(r) / 1.349)
  if (!is.finite(radial_scale) || radial_scale <= 0) radial_scale <- stats::sd(r)
  if (!is.finite(radial_scale) || radial_scale <= 0) radial_scale <- 0.25
  h_initial <- max(0.04, 1.06 * radial_scale * n^(-1 / 6))

  resultant <- Mod(mean(exp(1i * theta)))
  kappa_initial <- if (resultant < 0.53) {
    2 * resultant + resultant^3 + 5 * resultant^5 / 6
  } else if (resultant < 0.85) {
    -0.4 + 1.39 * resultant + 0.43 / (1 - resultant)
  } else {
    1 / (resultant^3 - 4 * resultant^2 + 3 * resultant)
  }
  kappa_initial <- min(max(kappa_initial, 0.15), 50)

  objective <- function(par) {
    h <- exp(par[[1L]])
    kappa <- exp(par[[2L]])
    total <- 0
    for (i in seq_len(n)) {
      other <- seq_len(n) != i
      terms <- log_von_mises_kernel(theta[[i]] - theta[other], kappa) +
        log_reflected_gaussian_kernel(r[[i]], r[other], h)
      total <- total + log_sum_exp(terms) - log(n - 1)
    }
    -total
  }

  fit <- stats::optim(
    log(c(h_initial, kappa_initial)), objective, method = "L-BFGS-B",
    lower = log(c(max(0.015, radial_scale / 20), 0.05)),
    upper = log(c(max(0.25, 2.5 * radial_scale), 100))
  )
  if (fit$convergence != 0L || any(!is.finite(fit$par))) {
    warning("Reflected-r bandwidth optimization failed; using initial values.")
    return(c(h_r = h_initial, kappa_theta = kappa_initial))
  }
  c(h_r = exp(fit$par[[1L]]), kappa_theta = exp(fit$par[[2L]]))
}

hvmf_r_cylinder_density <- function(r, theta, fit) {
  hvmf_hyperbolic_density(r, theta, fit) * sinh(r)
}

kde_r_cylinder_density <- function(r_grid, theta_grid, r, theta, bandwidths) {
  density <- numeric(length(r_grid))
  for (j in seq_along(r)) {
    density <- density + exp(
      log_von_mises_kernel(theta_grid - theta[[j]], bandwidths[["kappa_theta"]]) +
        log_reflected_gaussian_kernel(r_grid, r[[j]], bandwidths[["h_r"]])
    )
  }
  density / length(r)
}

reflected_kernel_mass_error <- function(r, bandwidths) {
  centers <- unique(stats::quantile(r, c(0, 0.5, 1), names = FALSE))
  radial_errors <- vapply(centers, function(center) {
    mass <- stats::integrate(
      function(x) exp(log_reflected_gaussian_kernel(x, center, bandwidths[["h_r"]])),
      0, Inf, rel.tol = 1e-11
    )$value
    abs(mass - 1)
  }, numeric(1))
  circular_mass <- stats::integrate(
    function(delta) exp(log_von_mises_kernel(delta, bandwidths[["kappa_theta"]])),
    -pi, pi, rel.tol = 1e-11
  )$value
  max(c(radial_errors, abs(circular_mass - 1)))
}

make_r_density_grid <- function(case, n_theta = 241L, n_r = 181L, r_max) {
  theta_values <- seq(0, 2 * pi, length.out = n_theta)
  r_values <- seq(0, r_max, length.out = n_r)
  grid <- expand.grid(theta = theta_values, r = r_values)
  fit <- hvmf_mle_h2(as.matrix(case[, c("x0", "x1", "x2")]))
  bandwidths <- select_r_bandwidths(case$r, case$theta)
  kernel_mass_error <- reflected_kernel_mass_error(case$r, bandwidths)
  if (!is.finite(kernel_mass_error) || kernel_mass_error > 1e-9) {
    stop("Reflected product-kernel normalization failed.", call. = FALSE)
  }

  parametric <- hvmf_r_cylinder_density(grid$r, grid$theta, fit)
  nonparametric <- kde_r_cylinder_density(
    grid$r, grid$theta, case$r, case$theta, bandwidths
  )
  parametric_matrix <- t(matrix(parametric, nrow = n_theta, ncol = n_r))
  nonparametric_matrix <- t(matrix(nonparametric, nrow = n_theta, ncol = n_r))
  seam_error <- max(
    abs(parametric_matrix[, 1L] - parametric_matrix[, n_theta]),
    abs(nonparametric_matrix[, 1L] - nonparametric_matrix[, n_theta])
  )
  if (seam_error > 1e-10) stop("Cylinder seam validation failed.", call. = FALSE)

  grid$theta_deg <- grid$theta * 180 / pi
  grid$window <- "nov_dec"
  grid$window_label <- "November + December"
  grid_parametric <- transform(
    grid, estimator = "Fitted HvMF", density = parametric,
    hdr_content = density_to_hdr_content(parametric)
  )
  grid_nonparametric <- transform(
    grid, estimator = "Nonparametric KDE", density = nonparametric,
    hdr_content = density_to_hdr_content(nonparametric)
  )
  list(
    grid = rbind(grid_parametric, grid_nonparametric),
    fit = fit,
    bandwidths = bandwidths,
    full_parametric_mass = full_hvmf_log_cylinder_mass(fit),
    full_nonparametric_mass = 1,
    parametric_mass_on_plot = trapezoid_2d_mass(parametric_matrix, theta_values, r_values),
    nonparametric_mass_on_plot = trapezoid_2d_mass(nonparametric_matrix, theta_values, r_values),
    kernel_mass_error = kernel_mass_error,
    seam_error = seam_error
  )
}

plot_r_pattern <- function(case, pattern_id, output_dir) {
  r_max <- max(3.5, 1.08 * max(case$r))
  object <- make_r_density_grid(case, r_max = r_max)
  grid_data <- object$grid
  grid_data$estimator <- factor(
    grid_data$estimator, levels = c("Fitted HvMF", "Nonparametric KDE")
  )
  points <- case
  palette <- grDevices::colorRampPalette(
    c("#FFFFFF", "#E8F3F8", "#B9DDE7", "#72B7C5", "#2A788E", "#234B70")
  )(256L)

  p <- ggplot2::ggplot(grid_data, ggplot2::aes(theta_deg, r)) +
    ggplot2::geom_raster(ggplot2::aes(fill = density), interpolate = FALSE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = hdr_content), breaks = c(0.25, 0.50, 0.75, 0.90),
      colour = "#17324D", linewidth = 0.28, alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = points, ggplot2::aes(theta_deg, r), inherit.aes = FALSE,
      shape = 21, size = 1.25, stroke = 0.28,
      colour = "#5A1A1A", fill = "#FFF7F2", alpha = 0.88
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(estimator)) +
    ggplot2::scale_fill_gradientn(
      colours = palette, limits = c(0, max(grid_data$density)),
      name = expression("Density w.r.t. " * dr * dtheta)
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 90, 180, 270, 360),
      labels = c("0°", "90°", "180°", "270°", "360°"), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, 360), ylim = c(0, r_max), expand = FALSE) +
    ggplot2::labs(
      x = expression("Wind direction " * theta * " (0° and 360° identified)"),
      y = expression("Scaled speed " * r == s / bar(s)),
      title = sprintf("Risø %d m wind data on the unwrapped cylinder", plot_height_m),
      subtitle = sprintf(
        "%s: days %s; November-December 1996-2003%s; HDR contours: 25%%, 50%%, 75%%, 90%%",
        pattern_id, paste(day_patterns[[pattern_id]], collapse = ", "), plot_bandwidth_note
      )
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.8, "lines"),
      strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey35"),
      strip.text = ggplot2::element_text(face = "bold"),
      aspect.ratio = 1,
      legend.position = "right",
      plot.title.position = "plot"
    )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(
    output_dir, sprintf("risoe_%dm_r_cylinder_density_contours_%s", plot_height_m, pattern_id)
  )
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 8.3, height = 4.3, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 8.3, height = 4.3, units = "in", dpi = 320)

  summary <- data.frame(
    window = "nov_dec", pattern = pattern_id, n = nrow(case),
    speed_mean = case$speed_mean[[1L]], kappa_hat = object$fit$kappa,
    mu_r_hat = object$fit$chi, mu_theta_deg_hat = object$fit$theta_deg,
    kde_h_r = object$bandwidths[["h_r"]],
    kde_kappa_theta = object$bandwidths[["kappa_theta"]],
    r_plot_max = r_max, full_parametric_mass = object$full_parametric_mass,
    full_nonparametric_mass = object$full_nonparametric_mass,
    parametric_mass_on_plot = object$parametric_mass_on_plot,
    nonparametric_mass_on_plot = object$nonparametric_mass_on_plot,
    kernel_mass_error = object$kernel_mass_error, seam_error = object$seam_error
  )
  utils::write.csv(summary, paste0(stem, "_summary.csv"), row.names = FALSE)
  summary
}

run_nov_dec_b5000_r_density_plots <- function(
    patterns = names(day_patterns),
    output_dir = file.path(results_dir, "r_cylinder_density_contours")) {
  validate_plot_samples()
  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  rows <- lapply(patterns, function(pattern_id) {
    case <- build_case(selected, "nov_dec", pattern_id)
    plot_r_pattern(case, pattern_id, output_dir)
  })
  summary <- do.call(rbind, rows)
  utils::write.csv(
    summary, file.path(
      output_dir, sprintf("risoe_%dm_r_cylinder_density_contours_summary.csv", plot_height_m)
    ),
    row.names = FALSE
  )
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  requested_patterns <- if (length(args) == 0L) names(day_patterns) else args
  result <- run_nov_dec_b5000_r_density_plots(patterns = requested_patterns)
  print(result, row.names = FALSE, digits = 5)
}
