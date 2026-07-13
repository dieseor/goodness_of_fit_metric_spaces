#!/usr/bin/env Rscript

# Log-cylinder-coordinate density plots for the Risoe 125 m wind samples.
#
# The hyperbolic polar chart is
#   x(r, theta) = (cosh(r), sinh(r) cos(theta), sinh(r) sin(theta)),
# with r = speed / mean(speed), followed by z = log(r).  Since the Riemannian
# area element on H^2 is sinh(r) dr dtheta and dr = exp(z) dz, every displayed
# density is with respect to dz dtheta on R x S^1.  Thus, the fitted HvMF
# density includes the full Jacobian sinh(exp(z)) exp(z).  The nonparametric
# estimate is a directional-linear product KDE: a von Mises kernel in theta
# and an ordinary Gaussian kernel in z.

script_path <- function() {
  file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_args) == 1L) {
    return(normalizePath(sub("^--file=", "", file_args), mustWork = TRUE))
  }
  normalizePath(file.path("real_data", "wind", "plot_risoe_cylinder_density_contours.R"),
                mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(script_path()), "..", ".."), mustWork = TRUE)
source(file.path(repo_root, "real_data", "wind", "preprocess_risoe_modern_hvmf.R"))
source(file.path(repo_root, "utils.R"))

required_packages <- c("ggplot2", "ncdf4")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

select_noon_all_months <- function(df, fixed_tz = "UTC") {
  date <- as.Date(df$datetime, tz = fixed_tz)
  noon <- as.POSIXct(paste(date, "12:00:00"), tz = fixed_tz)
  distance_to_noon <- abs(as.numeric(difftime(df$datetime, noon, units = "secs")))
  ord <- order(date, distance_to_noon, as.numeric(df$datetime))
  ordered <- df[ord, , drop = FALSE]
  selected <- ordered[!duplicated(as.Date(ordered$datetime, tz = fixed_tz)), , drop = FALSE]
  selected[order(selected$datetime), , drop = FALSE]
}

month_windows <- list(
  nov_dec = c(11L, 12L),
  feb = 2L,
  jan_feb = c(1L, 2L)
)

window_labels <- c(
  nov_dec = "November + December",
  feb = "February",
  jan_feb = "January + February"
)

# This reproduces the sampling convention used by the current screening code.
# In particular, start3 includes day 30 to match the stored Jensen-style set12.
day_patterns <- list(
  start1 = c(1L, 5L, 9L, 13L, 17L, 21L, 25L, 29L),
  start2 = c(2L, 6L, 10L, 14L, 18L, 22L, 26L, 30L),
  start3 = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
  start4 = c(4L, 8L, 12L, 16L, 20L, 24L, 28L)
)

build_case <- function(selected, window_id, pattern_id, years = 1996:2003) {
  months <- month_windows[[window_id]]
  days <- day_patterns[[pattern_id]]
  keep <- selected$year %in% years & selected$month %in% months & selected$day %in% days
  rows <- selected[keep, , drop = FALSE]
  rows <- rows[as.Date(rows$datetime, tz = "UTC") < as.Date("2004-12-01"), , drop = FALSE]
  rows <- rows[is.finite(rows$ws125) & is.finite(rows$wd125) & rows$ws125 > 0, , drop = FALSE]
  rownames(rows) <- NULL

  if (nrow(rows) < 2L) {
    stop("Too few valid observations for ", window_id, "/", pattern_id, call. = FALSE)
  }
  if (anyDuplicated(as.Date(rows$datetime, tz = "UTC"))) {
    stop("More than one observation was selected on the same date.", call. = FALSE)
  }

  speed <- as.numeric(rows$ws125)
  theta <- ((as.numeric(rows$wd125) %% 360) + 360) %% 360 * pi / 180
  speed_mean <- mean(speed)
  r <- speed / speed_mean
  z <- log(r)
  out <- data.frame(
    window = window_id,
    pattern = pattern_id,
    datetime = rows$datetime,
    speed = speed,
    speed_mean = speed_mean,
    r = r,
    z = z,
    theta = theta,
    theta_deg = theta * 180 / pi,
    x0 = cosh(r),
    x1 = sinh(r) * cos(theta),
    x2 = sinh(r) * sin(theta)
  )
  norm_error <- max(abs(-out$x0^2 + out$x1^2 + out$x2^2 + 1))
  if (!is.finite(norm_error) || norm_error > 1e-10) {
    stop("Hyperboloid-coordinate validation failed.", call. = FALSE)
  }
  out
}

log_i0 <- function(kappa) {
  log(besselI(kappa, nu = 0, expon.scaled = TRUE)) + abs(kappa)
}

log_von_mises_kernel <- function(delta, kappa) {
  kappa * cos(delta) - log(2 * pi) - log_i0(kappa)
}

log_gaussian_kernel <- function(z, center, bandwidth) {
  stats::dnorm((z - center) / bandwidth, log = TRUE) - log(bandwidth)
}

log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

select_cylinder_bandwidths <- function(z, theta) {
  n <- length(z)
  linear_scale <- min(stats::sd(z), stats::IQR(z) / 1.349)
  if (!is.finite(linear_scale) || linear_scale <= 0) linear_scale <- stats::sd(z)
  if (!is.finite(linear_scale) || linear_scale <= 0) linear_scale <- 0.5
  g_initial <- max(0.04, 1.06 * linear_scale * n^(-1 / 6))

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
    g <- exp(par[[1L]])
    kappa <- exp(par[[2L]])
    total <- 0
    for (i in seq_len(n)) {
      other <- seq_len(n) != i
      log_terms <- log_von_mises_kernel(theta[[i]] - theta[other], kappa) +
        log_gaussian_kernel(z[[i]], z[other], g)
      total <- total + log_sum_exp(log_terms) - log(n - 1)
    }
    -total
  }

  fit <- stats::optim(
    par = log(c(g_initial, kappa_initial)),
    fn = objective,
    method = "L-BFGS-B",
    lower = log(c(max(0.015, linear_scale / 20), 0.05)),
    upper = log(c(max(0.25, 2.5 * linear_scale), 100))
  )
  if (fit$convergence != 0L || any(!is.finite(fit$par))) {
    warning("Bandwidth likelihood optimization did not converge; using initial values.")
    return(c(g_z = g_initial, kappa_theta = kappa_initial))
  }
  c(g_z = exp(fit$par[[1L]]), kappa_theta = exp(fit$par[[2L]]))
}

hvmf_hyperbolic_density <- function(r, theta, fit) {
  x0 <- cosh(r)
  x1 <- sinh(r) * cos(theta)
  x2 <- sinh(r) * sin(theta)
  inner <- -x0 * fit$mu[[1L]] + x1 * fit$mu[[2L]] + x2 * fit$mu[[3L]]
  normalizer <- fit$kappa * exp(fit$kappa) / (2 * pi)
  normalizer * exp(fit$kappa * inner)
}

hvmf_log_cylinder_density <- function(z, theta, fit) {
  r <- exp(z)
  hvmf_hyperbolic_density(r, theta, fit) * sinh(r) * r
}

kde_log_cylinder_density <- function(z_grid, theta_grid, z, theta, bandwidths) {
  n_grid <- length(z_grid)
  density <- numeric(n_grid)
  for (j in seq_along(z)) {
    density <- density +
      exp(log_von_mises_kernel(theta_grid - theta[[j]], bandwidths[["kappa_theta"]]) +
            log_gaussian_kernel(z_grid, z[[j]], bandwidths[["g_z"]]))
  }
  density / length(z)
}

product_kernel_mass_error <- function(bandwidths) {
  circular_mass <- stats::integrate(
    function(delta) exp(log_von_mises_kernel(delta, bandwidths[["kappa_theta"]])),
    lower = -pi, upper = pi, rel.tol = 1e-11
  )$value
  linear_mass <- stats::integrate(
    function(u) exp(log_gaussian_kernel(u, 0, bandwidths[["g_z"]])),
    lower = -Inf, upper = Inf, rel.tol = 1e-11
  )$value
  max(abs(c(circular_mass, linear_mass) - 1))
}

trapezoid_2d_mass <- function(density, theta_values, linear_values) {
  theta_step <- diff(theta_values)[[1L]]
  linear_step <- diff(linear_values)[[1L]]
  theta_weights <- rep(1, length(theta_values))
  theta_weights[c(1L, length(theta_weights))] <- 0.5
  linear_weights <- rep(1, length(linear_values))
  linear_weights[c(1L, length(linear_weights))] <- 0.5
  sum(density * outer(linear_weights, theta_weights)) * linear_step * theta_step
}

density_to_hdr_content <- function(density) {
  if (any(!is.finite(density)) || any(density < 0) || sum(density) <= 0) {
    stop("HDR conversion requires a finite, nonnegative density grid.", call. = FALSE)
  }
  ord <- order(density, decreasing = TRUE)
  cumulative_mass <- cumsum(density[ord]) / sum(density)
  content <- numeric(length(density))
  content[ord] <- cumulative_mass
  content
}

full_hvmf_log_cylinder_mass <- function(fit) {
  log_radial_density <- function(z) {
    r <- exp(z)
    bessel_argument <- fit$kappa * fit$sinh_chi * sinh(r)
    log(fit$kappa) + fit$kappa - fit$kappa * fit$mu[[1L]] * cosh(r) +
      log_i0(bessel_argument) + log(sinh(r)) + z
  }
  integrand <- function(z) {
    out <- numeric(length(z))
    finite_region <- is.finite(z) & z > -30 & z < log(30)
    out[finite_region] <- exp(log_radial_density(z[finite_region]))
    out
  }
  stats::integrate(integrand, lower = -30, upper = log(30),
                   rel.tol = 1e-10, subdivisions = 500L)$value
}

make_density_grid <- function(case, n_theta = 241L, n_z = 201L, z_limits) {
  theta_values <- seq(0, 2 * pi, length.out = n_theta)
  z_values <- seq(z_limits[[1L]], z_limits[[2L]], length.out = n_z)
  grid <- expand.grid(theta = theta_values, z = z_values)
  X <- as.matrix(case[, c("x0", "x1", "x2")])
  fit <- hvmf_mle_h2(X)
  bandwidths <- select_cylinder_bandwidths(case$z, case$theta)
  kernel_mass_error <- product_kernel_mass_error(bandwidths)
  if (!is.finite(kernel_mass_error) || kernel_mass_error > 1e-10) {
    stop("Directional-linear product-kernel normalization failed.", call. = FALSE)
  }

  parametric <- hvmf_log_cylinder_density(grid$z, grid$theta, fit)
  nonparametric <- kde_log_cylinder_density(
    grid$z, grid$theta, case$z, case$theta, bandwidths
  )
  # expand.grid() varies its first coordinate fastest. Transpose so rows are
  # log-speed positions and columns are angular positions.
  parametric_matrix <- t(matrix(parametric, nrow = n_theta, ncol = n_z))
  nonparametric_matrix <- t(matrix(nonparametric, nrow = n_theta, ncol = n_z))

  r_check <- exp(grid$z)
  transformed_check <- hvmf_hyperbolic_density(r_check, grid$theta, fit) *
    sinh(r_check) * r_check
  jacobian_error <- max(abs(parametric - transformed_check))
  if (!is.finite(jacobian_error) || jacobian_error > 1e-12) {
    stop("Log-cylinder Jacobian validation failed.", call. = FALSE)
  }

  full_parametric_mass <- full_hvmf_log_cylinder_mass(fit)
  if (!is.finite(full_parametric_mass) || abs(full_parametric_mass - 1) > 1e-8) {
    stop("The transformed fitted HvMF density does not integrate to one.", call. = FALSE)
  }

  seam_error_parametric <- max(abs(parametric_matrix[, 1L] - parametric_matrix[, n_theta]))
  seam_error_nonparametric <- max(abs(nonparametric_matrix[, 1L] - nonparametric_matrix[, n_theta]))
  if (max(seam_error_parametric, seam_error_nonparametric) > 1e-10) {
    stop("Cylinder seam validation failed.", call. = FALSE)
  }

  grid$theta_deg <- grid$theta * 180 / pi
  grid$r <- exp(grid$z)
  grid$window <- case$window[[1L]]
  grid$window_label <- unname(window_labels[[case$window[[1L]]]])
  grid_parametric <- transform(
    grid,
    estimator = "Fitted HvMF",
    density = parametric,
    hdr_content = density_to_hdr_content(parametric)
  )
  grid_nonparametric <- transform(
    grid,
    estimator = "Nonparametric KDE",
    density = nonparametric,
    hdr_content = density_to_hdr_content(nonparametric)
  )

  list(
    grid = rbind(grid_parametric, grid_nonparametric),
    fit = fit,
    bandwidths = bandwidths,
    full_parametric_mass = full_parametric_mass,
    full_nonparametric_mass = 1,
    parametric_mass_on_plot = trapezoid_2d_mass(parametric_matrix, theta_values, z_values),
    nonparametric_mass_on_plot = trapezoid_2d_mass(nonparametric_matrix, theta_values, z_values),
    jacobian_error = jacobian_error,
    kernel_mass_error = kernel_mass_error,
    seam_error_parametric = seam_error_parametric,
    seam_error_nonparametric = seam_error_nonparametric
  )
}

plot_pattern <- function(cases, pattern_id, output_dir) {
  all_z <- unlist(lapply(cases, `[[`, "z"), use.names = FALSE)
  z_range <- diff(range(all_z))
  z_limits <- c(min(all_z) - 0.12 * z_range, max(all_z) + 0.12 * z_range)
  density_objects <- lapply(cases, make_density_grid, z_limits = z_limits)
  grid <- do.call(rbind, lapply(density_objects, `[[`, "grid"))
  points <- do.call(rbind, cases)
  points$window_label <- unname(window_labels[points$window])

  window_order <- unname(window_labels[names(cases)])
  grid$window_label <- factor(grid$window_label, levels = window_order)
  points$window_label <- factor(points$window_label, levels = window_order)
  grid$estimator <- factor(grid$estimator, levels = c("Fitted HvMF", "Nonparametric KDE"))

  hdr_breaks <- c(0.25, 0.50, 0.75, 0.90)
  palette <- grDevices::colorRampPalette(
    c("#FFFFFF", "#E8F3F8", "#B9DDE7", "#72B7C5", "#2A788E", "#234B70")
  )(256L)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = theta_deg, y = z)) +
    ggplot2::geom_raster(ggplot2::aes(fill = density), interpolate = FALSE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = hdr_content), breaks = hdr_breaks,
      colour = "#17324D", linewidth = 0.28, alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = points,
      mapping = ggplot2::aes(x = theta_deg, y = z),
      inherit.aes = FALSE, shape = 21, size = 1.25, stroke = 0.28,
      colour = "#5A1A1A", fill = "#FFF7F2", alpha = 0.88
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(window_label), cols = ggplot2::vars(estimator)) +
    ggplot2::scale_fill_gradientn(
      colours = palette, limits = c(0, max(grid$density)),
      name = expression("Density w.r.t. " * dz * dtheta)
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 90, 180, 270, 360),
      labels = c("0°", "90°", "180°", "270°", "360°"), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = log(c(0.05, 0.1, 0.25, 0.5, 1, 2, 3)),
      labels = c("0.05", "0.10", "0.25", "0.50", "1", "2", "3"),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 360), ylim = z_limits, expand = FALSE) +
    ggplot2::labs(
      x = expression("Wind direction " * theta * " (0° and 360° identified)"),
      y = expression("Scaled speed " * r == s / bar(s) * " (log scale)"),
      title = "Risø 125 m wind data on the unwrapped log-speed cylinder",
      subtitle = sprintf(
        "%s: days %s; closest-to-noon observations, 1996-2003; HDR contours: 25%%, 50%%, 75%%, 90%%",
        pattern_id, paste(day_patterns[[pattern_id]], collapse = ", ")
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
  stem <- file.path(output_dir, paste0("risoe_125m_log_cylinder_density_contours_", pattern_id))
  figure_height <- 3 * length(cases) + 1.3
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 8.3, height = figure_height, units = "in")
  ggplot2::ggsave(
    paste0(stem, ".png"), p,
    width = 8.3, height = figure_height, units = "in", dpi = 320
  )

  summary <- do.call(rbind, Map(function(case, object) {
    data.frame(
      window = case$window[[1L]],
      pattern = pattern_id,
      n = nrow(case),
      speed_mean = case$speed_mean[[1L]],
      kappa_hat = object$fit$kappa,
      mu_r_hat = object$fit$chi,
      mu_theta_deg_hat = object$fit$theta_deg,
      kde_g_z = object$bandwidths[["g_z"]],
      kde_kappa_theta = object$bandwidths[["kappa_theta"]],
      kde_h_theta = 1 / sqrt(object$bandwidths[["kappa_theta"]]),
      z_plot_min = z_limits[[1L]],
      z_plot_max = z_limits[[2L]],
      r_plot_min = exp(z_limits[[1L]]),
      r_plot_max = exp(z_limits[[2L]]),
      full_parametric_mass = object$full_parametric_mass,
      full_nonparametric_mass = object$full_nonparametric_mass,
      parametric_mass_on_plot = object$parametric_mass_on_plot,
      nonparametric_mass_on_plot = object$nonparametric_mass_on_plot,
      jacobian_error = object$jacobian_error,
      kernel_mass_error = object$kernel_mass_error,
      seam_error_parametric = object$seam_error_parametric,
      seam_error_nonparametric = object$seam_error_nonparametric
    )
  }, cases, density_objects))
  utils::write.csv(summary, paste0(stem, "_summary.csv"), row.names = FALSE)
  summary
}

run_cylinder_plots <- function(
    input_nc = file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"),
    output_dir = file.path(repo_root, "real_data", "wind", "log_cylinder_density_contours"),
    patterns = names(day_patterns)) {
  invalid_patterns <- setdiff(patterns, names(day_patterns))
  if (length(invalid_patterns) > 0L) {
    stop("Unknown pattern(s): ", paste(invalid_patterns, collapse = ", "), call. = FALSE)
  }
  all_data <- load_risoe_concurrent(input_nc, fixed_tz = "UTC")
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  all_summaries <- lapply(patterns, function(pattern_id) {
    cases <- lapply(names(month_windows), function(window_id) {
      build_case(selected, window_id, pattern_id)
    })
    names(cases) <- names(month_windows)
    plot_pattern(cases, pattern_id, output_dir)
  })
  summary <- do.call(rbind, all_summaries)
  utils::write.csv(summary, file.path(output_dir, "risoe_125m_log_cylinder_density_contours_summary.csv"),
                   row.names = FALSE)
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  requested_patterns <- if (length(args) == 0L) names(day_patterns) else args
  result <- run_cylinder_plots(patterns = requested_patterns)
  print(result, row.names = FALSE, digits = 5)
}
