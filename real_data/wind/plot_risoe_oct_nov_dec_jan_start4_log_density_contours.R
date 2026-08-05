#!/usr/bin/env Rscript

# Fitted-HvMF and LCV-selected Gaussian/von-Mises KDE density contours in
# z = log(speed / mean(speed)) for Risoe 125 m, Oct--Nov--Dec--Jan, start4.
# The sample matches the 1996--2004 screening: 125 m observations before
# 2004-11-12 only; no missing values are imputed.

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

source(file.path(repo_root, "real_data", "wind", "plot_risoe_cylinder_density_contours.R"))

target_months <- c(10L, 11L, 12L, 1L)
target_days <- c(4L, 8L, 12L, 16L, 20L, 24L, 28L)
window_labels <- c(oct_nov_dec_jan = "October + November + December + January")

build_target_case <- function(selected) {
  date <- as.Date(selected$datetime, tz = "UTC")
  rows <- selected[
    selected$year >= 1996L & selected$year <= 2004L &
      selected$month %in% target_months & selected$day %in% target_days &
      date < as.Date("2004-11-12"),
    , drop = FALSE
  ]
  valid <- is.finite(rows$ws125) & is.finite(rows$wd125) & rows$ws125 > 0
  rows <- rows[valid, , drop = FALSE]
  if (nrow(rows) != 220L || anyDuplicated(as.Date(rows$datetime, tz = "UTC"))) {
    stop("The plotting sample does not match the 125 m Oct--Nov--Dec--Jan start4 screening sample.", call. = FALSE)
  }
  speed <- as.numeric(rows$ws125)
  theta <- ((as.numeric(rows$wd125) %% 360) + 360) %% 360 * pi / 180
  speed_mean <- mean(speed)
  r <- speed / speed_mean
  data.frame(
    window = "oct_nov_dec_jan", pattern = "start4", datetime = rows$datetime,
    year = rows$year, month = rows$month, day = rows$day, speed = speed,
    speed_mean = speed_mean, r = r, z = log(r), theta = theta,
    theta_deg = theta * 180 / pi, x0 = cosh(r), x1 = sinh(r) * cos(theta),
    x2 = sinh(r) * sin(theta), stringsAsFactors = FALSE
  )
}

plot_target_case <- function(case, output_dir) {
  plot_height_m <<- 125L
  z_range <- diff(range(case$z))
  z_limits <- c(min(case$z) - .12 * z_range, max(case$z) + .12 * z_range)
  object <- make_density_grid(case, z_limits = z_limits)
  grid_data <- object$grid
  grid_data$estimator <- factor(grid_data$estimator, levels = c("Fitted HvMF", "Nonparametric KDE"))
  palette <- grDevices::colorRampPalette(c("#FFFFFF", "#E8F3F8", "#B9DDE7", "#72B7C5", "#2A788E", "#234B70"))(256L)
  p <- ggplot2::ggplot(grid_data, ggplot2::aes(theta_deg, z)) +
    ggplot2::geom_raster(ggplot2::aes(fill = density), interpolate = FALSE) +
    ggplot2::geom_contour(ggplot2::aes(z = hdr_content), breaks = c(.25, .50, .75, .90), colour = "#17324D", linewidth = .28, alpha = .9) +
    ggplot2::geom_point(data = case, ggplot2::aes(theta_deg, z), inherit.aes = FALSE, shape = 21, size = 1.25, stroke = .28, colour = "#5A1A1A", fill = "#FFF7F2", alpha = .88) +
    ggplot2::facet_grid(cols = ggplot2::vars(estimator)) +
    ggplot2::scale_fill_gradientn(colours = palette, limits = c(0, max(grid_data$density)), name = expression("Density w.r.t. " * dz * dtheta)) +
    ggplot2::scale_x_continuous(breaks = c(0, 90, 180, 270, 360), labels = c("0°", "90°", "180°", "270°", "360°"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = log(c(.05, .1, .25, .5, 1, 2, 3)), labels = c("0.05", "0.10", "0.25", "0.50", "1", "2", "3"), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, 360), ylim = z_limits, expand = FALSE) +
    ggplot2::labs(
      x = expression("Wind direction " * theta * " (0° and 360° identified)"),
      y = expression("Scaled speed " * r == s / bar(s) * " (log scale)"),
      title = "Risø 125 m wind data on the unwrapped log-speed cylinder",
      subtitle = sprintf("start4: days %s;\nOctober--November--December--January, 1996--2004; 2004 through 11 November; KDE bandwidth selected by LCV;\nHDR contours: 25%%, 50%%, 75%%, 90%%", paste(target_days, collapse = ", "))
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.spacing = grid::unit(.8, "lines"), strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey35"), strip.text = ggplot2::element_text(face = "bold"), aspect.ratio = 1, legend.position = "right", plot.title.position = "plot")
  stem <- file.path(output_dir, "risoe_125m_oct_nov_dec_jan_start4_log_density_contours")
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 8.3, height = 4.8, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 8.3, height = 4.8, units = "in", dpi = 320)
  data.frame(
    height_m = 125L, months = "10,11,12,1", pattern = "start4", n = nrow(case), n_2004 = sum(case$year == 2004L),
    speed_mean = case$speed_mean[[1L]], kappa_hat = object$fit$kappa, theta_deg_hat = object$fit$theta_deg,
    kde_g_z = object$bandwidths[["g_z"]], kde_kappa_theta = object$bandwidths[["kappa_theta"]],
    kde_h_theta = 1 / sqrt(object$bandwidths[["kappa_theta"]]), z_plot_min = z_limits[[1L]], z_plot_max = z_limits[[2L]],
    full_parametric_mass = object$full_parametric_mass, full_nonparametric_mass = object$full_nonparametric_mass,
    parametric_mass_on_plot = object$parametric_mass_on_plot, nonparametric_mass_on_plot = object$nonparametric_mass_on_plot,
    jacobian_error = object$jacobian_error, kernel_mass_error = object$kernel_mass_error,
    seam_error_parametric = object$seam_error_parametric, seam_error_nonparametric = object$seam_error_nonparametric,
    stringsAsFactors = FALSE
  )
}

run_oct_nov_dec_jan_start4_log_density_contours <- function(
    output_dir = file.path(repo_root, "real_data", "wind", "month_diagnostics", "oct_nov_dec_jan_125m_start4_log_density_contours_1996_2004")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- select_noon_all_months(load_risoe_concurrent(file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"), fixed_tz = "UTC")
  summary <- plot_target_case(build_target_case(selected), output_dir)
  utils::write.csv(summary, file.path(output_dir, "risoe_125m_oct_nov_dec_jan_start4_log_density_contours_summary.csv"), row.names = FALSE)
  summary
}

if (sys.nframe() == 0L) print(run_oct_nov_dec_jan_start4_log_density_contours(), row.names = FALSE, digits = 6)
