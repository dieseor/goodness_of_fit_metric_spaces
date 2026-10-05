#!/usr/bin/env Rscript

# Fitted-HvMF and LCV-selected reflected Gaussian/von-Mises KDE contours in
# r = speed / mean(speed), for Risoe 77 m, May--June--July, start4.
# The sample matches the 1996--2004 screening; no observations are imputed.

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

source(file.path(repo_root, "real_data", "wind", "plot_risoe_77m_nov_dec_b5000_r_cylinder_density_contours.R"))

target_months <- c(5L, 6L, 7L)
target_days <- c(4L, 8L, 12L, 16L, 20L, 24L, 28L)

build_target_case <- function(selected) {
  date <- as.Date(selected$datetime, tz = "UTC")
  rows <- selected[
    selected$year >= 1996L & selected$year <= 2004L &
      selected$month %in% target_months & selected$day %in% target_days &
      date < as.Date("2004-11-12"),
    , drop = FALSE
  ]
  valid <- is.finite(rows$ws77) & is.finite(rows$wd77) & rows$ws77 > 0
  rows <- rows[valid, , drop = FALSE]
  if (nrow(rows) != 188 || anyDuplicated(as.Date(rows$datetime, tz = "UTC"))) {
    stop("The plotting sample does not match the 77 m May--June--July start4 screening sample.", call. = FALSE)
  }
  speed <- as.numeric(rows$ws77)
  theta <- ((as.numeric(rows$wd77) %% 360) + 360) %% 360 * pi / 180
  speed_mean <- mean(speed)
  r <- speed / speed_mean
  data.frame(
    datetime = rows$datetime, year = rows$year, month = rows$month, day = rows$day,
    speed = speed, speed_mean = speed_mean, r = r, theta = theta,
    theta_deg = theta * 180 / pi, x0 = cosh(r), x1 = sinh(r) * cos(theta),
    x2 = sinh(r) * sin(theta), stringsAsFactors = FALSE
  )
}

plot_target_case <- function(case, output_dir) {
  plot_height_m <<- 77L
  r_max <- max(3.5, 1.08 * max(case$r))
  object <- make_r_density_grid(case, r_max = r_max)
  grid_data <- object$grid
  grid_data$estimator <- factor(grid_data$estimator, levels = c("Fitted HvMF", "Nonparametric KDE"))
  palette <- grDevices::colorRampPalette(c("#FFFFFF", "#E8F3F8", "#B9DDE7", "#72B7C5", "#2A788E", "#234B70"))(256L)
  p <- ggplot2::ggplot(grid_data, ggplot2::aes(theta_deg, r)) +
    ggplot2::geom_raster(ggplot2::aes(fill = density), interpolate = FALSE) +
    ggplot2::geom_contour(ggplot2::aes(z = hdr_content), breaks = c(.25, .50, .75, .90), colour = "#17324D", linewidth = .28, alpha = .9) +
    ggplot2::geom_point(data = case, ggplot2::aes(theta_deg, r), inherit.aes = FALSE, shape = 21, size = 1.25, stroke = .28, colour = "#5A1A1A", fill = "#FFF7F2", alpha = .88) +
    ggplot2::facet_grid(cols = ggplot2::vars(estimator)) +
    ggplot2::scale_fill_gradientn(colours = palette, limits = c(0, max(grid_data$density)), name = expression("Density w.r.t. " * dr * dtheta)) +
    ggplot2::scale_x_continuous(breaks = c(0, 90, 180, 270, 360), labels = c("0°", "90°", "180°", "270°", "360°"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, 360), ylim = c(0, r_max), expand = FALSE) +
    ggplot2::labs(
      x = expression("Wind direction " * theta * " (0° and 360° identified)"),
      y = expression("Scaled speed " * r == s / bar(s)),
      title = "Risø 77 m wind data on the unwrapped cylinder",
      subtitle = sprintf("start4: days %s; May--June--July, 1996--2004; KDE bandwidth selected by LCV;\nHDR contours: 25%%, 50%%, 75%%, 90%%", paste(target_days, collapse = ", "))
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.spacing = grid::unit(.8, "lines"), strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey35"), strip.text = ggplot2::element_text(face = "bold"), aspect.ratio = 1, legend.position = "right", plot.title.position = "plot")
  stem <- file.path(output_dir, "risoe_77m_may_jun_jul_start4_r_density_contours")
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 8.3, height = 4.3, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 8.3, height = 4.3, units = "in", dpi = 320)
  data.frame(
    height_m = 77L, months = "5,6,7", pattern = "start4", n = nrow(case), n_2004 = sum(case$year == 2004L),
    speed_mean = case$speed_mean[[1L]], kappa_hat = object$fit$kappa, theta_deg_hat = object$fit$theta_deg,
    kde_h_r = object$bandwidths[["h_r"]], kde_kappa_theta = object$bandwidths[["kappa_theta"]], kde_h_theta = 1 / sqrt(object$bandwidths[["kappa_theta"]]),
    r_plot_max = r_max, full_parametric_mass = object$full_parametric_mass, full_nonparametric_mass = object$full_nonparametric_mass,
    parametric_mass_on_plot = object$parametric_mass_on_plot, nonparametric_mass_on_plot = object$nonparametric_mass_on_plot,
    kernel_mass_error = object$kernel_mass_error, seam_error = object$seam_error, stringsAsFactors = FALSE
  )
}

run_may_jun_jul_start4_r_density_contours <- function(
    output_dir = file.path(repo_root, "real_data", "wind", "month_diagnostics", "may_jun_jul_77m_start4_r_density_contours_1996_2004")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- select_noon_all_months(load_risoe_concurrent(file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"), fixed_tz = "UTC")
  summary <- plot_target_case(build_target_case(selected), output_dir)
  utils::write.csv(summary, file.path(output_dir, "risoe_77m_may_jun_jul_start4_r_density_contours_summary.csv"), row.names = FALSE)
  summary
}

if (sys.nframe() == 0L) print(run_may_jun_jul_start4_r_density_contours(), row.names = FALSE, digits = 6)
