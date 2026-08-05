#!/usr/bin/env Rscript

# Fitted-HvMF and nonparametric Gaussian/von-Mises KDE density contours in
# z = log(speed / mean(speed)) for the cross-year Nov--Dec--Jan, start4 sample.
# The sample is exactly the one used in the 1996--2004 cross-year screening.

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

# Validated HvMF density, Gaussian log-speed KDE, HDR, and mass checks.
source(file.path(repo_root, "real_data", "wind", "plot_risoe_cylinder_density_contours.R"))

cross_year_months <- c(11L, 12L, 1L)
cross_year_days <- c(4L, 8L, 12L, 16L, 20L, 24L, 28L)
window_labels <- c(nov_dec_jan = "November + December + January")

build_cross_year_log_case <- function(selected, height_m) {
  speed_col <- sprintf("ws%d", height_m)
  direction_col <- sprintf("wd%d", height_m)
  cutoff <- if (height_m == 77L) as.Date("2005-01-01") else as.Date("2004-11-12")
  date <- as.Date(selected$datetime, tz = "UTC")
  rows <- selected[
    selected$year >= 1996L & selected$year <= 2004L &
      selected$month %in% cross_year_months & selected$day %in% cross_year_days &
      date < cutoff,
    , drop = FALSE
  ]
  valid <- is.finite(rows[[speed_col]]) & is.finite(rows[[direction_col]]) & rows[[speed_col]] > 0
  rows <- rows[valid, , drop = FALSE]
  if (nrow(rows) < 2L || anyDuplicated(as.Date(rows$datetime, tz = "UTC"))) {
    stop(sprintf("Invalid selected sample at %d m.", height_m), call. = FALSE)
  }
  speed <- as.numeric(rows[[speed_col]])
  theta <- ((as.numeric(rows[[direction_col]]) %% 360) + 360) %% 360 * pi / 180
  speed_mean <- mean(speed)
  r <- speed / speed_mean
  data.frame(
    window = "nov_dec_jan", pattern = "start4", datetime = rows$datetime,
    year = rows$year, month = rows$month, day = rows$day, speed = speed,
    speed_mean = speed_mean, r = r, z = log(r), theta = theta,
    theta_deg = theta * 180 / pi, x0 = cosh(r), x1 = sinh(r) * cos(theta),
    x2 = sinh(r) * sin(theta), stringsAsFactors = FALSE
  )
}

plot_cross_year_log_case <- function(case, height_m, output_dir) {
  plot_height_m <<- height_m
  coverage_note <- if (height_m == 125L) "1996--2004; 2004 through 11 November" else "1996--2004"
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
      title = sprintf("Risø %d m wind data on the unwrapped log-speed cylinder", height_m),
      subtitle = sprintf("start4: days %s;\nNovember--December--January, %s; HDR contours: 25%%, 50%%, 75%%, 90%%", paste(cross_year_days, collapse = ", "), coverage_note)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.spacing = grid::unit(.8, "lines"), strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey35"), strip.text = ggplot2::element_text(face = "bold"), aspect.ratio = 1, legend.position = "right", plot.title.position = "plot")
  stem <- file.path(output_dir, sprintf("risoe_%dm_nov_dec_jan_start4_log_density_contours", height_m))
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 8.3, height = 4.3, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 8.3, height = 4.3, units = "in", dpi = 320)
  data.frame(
    height_m = height_m, months = "11,12,1", pattern = "start4", n = nrow(case), n_2004 = sum(case$year == 2004L),
    speed_mean = case$speed_mean[[1L]], kappa_hat = object$fit$kappa, theta_deg_hat = object$fit$theta_deg,
    kde_g_z = object$bandwidths[["g_z"]], kde_kappa_theta = object$bandwidths[["kappa_theta"]],
    z_plot_min = z_limits[[1L]], z_plot_max = z_limits[[2L]], full_parametric_mass = object$full_parametric_mass,
    full_nonparametric_mass = object$full_nonparametric_mass, parametric_mass_on_plot = object$parametric_mass_on_plot,
    nonparametric_mass_on_plot = object$nonparametric_mass_on_plot, jacobian_error = object$jacobian_error,
    kernel_mass_error = object$kernel_mass_error, seam_error_parametric = object$seam_error_parametric,
    seam_error_nonparametric = object$seam_error_nonparametric, stringsAsFactors = FALSE
  )
}

run_nov_dec_jan_start4_log_density_contours <- function(
    output_dir = file.path(repo_root, "real_data", "wind", "month_diagnostics", "nov_dec_jan_start4_log_density_contours_1996_2004")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- select_noon_all_months(load_risoe_concurrent(file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"), fixed_tz = "UTC")
  cases <- lapply(c(77L, 125L), function(height_m) build_cross_year_log_case(selected, height_m))
  names(cases) <- c("77", "125")
  if (!identical(vapply(cases, nrow, integer(1)), c(`77` = 169L, `125` = 157L))) stop("Plot samples do not match the screening samples.", call. = FALSE)
  summary <- do.call(rbind, lapply(c(77L, 125L), function(height_m) plot_cross_year_log_case(cases[[as.character(height_m)]], height_m, output_dir)))
  utils::write.csv(summary, file.path(output_dir, "risoe_nov_dec_jan_start4_log_density_contours_summary.csv"), row.names = FALSE)
  summary
}

if (sys.nframe() == 0L) print(run_nov_dec_jan_start4_log_density_contours(), row.names = FALSE, digits = 6)
