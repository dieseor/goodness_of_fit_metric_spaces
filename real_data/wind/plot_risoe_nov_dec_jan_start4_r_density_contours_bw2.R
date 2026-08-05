#!/usr/bin/env Rscript

# As plot_risoe_nov_dec_jan_start4_r_density_contours.R, but with the
# nonparametric reflected product KDE smoothed by a factor of two in both
# effective bandwidths: h_r is doubled and h_theta is doubled, equivalently
# kappa_theta is divided by four.  The fitted HvMF panel is unchanged.

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

source(file.path(repo_root, "real_data", "wind", "plot_risoe_nov_dec_jan_start4_r_density_contours.R"))

smoothing_factor <- 2
automatic_r_bandwidth_selector <- select_r_bandwidths
select_r_bandwidths <- function(r, theta) {
  bandwidths <- automatic_r_bandwidth_selector(r, theta)
  bandwidths[["h_r"]] <- smoothing_factor * bandwidths[["h_r"]]
  bandwidths[["kappa_theta"]] <- bandwidths[["kappa_theta"]] / smoothing_factor^2
  bandwidths
}

plot_cross_year_case_bw2 <- function(case, height_m, output_dir) {
  plot_height_m <<- height_m
  coverage_note <- if (height_m == 125L) "1996--2004; 2004 through 11 November" else "1996--2004"
  r_max <- max(3.5, 1.08 * max(case$r))
  object <- make_r_density_grid(case, r_max = r_max)
  grid_data <- object$grid
  grid_data$estimator <- as.character(grid_data$estimator)
  grid_data$estimator[grid_data$estimator == "Nonparametric KDE"] <- "Nonparametric KDE; bandwidth × 2"
  grid_data$estimator <- factor(grid_data$estimator, levels = c("Fitted HvMF", "Nonparametric KDE; bandwidth × 2"))
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
      title = sprintf("Risø %d m wind data on the unwrapped cylinder", height_m),
      subtitle = sprintf("start4: days %s;\nNovember--December--January, %s;\nnonparametric KDE bandwidth × 2; HDR contours: 25%%, 50%%, 75%%, 90%%", paste(cross_year_days, collapse = ", "), coverage_note)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.spacing = grid::unit(.8, "lines"), strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey35"), strip.text = ggplot2::element_text(face = "bold"), aspect.ratio = 1, legend.position = "right", plot.title.position = "plot")
  stem <- file.path(output_dir, sprintf("risoe_%dm_nov_dec_jan_start4_r_density_contours_bw2", height_m))
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 8.3, height = 4.3, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 8.3, height = 4.3, units = "in", dpi = 320)
  data.frame(
    height_m = height_m, months = "11,12,1", pattern = "start4", n = nrow(case), n_2004 = sum(case$year == 2004L),
    kde_smoothing_factor = smoothing_factor, speed_mean = case$speed_mean[[1L]], kappa_hat = object$fit$kappa,
    theta_deg_hat = object$fit$theta_deg, kde_h_r = object$bandwidths[["h_r"]], kde_kappa_theta = object$bandwidths[["kappa_theta"]],
    r_plot_max = r_max, full_parametric_mass = object$full_parametric_mass, full_nonparametric_mass = object$full_nonparametric_mass,
    parametric_mass_on_plot = object$parametric_mass_on_plot, nonparametric_mass_on_plot = object$nonparametric_mass_on_plot,
    kernel_mass_error = object$kernel_mass_error, seam_error = object$seam_error, stringsAsFactors = FALSE
  )
}

run_nov_dec_jan_start4_r_density_contours_bw2 <- function(
    output_dir = file.path(repo_root, "real_data", "wind", "month_diagnostics", "nov_dec_jan_start4_r_density_contours_1996_2004_bw2")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- select_noon_all_months(load_risoe_concurrent(file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"), fixed_tz = "UTC")
  summary <- do.call(rbind, lapply(c(77L, 125L), function(height_m) {
    plot_cross_year_case_bw2(build_cross_year_case(selected, height_m), height_m, output_dir)
  }))
  utils::write.csv(summary, file.path(output_dir, "risoe_nov_dec_jan_start4_r_density_contours_bw2_summary.csv"), row.names = FALSE)
  summary
}

if (sys.nframe() == 0L) print(run_nov_dec_jan_start4_r_density_contours_bw2(), row.names = FALSE, digits = 6)
