#!/usr/bin/env Rscript

# Separate, publication-ready density panels in the original r=s/mean(s)
# coordinate for the three Risoe samples reported in the paper. All panels
# share axes, fill scale, and HDR levels. Densities are with respect to
# dtheta dr; hence the fitted HvMF density includes the Jacobian sinh(r).

source(file.path(
  "real_data", "wind", "plot_risoe_125m_nov_dec_b5000_r_cylinder_density_contours.R"
))

paper_output_dir <- file.path(
  "/Users/Diego/Documents/LaTEX Github/tex_GOF_metric_spaces", "AoS", "img"
)
paper_days <- c(4L, 8L, 12L, 16L, 20L, 24L, 28L)

build_paper_case <- function(selected, height_m, months, dataset_id) {
  plot_height_m <<- as.integer(height_m)
  plot_speed_col <<- sprintf("ws%d", height_m)
  plot_direction_col <<- sprintf("wd%d", height_m)
  plot_clean_before <<- if (height_m == 77L) {
    as.Date("2007-08-01")
  } else {
    as.Date("2004-12-01")
  }
  month_windows[[dataset_id]] <<- as.integer(months)
  day_patterns[["paper_days"]] <<- paper_days
  build_case(selected, dataset_id, "paper_days", years = 1996:2003)
}

make_paper_density_object <- function(case, n_theta = 241L, n_r = 181L) {
  theta_values <- seq(0, 2 * pi, length.out = n_theta)
  r_values <- seq(0, 2.5, length.out = n_r)
  grid <- expand.grid(theta = theta_values, r = r_values)
  fit <- hvmf_mle_h2(as.matrix(case[, c("x0", "x1", "x2")]))
  bandwidths <- select_r_bandwidths(case$r, case$theta)
  parametric <- hvmf_r_cylinder_density(grid$r, grid$theta, fit)
  nonparametric <- kde_r_cylinder_density(
    grid$r, grid$theta, case$r, case$theta, bandwidths
  )
  list(
    case = case,
    grid = grid,
    fit = fit,
    bandwidths = bandwidths,
    parametric = parametric,
    nonparametric = nonparametric
  )
}

extract_plot_legend <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot + ggplot2::theme(legend.position = "right"))
  index <- grep("^guide-box", grob$layout$name)
  index <- index[vapply(grob$grobs[index], function(x) {
    !inherits(x, "zeroGrob")
  }, logical(1))]
  if (length(index) != 1L) {
    stop("Could not extract a unique density legend.", call. = FALSE)
  }
  grob$grobs[[index]]
}

draw_plot_with_external_legend <- function(plot, output_file) {
  panel_width <- 3.25
  legend_width <- 1.20
  figure_height <- 3.15
  combined <- gridExtra::arrangeGrob(
    ggplot2::ggplotGrob(plot + ggplot2::theme(legend.position = "none")),
    extract_plot_legend(plot),
    nrow = 1L,
    widths = grid::unit(c(panel_width, legend_width), "in")
  )
  ggplot2::ggsave(
    output_file, combined, width = panel_width + legend_width,
    height = figure_height, units = "in"
  )
  ggplot2::ggsave(
    sub("[.]pdf$", ".png", output_file), combined,
    width = panel_width + legend_width, height = figure_height,
    units = "in", dpi = 320
  )
}

plot_paper_panel <- function(object, estimator, fill_max, show_legend,
                             show_x_title, show_y_title, output_file) {
  density <- object[[estimator]]
  grid <- object$grid
  grid$density <- density
  grid$hdr_content <- density_to_hdr_content(density)
  grid$theta_deg <- grid$theta * 180 / pi
  palette <- grDevices::colorRampPalette(
    c("#FFFFFF", "#E8F3F8", "#B9DDE7", "#72B7C5", "#2A788E", "#234B70")
  )(256L)

  p <- ggplot2::ggplot(grid, ggplot2::aes(theta_deg, r)) +
    ggplot2::geom_raster(ggplot2::aes(fill = density), interpolate = FALSE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = hdr_content), breaks = c(0.25, 0.50, 0.75, 0.90),
      colour = "#17324D", linewidth = 0.32, alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = object$case, ggplot2::aes(theta_deg, r), inherit.aes = FALSE,
      shape = 21, size = 1.05, stroke = 0.25,
      colour = "#5A1A1A", fill = "#FFF7F2", alpha = 0.86
    ) +
    ggplot2::scale_fill_gradientn(
      colours = palette, limits = c(0, fill_max),
      name = "Density"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 90, 180, 270, 360),
      labels = c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 2.5, by = 0.5), expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 360), ylim = c(0, 2.5), expand = FALSE) +
    ggplot2::labs(x = expression(theta), y = expression(s/bar(s))) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      aspect.ratio = 1,
      legend.position = if (show_legend) "right" else "none",
      axis.title = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(
        size = 11, colour = if (show_x_title) "black" else "transparent"
      ),
      axis.title.y = ggplot2::element_text(
        size = 11, colour = if (show_y_title) "black" else "transparent"
      ),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  if (show_legend) {
    draw_plot_with_external_legend(p, output_file)
  } else {
    ggplot2::ggsave(output_file, p, width = 3.25, height = 3.15, units = "in")
    ggplot2::ggsave(
      sub("[.]pdf$", ".png", output_file), p,
      width = 3.25, height = 3.15, units = "in", dpi = 320
    )
  }
}

run_paper_density_panels <- function(output_dir = paper_output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  cases <- list(
    nov_dec_77m = build_paper_case(selected, 77L, c(11L, 12L), "nov_dec_77m"),
    nov_dec_125m = build_paper_case(selected, 125L, c(11L, 12L), "nov_dec_125m"),
    may_jun_125m = build_paper_case(selected, 125L, c(5L, 6L), "may_jun_125m")
  )
  expected_n <- c(nov_dec_77m = 100L, nov_dec_125m = 101L, may_jun_125m = 108L)
  observed_n <- vapply(cases, nrow, integer(1))
  if (!identical(observed_n, expected_n)) {
    stop("Paper-panel samples do not match the B=5000 analyses.", call. = FALSE)
  }
  objects <- lapply(cases, make_paper_density_object)
  fill_max <- max(vapply(objects, function(x) {
    max(c(x$parametric, x$nonparametric))
  }, numeric(1)))

  for (dataset_id in names(objects)) {
    for (estimator in c("parametric", "nonparametric")) {
      show_legend <- dataset_id == "may_jun_125m" && estimator == "nonparametric"
      show_x_title <- dataset_id == "nov_dec_125m"
      show_y_title <- dataset_id == "nov_dec_77m"
      output_file <- file.path(
        output_dir,
        sprintf("risoe_%s_%s_r_density.pdf", dataset_id, estimator)
      )
      plot_paper_panel(
        objects[[dataset_id]], estimator, fill_max, show_legend,
        show_x_title, show_y_title, output_file
      )
    }
  }

  summary <- do.call(rbind, lapply(names(objects), function(dataset_id) {
    object <- objects[[dataset_id]]
    data.frame(
      dataset = dataset_id,
      n = nrow(object$case),
      speed_mean = object$case$speed_mean[[1L]],
      kappa_hat = object$fit$kappa,
      theta_deg_hat = object$fit$theta_deg,
      kde_h_r = object$bandwidths[["h_r"]],
      kde_kappa_theta = object$bandwidths[["kappa_theta"]],
      common_fill_max = fill_max
    )
  }))
  invisible(summary)
}

if (sys.nframe() == 0L) {
  result <- run_paper_density_panels()
  print(result, row.names = FALSE, digits = 6)
}
