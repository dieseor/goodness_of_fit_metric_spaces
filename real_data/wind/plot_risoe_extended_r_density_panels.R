#!/usr/bin/env Rscript

# Separate, publication-ready density panels in the original r=s/mean(s)
# coordinate for the November--December--January and May--June--July start4
# samples at 77 m.
#
# The November--December--January KDEs use the manual bandwidths from
# plot_risoe_nov_dec_jan_start4_r_density_contours_manual_bw.R.
# The May--June--July KDE uses the bandwidth selected by the existing LCV code.
#
# All four panels share axes, HDR levels, point styling, dimensions, and
# typography with plot_risoe_paper_r_density_panels.R. The two
# November--December--January panels share one density scale, while the two
# May--June--July panels share a second density scale. Each scale is shown to
# the right of the corresponding nonparametric panel: November--December--
# January and May--June--July, both at 77 m.

script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[grepl("^--file=", script_argument)]
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_argument),
  winslash = "/",
  mustWork = TRUE
)
repo_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
setwd(repo_root)

# Load the functions used to construct the two November--December--January
# samples and the common HvMF/KDE density machinery.
source(file.path(
  repo_root,
  "real_data", "wind",
  "plot_risoe_nov_dec_jan_start4_r_density_contours.R"
))

# Preserve the cross-year constructor before loading the May--June--July file.
build_nov_dec_jan_case <- build_cross_year_case
nov_dec_jan_days <- cross_year_days

# Load the exact May--June--July sample constructor.
source(file.path(
  repo_root,
  "real_data", "wind",
  "plot_risoe_may_jun_jul_77m_start4_r_density_contours.R"
))

paper_output_dir <- file.path(
  "/Users/Diego/Documents/LaTEX Github/tex_GOF_metric_spaces",
  "AoS", "img"
)

# These are the manual bandwidths used in
# plot_risoe_nov_dec_jan_start4_r_density_contours_manual_bw.R.
manual_bandwidths <- list(
  `77` = c(h_r = 0.30, kappa_theta = 1.50)
)

make_density_object <- function(case, bandwidths, n_theta = 241L,
                                n_r = 181L) {
  theta_values <- seq(0, 2 * pi, length.out = n_theta)
  r_values <- seq(0, 2.5, length.out = n_r)
  grid <- expand.grid(theta = theta_values, r = r_values)

  fit <- hvmf_mle_h2(as.matrix(case[, c("x0", "x1", "x2")]))
  parametric <- hvmf_r_cylinder_density(
    grid$r, grid$theta, fit
  )
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
  grob <- ggplot2::ggplotGrob(
    plot + ggplot2::theme(legend.position = "right")
  )
  index <- grep("^guide-box", grob$layout$name)
  index <- index[vapply(
    grob$grobs[index],
    function(x) !inherits(x, "zeroGrob"),
    logical(1)
  )]
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
    ggplot2::ggplotGrob(
      plot + ggplot2::theme(legend.position = "none")
    ),
    extract_plot_legend(plot),
    nrow = 1L,
    widths = grid::unit(c(panel_width, legend_width), "in")
  )

  ggplot2::ggsave(
    output_file,
    combined,
    width = panel_width + legend_width,
    height = figure_height,
    units = "in"
  )
  ggplot2::ggsave(
    sub("[.]pdf$", ".png", output_file),
    combined,
    width = panel_width + legend_width,
    height = figure_height,
    units = "in",
    dpi = 320
  )
}

plot_density_panel <- function(object, estimator, fill_max, show_legend,
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
    ggplot2::geom_raster(
      ggplot2::aes(fill = density),
      interpolate = FALSE
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = hdr_content),
      breaks = c(0.25, 0.50, 0.75, 0.90),
      colour = "#17324D",
      linewidth = 0.32,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = object$case,
      ggplot2::aes(theta_deg, r),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.05,
      stroke = 0.25,
      colour = "#5A1A1A",
      fill = "#FFF7F2",
      alpha = 0.86
    ) +
    ggplot2::scale_fill_gradientn(
      colours = palette,
      limits = c(0, fill_max),
      name = "Density"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 90, 180, 270, 360),
      labels = c(
        "0",
        expression(pi / 2),
        expression(pi),
        expression(3 * pi / 2),
        expression(2 * pi)
      ),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 2.5, by = 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 360),
      ylim = c(0, 2.5),
      expand = FALSE
    ) +
    ggplot2::labs(
      x = expression(theta),
      y = expression(s / bar(s))
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      aspect.ratio = 1,
      legend.position = if (show_legend) "right" else "none",
      axis.title = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(
        size = 11,
        colour = if (show_x_title) "black" else "transparent"
      ),
      axis.title.y = ggplot2::element_text(
        size = 11,
        colour = if (show_y_title) "black" else "transparent"
      ),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  if (show_legend) {
    draw_plot_with_external_legend(p, output_file)
  } else {
    ggplot2::ggsave(
      output_file,
      p,
      width = 3.25,
      height = 3.15,
      units = "in"
    )
    ggplot2::ggsave(
      sub("[.]pdf$", ".png", output_file),
      p,
      width = 3.25,
      height = 3.15,
      units = "in",
      dpi = 320
    )
  }
}

run_extended_density_panels <- function(output_dir = paper_output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"),
    fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")

  cases <- list(
    nov_dec_jan_77m_start4 = build_nov_dec_jan_case(selected, 77L),
    may_jun_jul_77m_start4 = build_target_case(selected)
  )

  observed_n <- vapply(cases, nrow, integer(1))
  if (any(observed_n <= 0L)) {
    stop("At least one plotting sample is empty.", call. = FALSE)
  }
  if (observed_n[["may_jun_jul_77m_start4"]] != 188L) {
    stop(
      "The May--June--July 77 m sample does not match the screening sample (n=188).",
      call. = FALSE
    )
  }

  # The May--June--July bandwidth is selected by the original LCV routine.
  plot_height_m <<- 77L
  may_jun_jul_bandwidths <- select_r_bandwidths(
    cases$may_jun_jul_77m_start4$r,
    cases$may_jun_jul_77m_start4$theta
  )

  objects <- list(
    nov_dec_jan_77m_start4 = make_density_object(
      cases$nov_dec_jan_77m_start4,
      manual_bandwidths[["77"]]
    ),
    may_jun_jul_77m_start4 = make_density_object(
      cases$may_jun_jul_77m_start4,
      may_jun_jul_bandwidths
    )
  )

  # Two common fill scales are used: one for all November--December--
  # January panels and one for both May--June--July panels.
  nov_dec_jan_ids <- "nov_dec_jan_77m_start4"
  may_jun_jul_ids <- "may_jun_jul_77m_start4"

  nov_dec_jan_fill_max <- max(vapply(
    objects[nov_dec_jan_ids],
    function(x) max(c(x$parametric, x$nonparametric), na.rm = TRUE),
    numeric(1)
  ))
  may_jun_jul_fill_max <- max(vapply(
    objects[may_jun_jul_ids],
    function(x) max(c(x$parametric, x$nonparametric), na.rm = TRUE),
    numeric(1)
  ))

  for (dataset_id in names(objects)) {
    for (estimator in c("parametric", "nonparametric")) {
      # Show one legend for each density scale, always to the right of
      # the corresponding nonparametric panel.
      show_legend <-
        estimator == "nonparametric" &&
        dataset_id %in% c(
          "nov_dec_jan_77m_start4",
          "may_jun_jul_77m_start4"
        )

      fill_max <- if (dataset_id %in% nov_dec_jan_ids) {
        nov_dec_jan_fill_max
      } else {
        may_jun_jul_fill_max
      }

      # Same axis-title layout as plot_risoe_paper_r_density_panels.R:
      # y title on the first dataset, x title on the middle dataset.
      show_x_title <- dataset_id == "nov_dec_jan_77m_start4"
      show_y_title <- dataset_id == "nov_dec_jan_77m_start4"

      output_file <- file.path(
        output_dir,
        sprintf(
          "risoe_%s_%s_r_density.pdf",
          dataset_id,
          estimator
        )
      )

      plot_density_panel(
        objects[[dataset_id]],
        estimator,
        fill_max,
        show_legend,
        show_x_title,
        show_y_title,
        output_file
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
      mu1_hat = unname(object$fit$mu[[1L]]),
      mu2_hat = unname(object$fit$mu[[2L]]),
      mu3_hat = unname(object$fit$mu[[3L]]),
      kde_h_r = object$bandwidths[["h_r"]],
      kde_kappa_theta = object$bandwidths[["kappa_theta"]],
      density_scale_group = if (dataset_id %in% nov_dec_jan_ids) {
        "nov_dec_jan"
      } else {
        "may_jun_jul"
      },
      group_fill_max = if (dataset_id %in% nov_dec_jan_ids) {
        nov_dec_jan_fill_max
      } else {
        may_jun_jul_fill_max
      },
      stringsAsFactors = FALSE
    )
  }))

  utils::write.csv(
    summary,
    file.path(output_dir, "risoe_extended_r_density_panels_summary.csv"),
    row.names = FALSE
  )

  invisible(summary)
}

if (sys.nframe() == 0L) {
  result <- run_extended_density_panels()
  print(result, row.names = FALSE, digits = 6)
}
