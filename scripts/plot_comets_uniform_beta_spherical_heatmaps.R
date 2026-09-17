#!/usr/bin/env Rscript

# Heat maps of the fitted uniform-beta (UB) densities on S^2, overlaid with
# the observed comet orbital poles. As in the paper's other spherical plots,
# rear points and contours are shown with reduced opacity. The heat-map fill
# itself represents the density on the front surface, so that two density
# values are never blended at one projected location. The saved paper
# parameters are never refitted.

resolve_comets_ub_heatmap_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_comets_ub_heatmap_path(
  "scripts", "plot_comets_uniform_beta_spherical_contours.R"
))

comets_ub_heatmap_scale <- function(theta, reference_size = 200001L) {
  reference_size <- as.integer(reference_size)
  if (!is.finite(reference_size) || reference_size < 10001L) {
    stop("`reference_size` must be at least 10001.")
  }

  y_reference <- (seq_len(reference_size) - 0.5) / reference_size
  density_reference <- comets_ub_density_y(y_reference, theta) / (4 * pi)
  lower <- min(density_reference)
  upper <- max(density_reference)
  if (!is.finite(lower) || lower <= 0 || !is.finite(upper) || upper <= lower) {
    stop("Could not construct a finite nondegenerate density scale.")
  }
  list(
    lower = lower,
    upper = upper,
    log_lower = log10(lower),
    log_upper = log10(upper),
    reference_size = reference_size
  )
}

comets_ub_heatmap_raster <- function(theta, camera, palette, scale,
                                     n_pixels = 701L,
                                     surface = c("front", "rear"),
                                     alpha_floor = 0,
                                     alpha_ceiling = 1) {
  surface <- match.arg(surface)
  n_pixels <- as.integer(n_pixels)
  if (!is.finite(n_pixels) || n_pixels < 201L) {
    stop("`n_pixels` must be an integer of at least 201.")
  }

  x_seq <- seq(-1, 1, length.out = n_pixels)
  y_seq <- seq(1, -1, length.out = n_pixels)
  plane <- expand.grid(x = x_seq, y = y_seq)
  radius_squared <- plane$x^2 + plane$y^2
  inside <- radius_squared <= 1

  colors <- rep("#00000000", nrow(plane))
  depth <- sqrt(pmax(0, 1 - radius_squared[inside]))
  if (identical(surface, "rear")) depth <- -depth
  xyz <- tcrossprod(plane$x[inside], camera$right) +
    tcrossprod(plane$y[inside], camera$up) +
    tcrossprod(depth, camera$view)
  z <- pmin(pmax(drop(xyz %*% theta$mu), -1), 1)

  density <- comets_ub_density_y((z + 1) / 2, theta) / (4 * pi)
  color_value <- (log10(density) - scale$log_lower) /
    (scale$log_upper - scale$log_lower)
  color_value <- pmin(pmax(color_value, 0), 1)
  color_index <- pmin(
    length(palette),
    pmax(1L, 1L + floor(color_value * (length(palette) - 1L)))
  )
  alpha_value <- alpha_floor +
    (alpha_ceiling - alpha_floor) * sqrt(color_value)
  rgb_value <- grDevices::col2rgb(palette[color_index]) / 255
  colors[inside] <- grDevices::rgb(
    rgb_value[1L, ], rgb_value[2L, ], rgb_value[3L, ],
    alpha = alpha_value
  )

  as.raster(matrix(colors, nrow = n_pixels, ncol = n_pixels, byrow = TRUE))
}

comets_ub_draw_heatmap_graticule <- function(camera) {
  lon_dense <- seq(-180, 180, length.out = 721L)
  lat_dense <- seq(-90, 90, length.out = 361L)
  for (latitude in c(-60, -30, 0, 30, 60)) {
    comets_ub_draw_curve(
      comets_ub_lon_lat_xyz(lon_dense, rep(latitude, length(lon_dense))),
      camera, color = "#FFFFFF", lwd = 0.45, front = TRUE, alpha = 0.58
    )
  }
  for (longitude in seq(-180, 120, by = 60)) {
    comets_ub_draw_curve(
      comets_ub_lon_lat_xyz(rep(longitude, length(lat_dense)), lat_dense),
      camera, color = "#FFFFFF", lwd = 0.45, front = TRUE, alpha = 0.58
    )
  }
  invisible(NULL)
}

comets_ub_draw_heatmap_points <- function(x, camera, front,
                                          back_alpha = 0.16) {
  projected <- comets_ub_project(x, camera)
  keep <- if (isTRUE(front)) projected$depth >= 0 else projected$depth < 0
  if (any(keep)) {
    alpha <- if (isTRUE(front)) 0.88 else back_alpha
    graphics::points(
      projected$x[keep], projected$y[keep],
      pch = 21, cex = 0.37, lwd = 0.35,
      col = grDevices::adjustcolor("#5A1A1A", alpha.f = alpha),
      bg = grDevices::adjustcolor("#FFF7F2", alpha.f = alpha)
    )
  }
  invisible(NULL)
}

comets_ub_draw_heatmap_panel <- function(x, theta, camera, title,
                                         palette, scale, boundaries,
                                         hdr_colors, n_pixels,
                                         back_alpha = 0.16,
                                         rear_density_alpha = 0.28) {
  rear_heatmap <- comets_ub_heatmap_raster(
    theta = theta,
    camera = camera,
    palette = palette,
    scale = scale,
    n_pixels = n_pixels,
    surface = "rear",
    alpha_floor = 0,
    alpha_ceiling = rear_density_alpha
  )
  front_heatmap <- comets_ub_heatmap_raster(
    theta = theta,
    camera = camera,
    palette = palette,
    scale = scale,
    n_pixels = n_pixels,
    surface = "front",
    alpha_floor = 0.10,
    alpha_ceiling = 0.95
  )

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-1.035, 1.035), ylim = c(-1.035, 1.035),
    asp = 1, xaxs = "i", yaxs = "i"
  )
  graphics::rasterImage(rear_heatmap, -1, -1, 1, 1, interpolate = TRUE)
  graphics::rasterImage(front_heatmap, -1, -1, 1, 1, interpolate = TRUE)
  comets_ub_draw_graticule(
    camera, front = FALSE, alpha = 0.12
  )
  comets_ub_draw_heatmap_graticule(camera)
  curves <- lapply(seq_len(nrow(boundaries)), function(i) {
    list(
      level = boundaries$probability[[i]],
      xyz = comets_ub_circle_xyz(theta$mu, boundaries$z[[i]])
    )
  })
  probabilities <- sort(unique(boundaries$probability))
  label_specs <- sunspots_joint_contour_label_specs(
    curves = curves,
    thresholds = probabilities,
    hdr_levels = probabilities,
    camera = camera,
    cex = 0.54
  )
  for (front in c(FALSE, TRUE)) {
    sunspots_joint_draw_density_contours_inline(
      curves = curves,
      thresholds = probabilities,
      camera = camera,
      colors = hdr_colors,
      lty = 1,
      lwd = 0.75,
      front = front,
      back_alpha = back_alpha,
      label_specs = label_specs
    )
    comets_ub_draw_heatmap_points(
      x, camera, front = front, back_alpha = back_alpha
    )
  }
  sunspots_joint_draw_contour_labels(
    label_specs = label_specs,
    colors = hdr_colors,
    cex = 0.54
  )
  graphics::symbols(
    0, 0, circles = 1, inches = FALSE, add = TRUE,
    fg = "#111111", lwd = 0.85
  )
  graphics::title(main = title, cex.main = 0.82, line = 0.25)
  invisible(NULL)
}

comets_ub_draw_heatmap_legend <- function(palette, scale, dataset_label) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  color_bar <- as.raster(matrix(palette, nrow = 1L))
  x_left <- 0.12
  x_right <- 0.88
  graphics::rasterImage(color_bar, x_left, 0.53, x_right, 0.73, interpolate = TRUE)
  graphics::rect(x_left, 0.53, x_right, 0.73, border = "#333333", lwd = 0.55)

  log_ticks <- pretty(c(scale$log_lower, scale$log_upper), n = 6L)
  log_ticks <- log_ticks[
    log_ticks >= scale$log_lower - 1e-12 &
      log_ticks <= scale$log_upper + 1e-12
  ]
  density_ticks <- 10^log_ticks
  tick_x <- x_left + (x_right - x_left) *
    (log_ticks - scale$log_lower) / (scale$log_upper - scale$log_lower)
  graphics::segments(tick_x, 0.49, tick_x, 0.53, col = "#333333", lwd = 0.55)
  graphics::text(
    tick_x, 0.43,
    labels = format(signif(density_ticks, 3L), trim = TRUE),
    cex = 0.68
  )
  graphics::text(
    0.5, 0.14,
    labels = "log f(x); rear surface faded",
    cex = 0.61
  )
  invisible(NULL)
}

comets_ub_render_heatmaps <- function(samples, thetas, scales, boundaries,
                                      output_path, device,
                                      palette, hdr_colors, n_pixels,
                                      angle_from_mu, back_alpha,
                                      rear_density_alpha) {
  if (identical(device, "pdf")) {
    grDevices::pdf(
      output_path, width = 7.2, height = 4.18,
      onefile = FALSE, paper = "special", version = "1.4",
      useDingbats = FALSE
    )
  } else if (identical(device, "png")) {
    grDevices::png(output_path, width = 2300, height = 1335, res = 300, bg = "white")
  } else {
    stop("Unknown graphics device: ", device)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::layout(
    matrix(c(1L, 2L, 3L, 4L), nrow = 2L, byrow = TRUE),
    heights = c(1, 0.21)
  )
  common_camera <- comets_ub_camera(
    thetas$short$mu,
    angle_from_mu = angle_from_mu
  )
  graphics::par(mar = c(0, 0.15, 1.35, 0.15), oma = c(0, 0, 0, 0))
  comets_ub_draw_heatmap_panel(
    samples$short, thetas$short, common_camera,
    sprintf("Short-period comets (n = %d)", nrow(samples$short)),
    palette, scales$short, boundaries$short, hdr_colors, n_pixels,
    back_alpha, rear_density_alpha
  )
  comets_ub_draw_heatmap_panel(
    samples$long, thetas$long, common_camera,
    sprintf("Long-period comets (n = %d)", nrow(samples$long)),
    palette, scales$long, boundaries$long, hdr_colors, n_pixels,
    back_alpha, rear_density_alpha
  )
  graphics::par(mar = c(0, 0, 0, 0))
  comets_ub_draw_heatmap_legend(palette, scales$short, "Short period")
  comets_ub_draw_heatmap_legend(palette, scales$long, "Long period")
  invisible(output_path)
}

run_comets_uniform_beta_spherical_heatmaps <- function(
    fit_root = file.path(
      "real_data", "reruns", "paper_main_realdata_B1000_3cores_20260831_113532",
      "comets", "uniform_beta"
    ),
    output_dir = file.path(
      "real_data", "comets", "mixture", "uniform_beta_spherical_heatmaps"
    ),
    n_pixels = 701L,
    angle_from_mu = 35,
    back_alpha = 0.16,
    rear_density_alpha = 0.28,
    hdr_probabilities = c(0.50, 0.80, 0.95)) {

  parameter_paths <- list(
    short = file.path(fit_root, "01_short_period_uniform_beta_mixture", "theta_hat.csv"),
    long = file.path(fit_root, "02_long_period_uniform_beta_mixture", "theta_hat.csv")
  )
  thetas <- lapply(parameter_paths, comets_ub_read_theta)

  comet_data <- load_comets_real_data(finite_normals = "both")
  samples <- list(
    short = as.matrix(comet_data$short$normal),
    long = as.matrix(comet_data$long$normal)
  )
  if (!identical(vapply(samples, nrow, integer(1L)), c(short = 784L, long = 610L))) {
    stop("The filtered comet samples do not have the paper sample sizes 784 and 610.")
  }

  if (!identical(sort(unique(as.numeric(hdr_probabilities))), c(0.50, 0.80, 0.95))) {
    stop("The paper-matched contour style requires HDR probabilities 0.50, 0.80, and 0.95.")
  }
  hdr_probabilities <- c(0.50, 0.80, 0.95)
  hdr_colors <- c("#1f78b4", "#e31a1c", "#1B5E20")
  palette <- grDevices::colorRampPalette(c(
    "#FFFFFF", "#E8F3F8", "#B9DDE7", "#72B7C5", "#2A788E", "#234B70"
  ))(256L)
  scales <- lapply(thetas, comets_ub_heatmap_scale)
  boundaries <- lapply(
    thetas,
    comets_ub_hdr_boundaries,
    probabilities = hdr_probabilities
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir, "comets_uniform_beta_spherical_heatmaps.pdf")
  png_path <- file.path(output_dir, "comets_uniform_beta_spherical_heatmaps.png")

  comets_ub_render_heatmaps(
    samples, thetas, scales, boundaries, pdf_path, "pdf",
    palette, hdr_colors, n_pixels, angle_from_mu, back_alpha,
    rear_density_alpha
  )
  comets_ub_render_heatmaps(
    samples, thetas, scales, boundaries, png_path, "png",
    palette, hdr_colors, n_pixels, angle_from_mu, back_alpha,
    rear_density_alpha
  )

  theta_df <- do.call(rbind, lapply(names(thetas), function(dataset) {
    theta <- thetas[[dataset]]
    data.frame(
      dataset = dataset,
      n = nrow(samples[[dataset]]),
      mu_1 = theta$mu[[1L]],
      mu_2 = theta$mu[[2L]],
      mu_3 = theta$mu[[3L]],
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      parameter_file = parameter_paths[[dataset]],
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(theta_df, file.path(output_dir, "parameters_used.csv"), row.names = FALSE)
  scale_df <- do.call(rbind, lapply(names(scales), function(dataset) {
    scale <- scales[[dataset]]
    data.frame(
      dataset = dataset,
      density_lower = scale$lower,
      density_upper = scale$upper,
      log10_density_lower = scale$log_lower,
      log10_density_upper = scale$log_upper,
      reference_size = scale$reference_size,
      stringsAsFactors = FALSE
    )
  }))
  rownames(scale_df) <- NULL
  utils::write.csv(scale_df, file.path(output_dir, "density_color_scales.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    fit_root = fit_root,
    sphere_coverage = paste(
      "single transparent orthographic projection per sample; front and rear",
      "surface densities superposed, with the rear surface faded"
    ),
    camera_reference = "short-period fitted axis",
    angle_from_short_fitted_axis_degrees = as.numeric(angle_from_mu),
    back_alpha = as.numeric(back_alpha),
    rear_density_alpha = as.numeric(rear_density_alpha),
    heatmap_quantity = paste(
      "fitted UB surface density f(x) on both projected sphere surfaces;",
      "rear surface faded"
    ),
    heatmap_scaling = paste(
      "log10 density; separate full finite-grid range for each fitted model;",
      "no HDR-based clipping"
    ),
    overlaid_hdr_probabilities = paste(hdr_probabilities, collapse = ","),
    overlaid_hdr_colors = paste(hdr_colors, collapse = ","),
    style_reference = "paper Risoe density panels and sunspot sphere panels",
    nonparametric_estimator = "none",
    n_pixels_per_view = as.integer(n_pixels),
    palette = "#FFFFFF,#E8F3F8,#B9DDE7,#72B7C5,#2A788E,#234B70",
    stringsAsFactors = FALSE
  ), file.path(output_dir, "plot_metadata.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  invisible(list(pdf = pdf_path, png = png_path, parameters = theta_df))
}

parse_comets_ub_heatmap_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --fit_root=PATH --output_dir=PATH --n_pixels=INTEGER ",
        "--angle_from_mu=NUMBER --back_alpha=NUMBER ",
        "--rear_density_alpha=NUMBER ",
        "--hdr_probabilities=0.5,0.8,0.95\n"
      ))
      quit(save = "no", status = 0L)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L) stop("Invalid option: ", arg)
    key <- pieces[[1L]]
    value <- pieces[[2L]]
    if (key %in% c("fit_root", "output_dir")) out[[key]] <- value
    if (identical(key, "n_pixels")) out[[key]] <- as.integer(value)
    if (key %in% c(
      "angle_from_mu", "back_alpha", "rear_density_alpha"
    )) {
      out[[key]] <- as.numeric(value)
    }
    if (identical(key, "hdr_probabilities")) {
      out[[key]] <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
    }
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(
    run_comets_uniform_beta_spherical_heatmaps,
    parse_comets_ub_heatmap_args()
  )
}
