#!/usr/bin/env Rscript

# Filled parametric highest-density regions on S^2, following the visual idea
# of ks::plot.kde(..., display = "filled.contour", cont = c(25, 50, 75)).
# The bands use fitted parametric densities only. No KDE is plotted or fitted.

resolve_spherical_hdr_band_path <- function(...) {
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

source(resolve_spherical_hdr_band_path(
  "scripts", "plot_comets_uniform_beta_spherical_contours.R"
))

spherical_hdr_band_validate <- function(levels, alphas, rear_factor) {
  levels <- sort(unique(as.numeric(levels)))
  alphas <- as.numeric(alphas)
  rear_factor <- as.numeric(rear_factor)
  if (!identical(levels, c(0.25, 0.50, 0.75))) {
    stop("The ks-style bands require probability contents 0.25, 0.50, and 0.75.")
  }
  if (length(alphas) != length(levels) || any(!is.finite(alphas)) ||
      any(alphas <= 0 | alphas > 1) || any(diff(alphas) < 0)) {
    stop("`alphas` must contain three nondecreasing values in (0, 1].")
  }
  if (length(rear_factor) != 1L || !is.finite(rear_factor) ||
      rear_factor <= 0 || rear_factor >= 1) {
    stop("`rear_factor` must lie strictly inside (0, 1).")
  }
  list(levels = levels, alphas = alphas, rear_factor = rear_factor)
}

spherical_hdr_band_raster <- function(
    density_fun, thresholds, camera, surface = c("front", "rear"),
    base_color = c("#E2E6BD", "#EEAB65", "#8E063B"), alphas = rep(1, 3),
    rear_factor = 0.35, n_pixels = 701L) {
  surface <- match.arg(surface)
  thresholds <- as.numeric(thresholds)
  n_pixels <- as.integer(n_pixels)
  if (length(thresholds) != 3L || any(!is.finite(thresholds)) ||
      is.unsorted(rev(thresholds), strictly = TRUE)) {
    stop("Thresholds must be finite and strictly decreasing for 25%, 50%, 75% HDRs.")
  }
  if (!is.finite(n_pixels) || n_pixels < 201L) {
    stop("`n_pixels` must be an integer of at least 201.")
  }

  x_seq <- seq(-1, 1, length.out = n_pixels)
  y_seq <- seq(1, -1, length.out = n_pixels)
  plane <- expand.grid(x = x_seq, y = y_seq)
  radius_squared <- plane$x^2 + plane$y^2
  inside <- radius_squared <= 1
  depth <- sqrt(pmax(0, 1 - radius_squared[inside]))
  if (identical(surface, "rear")) depth <- -depth

  xyz <- tcrossprod(plane$x[inside], camera$right) +
    tcrossprod(plane$y[inside], camera$up) +
    tcrossprod(depth, camera$view)
  density <- as.numeric(density_fun(xyz))
  if (length(density) != nrow(xyz) || any(!is.finite(density))) {
    stop("The density function returned non-finite or incorrectly sized values.")
  }

  # thresholds[1:3] correspond respectively to the 25%, 50%, and 75% HDRs.
  band <- integer(length(density))
  band[density >= thresholds[[3L]]] <- 1L
  band[density >= thresholds[[2L]]] <- 2L
  band[density >= thresholds[[1L]]] <- 3L
  surface_alphas <- alphas * if (identical(surface, "rear")) rear_factor else 1
  band_colors <- c(
    "#00000000",
    vapply(
      seq_along(surface_alphas),
      function(i) grDevices::adjustcolor(base_color[[i]], alpha.f = surface_alphas[[i]]),
      character(1L)
    )
  )
  colors <- rep("#00000000", nrow(plane))
  colors[inside] <- band_colors[band + 1L]
  as.raster(matrix(colors, nrow = n_pixels, ncol = n_pixels, byrow = TRUE))
}

spherical_hdr_band_draw_panel <- function(
    x, density_fun, thresholds, camera,
    base_color = c("#E2E6BD", "#EEAB65", "#8E063B"), alphas = rep(1, 3),
    rear_factor = 0.35, back_alpha = 0.18, n_pixels = 701L) {
  rear_raster <- spherical_hdr_band_raster(
    density_fun, thresholds, camera, surface = "rear",
    base_color = base_color, alphas = alphas,
    rear_factor = rear_factor, n_pixels = n_pixels
  )
  front_raster <- spherical_hdr_band_raster(
    density_fun, thresholds, camera, surface = "front",
    base_color = base_color, alphas = alphas,
    rear_factor = rear_factor, n_pixels = n_pixels
  )

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08),
    asp = 1, xaxs = "i", yaxs = "i"
  )
  graphics::rasterImage(rear_raster, -1, -1, 1, 1, interpolate = FALSE)
  graphics::rasterImage(front_raster, -1, -1, 1, 1, interpolate = FALSE)
  sunspots_joint_draw_sphere_outline(camera, back_alpha = back_alpha)
  for (front in c(FALSE, TRUE)) {
    sunspots_joint_draw_window_scatter(
      x_window = x,
      camera = camera,
      front = front,
      back_alpha = back_alpha
    )
  }
  invisible(NULL)
}

spherical_hdr_band_render_one <- function(
    x, density_fun, thresholds, camera, output_path, device,
    base_color, alphas, rear_factor, back_alpha, n_pixels) {
  if (identical(device, "pdf")) {
    grDevices::pdf(
      output_path, width = 3.2, height = 3.2,
      onefile = FALSE, paper = "special", bg = "transparent",
      version = "1.4", useDingbats = FALSE
    )
  } else if (identical(device, "png")) {
    grDevices::png(
      output_path, width = 960, height = 960, res = 300, bg = "white"
    )
  } else {
    stop("Unknown graphics device: ", device)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  tryCatch(
    {
      graphics::par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA)
      spherical_hdr_band_draw_panel(
        x = x, density_fun = density_fun, thresholds = thresholds,
        camera = camera, base_color = base_color, alphas = alphas,
        rear_factor = rear_factor, back_alpha = back_alpha,
        n_pixels = n_pixels
      )
    },
    finally = {
      graphics::par(old_par)
      grDevices::dev.off()
    }
  )
  invisible(output_path)
}

run_sunspot_parametric_hdr_bands <- function(
    source_dir = file.path(
      "real_data", "sunspots", "output",
      "cycle23_joint_spatial_window_parametric_shared_centered_daily"
    ),
    output_dir = file.path(source_dir, "individual_spheres_hdr_bands"),
    levels = c(0.25, 0.50, 0.75),
    base_color = c("#E2E6BD", "#EEAB65", "#8E063B"),
    alphas = rep(1, 3),
    rear_factor = 0.35,
    back_alpha = 0.18,
    view_theta = 35,
    view_phi = 18,
    n_pixels = 701L,
    threshold_reference_size = 200001L) {
  style <- spherical_hdr_band_validate(levels, alphas, rear_factor)
  levels <- style$levels
  alphas <- style$alphas
  rear_factor <- style$rear_factor
  metadata <- utils::read.csv(
    file.path(source_dir, "cycle23_joint_spatial_window_metadata.csv"),
    stringsAsFactors = FALSE
  )
  if (nrow(metadata) != 1L) {
    stop("The active sunspot plotting metadata are malformed.")
  }
  loaded <- sunspots_joint_load_window_plot_inputs(
    input_csv = metadata$input_csv[[1L]],
    start_date = metadata$start_date[[1L]],
    end_date = metadata$end_date_exclusive[[1L]],
    dequantization_seed = metadata$dequantization_seed[[1L]],
    hemisphere_regression = metadata$hemisphere_regression[[1L]],
    control = list(),
    from_joint_output_dir = metadata$from_joint_output_dir[[1L]]
  )
  retained <- loaded$retained
  theta_hat <- loaded$fit$theta_hat
  x <- jp_normalize_unit_matrix(
    as.matrix(retained[, c("x1", "x2", "x3")]),
    arg_name = "`retained`", min_ncol = 3L
  )
  windows <- sunspots_joint_spatial_rank_windows(
    retained$s,
    lower_levels = seq(0, 0.8, by = 0.2),
    upper_levels = seq(0.2, 1, by = 0.2),
    center_levels = seq(0.1, 0.9, by = 0.2)
  )
  definitions <- windows$summary

  camera <- sunspots_joint_sphere_camera(theta = view_theta, phi = view_phi)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  threshold_rows <- vector("list", nrow(definitions))
  pdf_paths <- png_paths <- character(nrow(definitions))

  for (i in seq_len(nrow(definitions))) {
    center_s <- definitions$center_s[[i]]
    density_fun <- local({
      center_s_i <- center_s
      theta_i <- theta_hat
      function(x_eval) sunspots_joint_parametric_conditional_density(
        x_eval, center_s = center_s_i, theta_hat = theta_i
      )
    })
    z_reference <- -1 +
      (2 * seq_len(as.integer(threshold_reference_size)) - 1) /
        as.integer(threshold_reference_size)
    x_reference <- cbind(
      sqrt(pmax(0, 1 - z_reference^2)),
      0,
      z_reference
    )
    density_reference <- density_fun(x_reference)
    reference_weights <- rep(
      1 / length(density_reference), length(density_reference)
    )
    hdr <- sunspots_joint_hdr_thresholds(
      density_reference,
      area_weights = reference_weights,
      levels = levels
    )
    total_reference_mass <- sum(density_reference * reference_weights)
    probability_content <- vapply(
      hdr$threshold,
      function(threshold) {
        sum(
          density_reference[density_reference >= threshold] *
            reference_weights[density_reference >= threshold]
        ) / total_reference_mass
      },
      numeric(1L)
    )
    threshold_rows[[i]] <- data.frame(
      dataset = "sunspots", window_id = i,
      probability = hdr$level, density_threshold = hdr$threshold,
      probability_content = probability_content,
      stringsAsFactors = FALSE
    )
    x_window <- x[windows$windows[[i]]$index, , drop = FALSE]
    pdf_paths[[i]] <- file.path(
      output_dir, sprintf("cycle23_joint_spatial_window_%02d_hdr_bands.pdf", i)
    )
    png_paths[[i]] <- sub("[.]pdf$", ".png", pdf_paths[[i]])
    spherical_hdr_band_render_one(
      x_window, density_fun, hdr$threshold, camera, pdf_paths[[i]], "pdf",
      base_color, alphas, rear_factor, back_alpha, n_pixels
    )
    spherical_hdr_band_render_one(
      x_window, density_fun, hdr$threshold, camera, png_paths[[i]], "png",
      base_color, alphas, rear_factor, back_alpha, n_pixels
    )
  }
  names(pdf_paths) <- names(png_paths) <- sprintf("window_%02d", seq_len(5L))
  utils::write.csv(
    definitions, file.path(output_dir, "window_definitions.csv"), row.names = FALSE
  )
  utils::write.csv(
    do.call(rbind, threshold_rows),
    file.path(output_dir, "hdr_band_thresholds.csv"), row.names = FALSE
  )
  utils::write.csv(data.frame(
    source_dir = source_dir,
    source_fit = metadata$from_joint_output_dir[[1L]],
    hdr_probability_contents = paste(levels, collapse = ","),
    interpretation = "nested parametric highest-density regions",
    band_alpha_outer_to_inner = paste(alphas, collapse = ","),
    rear_surface_alpha_factor = rear_factor,
    base_color = paste(base_color, collapse = ","),
    boundary_lines = "none",
    view_theta = view_theta,
    view_phi = view_phi,
    back_alpha = back_alpha,
    n_pixels = as.integer(n_pixels),
    threshold_reference_size = as.integer(threshold_reference_size),
    threshold_integration = "uniform midpoint rule in z on [-1,1]",
    style_reference = "ks filled.contour idea adapted to S2",
    stringsAsFactors = FALSE
  ), file.path(output_dir, "plot_metadata.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(list(pdf = pdf_paths, png = png_paths))
}

run_comet_parametric_hdr_bands <- function(
    fit_root = file.path(
      "real_data", "reruns", "paper_main_realdata_B1000_3cores_20260831_113532",
      "comets", "uniform_beta"
    ),
    output_dir = file.path(
      "real_data", "comets", "mixture", "uniform_beta_spherical_hdr_bands"
    ),
    levels = c(0.25, 0.50, 0.75),
    base_color = c("#E2E6BD", "#EEAB65", "#8E063B"),
    alphas = rep(1, 3),
    rear_factor = 0.35,
    back_alpha = 0.18,
    view_theta = 35,
    view_phi = 18,
    n_pixels = 701L) {
  style <- spherical_hdr_band_validate(levels, alphas, rear_factor)
  levels <- style$levels
  alphas <- style$alphas
  rear_factor <- style$rear_factor
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
  camera <- sunspots_joint_sphere_camera(theta = view_theta, phi = view_phi)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  threshold_rows <- list()
  pdf_paths <- png_paths <- character(length(samples))
  names(pdf_paths) <- names(png_paths) <- names(samples)

  for (dataset in names(samples)) {
    theta <- thetas[[dataset]]
    hdr <- do.call(rbind, lapply(
      levels, comets_ub_hdr_boundaries_one, theta = theta
    ))
    thresholds <- vapply(
      levels,
      function(level) unique(hdr$density_threshold[hdr$probability == level]),
      numeric(1L)
    )
    threshold_rows[[dataset]] <- data.frame(
      dataset = paste0(dataset, "_period"),
      probability = levels,
      density_threshold = thresholds,
      probability_content = vapply(
        levels,
        function(level) unique(hdr$probability_content[hdr$probability == level]),
        numeric(1L)
      ),
      stringsAsFactors = FALSE
    )
    density_fun <- local({
      theta_i <- theta
      function(x_eval) d_sph_uniform_beta_mixture_s2(
        x_eval, theta_i$mu, theta_i$weight_uniform,
        theta_i$alpha, theta_i$beta
      )
    })
    pdf_paths[[dataset]] <- file.path(
      output_dir,
      sprintf("comets_uniform_beta_%s_period_hdr_bands.pdf", dataset)
    )
    png_paths[[dataset]] <- sub("[.]pdf$", ".png", pdf_paths[[dataset]])
    spherical_hdr_band_render_one(
      samples[[dataset]], density_fun, thresholds, camera,
      pdf_paths[[dataset]], "pdf", base_color, alphas,
      rear_factor, back_alpha, n_pixels
    )
    spherical_hdr_band_render_one(
      samples[[dataset]], density_fun, thresholds, camera,
      png_paths[[dataset]], "png", base_color, alphas,
      rear_factor, back_alpha, n_pixels
    )
  }
  utils::write.csv(
    do.call(rbind, threshold_rows),
    file.path(output_dir, "hdr_band_thresholds.csv"), row.names = FALSE
  )
  utils::write.csv(data.frame(
    fit_root = fit_root,
    parameter_files = paste(unlist(parameter_paths), collapse = " | "),
    hdr_probability_contents = paste(levels, collapse = ","),
    interpretation = "nested parametric highest-density regions",
    band_alpha_outer_to_inner = paste(alphas, collapse = ","),
    rear_surface_alpha_factor = rear_factor,
    base_color = paste(base_color, collapse = ","),
    boundary_lines = "none",
    view_theta = view_theta,
    view_phi = view_phi,
    back_alpha = back_alpha,
    n_pixels = as.integer(n_pixels),
    style_reference = "ks filled.contour idea adapted to S2",
    stringsAsFactors = FALSE
  ), file.path(output_dir, "plot_metadata.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(list(pdf = pdf_paths, png = png_paths))
}

run_all_spherical_parametric_hdr_bands <- function() {
  list(
    sunspots = run_sunspot_parametric_hdr_bands(),
    comets = run_comet_parametric_hdr_bands()
  )
}

if (sys.nframe() == 0L) {
  run_all_spherical_parametric_hdr_bands()
}
