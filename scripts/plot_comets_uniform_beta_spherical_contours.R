#!/usr/bin/env Rscript

# Plot the comet orbital-pole samples on S^2 together with highest-density
# contours from the fitted uniform-beta (UB) models. The fitted parameters are
# read from the saved paper run; this script does not refit either model.

resolve_comets_ub_contour_path <- function(...) {
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

source(resolve_comets_ub_contour_path("utils.R"))
source(resolve_comets_ub_contour_path("real_data", "comets", "utils_comets_data.R"))
source(resolve_comets_ub_contour_path(
  "real_data", "sunspots", "run_sunspots_cycle23_joint_spatial_window_kde_plots.R"
))

comets_ub_read_theta <- function(path) {
  if (!file.exists(path)) stop("Saved UB parameter file not found: ", path)
  theta_df <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("mu_1", "mu_2", "mu_3", "weight_uniform", "alpha", "beta")
  if (nrow(theta_df) != 1L || !all(required %in% names(theta_df))) {
    stop("Saved UB parameter file must have one row and columns: ", paste(required, collapse = ", "))
  }
  theta <- list(
    mu = as.numeric(theta_df[1L, c("mu_1", "mu_2", "mu_3")]),
    weight_uniform = as.numeric(theta_df$weight_uniform[[1L]]),
    alpha = as.numeric(theta_df$alpha[[1L]]),
    beta = as.numeric(theta_df$beta[[1L]])
  )
  uniform_beta_mixture_normalize_theta(theta, ambient_dim = 3L)
}

comets_ub_cdf_y <- function(y, theta) {
  uniform_beta_mixture_cdf_y(
    y = y,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
}

comets_ub_density_y <- function(y, theta) {
  uniform_beta_mixture_density_y(
    y = y,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
}

comets_ub_quantile_y <- function(probability, theta, tol = 1e-12) {
  probability <- as.numeric(probability)
  if (probability <= 0) return(0)
  if (probability >= 1) return(1)
  stats::uniroot(
    function(y) comets_ub_cdf_y(y, theta) - probability,
    interval = c(0, 1),
    tol = tol
  )$root
}

comets_ub_hdr_boundaries_one <- function(probability, theta, tol = 1e-11) {
  probability <- as.numeric(probability)
  if (length(probability) != 1L || !is.finite(probability) ||
      probability <= 0 || probability >= 1) {
    stop("HDR probabilities must lie strictly inside (0, 1).")
  }

  alpha <- theta$alpha
  beta <- theta$beta
  mixture_cdf <- function(y) comets_ub_cdf_y(y, theta)
  log_beta_density <- function(y) stats::dbeta(y, alpha, beta, log = TRUE)
  root_eps <- 1e-8
  endpoint_eps <- .Machine$double.eps

  if (alpha >= 1 && beta <= 1 && !(alpha == 1 && beta == 1)) {
    y_boundary <- comets_ub_quantile_y(1 - probability, theta, tol = tol)
    boundaries <- y_boundary
    region <- sprintf("[%.16g,1]", y_boundary)
    mass <- 1 - mixture_cdf(y_boundary)
  } else if (alpha <= 1 && beta >= 1 && !(alpha == 1 && beta == 1)) {
    y_boundary <- comets_ub_quantile_y(probability, theta, tol = tol)
    boundaries <- y_boundary
    region <- sprintf("[0,%.16g]", y_boundary)
    mass <- mixture_cdf(y_boundary)
  } else if (alpha < 1 && beta < 1) {
    y_minimum <- (alpha - 1) / (alpha + beta - 2)

    right_boundary <- function(y_left) {
      target <- log_beta_density(y_left)
      minimum <- log_beta_density(y_minimum)
      if (target <= minimum + 1e-12) return(y_minimum)
      stats::uniroot(
        function(y) log_beta_density(y) - target,
        interval = c(y_minimum, 1 - endpoint_eps),
        tol = tol
      )$root
    }

    hdr_mass <- function(y_left) {
      y_right <- right_boundary(y_left)
      mixture_cdf(y_left) + 1 - mixture_cdf(y_right)
    }

    y_left <- stats::uniroot(
      function(y) hdr_mass(y) - probability,
      interval = c(root_eps, y_minimum - root_eps),
      tol = tol
    )$root
    y_right <- right_boundary(y_left)
    boundaries <- c(y_left, y_right)
    region <- sprintf("[0,%.16g] union [%.16g,1]", y_left, y_right)
    mass <- mixture_cdf(y_left) + 1 - mixture_cdf(y_right)
  } else if (alpha > 1 && beta > 1) {
    y_mode <- (alpha - 1) / (alpha + beta - 2)

    right_boundary <- function(y_left) {
      target <- log_beta_density(y_left)
      maximum <- log_beta_density(y_mode)
      if (target >= maximum - 1e-12) return(y_mode)
      stats::uniroot(
        function(y) log_beta_density(y) - target,
        interval = c(y_mode, 1 - root_eps),
        tol = tol
      )$root
    }

    hdr_mass <- function(y_left) {
      y_right <- right_boundary(y_left)
      mixture_cdf(y_right) - mixture_cdf(y_left)
    }

    y_left <- stats::uniroot(
      function(y) hdr_mass(y) - probability,
      interval = c(root_eps, y_mode - root_eps),
      tol = tol
    )$root
    y_right <- right_boundary(y_left)
    boundaries <- c(y_left, y_right)
    region <- sprintf("[%.16g,%.16g]", y_left, y_right)
    mass <- mixture_cdf(y_right) - mixture_cdf(y_left)
  } else {
    stop("The uniform density has no unique nontrivial HDR contours.")
  }

  density_thresholds <- comets_ub_density_y(boundaries, theta) / (4 * pi)
  if (max(density_thresholds) - min(density_thresholds) >
      1e-7 * max(1, max(abs(density_thresholds)))) {
    stop("Numerically inconsistent UB HDR boundary densities.")
  }
  if (abs(mass - probability) > 1e-7) {
    stop("Numerically inconsistent UB HDR probability content.")
  }

  data.frame(
    probability = probability,
    boundary_id = seq_along(boundaries),
    y = boundaries,
    z = 2 * boundaries - 1,
    density_threshold = mean(density_thresholds),
    probability_content = mass,
    hdr_region_y = region,
    stringsAsFactors = FALSE
  )
}

comets_ub_hdr_boundaries <- function(theta, probabilities) {
  do.call(rbind, lapply(
    probabilities,
    comets_ub_hdr_boundaries_one,
    theta = theta
  ))
}

comets_ub_cross3 <- function(a, b) {
  c(
    a[[2L]] * b[[3L]] - a[[3L]] * b[[2L]],
    a[[3L]] * b[[1L]] - a[[1L]] * b[[3L]],
    a[[1L]] * b[[2L]] - a[[2L]] * b[[1L]]
  )
}

comets_ub_normalize_vector <- function(x) {
  x <- as.numeric(x)
  norm_x <- sqrt(sum(x^2))
  if (!is.finite(norm_x) || norm_x <= 0) stop("Cannot normalize this vector.")
  x / norm_x
}

comets_ub_circle_xyz <- function(mu, z, n_points = 721L) {
  mu <- comets_ub_normalize_vector(mu)
  z <- as.numeric(z)
  if (length(z) != 1L || !is.finite(z) || z < -1 || z > 1) {
    stop("A contour coordinate z must lie in [-1, 1].")
  }
  basis <- jp_orthonormal_complement(mu)
  angle <- seq(0, 2 * pi, length.out = as.integer(n_points))
  tangent <- tcrossprod(cos(angle), basis[, 1L]) +
    tcrossprod(sin(angle), basis[, 2L])
  tcrossprod(rep(z, length(angle)), mu) +
    sqrt(max(0, 1 - z^2)) * tangent
}

comets_ub_camera <- function(mu, angle_from_mu = 35) {
  mu <- comets_ub_normalize_vector(mu)
  basis <- jp_orthonormal_complement(mu)
  angle <- as.numeric(angle_from_mu) * pi / 180
  view <- comets_ub_normalize_vector(cos(angle) * mu + sin(angle) * basis[, 1L])
  up <- comets_ub_normalize_vector(mu - sum(mu * view) * view)
  right <- comets_ub_normalize_vector(comets_ub_cross3(up, view))
  list(view = view, right = right, up = up)
}

comets_ub_project <- function(x, camera) {
  x <- as.matrix(x)
  if (ncol(x) != 3L || any(!is.finite(x))) {
    stop("Points to project must form a finite matrix with three columns.")
  }
  data.frame(
    x = drop(x %*% camera$right),
    y = drop(x %*% camera$up),
    depth = drop(x %*% camera$view),
    stringsAsFactors = FALSE
  )
}

comets_ub_visibility_runs <- function(depth, front) {
  visible <- if (isTRUE(front)) depth >= 0 else depth < 0
  if (length(visible) < 2L) return(list())
  runs <- rle(visible)
  run_end <- cumsum(runs$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  keep <- which(runs$values & runs$lengths >= 2L)
  lapply(keep, function(i) run_start[[i]]:run_end[[i]])
}

comets_ub_draw_curve <- function(xyz, camera, color, lwd, front, alpha) {
  projected <- comets_ub_project(xyz, camera)
  color <- grDevices::adjustcolor(color, alpha.f = alpha)
  for (idx in comets_ub_visibility_runs(projected$depth, front = front)) {
    graphics::lines(projected$x[idx], projected$y[idx], col = color, lwd = lwd)
  }
  invisible(NULL)
}

comets_ub_lon_lat_xyz <- function(lon_deg, lat_deg) {
  lon <- as.numeric(lon_deg) * pi / 180
  lat <- as.numeric(lat_deg) * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

comets_ub_draw_graticule <- function(camera, front, alpha) {
  lon_dense <- seq(-180, 180, length.out = 721L)
  lat_dense <- seq(-90, 90, length.out = 361L)
  for (latitude in c(-60, -30, 0, 30, 60)) {
    comets_ub_draw_curve(
      comets_ub_lon_lat_xyz(lon_dense, rep(latitude, length(lon_dense))),
      camera, color = "#666666", lwd = 0.35, front = front, alpha = alpha
    )
  }
  for (longitude in seq(-180, 120, by = 60)) {
    comets_ub_draw_curve(
      comets_ub_lon_lat_xyz(rep(longitude, length(lat_dense)), lat_dense),
      camera, color = "#666666", lwd = 0.35, front = front, alpha = alpha
    )
  }
  invisible(NULL)
}

comets_ub_draw_points <- function(x, camera, front, back_alpha = 0.13) {
  projected <- comets_ub_project(x, camera)
  keep <- if (isTRUE(front)) projected$depth >= 0 else projected$depth < 0
  if (any(keep)) {
    graphics::points(
      projected$x[keep], projected$y[keep],
      pch = 16, cex = 0.30,
      col = if (isTRUE(front)) "#111111" else
        grDevices::adjustcolor("#111111", alpha.f = back_alpha)
    )
  }
  invisible(NULL)
}

comets_ub_draw_panel <- function(x, theta, boundaries, title,
                                 colors, camera,
                                 back_alpha = 0.18) {
  curves <- lapply(seq_len(nrow(boundaries)), function(i) {
    list(
      level = boundaries$probability[[i]],
      xyz = comets_ub_circle_xyz(theta$mu, boundaries$z[[i]])
    )
  })
  probabilities <- sort(unique(boundaries$probability))

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08),
    asp = 1, xaxs = "i", yaxs = "i"
  )
  sunspots_joint_draw_sphere_outline(camera, back_alpha = back_alpha)
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
      colors = colors,
      lty = 1,
      lwd = 0.75,
      front = front,
      back_alpha = back_alpha,
      label_specs = label_specs
    )
    sunspots_joint_draw_window_scatter(
      x_window = x,
      camera = camera,
      front = front,
      back_alpha = back_alpha
    )
  }
  sunspots_joint_draw_contour_labels(
    label_specs = label_specs,
    colors = colors,
    cex = 0.54
  )
  if (!is.null(title) && nzchar(title)) {
    graphics::title(main = title, cex.main = 0.78, line = 0.25)
  }
  invisible(NULL)
}

comets_ub_render_individual_sunspot_style <- function(
    samples, thetas, boundaries, output_dir, device,
    colors, view_theta, view_phi, back_alpha) {
  common_camera <- sunspots_joint_sphere_camera(
    theta = view_theta,
    phi = view_phi
  )
  output_paths <- character(length(samples))
  names(output_paths) <- names(samples)

  for (dataset in names(samples)) {
    output_paths[[dataset]] <- file.path(
      output_dir,
      sprintf("comets_uniform_beta_%s_period_contours.%s", dataset, device)
    )
    if (identical(device, "pdf")) {
      grDevices::pdf(
        file = output_paths[[dataset]],
        width = 3.2,
        height = 3.2,
        onefile = FALSE,
        paper = "special",
        bg = "transparent",
        version = "1.4",
        useDingbats = FALSE
      )
    } else if (identical(device, "png")) {
      grDevices::png(
        filename = output_paths[[dataset]],
        width = 960,
        height = 960,
        res = 300,
        bg = "white"
      )
    } else {
      stop("Unknown graphics device: ", device)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    tryCatch(
      {
        graphics::par(
          mar = c(0, 0, 0, 0),
          oma = c(0, 0, 0, 0),
          xpd = NA
        )
        comets_ub_draw_panel(
          samples[[dataset]],
          thetas[[dataset]],
          boundaries[[dataset]],
          title = NULL,
          colors = colors,
          camera = common_camera,
          back_alpha = back_alpha
        )
      },
      finally = {
        graphics::par(old_par)
        grDevices::dev.off()
      }
    )
  }
  invisible(output_paths)
}

comets_ub_draw_legend <- function(probabilities, colors) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend_labels <- c(
    sprintf("%.0f%% HDR", 100 * probabilities),
    "Observed orbital poles"
  )
  graphics::legend(
    "center",
    legend = legend_labels,
    col = c(colors, "#111111"),
    lty = c(rep(1, length(probabilities)), NA),
    lwd = c(rep(1.35, length(probabilities)), NA),
    pch = c(rep(NA_integer_, length(probabilities)), 16),
    pt.cex = 0.7,
    horiz = TRUE,
    bty = "n",
    cex = 0.76,
    seg.len = 1.8,
    x.intersp = 0.55
  )
  invisible(NULL)
}

comets_ub_render_combined <- function(samples, thetas, boundaries,
                                      output_path, device,
                                      probabilities, colors,
                                      angle_from_mu, back_alpha) {
  if (identical(device, "pdf")) {
    grDevices::pdf(
      output_path, width = 7.2, height = 3.85,
      onefile = FALSE, paper = "special", bg = "transparent",
      version = "1.4", useDingbats = FALSE
    )
  } else if (identical(device, "png")) {
    grDevices::png(output_path, width = 2300, height = 1230, res = 300, bg = "white")
  } else {
    stop("Unknown graphics device: ", device)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::layout(
    matrix(c(1L, 2L, 3L, 3L), nrow = 2L, byrow = TRUE),
    heights = c(1, 0.12)
  )
  common_camera <- comets_ub_camera(
    thetas$short$mu,
    angle_from_mu = angle_from_mu
  )
  graphics::par(mar = c(0, 0.25, 1.5, 0.25), oma = c(0, 0, 0, 0))
  comets_ub_draw_panel(
    samples$short, thetas$short, boundaries$short,
    sprintf("Short-period comets (n = %d)", nrow(samples$short)),
    colors, common_camera, back_alpha
  )
  comets_ub_draw_panel(
    samples$long, thetas$long, boundaries$long,
    sprintf("Long-period comets (n = %d)", nrow(samples$long)),
    colors, common_camera, back_alpha
  )
  graphics::par(mar = c(0, 0, 0, 0))
  comets_ub_draw_legend(probabilities, colors)
  invisible(output_path)
}

run_comets_uniform_beta_spherical_contours <- function(
    fit_root = file.path(
      "real_data", "reruns", "paper_main_realdata_B1000_3cores_20260831_113532",
      "comets", "uniform_beta"
    ),
    output_dir = file.path(
      "real_data", "comets", "mixture", "uniform_beta_spherical_contours"
    ),
    hdr_probabilities = c(0.50, 0.80, 0.95),
    view_theta = 35,
    view_phi = 18,
    back_alpha = 0.18) {

  hdr_probabilities <- sort(unique(as.numeric(hdr_probabilities)))
  if (length(hdr_probabilities) < 1L || any(!is.finite(hdr_probabilities)) ||
      any(hdr_probabilities <= 0 | hdr_probabilities >= 1)) {
    stop("`hdr_probabilities` must contain values strictly inside (0, 1).")
  }
  if (!identical(hdr_probabilities, c(0.50, 0.80, 0.95))) {
    stop("The paper-matched contour style requires HDR probabilities 0.50, 0.80, and 0.95.")
  }
  colors <- c("#1f78b4", "#e31a1c", "#1B5E20")

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

  boundaries <- lapply(thetas, comets_ub_hdr_boundaries, probabilities = hdr_probabilities)
  boundary_df <- do.call(rbind, lapply(names(boundaries), function(dataset) {
    data.frame(dataset = dataset, boundaries[[dataset]], stringsAsFactors = FALSE)
  }))
  rownames(boundary_df) <- NULL

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_paths <- comets_ub_render_individual_sunspot_style(
    samples, thetas, boundaries, output_dir, "pdf",
    colors, view_theta, view_phi, back_alpha
  )
  png_paths <- comets_ub_render_individual_sunspot_style(
    samples, thetas, boundaries, output_dir, "png",
    colors, view_theta, view_phi, back_alpha
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
  utils::write.csv(boundary_df, file.path(output_dir, "hdr_contour_boundaries.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    fit_root = fit_root,
    hdr_probabilities = paste(hdr_probabilities, collapse = ","),
    contour_definition = "highest-density regions under the fitted UB surface density",
    contour_colors = paste(colors, collapse = ","),
    style_reference = paste(
      "exact paper sunspot individual sphere-panel format, including camera"
    ),
    output_geometry = "one 3.2 by 3.2 inch transparent PDF per sample",
    internal_title = "none, matching the active sunspot PDFs",
    internal_legend = "none, matching the active sunspot PDFs",
    contour_line_width = 0.75,
    contour_label_cex = 0.54,
    point_cex = 0.3375,
    plot_geometry = paste(
      "common paper sunspot orthographic camera; rear hemisphere shown with",
      "reduced opacity"
    ),
    nonparametric_estimator = "none",
    camera_reference = "paper sunspot sphere panels",
    view_theta_degrees = as.numeric(view_theta),
    view_phi_degrees = as.numeric(view_phi),
    back_alpha = as.numeric(back_alpha),
    stringsAsFactors = FALSE
  ), file.path(output_dir, "plot_metadata.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  invisible(list(
    pdf = pdf_paths,
    png = png_paths,
    parameters = theta_df,
    hdr_boundaries = boundary_df
  ))
}

parse_comets_ub_contour_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --fit_root=PATH --output_dir=PATH ",
        "--hdr_probabilities=0.5,0.8,0.95 ",
        "--view_theta=NUMBER --view_phi=NUMBER --back_alpha=NUMBER\n"
      ))
      quit(save = "no", status = 0L)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L) stop("Invalid option: ", arg)
    key <- pieces[[1L]]
    value <- pieces[[2L]]
    if (key %in% c("fit_root", "output_dir")) out[[key]] <- value
    if (identical(key, "hdr_probabilities")) {
      out[[key]] <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
    }
    if (key %in% c("view_theta", "view_phi", "back_alpha")) {
      out[[key]] <- as.numeric(value)
    }
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(
    run_comets_uniform_beta_spherical_contours,
    parse_comets_ub_contour_args()
  )
}
