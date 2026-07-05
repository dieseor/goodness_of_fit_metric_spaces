`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

parse_named_args <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
  if (length(args) == 0L) {
    return(list())
  }

  output <- vector("list", length(args))
  names(output) <- rep("", length(args))

  for (i in seq_along(args)) {
    arg <- substring(args[[i]], 3L)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[i]] <- value
    names(output)[[i]] <- key
  }

  output
}

parse_numeric_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) {
    return(as.numeric(default))
  }
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

format_beta_value <- function(beta) {
  format(beta, trim = TRUE, scientific = FALSE)
}

vmf_normalizing_constant_s2 <- function(kappa) {
  kappa <- as.numeric(kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a single nonnegative finite number.")
  }
  if (kappa == 0) {
    return(1 / (4 * pi))
  }
  kappa / (4 * pi * sinh(kappa))
}

build_sphere_mesh <- function(n_longitude = 180L, n_latitude = 120L) {
  n_longitude <- as.integer(n_longitude)
  n_latitude <- as.integer(n_latitude)
  if (n_longitude < 10L || n_latitude < 10L) {
    stop("`n_longitude` and `n_latitude` must both be at least 10.")
  }

  longitude <- seq(-pi, pi, length.out = n_longitude)
  latitude <- seq(-pi / 2, pi / 2, length.out = n_latitude)

  x <- outer(longitude, latitude, function(lon, lat) cos(lat) * cos(lon))
  y <- outer(longitude, latitude, function(lon, lat) cos(lat) * sin(lon))
  z <- outer(longitude, latitude, function(lon, lat) sin(lat))

  list(x = x, y = y, z = z)
}

vmf_antipodal_density_s2 <- function(x1, beta, kappa = 2) {
  density_constant <- vmf_normalizing_constant_s2(kappa = kappa)
  density_constant * ((1 - beta) * exp(kappa * x1) + beta * exp(-kappa * x1))
}

map_values_to_colors <- function(values,
                                 limits,
                                 palette_values,
                                 alpha_f = 0.92) {
  if (length(limits) != 2L || !all(is.finite(limits)) || limits[[1L]] >= limits[[2L]]) {
    stop("`limits` must be an increasing pair of finite numbers.")
  }

  scaled <- (values - limits[[1L]]) / (limits[[2L]] - limits[[1L]])
  scaled <- pmin(pmax(scaled, 0), 1)
  idx <- 1L + floor(scaled * (length(palette_values) - 1L))
  if (length(alpha_f) == 1L) {
    alpha_f <- rep(alpha_f, length(idx))
  }
  rgb_mat <- grDevices::col2rgb(palette_values[idx], alpha = FALSE)
  grDevices::rgb(
    red = rgb_mat[1L, ],
    green = rgb_mat[2L, ],
    blue = rgb_mat[3L, ],
    alpha = round(255 * alpha_f),
    maxColorValue = 255
  )
}

rotate_xyz <- function(X, azim, elev) {
  Rz <- matrix(
    c(cos(azim), -sin(azim), 0,
      sin(azim), cos(azim), 0,
      0, 0, 1),
    nrow = 3,
    byrow = TRUE
  )
  Rx <- matrix(
    c(1, 0, 0,
      0, cos(elev), -sin(elev),
      0, sin(elev), cos(elev)),
    nrow = 3,
    byrow = TRUE
  )
  R <- Rx %*% Rz
  X %*% t(R)
}

project_mesh_orthographic <- function(mesh, azim, elev) {
  xyz <- cbind(c(mesh$x), c(mesh$y), c(mesh$z))
  xyz_rot <- rotate_xyz(xyz, azim = azim, elev = elev)
  dims <- dim(mesh$x)

  list(
    u = matrix(xyz_rot[, 1L], nrow = dims[[1L]], ncol = dims[[2L]]),
    v = matrix(xyz_rot[, 3L], nrow = dims[[1L]], ncol = dims[[2L]]),
    depth = matrix(xyz_rot[, 2L], nrow = dims[[1L]], ncol = dims[[2L]])
  )
}

draw_segmented_visibility_line <- function(line_xyz,
                                           azim,
                                           elev,
                                           front_color,
                                           back_color,
                                           lwd = 0.8) {
  rotated <- rotate_xyz(line_xyz, azim = azim, elev = elev)
  u <- rotated[, 1L]
  v <- rotated[, 3L]
  depth <- rotated[, 2L]

  for (k in seq_len(nrow(rotated) - 1L)) {
    d_mid <- 0.5 * (depth[[k]] + depth[[k + 1L]])
    seg_col <- if (d_mid >= 0) front_color else back_color
    graphics::segments(
      x0 = u[[k]],
      y0 = v[[k]],
      x1 = u[[k + 1L]],
      y1 = v[[k + 1L]],
      col = seg_col,
      lwd = lwd
    )
  }
}

draw_sphere_graticule <- function(azim,
                                  elev,
                                  front_color = grDevices::adjustcolor("grey20", alpha.f = 0.45),
                                  back_color = grDevices::adjustcolor("grey20", alpha.f = 0.14)) {
  meridians <- seq(-135, 135, by = 45) * pi / 180
  parallels <- seq(-60, 60, by = 30) * pi / 180
  latitude_grid <- seq(-pi / 2, pi / 2, length.out = 240L)
  longitude_grid <- seq(-pi, pi, length.out = 320L)

  for (lon in meridians) {
    line_xyz <- cbind(
      cos(latitude_grid) * cos(lon),
      cos(latitude_grid) * sin(lon),
      sin(latitude_grid)
    )
    draw_segmented_visibility_line(
      line_xyz = line_xyz,
      azim = azim,
      elev = elev,
      front_color = front_color,
      back_color = back_color,
      lwd = 0.7
    )
  }

  for (lat in parallels) {
    line_xyz <- cbind(
      cos(lat) * cos(longitude_grid),
      cos(lat) * sin(longitude_grid),
      rep(sin(lat), length(longitude_grid))
    )
    draw_segmented_visibility_line(
      line_xyz = line_xyz,
      azim = azim,
      elev = elev,
      front_color = front_color,
      back_color = back_color,
      lwd = 0.7
    )
  }
}

draw_density_sphere_panel <- function(beta,
                                      mesh,
                                      kappa,
                                      density_limits,
                                      palette_values,
                                      azim = 0,
                                      elev = 0,
                                      show_title = FALSE) {
  projected <- project_mesh_orthographic(mesh, azim = azim, elev = elev)

  x_mid <- 0.25 * (
    mesh$x[-1L, -1L] +
      mesh$x[-nrow(mesh$x), -1L] +
      mesh$x[-1L, -ncol(mesh$x)] +
      mesh$x[-nrow(mesh$x), -ncol(mesh$x)]
  )
  depth_mid <- 0.25 * (
    projected$depth[-1L, -1L] +
      projected$depth[-nrow(projected$depth), -1L] +
      projected$depth[-1L, -ncol(projected$depth)] +
      projected$depth[-nrow(projected$depth), -ncol(projected$depth)]
  )
  density_mid <- vmf_antipodal_density_s2(x_mid, beta = beta, kappa = kappa)
  alpha_mid <- 0.22 + 0.73 * ((depth_mid + 1) / 2)^0.9
  col_mid <- matrix(
    map_values_to_colors(
      values = c(density_mid),
      limits = density_limits,
      palette_values = palette_values,
      alpha_f = c(alpha_mid)
    ),
    nrow = nrow(density_mid),
    ncol = ncol(density_mid)
  )

  order_idx <- order(c(depth_mid), decreasing = FALSE)
  n_i <- nrow(density_mid)
  ii <- ((order_idx - 1L) %% n_i) + 1L
  jj <- ((order_idx - 1L) %/% n_i) + 1L

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-1.05, 1.05),
    ylim = c(-1.05, 1.05),
    asp = 1,
    xaxs = "i",
    yaxs = "i"
  )

  for (k in seq_along(order_idx)) {
    i <- ii[[k]]
    j <- jj[[k]]
    graphics::polygon(
      x = c(
        projected$u[i, j],
        projected$u[i + 1L, j],
        projected$u[i + 1L, j + 1L],
        projected$u[i, j + 1L]
      ),
      y = c(
        projected$v[i, j],
        projected$v[i + 1L, j],
        projected$v[i + 1L, j + 1L],
        projected$v[i, j + 1L]
      ),
      col = col_mid[i, j],
      border = NA
    )
  }

  draw_sphere_graticule(azim = azim, elev = elev)
  t <- seq(0, 2 * pi, length.out = 720L)
  graphics::lines(cos(t), sin(t), col = "grey25", lwd = 1.1)

  if (isTRUE(show_title)) {
    graphics::title(main = sprintf("beta = %s", format_beta_value(beta)), line = 0.2, cex.main = 1.2)
  }
}

draw_horizontal_colorbar <- function(density_limits,
                                     palette_values,
                                     label = "Density",
                                     n_breaks = 4L) {
  x_breaks <- seq(density_limits[[1L]], density_limits[[2L]], length.out = as.integer(n_breaks))
  x_seq <- seq(density_limits[[1L]], density_limits[[2L]], length.out = length(palette_values) + 1L)

  graphics::plot.new()
  graphics::plot.window(
    xlim = density_limits,
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i"
  )

  for (i in seq_len(length(palette_values))) {
    graphics::rect(
      xleft = x_seq[[i]],
      ybottom = 0.25,
      xright = x_seq[[i + 1L]],
      ytop = 0.7,
      col = palette_values[[i]],
      border = NA
    )
  }

  graphics::rect(
    xleft = density_limits[[1L]],
    ybottom = 0.25,
    xright = density_limits[[2L]],
    ytop = 0.7,
    border = "grey40",
    lwd = 0.8
  )
  graphics::axis(
    1,
    at = x_breaks,
    labels = format(round(x_breaks, 2), nsmall = 2),
    lwd = 0,
    lwd.ticks = 0.8,
    cex.axis = 0.95,
    line = -0.3
  )
  graphics::mtext(label, side = 3, line = 0.2, cex = 1.05)
}

save_vmf_antipodal_density_figures <- function(output_dir = file.path(
                                                 "simulation_results",
                                                 "power_mixtures_pilot",
                                                 "plots",
                                                 "theoretical_densities"
                                               ),
                                               file_stem = "vmf_s2_antipodal_mixture_density",
                                               beta_values = c(0, 0.25, 0.5, 0.75, 1),
                                               kappa = 2,
                                               n_longitude = 180L,
                                               n_latitude = 120L,
                                               width = 4.2,
                                               height = 4.2,
                                               dpi = 300) {
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Package `viridisLite` is required.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  beta_values <- as.numeric(beta_values)
  if (length(beta_values) == 0L || any(!is.finite(beta_values)) ||
      any(beta_values < 0 | beta_values > 1)) {
    stop("`beta_values` must contain finite values in [0, 1].")
  }

  density_constant <- vmf_normalizing_constant_s2(kappa = kappa)
  density_limits <- c(density_constant * exp(-kappa), density_constant * exp(kappa))
  palette_values <- viridisLite::viridis(256L, option = "D", begin = 0.15, end = 0.9)
  mesh <- build_sphere_mesh(n_longitude = n_longitude, n_latitude = n_latitude)
  output_paths <- character(length(beta_values))

  for (i in seq_along(beta_values)) {
    beta <- beta_values[[i]]
    beta_slug <- gsub(".", "p", format_beta_value(beta), fixed = TRUE)
    output_path <- file.path(
      output_dir,
      sprintf("%s_beta_%s.png", file_stem, beta_slug)
    )

    grDevices::png(
      filename = output_path,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "transparent",
      type = "cairo"
    )
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(mar = c(0, 0, 0, 0), family = "sans")
    draw_density_sphere_panel(
      beta = beta,
      mesh = mesh,
      kappa = kappa,
      density_limits = density_limits,
      palette_values = palette_values,
      azim = 0,
      elev = 0,
      show_title = FALSE
    )
    graphics::par(old_par)
    grDevices::dev.off()

    output_paths[[i]] <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  }

  unname(output_paths)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))

  output_paths <- save_vmf_antipodal_density_figures(
    output_dir = as.character(args$output_dir %||% file.path(
      "simulation_results",
      "power_mixtures_pilot",
      "plots",
      "theoretical_densities"
    )),
    file_stem = as.character(args$file_stem %||% "vmf_s2_antipodal_mixture_density"),
    beta_values = parse_numeric_csv(args$beta_values, c(0, 0.25, 0.5, 0.75, 1)),
    kappa = as.numeric(args$kappa %||% 2),
    n_longitude = as.integer(args$n_longitude %||% 180L),
    n_latitude = as.integer(args$n_latitude %||% 120L),
    width = as.numeric(args$width %||% 4.2),
    height = as.numeric(args$height %||% 4.2),
    dpi = as.integer(args$dpi %||% 300L)
  )

  cat(paste(output_paths, collapse = "\n"), "\n")
}
