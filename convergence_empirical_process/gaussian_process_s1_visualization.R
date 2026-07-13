# ============================================================================
# S1 VISUALIZATION AND DIAGNOSTICS FOR THE LIMIT GAUSSIAN PROCESS
# ============================================================================

utils_path_s1 <- if (file.exists("utils.R")) {
  "utils.R"
} else if (file.exists(file.path("..", "utils.R"))) {
  file.path("..", "utils.R")
} else {
  stop("Could not find utils.R in current directory or parent directory.")
}
source(utils_path_s1)

vmf_gp_path_s1 <- if (file.exists(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))) {
  file.path("convergence_empirical_process", "gaussian_process_vmf.R")
} else if (file.exists("gaussian_process_vmf.R")) {
  "gaussian_process_vmf.R"
} else {
  stop("Could not find gaussian_process_vmf.R in the expected locations.")
}
source(vmf_gp_path_s1)


# ============================================================================
# INTERNAL HELPERS
# ============================================================================

make_theta_colors <- function(theta_values,
                              reference_theta = 0,
                              center_color = "#FDE725",
                              edge_color = "#2C7FB8") {
  theta_values <- wrap_angle_2pi(as.numeric(theta_values))
  reference_theta <- wrap_angle_2pi(as.numeric(reference_theta)[1L])
  signed_diff <- atan2(
    sin(theta_values - reference_theta),
    cos(theta_values - reference_theta)
  )
  normalized_distance <- abs(signed_diff) / pi
  palette_fun <- grDevices::colorRamp(c(center_color, edge_color), space = "Lab")
  rgb_values <- palette_fun(normalized_distance)
  grDevices::rgb(
    red = rgb_values[, 1],
    green = rgb_values[, 2],
    blue = rgb_values[, 3],
    maxColorValue = 255
  )
}

make_theta_rainbow_colors <- function(theta_values, start = 0, end = 5 / 6) {
  n_theta <- length(theta_values)
  if (n_theta < 1L) {
    return(character(0))
  }
  grDevices::rainbow(n_theta, start = start, end = end)
}

normalize_s1_color_scheme <- function(color_scheme) {
  color_scheme <- as.character(color_scheme)[1L]
  if (identical(color_scheme, "symmetric")) {
    return("yellow_blue")
  }
  valid_schemes <- c("yellow_blue", "rainbow")
  if (!color_scheme %in% valid_schemes) {
    stop(
      sprintf(
        "`color_scheme` must be one of %s.",
        paste(sprintf('"%s"', c(valid_schemes, "symmetric")), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  color_scheme
}

assign_s1_theta_colors <- function(theta_values, mu, color_scheme) {
  color_scheme <- normalize_s1_color_scheme(color_scheme)
  if (identical(color_scheme, "rainbow")) {
    return(make_theta_rainbow_colors(theta_values))
  }
  make_theta_colors(
    theta_values,
    reference_theta = circle_angle_from_point(mu)
  )
}

format_pi_fraction_label <- function(theta, max_denominator = 360, tol = 1e-10) {
  theta <- wrap_angle_2pi(as.numeric(theta)[1L])
  if (abs(theta) < tol) {
    return("0")
  }

  ratio <- theta / pi
  denominators <- seq_len(max_denominator)
  numerators <- round(ratio * denominators)
  errors <- abs(ratio - numerators / denominators)
  best_idx <- which.min(errors)
  num <- as.integer(numerators[[best_idx]])
  den <- as.integer(denominators[[best_idx]])

  if (!is.finite(num) || !is.finite(den) || den <= 0 || errors[[best_idx]] > tol) {
    return(format(theta, digits = 3, trim = TRUE, scientific = FALSE))
  }

  if (den == 1L) {
    if (num == 1L) return("\u03c0")
    return(paste0(num, "\u03c0"))
  }
  if (num == 1L) {
    return(paste0("\u03c0/", den))
  }
  paste0(num, "\u03c0/", den)
}

format_radians_label <- function(theta_values, digits = 3) {
  vapply(
    theta_values,
    function(theta) {
      label <- format_pi_fraction_label(theta)
      if (nzchar(label)) {
        label
      } else {
        format(theta, digits = digits, trim = TRUE, scientific = FALSE)
      }
    },
    character(1)
  )
}

draw_single_gaussian_realization <- function(cov_matrix, seed = NULL, tol = 1e-10) {
  validate_covariance_matrix(
    cov_matrix,
    symmetry_tol = tol,
    psd_tol = tol,
    stop_on_failure = TRUE
  )
  if (!is.null(seed)) {
    set.seed(seed)
  }
  z_standard <- stats::rnorm(nrow(cov_matrix))
  chol_upper <- tryCatch(
    chol(cov_matrix),
    error = function(e) {
      stop(
        paste0(
          "Cholesky factorization failed for the covariance matrix. ",
          "This visualization now simulates from N(0, I) and applies the global Cholesky factor; ",
          "the supplied covariance must therefore be numerically positive definite. Original error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  realization_vec <- as.numeric(t(chol_upper) %*% z_standard)
  list(
    realization_vec = realization_vec,
    standard_normal_draw = z_standard,
    chol_upper = chol_upper
  )
}

default_s1_t_grid <- function(t_points) {
  seq(0.05, 1.95, length.out = t_points)
}

save_combined_s1_plot <- function(left_panel,
                                  right_panel,
                                  save_path,
                                  width = 12,
                                  height = 6,
                                  dpi = 300,
                                  bg = "transparent") {
  grDevices::png(
    filename = save_path,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = bg
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 1,
        ncol = 2,
        widths = grid::unit(c(1, 1), "null")
      )
    )
  )
  print(left_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(right_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}

save_s1_plot_bundle <- function(plot_result, save_result) {
  saveRDS(plot_result, file = save_result)
  invisible(save_result)
}

build_s1_panels_from_curve_data <- function(curve_data,
                                            circle_grid,
                                            ray_extension,
                                            right_x_label,
                                            right_y_label) {
  ray_data <- data.frame(
    theta = circle_grid$theta,
    x = circle_grid$x,
    y = circle_grid$y,
    xend = (1 + ray_extension) * circle_grid$x,
    yend = (1 + ray_extension) * circle_grid$y,
    color = circle_grid$color,
    stringsAsFactors = FALSE
  )
  circle_outline <- data.frame(
    x = cos(seq(0, 2 * pi, length.out = 600)),
    y = sin(seq(0, 2 * pi, length.out = 600))
  )

  left_panel <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = circle_outline,
      ggplot2::aes(x = x, y = y),
      color = "grey30",
      linewidth = 0.9
    ) +
    ggplot2::geom_segment(
      data = ray_data,
      ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend, color = color),
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = curve_data,
      ggplot2::aes(x = x_circle, y = y_circle, group = omega_id, color = color),
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = circle_grid,
      ggplot2::aes(x = x, y = y, color = color),
      size = 2,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void()

  legend_breaks <- unique(curve_data$theta_label)
  color_values <- stats::setNames(circle_grid$color, circle_grid$legend_label)

  right_panel <- ggplot2::ggplot(
    curve_data,
    ggplot2::aes(x = plot_x, y = g, group = omega_id, color = theta_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = legend_breaks,
      name = "Radians"
    ) +
    ggplot2::labs(
      x = right_x_label,
      y = right_y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  list(left_panel = left_panel, right_panel = right_panel)
}

write_s1_plot_files <- function(left_panel, right_panel, save_plot) {
  save_ext <- tools::file_ext(save_plot)
  save_base <- if (nzchar(save_ext)) {
    sub(paste0("\\.", save_ext, "$"), "", save_plot)
  } else {
    save_plot
  }
  combined_path <- paste0(save_base, ifelse(nzchar(save_ext), paste0(".", save_ext), ".png"))
  left_path <- paste0(save_base, "_left.", ifelse(nzchar(save_ext), save_ext, "png"))
  right_path <- paste0(save_base, "_right.", ifelse(nzchar(save_ext), save_ext, "png"))
  ggplot2::ggsave(left_path, left_panel, width = 6, height = 6, dpi = 300)
  ggplot2::ggsave(right_path, right_panel, width = 7, height = 5, dpi = 300, bg = "transparent")
  save_combined_s1_plot(left_panel, right_panel, combined_path)
  invisible(c(combined = combined_path, left = left_path, right = right_path))
}


# ============================================================================
# PUBLIC VISUALIZATION
# ============================================================================

#' Visualize one realization of the vMF limit Gaussian process on S^1
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param n_angles Number of angular directions on the circle
#' @param t_points Number of interior t-points used to discretize the process
#' @param omega_grid Optional matrix of points on S^1
#' @param t_grid Optional vector of chordal thresholds in (0, 2)
#' @param cov_method Either "exact_s1_simple" or "mc"
#' @param n_mc_samples Monte Carlo sample size for the mc covariance route
#' @param n_cores Number of cores for covariance construction
#' @param seed Optional seed
#' @param curve_points Number of probability-grid points for plotting
#' @param plot_domain Either "u" or "t" for the display parametrization
#' @param color_scheme Either "yellow_blue" or "rainbow". The legacy alias
#'   "symmetric" is accepted and mapped to "yellow_blue"; this affects only the
#'   palette, not the stochastic process.
#' @param quiet Whether to suppress console output during covariance construction
#' @param ray_extension Radial extension beyond the unit circle
#' @param lateral_scale Relative amplitude used for left-panel mini-curves
#' @param cdf_grid_size Deterministic angular grid size used for S^1 CDF computations
#' @param save_plot Optional file path stem for saving the separate panels
#' @param save_result Optional path to an .rds bundle storing the simulated
#'   realization and all plot data so the figure can be regenerated later
#' @return List containing the plot, data, and covariance diagnostics
visualize_limit_gaussian_s1_vmf <- function(mu = c(1, 0),
                                            kappa = 2,
                                            n_angles = 6,
                                            t_points = 160,
                                            omega_grid = NULL,
                                            t_grid = NULL,
                                            cov_method = c("exact_s1_simple", "mc"),
                                            n_mc_samples = 5000,
                                            n_cores = 1,
                                            seed = NULL,
                                            curve_points = 200,
                                            plot_domain = c("u", "t"),
                                            color_scheme = c("yellow_blue", "rainbow", "symmetric"),
                                            quiet = TRUE,
                                            ray_extension = 0.35,
                                            lateral_scale = 0.10,
                                            cdf_grid_size = 16385,
                                            save_plot = NULL,
                                            save_result = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to generate the S1 visualization.")
  }

  cov_method <- match.arg(cov_method)
  plot_domain <- match.arg(plot_domain)
  color_scheme <- normalize_s1_color_scheme(match.arg(color_scheme))
  mu <- as.numeric(mu)
  if (length(mu) != 2) {
    stop("`mu` must have length 2 for S^1 visualization.")
  }
  mu <- mu / sqrt(sum(mu^2))
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(omega_grid)) {
    circle_grid <- generate_circle_grid(n_angles, theta0 = circle_angle_from_point(mu))
    omega_grid <- as.matrix(circle_grid[, c("x", "y")])
  } else {
    omega_grid <- as.matrix(omega_grid)
    if (ncol(omega_grid) != 2) {
      stop("`omega_grid` must have exactly 2 columns for S^1 visualization.")
    }
    theta_values <- circle_angles_from_matrix(omega_grid)
    ordering <- order(theta_values)
    omega_grid <- omega_grid[ordering, , drop = FALSE]
    theta_values <- theta_values[ordering]
    circle_grid <- data.frame(
      theta = theta_values,
      x = omega_grid[, 1],
      y = omega_grid[, 2],
      label = format_radians_label(theta_values),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(t_grid)) {
    t_grid <- default_s1_t_grid(t_points)
  } else {
    t_grid <- as.numeric(t_grid)
    t_points <- length(t_grid)
  }
  if (any(t_grid <= 0) || any(t_grid >= 2)) {
    stop("`t_grid` must lie strictly inside (0, 2) for the visualization grid.")
  }

  cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  theta_values <- circle_grid$theta
  colors <- assign_s1_theta_colors(theta_values, mu, color_scheme)
  circle_grid$color <- colors
  circle_grid$legend_label <- format_radians_label(theta_values)

  cov_matrix <- if (isTRUE(quiet)) {
    capture.output({
      cov_matrix_local <- cov_vmf(
        omega_grid = omega_grid,
        t_grid = t_grid,
        mu = mu,
        kappa = kappa,
        distance_type = "chordal",
        n_mc_samples = n_mc_samples,
        n_cores = n_cores,
        seed = seed,
        h0 = "simple",
        unknown_param = NULL,
        cov_method = cov_method,
        cdf_grid_size = cdf_grid_size
      )
    })
    cov_matrix_local
  } else {
    cov_vmf(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = mu,
      kappa = kappa,
      distance_type = "chordal",
      n_mc_samples = n_mc_samples,
      n_cores = n_cores,
      seed = seed,
      h0 = "simple",
      unknown_param = NULL,
      cov_method = cov_method,
      cdf_grid_size = cdf_grid_size
    )
  }
  cov_diagnostics <- validate_covariance_matrix(
    cov_matrix,
    symmetry_tol = 1e-10,
    psd_tol = 1e-10,
    stop_on_failure = TRUE
  )

  gaussian_draw <- draw_single_gaussian_realization(cov_matrix, seed = seed, tol = 1e-10)
  realization_vec <- gaussian_draw$realization_vec
  g_inner <- matrix(realization_vec, nrow = nrow(omega_grid), ncol = length(t_grid))

  curve_data_list <- vector("list", nrow(omega_grid))
  if (plot_domain == "u") {
    u_grid <- seq(0, 1, length.out = curve_points)
    for (i in seq_len(nrow(omega_grid))) {
      omega <- omega_grid[i, ]
      t_from_u <- invert_distance_profile_vmf_s1_chordal(
        omega = omega,
        mu = mu,
        kappa = kappa,
        u_values = u_grid,
        cdf_object = cdf_object,
        cdf_grid_size = cdf_grid_size
      )
      g_aug <- c(0, g_inner[i, ], 0)
      t_aug <- c(0, t_grid, 2)
      g_u <- approx(
        x = t_aug,
        y = g_aug,
        xout = t_from_u,
        method = "linear",
        ties = "ordered",
        rule = 2
      )$y
      curve_data_list[[i]] <- data.frame(
        omega_id = i,
        theta = theta_values[i],
        theta_label = circle_grid$legend_label[i],
        color = colors[i],
        u = u_grid,
        t = t_from_u,
        g = g_u,
        plot_x = u_grid,
        radial_progress = u_grid,
        stringsAsFactors = FALSE
      )
    }
    right_x_label <- "u"
    right_y_label <- expression(G(omega, F^{-1}*"("*u*")"))
  } else {
    t_plot <- seq(0, 2, length.out = curve_points)
    for (i in seq_len(nrow(omega_grid))) {
      omega <- omega_grid[i, ]
      g_aug <- c(0, g_inner[i, ], 0)
      t_aug <- c(0, t_grid, 2)
      g_t <- approx(
        x = t_aug,
        y = g_aug,
        xout = t_plot,
        method = "linear",
        ties = "ordered",
        rule = 2
      )$y
      u_from_t <- theoretical_distance_profile_vmf_s1_chordal(
        omega = omega,
        mu = mu,
        kappa = kappa,
        t_values = t_plot,
        cdf_object = cdf_object,
        cdf_grid_size = cdf_grid_size
      )
      curve_data_list[[i]] <- data.frame(
        omega_id = i,
        theta = theta_values[i],
        theta_label = circle_grid$legend_label[i],
        color = colors[i],
        u = u_from_t,
        t = t_plot,
        g = g_t,
        plot_x = t_plot,
        radial_progress = t_plot / 2,
        stringsAsFactors = FALSE
      )
    }
    right_x_label <- "t"
    right_y_label <- expression(G(omega, t))
  }
  curve_data <- do.call(rbind, curve_data_list)

  max_abs_g <- max(abs(curve_data$g))
  amplitude_scale <- if (max_abs_g > 0) lateral_scale / max_abs_g else 0
  curve_data$x_circle <- (1 + ray_extension * curve_data$radial_progress) * cos(curve_data$theta) -
    amplitude_scale * curve_data$g * sin(curve_data$theta)
  curve_data$y_circle <- (1 + ray_extension * curve_data$radial_progress) * sin(curve_data$theta) +
    amplitude_scale * curve_data$g * cos(curve_data$theta)

  panels <- build_s1_panels_from_curve_data(
    curve_data = curve_data,
    circle_grid = circle_grid,
    ray_extension = ray_extension,
    right_x_label = right_x_label,
    right_y_label = right_y_label
  )
  left_panel <- panels$left_panel
  right_panel <- panels$right_panel

  result <- list(
    plot = list(left_panel = left_panel, right_panel = right_panel),
    left_panel = left_panel,
    right_panel = right_panel,
    curve_data = curve_data,
    left_panel_data = curve_data[, c("omega_id", "theta", "u", "g", "x_circle", "y_circle", "color")],
    cov_matrix = cov_matrix,
    cov_diagnostics = cov_diagnostics,
    realization_vec = realization_vec,
    standard_normal_draw = gaussian_draw$standard_normal_draw,
    chol_upper = gaussian_draw$chol_upper,
    omega_grid = omega_grid,
    t_grid = t_grid,
    plot_grid = curve_data$plot_x,
    plot_domain = plot_domain,
    color_scheme = color_scheme,
    cov_method = cov_method,
    mu = mu,
    kappa = kappa,
    circle_grid = circle_grid,
    ray_extension = ray_extension,
    lateral_scale = lateral_scale,
    right_x_label = right_x_label,
    right_y_label = right_y_label,
    cdf_grid_size = cdf_grid_size
  )

  if (!is.null(save_plot)) {
    write_s1_plot_files(left_panel, right_panel, save_plot)
  }
  if (!is.null(save_result)) {
    save_s1_plot_bundle(result, save_result)
  }

  result
}

#' Rebuild an S1 visualization from a stored result bundle
#' @param result A result list returned by visualize_limit_gaussian_s1_vmf, or a
#'   path to a saved .rds bundle produced via `save_result`
#' @param color_scheme Either "yellow_blue" or "rainbow". The legacy alias
#'   "symmetric" is accepted and mapped to "yellow_blue"
#' @param save_plot Optional file path stem for saving the regenerated panels
#' @return List containing the regenerated plot objects and updated plot data
replot_limit_gaussian_s1_vmf_from_result <- function(result,
                                                     color_scheme = c("yellow_blue", "rainbow", "symmetric"),
                                                     save_plot = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to generate the S1 visualization.")
  }

  if (is.character(result) && length(result) == 1L) {
    result <- readRDS(result)
  }
  if (!is.list(result) || is.null(result$curve_data) || is.null(result$circle_grid)) {
    stop("`result` must be a visualization result list or a valid saved .rds bundle.")
  }

  color_scheme <- normalize_s1_color_scheme(match.arg(color_scheme))
  theta_values <- result$circle_grid$theta
  colors <- assign_s1_theta_colors(theta_values, result$mu, color_scheme)

  circle_grid <- result$circle_grid
  circle_grid$color <- colors
  circle_grid$legend_label <- format_radians_label(theta_values)

  curve_data <- result$curve_data
  color_map <- stats::setNames(colors, theta_values)
  label_map <- stats::setNames(circle_grid$legend_label, theta_values)
  curve_data$color <- unname(color_map[as.character(curve_data$theta)])
  curve_data$theta_label <- unname(label_map[as.character(curve_data$theta)])

  panels <- build_s1_panels_from_curve_data(
    curve_data = curve_data,
    circle_grid = circle_grid,
    ray_extension = result$ray_extension,
    right_x_label = result$right_x_label,
    right_y_label = result$right_y_label
  )

  replotted_result <- result
  replotted_result$plot <- panels
  replotted_result$left_panel <- panels$left_panel
  replotted_result$right_panel <- panels$right_panel
  replotted_result$curve_data <- curve_data
  replotted_result$left_panel_data <- curve_data[, c("omega_id", "theta", "u", "g", "x_circle", "y_circle", "color")]
  replotted_result$circle_grid <- circle_grid
  replotted_result$color_scheme <- color_scheme

  if (!is.null(save_plot)) {
    write_s1_plot_files(panels$left_panel, panels$right_panel, save_plot)
  }

  replotted_result
}


# ============================================================================
# BENCHMARKING AND DIAGNOSTICS
# These functions generate validation artifacts. When saved to disk, they
# should go under tests/benchmark_outputs/, not under routine output/ folders.
# ============================================================================

#' Benchmark Monte Carlo and exact S1 covariance routes for vMF
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param scenarios List of benchmark scenarios
#' @param mc_sizes Monte Carlo sample sizes for the mc route
#' @param n_cores Number of cores for covariance construction
#' @param seed Optional seed
#' @param cdf_grid_size Deterministic angular grid size used for S^1 CDF computations
#' @param output_dir Optional directory for saving benchmark tables and plots
#' @param save_plots Whether to save benchmark plots when output_dir is provided
#' @return List with the metrics data frame and optional ggplot objects
benchmark_covariance_vmf_s1 <- function(mu = c(1, 0),
                                        kappa = 2,
                                        scenarios = list(
                                          list(name = "small", n_angles = 6, t_points = 25),
                                          list(name = "medium", n_angles = 8, t_points = 50)
                                        ),
                                        mc_sizes = c(1000, 5000, 10000),
                                        n_cores = 1,
                                        seed = 123,
                                        cdf_grid_size = 16385,
                                        output_dir = NULL,
                                        save_plots = FALSE) {
  metrics_rows <- list()
  row_idx <- 0

  for (scenario in scenarios) {
    circle_grid <- generate_circle_grid(scenario$n_angles)
    omega_grid <- as.matrix(circle_grid[, c("x", "y")])
    t_grid <- default_s1_t_grid(scenario$t_points)

    exact_start <- Sys.time()
    sigma_exact <- cov_vmf(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = mu,
      kappa = kappa,
      distance_type = "chordal",
      n_cores = n_cores,
      seed = seed,
      h0 = "simple",
      cov_method = "exact_s1_simple",
      cdf_grid_size = cdf_grid_size
    )
    exact_time <- as.numeric(difftime(Sys.time(), exact_start, units = "secs"))
    exact_diag <- validate_covariance_matrix(
      sigma_exact,
      symmetry_tol = 1e-10,
      psd_tol = 1e-10,
      stop_on_failure = FALSE
    )

    row_idx <- row_idx + 1
    metrics_rows[[row_idx]] <- data.frame(
      scenario = scenario$name,
      n_angles = scenario$n_angles,
      t_points = scenario$t_points,
      method = "exact_s1_simple",
      n_mc_samples = NA_integer_,
      time_sec = exact_time,
      max_abs_diff = 0,
      rel_fro_error = 0,
      symmetry_gap = exact_diag$symmetry_gap,
      min_eigenvalue = exact_diag$min_eigenvalue,
      stringsAsFactors = FALSE
    )

    for (mc_n in mc_sizes) {
      mc_start <- Sys.time()
      sigma_mc <- cov_vmf(
        omega_grid = omega_grid,
        t_grid = t_grid,
        mu = mu,
        kappa = kappa,
        distance_type = "chordal",
        n_mc_samples = mc_n,
        n_cores = n_cores,
        seed = seed,
        h0 = "simple",
        cov_method = "mc",
        cdf_grid_size = cdf_grid_size
      )
      mc_time <- as.numeric(difftime(Sys.time(), mc_start, units = "secs"))
      mc_diag <- validate_covariance_matrix(
        sigma_mc,
        symmetry_tol = 1e-10,
        psd_tol = 1e-10,
        stop_on_failure = FALSE
      )

      diff_mat <- sigma_mc - sigma_exact
      fro_exact <- sqrt(sum(sigma_exact^2))
      row_idx <- row_idx + 1
      metrics_rows[[row_idx]] <- data.frame(
        scenario = scenario$name,
        n_angles = scenario$n_angles,
        t_points = scenario$t_points,
        method = "mc",
        n_mc_samples = mc_n,
        time_sec = mc_time,
        max_abs_diff = max(abs(diff_mat)),
        rel_fro_error = sqrt(sum(diff_mat^2)) / ifelse(fro_exact > 0, fro_exact, 1),
        symmetry_gap = mc_diag$symmetry_gap,
        min_eigenvalue = mc_diag$min_eigenvalue,
        stringsAsFactors = FALSE
      )
    }
  }

  metrics <- do.call(rbind, metrics_rows)
  error_plot <- NULL
  time_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    mc_metrics <- subset(metrics, method == "mc")
    error_plot <- ggplot2::ggplot(
      mc_metrics,
      ggplot2::aes(x = n_mc_samples, y = rel_fro_error, color = scenario)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Monte Carlo sample size",
        y = "Relative Frobenius error",
        title = "MC vs Exact S1 Covariance Error"
      ) +
      ggplot2::theme_minimal()

    time_plot <- ggplot2::ggplot(
      metrics,
      ggplot2::aes(
        x = ifelse(is.na(n_mc_samples), 0, n_mc_samples),
        y = time_sec,
        color = method,
        shape = scenario
      )
    ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_line(data = subset(metrics, method == "mc")) +
      ggplot2::labs(
        x = "Monte Carlo sample size (0 denotes exact route)",
        y = "Time (seconds)",
        title = "Covariance Construction Time"
      ) +
      ggplot2::theme_minimal()
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(metrics, file.path(output_dir, "s1_covariance_benchmark.csv"), row.names = FALSE)
    if (isTRUE(save_plots) && !is.null(error_plot) && !is.null(time_plot)) {
      ggplot2::ggsave(file.path(output_dir, "s1_covariance_error.png"), error_plot, width = 7, height = 5, dpi = 300)
      ggplot2::ggsave(file.path(output_dir, "s1_covariance_time.png"), time_plot, width = 7, height = 5, dpi = 300)
    }
  }

  list(
    metrics = metrics,
    error_plot = error_plot,
    time_plot = time_plot
  )
}

#' Reference probability on S^1 by adaptive quadrature over the exact arc intersection
#' @param omega1 First point on S^1
#' @param t1 First chordal threshold
#' @param omega2 Second point on S^1
#' @param t2 Second chordal threshold
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param rel.tol Relative tolerance passed to integrate()
#' @param abs.tol Absolute tolerance passed to integrate()
#' @return Joint probability computed by adaptive quadrature
joint_probability_vmf_s1_chordal_reference <- function(omega1,
                                                       t1,
                                                       omega2,
                                                       t2,
                                                       mu,
                                                       kappa,
                                                       rel.tol = 1e-11,
                                                       abs.tol = 1e-13) {
  seg1 <- s1_event_segments_chordal(omega1, t1)
  seg2 <- s1_event_segments_chordal(omega2, t2)
  overlap <- s1_intersect_segments(seg1, seg2)
  if (length(overlap) == 0 || nrow(overlap) == 0) {
    return(0)
  }

  total_prob <- 0
  for (i in seq_len(nrow(overlap))) {
    total_prob <- total_prob + integrate(
      f = function(phi) vmf_s1_angle_density(phi, mu = mu, kappa = kappa),
      lower = overlap[i, 1],
      upper = overlap[i, 2],
      rel.tol = rel.tol,
      abs.tol = abs.tol
    )$value
  }
  total_prob
}

#' Reference marginal probability on S^1 by adaptive quadrature
#' @param omega Point on S^1
#' @param t Chordal threshold
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param rel.tol Relative tolerance passed to integrate()
#' @param abs.tol Absolute tolerance passed to integrate()
#' @return Marginal probability computed by adaptive quadrature
distance_profile_vmf_s1_chordal_reference <- function(omega,
                                                      t,
                                                      mu,
                                                      kappa,
                                                      rel.tol = 1e-11,
                                                      abs.tol = 1e-13) {
  segments <- s1_event_segments_chordal(omega, t)
  if (length(segments) == 0 || nrow(segments) == 0) {
    return(0)
  }
  total_prob <- 0
  for (i in seq_len(nrow(segments))) {
    total_prob <- total_prob + integrate(
      f = function(phi) vmf_s1_angle_density(phi, mu = mu, kappa = kappa),
      lower = segments[i, 1],
      upper = segments[i, 2],
      rel.tol = rel.tol,
      abs.tol = abs.tol
    )$value
  }
  total_prob
}

#' Build a covariance matrix reference on S^1 by adaptive quadrature
#' @param omega_grid Matrix of points on S^1
#' @param t_grid Vector of chordal thresholds
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param rel.tol Relative tolerance passed to integrate()
#' @param abs.tol Absolute tolerance passed to integrate()
#' @return Covariance matrix under the simple null
cov_vmf_s1_simple_reference <- function(omega_grid,
                                        t_grid,
                                        mu,
                                        kappa,
                                        rel.tol = 1e-11,
                                        abs.tol = 1e-13) {
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t

  F_matrix <- matrix(0, nrow = n_omega, ncol = n_t)
  for (i in seq_len(n_omega)) {
    for (k in seq_len(n_t)) {
      F_matrix[i, k] <- distance_profile_vmf_s1_chordal_reference(
        omega = omega_grid[i, ],
        t = t_grid[k],
        mu = mu,
        kappa = kappa,
        rel.tol = rel.tol,
        abs.tol = abs.tol
      )
    }
  }
  F_vec <- as.vector(F_matrix)

  cov_matrix <- matrix(0, nrow = n_total, ncol = n_total)
  for (i in seq_len(n_total)) {
    omega_i <- ((i - 1) %% n_omega) + 1
    t_i <- ((i - 1) %/% n_omega) + 1
    for (j in i:n_total) {
      omega_j <- ((j - 1) %% n_omega) + 1
      t_j <- ((j - 1) %/% n_omega) + 1
      p_joint <- joint_probability_vmf_s1_chordal_reference(
        omega1 = omega_grid[omega_i, ],
        t1 = t_grid[t_i],
        omega2 = omega_grid[omega_j, ],
        t2 = t_grid[t_j],
        mu = mu,
        kappa = kappa,
        rel.tol = rel.tol,
        abs.tol = abs.tol
      )
      cov_value <- p_joint - F_vec[i] * F_vec[j]
      cov_matrix[i, j] <- cov_value
      cov_matrix[j, i] <- cov_value
    }
  }

  cov_matrix
}

#' Benchmark joint probabilities on S^1 against an independent adaptive-quadrature reference
#' @param mu Mean direction on S^1
#' @param kappas Vector of concentration parameters to benchmark
#' @param cases Optional list of benchmark cases
#' @param mc_sizes Monte Carlo sample sizes for the MC estimator
#' @param seed Optional seed
#' @param cdf_grid_sizes Deterministic CDF grid sizes for the current "exact" route
#' @param output_dir Optional directory for saving tables and plots
#' @param save_plots Whether to save plots when output_dir is provided
#' @return List with the benchmark table and optional plots
benchmark_joint_probability_vmf_s1 <- function(mu = c(1, 0),
                                               kappas = c(2, 8),
                                               cases = NULL,
                                               mc_sizes = c(1000, 5000, 10000, 50000),
                                               seed = 123,
                                               cdf_grid_sizes = c(2049, 4097, 8193, 16385),
                                               output_dir = NULL,
                                               save_plots = FALSE) {
  if (is.null(cases)) {
    cases <- list(
      list(name = "same_center", theta1 = 0, theta2 = 0, t1 = 0.70, t2 = 1.05),
      list(name = "partial_overlap", theta1 = 0, theta2 = 0.85, t1 = 0.95, t2 = 1.05),
      list(name = "near_tangent", theta1 = 0, theta2 = 1.15, t1 = 0.60, t2 = 0.60),
      list(name = "wraparound", theta1 = 6.00, theta2 = 0.35, t1 = 0.90, t2 = 0.85)
    )
  }

  if (!is.null(seed)) set.seed(seed)
  results <- list()
  row_idx <- 0

  for (kappa in kappas) {
    for (case in cases) {
      omega1 <- c(cos(case$theta1), sin(case$theta1))
      omega2 <- c(cos(case$theta2), sin(case$theta2))

      ref_start <- Sys.time()
      reference_prob <- joint_probability_vmf_s1_chordal_reference(
        omega1 = omega1,
        t1 = case$t1,
        omega2 = omega2,
        t2 = case$t2,
        mu = mu,
        kappa = kappa
      )
      ref_time <- as.numeric(difftime(Sys.time(), ref_start, units = "secs"))

      row_idx <- row_idx + 1
      results[[row_idx]] <- data.frame(
        kappa = kappa,
        case = case$name,
        method = "reference_integrate",
        tuning = NA_character_,
        estimate = reference_prob,
        abs_error = 0,
        rel_error = 0,
        time_sec = ref_time,
        stringsAsFactors = FALSE
      )

      for (grid_size in cdf_grid_sizes) {
        det_start <- Sys.time()
        det_prob <- joint_probability_vmf_s1_chordal_exact(
          omega1 = omega1,
          t1 = case$t1,
          omega2 = omega2,
          t2 = case$t2,
          mu = mu,
          kappa = kappa,
          cdf_grid_size = grid_size
        )
        det_time <- as.numeric(difftime(Sys.time(), det_start, units = "secs"))
        row_idx <- row_idx + 1
        results[[row_idx]] <- data.frame(
          kappa = kappa,
          case = case$name,
          method = "deterministic_cdf",
          tuning = paste0("n_grid=", grid_size),
          estimate = det_prob,
          abs_error = abs(det_prob - reference_prob),
          rel_error = abs(det_prob - reference_prob) / ifelse(reference_prob > 0, reference_prob, 1),
          time_sec = det_time,
          stringsAsFactors = FALSE
        )
      }

      for (mc_n in mc_sizes) {
        if (!is.null(seed)) set.seed(seed)
        mc_start <- Sys.time()
        sample_points <- rotasym::r_vMF(n = mc_n, mu = mu, kappa = kappa)
        mc_prob <- mean(
          sqrt(rowSums((sample_points - matrix(omega1, nrow = mc_n, ncol = 2, byrow = TRUE))^2)) <= case$t1 &
            sqrt(rowSums((sample_points - matrix(omega2, nrow = mc_n, ncol = 2, byrow = TRUE))^2)) <= case$t2
        )
        mc_time <- as.numeric(difftime(Sys.time(), mc_start, units = "secs"))
        row_idx <- row_idx + 1
        results[[row_idx]] <- data.frame(
          kappa = kappa,
          case = case$name,
          method = "mc",
          tuning = paste0("n_mc=", mc_n),
          estimate = mc_prob,
          abs_error = abs(mc_prob - reference_prob),
          rel_error = abs(mc_prob - reference_prob) / ifelse(reference_prob > 0, reference_prob, 1),
          time_sec = mc_time,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  benchmark_df <- do.call(rbind, results)
  error_plot <- NULL
  time_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    comparable <- subset(benchmark_df, method != "reference_integrate")
    error_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = abs_error, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(case ~ kappa, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Absolute error vs adaptive quadrature",
        title = "Joint-Probability Benchmark On S1"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    time_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = time_sec, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(case ~ kappa, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Time (seconds)",
        title = "Joint-Probability Timing On S1"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      benchmark_df,
      file.path(output_dir, "s1_joint_probability_benchmark.csv"),
      row.names = FALSE
    )
    if (isTRUE(save_plots) && !is.null(error_plot) && !is.null(time_plot)) {
      ggplot2::ggsave(
        file.path(output_dir, "s1_joint_probability_error.png"),
        error_plot,
        width = 10,
        height = 7,
        dpi = 300
      )
      ggplot2::ggsave(
        file.path(output_dir, "s1_joint_probability_time.png"),
        time_plot,
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }

  list(
    metrics = benchmark_df,
    error_plot = error_plot,
    time_plot = time_plot
  )
}

#' Benchmark covariance matrices on S^1 against an external adaptive-quadrature reference
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param n_angles Number of fixed omega directions
#' @param t_points Number of t values per omega
#' @param mc_sizes Monte Carlo sample sizes for the MC route
#' @param seed Optional seed
#' @param cdf_grid_sizes Deterministic CDF grid sizes for the deterministic route
#' @param reference_rel_tol Relative tolerance for the reference integrals
#' @param reference_abs_tol Absolute tolerance for the reference integrals
#' @param output_dir Optional directory for saving tables and plots
#' @param save_plots Whether to save plots when output_dir is provided
#' @return List with the benchmark table, plots, and reference matrix diagnostics
benchmark_covariance_vmf_s1_reference <- function(mu = c(1, 0),
                                                  kappa = 2,
                                                  n_angles = 6,
                                                  t_points = 10,
                                                  mc_sizes = c(10000, 50000, 100000),
                                                  seed = 123,
                                                  cdf_grid_sizes = c(2049, 4097, 8193, 16385),
                                                  reference_rel_tol = 1e-11,
                                                  reference_abs_tol = 1e-13,
                                                  output_dir = NULL,
                                                  save_plots = FALSE) {
  circle_grid <- generate_circle_grid(n_angles)
  omega_grid <- as.matrix(circle_grid[, c("x", "y")])
  t_grid <- default_s1_t_grid(t_points)

  ref_start <- Sys.time()
  sigma_ref <- cov_vmf_s1_simple_reference(
    omega_grid = omega_grid,
    t_grid = t_grid,
    mu = mu,
    kappa = kappa,
    rel.tol = reference_rel_tol,
    abs.tol = reference_abs_tol
  )
  ref_time <- as.numeric(difftime(Sys.time(), ref_start, units = "secs"))
  ref_diag <- validate_covariance_matrix(
    sigma_ref,
    symmetry_tol = 1e-10,
    psd_tol = 1e-10,
    stop_on_failure = FALSE
  )

  rows <- list()
  row_idx <- 1
  rows[[row_idx]] <- data.frame(
    method = "reference_integrate",
    tuning = NA_character_,
    time_sec = ref_time,
    max_abs_diff = 0,
    rel_fro_error = 0,
    symmetry_gap = ref_diag$symmetry_gap,
    min_eigenvalue = ref_diag$min_eigenvalue,
    stringsAsFactors = FALSE
  )

  fro_ref <- sqrt(sum(sigma_ref^2))

  for (grid_size in cdf_grid_sizes) {
    det_start <- Sys.time()
    sigma_det <- cov_vmf(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = mu,
      kappa = kappa,
      distance_type = "chordal",
      n_cores = 1,
      seed = seed,
      h0 = "simple",
      cov_method = "exact_s1_simple",
      cdf_grid_size = grid_size
    )
    det_time <- as.numeric(difftime(Sys.time(), det_start, units = "secs"))
    det_diag <- validate_covariance_matrix(
      sigma_det,
      symmetry_tol = 1e-10,
      psd_tol = 1e-10,
      stop_on_failure = FALSE
    )
    diff_mat <- sigma_det - sigma_ref
    row_idx <- row_idx + 1
    rows[[row_idx]] <- data.frame(
      method = "deterministic_cdf",
      tuning = paste0("n_grid=", grid_size),
      time_sec = det_time,
      max_abs_diff = max(abs(diff_mat)),
      rel_fro_error = sqrt(sum(diff_mat^2)) / ifelse(fro_ref > 0, fro_ref, 1),
      symmetry_gap = det_diag$symmetry_gap,
      min_eigenvalue = det_diag$min_eigenvalue,
      stringsAsFactors = FALSE
    )
  }

  for (mc_n in mc_sizes) {
    mc_start <- Sys.time()
    sigma_mc <- cov_vmf(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = mu,
      kappa = kappa,
      distance_type = "chordal",
      n_mc_samples = mc_n,
      n_cores = 1,
      seed = seed,
      h0 = "simple",
      cov_method = "mc"
    )
    mc_time <- as.numeric(difftime(Sys.time(), mc_start, units = "secs"))
    mc_diag <- validate_covariance_matrix(
      sigma_mc,
      symmetry_tol = 1e-10,
      psd_tol = 1e-10,
      stop_on_failure = FALSE
    )
    diff_mat <- sigma_mc - sigma_ref
    row_idx <- row_idx + 1
    rows[[row_idx]] <- data.frame(
      method = "mc",
      tuning = paste0("n_mc=", mc_n),
      time_sec = mc_time,
      max_abs_diff = max(abs(diff_mat)),
      rel_fro_error = sqrt(sum(diff_mat^2)) / ifelse(fro_ref > 0, fro_ref, 1),
      symmetry_gap = mc_diag$symmetry_gap,
      min_eigenvalue = mc_diag$min_eigenvalue,
      stringsAsFactors = FALSE
    )
  }

  metrics <- do.call(rbind, rows)
  error_plot <- NULL
  time_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    comparable <- subset(metrics, method != "reference_integrate")
    error_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = rel_fro_error, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Relative Frobenius error vs reference_integrate",
        title = sprintf("Covariance Benchmark On S1 (%d angles x %d t's)", n_angles, t_points)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    time_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = time_sec, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Time (seconds)",
        title = sprintf("Covariance Timing On S1 (%d angles x %d t's)", n_angles, t_points)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      metrics,
      file.path(output_dir, sprintf("s1_covariance_reference_benchmark_angles%d_t%d.csv", n_angles, t_points)),
      row.names = FALSE
    )
    if (isTRUE(save_plots) && !is.null(error_plot) && !is.null(time_plot)) {
      ggplot2::ggsave(
        file.path(output_dir, sprintf("s1_covariance_reference_error_angles%d_t%d.png", n_angles, t_points)),
        error_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
      ggplot2::ggsave(
        file.path(output_dir, sprintf("s1_covariance_reference_time_angles%d_t%d.png", n_angles, t_points)),
        time_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
    }
  }

  list(
    metrics = metrics,
    error_plot = error_plot,
    time_plot = time_plot,
    sigma_reference = sigma_ref,
    reference_diagnostics = ref_diag,
    omega_grid = omega_grid,
    t_grid = t_grid
  )
}

cat("S1 visualization and benchmark functions loaded successfully!\n")
