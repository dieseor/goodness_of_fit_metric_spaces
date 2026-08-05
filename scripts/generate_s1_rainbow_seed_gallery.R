#!/usr/bin/env Rscript

# Generate several S^1 Gaussian-process realizations from one cached exact
# covariance factorization.  Within each PNG all directions share one
# realization; different PNGs use different seeds.

resolve_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  existing <- candidates[file.exists(candidates) | dir.exists(candidates)]
  if (!length(existing)) {
    stop(sprintf("Could not resolve path: %s", file.path(...)))
  }
  existing[[1L]]
}

parse_arguments <- function(args) {
  values <- list(
    seeds = seq_len(10L),
    n_cores = 10L,
    output_dir = file.path(
      "output", "convergence", "gaussian_process", "s1",
      "seed_gallery_angles16_t200_tdomain_rainbow"
    ),
    cache_result = file.path(
      "output", "convergence", "gaussian_process", "s1",
      "limit_gaussian_s1_exact_kappa2_angles16_t200_tdomain_result.rds"
    )
  )

  for (argument in args) {
    if (startsWith(argument, "--seeds=")) {
      values$seeds <- as.integer(strsplit(sub("^--seeds=", "", argument), ",", fixed = TRUE)[[1L]])
    } else if (startsWith(argument, "--n-cores=")) {
      values$n_cores <- as.integer(sub("^--n-cores=", "", argument))
    } else if (startsWith(argument, "--output-dir=")) {
      values$output_dir <- sub("^--output-dir=", "", argument)
    } else if (startsWith(argument, "--cache-result=")) {
      values$cache_result <- sub("^--cache-result=", "", argument)
    } else {
      stop(sprintf("Unknown argument: %s", argument))
    }
  }

  values$seeds <- unique(values$seeds)
  if (length(values$seeds) != 10L || any(!is.finite(values$seeds))) {
    stop("`--seeds` must contain exactly 10 distinct integer seeds.")
  }
  if (length(values$n_cores) != 1L || !is.finite(values$n_cores) || values$n_cores < 1L) {
    stop("`--n-cores` must be a positive integer.")
  }
  values$n_cores <- min(values$n_cores, length(values$seeds))
  values
}

validate_cache <- function(result) {
  required <- c("cov_matrix", "chol_upper", "curve_data", "circle_grid", "omega_grid", "t_grid")
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop(sprintf("The cached result is missing: %s", paste(missing, collapse = ", ")))
  }
  if (nrow(result$omega_grid) != 16L || length(result$t_grid) != 200L || result$plot_domain != "t") {
    stop("The cached result must use 16 directions, 200 t-points, and `plot_domain = 't'`.")
  }
  if (!identical(result$color_scheme, "rainbow")) {
    stop("The cached result must use `color_scheme = 'rainbow'`.")
  }
  invisible(result)
}

draw_realization <- function(seed, chol_upper) {
  set.seed(seed)
  standard_normal_draw <- stats::rnorm(nrow(chol_upper))
  list(
    seed = seed,
    standard_normal_draw = standard_normal_draw,
    realization_vec = as.numeric(t(chol_upper) %*% standard_normal_draw)
  )
}

build_panels_for_realization <- function(base_result, realization_vec) {
  n_omega <- nrow(base_result$omega_grid)
  g_inner <- matrix(realization_vec, nrow = n_omega, ncol = length(base_result$t_grid))
  curve_data <- base_result$curve_data

  for (i in seq_len(n_omega)) {
    row_indices <- which(curve_data$omega_id == i)
    curve_data$g[row_indices] <- stats::approx(
      x = c(0, base_result$t_grid, 2),
      y = c(0, g_inner[i, ], 0),
      xout = curve_data$t[row_indices],
      method = "linear",
      ties = "ordered",
      rule = 2
    )$y
  }

  max_abs_g <- max(abs(curve_data$g))
  amplitude_scale <- if (max_abs_g > 0) base_result$lateral_scale / max_abs_g else 0
  curve_data$x_circle <- (1 + base_result$ray_extension * curve_data$radial_progress) * cos(curve_data$theta) -
    amplitude_scale * curve_data$g * sin(curve_data$theta)
  curve_data$y_circle <- (1 + base_result$ray_extension * curve_data$radial_progress) * sin(curve_data$theta) +
    amplitude_scale * curve_data$g * cos(curve_data$theta)

  panels <- build_s1_panels_from_curve_data(
    curve_data = curve_data,
    circle_grid = base_result$circle_grid,
    ray_extension = base_result$ray_extension,
    right_x_label = base_result$right_x_label,
    right_y_label = base_result$right_y_label
  )
  list(curve_data = curve_data, panels = panels)
}

main <- function() {
  options(warn = 1)
  arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
  visualization_path <- resolve_path("convergence_empirical_process", "gaussian_process_s1_visualization.R")
  source(visualization_path)

  cache_path <- resolve_path(arguments$cache_result)
  base_result <- readRDS(cache_path)
  validate_cache(base_result)
  dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Forking shares the read-only Cholesky factor on macOS/Linux.  Limit BLAS
  # threads so the ten workers do not oversubscribe the machine.
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  seeds <- arguments$seeds
  draws <- if (.Platform$OS.type == "windows" || arguments$n_cores == 1L) {
    lapply(seeds, draw_realization, chol_upper = base_result$chol_upper)
  } else {
    parallel::mclapply(
      seeds,
      draw_realization,
      chol_upper = base_result$chol_upper,
      mc.cores = arguments$n_cores,
      mc.preschedule = FALSE
    )
  }

  output_rows <- vector("list", length(draws))
  for (i in seq_along(draws)) {
    draw <- draws[[i]]
    plot_data <- build_panels_for_realization(base_result, draw$realization_vec)
    plot_stem <- file.path(
      arguments$output_dir,
      sprintf("limit_gaussian_s1_exact_kappa2_angles16_t200_tdomain_rainbow_seed%04d.png", draw$seed)
    )
    paths <- write_s1_plot_files(plot_data$panels$left_panel, plot_data$panels$right_panel, plot_stem)
    output_rows[[i]] <- data.frame(
      seed = draw$seed,
      combined = unname(paths[["combined"]]),
      left_panel = unname(paths[["left"]]),
      right_panel = unname(paths[["right"]]),
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, output_rows)
  manifest_path <- file.path(arguments$output_dir, "manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE)
  message(sprintf("Generated %d plots in %s", nrow(manifest), normalizePath(arguments$output_dir)))
  message(sprintf("Manifest: %s", normalizePath(manifest_path)))
}

main()
