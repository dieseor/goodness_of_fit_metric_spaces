sc_parallel_lapply <- function(X, FUN, n_cores = 1L, ...) {
  n_cores <- as.integer(n_cores[[1L]])
  if (!is.finite(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }

  if (n_cores <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }

  if (.Platform$OS.type == "windows") {
    warning("`n_cores > 1` requested on Windows; falling back to sequential lapply().", call. = FALSE)
    return(lapply(X, FUN, ...))
  }

  parallel::mclapply(X, FUN, ..., mc.cores = n_cores)
}

sc_axis_projected_cdf <- get("spherical_cauchy_axis_projected_cdf", mode = "function")
sc_axis_projected_density <- get("spherical_cauchy_axis_projected_density", mode = "function")

normalize_sc_s2_data <- function(data, tol = 1e-8) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) != 3L) {
    stop("`data` must be a non-empty n x 3 matrix.")
  }
  if (any(!is.finite(data))) {
    stop("`data` must contain only finite values.")
  }

  row_norms <- sqrt(rowSums(data^2))
  if (any(row_norms <= tol)) {
    stop("`data` contains rows with near-zero norm.")
  }

  data_unit <- data / row_norms
  if (any(abs(sqrt(rowSums(data_unit^2)) - 1) > 10 * tol)) {
    stop("Failed to normalize `data` onto S^2.")
  }

  data_unit
}

compute_sc_projection_direction <- function(data,
                                            resultant_tol = 1e-8,
                                            fallback = c(0, 0, 1),
                                            warn_on_fallback = TRUE) {
  data <- normalize_sc_s2_data(data)
  fallback <- as.numeric(fallback)

  if (length(fallback) != 3L || any(!is.finite(fallback))) {
    stop("`fallback` must be a finite vector of length 3.")
  }

  fallback_norm <- sqrt(sum(fallback^2))
  if (fallback_norm <= 0) {
    stop("`fallback` must have positive norm.")
  }
  fallback <- fallback / fallback_norm

  resultant <- colMeans(data)
  resultant_norm <- sqrt(sum(resultant^2))

  if (is.finite(resultant_norm) && resultant_norm > resultant_tol) {
    direction <- resultant / resultant_norm
  } else {
    direction <- fallback
    if (isTRUE(warn_on_fallback)) {
      warning(
        sprintf(
          "Resultant norm %.3e is not above tolerance %.3e; using fallback direction.",
          resultant_norm,
          resultant_tol
        ),
        call. = FALSE
      )
    }
  }

  list(
    direction = direction,
    resultant = resultant,
    resultant_norm = resultant_norm
  )
}

project_sc_data <- function(data, direction) {
  data <- normalize_sc_s2_data(data)
  direction <- as.numeric(direction)

  if (length(direction) != 3L || any(!is.finite(direction))) {
    stop("`direction` must be a finite vector of length 3.")
  }

  direction <- direction / sqrt(sum(direction^2))
  pmin(pmax(as.numeric(data %*% direction), -1), 1)
}

fit_projected_sc_rho <- function(z,
                                 rho_lower = 1e-8,
                                 rho_upper = 0.999,
                                 tol = 1e-8) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  z <- pmin(pmax(z, -1), 1)

  if (length(z) == 0L) {
    stop("`z` must contain at least one finite projected value.")
  }
  if (!is.finite(rho_lower) || !is.finite(rho_upper) || rho_lower < 0 || rho_upper >= 1 || rho_lower >= rho_upper) {
    stop("Invalid rho optimization interval.")
  }

  neg_loglik <- function(rho) {
    -sum(sc_axis_projected_density(z, rho = rho, log = TRUE))
  }

  opt <- stats::optimize(
    f = neg_loglik,
    interval = c(rho_lower, rho_upper),
    tol = tol
  )

  data.frame(
    rho = opt$minimum,
    loglik_projected = -opt$objective,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

score_projected_sc_grid <- function(z,
                                    rho_grid,
                                    n_cores = 1L) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  z <- pmin(pmax(z, -1), 1)

  rho_grid <- unique(as.numeric(rho_grid))
  rho_grid <- rho_grid[is.finite(rho_grid) & rho_grid >= 0 & rho_grid < 1]
  if (length(rho_grid) == 0L) {
    stop("`rho_grid` must contain at least one finite value in [0, 1).")
  }

  empirical_cdf <- stats::ecdf(z)
  z_grid <- seq(-1, 1, length.out = 401L)

  fits <- do.call(
    rbind,
    sc_parallel_lapply(rho_grid, function(rho) {
      loglik <- sum(sc_axis_projected_density(z, rho = rho, log = TRUE))
      ks_distance <- max(abs(empirical_cdf(z_grid) - sc_axis_projected_cdf(z_grid, rho = rho)))
      data.frame(
        rho = rho,
        loglik_projected = loglik,
        status = "ok",
        ks_distance = ks_distance,
        stringsAsFactors = FALSE
      )
    }, n_cores = n_cores)
  )

  fits[order(-fits$loglik_projected, fits$rho), , drop = FALSE]
}

plot_sc_projected_diagnostic <- function(z,
                                         cdf_grid,
                                         fits,
                                         sample_label,
                                         output_path = NULL,
                                         legend_position = "topleft",
                                         hat_curve = NULL,
                                         width = 1200,
                                         height = 900,
                                         res = 140) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  empirical_cdf <- stats::ecdf(z)

  if (!is.null(output_path)) {
    grDevices::png(output_path, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  line_colors <- grDevices::hcl.colors(max(nrow(fits), 3L), palette = "Dark 3")[seq_len(nrow(fits))]

  plot(
    NA,
    NA,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    xlab = "Projected coordinate z = <x, omega>",
    ylab = "CDF",
    main = sprintf("Spherical Cauchy projected diagnostic (%s)", sample_label)
  )

  z_emp <- sort(unique(c(-1, z, 1)))
  lines(
    z_emp,
    empirical_cdf(z_emp),
    type = "s",
    lwd = 2,
    col = "black"
  )
  grid(col = "#d9d9d9")

  for (i in seq_len(nrow(fits))) {
    fit_i <- fits[i, , drop = FALSE]
    cdf_i <- cdf_grid[cdf_grid$rho == fit_i$rho, , drop = FALSE]
    lines(cdf_i$z_grid, cdf_i$cdf, col = line_colors[i], lwd = 2)
  }

  if (!is.null(hat_curve)) {
    lines(
      hat_curve$z_grid,
      hat_curve$cdf,
      col = hat_curve$color,
      lwd = hat_curve$lwd,
      lty = hat_curve$lty
    )
  }

  legend_labels <- sprintf(
    "rho=%s, ll=%.3f",
    format(fits$rho, trim = TRUE, scientific = FALSE),
    fits$loglik_projected
  )

  if (!is.null(hat_curve)) {
    legend_labels <- c(legend_labels, hat_curve$label)
    legend_colors <- c("black", line_colors, hat_curve$color)
    legend_lty <- c(1, rep(1, nrow(fits)), hat_curve$lty)
    legend_lwd <- c(2, rep(2, nrow(fits)), hat_curve$lwd)
  } else {
    legend_colors <- c("black", line_colors)
    legend_lty <- c(1, rep(1, nrow(fits)))
    legend_lwd <- c(2, rep(2, nrow(fits)))
  }

  legend(
    legend_position,
    legend = c("ECDF", legend_labels),
    col = legend_colors,
    lty = legend_lty,
    lwd = legend_lwd,
    bty = "n",
    cex = 0.9
  )
}

run_sc_projected_poor_mans_diagnostic <- function(data,
                                                  rho_grid = c(0, 0.25, 0.5, 0.7, 0.85, 0.93),
                                                  output_dir = file.path("output", "spherical_cauchy_projected_poor_mans_diagnostic"),
                                                  sample_label = "sample",
                                                  resultant_tol = 1e-8,
                                                  fallback = c(0, 0, 1),
                                                  z_grid_size = 401L,
                                                  n_cores = 1L,
                                                  add_rho_hat = TRUE,
                                                  rho_fit_lower = 1e-8,
                                                  rho_fit_upper = 0.999,
                                                  legend_position = "topleft",
                                                  save_plot = TRUE,
                                                  verbose = TRUE) {
  data <- normalize_sc_s2_data(data)
  rho_grid <- unique(as.numeric(rho_grid))
  rho_grid <- rho_grid[is.finite(rho_grid) & rho_grid >= 0 & rho_grid < 1]

  if (length(rho_grid) == 0L) {
    stop("`rho_grid` must contain at least one finite value in [0, 1).")
  }
  if (length(sample_label) != 1L || !nzchar(sample_label)) {
    stop("`sample_label` must be a non-empty string.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  direction_info <- compute_sc_projection_direction(
    data = data,
    resultant_tol = resultant_tol,
    fallback = fallback
  )
  z <- project_sc_data(data = data, direction = direction_info$direction)
  fits <- score_projected_sc_grid(z = z, rho_grid = rho_grid, n_cores = n_cores)

  z_grid <- seq(-1, 1, length.out = max(as.integer(z_grid_size), 101L))
  hat_curve <- NULL
  if (isTRUE(add_rho_hat)) {
    mle_fit <- fit_projected_sc_rho(
      z = z,
      rho_lower = rho_fit_lower,
      rho_upper = rho_fit_upper
    )

    hat_curve <- list(
      z_grid = z_grid,
      cdf = sc_axis_projected_cdf(z_grid, rho = mle_fit$rho),
      label = sprintf("rho_hat=%.4f", mle_fit$rho),
      color = "#b30000",
      lwd = 3,
      lty = 2,
      rho_hat = mle_fit$rho,
      loglik_hat = mle_fit$loglik_projected
    )

    fits <- rbind(
      data.frame(
        rho = mle_fit$rho,
        loglik_projected = mle_fit$loglik_projected,
        status = "mle",
        ks_distance = max(abs(stats::ecdf(z)(z_grid) - sc_axis_projected_cdf(z_grid, rho = mle_fit$rho))),
        stringsAsFactors = FALSE
      ),
      fits
    )
    fits <- fits[order(-fits$loglik_projected, fits$rho), , drop = FALSE]
    rownames(fits) <- NULL
  }

  cdf_grid <- do.call(
    rbind,
    lapply(seq_len(nrow(fits)), function(i) {
      data.frame(
        rho = fits$rho[i],
        z_grid = z_grid,
        cdf = sc_axis_projected_cdf(z_grid, rho = fits$rho[i]),
        stringsAsFactors = FALSE
      )
    })
  )

  projection_df <- data.frame(
    index = seq_along(z),
    z = z,
    stringsAsFactors = FALSE
  )

  direction_df <- data.frame(
    sample_label = sample_label,
    n = nrow(data),
    resultant_norm = direction_info$resultant_norm,
    resultant_1 = direction_info$resultant[1],
    resultant_2 = direction_info$resultant[2],
    resultant_3 = direction_info$resultant[3],
    direction_1 = direction_info$direction[1],
    direction_2 = direction_info$direction[2],
    direction_3 = direction_info$direction[3],
    stringsAsFactors = FALSE
  )

  fits_path <- file.path(output_dir, sprintf("%s_projected_rho_fits.csv", sample_label))
  cdf_path <- file.path(output_dir, sprintf("%s_projected_cdf_grid.csv", sample_label))
  z_path <- file.path(output_dir, sprintf("%s_projected_data.csv", sample_label))
  direction_path <- file.path(output_dir, sprintf("%s_projection_direction.csv", sample_label))
  plot_path <- file.path(output_dir, sprintf("%s_projected_ecdf_overlay.png", sample_label))

  utils::write.csv(fits, fits_path, row.names = FALSE)
  utils::write.csv(cdf_grid, cdf_path, row.names = FALSE)
  utils::write.csv(projection_df, z_path, row.names = FALSE)
  utils::write.csv(direction_df, direction_path, row.names = FALSE)

  if (isTRUE(save_plot)) {
    plot_sc_projected_diagnostic(
      z = z,
      cdf_grid = cdf_grid,
      fits = fits,
      sample_label = sample_label,
      output_path = plot_path,
      legend_position = legend_position,
      hat_curve = hat_curve
    )
  }

  if (isTRUE(verbose)) {
    message(sprintf("Projected spherical Cauchy diagnostic written to %s", output_dir))
  }

  list(
    sample_label = sample_label,
    n = nrow(data),
    direction_info = direction_info,
    z = z,
    fits = fits,
    cdf_grid = cdf_grid,
    hat_curve = hat_curve,
    output_paths = list(
      fits = fits_path,
      cdf_grid = cdf_path,
      projected_data = z_path,
      direction = direction_path,
      plot = if (isTRUE(save_plot)) plot_path else NA_character_
    )
  )
}