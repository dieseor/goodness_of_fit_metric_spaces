jp_logspace_sub <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  max_ab <- pmax(a, b)
  min_ab <- pmin(a, b)
  max_ab + log1p(-exp(min_ab - max_ab))
}

jp_parallel_lapply <- function(X, FUN, n_cores = 1L, ...) {
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

normalize_jp_s2_data <- function(data, tol = 1e-8) {
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

compute_jp_projection_direction <- function(data,
                                            resultant_tol = 1e-8,
                                            fallback = c(0, 0, 1),
                                            warn_on_fallback = TRUE) {
  data <- normalize_jp_s2_data(data)
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

project_jp_data <- function(data, direction) {
  data <- normalize_jp_s2_data(data)
  direction <- as.numeric(direction)

  if (length(direction) != 3L || any(!is.finite(direction))) {
    stop("`direction` must be a finite vector of length 3.")
  }

  direction <- direction / sqrt(sum(direction^2))
  pmin(pmax(as.numeric(data %*% direction), -1), 1)
}

jp_projected_log_integral_z <- function(alpha,
                                        beta,
                                        alpha_tol = 1e-12,
                                        beta_minus_one_tol = 1e-10) {
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)

  if (length(alpha) != 1L || length(beta) != 1L || !is.finite(alpha) || !is.finite(beta)) {
    stop("`alpha` and `beta` must be finite scalars.")
  }
  if (abs(alpha) >= 1) {
    stop("`alpha` must lie in (-1, 1).")
  }

  if (abs(alpha) <= alpha_tol) {
    return(log(2))
  }

  if (abs(beta + 1) <= beta_minus_one_tol) {
    log_ratio <- log1p(alpha) - log1p(-alpha)
    return(log(abs(log_ratio)) - log(abs(alpha)))
  }

  gamma <- beta + 1
  a_term <- gamma * log1p(alpha)
  b_term <- gamma * log1p(-alpha)

  jp_logspace_sub(a_term, b_term) - log(abs(alpha * gamma))
}

jp_projected_log_integral_vmf_z <- function(kappa, kappa_tol = 1e-12) {
  kappa <- as.numeric(kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite nonnegative scalar.")
  }

  if (abs(kappa) <= kappa_tol) {
    return(log(2))
  }

  log(expm1(2 * kappa)) - kappa - log(kappa)
}

jp_projected_log_density_z <- function(z,
                                       kappa,
                                       psi,
                                       alpha_tol = 1e-12,
                                       kappa_tol = 1e-12) {
  z <- as.numeric(z)
  out <- rep(-Inf, length(z))
  valid <- is.finite(z) & z >= -1 & z <= 1

  if (!any(valid)) {
    return(out)
  }
  if (length(kappa) != 1L || length(psi) != 1L || !is.finite(kappa) || !is.finite(psi) || kappa < 0) {
    stop("`kappa` and `psi` must be finite scalars with `kappa >= 0`.")
  }

  if (abs(psi) <= alpha_tol) {
    if (abs(kappa) <= kappa_tol) {
      out[valid] <- -log(2)
      return(out)
    }

    log_I0 <- jp_projected_log_integral_vmf_z(kappa)
    out[valid] <- kappa * z[valid] - log_I0
    return(out)
  }

  alpha <- tanh(kappa * psi)
  beta <- 1 / psi

  if (abs(alpha) <= alpha_tol) {
    out[valid] <- -log(2)
    return(out)
  }

  term <- 1 + alpha * z[valid]
  if (any(term <= 0) || any(!is.finite(term))) {
    return(out)
  }

  log_I <- jp_projected_log_integral_z(alpha = alpha, beta = beta)
  out[valid] <- beta * log(term) - log_I
  out
}

jp_projected_density_z <- function(z, kappa, psi, log = FALSE) {
  log_density <- jp_projected_log_density_z(z = z, kappa = kappa, psi = psi)
  if (isTRUE(log)) log_density else exp(log_density)
}

jp_projected_cdf_z <- function(z,
                               kappa,
                               psi,
                               alpha_tol = 1e-12,
                               kappa_tol = 1e-12,
                               beta_minus_one_tol = 1e-10) {
  z <- as.numeric(z)
  out <- rep(NA_real_, length(z))

  below <- is.finite(z) & z < -1
  above <- is.finite(z) & z > 1
  valid <- is.finite(z) & z >= -1 & z <= 1

  out[below] <- 0
  out[above] <- 1

  if (!any(valid)) {
    return(out)
  }
  if (length(kappa) != 1L || length(psi) != 1L || !is.finite(kappa) || !is.finite(psi) || kappa < 0) {
    stop("`kappa` and `psi` must be finite scalars with `kappa >= 0`.")
  }

  z_valid <- z[valid]

  if (abs(psi) <= alpha_tol) {
    if (abs(kappa) <= kappa_tol) {
      out[valid] <- (z_valid + 1) / 2
      return(pmin(pmax(out, 0), 1))
    }

    denominator <- expm1(2 * kappa)
    numerator <- expm1(kappa * (z_valid + 1))
    out[valid] <- numerator / denominator
    return(pmin(pmax(out, 0), 1))
  }

  alpha <- tanh(kappa * psi)
  beta <- 1 / psi

  if (abs(alpha) <= alpha_tol) {
    out[valid] <- (z_valid + 1) / 2
    return(pmin(pmax(out, 0), 1))
  }

  if (abs(beta + 1) <= beta_minus_one_tol) {
    numerator <- log1p(alpha * z_valid) - log1p(-alpha)
    denominator <- log1p(alpha) - log1p(-alpha)
    out[valid] <- numerator / denominator
    return(pmin(pmax(out, 0), 1))
  }

  gamma <- beta + 1
  lower_term <- gamma * log1p(-alpha)
  upper_term <- gamma * log1p(alpha)
  z_term <- gamma * log1p(alpha * z_valid)
  log_num <- jp_logspace_sub(z_term, lower_term)
  log_den <- jp_logspace_sub(upper_term, lower_term)
  out[valid] <- exp(log_num - log_den)
  pmin(pmax(out, 0), 1)
}

fit_projected_kappa_given_psi <- function(z,
                                          psi,
                                          kappa_upper_vmf = 100,
                                          alpha_tol = 1e-12) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]

  if (length(z) == 0L) {
    stop("`z` must contain at least one finite projection.")
  }
  if (any(z < -1 | z > 1)) {
    stop("`z` must lie in [-1, 1].")
  }
  if (length(psi) != 1L || !is.finite(psi)) {
    stop("`psi` must be a finite scalar.")
  }

  kappa_upper <- if (abs(psi) <= alpha_tol) kappa_upper_vmf else 6 / abs(psi)
  kappa_upper <- as.numeric(kappa_upper)

  if (!is.finite(kappa_upper) || kappa_upper <= 0) {
    stop("Computed `kappa_upper` is not positive and finite.")
  }

  objective <- function(kappa) {
    sum(jp_projected_log_density_z(z = z, kappa = kappa, psi = psi))
  }

  opt <- try(
    stats::optimize(f = objective, interval = c(0, kappa_upper), maximum = TRUE),
    silent = TRUE
  )

  if (inherits(opt, "try-error") || !is.list(opt) || !is.finite(opt$objective) || !is.finite(opt$maximum)) {
    return(data.frame(
      psi = psi,
      kappa_hat = NA_real_,
      loglik_projected = NA_real_,
      status = "optimize_failed",
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    psi = psi,
    kappa_hat = opt$maximum,
    loglik_projected = opt$objective,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

plot_jp_projected_diagnostic <- function(z,
                                         cdf_grid,
                                         fits,
                                         sample_label,
                                         output_path,
                                         legend_position = "topleft",
                                         hat_curve = NULL) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  palette_values <- grDevices::hcl.colors(max(nrow(fits), 3L), palette = "Dark 3")
  line_colors <- palette_values[seq_len(nrow(fits))]

  grDevices::png(filename = output_path, width = 1200, height = 900, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(
    stats::ecdf(z),
    do.points = FALSE,
    verticals = TRUE,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    main = sprintf("JP projected poor-man diagnostic: %s", sample_label),
    xlab = "z = m^T x",
    ylab = "CDF",
    lwd = 2,
    col = "black"
  )
  grid(col = "#d9d9d9")

  for (i in seq_len(nrow(fits))) {
    fit_i <- fits[i, , drop = FALSE]
    cdf_i <- cdf_grid[cdf_grid$psi == fit_i$psi, , drop = FALSE]
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
    "psi=%s, kappa=%.3f, ll=%.3f",
    format(fits$psi, trim = TRUE, scientific = FALSE),
    fits$kappa_hat,
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

run_jp_projected_poor_mans_diagnostic <- function(data,
                                                  psi_grid = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
                                                  output_dir = file.path("output", "distance_profiles", "jp", "projected_poor_mans"),
                                                  sample_label = "sample",
                                                  resultant_tol = 1e-8,
                                                  fallback = c(0, 0, 1),
                                                  kappa_upper_vmf = 100,
                                                  z_grid_size = 401L,
                                                  n_cores = 1L,
                                                  hat_params = NULL,
                                                  legend_position = "topleft",
                                                  save_plot = TRUE,
                                                  verbose = TRUE) {
  data <- normalize_jp_s2_data(data)
  psi_grid <- unique(as.numeric(psi_grid))
  psi_grid <- psi_grid[is.finite(psi_grid)]

  if (length(psi_grid) == 0L) {
    stop("`psi_grid` must contain at least one finite value.")
  }
  if (length(sample_label) != 1L || !nzchar(sample_label)) {
    stop("`sample_label` must be a non-empty string.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  direction_info <- compute_jp_projection_direction(
    data = data,
    resultant_tol = resultant_tol,
    fallback = fallback
  )
  z <- project_jp_data(data = data, direction = direction_info$direction)
  empirical_cdf <- stats::ecdf(z)

  fits <- do.call(
    rbind,
    jp_parallel_lapply(psi_grid, function(psi) {
      fit_projected_kappa_given_psi(
        z = z,
        psi = psi,
        kappa_upper_vmf = kappa_upper_vmf
      )
    }, n_cores = n_cores)
  )

  fits$ks_distance <- vapply(seq_len(nrow(fits)), function(i) {
    if (!is.finite(fits$kappa_hat[i])) {
      return(NA_real_)
    }

    z_grid <- seq(-1, 1, length.out = max(as.integer(z_grid_size), 101L))
    max(abs(empirical_cdf(z_grid) - jp_projected_cdf_z(z_grid, kappa = fits$kappa_hat[i], psi = fits$psi[i])))
  }, numeric(1))

  fits <- fits[order(-fits$loglik_projected, fits$psi), , drop = FALSE]
  rownames(fits) <- NULL

  z_grid <- seq(-1, 1, length.out = max(as.integer(z_grid_size), 101L))
  hat_curve <- NULL
  if (!is.null(hat_params)) {
    if (is.list(hat_params)) {
      kappa_hat <- as.numeric(hat_params$kappa)
      psi_hat <- as.numeric(hat_params$psi)
      hat_label <- hat_params$label
      hat_color <- hat_params$color
    } else if (is.numeric(hat_params) && length(hat_params) >= 2L) {
      kappa_hat <- as.numeric(hat_params[[1L]])
      psi_hat <- as.numeric(hat_params[[2L]])
      hat_label <- NULL
      hat_color <- NULL
    } else {
      stop("`hat_params` must be a list with `kappa` and `psi`, or a numeric vector of length >= 2.")
    }

    if (!is.finite(kappa_hat) || !is.finite(psi_hat) || kappa_hat < 0) {
      stop("`hat_params` contains invalid `kappa`/`psi` values.")
    }

    if (is.null(hat_label) || !nzchar(as.character(hat_label)[[1L]])) {
      hat_label <- sprintf("hat: psi=%.4f, kappa=%.3f", psi_hat, kappa_hat)
    }
    if (is.null(hat_color) || !nzchar(as.character(hat_color)[[1L]])) {
      hat_color <- "#b30000"
    }

    hat_curve <- list(
      z_grid = z_grid,
      cdf = jp_projected_cdf_z(z_grid, kappa = kappa_hat, psi = psi_hat),
      label = as.character(hat_label)[[1L]],
      color = as.character(hat_color)[[1L]],
      lwd = 3,
      lty = 2,
      kappa_hat = kappa_hat,
      psi_hat = psi_hat
    )
  }

  cdf_grid <- do.call(
    rbind,
    lapply(seq_len(nrow(fits)), function(i) {
      data.frame(
        psi = fits$psi[i],
        kappa_hat = fits$kappa_hat[i],
        z_grid = z_grid,
        cdf = jp_projected_cdf_z(z_grid, kappa = fits$kappa_hat[i], psi = fits$psi[i]),
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

  fits_path <- file.path(output_dir, sprintf("%s_projected_kappa_fits.csv", sample_label))
  cdf_path <- file.path(output_dir, sprintf("%s_projected_cdf_grid.csv", sample_label))
  z_path <- file.path(output_dir, sprintf("%s_projected_data.csv", sample_label))
  direction_path <- file.path(output_dir, sprintf("%s_projection_direction.csv", sample_label))
  plot_path <- file.path(output_dir, sprintf("%s_projected_ecdf_overlay.png", sample_label))

  utils::write.csv(fits, fits_path, row.names = FALSE)
  utils::write.csv(cdf_grid, cdf_path, row.names = FALSE)
  utils::write.csv(projection_df, z_path, row.names = FALSE)
  utils::write.csv(direction_df, direction_path, row.names = FALSE)

  if (isTRUE(save_plot)) {
    plot_jp_projected_diagnostic(
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
    message(sprintf("Projected JP diagnostic written to %s", output_dir))
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

run_jp_projected_poor_mans_diagnostic_h0 <- function(H0,
                                                     n_values,
                                                     psi_grid = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
                                                     output_dir = file.path("output", "distance_profiles", "jp", "projected_poor_mans", "h0"),
                                                     n_cores = 1L,
                                                     n_cores_inner = 1L,
                                                     seed = NULL,
                                                     ...) {
  if (!is.function(H0)) {
    stop("`H0` must be a function that takes `n` and returns an n x 3 sample on S^2.")
  }

  n_values <- as.integer(n_values)
  if (length(n_values) == 0L || any(!is.finite(n_values)) || any(n_values < 1L)) {
    stop("`n_values` must contain positive integers.")
  }
  n_cores <- as.integer(n_cores[[1L]])
  n_cores_inner <- as.integer(n_cores_inner[[1L]])
  if (!is.finite(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }
  if (!is.finite(n_cores_inner) || n_cores_inner < 1L) {
    stop("`n_cores_inner` must be a positive integer.")
  }
  if (n_cores > 1L && n_cores_inner > 1L) {
    stop("Avoid nested parallelism: use either `n_cores > 1` across `n_values` or `n_cores_inner > 1` within each sample, but not both.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  results <- vector("list", length(n_values))
  result_list <- jp_parallel_lapply(n_values, function(n_i) {
    sample_i <- H0(n_i)
    run_jp_projected_poor_mans_diagnostic(
      data = sample_i,
      psi_grid = psi_grid,
      output_dir = output_dir,
      sample_label = paste0("n_", n_i),
      n_cores = n_cores_inner,
      ...
    )
  }, n_cores = n_cores)
  names(result_list) <- paste0("n_", n_values)
  results[] <- result_list

  summary_df <- do.call(
    rbind,
    lapply(results, function(res) {
      best_fit <- res$fits[1, , drop = FALSE]
      data.frame(
        sample_label = res$sample_label,
        n = res$n,
        best_psi = best_fit$psi,
        best_kappa_hat = best_fit$kappa_hat,
        best_loglik_projected = best_fit$loglik_projected,
        best_ks_distance = best_fit$ks_distance,
        stringsAsFactors = FALSE
      )
    })
  )

  summary_path <- file.path(output_dir, "projected_poor_mans_summary.csv")
  utils::write.csv(summary_df, summary_path, row.names = FALSE)

  list(
    results = results,
    summary = summary_df,
    summary_path = summary_path
  )
}

# Example for a concrete sample:
# source("distance_profiles/jp_projected_poor_mans_diagnostic.R")
# source("utils.R")
# x <- r_sph_jp(n = 200, mu = c(0, 0, 1), kappa = 1, psi = 0.5)
# res <- run_jp_projected_poor_mans_diagnostic(
#   data = x,
#   psi_grid = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
#   output_dir = file.path("output", "jp_projected_poor_mans_example"),
#   sample_label = "jp_sample",
#   n_cores = 12
# )

# Example under H0 for several n:
# source("distance_profiles/jp_projected_poor_mans_diagnostic.R")
# source("utils.R")
# H0 <- function(n) r_sph_jp(n = n, mu = c(0, 0, 1), kappa = 1, psi = 0.5)
# res_h0 <- run_jp_projected_poor_mans_diagnostic_h0(
#   H0 = H0,
#   n_values = c(50, 100, 200),
#   psi_grid = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
#   output_dir = file.path("output", "jp_projected_poor_mans_h0_example"),
#   n_cores = 3,
#   n_cores_inner = 1
# )