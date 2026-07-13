suppressPackageStartupMessages({
  source(file.path("utils.R"))
  source(file.path("real_data", "comets", "utils_comets_data.R"))
  source(file.path("scripts", "path_helpers.R"))
  source(file.path("bootstrap", "cardioid_model_spec.R"))
})

normalize_projection_vector <- function(x,
                                        arg_name = "`x`",
                                        tol = 1e-12) {
  x <- as.numeric(x)
  if (length(x) != 3L || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite vector of length 3.", arg_name))
  }
  norm_x <- sqrt(sum(x^2))
  if (!is.finite(norm_x) || norm_x <= tol) {
    stop(sprintf("%s must have norm strictly larger than %.3e.", arg_name, tol))
  }
  x / norm_x
}

normalize_projection_matrix <- function(x) {
  jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
}

monotone_clip_cdf <- function(z_grid, cdf_values) {
  ord <- order(z_grid)
  z_grid <- as.numeric(z_grid)[ord]
  cdf_values <- as.numeric(cdf_values)[ord]
  cdf_values <- pmin(pmax(cdf_values, 0), 1)
  cdf_values <- cummax(cdf_values)
  cdf_values
}

rbind_fill_projected_summary <- function(rows) {
  if (length(rows) == 0L) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned_rows <- lapply(rows, function(row) {
    missing_names <- setdiff(all_names, names(row))
    for (missing_name in missing_names) {
      row[[missing_name]] <- NA_real_
    }
    row[, all_names, drop = FALSE]
  })

  do.call(rbind, aligned_rows)
}

canonicalize_eigen_axis_sign <- function(axis,
                                         sample_mean_direction,
                                         tol = 1e-12) {
  axis <- normalize_projection_vector(axis, arg_name = "`axis`", tol = tol)
  sample_mean_direction <- normalize_projection_vector(
    sample_mean_direction,
    arg_name = "`sample_mean_direction`",
    tol = tol
  )
  alignment <- sum(axis * sample_mean_direction)

  if (is.finite(alignment) && abs(alignment) > tol) {
    if (alignment < 0) {
      axis <- -axis
      alignment <- -alignment
    }
    sign_rule <- "aligned_with_sample_mean"
  } else if (axis[[3L]] < 0) {
    axis <- -axis
    sign_rule <- "nonnegative_third_coordinate"
  } else {
    sign_rule <- "nonnegative_third_coordinate"
  }

  list(
    axis = axis,
    alignment_with_sample_mean = alignment,
    sign_rule = sign_rule
  )
}

compute_director_projection_axis <- function(data_matrix,
                                             dataset_label,
                                             tol = 1e-12) {
  data_matrix <- normalize_projection_matrix(data_matrix)
  sample_mean <- colMeans(data_matrix)
  mean_norm <- sqrt(sum(sample_mean^2))
  if (!is.finite(mean_norm) || mean_norm <= tol) {
    stop(sprintf("The sample mean is numerically zero for dataset `%s`.", dataset_label))
  }
  sample_mean_direction <- sample_mean / mean_norm

  scatter_matrix <- crossprod(data_matrix)
  eig <- eigen(scatter_matrix, symmetric = TRUE)
  lead_axis_signed <- canonicalize_eigen_axis_sign(
    axis = eig$vectors[, 1L],
    sample_mean_direction = sample_mean_direction,
    tol = tol
  )

  if (identical(dataset_label, "short_period")) {
    axis <- sample_mean_direction
    axis_type <- "sample_mean_direction"
    sign_rule <- "sample_mean_direction"
    alignment <- 1
  } else if (identical(dataset_label, "long_period")) {
    axis <- lead_axis_signed$axis
    axis_type <- "scatter_first_eigenvector"
    sign_rule <- lead_axis_signed$sign_rule
    alignment <- lead_axis_signed$alignment_with_sample_mean
  } else {
    stop("`dataset_label` must be either `short_period` or `long_period`.")
  }

  list(
    axis = axis,
    axis_type = axis_type,
    sign_rule = sign_rule,
    alignment_with_sample_mean = alignment,
    sample_mean = sample_mean,
    sample_mean_direction = sample_mean_direction,
    sample_mean_norm = mean_norm,
    scatter_matrix = scatter_matrix,
    scatter_eigenvalues = eig$values,
    scatter_first_eigenvector = lead_axis_signed$axis
  )
}

project_data_onto_axis <- function(data_matrix, axis) {
  data_matrix <- normalize_projection_matrix(data_matrix)
  axis <- normalize_projection_vector(axis, arg_name = "`axis`")
  pmin(pmax(as.numeric(data_matrix %*% axis), -1), 1)
}

build_projected_fit <- function(z_observed,
                                cdf_fun,
                                grid_size = 1001L) {
  z_observed <- as.numeric(z_observed)
  z_observed <- z_observed[is.finite(z_observed)]
  z_observed <- pmin(pmax(z_observed, -1), 1)
  if (length(z_observed) == 0L) {
    stop("`z_observed` must contain at least one finite projected value.")
  }

  ecdf_z <- stats::ecdf(z_observed)
  z_grid <- sort(unique(c(seq(-1, 1, length.out = as.integer(grid_size)), z_observed)))
  fitted_cdf_raw <- as.numeric(cdf_fun(z_grid))
  fitted_cdf <- monotone_clip_cdf(z_grid = z_grid, cdf_values = fitted_cdf_raw)
  empirical_cdf <- ecdf_z(z_grid)

  list(
    z_grid = z_grid,
    fitted_cdf = fitted_cdf,
    empirical_cdf = empirical_cdf,
    ks = max(abs(empirical_cdf - fitted_cdf)),
    cvm = mean((empirical_cdf - fitted_cdf)^2)
  )
}

make_axis_metadata_row <- function(dataset_label, axis_info) {
  data.frame(
    dataset = dataset_label,
    axis_type = axis_info$axis_type,
    sign_rule = axis_info$sign_rule,
    axis_1 = axis_info$axis[[1L]],
    axis_2 = axis_info$axis[[2L]],
    axis_3 = axis_info$axis[[3L]],
    sample_mean_1 = axis_info$sample_mean[[1L]],
    sample_mean_2 = axis_info$sample_mean[[2L]],
    sample_mean_3 = axis_info$sample_mean[[3L]],
    sample_mean_direction_1 = axis_info$sample_mean_direction[[1L]],
    sample_mean_direction_2 = axis_info$sample_mean_direction[[2L]],
    sample_mean_direction_3 = axis_info$sample_mean_direction[[3L]],
    sample_mean_norm = axis_info$sample_mean_norm,
    alignment_with_sample_mean = axis_info$alignment_with_sample_mean,
    scatter_eigenvalue_1 = axis_info$scatter_eigenvalues[[1L]],
    scatter_eigenvalue_2 = axis_info$scatter_eigenvalues[[2L]],
    scatter_eigenvalue_3 = axis_info$scatter_eigenvalues[[3L]],
    scatter_first_eigenvector_1 = axis_info$scatter_first_eigenvector[[1L]],
    scatter_first_eigenvector_2 = axis_info$scatter_first_eigenvector[[2L]],
    scatter_first_eigenvector_3 = axis_info$scatter_first_eigenvector[[3L]],
    stringsAsFactors = FALSE
  )
}

plot_projected_ecdf_overlay <- function(z_observed,
                                        fit,
                                        output_path,
                                        dataset_label,
                                        model_label,
                                        axis_type) {
  z_observed <- as.numeric(z_observed)
  ecdf_z <- stats::ecdf(z_observed)

  grDevices::png(output_path, width = 1200, height = 900, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(
    NA,
    NA,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    xlab = "Projected coordinate z = <x, gamma>",
    ylab = "CDF",
    main = sprintf("%s projected ECDF overlay (%s)", dataset_label, model_label),
    sub = sprintf("Axis: %s", axis_type)
  )
  grid(col = "#d9d9d9")
  z_emp <- sort(unique(c(-1, z_observed, 1)))
  lines(z_emp, ecdf_z(z_emp), type = "s", lwd = 2.2, col = "black")
  lines(fit$z_grid, fit$fitted_cdf, lwd = 2.2, col = "#c23b23")
  legend(
    "topleft",
    legend = c("ECDF", "Fitted CDF"),
    col = c("black", "#c23b23"),
    lty = c(1, 1),
    lwd = c(2.2, 2.2),
    bty = "n"
  )
}

load_saved_spherical_cauchy_theta <- function(dataset_label,
                                              tol = 1e-10) {
  summary_path <- canonical_comets_spherical_cauchy_dir(
    "paper_results_B5000_sampleks",
    "fast",
    "comets_spherical_cauchy_short_long_fast_summary.csv"
  )
  summary_df <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  rows <- summary_df[summary_df$dataset == dataset_label, , drop = FALSE]
  if (nrow(rows) == 0L) {
    stop(sprintf("No spherical Cauchy rows found for dataset `%s`.", dataset_label))
  }

  mu_mat <- as.matrix(rows[, c("mu_hat_1", "mu_hat_2", "mu_hat_3"), drop = FALSE])
  rho_vals <- as.numeric(rows$rho_hat)
  if (max(abs(mu_mat - rep(mu_mat[1L, ], each = nrow(mu_mat)))) > tol ||
      max(abs(rho_vals - rho_vals[[1L]])) > tol) {
    stop(sprintf("Stored spherical Cauchy parameters disagree across statistics for `%s`.", dataset_label))
  }

  list(
    theta = list(
      mu = normalize_projection_vector(mu_mat[1L, ], arg_name = "`mu_hat`"),
      rho = rho_vals[[1L]]
    ),
    source_path = summary_path
  )
}

load_saved_small_circle_theta <- function(dataset_label,
                                          tol = 1e-10) {
  run_name_ks <- if (identical(dataset_label, "short_period")) {
    "short_comets_B5000_sampleks"
  } else {
    "long_comets_B5000_sampleks"
  }
  run_name_cvm <- if (identical(dataset_label, "short_period")) {
    "short_comets_B5000_samplecvm"
  } else {
    "long_comets_B5000_samplecvm"
  }

  ks_path <- canonical_comets_small_circle_dir(run_name_ks, "fast", "stage_01_M5000_B5000.rds")
  cvm_path <- canonical_comets_small_circle_dir(run_name_cvm, "fast", "stage_01_M5000_B5000.rds")
  ks_result <- readRDS(ks_path)
  cvm_result <- readRDS(cvm_path)
  ks_theta <- ks_result$observed$theta_hat
  cvm_theta <- cvm_result$observed$theta_hat

  mu_ks <- normalize_projection_vector(ks_theta$mu, arg_name = "`small_circle mu`")
  mu_cvm <- normalize_projection_vector(cvm_theta$mu, arg_name = "`small_circle mu`")
  if (max(abs(mu_ks - mu_cvm)) > tol ||
      abs(as.numeric(ks_theta$kappa) - as.numeric(cvm_theta$kappa)) > tol ||
      abs(as.numeric(ks_theta$nu) - as.numeric(cvm_theta$nu)) > tol) {
    stop(sprintf(
      "Stored small-circle parameters disagree between KS and CvM for dataset `%s`.",
      dataset_label
    ))
  }

  list(
    theta = list(
      mu = mu_ks,
      kappa = as.numeric(ks_theta$kappa),
      nu = as.numeric(ks_theta$nu),
      ambient_dim = as.integer(ks_theta$ambient_dim %||% 3L)
    ),
    source_path = paste(c(ks_path, cvm_path), collapse = ";")
  )
}

load_saved_theta_csv <- function(path,
                                 family) {
  theta_df <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(theta_df) != 1L) {
    stop(sprintf("Expected a single theta row for `%s` at `%s`.", family, path))
  }
  theta_df[[1L]]
  list(theta_df = theta_df, source_path = path)
}

load_saved_uniform_beta_theta <- function(dataset_label) {
  subdir <- if (identical(dataset_label, "short_period")) {
    "01_short_period_uniform_beta_mixture"
  } else {
    "02_long_period_uniform_beta_mixture"
  }
  theta_path <- canonical_comets_mixture_dir(
    "uniform_beta_mixture_short_long_B5000",
    "fast",
    subdir,
    "theta_hat.csv"
  )
  theta_df <- utils::read.csv(theta_path, stringsAsFactors = FALSE)
  list(
    theta = list(
      mu = normalize_projection_vector(theta_df[1L, c("mu_1", "mu_2", "mu_3")], arg_name = "`mu`"),
      weight_uniform = theta_df$weight_uniform[[1L]],
      alpha = theta_df$alpha[[1L]],
      beta = theta_df$beta[[1L]]
    ),
    source_path = theta_path
  )
}

load_saved_beta_mixture2_theta <- function(dataset_label) {
  subdir <- if (identical(dataset_label, "short_period")) {
    "01_short_period_beta_mixture2"
  } else {
    "02_long_period_beta_mixture2"
  }
  theta_path <- canonical_comets_mixture_dir(
    "beta_mixture2_short_long_B1000",
    "fast",
    subdir,
    "theta_hat.csv"
  )
  theta_df <- utils::read.csv(theta_path, stringsAsFactors = FALSE)
  list(
    theta = list(
      mu = normalize_projection_vector(theta_df[1L, c("mu_1", "mu_2", "mu_3")], arg_name = "`mu`"),
      weight1 = theta_df$weight1[[1L]],
      alpha1 = theta_df$alpha1[[1L]] %||% theta_df$param1[[1L]],
      beta1 = theta_df$beta1[[1L]] %||% theta_df$param2[[1L]],
      alpha2 = theta_df$alpha2[[1L]] %||% theta_df$param3[[1L]],
      beta2 = theta_df$beta2[[1L]] %||% theta_df$param4[[1L]]
    ),
    source_path = theta_path
  )
}

extract_cardioid_theta_map <- function(stage_bundle,
                                       model_ids = paste0("C", 1:4)) {
  stats::setNames(
    lapply(model_ids, function(model_id) {
      theta_hat <- stage_bundle$results[[model_id]]$observed$theta_hat
      list(
        mu = normalize_projection_vector(theta_hat$mu, arg_name = sprintf("`theta_hat$mu` (%s)", model_id)),
        rho = as.numeric(theta_hat$rho),
        k = as.integer(theta_hat$k),
        p = as.integer(theta_hat$p)
      )
    }),
    model_ids
  )
}

assert_cardioid_theta_maps_equal <- function(lhs,
                                             rhs,
                                             dataset_label,
                                             tol = 1e-10) {
  for (model_id in names(lhs)) {
    lhs_theta <- lhs[[model_id]]
    rhs_theta <- rhs[[model_id]]
    if (lhs_theta$k != rhs_theta$k ||
        lhs_theta$p != rhs_theta$p ||
        max(abs(lhs_theta$mu - rhs_theta$mu)) > tol ||
        abs(lhs_theta$rho - rhs_theta$rho) > tol) {
      stop(sprintf(
        "Stored cardioid parameters differ between statistics for dataset `%s`, model `%s`.",
        dataset_label,
        model_id
      ))
    }
  }
  invisible(TRUE)
}

load_saved_cardioid_theta <- function(dataset_label,
                                      model_id) {
  pipeline_path <- canonical_comets_cardioid_dir(
    "paper_results_B5000_sampleks",
    "fast",
    "pipeline_result.rds"
  )
  pipeline_result <- readRDS(pipeline_path)
  manifest <- pipeline_result$manifest
  stage_pattern <- if (identical(dataset_label, "short_period")) {
    "Short-period cardioid"
  } else {
    "Oort cardioid"
  }
  stage_dirs <- manifest$stage_dir[grepl(stage_pattern, manifest$stage_label, fixed = TRUE)]
  if (length(stage_dirs) != 2L) {
    stop(sprintf("Expected two cardioid stage directories for dataset `%s`.", dataset_label))
  }

  stage_bundles <- lapply(stage_dirs, function(stage_dir) {
    readRDS(file.path(stage_dir, "stage_bundle.rds"))
  })
  theta_maps <- lapply(stage_bundles, extract_cardioid_theta_map)
  assert_cardioid_theta_maps_equal(
    lhs = theta_maps[[1L]],
    rhs = theta_maps[[2L]],
    dataset_label = dataset_label
  )

  list(
    theta = theta_maps[[1L]][[model_id]],
    source_path = paste(file.path(stage_dirs, "stage_bundle.rds"), collapse = ";")
  )
}

cardioid_projected_cdf_on_axis <- function(z,
                                           theta,
                                           gamma) {
  z <- pmin(pmax(as.numeric(z), -1), 1)
  gamma <- normalize_projection_vector(gamma, arg_name = "`gamma`")
  p_proj_car_gamma(
    x = z,
    rho = theta$rho,
    k = theta$k,
    p = theta$p %||% length(theta$mu),
    mu = theta$mu,
    gamma = gamma
  )
}

spherical_cauchy_projected_cdf_on_axis <- function(z,
                                                   theta,
                                                   gamma) {
  spherical_cauchy_projected_cdf(
    x = z,
    omega = normalize_projection_vector(gamma, arg_name = "`gamma`"),
    mu = theta$mu,
    rho = theta$rho,
    warn = FALSE
  )
}

small_circle_projected_cdf_on_axis <- function(z,
                                               theta,
                                               gamma,
                                               method = c("integral", "legendre"),
                                               quad_n = 2048L,
                                               l_max = 200L,
                                               tol = 1e-10) {
  method <- match.arg(method)
  z <- pmin(pmax(as.numeric(z), -1), 1)
  t_values <- acos(z)
  tail_profile <- distance_profile_small_circle(
    omega = normalize_projection_vector(gamma, arg_name = "`gamma`"),
    t_values = t_values,
    mu = theta$mu,
    kappa = theta$kappa,
    nu = theta$nu,
    distance_type = "geodesic",
    method = method,
    quad_n = as.integer(quad_n),
    l_max = as.integer(l_max),
    tol = tol,
    validate_against_integral = FALSE
  )
  monotone_clip_cdf(z_grid = z, cdf_values = 1 - tail_profile)
}

beta_mixture2_projected_cdf_on_axis <- function(z,
                                                theta,
                                                gamma,
                                                method = c("integral", "legendre"),
                                                quad_n = 2048L,
                                                l_max = 200L,
                                                tol = 1e-10) {
  method <- match.arg(method)
  z <- pmin(pmax(as.numeric(z), -1), 1)
  t_values <- acos(z)
  tail_profile <- distance_profile_beta_mixture2(
    omega = normalize_projection_vector(gamma, arg_name = "`gamma`"),
    t_values = t_values,
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    distance_type = "geodesic",
    method = method,
    quad_n = as.integer(quad_n),
    l_max = as.integer(l_max),
    tol = tol,
    validate_against_integral = FALSE
  )
  monotone_clip_cdf(z_grid = z, cdf_values = 1 - tail_profile)
}

uniform_beta_mixture_projected_cdf_on_axis <- function(z,
                                                       theta,
                                                       gamma,
                                                       method = c("integral", "legendre"),
                                                       quad_n = 2048L,
                                                       l_max = 200L,
                                                       tol = 1e-10) {
  method <- match.arg(method)
  z <- pmin(pmax(as.numeric(z), -1), 1)
  t_values <- acos(z)
  tail_profile <- distance_profile_uniform_beta_mixture(
    omega = normalize_projection_vector(gamma, arg_name = "`gamma`"),
    t_values = t_values,
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta,
    distance_type = "geodesic",
    method = method,
    quad_n = as.integer(quad_n),
    l_max = as.integer(l_max),
    tol = tol,
    validate_against_integral = FALSE
  )
  monotone_clip_cdf(z_grid = z, cdf_values = 1 - tail_profile)
}

projected_fit_summary_row <- function(dataset_label,
                                      model_name,
                                      axis_info,
                                      source_path,
                                      fit,
                                      theta_fields) {
  cbind(
    data.frame(
      dataset = dataset_label,
      model = model_name,
      axis_type = axis_info$axis_type,
      sign_rule = axis_info$sign_rule,
      parameter_source = source_path,
      projected_ks = fit$ks,
      projected_cvm = fit$cvm,
      stringsAsFactors = FALSE
    ),
    theta_fields
  )
}

theta_fields_cardioid <- function(theta) {
  data.frame(
    mu_1 = theta$mu[[1L]],
    mu_2 = theta$mu[[2L]],
    mu_3 = theta$mu[[3L]],
    rho = theta$rho,
    k = theta$k,
    stringsAsFactors = FALSE
  )
}

theta_fields_spherical_cauchy <- function(theta) {
  data.frame(
    mu_1 = theta$mu[[1L]],
    mu_2 = theta$mu[[2L]],
    mu_3 = theta$mu[[3L]],
    rho = theta$rho,
    stringsAsFactors = FALSE
  )
}

theta_fields_small_circle <- function(theta) {
  data.frame(
    mu_1 = theta$mu[[1L]],
    mu_2 = theta$mu[[2L]],
    mu_3 = theta$mu[[3L]],
    kappa = theta$kappa,
    nu = theta$nu,
    stringsAsFactors = FALSE
  )
}

theta_fields_beta_mixture2 <- function(theta) {
  data.frame(
    mu_1 = theta$mu[[1L]],
    mu_2 = theta$mu[[2L]],
    mu_3 = theta$mu[[3L]],
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2,
    stringsAsFactors = FALSE
  )
}

theta_fields_uniform_beta_mixture <- function(theta) {
  data.frame(
    mu_1 = theta$mu[[1L]],
    mu_2 = theta$mu[[2L]],
    mu_3 = theta$mu[[3L]],
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta,
    stringsAsFactors = FALSE
  )
}

run_single_comets_projected_axis_diagnostic <- function(data_matrix,
                                                        dataset_label,
                                                        model_name,
                                                        output_dir,
                                                        axis_info,
                                                        grid_size = 1001L,
                                                        rotational_method = c("integral", "legendre"),
                                                        rotational_quad_n = 2048L,
                                                        rotational_l_max = 200L,
                                                        rotational_tol = 1e-10,
                                                        save_plot = TRUE) {
  rotational_method <- match.arg(rotational_method)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  model_bundle <- switch(
    model_name,
    C1 = load_saved_cardioid_theta(dataset_label = dataset_label, model_id = "C1"),
    C2 = load_saved_cardioid_theta(dataset_label = dataset_label, model_id = "C2"),
    C3 = load_saved_cardioid_theta(dataset_label = dataset_label, model_id = "C3"),
    C4 = load_saved_cardioid_theta(dataset_label = dataset_label, model_id = "C4"),
    small_circle = load_saved_small_circle_theta(dataset_label = dataset_label),
    beta_mixture2 = load_saved_beta_mixture2_theta(dataset_label = dataset_label),
    uniform_beta_mixture = load_saved_uniform_beta_theta(dataset_label = dataset_label),
    stop(sprintf("Unsupported model `%s`.", model_name))
  )
  theta <- model_bundle$theta

  cdf_fun <- switch(
    model_name,
    C1 = function(z_grid) cardioid_projected_cdf_on_axis(z_grid, theta = theta, gamma = axis_info$axis),
    C2 = function(z_grid) cardioid_projected_cdf_on_axis(z_grid, theta = theta, gamma = axis_info$axis),
    C3 = function(z_grid) cardioid_projected_cdf_on_axis(z_grid, theta = theta, gamma = axis_info$axis),
    C4 = function(z_grid) cardioid_projected_cdf_on_axis(z_grid, theta = theta, gamma = axis_info$axis),
    small_circle = function(z_grid) small_circle_projected_cdf_on_axis(
      z_grid,
      theta = theta,
      gamma = axis_info$axis,
      method = rotational_method,
      quad_n = rotational_quad_n,
      l_max = rotational_l_max,
      tol = rotational_tol
    ),
    beta_mixture2 = function(z_grid) beta_mixture2_projected_cdf_on_axis(
      z_grid,
      theta = theta,
      gamma = axis_info$axis,
      method = rotational_method,
      quad_n = rotational_quad_n,
      l_max = rotational_l_max,
      tol = rotational_tol
    ),
    uniform_beta_mixture = function(z_grid) uniform_beta_mixture_projected_cdf_on_axis(
      z_grid,
      theta = theta,
      gamma = axis_info$axis,
      method = rotational_method,
      quad_n = rotational_quad_n,
      l_max = rotational_l_max,
      tol = rotational_tol
    )
  )

  theta_fields <- switch(
    model_name,
    C1 = theta_fields_cardioid(theta),
    C2 = theta_fields_cardioid(theta),
    C3 = theta_fields_cardioid(theta),
    C4 = theta_fields_cardioid(theta),
    small_circle = theta_fields_small_circle(theta),
    beta_mixture2 = theta_fields_beta_mixture2(theta),
    uniform_beta_mixture = theta_fields_uniform_beta_mixture(theta)
  )

  z_observed <- project_data_onto_axis(data_matrix, axis_info$axis)
  fit <- build_projected_fit(
    z_observed = z_observed,
    cdf_fun = cdf_fun,
    grid_size = as.integer(grid_size)
  )

  axis_df <- make_axis_metadata_row(dataset_label = dataset_label, axis_info = axis_info)
  summary_df <- projected_fit_summary_row(
    dataset_label = dataset_label,
    model_name = model_name,
    axis_info = axis_info,
    source_path = model_bundle$source_path,
    fit = fit,
    theta_fields = theta_fields
  )
  cdf_grid_df <- data.frame(
    z_grid = fit$z_grid,
    fitted_cdf = fit$fitted_cdf,
    empirical_cdf = fit$empirical_cdf,
    abs_diff = abs(fit$empirical_cdf - fit$fitted_cdf),
    stringsAsFactors = FALSE
  )
  projected_data_df <- data.frame(
    z = z_observed,
    stringsAsFactors = FALSE
  )

  utils::write.csv(summary_df, file = file.path(output_dir, "fit_summary.csv"), row.names = FALSE)
  utils::write.csv(axis_df, file = file.path(output_dir, "projection_direction.csv"), row.names = FALSE)
  utils::write.csv(cdf_grid_df, file = file.path(output_dir, "projected_cdf_grid.csv"), row.names = FALSE)
  utils::write.csv(projected_data_df, file = file.path(output_dir, "projected_data.csv"), row.names = FALSE)

  if (isTRUE(save_plot)) {
    plot_projected_ecdf_overlay(
      z_observed = z_observed,
      fit = fit,
      output_path = file.path(output_dir, "projected_ecdf_overlay.png"),
      dataset_label = dataset_label,
      model_label = model_name,
      axis_type = axis_info$axis_type
    )
  }

  summary_df
}

run_comets_projected_axis_diagnostics <- function(
    output_root = file.path("real_data", "comets", "projected_axis_diagnostics", "director_axes"),
    datasets = c("short", "long"),
    models = c("C1", "C2", "C3", "C4", "small_circle", "beta_mixture2", "uniform_beta_mixture"),
    grid_size = 1001L,
    rotational_method = c("integral", "legendre"),
    rotational_quad_n = 2048L,
    rotational_l_max = 200L,
    rotational_tol = 1e-10,
    save_plot = TRUE) {
  rotational_method <- match.arg(rotational_method)
  datasets <- unique(as.character(datasets))
  models <- unique(as.character(models))

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  comets_data <- load_comets_real_data(finite_normals = "both")
  summary_rows <- list()
  counter <- 0L

  for (dataset_name in datasets) {
    if (!dataset_name %in% c("short", "long")) {
      stop("`datasets` must be a subset of c('short', 'long').")
    }

    data_matrix <- if (identical(dataset_name, "short")) {
      as.matrix(comets_data$short$normal)
    } else {
      as.matrix(comets_data$long$normal)
    }
    dataset_label <- if (identical(dataset_name, "short")) "short_period" else "long_period"
    axis_info <- compute_director_projection_axis(data_matrix = data_matrix, dataset_label = dataset_label)

    for (model_name in models) {
      counter <- counter + 1L
      message(sprintf("[Comets projected axis] %s / %s", dataset_label, model_name))
      output_dir <- file.path(output_root, sprintf("%02d_%s_%s", counter, dataset_label, model_name))
      summary_rows[[length(summary_rows) + 1L]] <- run_single_comets_projected_axis_diagnostic(
        data_matrix = data_matrix,
        dataset_label = dataset_label,
        model_name = model_name,
        output_dir = output_dir,
        axis_info = axis_info,
        grid_size = as.integer(grid_size),
        rotational_method = rotational_method,
        rotational_quad_n = as.integer(rotational_quad_n),
        rotational_l_max = as.integer(rotational_l_max),
        rotational_tol = rotational_tol,
        save_plot = isTRUE(save_plot)
      )
    }
  }

  summary_df <- rbind_fill_projected_summary(summary_rows)
  utils::write.csv(summary_df, file = file.path(output_root, "projected_axis_diagnostics_summary.csv"), row.names = FALSE)
  invisible(summary_df)
}

if (sys.nframe() == 0L) {
  run_comets_projected_axis_diagnostics()
}
