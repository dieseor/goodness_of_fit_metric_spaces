#!/usr/bin/env Rscript

resolve_skyelavas_sensitivity_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

utils_path <- resolve_skyelavas_sensitivity_path(
  "real_data",
  "logistic_gaussian",
  "utils_logistic_gaussian_screening.R"
)
bootstrap_path <- resolve_skyelavas_sensitivity_path(
  "bootstrap",
  "multiplier_bootstrap.R"
)

source(utils_path)
source(bootstrap_path)

if (!requireNamespace("goftest", quietly = TRUE)) {
  stop("The package `goftest` is required for the Aitchison-style sensitivity analysis.")
}

parse_screening_args <- function(args) {
  output <- list()

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pieces[[1]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[key]] <- value
  }

  output
}

parse_integer_arg <- function(x, default) {
  if (is.null(x)) {
    return(as.integer(default))
  }
  as.integer(x)
}

parse_numeric_arg <- function(x, default) {
  if (is.null(x)) {
    return(as.numeric(default))
  }
  as.numeric(x)
}

parse_optional_numeric_arg <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  as.numeric(x)
}

collapse_numeric_vector <- function(x, digits = 10) {
  paste(formatC(as.numeric(x), digits = digits, format = "fg", flag = "#"), collapse = ";")
}

collapse_character_vector <- function(x) {
  paste(as.character(x), collapse = ";")
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  min(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(x)
}

choose_omega_grid_for_subset <- function(x_closed,
                                         omega_grid_type = c("sample_points", "fixed_simplex_lattice"),
                                         max_centers = 100L,
                                         seed = 123L,
                                         boundary_epsilon = NULL) {
  omega_grid_type <- match.arg(omega_grid_type)
  if (identical(omega_grid_type, "fixed_simplex_lattice")) {
    return(build_fixed_simplex_omega_grid(
      ambient_dim = ncol(x_closed),
      max_centers = max_centers,
      boundary_epsilon = boundary_epsilon
    ))
  }

  omega <- as.matrix(x_closed)
  rownames(omega) <- rownames(x_closed)
  list(
    omega = omega,
    lattice_level = NA_integer_,
    n_centers = nrow(omega),
    construction = "sample_points_exact"
  )
}

make_top_ks_discrepancy_rows <- function(excluded_sample_id,
                                         omega_grid_type,
                                         centers,
                                         t_grid,
                                         empirical_profile,
                                         theoretical_profile,
                                         distance_observed,
                                         sample_ids,
                                         top_k = 50L,
                                         scaled_process_matrix = NULL,
                                         localization_source = "external_profiles") {
  n_obs <- length(sample_ids)
  n_centers <- nrow(centers$omega)
  n_thresholds <- length(t_grid)

  orient_profile <- function(profile, name) {
    profile <- as.matrix(profile)
    if (identical(dim(profile), c(n_centers, n_thresholds))) {
      return(profile)
    }
    if (identical(dim(profile), c(n_thresholds, n_centers))) {
      return(t(profile))
    }
    stop(sprintf(
      "%s has incompatible dimensions %s; expected centers x thresholds = %d x %d or thresholds x centers = %d x %d.",
      name,
      paste(dim(profile), collapse = " x "),
      n_centers,
      n_thresholds,
      n_thresholds,
      n_centers
    ))
  }

  empirical_profile_ct <- orient_profile(empirical_profile, "empirical_profile")
  theoretical_profile_ct <- orient_profile(theoretical_profile, "theoretical_profile")

  if (!is.null(scaled_process_matrix)) {
    scaled_signed_matrix <- orient_profile(scaled_process_matrix, "scaled_process_matrix")
    signed_matrix <- scaled_signed_matrix / sqrt(n_obs)
  } else {
    signed_matrix <- empirical_profile_ct - theoretical_profile_ct
    scaled_signed_matrix <- sqrt(n_obs) * signed_matrix
  }
  abs_matrix <- abs(signed_matrix)
  scaled_abs_matrix <- abs(scaled_signed_matrix)

  ordering <- order(scaled_abs_matrix, decreasing = TRUE)
  keep <- ordering[seq_len(min(as.integer(top_k), length(ordering)))]
  locations <- arrayInd(keep, .dim = dim(scaled_abs_matrix))

  rows <- vector("list", nrow(locations))
  for (i in seq_len(nrow(locations))) {
    center_idx <- locations[i, 1]
    threshold_idx <- locations[i, 2]
    t_value <- t_grid[[threshold_idx]]
    center_distances <- distance_observed[, center_idx]
    nearest_center_idx <- which.min(center_distances)
    nearest_center_sample_id <- sample_ids[[nearest_center_idx]]
    distance_center_to_nearest_sample <- center_distances[[nearest_center_idx]]

    threshold_gap <- abs(center_distances - t_value)
    threshold_min <- min(threshold_gap)
    closest_threshold_indices <- which(threshold_gap == threshold_min)
    closest_threshold_ids <- sample_ids[closest_threshold_indices]
    closest_threshold_distances <- center_distances[closest_threshold_indices]

    rows[[i]] <- data.frame(
      excluded_sample_id = excluded_sample_id,
      omega_grid_type = omega_grid_type,
      rank = i,
      center_index = center_idx,
      threshold_index = threshold_idx,
      t = t_value,
      empirical_profile = empirical_profile_ct[center_idx, threshold_idx],
      fitted_profile = theoretical_profile_ct[center_idx, threshold_idx],
      signed_discrepancy = signed_matrix[center_idx, threshold_idx],
      absolute_discrepancy = abs_matrix[center_idx, threshold_idx],
      scaled_signed_discrepancy = scaled_signed_matrix[center_idx, threshold_idx],
      scaled_absolute_discrepancy = scaled_abs_matrix[center_idx, threshold_idx],
      center_coordinates = collapse_numeric_vector(centers$omega[center_idx, ]),
      nearest_sample_id_to_center = nearest_center_sample_id,
      distance_center_to_nearest_sample = distance_center_to_nearest_sample,
      sample_ids_closest_to_threshold = collapse_character_vector(closest_threshold_ids),
      distances_closest_to_threshold = collapse_numeric_vector(closest_threshold_distances),
      distance_to_threshold_min = threshold_min,
      profile_orientation = "centers_x_thresholds",
      localization_source = localization_source,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

run_single_distance_profile_analysis <- function(x_comp,
                                                 sample_ids,
                                                 excluded_sample_id,
                                                 omega_grid_type,
                                                 B,
                                                 n_cores,
                                                 seed,
                                                 max_centers,
                                                 n_t,
                                                 t_grid_tail_prob,
                                                 boundary_epsilon,
                                                 quadform_method) {
  started_proc <- proc.time()[["elapsed"]]

  x_closed <- close_composition_rows(x_comp)
  fit <- fit_logistic_gaussian_plugin(x_closed = x_closed, ridge = 1e-8)
  centers <- choose_omega_grid_for_subset(
    x_closed = x_closed,
    omega_grid_type = omega_grid_type,
    max_centers = max_centers,
    seed = seed,
    boundary_epsilon = boundary_epsilon
  )
  z_centers <- ilr_transform_closed(centers$omega, V = fit$ilr_basis)
  distance_observed <- distance_matrix_from_ilr(fit$Z, z_centers)
  t_grid_info <- build_fixed_t_grid_logistic_gaussian(
    fit = fit,
    omega_grid = centers$omega,
    n_t = n_t,
    tail_prob = t_grid_tail_prob
  )
  t_grid <- t_grid_info$t_grid

  control <- normalize_logistic_gaussian_screening_control(
    list(mvnormal_quadform_method = quadform_method)
  )

  ## For localization of the observed KS maximum, use the same fitted-profile
  ## route as the standard screening code before the composite bootstrap call.
  ## In particular, do not force the HBE approximation here: for extreme
  ## artificial lattice centers it may return degenerate tail values and then
  ## the stored top discrepancies no longer correspond to the reported KS
  ## statistic.
  theoretical_profile <- evaluate_fitted_logistic_gaussian_profile(
    omega_grid = centers$omega,
    t_grid = t_grid,
    fit = fit
  )
  empirical_profile <- compute_profile_matrix_from_distances(distance_observed, t_grid)
  observed_stats <- compute_screening_statistics_from_profiles(
    empirical_profile = empirical_profile,
    theoretical_profile = theoretical_profile,
    n_obs = nrow(x_closed)
  )

  external_ks_statistic <- find_named_numeric_value(
    observed_stats,
    c("ks", "ks_statistic", "KS", "Kolmogorov", "kolmogorov_smirnov")
  )
  external_cvm_statistic <- find_named_numeric_value(
    observed_stats,
    c("cvm", "cvm_statistic", "CvM", "cramer_von_mises")
  )

  bootstrap_result <- multiplier_bootstrap_logistic_gaussian(
    data = x_closed,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = n_cores,
    seed = seed,
    control = control,
    unknown_param = "both"
  )

  ks_extracted <- extract_metric_result(bootstrap_result, "ks")
  cvm_extracted <- extract_metric_result(bootstrap_result, "cvm")
  diagnostics_mardia <- mardia_multivariate_normality(fit$Z)
  shapiro_min_p <- safe_min(compute_marginal_shapiro_pvalues(fit$Z))

  localization_empirical_profile <- empirical_profile
  localization_theoretical_profile <- theoretical_profile
  localization_scaled_process_matrix <- NULL
  localization_source <- "external_profiles"

  if (!is.null(bootstrap_result$observed$ks)) {
    observed_ks <- bootstrap_result$observed$ks
    if (!is.null(observed_ks$empirical_profile)) {
      localization_empirical_profile <- observed_ks$empirical_profile
      localization_source <- "bootstrap_result$observed$ks$profiles"
    }
    if (!is.null(observed_ks$theoretical_profile)) {
      localization_theoretical_profile <- observed_ks$theoretical_profile
      localization_source <- "bootstrap_result$observed$ks$profiles"
    }
    if (!is.null(observed_ks$process_matrix)) {
      localization_scaled_process_matrix <- observed_ks$process_matrix
      localization_source <- "bootstrap_result$observed$ks$process_matrix"
    } else if (!is.null(observed_ks$process)) {
      localization_scaled_process_matrix <- observed_ks$process
      localization_source <- "bootstrap_result$observed$ks$process"
    }
  }

  top_discrepancies <- make_top_ks_discrepancy_rows(
    excluded_sample_id = excluded_sample_id,
    omega_grid_type = omega_grid_type,
    centers = centers,
    t_grid = t_grid,
    empirical_profile = localization_empirical_profile,
    theoretical_profile = localization_theoretical_profile,
    distance_observed = distance_observed,
    sample_ids = sample_ids,
    top_k = 50L,
    scaled_process_matrix = localization_scaled_process_matrix,
    localization_source = localization_source
  )

  ks_top_scaled_max <- top_discrepancies$scaled_absolute_discrepancy[[1L]]
  ks_summary_abs_diff <- abs(ks_top_scaled_max - ks_extracted$statistic)
  ks_external_abs_diff <- abs(ks_top_scaled_max - external_ks_statistic)
  ks_top_matches_summary <- is.finite(ks_summary_abs_diff) && ks_summary_abs_diff <= 1e-8
  ks_top_matches_external <- is.finite(ks_external_abs_diff) && ks_external_abs_diff <= 1e-8
  if (!ks_top_matches_summary || !ks_top_matches_external) {
    warning(sprintf(
      "KS localization check for excluded_sample_id=%s, omega_grid_type=%s: top scaled discrepancy %.12g, extracted KS %.12g, external-profile KS %.12g.",
      excluded_sample_id,
      omega_grid_type,
      ks_top_scaled_max,
      ks_extracted$statistic,
      external_ks_statistic
    ))
  }

  elapsed_seconds <- unname(proc.time()[["elapsed"]] - started_proc)
  summary_row <- data.frame(
    excluded_sample_id = excluded_sample_id,
    omega_grid_type = omega_grid_type,
    n = nrow(x_closed),
    D = ncol(x_closed),
    n_centers = centers$n_centers,
    omega_grid_construction = centers$construction,
    boundary_epsilon = if (identical(omega_grid_type, "fixed_simplex_lattice")) centers$boundary_epsilon %||% boundary_epsilon else NA_real_,
    n_t = length(t_grid),
    t_grid_tail_prob = t_grid_tail_prob,
    t_grid_max = max(t_grid),
    ks_statistic = ks_extracted$statistic,
    ks_top_scaled_max = ks_top_scaled_max,
    ks_top_summary_abs_diff = ks_summary_abs_diff,
    ks_top_matches_summary = ks_top_matches_summary,
    external_ks_statistic = external_ks_statistic,
    external_cvm_statistic = external_cvm_statistic,
    ks_top_external_abs_diff = ks_external_abs_diff,
    ks_top_matches_external = ks_top_matches_external,
    ks_statistic_node_source = ks_extracted$statistic_node_source %||% NA_character_,
    ks_statistic_node_names = collapse_character_vector(ks_extracted$statistic_node_names %||% character(0)),
    cvm_statistic_node_source = cvm_extracted$statistic_node_source %||% NA_character_,
    cvm_statistic_node_names = collapse_character_vector(cvm_extracted$statistic_node_names %||% character(0)),
    ks_pvalue = ks_extracted$p_value,
    cvm_statistic = cvm_extracted$statistic,
    cvm_pvalue = cvm_extracted$p_value,
    mardia_skew_pvalue = diagnostics_mardia$skewness_p_value,
    mardia_kurtosis_pvalue = diagnostics_mardia$kurtosis_p_value,
    shapiro_min_pvalue = shapiro_min_p,
    eigenvalues = collapse_numeric_vector(fit$eigenvalues),
    condition_number = fit$condition_number,
    elapsed_seconds = elapsed_seconds,
    localization_source = localization_source,
    stringsAsFactors = FALSE
  )

  list(
    summary_row = summary_row,
    top_discrepancies = top_discrepancies
  )
}

alr_transform_last <- function(x_closed) {
  x_closed <- as.matrix(x_closed)
  if (ncol(x_closed) < 2L) {
    stop("ALR transform requires at least two compositional parts.")
  }
  log(x_closed[, -ncol(x_closed), drop = FALSE] / x_closed[, ncol(x_closed)])
}

whiten_alr_data <- function(y, ridge = 1e-8, tol = 1e-10) {
  y <- as.matrix(y)
  storage.mode(y) <- "double"
  mu_hat <- colMeans(y)
  centered <- sweep(y, 2L, mu_hat, FUN = "-")
  sigma_hat <- stats::cov(y)
  sigma_hat <- 0.5 * (sigma_hat + t(sigma_hat))
  eig <- eigen(sigma_hat, symmetric = TRUE)
  min_eig <- min(eig$values)
  ridge_added <- 0
  if (min_eig <= tol) {
    ridge_added <- max(ridge, tol - min_eig + ridge)
    sigma_hat <- sigma_hat + ridge_added * diag(ncol(sigma_hat))
    eig <- eigen(sigma_hat, symmetric = TRUE)
  }
  inv_sqrt <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, tol)), nrow = ncol(y)) %*% t(eig$vectors)
  z <- centered %*% inv_sqrt
  list(
    z = z,
    mu_hat = mu_hat,
    Sigma_hat = sigma_hat,
    eigenvalues = eig$values,
    ridge_added = ridge_added,
    condition_number = max(eig$values) / min(eig$values)
  )
}

run_aitchison_style_tests <- function(x_comp, excluded_sample_id) {
  x_closed <- close_composition_rows(x_comp)
  y <- alr_transform_last(x_closed)
  whitened <- whiten_alr_data(y)
  z <- as.matrix(whitened$z)
  q_dim <- ncol(z)

  test_rows <- list()
  append_test <- function(variable, family, x, null_cdf, null_name) {
    ks_p <- stats::ks.test(x, null_cdf)$p.value
    cvm_p <- goftest::cvm.test(x, null = null_cdf, nullname = null_name)$p.value
    test_rows[[length(test_rows) + 1L]] <<- data.frame(
      excluded_sample_id = excluded_sample_id,
      family = family,
      variable = variable,
      null_distribution = null_name,
      ks_pvalue = ks_p,
      cvm_pvalue = cvm_p,
      aitchison_style_approx = TRUE,
      stephens_correction = FALSE,
      stringsAsFactors = FALSE
    )
  }

  for (j in seq_len(q_dim)) {
    append_test(
      variable = sprintf("marginal_z%d", j),
      family = "marginal",
      x = z[, j],
      null_cdf = "pnorm",
      null_name = "N(0,1)"
    )
  }

  for (j in seq_len(q_dim - 1L)) {
    for (k in seq.int(j + 1L, q_dim)) {
      u <- (atan2(z[, k], z[, j]) + pi) / (2 * pi)
      append_test(
        variable = sprintf("angle_z%d_z%d", j, k),
        family = "angle",
        x = u,
        null_cdf = "punif",
        null_name = "U(0,1)"
      )
    }
  }

  radius_sq <- rowSums(z^2)
  append_test(
    variable = sprintf("radius_sq_chisq_%d", q_dim),
    family = "radius",
    x = radius_sq,
    null_cdf = function(q) stats::pchisq(q, df = q_dim),
    null_name = sprintf("Chi-square(%d)", q_dim)
  )

  tests_df <- do.call(rbind, test_rows)
  min_ks_idx <- which.min(tests_df$ks_pvalue)
  min_cvm_idx <- which.min(tests_df$cvm_pvalue)

  summary_row <- data.frame(
    excluded_sample_id = excluded_sample_id,
    n = nrow(x_closed),
    D = ncol(x_closed),
    alr_dim = q_dim,
    n_univariate_quantities = nrow(tests_df),
    n_tests_total = 2L * nrow(tests_df),
    aitchison_style_approx = TRUE,
    stephens_correction = FALSE,
    ks_rejections_5pct = sum(tests_df$ks_pvalue <= 0.05),
    cvm_rejections_5pct = sum(tests_df$cvm_pvalue <= 0.05),
    ks_rejections_marginal_5pct = sum(tests_df$family == "marginal" & tests_df$ks_pvalue <= 0.05),
    ks_rejections_angle_5pct = sum(tests_df$family == "angle" & tests_df$ks_pvalue <= 0.05),
    ks_rejections_radius_5pct = sum(tests_df$family == "radius" & tests_df$ks_pvalue <= 0.05),
    cvm_rejections_marginal_5pct = sum(tests_df$family == "marginal" & tests_df$cvm_pvalue <= 0.05),
    cvm_rejections_angle_5pct = sum(tests_df$family == "angle" & tests_df$cvm_pvalue <= 0.05),
    cvm_rejections_radius_5pct = sum(tests_df$family == "radius" & tests_df$cvm_pvalue <= 0.05),
    ks_min_pvalue = tests_df$ks_pvalue[[min_ks_idx]],
    ks_min_family = tests_df$family[[min_ks_idx]],
    ks_min_variable = tests_df$variable[[min_ks_idx]],
    cvm_min_pvalue = tests_df$cvm_pvalue[[min_cvm_idx]],
    cvm_min_family = tests_df$family[[min_cvm_idx]],
    cvm_min_variable = tests_df$variable[[min_cvm_idx]],
    eigenvalues = collapse_numeric_vector(whitened$eigenvalues),
    condition_number = whitened$condition_number,
    ridge_added = whitened$ridge_added,
    stringsAsFactors = FALSE
  )

  list(
    summary_row = summary_row,
    individual_tests = tests_df
  )
}

run_skyelavas_sensitivity_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- parse_screening_args(args)
  B <- parse_integer_arg(args$B, 1000L)
  n_cores <- parse_integer_arg(args$n_cores, 12L)
  seed <- parse_integer_arg(args$seed, 123L)
  max_centers <- parse_integer_arg(args$max_centers, 100L)
  n_t <- parse_integer_arg(args$n_t, 60L)
  t_grid_tail_prob <- parse_numeric_arg(args$t_grid_tail_prob, 1e-8)
  boundary_epsilon <- parse_optional_numeric_arg(args$boundary_epsilon)
  # EN DUDA (2026-07-26): shared exact dispatcher is the production default.
  quadform_method <- as.character(args$quadform_method %||% "auto")
  excluded_sample_ids_arg <- as.character(args$excluded_sample_ids %||% "982")
  excluded_sample_ids_arg <- trimws(strsplit(excluded_sample_ids_arg, ",", fixed = TRUE)[[1]])
  run_all_exclusions <- length(excluded_sample_ids_arg) == 1L && identical(tolower(excluded_sample_ids_arg), "all")
  out_dir <- as.character(args$out_dir %||% file.path(
    "output", "calibration", "bootstrap", "logistic_gaussian", "skyelavas_sensitivity"
  ))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  base_data <- prepare_composition_dataset("SkyeLavas")
  sample_ids_full <- as.character(base_data$X_raw$sample_id)
  x_comp_full <- as.matrix(base_data$X_comp)
  rownames(x_comp_full) <- sample_ids_full

  cat("\nSkyeLavas sensitivity analysis\n")
  cat("Full basalt sample size:", length(sample_ids_full), "\n")
  cat("Full basalt sample ids:", paste(sample_ids_full, collapse = ", "), "\n")
  if ("982" %in% sample_ids_full) {
    cat("Excluding sample 982 gives the current SkyeLavasAitchison32 reconstruction.\n")
  } else {
    warning("Sample 982 was not found in SkyeLavas; check the local data constructor.")
  }
  if (run_all_exclusions) {
    excluded_sample_ids <- sample_ids_full
    cat("Exclusions to run: all 33 leave-one-out exclusions.\n")
  } else {
    excluded_sample_ids <- excluded_sample_ids_arg
    missing_exclusions <- setdiff(excluded_sample_ids, sample_ids_full)
    if (length(missing_exclusions) > 0L) {
      stop(sprintf(
        "The following requested excluded_sample_ids are not in SkyeLavas: %s",
        paste(missing_exclusions, collapse = ", ")
      ))
    }
    cat("Exclusions to run:", paste(excluded_sample_ids, collapse = ", "), "\n")
  }
  cat("B:", B, "\n")
  cat("n_cores:", n_cores, "\n")
  cat("seed:", seed, "\n")
  cat("max_centers:", max_centers, "\n")
  cat("n_t:", n_t, "\n")
  cat("t_grid_tail_prob:", format(t_grid_tail_prob, scientific = TRUE), "\n")
  cat("boundary_epsilon:", if (is.null(boundary_epsilon)) "default(0.5/D)" else format(boundary_epsilon, scientific = TRUE), "\n")
  cat("quadform_method:", quadform_method, "\n")
  cat("out_dir:", out_dir, "\n\n")

  omega_grid_types <- c("fixed_simplex_lattice", "sample_points")
  dp_rows <- list()
  aitchison_rows <- list()
  aitchison_individual_rows <- list()
  discrepancy_rows <- list()

  for (i in seq_along(excluded_sample_ids)) {
    excluded_sample_id <- excluded_sample_ids[[i]]
    keep <- sample_ids_full != excluded_sample_id
    x_comp_subset <- x_comp_full[keep, , drop = FALSE]
    sample_ids_subset <- sample_ids_full[keep]

    cat(sprintf("[%s] leave-one-out excluding %s\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                excluded_sample_id))

    aitchison_analysis <- run_aitchison_style_tests(
      x_comp = x_comp_subset,
      excluded_sample_id = excluded_sample_id
    )
    aitchison_rows[[i]] <- aitchison_analysis$summary_row
    aitchison_individual_rows[[i]] <- aitchison_analysis$individual_tests

    for (grid_type in omega_grid_types) {
      analysis_seed <- as.integer(seed) + 1000L * (i - 1L) + if (identical(grid_type, "sample_points")) 500L else 0L
      analysis <- run_single_distance_profile_analysis(
        x_comp = x_comp_subset,
        sample_ids = sample_ids_subset,
        excluded_sample_id = excluded_sample_id,
        omega_grid_type = grid_type,
        B = B,
        n_cores = n_cores,
        seed = analysis_seed,
        max_centers = max_centers,
        n_t = n_t,
        t_grid_tail_prob = t_grid_tail_prob,
        boundary_epsilon = boundary_epsilon,
        quadform_method = quadform_method
      )
      dp_rows[[length(dp_rows) + 1L]] <- analysis$summary_row
      discrepancy_rows[[length(discrepancy_rows) + 1L]] <- analysis$top_discrepancies
    }
  }

  dp_df <- do.call(rbind, dp_rows)
  aitchison_df <- do.call(rbind, aitchison_rows)
  aitchison_individual_df <- do.call(rbind, aitchison_individual_rows)
  discrepancies_df <- do.call(rbind, discrepancy_rows)

  dp_path <- file.path(out_dir, "dp_leave_one_out.csv")
  aitchison_path <- file.path(out_dir, "aitchison_style_leave_one_out.csv")
  aitchison_individual_path <- file.path(out_dir, "aitchison_style_individual_tests.csv")
  discrepancies_path <- file.path(out_dir, "ks_top_discrepancies.csv")

  utils::write.csv(dp_df, file = dp_path, row.names = FALSE)
  utils::write.csv(aitchison_df, file = aitchison_path, row.names = FALSE)
  utils::write.csv(aitchison_individual_df, file = aitchison_individual_path, row.names = FALSE)
  utils::write.csv(discrepancies_df, file = discrepancies_path, row.names = FALSE)

  cat("dp_leave_one_out.csv:", dp_path, "\n")
  cat("aitchison_style_leave_one_out.csv:", aitchison_path, "\n")
  cat("aitchison_style_individual_tests.csv:", aitchison_individual_path, "\n")
  cat("ks_top_discrepancies.csv:", discrepancies_path, "\n")
}

if (sys.nframe() == 0L) {
  run_skyelavas_sensitivity_cli()
}
