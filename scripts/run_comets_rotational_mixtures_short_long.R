resolve_rotmix_comets_path <- function(...) {
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

model_specs_script_path_rotmix_comets <- resolve_rotmix_comets_path("bootstrap", "model_specs.R")
multiplier_bootstrap_script_path_rotmix_comets <- resolve_rotmix_comets_path("bootstrap", "multiplier_bootstrap.R")
utils_script_path_rotmix_comets <- resolve_rotmix_comets_path("utils.R")
comets_utils_script_path_rotmix_comets <- resolve_rotmix_comets_path(
  "real_data",
  "comets",
  "utils_comets_data.R"
)

source(utils_script_path_rotmix_comets)
source(model_specs_script_path_rotmix_comets)
source(multiplier_bootstrap_script_path_rotmix_comets)
source(comets_utils_script_path_rotmix_comets)

make_beta_mixture2_spec <- get("make_beta_mixture2_spec", mode = "function")
make_logitnormal_mixture2_spec <- get("make_logitnormal_mixture2_spec", mode = "function")
generate_canonical_lattice <- get("generate_canonical_lattice", mode = "function")
multiplier_bootstrap_gof <- get("multiplier_bootstrap_gof", mode = "function")

timestamp_tag_rotmix_comets <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

write_lines_if_possible_rotmix_comets <- function(lines, path) {
  writeLines(as.character(lines), con = path)
  invisible(path)
}

load_comets_distance_profile_data_rotmix <- function() {
  load_comets_real_data(finite_normals = "both")
}

rotmix_model_n_parameters <- function(model_name) {
  if (!model_name %in% c("beta_mixture2", "logitnormal_mixture2")) {
    stop(sprintf("Unsupported model name: %s", model_name))
  }
  6L
}

rotmix_model_density_log <- function(model_name, data_matrix, theta) {
  if (identical(model_name, "beta_mixture2")) {
    return(d_sph_beta_mixture2_s2(
      x = data_matrix,
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      log = TRUE
    ))
  }

  if (identical(model_name, "logitnormal_mixture2")) {
    return(d_sph_logitnormal_mixture2_s2(
      x = data_matrix,
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      log = TRUE
    ))
  }

  stop(sprintf("Unsupported model name: %s", model_name))
}

rotmix_projected_density_z <- function(model_name, z, theta) {
  if (identical(model_name, "beta_mixture2")) {
    return(beta_mixture2_density_gz(
      z = z,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    ))
  }

  if (identical(model_name, "logitnormal_mixture2")) {
    return(logitnormal_mixture2_density_gz(
      z = z,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2
    ))
  }

  stop(sprintf("Unsupported model name: %s", model_name))
}

rotmix_projected_cdf_z <- function(model_name, z, theta) {
  y <- (as.numeric(z) + 1) / 2

  if (identical(model_name, "beta_mixture2")) {
    return(beta_mixture2_cdf_y(
      y = y,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    ))
  }

  if (identical(model_name, "logitnormal_mixture2")) {
    return(logitnormal_mixture2_cdf_y(
      y = y,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2
    ))
  }

  stop(sprintf("Unsupported model name: %s", model_name))
}

rotmix_fit_model_theta <- function(model_name, data_matrix, control = list()) {
  if (identical(model_name, "beta_mixture2")) {
    return(fit_beta_mixture2_theta(
      data = data_matrix,
      null = list(type = "composite"),
      control = control
    ))
  }

  if (identical(model_name, "logitnormal_mixture2")) {
    return(fit_logitnormal_mixture2_theta(
      data = data_matrix,
      null = list(type = "composite"),
      control = control
    ))
  }

  stop(sprintf("Unsupported model name: %s", model_name))
}

rotmix_make_spec <- function(model_name,
                             distance_type = "geodesic") {
  if (identical(model_name, "beta_mixture2")) {
    return(make_beta_mixture2_spec(distance_type = distance_type))
  }

  if (identical(model_name, "logitnormal_mixture2")) {
    return(make_logitnormal_mixture2_spec(distance_type = distance_type))
  }

  stop(sprintf("Unsupported model name: %s", model_name))
}

rotmix_projected_gof <- function(z,
                                 cdf_fun,
                                 grid_size = 1001L) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  if (length(z) == 0L) {
    stop("`z` must contain at least one finite projected value.")
  }

  ecdf_z <- stats::ecdf(z)
  z_grid <- sort(unique(c(seq(-1, 1, length.out = as.integer(grid_size)), z)))
  fitted_cdf <- cdf_fun(z_grid)
  empirical_cdf <- ecdf_z(z_grid)

  list(
    ks = max(abs(empirical_cdf - fitted_cdf)),
    cvm = mean((empirical_cdf - fitted_cdf)^2),
    z_grid = z_grid,
    empirical_cdf = empirical_cdf,
    fitted_cdf = fitted_cdf
  )
}

plot_rotmix_projected_diagnostics <- function(z,
                                              density_grid_df,
                                              cdf_grid_df,
                                              output_dir,
                                              dataset_label,
                                              model_label) {
  hist_path <- file.path(output_dir, "projected_hist_density.png")
  cdf_path <- file.path(output_dir, "projected_ecdf_cdf.png")

  grDevices::png(hist_path, width = 1200, height = 900, res = 140)
  hist(
    z,
    breaks = "FD",
    freq = FALSE,
    col = "#d9e6f2",
    border = "white",
    main = sprintf("%s: projected density (%s)", dataset_label, model_label),
    xlab = "z = <mu_hat, x>"
  )
  lines(density_grid_df$z_grid, density_grid_df$density, lwd = 2, col = "#c23b23")
  grid(col = "#d9d9d9")
  grDevices::dev.off()

  grDevices::png(cdf_path, width = 1200, height = 900, res = 140)
  plot(
    NA,
    NA,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    xlab = "z = <mu_hat, x>",
    ylab = "CDF",
    main = sprintf("%s: projected CDF (%s)", dataset_label, model_label)
  )
  z_emp <- sort(unique(c(-1, z, 1)))
  lines(z_emp, stats::ecdf(z)(z_emp), type = "s", lwd = 2, col = "black")
  lines(cdf_grid_df$z_grid, cdf_grid_df$cdf, lwd = 2, col = "#1f78b4")
  legend(
    "topleft",
    legend = c("ECDF", "Fitted CDF"),
    col = c("black", "#1f78b4"),
    lwd = c(2, 2),
    lty = c(1, 1),
    bty = "n"
  )
  grid(col = "#d9d9d9")
  grDevices::dev.off()

  list(hist_path = hist_path, cdf_path = cdf_path)
}

rotmix_theta_summary_row <- function(model_name, theta) {
  if (identical(model_name, "beta_mixture2")) {
    return(data.frame(
      model = model_name,
      mu_1 = theta$mu[[1L]],
      mu_2 = theta$mu[[2L]],
      mu_3 = theta$mu[[3L]],
      weight1 = theta$weight1,
      param1 = theta$alpha1,
      param2 = theta$beta1,
      param3 = theta$alpha2,
      param4 = theta$beta2,
      min_shape_hat = min(theta$alpha1, theta$beta1, theta$alpha2, theta$beta2),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    model = model_name,
    mu_1 = theta$mu[[1L]],
    mu_2 = theta$mu[[2L]],
    mu_3 = theta$mu[[3L]],
    weight1 = theta$weight1,
    param1 = theta$mean1,
    param2 = theta$sd1,
    param3 = theta$mean2,
    param4 = theta$sd2,
    stringsAsFactors = FALSE
  )
}

run_single_rotmix_comet_fit <- function(data_matrix,
                                        dataset_label,
                                        model_name,
                                        output_dir,
                                        B = 500L,
                                        statistics = c("ks", "cvm"),
                                        n_cores = 6L,
                                        seed = 20260601L,
                                        M_value = 60L,
                                        ks_t_points = 200L,
                                        bootstrap_method = "reestimated",
                                        distance_type = "geodesic",
                                        control = list()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  theta_hat <- rotmix_fit_model_theta(
    model_name = model_name,
    data_matrix = data_matrix,
    control = control
  )

  z <- pmin(pmax(as.numeric(data_matrix %*% theta_hat$mu), -1), 1)
  loglik <- sum(rotmix_model_density_log(model_name, data_matrix, theta_hat))
  n <- nrow(data_matrix)
  k <- rotmix_model_n_parameters(model_name)
  aic <- 2 * k - 2 * loglik
  bic <- log(n) * k - 2 * loglik

  projected_fit <- rotmix_projected_gof(
    z = z,
    cdf_fun = function(z_grid) rotmix_projected_cdf_z(model_name, z_grid, theta_hat)
  )

  z_grid <- seq(-0.999, 0.999, length.out = 1001L)
  density_grid_df <- data.frame(
    z_grid = z_grid,
    density = rotmix_projected_density_z(model_name, z_grid, theta_hat),
    stringsAsFactors = FALSE
  )
  cdf_grid_df <- data.frame(
    z_grid = projected_fit$z_grid,
    cdf = projected_fit$fitted_cdf,
    empirical_cdf = projected_fit$empirical_cdf,
    stringsAsFactors = FALSE
  )

  spec <- rotmix_make_spec(model_name, distance_type = distance_type)
  ks_grid <- list(
    omega_grid = generate_canonical_lattice(as.integer(M_value), dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = as.integer(ks_t_points))
  )

  run_one_statistic <- function(statistic_name, statistic_seed) {
    gc(verbose = FALSE)
    multiplier_bootstrap_gof(
      data = data_matrix,
      spec = spec,
      null = list(type = "composite"),
      statistics = statistic_name,
      ks_grid = if (identical(statistic_name, "ks")) ks_grid else NULL,
      B = as.integer(B),
      alpha = 0.05,
      n_cores = as.integer(n_cores),
      seed = as.integer(statistic_seed),
      bootstrap_method = bootstrap_method,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = FALSE,
        bootstrap_thetas = FALSE
      ),
      control = control,
      observed_theta_hat = theta_hat
    )
  }

  statistics <- intersect(as.character(statistics), c("ks", "cvm"))
  if (length(statistics) == 0L) {
    stop("`statistics` must contain at least one of 'ks' or 'cvm'.")
  }

  gof_parts <- lapply(seq_along(statistics), function(j) {
    run_one_statistic(statistics[[j]], as.integer(seed) + 1000L * j)
  })
  names(gof_parts) <- statistics

  extract_stat_inference <- function(part, statistic_name) {
    inference <- part$inference
    if (is.list(inference) && !is.null(inference[[statistic_name]])) {
      return(inference[[statistic_name]])
    }
    inference
  }

  gof_result <- list(
    inference = stats::setNames(
      lapply(seq_along(statistics), function(j) {
        extract_stat_inference(gof_parts[[j]], statistics[[j]])
      }),
      statistics
    ),
    diagnostics = list(
      elapsed_seconds = sum(vapply(
        gof_parts,
        function(part) as.numeric(part$diagnostics$elapsed_seconds %||% NA_real_),
        numeric(1)
      ), na.rm = TRUE),
      split_statistics = TRUE,
      statistics = statistics,
      per_statistic_elapsed_seconds = vapply(
        gof_parts,
        function(part) as.numeric(part$diagnostics$elapsed_seconds %||% NA_real_),
        numeric(1)
      ),
      bootstrap_method = unique(vapply(
        gof_parts,
        function(part) as.character(part$diagnostics$bootstrap_method %||% NA_character_),
        character(1)
      )),
      effective_bootstrap_method = unique(vapply(
        gof_parts,
        function(part) as.character(part$diagnostics$effective_bootstrap_method %||% NA_character_),
        character(1)
      )),
      fallback_reason = unique(vapply(
        gof_parts,
        function(part) as.character(part$diagnostics$fallback_reason %||% NA_character_),
        character(1)
      ))
    )
  )
  rm(gof_parts)
  gc(verbose = FALSE)

  get_gof_inference_value <- function(statistic_name, candidate_names) {
    inference <- gof_result$inference[[statistic_name]]
    if (is.null(inference)) {
      return(NA_real_)
    }
    for (candidate_name in candidate_names) {
      value <- NULL
      if (is.list(inference) && !is.null(inference[[candidate_name]])) {
        value <- inference[[candidate_name]]
      } else if (is.data.frame(inference) && candidate_name %in% names(inference)) {
        value <- inference[[candidate_name]][[1L]]
      }
      if (!is.null(value) && length(value) >= 1L && is.finite(as.numeric(value[[1L]]))) {
        return(as.numeric(value[[1L]]))
      }
    }
    NA_real_
  }

  theta_df <- rotmix_theta_summary_row(model_name, theta_hat)
  summary_df <- data.frame(
    dataset = dataset_label,
    model = model_name,
    n = n,
    loglik = loglik,
    AIC = aic,
    BIC = bic,
    projected_ks = projected_fit$ks,
    projected_cvm = projected_fit$cvm,
    gof_ks_observed = get_gof_inference_value("ks", c("observed", "observed_statistic", "statistic", "Tn")),
    gof_ks_p_value = get_gof_inference_value("ks", c("p_value", "p.value", "pvalue", "p")),
    gof_cvm_observed = get_gof_inference_value("cvm", c("observed", "observed_statistic", "statistic", "Tn")),
    gof_cvm_p_value = get_gof_inference_value("cvm", c("p_value", "p.value", "pvalue", "p")),
    B = as.integer(B),
    M = as.integer(M_value),
    n_cores = as.integer(n_cores),
    bootstrap_method = paste(gof_result$diagnostics$bootstrap_method, collapse = ";"),
    effective_bootstrap_method = paste(gof_result$diagnostics$effective_bootstrap_method, collapse = ";"),
    fallback_reason = paste(gof_result$diagnostics$fallback_reason[!is.na(gof_result$diagnostics$fallback_reason) &
      nzchar(gof_result$diagnostics$fallback_reason)], collapse = ";"),
    elapsed_seconds = gof_result$diagnostics$elapsed_seconds,
    stringsAsFactors = FALSE
  )
  summary_df <- cbind(summary_df, theta_df[, setdiff(names(theta_df), "model"), drop = FALSE])

  utils::write.csv(summary_df, file = file.path(output_dir, "summary.csv"), row.names = FALSE)
  utils::write.csv(theta_df, file = file.path(output_dir, "theta_hat.csv"), row.names = FALSE)
  utils::write.csv(density_grid_df, file = file.path(output_dir, "projected_density_grid.csv"), row.names = FALSE)
  utils::write.csv(cdf_grid_df, file = file.path(output_dir, "projected_cdf_grid.csv"), row.names = FALSE)
  utils::write.csv(data.frame(z = z), file = file.path(output_dir, "projected_data.csv"), row.names = FALSE)

  plot_paths <- plot_rotmix_projected_diagnostics(
    z = z,
    density_grid_df = density_grid_df,
    cdf_grid_df = cdf_grid_df,
    output_dir = output_dir,
    dataset_label = dataset_label,
    model_label = model_name
  )
  rm(gof_result, projected_fit, density_grid_df, cdf_grid_df, z)
  gc(verbose = FALSE)

  summary_df
}

run_comets_mixtures_short_long <- function(
  output_root = file.path(
    "output",
    "comets",
    paste0("rotational_mixtures_", timestamp_tag_rotmix_comets())
  ),
  datasets = c("short", "long"),
  models = c("beta_mixture2", "logitnormal_mixture2"),
  B = 500L,
  statistics = c("ks", "cvm"),
  n_cores = 6L,
  seed = 20260601L,
  M_value = 60L,
  ks_t_points = 200L,
  bootstrap_method = "reestimated",
  distance_type = "geodesic",
  control_beta = list(
    beta_mixture2_profile_method = "legendre",
    beta_mixture2_quad_n = 100L,
    beta_mixture2_optim_control = list(maxit = 350L, reltol = 1e-9)
  ),
  control_logit = list(
    logitnormal_mixture2_profile_method = "legendre",
    logitnormal_mixture2_quad_n = 400L,
    logitnormal_mixture2_optim_control = list(maxit = 350L, reltol = 1e-9)
  )
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  comets_data <- load_comets_distance_profile_data_rotmix()
  write_lines_if_possible_rotmix_comets(
    capture.output(sessionInfo()),
    file.path(output_root, "sessionInfo.txt")
  )

  summary_rows <- list()
  counter <- 0L

  for (dataset_name in datasets) {
    data_matrix <- if (identical(dataset_name, "short")) {
      as.matrix(comets_data$short$normal)
    } else {
      as.matrix(comets_data$long$normal)
    }
    dataset_label <- if (identical(dataset_name, "short")) "short_period" else "long_period"

    for (model_name in models) {
      counter <- counter + 1L
      model_output_dir <- file.path(output_root, sprintf("%02d_%s_%s", counter, dataset_label, model_name))
      control <- if (identical(model_name, "beta_mixture2")) control_beta else control_logit

      message(sprintf(
        "[Comets] %s / %s with B = %d, n_cores = %d, bootstrap_method = %s",
        dataset_label,
        model_name,
        as.integer(B),
        as.integer(n_cores),
        bootstrap_method
      ))

      summary_rows[[length(summary_rows) + 1L]] <- run_single_rotmix_comet_fit(
        data_matrix = data_matrix,
        dataset_label = dataset_label,
        model_name = model_name,
        output_dir = model_output_dir,
        B = as.integer(B),
        statistics = statistics,
        n_cores = as.integer(n_cores),
        seed = as.integer(seed) + counter,
        M_value = as.integer(M_value),
        ks_t_points = as.integer(ks_t_points),
        bootstrap_method = bootstrap_method,
        distance_type = distance_type,
        control = control
      )
    }
  }

  summary_all <- do.call(rbind, summary_rows)
  utils::write.csv(summary_all, file = file.path(output_root, "comets_rotational_mixtures_summary.csv"), row.names = FALSE)

  list(
    output_root = output_root,
    summary = summary_all
  )
}

parse_named_args_rotmix_comets <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
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

if (sys.nframe() == 0L) {
  args <- parse_named_args_rotmix_comets(commandArgs(trailingOnly = TRUE))

  result <- run_comets_mixtures_short_long(
    output_root = args$output_root %||% file.path(
      "output",
      "comets",
      paste0("rotational_mixtures_", timestamp_tag_rotmix_comets())
    ),
    datasets = if (!is.null(args$datasets)) {
      strsplit(tolower(args$datasets), ",", fixed = TRUE)[[1L]]
    } else {
      c("short", "long")
    },
    models = if (!is.null(args$models)) {
      strsplit(tolower(args$models), ",", fixed = TRUE)[[1L]]
    } else {
      c("beta_mixture2", "logitnormal_mixture2")
    },
    B = if (!is.null(args$B)) as.integer(args$B) else 500L,
    statistics = if (!is.null(args$statistics)) {
      strsplit(tolower(args$statistics), ",", fixed = TRUE)[[1L]]
    } else {
      c("ks", "cvm")
    },
    n_cores = if (!is.null(args$n_cores)) as.integer(args$n_cores) else 6L,
    seed = if (!is.null(args$seed)) as.integer(args$seed) else 20260601L,
    M_value = if (!is.null(args$M)) as.integer(args$M) else 60L,
    ks_t_points = if (!is.null(args$ks_t_points)) as.integer(args$ks_t_points) else 200L,
    bootstrap_method = args$bootstrap_method %||% "reestimated",
    distance_type = args$distance_type %||% "geodesic",
    control_beta = list(
      beta_mixture2_profile_method = "legendre",
      beta_mixture2_quad_n = 100L,
      beta_mixture2_optim_control = list(maxit = 350L, reltol = 1e-9),
      beta_mixture2_fast_shape_regular_eps = if (!is.null(args$beta_mixture2_fast_shape_regular_eps)) as.numeric(args$beta_mixture2_fast_shape_regular_eps) else 0,
      fast_bootstrap_chunk_size = if (!is.null(args$fast_bootstrap_chunk_size)) as.integer(args$fast_bootstrap_chunk_size) else NULL
    )
  )

  message(sprintf("Summary CSV: %s", file.path(result$output_root, "comets_rotational_mixtures_summary.csv")))
}
