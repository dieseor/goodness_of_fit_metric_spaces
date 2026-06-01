resolve_comets_jp_observed_method_path <- function(...) {
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

benchmark_runner_path_jp_observed <- resolve_comets_jp_observed_method_path(
  "scripts",
  "run_comets_distance_profile_jp_short_benchmark.R"
)
profile_refine_helper_path_jp_observed <- resolve_comets_jp_observed_method_path(
  "distance_profiles",
  "jp_profile_refine_observed_helpers.R"
)

source(benchmark_runner_path_jp_observed)
source(profile_refine_helper_path_jp_observed)

jp_mle_s2_weighted <- get("jp_mle_s2_weighted", mode = "function")
fit_jp_theta_profile_refine_observed_details <- get("fit_jp_theta_profile_refine_observed_details", mode = "function")
make_jp_ks_grid_comets <- get("make_jp_ks_grid_comets", mode = "function")
multiplier_bootstrap_gof <- get("multiplier_bootstrap_gof", mode = "function")
make_jp_spec <- get("make_jp_spec", mode = "function")
write_stage_bundle_jp_comets <- get("write_stage_bundle_jp_comets", mode = "function")
timestamp_tag_jp_comets <- get("timestamp_tag_jp_comets", mode = "function")
load_comets_distance_profile_data_jp <- get("load_comets_distance_profile_data_jp", mode = "function")
write_lines_if_possible_jp_comets <- get("write_lines_if_possible_jp_comets", mode = "function")
append_manifest_row_jp_comets <- get("append_manifest_row_jp_comets", mode = "function")

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "jp_mle_s2_weighted",
    "fit_jp_theta_profile_refine_observed_details",
    "make_jp_ks_grid_comets",
    "multiplier_bootstrap_gof",
    "make_jp_spec",
    "write_stage_bundle_jp_comets",
    "timestamp_tag_jp_comets",
    "load_comets_distance_profile_data_jp",
    "write_lines_if_possible_jp_comets",
    "append_manifest_row_jp_comets"
  ))
}

summarize_theta_star_ranges_jp_observed <- function(theta_star) {
  if (is.null(theta_star) || length(theta_star) == 0L) {
    return(list(
      kappa_star_min = NA_real_,
      kappa_star_max = NA_real_,
      psi_star_min = NA_real_,
      psi_star_max = NA_real_
    ))
  }

  list(
    kappa_star_min = min(vapply(theta_star, `[[`, numeric(1), "kappa")),
    kappa_star_max = max(vapply(theta_star, `[[`, numeric(1), "kappa")),
    psi_star_min = min(vapply(theta_star, `[[`, numeric(1), "psi")),
    psi_star_max = max(vapply(theta_star, `[[`, numeric(1), "psi"))
  )
}

summarize_bootstrap_quality_jp_observed <- function(result,
                                                    statistic) {
  bootstrap_values <- result$bootstrap$statistics[[statistic]]
  if (is.null(bootstrap_values) || length(bootstrap_values) == 0L) {
    return(list(
      stat_q50 = NA_real_,
      stat_q90 = NA_real_,
      stat_q95 = NA_real_,
      stat_q99 = NA_real_,
      bootstrap_failures = NA_integer_
    ))
  }

  stat_quantiles <- stats::quantile(
    bootstrap_values,
    probs = c(0.5, 0.9, 0.95, 0.99),
    na.rm = TRUE,
    names = FALSE
  )

  theta_star <- result$bootstrap$theta_star
  theta_failures <- if (is.null(theta_star)) {
    0L
  } else {
    sum(!vapply(theta_star, function(theta) {
      is.list(theta) &&
        all(is.finite(as.numeric(theta$mu))) &&
        is.finite(as.numeric(theta$kappa)) &&
        is.finite(as.numeric(theta$psi))
    }, logical(1)))
  }

  list(
    stat_q50 = as.numeric(stat_quantiles[[1L]]),
    stat_q90 = as.numeric(stat_quantiles[[2L]]),
    stat_q95 = as.numeric(stat_quantiles[[3L]]),
    stat_q99 = as.numeric(stat_quantiles[[4L]]),
    bootstrap_failures = as.integer(sum(!is.finite(bootstrap_values)) + theta_failures)
  )
}

fit_jp_observed_standard_bfgs <- function(data_matrix,
                                          control) {
  t0 <- proc.time()[["elapsed"]]
  theta <- jp_mle_s2_weighted(
    data = data_matrix,
    weights = NULL,
    control = control
  )
  list(
    theta = theta,
    elapsed_seconds = proc.time()[["elapsed"]] - t0,
    method_label = "standard_bfgs",
    details = list(convergence_status = "ok")
  )
}

fit_jp_observed_standard_nelder_mead <- function(data_matrix,
                                                 control) {
  t0 <- proc.time()[["elapsed"]]
  theta <- jp_mle_s2_weighted(
    data = data_matrix,
    weights = NULL,
    control = control
  )
  list(
    theta = theta,
    elapsed_seconds = proc.time()[["elapsed"]] - t0,
    method_label = "standard_nelder_mead",
    details = list(convergence_status = "ok")
  )
}

fit_jp_observed_profile_refine <- function(data_matrix,
                                           control) {
  fit <- fit_jp_theta_profile_refine_observed_details(
    data = data_matrix,
    weights = NULL,
    control = control
  )
  list(
    theta = fit$theta,
    elapsed_seconds = fit$elapsed_seconds,
    method_label = "profile_refine_observed",
    details = fit
  )
}

run_single_observed_method_stage_jp <- function(data_matrix,
                                                dataset_label,
                                                observed_method,
                                                statistic,
                                                B,
                                                M_value,
                                                n_cores,
                                                seed,
                                                stage_dir,
                                                distance_type = "geodesic",
                                                ks_t_points = 250L,
                                                bootstrap_control = list(
                                                  jp_profile_method = "tabulated",
                                                  jp_profile_n_u = 1025L,
                                                  jp_profile_n_delta = 257L,
                                                  jp_vmf_switch_abs_kappa_psi = 1e-3,
                                                  jp_mle_max_abs_kappa_psi = 6
                                                ),
                                                observed_control = list()) {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  stage_start <- Sys.time()

  observed_fit <- if (identical(observed_method, "standard_bfgs")) {
    fit_jp_observed_standard_bfgs(
      data_matrix = data_matrix,
      control = bootstrap_control
    )
  } else if (identical(observed_method, "standard_nelder_mead")) {
    nm_control <- bootstrap_control
    nm_control$jp_mle_method <- "Nelder-Mead"
    fit_jp_observed_standard_nelder_mead(
      data_matrix = data_matrix,
      control = nm_control
    )
  } else if (identical(observed_method, "profile_refine_observed_bootstrap_local")) {
    profile_control <- observed_control
    fit_jp_observed_profile_refine(
      data_matrix = data_matrix,
      control = profile_control
    )
  } else {
    stop(sprintf("Unsupported observed method: %s", observed_method))
  }

  engine_control <- bootstrap_control

  ks_grid <- if (identical(statistic, "ks")) {
    make_jp_ks_grid_comets(M_value = M_value, ks_t_points = ks_t_points)
  } else {
    NULL
  }

  result <- multiplier_bootstrap_gof(
    data = data_matrix,
    spec = make_jp_spec(distance_type = distance_type),
    null = list(type = "composite"),
    statistics = statistic,
    ks_grid = ks_grid,
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    observed_theta_hat = observed_fit$theta,
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = engine_control
  )

  theta_ranges <- summarize_theta_star_ranges_jp_observed(result$bootstrap$theta_star)
  bootstrap_quality <- summarize_bootstrap_quality_jp_observed(result, statistic)
  inference <- result$inference[[statistic]]
  theta_hat <- result$observed$theta_hat
  total_elapsed_seconds <- as.numeric(difftime(Sys.time(), stage_start, units = "secs"))

  summary_df <- data.frame(
    dataset = dataset_label,
    method_observed = observed_method,
    method_bootstrap = "standard_local_fast",
    statistic = statistic,
    B = as.integer(B),
    M = as.integer(M_value),
    n = result$diagnostics$n,
    n_cores = result$diagnostics$n_cores,
    observed_elapsed_seconds = observed_fit$elapsed_seconds,
    total_elapsed_seconds = total_elapsed_seconds,
    observed_statistic = inference$observed,
    critical_value = inference$critical_value,
    p_value = inference$p_value,
    reject = inference$reject,
    mu_hat_1 = theta_hat$mu[[1L]],
    mu_hat_2 = theta_hat$mu[[2L]],
    mu_hat_3 = theta_hat$mu[[3L]],
    kappa_hat = theta_hat$kappa,
    psi_hat = theta_hat$psi,
    bootstrap_stat_q50 = bootstrap_quality$stat_q50,
    bootstrap_stat_q90 = bootstrap_quality$stat_q90,
    bootstrap_stat_q95 = bootstrap_quality$stat_q95,
    bootstrap_stat_q99 = bootstrap_quality$stat_q99,
    kappa_star_min = theta_ranges$kappa_star_min,
    kappa_star_max = theta_ranges$kappa_star_max,
    psi_star_min = theta_ranges$psi_star_min,
    psi_star_max = theta_ranges$psi_star_max,
    bootstrap_failures = bootstrap_quality$bootstrap_failures,
    output_dir = stage_dir,
    stringsAsFactors = FALSE
  )

  bundle <- list(
    dataset_label = dataset_label,
    observed_method = observed_method,
    statistic = statistic,
    observed_fit = observed_fit,
    result = result,
    summary = summary_df,
    config = list(
      B = as.integer(B),
      M = as.integer(M_value),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      distance_type = distance_type,
      ks_t_points = as.integer(ks_t_points)
    ),
    diagnostics = list(
      observed_elapsed_seconds = observed_fit$elapsed_seconds,
      total_elapsed_seconds = total_elapsed_seconds,
      observed_method = observed_method,
      bootstrap_method = "standard_local_fast"
    )
  )

  write_stage_bundle_jp_comets(
    stage_dir = stage_dir,
    bundle = bundle,
    summary_df = summary_df,
    metadata = list(
      dataset_label = dataset_label,
      observed_method = observed_method,
      statistic = statistic,
      B = as.integer(B),
      M = as.integer(M_value),
      n = nrow(data_matrix),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      distance_type = distance_type,
      ks_t_points = as.integer(ks_t_points),
      observed_elapsed_seconds = observed_fit$elapsed_seconds,
      total_elapsed_seconds = total_elapsed_seconds,
      bootstrap_failures = bootstrap_quality$bootstrap_failures,
      completed_at = Sys.time()
    ),
    checkpoint_name = "stage_bundle"
  )

  utils::write.csv(
    summary_df,
    file = file.path(stage_dir, "timing.csv"),
    row.names = FALSE
  )

  bundle
}

build_method_comparison_table_jp <- function(summary_df) {
  standard_df <- summary_df[summary_df$method_observed == "standard_nelder_mead", , drop = FALSE]
  profile_df <- summary_df[summary_df$method_observed == "profile_refine_observed_bootstrap_local", , drop = FALSE]

  if (nrow(standard_df) == 0L || nrow(profile_df) == 0L) {
    return(data.frame())
  }

  merged <- merge(
    standard_df,
    profile_df,
    by = c("dataset", "statistic", "B", "M", "n", "n_cores"),
    suffixes = c("_nm", "_profile"),
    all = FALSE
  )

  merged$observed_time_ratio_profile_over_nm <- merged$observed_elapsed_seconds_profile / merged$observed_elapsed_seconds_nm
  merged$total_time_ratio_profile_over_nm <- merged$total_elapsed_seconds_profile / merged$total_elapsed_seconds_nm
  merged
}

run_comets_distance_profile_jp_short_observed_method_benchmark <- function(output_root = NULL,
                                                                           observed_methods = c("standard_bfgs", "standard_nelder_mead", "profile_refine_observed_bootstrap_local"),
                                                                           B_values = c(50L, 100L, 200L, 500L),
                                                                           statistic = "ks",
                                                                           n_cores = 12L,
                                                                           ks_t_points = 250L,
                                                                           base_seed = 20260530L,
                                                                           distance_type = "geodesic",
                                                                           bootstrap_control = list(
                                                                             jp_profile_method = "tabulated",
                                                                             jp_profile_n_u = 1025L,
                                                                             jp_profile_n_delta = 257L,
                                                                             jp_vmf_switch_abs_kappa_psi = 1e-3,
                                                                             jp_mle_max_abs_kappa_psi = 6
                                                                           ),
                                                                           profile_refine_control = list(
                                                                             max_abs_kappa_psi = 6,
                                                                             kappa_grid = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0),
                                                                             psi_abs_grid = c(0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0),
                                                                             method_mu = "Nelder-Mead",
                                                                             maxit_mu = 150L,
                                                                             reltol_mu = 1e-8,
                                                                             top_k_per_branch = 3L,
                                                                             psi_min = 1e-3,
                                                                             maxit_refine = 500L,
                                                                             reltol_refine = 1e-10
                                                                           )) {
  if (is.null(output_root)) {
    output_root <- file.path(
      "output",
      "comets_distance_profile_jp",
      paste0("run_", timestamp_tag_jp_comets(), "_short_", statistic, "_observed_method_benchmark")
    )
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  B_values <- unique(as.integer(B_values))
  B_values <- B_values[is.finite(B_values) & B_values > 0L]
  if (length(B_values) == 0L) {
    stop("`B_values` must contain at least one positive integer.")
  }

  comets_data <- load_comets_distance_profile_data_jp()
  dataset_summary <- data.frame(
    dataset = "short_period",
    n = nrow(comets_data$short),
    ambient_dim = ncol(comets_data$short$normal),
    stringsAsFactors = FALSE
  )
  utils::write.csv(dataset_summary, file = file.path(output_root, "dataset_summary.csv"), row.names = FALSE)
  write_lines_if_possible_jp_comets(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))

  summary_rows <- list()
  manifest_rows <- list()
  stage_results <- list()
  idx <- 1L

  for (observed_method in observed_methods) {
    method_dir <- file.path(output_root, observed_method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

    for (i in seq_along(B_values)) {
      B_value <- B_values[[i]]
      stage_id <- sprintf("%02d", i)
      stage_name <- sprintf("short_%s_M%d_B%d", statistic, B_value, B_value)
      stage_dir <- file.path(method_dir, paste0(stage_id, "_", stage_name))
      stage_start <- Sys.time()

      message(sprintf(
        "[Short-period JP %s | %s] %d/%d: M = B = %d with %d cores",
        toupper(statistic),
        observed_method,
        i,
        length(B_values),
        B_value,
        as.integer(n_cores)
      ))

      bundle <- run_single_observed_method_stage_jp(
        data_matrix = comets_data$short$normal,
        dataset_label = sprintf("Short-period JP %s", toupper(statistic)),
        observed_method = observed_method,
        statistic = statistic,
        B = B_value,
        M_value = B_value,
        n_cores = n_cores,
        seed = base_seed + idx,
        stage_dir = stage_dir,
        distance_type = distance_type,
        ks_t_points = ks_t_points,
        bootstrap_control = bootstrap_control,
        observed_control = profile_refine_control
      )
      stage_results[[paste(observed_method, stage_name, sep = "::")]] <- bundle
      summary_rows[[length(summary_rows) + 1L]] <- bundle$summary

      manifest_rows <- append_manifest_row_jp_comets(
        manifest_rows = manifest_rows,
        stage_id = paste(observed_method, stage_id, sep = "::"),
        stage_label = sprintf("%s | Short-period JP %s M=B=%d", observed_method, toupper(statistic), B_value),
        status = "completed",
        stage_dir = stage_dir,
        started_at = stage_start,
        finished_at = Sys.time()
      )

      summary_df <- do.call(rbind, summary_rows)
      utils::write.csv(summary_df, file = file.path(output_root, "benchmark_summary.csv"), row.names = FALSE)
      utils::write.csv(do.call(rbind, manifest_rows), file = file.path(output_root, "pipeline_manifest.csv"), row.names = FALSE)
      comparison_df <- build_method_comparison_table_jp(summary_df)
      if (nrow(comparison_df) > 0L) {
        utils::write.csv(comparison_df, file = file.path(output_root, "method_comparison.csv"), row.names = FALSE)
      }

      idx <- idx + 1L
    }
  }

  pipeline_result <- list(
    output_root = output_root,
    dataset_summary = dataset_summary,
    stages = stage_results,
    manifest = do.call(rbind, manifest_rows),
    benchmark_summary = do.call(rbind, summary_rows)
  )
  saveRDS(pipeline_result, file = file.path(output_root, "pipeline_result.rds"))
  pipeline_result
}

parse_named_args_jp_observed_methods <- function(args) {
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

if (sys.nframe() == 0L) {
  args <- parse_named_args_jp_observed_methods(commandArgs(trailingOnly = TRUE))

  B_values <- if (!is.null(args$B_values)) {
    as.integer(strsplit(args$B_values, ",", fixed = TRUE)[[1L]])
  } else {
    c(50L, 100L, 200L, 500L)
  }

  observed_methods <- if (!is.null(args$observed_methods)) {
    strsplit(args$observed_methods, ",", fixed = TRUE)[[1L]]
  } else {
    c("standard_bfgs", "standard_nelder_mead", "profile_refine_observed_bootstrap_local")
  }

  output_root <- if (!is.null(args$output_root)) args$output_root else NULL
  statistic <- if (!is.null(args$statistic)) tolower(args$statistic) else "ks"
  n_cores <- if (!is.null(args$n_cores)) as.integer(args$n_cores) else 12L
  ks_t_points <- if (!is.null(args$ks_t_points)) as.integer(args$ks_t_points) else 250L
  base_seed <- if (!is.null(args$seed)) as.integer(args$seed) else 20260530L
  distance_type <- if (!is.null(args$distance_type)) args$distance_type else "geodesic"

  run_comets_distance_profile_jp_short_observed_method_benchmark(
    output_root = output_root,
    observed_methods = observed_methods,
    B_values = B_values,
    statistic = statistic,
    n_cores = n_cores,
    ks_t_points = ks_t_points,
    base_seed = base_seed,
    distance_type = distance_type
  )
}