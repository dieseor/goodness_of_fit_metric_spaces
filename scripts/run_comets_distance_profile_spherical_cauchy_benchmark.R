resolve_comets_spherical_cauchy_path <- function(...) {
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

model_specs_script_path_sc_comets <- resolve_comets_spherical_cauchy_path(
  "bootstrap",
  "model_specs.R"
)
multiplier_bootstrap_script_path_sc_comets <- resolve_comets_spherical_cauchy_path(
  "bootstrap",
  "multiplier_bootstrap.R"
)
utils_script_path_sc_comets <- resolve_comets_spherical_cauchy_path("utils.R")
comets_utils_script_path_sc_comets <- resolve_comets_spherical_cauchy_path(
  "real_data",
  "comets",
  "utils_comets_data.R"
)

source(model_specs_script_path_sc_comets)
source(multiplier_bootstrap_script_path_sc_comets)
source(utils_script_path_sc_comets)
source(comets_utils_script_path_sc_comets)
source(resolve_comets_spherical_cauchy_path("scripts", "path_helpers.R"))

make_spherical_cauchy_spec <- get("make_spherical_cauchy_spec", mode = "function")
multiplier_bootstrap_gof <- get("multiplier_bootstrap_gof", mode = "function")
generate_canonical_lattice <- get("generate_canonical_lattice", mode = "function")

timestamp_tag_sc_comets <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

write_lines_if_possible_sc_comets <- function(lines, path) {
  writeLines(as.character(lines), con = path)
  invisible(path)
}

load_existing_sc_benchmark_summary <- function(output_root) {
  summary_path <- file.path(output_root, "benchmark_summary.csv")
  if (!file.exists(summary_path)) {
    return(NULL)
  }

  utils::read.csv(summary_path, stringsAsFactors = FALSE)
}

write_stage_bundle_sc_comets <- function(stage_dir,
                                         bundle,
                                         summary_df = NULL,
                                         metadata = NULL,
                                         checkpoint_name = "stage_bundle") {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  bundle_rds <- file.path(stage_dir, paste0(checkpoint_name, ".rds"))
  bundle_rdata <- file.path(stage_dir, paste0(checkpoint_name, ".RData"))
  saveRDS(bundle, file = bundle_rds)

  stage_bundle <- bundle
  save(stage_bundle, file = bundle_rdata)

  if (!is.null(summary_df)) {
    utils::write.csv(
      summary_df,
      file = file.path(stage_dir, "pvalues.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(metadata)) {
    saveRDS(metadata, file = file.path(stage_dir, "metadata.rds"))
    write_lines_if_possible_sc_comets(
      capture.output(str(metadata, max.level = 2)),
      file.path(stage_dir, "metadata.txt")
    )
  }

  invisible(
    list(
      bundle_rds = bundle_rds,
      bundle_rdata = bundle_rdata
    )
  )
}

load_comets_distance_profile_data_sc <- function() {
  load_comets_real_data(
    finite_normals = "both",
    warn_incomplete = TRUE,
    warn_nonfinite = TRUE
  )
}

append_manifest_row_sc_comets <- function(manifest_rows,
                                          stage_id,
                                          stage_label,
                                          status,
                                          stage_dir,
                                          started_at,
                                          finished_at) {
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    stage_id = stage_id,
    stage_label = stage_label,
    status = status,
    stage_dir = stage_dir,
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
    finished_at = format(finished_at, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_seconds = as.numeric(difftime(finished_at, started_at, units = "secs")),
    stringsAsFactors = FALSE
  )
  manifest_rows
}

summarize_sc_model_result <- function(result,
                                      dataset_label,
                                      statistic,
                                      M_value) {
  theta_hat <- result$observed$theta_hat
  inference <- result$inference[[statistic]]
  theta_star <- result$bootstrap$theta_star

  rho_star_range <- if (is.null(theta_star) || length(theta_star) == 0L) {
    c(NA_real_, NA_real_)
  } else {
    range(vapply(theta_star, `[[`, numeric(1), "rho"))
  }

  data.frame(
    dataset = dataset_label,
    model = "spherical_cauchy_composite",
    null_type = "composite",
    statistic = statistic,
    spec_name = result$diagnostics$spec_name,
    engine = result$diagnostics$engine,
    method = result$diagnostics$method,
    weighted_mle = isTRUE(result$diagnostics$weighted_mle),
    observed_statistic = inference$observed,
    critical_value = inference$critical_value,
    p_value = inference$p_value,
    reject = inference$reject,
    mu_hat_1 = theta_hat$mu[[1L]],
    mu_hat_2 = theta_hat$mu[[2L]],
    mu_hat_3 = theta_hat$mu[[3L]],
    rho_hat = theta_hat$rho,
    rho_star_min = rho_star_range[[1L]],
    rho_star_max = rho_star_range[[2L]],
    B = result$bootstrap$B,
    M = NA_integer_,
    ks_grid_mode = if (identical(statistic, "ks")) "sample_points_unique_distances" else NA_character_,
    n = result$diagnostics$n,
    n_cores = result$diagnostics$n_cores,
    bootstrap_method = result$diagnostics$bootstrap_method %||% NA_character_,
    effective_bootstrap_method = result$diagnostics$effective_bootstrap_method %||% NA_character_,
    elapsed_seconds = result$diagnostics$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

make_sc_ks_grid_comets <- function(M_value,
                                   ks_t_points = 250L) {
  make_sample_unique_distance_ks_grid()
}

run_single_sc_comet_model <- function(data_matrix,
                                      statistic,
                                      B,
                                      M_value,
                                      n_cores,
                                      seed,
                                      bootstrap_method = "reestimated",
                                      distance_type = "geodesic",
                                      ks_t_points = 250L,
                                      control = list(
                                        spherical_cauchy_maxit = 500L,
                                        spherical_cauchy_reltol = 1e-10,
                                        spherical_cauchy_optim_method = "BFGS",
                                        spherical_cauchy_use_gradient = TRUE,
                                        spherical_cauchy_profile_tol = 1e-10,
                                        spherical_cauchy_profile_warn = FALSE
                                      )) {
  spec <- make_spherical_cauchy_spec(distance_type = distance_type)

  ks_grid <- if (identical(statistic, "ks")) {
    make_sc_ks_grid_comets(M_value = M_value, ks_t_points = ks_t_points)
  } else {
    NULL
  }

  multiplier_bootstrap_gof(
    data = data_matrix,
    spec = spec,
    null = list(type = "composite"),
    statistics = statistic,
    ks_grid = ks_grid,
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    bootstrap_method = bootstrap_method,
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = control
  )
}

run_sc_comet_stage <- function(data_matrix,
                               dataset_label,
                               statistic,
                               stage_dir,
                               B,
                               M_value,
                               n_cores,
                               seed,
                               bootstrap_method = "reestimated",
                               distance_type = "geodesic",
                               ks_t_points = 250L,
                               control = list(
                                 spherical_cauchy_maxit = 500L,
                                 spherical_cauchy_reltol = 1e-10,
                                 spherical_cauchy_optim_method = "BFGS",
                                 spherical_cauchy_use_gradient = TRUE,
                                 spherical_cauchy_profile_tol = 1e-10,
                                 spherical_cauchy_profile_warn = FALSE
                               )) {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  stage_start <- Sys.time()

  result <- run_single_sc_comet_model(
    data_matrix = data_matrix,
    statistic = statistic,
    B = B,
    M_value = M_value,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    distance_type = distance_type,
    ks_t_points = ks_t_points,
    control = control
  )

  summary_df <- summarize_sc_model_result(
    result = result,
    dataset_label = dataset_label,
    statistic = statistic,
    M_value = M_value
  )

  elapsed_seconds <- as.numeric(difftime(Sys.time(), stage_start, units = "secs"))

  final_bundle <- list(
    dataset_label = dataset_label,
    statistic = statistic,
    result = result,
    summary = summary_df,
    config = list(
      B = as.integer(B),
      M = NA_integer_,
      ks_grid_mode = if (identical(statistic, "ks")) "sample_points_unique_distances" else NA_character_,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      bootstrap_method = bootstrap_method,
      distance_type = distance_type,
      ks_t_points = as.integer(ks_t_points)
    ),
    diagnostics = list(
      n = nrow(data_matrix),
      p = ncol(data_matrix),
      elapsed_seconds = elapsed_seconds,
      engine = "multiplier_bootstrap_gof",
      method = "distance_profiles",
      weighted_mle = TRUE
    )
  )

  write_stage_bundle_sc_comets(
    stage_dir = stage_dir,
    bundle = final_bundle,
    summary_df = summary_df,
    metadata = list(
      dataset_label = dataset_label,
      statistic = statistic,
      B = as.integer(B),
      M = NA_integer_,
      ks_grid_mode = if (identical(statistic, "ks")) "sample_points_unique_distances" else NA_character_,
      n = nrow(data_matrix),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      bootstrap_method = bootstrap_method,
      distance_type = distance_type,
      ks_t_points = as.integer(ks_t_points),
      engine = "multiplier_bootstrap_gof",
      method = "distance_profiles",
      weighted_mle = TRUE,
      completed_at = Sys.time(),
      elapsed_seconds = elapsed_seconds
    ),
    checkpoint_name = "stage_bundle"
  )

  utils::write.csv(
    summary_df,
    file = file.path(stage_dir, "timing.csv"),
    row.names = FALSE
  )

  final_bundle
}

run_comets_distance_profile_spherical_cauchy_benchmark <- function(output_root = NULL,
                                                                    dataset = c("short", "long"),
                                                                    B_values = c(10L, 30L, 50L, 100L, 200L, 500L, 1000L),
                                                                    statistic = "ks",
                                                                    n_cores = 12L,
                                                                    ks_t_points = 250L,
                                                                    base_seed = 20260529L,
                                                                    bootstrap_method = "reestimated",
                                                                    distance_type = "geodesic",
                                                                    control = list(
                                                                      spherical_cauchy_maxit = 500L,
                                                                      spherical_cauchy_reltol = 1e-10,
                                                                      spherical_cauchy_optim_method = "BFGS",
                                                                      spherical_cauchy_use_gradient = TRUE,
                                                                      spherical_cauchy_profile_tol = 1e-10,
                                                                      spherical_cauchy_profile_warn = FALSE
                                                                    )) {
  dataset <- match.arg(dataset)

  if (is.null(output_root)) {
    speed_dir <- if (identical(bootstrap_method, "fast_multiplier")) "fast" else "slow"
    run_name <- sprintf("%s_period_%s_benchmark", dataset, statistic)
    output_root <- canonical_comets_spherical_cauchy_dir(run_name, speed_dir)
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  B_values <- unique(as.integer(B_values))
  B_values <- B_values[is.finite(B_values) & B_values > 0L]
  if (length(B_values) == 0L) {
    stop("`B_values` must contain at least one positive integer.")
  }

  comets_data <- load_comets_distance_profile_data_sc()
  data_matrix <- if (identical(dataset, "short")) comets_data$short$normal else comets_data$long$normal
  dataset_label <- if (identical(dataset, "short")) "short_period" else "long_period"

  write_comets_dataset_summary(
    path = file.path(output_root, "dataset_summary.csv"),
    dataset_label = dataset_label,
    data_matrix = data_matrix
  )

  write_lines_if_possible_sc_comets(
    capture.output(sessionInfo()),
    file.path(output_root, "sessionInfo.txt")
  )

  saveRDS(
    list(
      output_root = output_root,
      benchmark_B_values = B_values,
      benchmark_M_values = B_values,
      statistic = statistic,
      n_cores = as.integer(n_cores),
      ks_t_points = as.integer(ks_t_points),
      base_seed = as.integer(base_seed),
      bootstrap_method = bootstrap_method,
      distance_type = distance_type,
      dataset = dataset_label,
      model = "spherical_cauchy_composite"
    ),
    file = file.path(output_root, "run_config.rds")
  )

  existing_summary <- load_existing_sc_benchmark_summary(output_root)
  manifest_rows <- list()
  stage_results <- list()
  summary_rows <- list()

  for (i in seq_along(B_values)) {
    B_value <- B_values[[i]]
    stage_id <- sprintf("%02d", i)
    stage_name <- sprintf("%s_%s_M%d_B%d", dataset, statistic, B_value, B_value)
    stage_dir <- file.path(output_root, paste0(stage_id, "_", stage_name))
    stage_file <- file.path(stage_dir, "stage_bundle.rds")

    if (file.exists(stage_file)) {
      existing_row <- NULL
      if (!is.null(existing_summary)) {
        row_idx <- which(
          existing_summary$dataset == dataset_label &
            existing_summary$statistic == statistic &
            existing_summary$B == B_value
        )
        if (length(row_idx) >= 1L) {
          existing_row <- existing_summary[row_idx[[1L]], , drop = FALSE]
        }
      }

      if (is.null(existing_row)) {
        bundle_existing <- readRDS(stage_file)
        existing_row <- bundle_existing$summary
      }

      summary_rows[[length(summary_rows) + 1L]] <- existing_row
      next
    }

    stage_start <- Sys.time()

    message(sprintf(
      "[%s %s] %d/%d: M = B = %d with %d cores",
      if (identical(dataset, "short")) "Short-period SC" else "Long-period SC",
      toupper(statistic),
      i,
      length(B_values),
      B_value,
      as.integer(n_cores)
    ))

    stage_results[[stage_name]] <- run_sc_comet_stage(
      data_matrix = data_matrix,
      dataset_label = dataset_label,
      statistic = statistic,
      stage_dir = stage_dir,
      B = B_value,
      M_value = B_value,
      n_cores = n_cores,
      seed = base_seed + i,
      bootstrap_method = bootstrap_method,
      distance_type = distance_type,
      ks_t_points = ks_t_points,
      control = control
    )

    manifest_rows <- append_manifest_row_sc_comets(
      manifest_rows = manifest_rows,
      stage_id = stage_id,
      stage_label = sprintf("%s spherical_cauchy %s M=B=%d", dataset_label, toupper(statistic), B_value),
      status = "completed",
      stage_dir = stage_dir,
      started_at = stage_start,
      finished_at = Sys.time()
    )

    utils::write.csv(
      do.call(rbind, manifest_rows),
      file = file.path(output_root, "pipeline_manifest.csv"),
      row.names = FALSE
    )

    summary_rows[[length(summary_rows) + 1L]] <- stage_results[[stage_name]]$summary
    benchmark_summary <- do.call(rbind, summary_rows)
    utils::write.csv(
      benchmark_summary,
      file = file.path(output_root, "benchmark_summary.csv"),
      row.names = FALSE
    )
  }

  pipeline_result <- list(
    output_root = output_root,
    dataset_summary = data.frame(
      dataset = dataset_label,
      n = nrow(data_matrix),
      ambient_dim = ncol(data_matrix),
      stringsAsFactors = FALSE
    ),
    stages = stage_results,
    manifest = do.call(rbind, manifest_rows),
    engine = "multiplier_bootstrap_gof",
    method = "distance_profiles"
  )

  saveRDS(pipeline_result, file = file.path(output_root, "pipeline_result.rds"))
  pipeline_result
}

parse_named_args_sc_comets <- function(args) {
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
  args <- parse_named_args_sc_comets(commandArgs(trailingOnly = TRUE))

  dataset <- if (!is.null(args$dataset)) tolower(args$dataset) else "short"
  B_values <- if (!is.null(args$B_values)) {
    as.integer(strsplit(args$B_values, ",", fixed = TRUE)[[1L]])
  } else {
    c(10L, 30L, 50L, 100L, 200L, 500L, 1000L)
  }

  output_root <- if (!is.null(args$output_root)) args$output_root else NULL
  statistic <- if (!is.null(args$statistic)) tolower(args$statistic) else "ks"
  n_cores <- if (!is.null(args$n_cores)) as.integer(args$n_cores) else 12L
  ks_t_points <- if (!is.null(args$ks_t_points)) as.integer(args$ks_t_points) else 250L
  base_seed <- if (!is.null(args$seed)) as.integer(args$seed) else 20260529L
  bootstrap_method <- if (!is.null(args$bootstrap_method)) args$bootstrap_method else "reestimated"
  distance_type <- if (!is.null(args$distance_type)) args$distance_type else "geodesic"

  run_comets_distance_profile_spherical_cauchy_benchmark(
    output_root = output_root,
    dataset = dataset,
    B_values = B_values,
    statistic = statistic,
    n_cores = n_cores,
    ks_t_points = ks_t_points,
    base_seed = base_seed,
    bootstrap_method = bootstrap_method,
    distance_type = distance_type
  )
}
