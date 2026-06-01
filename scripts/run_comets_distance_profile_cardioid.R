resolve_comets_cardioid_path <- function(...) {
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

cardioid_model_spec_script_path <- resolve_comets_cardioid_path(
  "bootstrap",
  "cardioid_model_spec.R"
)
multiplier_bootstrap_script_path <- resolve_comets_cardioid_path(
  "bootstrap",
  "multiplier_bootstrap.R"
)
utils_script_path_cardioid <- resolve_comets_cardioid_path("utils.R")

source(cardioid_model_spec_script_path)
source(multiplier_bootstrap_script_path)
source(utils_script_path_cardioid)

timestamp_tag_cardioid_comets <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

write_lines_if_possible_cardioid <- function(lines, path) {
  writeLines(as.character(lines), con = path)
  invisible(path)
}

write_stage_bundle_cardioid <- function(stage_dir,
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
    write_lines_if_possible_cardioid(
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

load_comets_distance_profile_data <- function() {
  if (!requireNamespace("sphunif", quietly = TRUE)) {
    stop("Package `sphunif` is required for the comet analyses.")
  }

  data("comets", package = "sphunif")
  comets$normal <- cbind(
    sin(comets$i) * sin(comets$om),
    -sin(comets$i) * cos(comets$om),
    cos(comets$i)
  )

  valid_rows <-
    !is.na(comets$class) &
    is.finite(comets$per_y) &
    is.finite(comets$i) &
    is.finite(comets$om) &
    !is.na(comets$frag)

  dropped_incomplete <- sum(!valid_rows)
  if (dropped_incomplete > 0L) {
    warning(
      sprintf(
        "Dropping %d comet rows with incomplete orbital elements before cardioid comet filtering.",
        dropped_incomplete
      ),
      call. = FALSE
    )
  }

  comets_valid <- comets[valid_rows, , drop = FALSE]

  comets_oort <- subset(
    x = comets_valid,
    subset = !(class %in% c("HYP", "PAR")) & per_y >= 200
  )
  comets_oort <- comets_oort[!comets_oort$frag, ]

  comets_short <- subset(
    x = comets_valid,
    subset = !(class %in% c("HYP", "PAR")) & per_y < 200
  )
  comets_short <- comets_short[!comets_short$frag, ]

  list(
    raw = comets,
    oort = comets_oort,
    short = comets_short
  )
}

make_comet_cardioid_models <- function(model_ids = NULL) {
  models <- list(
    list(
      id = "Uniform",
      label = "Uniform",
      k = 1L,
      null = list(
        type = "simple",
        theta = list(mu = c(0, 0, 1), rho = 0, k = 1L)
      )
    ),
    list(
      id = "C1",
      label = "C1",
      k = 1L,
      null = list(type = "composite")
    ),
    list(
      id = "C2",
      label = "C2",
      k = 2L,
      null = list(type = "composite")
    ),
    list(
      id = "C3",
      label = "C3",
      k = 3L,
      null = list(type = "composite")
    ),
    list(
      id = "C4",
      label = "C4",
      k = 4L,
      null = list(type = "composite")
    )
  )

  if (is.null(model_ids)) {
    return(models)
  }

  keep_ids <- as.character(model_ids)
  output <- models[vapply(models, function(model) model$id %in% keep_ids, logical(1))]
  if (length(output) == 0L) {
    stop("`model_ids` did not match any supported cardioid comet models.")
  }

  output
}

append_manifest_row_cardioid <- function(manifest_rows,
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

summarize_cardioid_model_result <- function(result,
                                            model,
                                            dataset_label,
                                            statistic) {
  theta_hat <- result$observed$theta_hat
  inference <- result$inference[[statistic]]
  theta_star <- result$bootstrap$theta_star

  rho_star_range <- if (is.null(theta_star) || length(theta_star) == 0L) {
    c(NA_real_, NA_real_)
  } else {
    rho_values <- vapply(theta_star, `[[`, numeric(1), "rho")
    range(rho_values)
  }

  mu1_star_range <- if (is.null(theta_star) || length(theta_star) == 0L) {
    c(NA_real_, NA_real_)
  } else {
    mu1_values <- vapply(theta_star, function(theta) theta$mu[[1L]], numeric(1))
    range(mu1_values)
  }

  data.frame(
    dataset = dataset_label,
    model = model$label,
    k = as.integer(model$k),
    null_type = model$null$type,
    statistic = statistic,
    spec_name = result$diagnostics$spec_name,
    engine = result$diagnostics$engine,
    method = result$diagnostics$method,
    weighted_mle = isTRUE(result$diagnostics$weighted_mle),
    observed_statistic = inference$observed,
    critical_value = inference$critical_value,
    p_value = inference$p_value,
    reject = inference$reject,
    rho_hat = theta_hat$rho,
    mu_hat_1 = theta_hat$mu[[1L]],
    mu_hat_2 = theta_hat$mu[[2L]],
    mu_hat_3 = theta_hat$mu[[3L]],
    rho_star_min = rho_star_range[[1L]],
    rho_star_max = rho_star_range[[2L]],
    mu_star_1_min = mu1_star_range[[1L]],
    mu_star_1_max = mu1_star_range[[2L]],
    B = result$bootstrap$B,
    n = result$diagnostics$n,
    n_cores = result$diagnostics$n_cores,
    elapsed_seconds = result$diagnostics$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

run_single_cardioid_comet_model <- function(data_matrix,
                                            model,
                                            statistic,
                                            B,
                                            n_cores,
                                            seed,
                                            ks_grid = NULL,
                                            distance_type = "geodesic",
                                            control = list()) {
  spec <- make_cardioid_spec(
    k = model$k,
    distance_type = distance_type,
    unknown_param = "both"
  )

  multiplier_bootstrap_gof(
    data = data_matrix,
    spec = spec,
    null = model$null,
    statistics = statistic,
    ks_grid = ks_grid,
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = identical(model$null$type, "composite")
    ),
    control = control
  )
}

run_cardioid_comet_stage <- function(data_matrix,
                                     dataset_label,
                                     statistic,
                                     stage_dir,
                                     models,
                                     B,
                                     n_cores,
                                     seed,
                                     ks_grid = NULL,
                                     distance_type = "geodesic",
                                     control = list()) {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  stage_start <- Sys.time()
  stage_pb <- utils::txtProgressBar(min = 0, max = length(models), style = 3)
  on.exit(close(stage_pb), add = TRUE)

  results_by_model <- vector("list", length(models))
  names(results_by_model) <- vapply(models, `[[`, character(1), "id")
  summary_rows <- vector("list", length(models))

  for (i in seq_along(models)) {
    model <- models[[i]]
    message(sprintf(
      "[%s][%s] Running %d/%d: %s",
      dataset_label,
      statistic,
      i,
      length(models),
      model$label
    ))

    result <- run_single_cardioid_comet_model(
      data_matrix = data_matrix,
      model = model,
      statistic = statistic,
      B = B,
      n_cores = n_cores,
      seed = seed + i,
      ks_grid = ks_grid,
      distance_type = distance_type,
      control = control
    )

    results_by_model[[i]] <- result
    summary_rows[[i]] <- summarize_cardioid_model_result(
      result = result,
      model = model,
      dataset_label = dataset_label,
      statistic = statistic
    )

    partial_bundle <- list(
      dataset_label = dataset_label,
      statistic = statistic,
      completed_models = results_by_model[seq_len(i)],
      summary = do.call(rbind, summary_rows[seq_len(i)]),
      config = list(
        B = as.integer(B),
        n_cores = as.integer(n_cores),
        seed = as.integer(seed),
        distance_type = distance_type,
        ks_grid = ks_grid
      ),
      engine = "multiplier_bootstrap_gof",
      method = "distance_profiles"
    )

    write_stage_bundle_cardioid(
      stage_dir = stage_dir,
      bundle = partial_bundle,
      summary_df = partial_bundle$summary,
      metadata = list(
        dataset_label = dataset_label,
        statistic = statistic,
        status = "partial",
        completed_models = i,
        total_models = length(models),
        engine = "multiplier_bootstrap_gof",
        method = "distance_profiles",
        updated_at = Sys.time()
      ),
      checkpoint_name = "stage_checkpoint"
    )

    utils::setTxtProgressBar(stage_pb, i)
  }

  summary_df <- do.call(rbind, summary_rows)
  elapsed_seconds <- as.numeric(difftime(Sys.time(), stage_start, units = "secs"))

  final_bundle <- list(
    dataset_label = dataset_label,
    statistic = statistic,
    results = results_by_model,
    summary = summary_df,
    config = list(
      B = as.integer(B),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      distance_type = distance_type,
      ks_grid = ks_grid
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

  write_stage_bundle_cardioid(
    stage_dir = stage_dir,
    bundle = final_bundle,
    summary_df = summary_df,
    metadata = list(
      dataset_label = dataset_label,
      statistic = statistic,
      B = as.integer(B),
      n = nrow(data_matrix),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      distance_type = distance_type,
      ks_grid = ks_grid,
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

run_comets_distance_profile_cardioid <- function(output_root = NULL,
                                                 stages = c("oort_cvm", "short_cvm", "short_ks"),
                                                 model_ids = NULL,
                                                 cvm_B = 1000L,
                                                 ks_B = 200L,
                                                 ks_omega_points = 200L,
                                                 ks_t_points = 50L,
                                                 n_cores = 10L,
                                                 base_seed = 20260524L,
                                                 distance_type = "geodesic",
                                                 control = list(
                                                   cardioid_optim_control = list(maxit = 1000)
                                                 )) {
  if (is.null(output_root)) {
    output_root <- file.path(
      "output",
      "comets_distance_profile_cardioid",
      paste0("run_", timestamp_tag_cardioid_comets())
    )
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  comets_data <- load_comets_distance_profile_data()
  models <- make_comet_cardioid_models(model_ids = model_ids)
  dataset_summary <- data.frame(
    dataset = c("oort", "short_period"),
    n = c(nrow(comets_data$oort), nrow(comets_data$short)),
    ambient_dim = c(ncol(comets_data$oort$normal), ncol(comets_data$short$normal)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    dataset_summary,
    file = file.path(output_root, "dataset_summary.csv"),
    row.names = FALSE
  )

  write_lines_if_possible_cardioid(
    capture.output(sessionInfo()),
    file.path(output_root, "sessionInfo.txt")
  )

  saveRDS(
    list(
      output_root = output_root,
      stages = as.character(stages),
      model_ids = vapply(models, `[[`, character(1), "id"),
      cvm_B = as.integer(cvm_B),
      ks_B = as.integer(ks_B),
      ks_omega_points = as.integer(ks_omega_points),
      ks_t_points = as.integer(ks_t_points),
      n_cores = as.integer(n_cores),
      base_seed = as.integer(base_seed),
      distance_type = distance_type
    ),
    file = file.path(output_root, "run_config.rds")
  )

  ks_grid <- list(
    omega_grid = generate_canonical_lattice(as.integer(ks_omega_points), dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = as.integer(ks_t_points))
  )

  stage_definitions <- list(
    oort_cvm = list(
      stage_id = "01",
      label = "Oort cardioid CvM",
      dataset = comets_data$oort$normal,
      statistic = "cvm",
      B = cvm_B,
      ks_grid = NULL
    ),
    oort_ks = list(
      stage_id = "04",
      label = "Oort cardioid KS",
      dataset = comets_data$oort$normal,
      statistic = "ks",
      B = ks_B,
      ks_grid = ks_grid
    ),
    short_cvm = list(
      stage_id = "02",
      label = "Short-period cardioid CvM",
      dataset = comets_data$short$normal,
      statistic = "cvm",
      B = cvm_B,
      ks_grid = NULL
    ),
    short_ks = list(
      stage_id = "03",
      label = "Short-period cardioid KS",
      dataset = comets_data$short$normal,
      statistic = "ks",
      B = ks_B,
      ks_grid = ks_grid
    )
  )

  requested_stages <- as.character(stages)
  if (!all(requested_stages %in% names(stage_definitions))) {
    stop("`stages` contains unsupported stage identifiers.")
  }

  manifest_rows <- list()
  stage_results <- list()

  for (stage_name in requested_stages) {
    stage <- stage_definitions[[stage_name]]
    stage_dir <- file.path(output_root, paste0(stage$stage_id, "_", stage_name))
    stage_start <- Sys.time()

    stage_results[[stage_name]] <- run_cardioid_comet_stage(
      data_matrix = stage$dataset,
      dataset_label = stage$label,
      statistic = stage$statistic,
      stage_dir = stage_dir,
      models = models,
      B = stage$B,
      n_cores = n_cores,
      seed = base_seed + as.integer(stage$stage_id),
      ks_grid = stage$ks_grid,
      distance_type = distance_type,
      control = control
    )

    manifest_rows <- append_manifest_row_cardioid(
      manifest_rows = manifest_rows,
      stage_id = stage$stage_id,
      stage_label = stage$label,
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
  }

  pipeline_result <- list(
    output_root = output_root,
    dataset_summary = dataset_summary,
    stages = stage_results,
    manifest = do.call(rbind, manifest_rows),
    engine = "multiplier_bootstrap_gof",
    method = "distance_profiles"
  )

  saveRDS(pipeline_result, file = file.path(output_root, "pipeline_result.rds"))
  pipeline_result
}
