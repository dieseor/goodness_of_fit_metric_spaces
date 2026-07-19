resolve_comets_watson_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_comets_watson_path("bootstrap", "model_specs.R"))
source(resolve_comets_watson_path("bootstrap", "multiplier_bootstrap.R"))
source(resolve_comets_watson_path("real_data", "comets", "utils_comets_data.R"))
source(resolve_comets_watson_path("scripts", "path_helpers.R"))

summarize_watson_comet_result <- function(result, dataset_label, statistic) {
  theta <- result$observed$theta_hat
  inference <- result$inference[[statistic]]
  data.frame(
    dataset = dataset_label,
    model = "watson_composite",
    null_type = "composite",
    statistic = statistic,
    spec_name = result$diagnostics$spec_name,
    bootstrap_method = result$diagnostics$bootstrap_method,
    effective_bootstrap_method = result$diagnostics$effective_bootstrap_method,
    observed_statistic = inference$observed,
    critical_value = inference$critical_value,
    p_value = inference$p_value,
    reject = inference$reject,
    mu_hat_1 = theta$mu[[1L]],
    mu_hat_2 = theta$mu[[2L]],
    mu_hat_3 = theta$mu[[3L]],
    kappa_hat = theta$kappa,
    B = result$bootstrap$B,
    n = result$diagnostics$n,
    n_cores = result$diagnostics$n_cores,
    elapsed_seconds = result$diagnostics$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

run_comets_distance_profile_watson_benchmark <- function(output_root = NULL,
                                                          dataset = c("short", "long"),
                                                          B = 5000L,
                                                          statistic = c("ks", "cvm"),
                                                          n_cores = 8L,
                                                          seed = 20260714L,
                                                          distance_type = "geodesic",
                                                          control = list(
                                                            watson_L_max = 200L,
                                                            watson_quad_n = 400L,
                                                            watson_tol = 1e-10,
                                                            derivative_mc_size = 1000L,
                                                            derivative_mc_seed = 20260714L
                                                          )) {
  dataset <- match.arg(dataset)
  statistic <- match.arg(statistic)
  B <- as.integer(B)
  if (!is.finite(B) || B <= 0L) stop("`B` must be a strictly positive integer.")
  if (is.null(output_root)) {
    output_root <- canonical_comets_watson_dir(sprintf("%s_comets_B%d_sample%s", dataset, B, statistic), "fast")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  comets <- load_comets_real_data(finite_normals = "none")
  x <- if (identical(dataset, "short")) comets$short$normal else comets$long$normal
  dataset_label <- if (identical(dataset, "short")) "short_period" else "long_period"
  ks_grid <- if (identical(statistic, "ks")) make_sample_unique_distance_ks_grid() else NULL
  result <- multiplier_bootstrap_watson(
    data = x,
    null = list(type = "composite"),
    statistics = statistic,
    ks_grid = ks_grid,
    B = B,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = TRUE, bootstrap_statistics = TRUE, bootstrap_thetas = TRUE),
    control = utils::modifyList(control, list(derivative_mc_seed = control$derivative_mc_seed %||% as.integer(seed))),
    distance_type = distance_type
  )
  summary <- summarize_watson_comet_result(result, dataset_label, statistic)
  utils::write.csv(summary, file.path(output_root, "benchmark_summary.csv"), row.names = FALSE)
  write_comets_dataset_summary(file.path(output_root, "dataset_summary.csv"), dataset_label, x)
  writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
  saveRDS(result, file.path(output_root, sprintf("stage_01_M%d_B%d.rds", B, B)))
  invisible(summary)
}
