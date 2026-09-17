#!/usr/bin/env Rscript

# Diagnostic only: evaluates selected pre-existing UB fits with the production
# re-estimated multiplier bootstrap.  It does not alter production code/results.
source("utils.R")
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "uniform_beta_mixture_model_spec.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("real_data", "comets", "utils_comets_data.R"))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
B <- as.integer(value("--B", "200")); n_cores <- min(8L, as.integer(value("--cores", "8")))
output_root <- value("--output", file.path("real_data", "comets", "diagnostics", "uniform_beta_gof_basin_sensitivity_B200_20260916"))
requested_cases <- strsplit(value("--cases", ""), ",", fixed = TRUE)[[1L]]
requested_cases <- trimws(requested_cases[nzchar(trimws(requested_cases))])
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_theta <- function(path) {
  z <- utils::read.csv(path, stringsAsFactors = FALSE)[1L, ]
  list(mu = as.numeric(z[c("mu_1", "mu_2", "mu_3")]), weight_uniform = z$weight_uniform, alpha = z$alpha, beta = z$beta)
}
nearest_and_clipping <- function(theta, x) {
  dot <- pmin(1, pmax(-1, as.numeric(x %*% theta$mu))); y <- (dot + 1) / 2
  eps <- 1e-12
  endpoint_y_clip <- any(y <= eps * (1 + 1e-6) | y >= 1 - eps * (1 + 1e-6))
  parameter_clip <- theta$weight_uniform <= .01 * (1 + 1e-8) || theta$weight_uniform >= .99 * (1 - 1e-8) ||
    theta$alpha <= .05 * (1 + 1e-8) || theta$alpha >= 1000 * (1 - 1e-8) ||
    theta$beta <= .05 * (1 + 1e-8) || theta$beta >= 1000 * (1 - 1e-8)
  c(nearest_mu = min(acos(dot)), nearest_antimu = min(acos(-dot)), endpoint_y_clip = endpoint_y_clip,
    parameter_clip = parameter_clip, clipping_active = endpoint_y_clip || parameter_clip)
}
extract_inference <- function(result, statistic) {
  z <- result$inference; if (is.list(z) && !is.null(z[[statistic]])) z <- z[[statistic]]
  get <- function(keys) for (key in keys) if (!is.null(z[[key]])) return(as.numeric(z[[key]][[1L]]))
  c(observed = get(c("observed", "observed_statistic", "statistic", "Tn")),
    p_value = get(c("p_value", "p.value", "pvalue", "p")),
    critical_value = get(c("critical_value", "critical", "quantile")))
}

rerun_root <- file.path("real_data", "reruns", "paper_main_realdata_B1000_3cores_20260831_113532", "comets", "uniform_beta")
audit_root <- file.path("real_data", "comets", "diagnostics", "uniform_beta_mle_repeated_1000_20260916")
reported_short <- read_theta(file.path(rerun_root, "01_short_period_uniform_beta_mixture", "theta_hat.csv"))
reported_long <- read_theta(file.path(rerun_root, "02_long_period_uniform_beta_mixture", "theta_hat.csv"))
representatives <- utils::read.csv(file.path(audit_root, "long_representative_density_basins.csv"), stringsAsFactors = FALSE)
pick_medoid <- function(basin) representatives[representatives$role == "medoid" & representatives$density_basin == basin, ][1L, ]
theta_from_rep <- function(z) list(mu = as.numeric(z[c("mu1_canonical", "mu2_canonical", "mu3_canonical")]),
  weight_uniform = z$weight_uniform, alpha = z$alpha_canonical, beta = z$beta_canonical)

# Basins 2 and 4 are genuinely distinct from the reported fit and their medoids
# are neither clipped nor within 1e-5 radians of a sample point or antipode.
cases <- list(
  short_reported = list(dataset = "short", basin = "reported", theta = reported_short),
  long_reported = list(dataset = "long", basin = "reported", theta = reported_long),
  long_basin2 = list(dataset = "long", basin = "density_basin_2_medoid", theta = theta_from_rep(pick_medoid(2L))),
  long_basin4 = list(dataset = "long", basin = "density_basin_4_medoid", theta = theta_from_rep(pick_medoid(4L)))
)
if (length(requested_cases)) {
  unknown_cases <- setdiff(requested_cases, names(cases))
  if (length(unknown_cases)) stop("Unknown case(s): ", paste(unknown_cases, collapse = ", "))
  cases <- cases[requested_cases]
}
comets <- load_comets_real_data(finite_normals = "both")
samples <- list(short = as.matrix(comets$short$normal), long = as.matrix(comets$long$normal))
control <- list(
  uniform_beta_mixture_profile_method = "legendre", uniform_beta_mixture_quad_n = 100L,
  uniform_beta_mixture_optim_control = list(maxit = 350L, reltol = 1e-9),
  uniform_beta_mixture_bootstrap_n_starts = 1L,
  uniform_beta_mixture_bootstrap_optim_control = list(maxit = 80L, reltol = 1e-6),
  progress_bar = TRUE, reestimated_bootstrap_chunk_size = max(1L, ceiling(B / (8L * n_cores)))
)
spec <- make_uniform_beta_mixture_spec(distance_type = "geodesic")
summary_rows <- list(); selected_rows <- list()
for (case_name in names(cases)) {
  case <- cases[[case_name]]; x <- samples[[case$dataset]]; base_seed <- if (case$dataset == "short") 20260602L else 20260603L
  nd <- nearest_and_clipping(case$theta, x)
  selected_rows[[case_name]] <- data.frame(case = case_name, dataset = case$dataset, basin = case$basin,
    mu_1 = case$theta$mu[[1L]], mu_2 = case$theta$mu[[2L]], mu_3 = case$theta$mu[[3L]],
    weight_uniform = case$theta$weight_uniform, alpha = case$theta$alpha, beta = case$theta$beta,
    t(nd), stringsAsFactors = FALSE)
  seed <- base_seed + 1000L
  control$progress_label <- sprintf("UB joint basin diagnostic: %s", case_name)
  message("[", case_name, "] KS+CvM jointly, B=", B, ", seed=", seed)
  result <- multiplier_bootstrap_gof(data = x, spec = spec, null = list(type = "composite"), statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = .05, n_cores = n_cores, seed = seed,
    bootstrap_method = "reestimated", keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = TRUE),
    control = control, observed_theta_hat = case$theta)
  saveRDS(result, file.path(output_root, sprintf("%s_joint.rds", case_name)))
  reps <- do.call(rbind, lapply(seq_along(result$bootstrap$theta_star), function(b) {
    th <- result$bootstrap$theta_star[[b]]; q <- nearest_and_clipping(th, x)
    data.frame(case = case_name, dataset = case$dataset, basin = case$basin, replicate = b,
      mu_1 = th$mu[[1L]], mu_2 = th$mu[[2L]], mu_3 = th$mu[[3L]], weight_uniform = th$weight_uniform,
      alpha = th$alpha, beta = th$beta, t(q), bootstrap_ks = result$bootstrap$statistics$ks[[b]],
      bootstrap_cvm = result$bootstrap$statistics$cvm[[b]])
  }))
  utils::write.csv(reps, file.path(output_root, sprintf("%s_joint_replicates.csv", case_name)), row.names = FALSE)
  for (statistic in c("ks", "cvm")) {
    stats_star <- reps[[paste0("bootstrap_", statistic)]]
    inf <- extract_inference(result, statistic); keep <- !reps$clipping_active
    p_sensitivity <- (1 + sum(stats_star[keep] >= inf[["observed"]])) / (1 + sum(keep))
    endpoint <- stats_star[!keep]; interior <- stats_star[keep]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      case = case_name, dataset = case$dataset, basin = case$basin, statistic = statistic, B = B, seed = seed,
      observed = inf[["observed"]], p_value_all = inf[["p_value"]], p_value_excluding_clipped = p_sensitivity,
      critical_value_all = inf[["critical_value"]], n_clipped = sum(!keep), n_interior = sum(keep),
      endpoint_tail_count = sum(endpoint >= inf[["observed"]]), interior_tail_count = sum(interior >= inf[["observed"]]),
      endpoint_median = if (length(endpoint)) median(endpoint) else NA_real_, endpoint_max = if (length(endpoint)) max(endpoint) else NA_real_,
      interior_median = if (length(interior)) median(interior) else NA_real_, interior_max = if (length(interior)) max(interior) else NA_real_,
      upper_tail_from_endpoint = if (sum(stats_star >= inf[["observed"]])) sum(endpoint >= inf[["observed"]]) / sum(stats_star >= inf[["observed"]]) else NA_real_,
      elapsed_seconds = result$diagnostics$elapsed_seconds, effective_bootstrap_method = result$diagnostics$effective_bootstrap_method,
      stringsAsFactors = FALSE)
    utils::write.csv(do.call(rbind, summary_rows), file.path(output_root, "summary.csv"), row.names = FALSE)
  }
}
utils::write.csv(do.call(rbind, selected_rows), file.path(output_root, "selected_fits.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
writeLines(c(
  sprintf("B=%d; n_cores=%d", B, n_cores),
  sprintf("cases=%s", paste(names(cases), collapse = ",")),
  "Each call jointly requests KS and CvM, sharing multipliers, reestimated theta_star, and the full sample-based profile kernel.",
  "Multipliers are normalized Exp(1); each bootstrap refit has one warm start, maxit=80, reltol=1e-6.",
  "KS uses make_sample_unique_distance_ks_grid().",
  "Clipped means the production endpoint-y or parameter clip is active. The sensitivity p-value removes only these rows and uses the same +1 correction on the retained B."
), file.path(output_root, "manifest.txt"))
message("Wrote UB GOF basin sensitivity diagnostic to: ", normalizePath(output_root))
