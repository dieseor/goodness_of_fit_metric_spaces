## B=200 diagnostic run of the production re-estimated multiplier bootstrap.

source(file.path("utils.R"))
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "uniform_beta_mixture_model_spec.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("real_data", "comets", "utils_comets_data.R"))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) substring(hit[[1L]], nchar(name) + 4L) else default
}
B <- as.integer(arg("B", "200"))
n_cores <- as.integer(arg("n_cores", "3"))
output_root <- arg("output_root", file.path(
  "real_data", "comets", "diagnostics", "uniform_beta_bootstrap_gof_B200_20260915"
))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

fit_root <- file.path(
  "real_data", "reruns", "paper_main_realdata_B1000_3cores_20260831_113532",
  "comets", "uniform_beta"
)
read_theta <- function(path) {
  z <- utils::read.csv(path, stringsAsFactors = FALSE)[1L, ]
  list(mu = as.numeric(z[c("mu_1", "mu_2", "mu_3")]),
       weight_uniform = z$weight_uniform, alpha = z$alpha, beta = z$beta)
}
nearest <- function(mu, x) {
  z <- pmin(pmax(as.numeric(x %*% mu), -1), 1)
  c(min(acos(z)), min(acos(-z)))
}
extract_inference <- function(result, statistic) {
  z <- result$inference
  if (is.list(z) && !is.null(z[[statistic]])) z <- z[[statistic]]
  get <- function(names) {
    for (name in names) if (!is.null(z[[name]])) return(as.numeric(z[[name]][[1L]]))
    NA_real_
  }
  c(observed = get(c("observed", "observed_statistic", "statistic", "Tn")),
    p_value = get(c("p_value", "p.value", "pvalue", "p")),
    critical_value = get(c("critical_value", "critical", "quantile")))
}

comets <- load_comets_real_data(finite_normals = "both")
samples <- list(short = as.matrix(comets$short$normal), long = as.matrix(comets$long$normal))
fit_dirs <- c(short = "01_short_period_uniform_beta_mixture",
              long = "02_long_period_uniform_beta_mixture")
base_case_seed <- c(short = 20260602L, long = 20260603L)
control <- list(
  uniform_beta_mixture_profile_method = "legendre",
  uniform_beta_mixture_quad_n = 100L,
  uniform_beta_mixture_optim_control = list(maxit = 350L, reltol = 1e-9),
  uniform_beta_mixture_bootstrap_n_starts = 1L,
  uniform_beta_mixture_bootstrap_optim_control = list(maxit = 80L, reltol = 1e-6),
  progress_bar = TRUE,
  reestimated_bootstrap_chunk_size = max(1L, ceiling(B / (8L * n_cores)))
)
spec <- make_uniform_beta_mixture_spec(distance_type = "geodesic")
summaries <- list()

for (dataset in names(samples)) {
  x <- samples[[dataset]]
  theta <- read_theta(file.path(fit_root, fit_dirs[[dataset]], "theta_hat.csv"))
  for (j in seq_along(c("ks", "cvm"))) {
    statistic <- c("ks", "cvm")[[j]]
    seed <- base_case_seed[[dataset]] + 1000L * j
    control$progress_label <- paste(dataset, statistic, "UB audit")
    message("[bootstrap GOF] ", dataset, " / ", statistic, " / B=", B, " / seed=", seed)
    result <- multiplier_bootstrap_gof(
      data = x, spec = spec, null = list(type = "composite"),
      statistics = statistic,
      ks_grid = if (statistic == "ks") make_sample_unique_distance_ks_grid() else NULL,
      B = B, alpha = 0.05, n_cores = n_cores, seed = seed,
      bootstrap_method = "reestimated",
      keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = TRUE),
      control = control, observed_theta_hat = theta
    )
    prefix <- paste(dataset, statistic, sep = "_")
    saveRDS(result, file.path(output_root, paste0(prefix, ".rds")))
    theta_rows <- do.call(rbind, lapply(seq_along(result$bootstrap$theta_star), function(b) {
      th <- result$bootstrap$theta_star[[b]]
      nd <- nearest(th$mu, x)
      data.frame(
        dataset = dataset, statistic = statistic, replicate = b,
        mu_1 = th$mu[[1L]], mu_2 = th$mu[[2L]], mu_3 = th$mu[[3L]],
        weight_uniform = th$weight_uniform, alpha = th$alpha, beta = th$beta,
        distance_to_nearest = nd[[1L]], antipode_distance_to_nearest = nd[[2L]],
        bootstrap_statistic = result$bootstrap$statistics[[statistic]][[b]]
      )
    }))
    utils::write.csv(theta_rows, file.path(output_root, paste0(prefix, "_replicates.csv")), row.names = FALSE)
    inf <- extract_inference(result, statistic)
    summaries[[length(summaries) + 1L]] <- data.frame(
      dataset = dataset, statistic = statistic, B = B, seed = seed,
      observed = inf[["observed"]], p_value = inf[["p_value"]],
      critical_value = inf[["critical_value"]],
      bootstrap_min = min(theta_rows$bootstrap_statistic),
      bootstrap_median = median(theta_rows$bootstrap_statistic),
      bootstrap_max = max(theta_rows$bootstrap_statistic),
      singular_count = sum(pmin(theta_rows$distance_to_nearest,
                                theta_rows$antipode_distance_to_nearest) <= 2.1e-6),
      elapsed_seconds = result$diagnostics$elapsed_seconds,
      effective_bootstrap_method = result$diagnostics$effective_bootstrap_method
    )
    utils::write.csv(do.call(rbind, summaries), file.path(output_root, "summary.csv"), row.names = FALSE)
  }
}

writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
writeLines(c(
  paste("B:", B), paste("n_cores:", n_cores),
  "multipliers: Exp(1), normalized by replicate row mean",
  "method: production reestimated multiplier bootstrap",
  "fit: saved paper theta; warm start only; maxit=80; reltol=1e-6",
  "KS grid: sample points and unique distances",
  "seeds: paper case seed + 1000 for KS and +2000 for CvM"
), file.path(output_root, "manifest.txt"))
message("Bootstrap GOF audit artifacts: ", output_root)
