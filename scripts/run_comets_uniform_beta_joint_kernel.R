#!/usr/bin/env Rscript

# Standalone rerun of the comet UB experiment with a joint re-estimated
# KS/CvM kernel.  It leaves the paper runner and its saved outputs untouched.
source(file.path("scripts", "run_comets_rotational_mixtures_short_long.R"))

args <- parse_named_args_rotmix_comets(commandArgs(trailingOnly = TRUE))
B <- as.integer(args$B %||% 1000L)
n_cores <- min(8L, as.integer(args$n_cores %||% 8L))
output_root <- args$output_root %||% file.path(
  "real_data", "reruns", sprintf("comets_uniform_beta_joint_kernel_B%d", B)
)
if (!is.finite(B) || B < 1L || !is.finite(n_cores) || n_cores < 1L) stop("`B` and `n_cores` must be positive.")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

control <- list(
  uniform_beta_mixture_profile_method = "legendre",
  uniform_beta_mixture_quad_n = 100L,
  uniform_beta_mixture_optim_control = list(maxit = 350L, reltol = 1e-9),
  reestimated_fuse_ks_cvm = TRUE,
  progress_bar = TRUE
)
comets <- load_comets_distance_profile_data_rotmix()
datasets <- list(short = as.matrix(comets$short$normal), long = as.matrix(comets$long$normal))
# Historical per-dataset base seeds; the joint run uses the former KS stream so
# it is directly comparable to the saved KS experiment while sharing it with CvM.
joint_seeds <- c(short = 20261711L, long = 20261712L)
spec <- rotmix_make_spec("uniform_beta_mixture", distance_type = "geodesic")
rows <- list()

for (dataset in names(datasets)) {
  x <- datasets[[dataset]]; label <- paste0(dataset, "_period")
  out <- file.path(output_root, sprintf("%s_uniform_beta_mixture", label))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  theta_hat <- rotmix_fit_model_theta("uniform_beta_mixture", x, control = control)
  result <- multiplier_bootstrap_gof(
    data = x, spec = spec, null = list(type = "composite"), statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = .05, n_cores = n_cores,
    seed = joint_seeds[[dataset]], bootstrap_method = "reestimated",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
    control = control, observed_theta_hat = theta_hat
  )
  utils::write.csv(rotmix_theta_export_row("uniform_beta_mixture", theta_hat), file.path(out, "theta_hat.csv"), row.names = FALSE)
  saveRDS(result, file.path(out, "gof_joint_ks_cvm.rds"))
  inference <- result$inference
  rows[[dataset]] <- data.frame(
    dataset = label, model = "uniform_beta_mixture", n = nrow(x), B = B, n_cores = n_cores,
    seed = joint_seeds[[dataset]], gof_ks_observed = inference$ks$observed, gof_ks_p_value = inference$ks$p_value,
    gof_cvm_observed = inference$cvm$observed, gof_cvm_p_value = inference$cvm$p_value,
    bootstrap_method = result$diagnostics$effective_bootstrap_method,
    shared_sample_ks_cvm_cache = result$diagnostics$shared_sample_ks_cvm_cache,
    reestimated_fuse_ks_cvm_effective = result$diagnostics$reestimated_fuse_ks_cvm_effective,
    elapsed_seconds = result$diagnostics$elapsed_seconds,
    mu_1 = theta_hat$mu[[1L]], mu_2 = theta_hat$mu[[2L]], mu_3 = theta_hat$mu[[3L]],
    weight_uniform = theta_hat$weight_uniform, alpha = theta_hat$alpha, beta = theta_hat$beta
  )
  utils::write.csv(rows[[dataset]], file.path(out, "summary.csv"), row.names = FALSE)
}
utils::write.csv(do.call(rbind, rows), file.path(output_root, "comets_uniform_beta_joint_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
writeLines(c(
  "Standalone joint-kernel rerun; official paper outputs were not changed.",
  "Observed fit: current production UB MLE with the paper's control values.",
  "Bootstrap controls are inherited unchanged from the paper runner; only KS/CvM are now requested jointly.",
  "KS and CvM are requested together, sharing the full sample-based profile kernel and theta_star per replicate.",
  "The joint seed for each dataset is the saved paper KS seed (short 20261711; long 20261712)."
), file.path(output_root, "manifest.txt"))
message("Joint-kernel comet UB rerun written to: ", normalizePath(output_root))
