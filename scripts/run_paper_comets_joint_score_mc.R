#!/usr/bin/env Rscript

# Joint KS/CvM rerun for the C2 and small-circle fits in the paper.
# Run from the repository root. Historical results and fits are preserved.
Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("bootstrap/multiplier_bootstrap.R")
source("real_data/comets/utils_comets_data.R")

paper_comets_joint_control <- function(derivative_seed) {
  list(
    derivative_method = "score_mc",
    derivative_mc_size = 10000L,
    derivative_mc_seed = as.integer(derivative_seed),
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fast_multiplier_fuse_ks_cvm = TRUE,
    fast_multiplier_cache_corrections = "auto",
    fast_multiplier_cvm_block_size = 50L,
    small_circle_profile_method = "legendre",
    small_circle_L_max = 200L,
    small_circle_quad_n = 400L,
    small_circle_tol = 1e-10
  )
}

run_paper_comets_joint_case <- function(data, model, theta_hat, B = 1000L,
                                      n_cores = 3L, seed,
                                      derivative_seed = seed + 50000L) {
  model <- match.arg(model, c("C2", "SC"))
  spec <- if (model == "C2") {
    make_cardioid_spec(k = 2L, distance_type = "geodesic", unknown_param = "both")
  } else {
    make_small_circle_spec(distance_type = "geodesic")
  }
  result <- multiplier_bootstrap_gof(
    data = data, spec = spec, null = list(type = "composite"),
    observed_theta_hat = theta_hat,
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = as.integer(B), alpha = 0.05, n_cores = as.integer(n_cores),
    seed = as.integer(seed), bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE,
                bootstrap_thetas = FALSE),
    control = paper_comets_joint_control(derivative_seed)
  )
  d <- result$diagnostics
  checks <- c(
    fast = identical(d$effective_bootstrap_method, "fast_multiplier"),
    no_fallback = identical(d$fallback_to_reestimated, FALSE),
    cpp = identical(d$fast_multiplier_backend_effective, "cpp"),
    kernel = identical(d$fast_multiplier_cpp_kernel_effective, "contiguous_double"),
    fused = isTRUE(d$fast_multiplier_fuse_ks_cvm_effective),
    shared = isTRUE(d$shared_sample_ks_cvm_cache),
    score_mc = identical(d$derivative_method_effective, "score_mc"),
    mc_size = identical(as.integer(d$derivative_mc_size), 10000L),
    cores = identical(as.integer(d$n_cores), as.integer(n_cores))
  )
  if (!all(checks)) {
    stop("Joint comet configuration was not effective: ",
         paste(names(checks)[!checks], collapse = ", "))
  }
  result
}

run_paper_comets_joint <- function(output_dir, B = 1000L, n_cores = 3L,
    reference_dir = "real_data/reruns/paper_main_realdata_B1000_3cores_20260831_113532") {
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    stop("Use a new output directory to preserve earlier results.")
  }
  if (B < 1L || n_cores < 1L) stop("B and n_cores must be positive.")
  design <- data.frame(
    model = c("C2", "C2", "SC", "SC"),
    period = c("short", "long", "short", "long"),
    seed = c(20260530L, 20260531L, 20260532L, 20260532L),
    fit_source = file.path(reference_dir, c(
      "comets/cardioid_C2/03_short_ks/stage_bundle.rds",
      "comets/cardioid_C2/04_oort_ks/stage_bundle.rds",
      "comets/small_circle/short_ks/stage_01_M1000_B1000.rds",
      "comets/small_circle/long_ks/stage_01_M1000_B1000.rds"
    )), stringsAsFactors = FALSE
  )
  if (!all(file.exists(design$fit_source))) stop("A reference fit is missing.")
  data <- load_comets_real_data(finite_normals = "both")
  stopifnot(nrow(data$short$normal) == 784L, nrow(data$long$normal) == 610L)
  design$derivative_seed <- design$seed + 50000L
  design$B <- as.integer(B)
  design$n_cores <- as.integer(n_cores)
  design$derivative_mc_size <- 10000L
  design$fit_source_md5 <- unname(tools::md5sum(design$fit_source))
  dir.create(output_dir, recursive = TRUE)
  utils::write.csv(design, file.path(output_dir, "manifest.csv"), row.names = FALSE)
  saveRDS(lapply(design$derivative_seed, paper_comets_joint_control),
          file.path(output_dir, "controls.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  summaries <- list()
  for (i in seq_len(nrow(design))) {
    job <- design[i, ]
    reference <- readRDS(job$fit_source)
    if (job$model == "C2") reference <- reference$results$C2
    cat(sprintf("Running %s %s: joint KS/CvM, B=%d, Nderiv=10000, cores=%d\n",
                job$model, job$period, B, n_cores))
    result <- run_paper_comets_joint_case(
      data = data[[job$period]]$normal, model = job$model,
      theta_hat = reference$observed$theta_hat, B = B, n_cores = n_cores,
      seed = job$seed, derivative_seed = job$derivative_seed
    )
    saveRDS(result, file.path(output_dir, paste0(job$model, "_", job$period, ".rds")))
    summaries[[i]] <- data.frame(
      model = job$model, period = job$period, n = result$diagnostics$n,
      B = B, derivative_mc_size = 10000L, n_cores = n_cores,
      ks_pvalue = result$inference$ks$p_value,
      cvm_pvalue = result$inference$cvm$p_value,
      ks_reject = result$inference$ks$reject,
      cvm_reject = result$inference$cvm$reject,
      backend = result$diagnostics$fast_multiplier_backend_effective,
      kernel = result$diagnostics$fast_multiplier_cpp_kernel_effective,
      fused = result$diagnostics$fast_multiplier_fuse_ks_cvm_effective,
      shared = result$diagnostics$shared_sample_ks_cvm_cache
    )
    utils::write.csv(do.call(rbind, summaries), file.path(output_dir, "summary.csv"),
                     row.names = FALSE)
  }
  invisible(do.call(rbind, summaries))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  cli <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) stop("Use --name=value arguments.")
    pair <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    cli[[pair[1L]]] <- paste(pair[-1L], collapse = "=")
  }
  unknown <- setdiff(names(cli), c("output_dir", "B", "n_cores", "reference_dir"))
  if (length(unknown)) stop("Unknown arguments: ", paste(unknown, collapse = ", "))
  requested_B <- as.integer(cli$B %||% 1000L)
  run_paper_comets_joint(
    output_dir = cli$output_dir %||% file.path("real_data/reruns",
      paste0("paper_comets_joint_B", requested_B, "_Nderiv10000_",
             format(Sys.time(), "%Y%m%d_%H%M%S"))),
    B = requested_B, n_cores = as.integer(cli$n_cores %||% 3L),
    reference_dir = cli$reference_dir %||%
      "real_data/reruns/paper_main_realdata_B1000_3cores_20260831_113532"
  )
}
