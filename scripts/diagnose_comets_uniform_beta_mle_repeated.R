#!/usr/bin/env Rscript

# Diagnostic only: repeatedly calls the production UB MLE routine with one
# externally supplied warm start.  It does not modify the production pipeline.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) default else args[[hit + 1L]]
}
n_fits <- as.integer(arg_value("--n-fits", 100L))
n_cores <- min(8L, as.integer(arg_value("--cores", 8L)))
seed <- as.integer(arg_value("--seed", 20260916L))
output_root <- arg_value("--output", NULL)
if (!is.finite(n_fits) || n_fits < 1L) stop("`--n-fits` must be positive.")
if (!is.finite(n_cores) || n_cores < 1L) stop("`--cores` must be positive.")

source("utils.R")
source(file.path("real_data", "comets", "utils_comets_data.R"))

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (is.null(output_root)) {
  output_root <- file.path("real_data", "comets", "diagnostics",
    sprintf("uniform_beta_mle_repeated_%dfits_%s", n_fits, tag))
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# These are exactly the controls passed by the comet application.  In
# particular, no bounds, clipping, method, or tolerance is supplied here
# beyond the application control list.
production_control <- list(
  uniform_beta_mixture_profile_method = "legendre",
  uniform_beta_mixture_quad_n = 100L,
  uniform_beta_mixture_optim_control = list(maxit = 350L, reltol = 1e-9)
)

read_saved_theta <- function(which) {
  root <- file.path("real_data", "reruns",
    "paper_main_realdata_B1000_3cores_20260831_113532", "comets",
    "uniform_beta")
  dir_name <- if (identical(which, "short")) {
    "01_short_period_uniform_beta_mixture"
  } else "02_long_period_uniform_beta_mixture"
  z <- utils::read.csv(file.path(root, dir_name, "theta_hat.csv"), check.names = FALSE)
  list(mu = as.numeric(z[1L, c("mu_1", "mu_2", "mu_3")]),
       weight_uniform = z$weight_uniform[[1L]], alpha = z$alpha[[1L]], beta = z$beta[[1L]])
}

geodesic_nearest <- function(mu, x) {
  acos(pmin(1, pmax(-1, as.numeric(x %*% mu))))
}

make_start_table <- function(x, theta_saved, n, seed) {
  set.seed(seed)
  m <- nrow(x)
  dirs <- matrix(NA_real_, nrow = n, ncol = 3L)
  # Exact saved estimate first; the remainder deliberately cover general,
  # data-centred, and antipodal directions without changing the MLE routine.
  dirs[1L, ] <- theta_saved$mu
  if (n >= 2L) {
    for (i in seq.int(2L, n)) {
      mode <- (i - 2L) %% 10L
      if (mode <= 4L) {
        u <- stats::rnorm(3L); dirs[i, ] <- u / sqrt(sum(u^2))
      } else if (mode <= 7L) {
        dirs[i, ] <- x[sample.int(m, 1L), ]
      } else {
        dirs[i, ] <- -x[sample.int(m, 1L), ]
      }
    }
  }
  w_grid <- c(0.02, 0.05, 0.12, 0.30, 0.55, 0.80, 0.95)
  shape_grid <- c(0.08, 0.20, 0.50, 0.75, 1.20, 2, 5, 20, 40)
  out <- data.frame(
    start_id = seq_len(n),
    start_mu1 = dirs[, 1L], start_mu2 = dirs[, 2L], start_mu3 = dirs[, 3L],
    start_weight_uniform = w_grid[(seq_len(n) - 1L) %% length(w_grid) + 1L],
    start_alpha = shape_grid[(seq_len(n) - 1L) %% length(shape_grid) + 1L],
    start_beta = shape_grid[(3L * seq_len(n) - 1L) %% length(shape_grid) + 1L]
  )
  out[1L, c("start_weight_uniform", "start_alpha", "start_beta")] <-
    unlist(theta_saved[c("weight_uniform", "alpha", "beta")])
  out
}

one_fit <- function(row, x, dataset_name) {
  start_theta <- list(
    mu = as.numeric(row[c("start_mu1", "start_mu2", "start_mu3")]),
    weight_uniform = row$start_weight_uniform,
    alpha = row$start_alpha,
    beta = row$start_beta
  )
  control <- production_control
  control$uniform_beta_mixture_start_theta <- start_theta
  control$uniform_beta_mixture_warm_start_only <- TRUE
  warnings <- character()
  elapsed <- system.time({
    fit <- withCallingHandlers(
      try(uniform_beta_mixture_mle_s2_weighted(x = x, control = control), silent = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })[["elapsed"]]
  if (inherits(fit, "try-error")) {
    return(cbind(row, data.frame(dataset = dataset_name, error = as.character(fit),
      elapsed_seconds = elapsed, stringsAsFactors = FALSE)))
  }
  theta <- fit[c("mu", "weight_uniform", "alpha", "beta")]
  y_raw <- (pmin(1, pmax(-1, as.numeric(x %*% theta$mu))) + 1) / 2
  eps <- 1e-12
  loglik <- sum(d_sph_uniform_beta_mixture_s2(
    x, theta$mu, theta$weight_uniform, theta$alpha, theta$beta, log = TRUE, eps = eps
  ))
  nearest_mu <- geodesic_nearest(theta$mu, x)
  nearest_antimu <- geodesic_nearest(-theta$mu, x)
  endpoint_y_clip <- any(y_raw <= eps * (1 + 1e-6) | y_raw >= 1 - eps * (1 + 1e-6))
  parameter_clip <- theta$weight_uniform <= 0.01 * (1 + 1e-8) ||
    theta$weight_uniform >= 0.99 * (1 - 1e-8) || theta$alpha <= 0.05 * (1 + 1e-8) ||
    theta$alpha >= 1000 * (1 - 1e-8) || theta$beta <= 0.05 * (1 + 1e-8) ||
    theta$beta >= 1000 * (1 - 1e-8)
  counts <- fit$opt$counts
  returned_start <- fit$start_theta
  changed_start <- max(abs(c(returned_start$mu, returned_start$weight_uniform,
    returned_start$alpha, returned_start$beta) - c(start_theta$mu, start_theta$weight_uniform,
    start_theta$alpha, start_theta$beta))) > 1e-10
  cbind(row, data.frame(
    dataset = dataset_name,
    mu1 = theta$mu[[1L]], mu2 = theta$mu[[2L]], mu3 = theta$mu[[3L]],
    weight_uniform = theta$weight_uniform, alpha = theta$alpha, beta = theta$beta,
    loglik = loglik, AIC = -2 * loglik + 10, BIC = -2 * loglik + 5 * log(nrow(x)),
    convergence = fit$opt$convergence, optim_message = fit$opt$message %||% "",
    function_evaluations = unname(counts[["function"]] %||% NA_integer_),
    gradient_evaluations = unname(counts[["gradient"]] %||% NA_integer_),
    gradient_norm = NA_real_, gradient_norm_note = "not returned by production optim call",
    elapsed_seconds = elapsed,
    nearest_mu = min(nearest_mu), nearest_antimu = min(nearest_antimu),
    endpoint_y_clip = endpoint_y_clip,
    parameter_clip = parameter_clip,
    clipping_active = endpoint_y_clip | parameter_clip,
    returned_start_changed = changed_start,
    warning = paste(unique(warnings), collapse = " | "), error = "",
    stringsAsFactors = FALSE
  ))
}

fibonacci_sphere <- function(n) {
  i <- seq_len(n) - 0.5
  z <- 1 - 2 * i / n
  phi <- pi * (3 - sqrt(5)) * i
  cbind(sqrt(pmax(0, 1 - z^2)) * cos(phi), sqrt(pmax(0, 1 - z^2)) * sin(phi), z)
}

analyse_density_basins <- function(fits, label) {
  ok <- fits$error == "" & is.finite(fits$loglik)
  interior <- fits[ok & !fits$clipping_active, , drop = FALSE]
  if (!nrow(interior)) return(invisible(NULL))
  grid <- fibonacci_sphere(512L)
  sqrt_density <- t(vapply(seq_len(nrow(interior)), function(i) {
    z <- interior[i, ]
    sqrt((4 * pi / nrow(grid)) / 2 * d_sph_uniform_beta_mixture_s2(
      grid, c(z$mu1, z$mu2, z$mu3), z$weight_uniform, z$alpha, z$beta, log = FALSE
    ))
  }, numeric(nrow(grid))))
  gram <- tcrossprod(sqrt_density)
  norms <- rowSums(sqrt_density^2)
  h2 <- pmax(0, outer(norms, norms, "+") - 2 * gram)
  h <- matrix(sqrt(h2), nrow = nrow(interior), ncol = nrow(interior))
  # Components join only models with fitted-density Hellinger distance <= .025.
  # This deliberately does not use a parameter-distance clustering rule.
  parent <- seq_len(nrow(interior))
  root <- function(a) { while (parent[a] != a) { parent[a] <<- parent[parent[a]]; a <- parent[a] }; a }
  if (nrow(interior) >= 2L) {
    for (i in seq_len(nrow(interior) - 1L)) for (j in seq.int(i + 1L, nrow(interior))) {
      if (h[i, j] <= 0.025) { ri <- root(i); rj <- root(j); if (ri != rj) parent[rj] <- ri }
    }
  }
  component <- vapply(seq_len(nrow(interior)), root, integer(1))
  component <- match(component, unique(component))
  interior$density_component <- component
  utils::write.csv(interior, file.path(output_root, sprintf("%s_interior_fits_with_density_components.csv", label)), row.names = FALSE)
  summary <- do.call(rbind, lapply(split(interior, component), function(z) {
    data.frame(component = z$density_component[[1L]], n = nrow(z), fraction_all_starts = nrow(z) / n_fits,
      loglik_median = stats::median(z$loglik), loglik_max = max(z$loglik),
      mu1_sd = stats::sd(z$mu1), mu2_sd = stats::sd(z$mu2), mu3_sd = stats::sd(z$mu3),
      weight_sd = stats::sd(z$weight_uniform), alpha_sd = stats::sd(z$alpha), beta_sd = stats::sd(z$beta))
  }))
  summary <- summary[order(-summary$n, -summary$loglik_max), ]
  utils::write.csv(summary, file.path(output_root, sprintf("%s_density_basin_summary.csv", label)), row.names = FALSE)
  utils::write.csv(as.data.frame(as.table(h)), file.path(output_root, sprintf("%s_interior_pairwise_hellinger.csv", label)), row.names = FALSE)
}

comets <- load_comets_real_data(finite_normals = "both")
datasets <- list(short = as.matrix(comets$short$normal), long = as.matrix(comets$long$normal))
all_fits <- list()
for (label in names(datasets)) {
  x <- datasets[[label]]
  saved <- read_saved_theta(label)
  # First, reproduce the normal application call with no diagnostic warm start.
  reported_run <- uniform_beta_mixture_mle_s2_weighted(x = x, control = production_control)
  reported_theta <- reported_run[c("mu", "weight_uniform", "alpha", "beta")]
  starts <- make_start_table(x, saved, n_fits, seed + if (label == "short") 1L else 2L)
  fit_rows <- parallel::mclapply(seq_len(nrow(starts)), function(i) one_fit(starts[i, , drop = FALSE], x, label),
    mc.cores = n_cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  fits <- do.call(rbind, fit_rows)
  utils::write.csv(fits, file.path(output_root, sprintf("%s_all_fits.csv", label)), row.names = FALSE)
  utils::write.csv(data.frame(dataset = label, n = nrow(x), mu1 = reported_theta$mu[[1L]],
    mu2 = reported_theta$mu[[2L]], mu3 = reported_theta$mu[[3L]],
    weight_uniform = reported_theta$weight_uniform, alpha = reported_theta$alpha, beta = reported_theta$beta,
    projected_loglik = reported_run$loglik,
    surface_loglik = sum(d_sph_uniform_beta_mixture_s2(x, reported_theta$mu, reported_theta$weight_uniform,
      reported_theta$alpha, reported_theta$beta, log = TRUE)), convergence = reported_run$opt$convergence),
    file.path(output_root, sprintf("%s_direct_production_reproduction.csv", label)), row.names = FALSE)
  analyse_density_basins(fits, label)
  all_fits[[label]] <- fits
}

writeLines(c(
  "Repeated Uniform-Beta MLE diagnostic (separate from production outputs).",
  sprintf("n_fits per dataset: %d; cores: %d; seed: %d", n_fits, n_cores, seed),
  "Every fit calls uniform_beta_mixture_mle_s2_weighted() directly.",
  "The only added controls provide one external warm start and request warm_start_only=TRUE.",
  "If that direct production routine returns non-convergence, its own documented fallback to its ordinary candidate starts remains active.",
  "Density components use a 512-point sphere grid and Hellinger distance threshold 0.025; parameter distances are not used for grouping."
), file.path(output_root, "README.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message("Wrote diagnostic results to: ", normalizePath(output_root))
