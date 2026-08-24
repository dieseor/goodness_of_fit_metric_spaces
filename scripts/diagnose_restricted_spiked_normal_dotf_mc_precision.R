#!/usr/bin/env Rscript

# Precision/timing diagnostic for the score-Monte-Carlo estimator of dot F in
# the restricted mean-aligned single-spiked normal model.  This script is
# intentionally isolated: it neither runs GOF nor changes any common pipeline.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  value <- args[startsWith(args, key)]
  if (!length(value)) return(default)
  sub(key, "", value[[1L]], fixed = TRUE)
}
parse_csv_integer <- function(value, name) {
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (!length(out) || any(!is.finite(out)) || any(out <= 0L)) {
    stop("`--", name, "` must contain positive integers.")
  }
  out
}

cores <- as.integer(arg_value("cores", "1"))
mc_sizes <- parse_csv_integer(arg_value("mc_sizes", "1000,10000"), "mc_sizes")
replications <- as.integer(arg_value("replications", "20"))
reference_mc_size <- as.integer(arg_value("reference_mc_size", "200000"))
n_centers <- as.integer(arg_value("n_centers", "12"))
radius_quantiles <- as.numeric(strsplit(arg_value("radius_quantiles", "0.1,0.3,0.5,0.7,0.9"), ",", fixed = TRUE)[[1L]])
radius_pilot_size <- as.integer(arg_value("radius_pilot_size", "20000"))
lambda <- as.numeric(arg_value("lambda", "2"))
seed <- as.integer(arg_value("seed", "20260834"))
output_dir <- arg_value(
  "output_dir",
  "simulation_results/restricted_spiked_normal_dotf_mc_precision_lambda2"
)

if (any(!is.finite(c(cores, mc_sizes, replications, reference_mc_size, n_centers,
                     radius_quantiles, radius_pilot_size, lambda, seed))) ||
    cores < 1L || replications < 1L || reference_mc_size < 1L ||
    n_centers < 1L || radius_pilot_size < 1L || lambda <= 0 ||
    any(radius_quantiles <= 0 | radius_quantiles >= 1)) {
  stop("Invalid diagnostic settings.")
}

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theta_from_design <- function(d, design) {
  switch(
    design,
    axis_075 = c(0.75, rep(0, d - 1L)),
    axis_100 = c(1, rep(0, d - 1L)),
    diagonal_100 = rep(1 / sqrt(d), d),
    diagonal_150 = rep(1.5 / sqrt(d), d),
    stop("Unknown theta design.")
  )
}

config <- do.call(rbind, lapply(c(2L, 5L), function(d) {
  data.frame(
    config_id = sprintf("d%d_%s", d, c("axis_075", "axis_100", "diagonal_100", "diagonal_150")),
    d = d,
    theta_design = c("axis_075", "axis_100", "diagonal_100", "diagonal_150"),
    stringsAsFactors = FALSE
  )
}))
config$theta <- vapply(seq_len(nrow(config)), function(i) {
  paste(theta_from_design(config$d[[i]], config$theta_design[[i]]), collapse = ";")
}, character(1L))
utils::write.csv(config, file.path(output_dir, "configurations.csv"), row.names = FALSE)
max_parameter_dimension <- max(config$d) + 1L

build_evaluation_grid <- function(row, grid_seed) {
  set.seed(grid_seed)
  d <- as.integer(row$d)
  theta <- theta_from_design(d, as.character(row$theta_design))
  centers <- rrestricted_spiked_normal(n_centers, theta, lambda)
  radius_pilot <- rrestricted_spiked_normal(radius_pilot_size, theta, lambda)
  distance_matrix <- sqrt(pmax(
    rowSums(radius_pilot^2) + rep(rowSums(centers^2), each = nrow(radius_pilot)) -
      2 * radius_pilot %*% t(centers),
    0
  ))
  # Each centre receives thresholds from its own null distance distribution.
  radii <- vapply(seq_len(n_centers), function(j) {
    stats::quantile(distance_matrix[, j], probs = radius_quantiles, names = FALSE, type = 8)
  }, numeric(length(radius_quantiles)))
  radii <- t(matrix(
    radii,
    nrow = length(radius_quantiles),
    ncol = n_centers
  ))
  grid <- expand.grid(center_id = seq_len(n_centers), radius_id = seq_along(radius_quantiles))
  grid$omega <- vapply(seq_len(nrow(grid)), function(i) {
    paste(centers[grid$center_id[[i]], ], collapse = ";")
  }, character(1L))
  grid$t <- vapply(seq_len(nrow(grid)), function(i) {
    radii[grid$center_id[[i]], grid$radius_id[[i]]]
  }, numeric(1L))
  grid$target_quantile <- radius_quantiles[grid$radius_id]
  grid
}

estimate_dotf <- function(theta, grid, n_aux, simulation_seed) {
  set.seed(simulation_seed)
  started <- proc.time()[["elapsed"]]
  aux <- rrestricted_spiked_normal(n_aux, theta, lambda)
  psi <- restricted_spiked_normal_score_matrix(aux, list(theta = theta, lambda = lambda))
  omega <- do.call(rbind, strsplit(grid$omega, ";", fixed = TRUE))
  storage.mode(omega) <- "double"
  distance_matrix <- sqrt(pmax(
    rowSums(aux^2) + rep(rowSums(omega^2), each = nrow(aux)) - 2 * aux %*% t(omega),
    0
  ))
  dotf <- matrix(NA_real_, nrow = nrow(grid), ncol = ncol(psi))
  for (j in seq_len(nrow(grid))) {
    dotf[j, ] <- colMeans(psi * (distance_matrix[, j] <= grid$t[[j]]))
  }
  list(dotf = dotf, elapsed_seconds = proc.time()[["elapsed"]] - started)
}

grid_path <- file.path(output_dir, "evaluation_grid.csv")
grid_list <- lapply(seq_len(nrow(config)), function(i) {
  row <- config[i, , drop = FALSE]
  grid <- build_evaluation_grid(row, seed + 100003L * i)
  grid$config_id <- row$config_id
  grid$d <- row$d
  grid$theta_design <- row$theta_design
  grid
})
grid_all <- do.call(rbind, grid_list)
utils::write.csv(grid_all, grid_path, row.names = FALSE)

reference_path <- file.path(output_dir, "dotf_reference.csv")
reference <- if (file.exists(reference_path)) {
  utils::read.csv(reference_path, stringsAsFactors = FALSE)
} else {
  reference_rows <- parallel::mclapply(seq_len(nrow(config)), function(i) {
    row <- config[i, , drop = FALSE]
    theta <- theta_from_design(row$d, row$theta_design)
    grid <- grid_all[grid_all$config_id == row$config_id, , drop = FALSE]
    estimate <- estimate_dotf(theta, grid, reference_mc_size, seed + 2000003L + 100003L * i)
    out <- grid
    out$reference_elapsed_seconds <- estimate$elapsed_seconds
    for (j in seq_len(max_parameter_dimension)) {
      out[[paste0("dotf_", j)]] <- if (j <= ncol(estimate$dotf)) {
        estimate$dotf[, j]
      } else {
        NA_real_
      }
    }
    out
  }, mc.cores = min(cores, nrow(config)), mc.preschedule = FALSE)
  reference <- do.call(rbind, reference_rows)
  utils::write.csv(reference, reference_path, row.names = FALSE)
  reference
}

jobs <- expand.grid(
  config_id = config$config_id,
  mc_size = mc_sizes,
  replication = seq_len(replications),
  stringsAsFactors = FALSE
)
jobs$job_id <- seq_len(nrow(jobs))
results_path <- file.path(output_dir, "dotf_estimation_runs.csv")
existing <- if (file.exists(results_path)) utils::read.csv(results_path, stringsAsFactors = FALSE) else data.frame()
point_results_path <- file.path(output_dir, "dotf_point_errors.csv")
existing_points <- if (file.exists(point_results_path)) {
  utils::read.csv(point_results_path, stringsAsFactors = FALSE)
} else data.frame()
done <- if (nrow(existing)) as.integer(existing$job_id) else integer()
pending <- jobs[!jobs$job_id %in% done, , drop = FALSE]
message(sprintf("Restricted-spiked dot-F MC diagnostic: %d pending / %d runs; cores=%d.",
                nrow(pending), nrow(jobs), cores))

run_job <- function(job) {
  config_row <- config[config$config_id == job$config_id, , drop = FALSE]
  theta <- theta_from_design(config_row$d[[1L]], config_row$theta_design[[1L]])
  grid <- grid_all[grid_all$config_id == job$config_id, , drop = FALSE]
  ref <- reference[reference$config_id == job$config_id, , drop = FALSE]
  reference_matrix <- as.matrix(ref[paste0("dotf_", seq_len(length(theta) + 1L))])
  started <- proc.time()[["elapsed"]]
  tryCatch({
    estimate <- estimate_dotf(theta, grid, as.integer(job$mc_size),
                              seed + 5000003L + 100003L * as.integer(job$job_id))
    error <- estimate$dotf - reference_matrix
    reference_rms <- sqrt(mean(reference_matrix^2))
    point_ref_norm <- sqrt(rowSums(reference_matrix^2))
    point_error_norm <- sqrt(rowSums(error^2))
    summary_row <- data.frame(
      job_id = as.integer(job$job_id), config_id = as.character(job$config_id),
      d = config_row$d[[1L]], theta_design = config_row$theta_design[[1L]],
      theta_norm = sqrt(sum(theta^2)), lambda = lambda,
      mc_size = as.integer(job$mc_size), replication = as.integer(job$replication),
      status = "ok", error_message = NA_character_,
      elapsed_seconds = estimate$elapsed_seconds,
      coordinate_rmse = sqrt(mean(error^2)),
      relative_rmse = sqrt(mean(error^2)) / reference_rms,
      median_point_l2_error = stats::median(point_error_norm),
      median_point_relative_l2_error = stats::median(point_error_norm / pmax(point_ref_norm, 1e-12)),
      max_point_l2_error = max(point_error_norm),
      stringsAsFactors = FALSE
    )
    point_rows <- data.frame(
      job_id = as.integer(job$job_id), config_id = as.character(job$config_id),
      d = config_row$d[[1L]], theta_design = config_row$theta_design[[1L]],
      mc_size = as.integer(job$mc_size), replication = as.integer(job$replication),
      point_id = seq_len(nrow(grid)), center_id = grid$center_id,
      radius_id = grid$radius_id, target_quantile = grid$target_quantile, t = grid$t,
      reference_l2_norm = point_ref_norm, l2_error = point_error_norm,
      relative_l2_error = point_error_norm / pmax(point_ref_norm, 1e-12),
      coordinate_rmse = sqrt(rowMeans(error^2)),
      coordinate_max_abs_error = apply(abs(error), 1L, max),
      stringsAsFactors = FALSE
    )
    list(summary = summary_row, points = point_rows)
  }, error = function(error) {
    list(summary = data.frame(
      job_id = as.integer(job$job_id), config_id = as.character(job$config_id),
      d = config_row$d[[1L]], theta_design = config_row$theta_design[[1L]],
      theta_norm = sqrt(sum(theta^2)), lambda = lambda,
      mc_size = as.integer(job$mc_size), replication = as.integer(job$replication),
      status = "error", error_message = conditionMessage(error), elapsed_seconds = proc.time()[["elapsed"]] - started,
      coordinate_rmse = NA_real_, relative_rmse = NA_real_, median_point_l2_error = NA_real_,
      median_point_relative_l2_error = NA_real_, max_point_l2_error = NA_real_, stringsAsFactors = FALSE
    ), points = NULL)
  })
}

started_all <- proc.time()[["elapsed"]]
if (nrow(pending)) {
  for (begin in seq.int(1L, nrow(pending), by = cores)) {
    end <- min(nrow(pending), begin + cores - 1L)
    batch <- pending[begin:end, , drop = FALSE]
    outcomes <- parallel::mclapply(seq_len(nrow(batch)), function(i) run_job(batch[i, , drop = FALSE]),
                                   mc.cores = min(cores, nrow(batch)), mc.preschedule = FALSE)
    existing <- rbind(existing, do.call(rbind, lapply(outcomes, `[[`, "summary")))
    point_rows <- Filter(Negate(is.null), lapply(outcomes, `[[`, "points"))
    if (length(point_rows)) {
      existing_points <- rbind(existing_points, do.call(rbind, point_rows))
      utils::write.csv(existing_points, point_results_path, row.names = FALSE)
    }
    utils::write.csv(existing, results_path, row.names = FALSE)
    elapsed <- proc.time()[["elapsed"]] - started_all
    message(sprintf("completed %d/%d runs (this run %d/%d); elapsed %.1fs; ETA %.1fs",
                    nrow(existing), nrow(jobs), end, nrow(pending), elapsed,
                    elapsed / end * (nrow(pending) - end)))
  }
}

ok <- existing[existing$status == "ok", , drop = FALSE]
summary <- do.call(rbind, lapply(split(ok, interaction(ok$config_id, ok$mc_size, drop = TRUE)), function(z) {
  data.frame(
    config_id = z$config_id[[1L]], d = z$d[[1L]], theta_design = z$theta_design[[1L]],
    theta_norm = z$theta_norm[[1L]], lambda = z$lambda[[1L]], mc_size = z$mc_size[[1L]],
    completed = nrow(z), failures = sum(existing$status != "ok" &
      existing$config_id == z$config_id[[1L]] & existing$mc_size == z$mc_size[[1L]]),
    mean_coordinate_rmse = mean(z$coordinate_rmse),
    median_coordinate_rmse = stats::median(z$coordinate_rmse),
    mean_relative_rmse = mean(z$relative_rmse),
    median_relative_rmse = stats::median(z$relative_rmse),
    mean_median_point_l2_error = mean(z$median_point_l2_error),
    mean_median_point_relative_l2_error = mean(z$median_point_relative_l2_error),
    mean_max_point_l2_error = mean(z$max_point_l2_error),
    mean_elapsed_seconds = mean(z$elapsed_seconds),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(summary, file.path(output_dir, "dotf_precision_summary.csv"), row.names = FALSE)
writeLines(c(
  "restricted single-spiked normal dot-F score-MC precision diagnostic",
  sprintf("lambda=%g; cores=%d; seed=%d", lambda, cores, seed),
  sprintf("candidate auxiliary sizes: %s; replications each: %d", paste(mc_sizes, collapse = ","), replications),
  sprintf("reference auxiliary size: %d", reference_mc_size),
  sprintf("evaluation grid: %d centres x %d conditional-distance quantiles = %d points/configuration",
          n_centers, length(radius_quantiles), n_centers * length(radius_quantiles)),
  "configurations: d=2,5 and theta=.75e1, e1, 1.5*1/sqrt(d), 1/sqrt(d)",
  "target: E[1{||X-omega||<=t} psi_(theta,lambda)(X)] at the true restricted null parameter",
  "dotf_point_errors.csv stores every evaluation-point error, not only its aggregate.",
  "This diagnoses auxiliary Monte Carlo error only; it does not run GOF or assess plug-in theta error."
), file.path(output_dir, "manifest.txt"))
message("Summary: ", normalizePath(file.path(output_dir, "dotf_precision_summary.csv"), mustWork = FALSE))
