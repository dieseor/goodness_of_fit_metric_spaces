resolve_c3_long_seed_sensitivity_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_c3_long_seed_sensitivity_path("scripts", "run_comets_distance_profile_cardioid.R"))

parse_c3_long_seed_sensitivity_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else TRUE
  }
  out
}

run_comets_cardioid_c3_long_seed_sensitivity <- function(
    output_root = "real_data/comets/cardioid/c3_long_seed_sensitivity_B5000/fast",
    seeds = 20260715L:20260719L,
    B = 5000L,
    n_cores = 1L) {
  seeds <- as.integer(seeds)
  if (length(seeds) != 5L || any(!is.finite(seeds))) {
    stop("`seeds` must contain exactly five finite integer seeds.")
  }
  if (as.integer(n_cores) != 1L) stop("This sensitivity runner is intentionally restricted to one core.")
  if (!is.finite(B) || B <= 0L) stop("`B` must be a positive integer.")

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  summaries <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    seed <- seeds[[i]]
    seed_root <- file.path(output_root, sprintf("seed_%d", seed))
    completed_path <- file.path(seed_root, "pipeline_result.rds")
    if (file.exists(completed_path)) {
      message(sprintf("[C3 long seed sensitivity] %d/%d: seed=%d (already completed)", i, length(seeds), seed))
      result <- readRDS(completed_path)
    } else {
      message(sprintf("[C3 long seed sensitivity] %d/%d: seed=%d", i, length(seeds), seed))
      result <- run_comets_distance_profile_cardioid(
        output_root = seed_root,
        stages = c("oort_cvm", "oort_ks"),
        model_ids = "C3",
        cvm_B = as.integer(B),
        ks_B = as.integer(B),
        ks_grid_mode = "sample_points_unique_distances",
        n_cores = 1L,
        base_seed = seed,
        bootstrap_method = "fast_multiplier",
        distance_type = "geodesic",
        control = list(cardioid_optim_control = list(maxit = 1000))
      )
    }
    summaries[[i]] <- do.call(rbind, lapply(result$stages, `[[`, "summary"))
    summaries[[i]]$base_seed <- seed
  }

  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, file.path(output_root, "c3_long_seed_sensitivity_summary.csv"), row.names = FALSE)
  saveRDS(
    list(output_root = output_root, seeds = seeds, B = as.integer(B), n_cores = 1L, summary = summary),
    file.path(output_root, "c3_long_seed_sensitivity_run.rds")
  )
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- parse_c3_long_seed_sensitivity_args(commandArgs(trailingOnly = TRUE))
  seeds <- if (is.null(args$seeds)) 20260715L:20260719L else as.integer(strsplit(args$seeds, ",", fixed = TRUE)[[1L]])
  run_comets_cardioid_c3_long_seed_sensitivity(
    output_root = args$output_root %||% "real_data/comets/cardioid/c3_long_seed_sensitivity_B5000/fast",
    seeds = seeds,
    B = as.integer(args$B %||% 5000L),
    n_cores = as.integer(args$n_cores %||% 1L)
  )
}
