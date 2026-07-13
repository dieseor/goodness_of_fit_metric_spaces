resolve_comets_sc_fast_path <- function(...) {
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

source(resolve_comets_sc_fast_path(
  "scripts",
  "run_comets_distance_profile_spherical_cauchy_benchmark.R"
))
source(resolve_comets_sc_fast_path("scripts", "path_helpers.R"))

parse_named_args_sc_fast <- function(args) {
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

run_comets_spherical_cauchy_short_long_fast <- function(output_root = NULL,
                                                        B = 5000L,
                                                        n_cores = 8L,
                                                        ks_t_points = 250L,
                                                        base_seed = 20260708L,
                                                        distance_type = "geodesic",
                                                        control = list(
                                                          spherical_cauchy_maxit = 500L,
                                                          spherical_cauchy_reltol = 1e-10,
                                                          spherical_cauchy_optim_method = "BFGS",
                                                          spherical_cauchy_use_gradient = TRUE,
                                                          spherical_cauchy_profile_tol = 1e-10,
                                                          spherical_cauchy_profile_warn = FALSE
                                                        )) {
  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  ks_t_points <- as.integer(ks_t_points)
  base_seed <- as.integer(base_seed)

  if (!is.finite(B) || B <= 0L) {
    stop("`B` must be a strictly positive integer.")
  }
  if (!is.finite(n_cores) || n_cores <= 0L) {
    stop("`n_cores` must be a strictly positive integer.")
  }
  if (!is.finite(ks_t_points) || ks_t_points <= 0L) {
    stop("`ks_t_points` must be a strictly positive integer.")
  }

  if (is.null(output_root)) {
    output_root <- canonical_comets_spherical_cauchy_dir(
      sprintf("paper_results_B%d_sampleks", B),
      "fast"
    )
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  job_grid <- data.frame(
    dataset = c("short", "short", "long", "long"),
    statistic = c("cvm", "ks", "cvm", "ks"),
    stage_id = sprintf("%02d", 1:4),
    stringsAsFactors = FALSE
  )

  summary_rows <- vector("list", nrow(job_grid))

  for (i in seq_len(nrow(job_grid))) {
    stage_dir <- file.path(
      output_root,
      sprintf("%s_%s_%s", job_grid$stage_id[[i]], job_grid$dataset[[i]], job_grid$statistic[[i]])
    )

    message(sprintf(
      "[spherical_cauchy fast comets] %d/%d: dataset=%s statistic=%s B=%d",
      i,
      nrow(job_grid),
      job_grid$dataset[[i]],
      toupper(job_grid$statistic[[i]]),
      B
    ))

    run_comets_distance_profile_spherical_cauchy_benchmark(
      output_root = stage_dir,
      dataset = job_grid$dataset[[i]],
      B_values = B,
      statistic = job_grid$statistic[[i]],
      n_cores = n_cores,
      ks_t_points = ks_t_points,
      base_seed = base_seed + 100L * i,
      bootstrap_method = "fast_multiplier",
      distance_type = distance_type,
      control = control
    )

    summary_path <- file.path(stage_dir, "benchmark_summary.csv")
    if (!file.exists(summary_path)) {
      stop(sprintf("Expected benchmark summary was not created: %s", summary_path))
    }
    summary_rows[[i]] <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  }

  combined_summary <- do.call(rbind, summary_rows)
  utils::write.csv(
    combined_summary,
    file = file.path(output_root, "comets_spherical_cauchy_short_long_fast_summary.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      output_root = output_root,
      B = B,
      n_cores = n_cores,
      ks_t_points = ks_t_points,
      base_seed = base_seed,
      distance_type = distance_type,
      bootstrap_method = "fast_multiplier",
      jobs = job_grid,
      summary = combined_summary
    ),
    file = file.path(output_root, "comets_spherical_cauchy_short_long_fast_run.rds")
  )

  invisible(combined_summary)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_sc_fast(commandArgs(trailingOnly = TRUE))

  run_comets_spherical_cauchy_short_long_fast(
    output_root = args$output_root %||% NULL,
    B = as.integer(args$B %||% 5000L),
    n_cores = as.integer(args$n_cores %||% 8L),
    ks_t_points = as.integer(args$ks_t_points %||% 250L),
    base_seed = as.integer(args$seed %||% 20260708L),
    distance_type = as.character(args$distance_type %||% "geodesic")
  )
}
