resolve_comets_watson_fast_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_comets_watson_fast_path("scripts", "run_comets_distance_profile_watson_benchmark.R"))
source(resolve_comets_watson_fast_path("scripts", "path_helpers.R"))

run_comets_watson_short_long_fast <- function(output_root = NULL,
                                              B = 5000L,
                                              n_cores = 8L,
                                              base_seed = 20260714L,
                                              distance_type = "geodesic",
                                              control = list(
                                                watson_L_max = 200L,
                                                watson_quad_n = 400L,
                                                watson_tol = 1e-10,
                                                derivative_mc_size = 1000L
                                              )) {
  if (is.null(output_root)) output_root <- canonical_comets_watson_dir(sprintf("paper_results_B%d_sampleks", as.integer(B)), "fast")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  jobs <- data.frame(dataset = rep(c("short", "long"), each = 2L), statistic = rep(c("ks", "cvm"), 2L), stringsAsFactors = FALSE)
  summaries <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    stage_dir <- file.path(output_root, sprintf("%02d_%s_%s", i, jobs$dataset[[i]], jobs$statistic[[i]]))
    message(sprintf("[watson fast comets] %d/%d: dataset=%s statistic=%s B=%d", i, nrow(jobs), jobs$dataset[[i]], toupper(jobs$statistic[[i]]), B))
    summaries[[i]] <- run_comets_distance_profile_watson_benchmark(
      output_root = stage_dir, dataset = jobs$dataset[[i]], statistic = jobs$statistic[[i]],
      B = B, n_cores = n_cores,
      seed = as.integer(base_seed + 100L * i), distance_type = distance_type,
      control = utils::modifyList(control, list(derivative_mc_seed = as.integer(base_seed + 1000L * i)))
    )
  }
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, file.path(output_root, "comets_watson_short_long_fast_summary.csv"), row.names = FALSE)
  saveRDS(list(output_root = output_root, B = B, n_cores = n_cores, base_seed = base_seed, jobs = jobs, summary = summary), file.path(output_root, "comets_watson_short_long_fast_run.rds"))
  invisible(summary)
}

parse_watson_comets_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else TRUE
  }
  out
}

if (sys.nframe() == 0L) {
  args <- parse_watson_comets_args(commandArgs(trailingOnly = TRUE))
  run_comets_watson_short_long_fast(
    output_root = args$output_root %||% NULL,
    B = as.integer(args$B %||% 5000L),
    n_cores = as.integer(args$n_cores %||% 8L),
    base_seed = as.integer(args$seed %||% 20260714L),
    distance_type = as.character(args$distance_type %||% "geodesic")
  )
}
