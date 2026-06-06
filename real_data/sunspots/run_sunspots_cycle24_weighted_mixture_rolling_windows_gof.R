#!/usr/bin/env Rscript

resolve_cycle24_weighted_rolling_wrapper_path <- function(...) {
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

general_runner_path_cycle24_weighted_rolling <- resolve_cycle24_weighted_rolling_wrapper_path(
  "real_data", "sunspots", "run_sunspots_weighted_mixture_rolling_windows_gof.R"
)
source(general_runner_path_cycle24_weighted_rolling)

run_sunspots_cycle24_weighted_mixture_rolling_windows_gof <- function(...) {
  args <- list(...)
  if (is.null(args$cycles)) {
    args$cycles <- 24L
  }
  if (is.null(args$input_csv)) {
    args$input_csv <- file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv")
  }
  if (is.null(args$output_dir)) {
    args$output_dir <- file.path("real_data", "sunspots", "output", "cycle24_weighted_mixture_rolling_windows_10cr_B1000")
  }
  do.call(run_sunspots_weighted_mixture_rolling_windows_gof, args)
}

if (sys.nframe() == 0L) {
  cli_args <- parse_sunspots_weighted_rolling_args()
  if (is.null(cli_args$cycles)) {
    cli_args$cycles <- 24L
  }
  if (is.null(cli_args$input_csv)) {
    cli_args$input_csv <- file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv")
  }
  if (is.null(cli_args$output_dir)) {
    cli_args$output_dir <- file.path("real_data", "sunspots", "output", "cycle24_weighted_mixture_rolling_windows_10cr_B1000")
  }
  output <- do.call(run_sunspots_weighted_mixture_rolling_windows_gof, cli_args)
  print(output$summary_by_cycle)
  cat(sprintf("Total elapsed seconds: %.3f\n", output$elapsed_total))
}
