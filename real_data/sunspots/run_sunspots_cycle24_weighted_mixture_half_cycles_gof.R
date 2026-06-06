#!/usr/bin/env Rscript

resolve_cycle24_half_cycles_path <- function(...) {
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

runner_path_cycle24_half_cycles <- resolve_cycle24_half_cycles_path(
  "real_data", "sunspots", "run_sunspots_cycle24_small_circle_weighted_mixture_gof.R"
)
source(runner_path_cycle24_half_cycles)

split_cycle24_half_datasets <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output")) {
  df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  if (!"date" %in% names(df)) {
    stop("The cycle-24 CSV must contain a `date` column to split the cycle temporally.")
  }

  dates <- as.POSIXct(df$date, tz = "UTC")
  if (any(!is.finite(as.numeric(dates)))) {
    stop("Could not parse all dates in the cycle-24 CSV.")
  }

  ordering <- order(dates)
  df <- df[ordering, , drop = FALSE]
  n <- nrow(df)
  n_first <- floor(n / 2)
  n_second <- n - n_first
  if (n_first <= 0L || n_second <= 0L) {
    stop("Need at least two observations to split the cycle into two halves.")
  }

  first_half <- df[seq_len(n_first), , drop = FALSE]
  second_half <- df[(n_first + 1L):n, , drop = FALSE]

  first_path <- file.path(output_dir, "sunspots_cycle24_s2_first_half.csv")
  second_path <- file.path(output_dir, "sunspots_cycle24_s2_second_half.csv")
  utils::write.csv(first_half, first_path, row.names = FALSE)
  utils::write.csv(second_half, second_path, row.names = FALSE)

  list(
    first_path = first_path,
    second_path = second_path,
    split_summary = data.frame(
      n_total = n,
      n_first_half = n_first,
      n_second_half = n_second,
      first_date = df$date[[1L]],
      first_half_end = first_half$date[[n_first]],
      second_half_start = second_half$date[[1L]],
      last_date = df$date[[n]],
      stringsAsFactors = FALSE
    )
  )
}

run_sunspots_cycle24_weighted_mixture_half_cycles_gof <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output"),
    B = 500L,
    n_cores = 6L,
    seed_first_half = 123L,
    seed_second_half = 456L) {
  split_output <- split_cycle24_half_datasets(
    input_csv = input_csv,
    output_dir = output_dir
  )

  split_summary_path <- file.path(output_dir, "sunspots_cycle24_half_split_summary.csv")
  utils::write.csv(split_output$split_summary, split_summary_path, row.names = FALSE)

  first_dir <- file.path(output_dir, "cycle24_small_circle_weighted_mixture_first_half")
  second_dir <- file.path(output_dir, "cycle24_small_circle_weighted_mixture_second_half")

  first_result <- run_sunspots_cycle24_small_circle_weighted_mixture_gof(
    input_csv = split_output$first_path,
    output_dir = first_dir,
    dataset_label = "cycle24_first_half",
    statistics = "ks",
    B = as.integer(B),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed_first_half)
  )

  second_result <- run_sunspots_cycle24_small_circle_weighted_mixture_gof(
    input_csv = split_output$second_path,
    output_dir = second_dir,
    dataset_label = "cycle24_second_half",
    statistics = "ks",
    B = as.integer(B),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed_second_half)
  )

  combined <- rbind(first_result, second_result)
  utils::write.csv(
    combined,
    file.path(output_dir, "sunspots_cycle24_weighted_mixture_half_cycles_gof_results.csv"),
    row.names = FALSE
  )

  invisible(list(
    split_summary = split_output$split_summary,
    first_result = first_result,
    second_result = second_result,
    combined = combined
  ))
}

if (sys.nframe() == 0L) {
  output <- run_sunspots_cycle24_weighted_mixture_half_cycles_gof()
  print(output$split_summary)
  print(output$combined)
}
