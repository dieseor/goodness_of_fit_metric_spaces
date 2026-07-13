#!/usr/bin/env Rscript

# Fast composite-HvMF KS (sample points / unique distances) and CvM analysis
# for the preselected poorly fitting contiguous month pairs at Risoe 125 m.
# Default candidates: May-June (clearest bimodality) and June-July (strongest
# robust CvM rejection in the stored B=1000 screening).

source(file.path("real_data", "wind", "run_risoe_125m_screening_ks_cvm.R"))

candidate_windows <- list(
  may_jun = c(5L, 6L),
  jun_jul = c(6L, 7L)
)

strict_step4_patterns <- lapply(seq_len(4L), function(start_day) {
  seq.int(start_day, 30L, by = 4L)
})
names(strict_step4_patterns) <- paste0("start", seq_len(4L))

make_candidate_configs <- function() {
  configs <- list()
  index <- 1L
  for (window_id in names(candidate_windows)) {
    for (pattern_id in names(strict_step4_patterns)) {
      configs[[index]] <- list(
        dataset_id = sprintf("risoe_clean_125m_%s_%s_1996_2003", window_id, pattern_id),
        window = window_id,
        months = candidate_windows[[window_id]],
        pattern = pattern_id,
        day_pattern = strict_step4_patterns[[pattern_id]]
      )
      index <- index + 1L
    }
  }
  configs
}

run_bad_pair_candidates_fast_ks_cvm <- function(
    B = 5000L,
    n_cores = 12L,
    base_seed = 2026071300L,
    input_nc = file.path("real_data", "wind", "risoe_m_all.nc"),
    output_dir = file.path(
      "real_data", "wind", "month_diagnostics",
      sprintf("bad_pair_candidates_125m_ks_cvm_b%d", as.integer(B)),
      "fast_sampleks_strict_step4"
    )) {
  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  base_seed <- as.integer(base_seed)
  if (B < 1L || n_cores < 1L) {
    stop("`B` and `n_cores` must be positive integers.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_df <- load_risoe_concurrent(input_nc, fixed_tz = "UTC")
  selected_df <- select_noon_all_months(all_df, fixed_tz = "UTC")
  configs <- make_candidate_configs()

  rows <- lapply(seq_along(configs), function(i) {
    config <- configs[[i]]
    message(sprintf(
      "Running %s/%s: days %s; B=%d; seed=%d",
      config$window,
      config$pattern,
      paste(config$day_pattern, collapse = ","),
      B,
      base_seed + i
    ))
    run_single_risoe_125m_case(
      config = config,
      selected_df = selected_df,
      B = B,
      n_cores = n_cores,
      bootstrap_method = "reestimated",
      fixed_tz = "UTC",
      ks_grid_mode = "sample_points_unique_distances",
      seed = base_seed + i
    )
  })

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  results$months <- vapply(
    results$window,
    function(window_id) paste(candidate_windows[[window_id]], collapse = ","),
    character(1)
  )
  results$day_pattern <- vapply(
    results$pattern,
    function(pattern_id) paste(strict_step4_patterns[[pattern_id]], collapse = ","),
    character(1)
  )

  output_csv <- file.path(
    output_dir,
    sprintf("risoe_125m_bad_pair_candidates_fast_sampleks_ks_cvm_b%d.csv", B)
  )
  utils::write.csv(results, output_csv, row.names = FALSE)
  message("Results written to: ", normalizePath(output_csv, mustWork = TRUE))
  invisible(results)
}

if (sys.nframe() == 0L) {
  B <- as.integer(Sys.getenv("B", unset = "5000"))
  n_cores <- as.integer(Sys.getenv("N_CORES", unset = "12"))
  results <- run_bad_pair_candidates_fast_ks_cvm(B = B, n_cores = n_cores)
  print(results, row.names = FALSE, digits = 6)
}
