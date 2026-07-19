#!/usr/bin/env Rscript

# Fast composite-HvMF KS (sample points / unique distances) and CvM analysis
# for the cleaned Risoe 77 m November-December data, 1996-2003.

source(file.path("real_data", "wind", "run_risoe_125m_screening_ks_cvm.R"))

strict_step4_patterns <- lapply(seq_len(4L), function(start_day) {
  seq.int(start_day, 30L, by = 4L)
})
names(strict_step4_patterns) <- paste0("start", seq_len(4L))

make_nov_dec_77m_configs <- function() {
  lapply(names(strict_step4_patterns), function(pattern_id) {
    list(
      dataset_id = sprintf("risoe_clean_77m_nov_dec_%s_1996_2003", pattern_id),
      window = "nov_dec",
      months = c(11L, 12L),
      pattern = pattern_id,
      day_pattern = strict_step4_patterns[[pattern_id]]
    )
  })
}

run_nov_dec_77m_fast_ks_cvm <- function(
    B = 5000L,
    n_cores = 12L,
    base_seed = 2026071300L,
    input_nc = file.path("real_data", "wind", "risoe_m_all.nc"),
    output_dir = file.path(
      "real_data", "wind", "month_diagnostics",
      sprintf("nov_dec_77m_ks_cvm_b%d", as.integer(B)),
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
  configs <- make_nov_dec_77m_configs()

  rows <- lapply(seq_along(configs), function(i) {
    config <- configs[[i]]
    message(sprintf(
      "Running %s: days %s; B=%d; seed=%d",
      config$pattern, paste(config$day_pattern, collapse = ","),
      B, base_seed + i
    ))
    run_single_risoe_case(
      config = config,
      selected_df = selected_df,
      speed_col = "ws77",
      direction_col = "wd77",
      height_m = 77L,
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
  results$height_m <- 77L
  results$day_pattern <- vapply(
    results$pattern,
    function(pattern_id) paste(strict_step4_patterns[[pattern_id]], collapse = ","),
    character(1)
  )

  output_csv <- file.path(
    output_dir,
    sprintf("risoe_77m_nov_dec_fast_sampleks_ks_cvm_b%d.csv", B)
  )
  utils::write.csv(results, output_csv, row.names = FALSE)
  message("Results written to: ", normalizePath(output_csv, mustWork = TRUE))
  invisible(results)
}

if (sys.nframe() == 0L) {
  B <- as.integer(Sys.getenv("B", unset = "5000"))
  n_cores <- as.integer(Sys.getenv("N_CORES", unset = "12"))
  results <- run_nov_dec_77m_fast_ks_cvm(B = B, n_cores = n_cores)
  print(results, row.names = FALSE, digits = 6)
}
