#!/usr/bin/env Rscript

# Log-cylinder density plots matched exactly to the strict-step-four
# November-December B=5000 KS/CvM analysis.

source(file.path("real_data", "wind", "plot_risoe_cylinder_density_contours.R"))

results_dir <- file.path(
  repo_root, "real_data", "wind", "month_diagnostics",
  "nov_dec_125m_ks_cvm_b5000", "fast_sampleks_strict_step4"
)
results_csv <- file.path(
  results_dir, "risoe_125m_nov_dec_fast_sampleks_ks_cvm_b5000.csv"
)

if (!file.exists(results_csv)) {
  stop("B=5000 November-December results CSV not found: ", results_csv, call. = FALSE)
}

screening_results <- utils::read.csv(results_csv, stringsAsFactors = FALSE)
screening_results <- screening_results[screening_results$window == "nov_dec", , drop = FALSE]

month_windows <- list(nov_dec = c(11L, 12L))
window_labels <- c(nov_dec = "November + December")

pattern_rows <- screening_results[!duplicated(screening_results$pattern), , drop = FALSE]
day_patterns <- lapply(pattern_rows$day_pattern, function(x) {
  as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
})
names(day_patterns) <- pattern_rows$pattern
day_patterns <- day_patterns[paste0("start", seq_len(4L))]

validate_plot_samples <- function() {
  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  for (pattern_id in names(day_patterns)) {
    case <- build_case(selected, "nov_dec", pattern_id)
    expected_n <- unique(screening_results$n[screening_results$pattern == pattern_id])
    if (length(expected_n) != 1L || nrow(case) != expected_n) {
      stop(
        sprintf(
          "Sample mismatch for %s: plot n=%d, B=5000 result n=%s.",
          pattern_id, nrow(case), paste(expected_n, collapse = ",")
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

run_nov_dec_b5000_density_plots <- function(
    patterns = names(day_patterns),
    output_dir = file.path(results_dir, "log_cylinder_density_contours")) {
  validate_plot_samples()
  run_cylinder_plots(output_dir = output_dir, patterns = patterns)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  requested_patterns <- if (length(args) == 0L) names(day_patterns) else args
  result <- run_nov_dec_b5000_density_plots(patterns = requested_patterns)
  print(result, row.names = FALSE, digits = 5)
}
