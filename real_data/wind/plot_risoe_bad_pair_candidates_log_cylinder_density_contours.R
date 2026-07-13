#!/usr/bin/env Rscript

# Log-cylinder parametric/nonparametric density comparison for the two
# preselected poorly fitting contiguous month pairs at Risoe 125 m.

source(file.path("real_data", "wind", "plot_risoe_cylinder_density_contours.R"))

month_windows <- list(
  may_jun = c(5L, 6L),
  jun_jul = c(6L, 7L)
)

window_labels <- c(
  may_jun = "May + June",
  jun_jul = "June + July"
)

# Match exactly the strict step-four samples used in the stored B=1000
# screening and in the B=5000 rerun script.
day_patterns <- lapply(seq_len(4L), function(start_day) {
  seq.int(start_day, 30L, by = 4L)
})
names(day_patterns) <- paste0("start", seq_len(4L))

run_bad_pair_density_plots <- function(
    patterns = names(day_patterns),
    output_dir = file.path(
      repo_root, "real_data", "wind", "bad_pair_log_cylinder_density_contours"
    )) {
  run_cylinder_plots(output_dir = output_dir, patterns = patterns)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  requested_patterns <- if (length(args) == 0L) names(day_patterns) else args
  result <- run_bad_pair_density_plots(patterns = requested_patterns)
  print(result, row.names = FALSE, digits = 5)
}
