#!/usr/bin/env Rscript

# Sunspots GOF analysis on S^2.
#
# Current target:
#   - use the sunspot-birth data from rotasym::sunspots_births;
#   - export clean S^2 files for solar cycles 21, 22, 23, and 24;
#   - represent each sunspot birth by its full spherical position on S^2;
#   - save per-cycle CSV files plus a joint 21--23 file for rolling-window
#     analyses comparable to Fernandez-de-Marcos and Garcia-Portugues;
#   - leave skeleton output files consumed by the GOF runners.

suppressPackageStartupMessages({
  library(rotasym)
  library(dplyr)
})

output_dir <- file.path("real_data", "sunspots", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cycle_targets <- 21:24

# In rotasym::sunspots_births, theta is the longitude and phi is the latitude,
# both in radians. The S^2 embedding is
#   X = (cos(phi) cos(theta), cos(phi) sin(theta), sin(phi)).
sunspots_raw <- rotasym::sunspots_births

sunspots_cycles_s2 <- sunspots_raw %>%
  filter(cycle %in% cycle_targets) %>%
  mutate(
    hemisphere = case_when(
      phi > 0 ~ "N",
      phi < 0 ~ "S",
      TRUE ~ "Equator"
    ),
    x1 = cos(phi) * cos(theta),
    x2 = cos(phi) * sin(theta),
    x3 = sin(phi),
    norm_s2 = sqrt(x1^2 + x2^2 + x3^2)
  )

for (cycle_target in cycle_targets) {
  sunspots_cycle_s2 <- sunspots_cycles_s2 %>%
    filter(cycle == cycle_target)

  utils::write.csv(
    sunspots_cycle_s2,
    file.path(output_dir, sprintf("sunspots_cycle%d_s2_all.csv", cycle_target)),
    row.names = FALSE
  )

  sunspots_cycle_s2 %>%
    filter(hemisphere == "N") %>%
    utils::write.csv(file.path(output_dir, sprintf("sunspots_cycle%d_s2_north.csv", cycle_target)), row.names = FALSE)

  sunspots_cycle_s2 %>%
    filter(hemisphere == "S") %>%
    utils::write.csv(file.path(output_dir, sprintf("sunspots_cycle%d_s2_south.csv", cycle_target)), row.names = FALSE)
}

sunspots_cycles21_23_s2 <- sunspots_cycles_s2 %>%
  filter(cycle %in% 21:23)
utils::write.csv(
  sunspots_cycles21_23_s2,
  file.path(output_dir, "sunspots_cycles21_23_s2_all.csv"),
  row.names = FALSE
)

summary_table <- sunspots_cycles_s2 %>%
  group_by(cycle) %>%
  summarise(
    n = n(),
    n_north = sum(hemisphere == "N"),
    n_south = sum(hemisphere == "S"),
    n_equator = sum(hemisphere == "Equator"),
    min_norm_s2 = min(norm_s2),
    max_norm_s2 = max(norm_s2),
    max_abs_norm_error = max(abs(norm_s2 - 1))
  )

utils::write.csv(
  summary_table,
  file.path(output_dir, "sunspots_cycles21_24_s2_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  summary_table %>% filter(cycle == 24L),
  file.path(output_dir, "sunspots_cycle24_s2_summary.csv"),
  row.names = FALSE
)

# Skeleton of the final GOF output table for the symmetric two-small-circles model.
gof_results_skeleton <- tibble::tibble(
  cycle = integer(),
  dataset = character(),
  model = character(),
  mixture_components = integer(),
  statistic_type = character(),
  n = integer(),
  loglik = double(),
  aic = double(),
  bic = double(),
  test_statistic = double(),
  p_value_raw = double(),
  mu_1 = double(),
  mu_2 = double(),
  mu_3 = double(),
  kappa_hat = double(),
  nu_hat = double()
)

utils::write.csv(
  gof_results_skeleton,
  file.path(output_dir, "sunspots_cycle24_small_circle_symmetric_mixture_gof_results.csv"),
  row.names = FALSE
)

print(summary_table)
message("Saved cycle-21--24 S^2 sunspot files in: ", output_dir)
