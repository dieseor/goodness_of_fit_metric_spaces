source("bootstrap/calibration_study.R")

output_dir <- file.path("output", "calibration", "bootstrap", "hvmf_composite_cvm_full")

res <- run_full_hvmf_composite_cvm_calibration_study(
  output_dir = output_dir,
  n_cores_outer = 10,
  M_outer = 200,
  B = 200,
  seed = 456,
  show_progress = TRUE,
  verbose = FALSE
)

saveRDS(res, file.path(output_dir, "run_object_M200_B200.rds"))
