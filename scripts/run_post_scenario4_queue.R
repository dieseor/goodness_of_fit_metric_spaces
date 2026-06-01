is_scenario4_running <- function() {
  status <- suppressWarnings(system(
    "pgrep -f 'run_spherical_cauchy_calibration_m500_b500_sequential.R' >/dev/null 2>&1",
    intern = FALSE,
    ignore.stdout = TRUE,
    ignore.stderr = TRUE
  ))
  identical(as.integer(status), 0L)
}

run_step <- function(step_id, description, command, args, log_dir) {
  log_path <- file.path(log_dir, sprintf("%02d_%s.log", step_id, gsub("[^a-zA-Z0-9_]+", "_", description)))
  cat(sprintf("[%s] Starting step %02d: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step_id, description))

  exit_code <- suppressWarnings(system2(
    command = command,
    args = args,
    stdout = log_path,
    stderr = log_path
  ))

  if (!identical(as.integer(exit_code), 0L)) {
    stop(sprintf("Step %02d failed (exit code %d). Check log: %s", step_id, as.integer(exit_code), log_path))
  }

  cat(sprintf("[%s] Completed step %02d: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step_id, description))
  invisible(log_path)
}

main <- function() {
  log_dir <- file.path("output", "queued_runs", "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("[%s] Waiting for spherical Cauchy scenario 4 to finish...\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  while (is_scenario4_running()) {
    cat(sprintf("[%s] Scenario 4 still running; re-checking in 120 seconds.\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    Sys.sleep(120)
  }

  cat(sprintf("[%s] Scenario 4 finished. Launching queued tasks.\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  run_step(
    step_id = 1L,
    description = "sc_short",
    command = "Rscript",
    args = c(
      "scripts/run_comets_distance_profile_spherical_cauchy_benchmark.R",
      "--dataset=short",
      "--output_root=output/comets_distance_profile_spherical_cauchy/queued_short_ks_benchmark",
      "--statistic=ks",
      "--n_cores=12",
      "--B_values=10,30,50,100,200,500,1000"
    ),
    log_dir = log_dir
  )

  run_step(
    step_id = 2L,
    description = "jp_long",
    command = "Rscript",
    args = c(
      "scripts/run_comets_distance_profile_jp_long_benchmark.R",
      "--output_root=output/comets_distance_profile_jp/queued_long_ks_benchmark",
      "--statistic=ks",
      "--n_cores=12",
      "--B_values=10,30,50,100,200,500,1000"
    ),
    log_dir = log_dir
  )

  run_step(
    step_id = 3L,
    description = "sc_long",
    command = "Rscript",
    args = c(
      "scripts/run_comets_distance_profile_spherical_cauchy_benchmark.R",
      "--dataset=long",
      "--output_root=output/comets_distance_profile_spherical_cauchy/queued_long_ks_benchmark",
      "--statistic=ks",
      "--n_cores=12",
      "--B_values=10,30,50,100,200,500,1000"
    ),
    log_dir = log_dir
  )

  run_step(
    step_id = 4L,
    description = "vmf_composite_n50_n100",
    command = "Rscript",
    args = c(
      "scripts/run_vmf_composite_calibration_n50_n100.R",
      "--output_root=output/bootstrap_calibration/vmf_composite_M1000_B1000_n50_100_queue",
      "--n_values=50,100",
      "--M_outer=1000",
      "--B=1000",
      "--n_cores_outer=12"
    ),
    log_dir = log_dir
  )

  cat(sprintf("[%s] All queued tasks completed successfully.\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

if (sys.nframe() == 0L) {
  main()
}
