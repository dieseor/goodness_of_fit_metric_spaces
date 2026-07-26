#!/usr/bin/env Rscript

# C3 pilot for the Section 6 normal and logistic-Gaussian experiments.
#
# This script sources the production runner and calls run_section6_family()
# directly.  It therefore retains the production scenario catalogue, composite
# fit, sample-point KS and CvM statistics, fast multiplier bootstrap, and fused
# C++ kernel.  Only numerical budgets/design sizes are reduced:
#
#   final: dimensions = c(2, 10), n = c(50, 100, 200, 400),
#          beta = c(0, 0.5, 1), M = 1000, B = 5000,
#          derivative_mc_size = 1000
#   pilot: dimensions = c(2, 10), n = 50, beta = c(0, 1), M = 1, B = 100,
#          derivative_mc_size = 200
#
# Both production scenarios in each family are retained.  With these settings
# there are 2 scenarios x 2 dimensions x 2 beta values x 1 replication =
# 8 jobs per family.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

required_files <- c(
  file.path("scripts", "run_section6_new_scenarios.R"),
  file.path("bootstrap", "multiplier_bootstrap.R"),
  file.path("cpp", "distance_profile_backend.cpp")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    sprintf(
      "Run this pilot from the repository root. Missing: %s",
      paste(missing_files, collapse = ", ")
    ),
    call. = FALSE
  )
}

source(file.path("scripts", "run_section6_new_scenarios.R"))

pilot_parameters <- list(
  M = 1L,
  B = 100L,
  dimensions = c(2L, 10L),
  n_values = 50L,
  beta_values = c(0, 1),
  derivative_mc_size = 200L,
  cvm_block_size = 50L,
  checkpoint_results = 2L
)

normal_seed <- 20260728L
lg_seed <- 20260727L

allocated_cores <- suppressWarnings(as.integer(Sys.getenv(
  "SLURM_CPUS_PER_TASK",
  unset = "1"
)))
if (!is.finite(allocated_cores) || allocated_cores < 1L) {
  allocated_cores <- 1L
}
pilot_cores <- min(4L, allocated_cores)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Europe/Madrid")
job_tag <- Sys.getenv(
  "SLURM_JOB_ID",
  unset = sprintf("local_%d", Sys.getpid())
)
job_tag <- gsub("[^A-Za-z0-9_.-]", "_", job_tag)
pilot_root <- file.path(
  "results",
  "c3_pilot",
  sprintf("run_%s_job_%s", timestamp, job_tag)
)
if (file.exists(pilot_root)) {
  stop(
    sprintf("Refusing to overwrite existing pilot directory: %s", pilot_root),
    call. = FALSE
  )
}
if (!dir.create(pilot_root, recursive = TRUE, showWarnings = FALSE)) {
  stop(sprintf("Could not create pilot directory: %s", pilot_root), call. = FALSE)
}

normal_output <- file.path(pilot_root, "normal")
lg_output <- file.path(pilot_root, "logistic_gaussian")
pilot_summary_path <- file.path(pilot_root, "pilot_run_summary.csv")

log_message <- function(...) {
  cat(sprintf(
    "[%s] %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Madrid"),
    sprintf(...)
  ))
  flush.console()
}

expected_replication_rows <- function(family) {
  design <- make_section6_design(
    family = family,
    dimensions = pilot_parameters$dimensions,
    n_values = pilot_parameters$n_values,
    beta_values = pilot_parameters$beta_values
  )
  as.integer(nrow(design) * pilot_parameters$M)
}

empty_summary_row <- function(family, output_dir, seed) {
  data.frame(
    family = family,
    scenarios = paste(section6_family_scenarios(family), collapse = ","),
    status = "not_run",
    error_message = NA_character_,
    started_at = NA_character_,
    finished_at = NA_character_,
    user_seconds = NA_real_,
    system_seconds = NA_real_,
    elapsed_seconds = NA_real_,
    M = pilot_parameters$M,
    B = pilot_parameters$B,
    dimensions = paste(pilot_parameters$dimensions, collapse = ","),
    n_values = paste(pilot_parameters$n_values, collapse = ","),
    beta_values = paste(pilot_parameters$beta_values, collapse = ","),
    derivative_mc_size = pilot_parameters$derivative_mc_size,
    cvm_block_size = pilot_parameters$cvm_block_size,
    checkpoint_results = pilot_parameters$checkpoint_results,
    cores = pilot_cores,
    base_seed = seed,
    expected_replication_rows = expected_replication_rows(family),
    result_rows = NA_integer_,
    ok_rows = NA_integer_,
    nonconforming_rows = NA_integer_,
    error_rows = NA_integer_,
    output_dir = output_dir,
    raw_results_path = file.path(output_dir, "raw_results.csv"),
    family_summary_path = file.path(output_dir, "summary.csv"),
    stringsAsFactors = FALSE
  )
}

family_rows <- list(
  normal = empty_summary_row("normal", normal_output, normal_seed),
  lg = empty_summary_row("lg", lg_output, lg_seed)
)

validate_family_result <- function(family, fit) {
  if (!is.list(fit) || !is.data.frame(fit$results)) {
    stop(sprintf("%s runner returned no results data frame.", family))
  }
  results <- fit$results
  expected_rows <- expected_replication_rows(family)
  if (nrow(results) != expected_rows) {
    stop(sprintf(
      "%s pilot returned %d rows; expected %d.",
      family, nrow(results), expected_rows
    ))
  }
  bad <- is.na(results$status) | results$status != "ok"
  if (any(bad)) {
    details <- paste(
      sprintf(
        "%s(rep=%s): %s%s",
        results$scenario[bad],
        results$rep[bad],
        results$status[bad],
        ifelse(
          is.na(results$error_message[bad]),
          "",
          paste0(" - ", results$error_message[bad])
        )
      ),
      collapse = "; "
    )
    stop(sprintf("%s pilot has failed/nonconforming rows: %s", family, details))
  }
  invisible(results)
}

run_family_pilot <- function(family, output_dir, seed) {
  started <- Sys.time()
  started_elapsed <- proc.time()[["elapsed"]]
  timing <- NULL
  fit <- NULL
  error <- NULL

  log_message(
    "START family=%s output=%s cores=%d M=%d B=%d",
    family, output_dir, pilot_cores, pilot_parameters$M, pilot_parameters$B
  )

  fit <- tryCatch({
    timing <- system.time({
      value <- run_section6_family(
        family = family,
        output_dir = output_dir,
        M = pilot_parameters$M,
        B = pilot_parameters$B,
        dimensions = pilot_parameters$dimensions,
        n_values = pilot_parameters$n_values,
        beta_values = pilot_parameters$beta_values,
        cores = pilot_cores,
        base_seed = seed,
        derivative_mc_size = pilot_parameters$derivative_mc_size,
        cvm_block_size = pilot_parameters$cvm_block_size,
        checkpoint_results = pilot_parameters$checkpoint_results,
        show_progress = TRUE
      )
    })
    validate_family_result(family, value)
    value
  }, error = function(e) {
    error <<- e
    NULL
  })

  finished <- Sys.time()
  elapsed_fallback <- proc.time()[["elapsed"]] - started_elapsed
  results <- if (is.list(fit) && is.data.frame(fit$results)) {
    fit$results
  } else if (file.exists(file.path(output_dir, "raw_results.csv"))) {
    utils::read.csv(
      file.path(output_dir, "raw_results.csv"),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  row <- empty_summary_row(family, output_dir, seed)
  row$status <- if (is.null(error)) "ok" else "error"
  row$error_message <- if (is.null(error)) NA_character_ else conditionMessage(error)
  row$started_at <- format(started, "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Madrid")
  row$finished_at <- format(finished, "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Madrid")
  row$user_seconds <- if (is.null(timing)) NA_real_ else unname(timing[["user.self"]])
  row$system_seconds <- if (is.null(timing)) NA_real_ else unname(timing[["sys.self"]])
  row$elapsed_seconds <- if (is.null(timing)) elapsed_fallback else unname(timing[["elapsed"]])
  if (!is.null(results)) {
    row$result_rows <- nrow(results)
    row$ok_rows <- sum(results$status == "ok")
    row$nonconforming_rows <- sum(results$status == "nonconforming")
    row$error_rows <- sum(results$status == "error")
  }

  if (is.null(error)) {
    log_message("END family=%s status=ok elapsed_seconds=%.3f", family, row$elapsed_seconds)
  } else {
    log_message(
      "END family=%s status=error elapsed_seconds=%.3f message=%s",
      family, row$elapsed_seconds, conditionMessage(error)
    )
  }

  list(row = row, error = error)
}

write_pilot_summary <- function(rows, total_started, total_elapsed, status, error_message = NA_character_) {
  total_row <- data.frame(
    family = "total",
    scenarios = NA_character_,
    status = status,
    error_message = error_message,
    started_at = format(total_started, "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Madrid"),
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Madrid"),
    user_seconds = sum(vapply(rows, function(x) x$user_seconds, numeric(1)), na.rm = TRUE),
    system_seconds = sum(vapply(rows, function(x) x$system_seconds, numeric(1)), na.rm = TRUE),
    elapsed_seconds = total_elapsed,
    M = NA_integer_,
    B = NA_integer_,
    dimensions = NA_character_,
    n_values = NA_character_,
    beta_values = NA_character_,
    derivative_mc_size = NA_integer_,
    cvm_block_size = NA_integer_,
    checkpoint_results = NA_integer_,
    cores = pilot_cores,
    base_seed = NA_integer_,
    expected_replication_rows = sum(vapply(rows, function(x) x$expected_replication_rows, integer(1))),
    result_rows = sum(vapply(rows, function(x) x$result_rows, integer(1)), na.rm = TRUE),
    ok_rows = sum(vapply(rows, function(x) x$ok_rows, integer(1)), na.rm = TRUE),
    nonconforming_rows = sum(vapply(rows, function(x) x$nonconforming_rows, integer(1)), na.rm = TRUE),
    error_rows = sum(vapply(rows, function(x) x$error_rows, integer(1)), na.rm = TRUE),
    output_dir = pilot_root,
    raw_results_path = NA_character_,
    family_summary_path = pilot_summary_path,
    stringsAsFactors = FALSE
  )
  summary <- do.call(rbind, c(rows, list(total = total_row)))
  section6_write_atomic_csv(summary, pilot_summary_path)
  invisible(summary)
}

utils::capture.output(
  {
    cat(sprintf("pilot_root: %s\n", normalizePath(pilot_root)))
    cat(sprintf("SLURM_JOB_ID: %s\n", Sys.getenv("SLURM_JOB_ID", unset = "not set")))
    cat(sprintf("SLURM_CPUS_PER_TASK: %s\n", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "not set")))
    cat(sprintf("pilot_cores: %d\n", pilot_cores))
    cat(sprintf("R_LIBS_USER: %s\n", Sys.getenv("R_LIBS_USER", unset = "not set")))
    cat(sprintf(
      "thread_limits: OMP=%s OPENBLAS=%s MKL=%s\n",
      Sys.getenv("OMP_NUM_THREADS", unset = "not set"),
      Sys.getenv("OPENBLAS_NUM_THREADS", unset = "not set"),
      Sys.getenv("MKL_NUM_THREADS", unset = "not set")
    ))
    print(sessionInfo())
  },
  file = file.path(pilot_root, "session_info.txt")
)

total_started <- Sys.time()
total_started_elapsed <- proc.time()[["elapsed"]]
log_message("PILOT START root=%s", pilot_root)

normal_run <- run_family_pilot("normal", normal_output, normal_seed)
family_rows$normal <- normal_run$row
if (!is.null(normal_run$error)) {
  total_elapsed <- proc.time()[["elapsed"]] - total_started_elapsed
  family_rows$lg$error_message <- "Not run because the normal pilot failed."
  write_pilot_summary(
    family_rows,
    total_started,
    total_elapsed,
    status = "error",
    error_message = conditionMessage(normal_run$error)
  )
  stop(
    sprintf(
      "C3 normal pilot failed. See %s and %s.",
      pilot_summary_path, normal_output
    ),
    call. = FALSE
  )
}

lg_run <- run_family_pilot("lg", lg_output, lg_seed)
family_rows$lg <- lg_run$row
total_elapsed <- proc.time()[["elapsed"]] - total_started_elapsed
if (!is.null(lg_run$error)) {
  write_pilot_summary(
    family_rows,
    total_started,
    total_elapsed,
    status = "error",
    error_message = conditionMessage(lg_run$error)
  )
  stop(
    sprintf(
      "C3 logistic-Gaussian pilot failed. See %s and %s.",
      pilot_summary_path, lg_output
    ),
    call. = FALSE
  )
}

summary <- write_pilot_summary(
  family_rows,
  total_started,
  total_elapsed,
  status = "ok"
)
log_message("PILOT END status=ok elapsed_seconds=%.3f", total_elapsed)
log_message("SUMMARY %s", pilot_summary_path)
print(summary, row.names = FALSE)
