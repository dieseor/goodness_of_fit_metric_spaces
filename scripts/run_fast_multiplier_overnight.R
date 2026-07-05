resolve_fast_overnight_path <- function(...) {
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

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

source(resolve_fast_overnight_path("scripts", "path_helpers.R"))

parse_fast_overnight_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

run_logged_step <- function(step_index,
                            n_steps,
                            label,
                            command,
                            args,
                            env = character(0),
                            log_dir) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_]+", "_", tolower(label))
  log_path <- file.path(log_dir, sprintf("%02d_%s.log", step_index, safe_label))
  border <- paste(rep("=", 72), collapse = "")
  message(border)
  message(sprintf("[step %d/%d] %s", step_index, n_steps, label))
  message(sprintf("[step %d/%d] log: %s", step_index, n_steps, log_path))
  message(border)
  start_time <- Sys.time()
  status <- system2(command, args = args, stdout = log_path, stderr = log_path, env = env)
  elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (!identical(status, 0L)) {
    stop(sprintf("Step failed: %s. See %s", label, log_path))
  }
  message(sprintf("[step %d/%d] completed in %.1f seconds", step_index, n_steps, elapsed_seconds))
  invisible(list(log_path = log_path, elapsed_seconds = elapsed_seconds))
}

run_fast_multiplier_overnight <- function(output_root = file.path("output"),
                                          n_cores = 12L,
                                          derivative_mc_size = 1000L,
                                          derivative_mc_seed = 20260613L,
                                          validation_B = 1000L,
                                          validation_M_outer = 5L) {
  scripts_dir <- dirname(resolve_fast_overnight_path("scripts", "run_fast_multiplier_overnight.R"))
  log_dir <- canonical_fast_multiplier_logs_dir("overnight")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  steps <- list(
    list(
      label = "Fast multiplier regression tests",
      command = "Rscript",
      args = c(
        file.path(scripts_dir, "run_fast_multiplier_tests.R")
      ),
      env = c(
        "RENV_CONFIG_SANDBOX_ENABLED=FALSE",
        "RENV_CONFIG_AUTO_SNAPSHOT=FALSE"
      )
    ),
    list(
      label = "Sequential fast validation and paper experiments",
      command = "Rscript",
      args = c(
        file.path(scripts_dir, "run_all_fast_multiplier_sequential.R"),
        sprintf("--n_cores=%d", as.integer(n_cores)),
        sprintf("--validation_B=%d", as.integer(validation_B)),
        sprintf("--validation_M_outer=%d", as.integer(validation_M_outer)),
        sprintf("--derivative_mc_size=%d", as.integer(derivative_mc_size)),
        sprintf("--derivative_mc_seed=%d", as.integer(derivative_mc_seed)),
        sprintf("--output_root=%s", output_root)
      ),
      env = c(
        "RENV_CONFIG_SANDBOX_ENABLED=FALSE",
        "RENV_CONFIG_AUTO_SNAPSHOT=FALSE"
      )
    )
  )

  progress_bar <- utils::txtProgressBar(min = 0, max = length(steps), style = 3)
  on.exit(close(progress_bar), add = TRUE)

  results <- vector("list", length(steps))
  for (i in seq_along(steps)) {
    results[[i]] <- run_logged_step(
      step_index = i,
      n_steps = length(steps),
      label = steps[[i]]$label,
      command = steps[[i]]$command,
      args = steps[[i]]$args,
      env = steps[[i]]$env %||% character(0),
      log_dir = log_dir
    )
    utils::setTxtProgressBar(progress_bar, i)
  }

  invisible(results)
}

if (sys.nframe() == 0L) {
  args <- parse_fast_overnight_args(commandArgs(trailingOnly = TRUE))
  run_fast_multiplier_overnight(
    output_root = args$output_root %||% file.path("output"),
    n_cores = as.integer(args$n_cores %||% 12L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L),
    derivative_mc_seed = as.integer(args$derivative_mc_seed %||% 20260613L),
    validation_B = as.integer(args$validation_B %||% 1000L),
    validation_M_outer = as.integer(args$validation_M_outer %||% 5L)
  )
}
