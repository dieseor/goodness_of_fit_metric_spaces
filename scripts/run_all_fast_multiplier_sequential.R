resolve_run_all_fast_path <- function(...) {
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

parse_run_all_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

find_orphan_rsock_workers <- function() {
  cmd <- paste(
    "ps -axo pid=,ppid=,command=",
    "| awk '",
    "$2 == 1 && $0 ~ /parallel:::\\.workRSOCK/ { print $1 }",
    "'"
  )
  out <- tryCatch(system(cmd, intern = TRUE), warning = function(e) character(0), error = function(e) character(0))
  pids <- suppressWarnings(as.integer(trimws(out)))
  pids[is.finite(pids) & pids > 0L]
}

kill_orphan_rsock_workers <- function(verbose = TRUE) {
  pids <- find_orphan_rsock_workers()
  if (length(pids) == 0L) {
    if (isTRUE(verbose)) message("[cleanup] no orphan RSOCK workers found")
    return(invisible(integer(0)))
  }

  if (isTRUE(verbose)) {
    message(sprintf("[cleanup] killing orphan RSOCK workers: %s", paste(pids, collapse = ", ")))
  }
  tools::pskill(pids, tools::SIGTERM)
  Sys.sleep(2)
  ps_out <- tryCatch(
    system(sprintf("ps -p %s -o pid=", paste(pids, collapse = ",")), intern = TRUE),
    warning = function(e) character(0),
    error = function(e) character(0)
  )
  still_alive <- suppressWarnings(as.integer(trimws(ps_out)))
  still_alive <- still_alive[is.finite(still_alive) & still_alive > 0L]
  if (length(still_alive) > 0L) {
    tools::pskill(still_alive, tools::SIGKILL)
  }
  invisible(pids)
}

run_task <- function(task_index,
                     n_tasks,
                     label,
                     command_args,
                     log_dir) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, sprintf("%02d_%s.log", task_index, gsub("[^A-Za-z0-9_]+", "_", tolower(label))))
  border <- paste(rep("=", 72), collapse = "")
  message(border)
  message(sprintf("[task %d/%d] %s", task_index, n_tasks, label))
  message(sprintf("[task %d/%d] log: %s", task_index, n_tasks, log_path))
  message(border)
  start <- Sys.time()
  full_cmd <- sprintf(
    "env RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_AUTO_SNAPSHOT=FALSE Rscript %s 2>&1 | tee %s",
    paste(shQuote(command_args), collapse = " "),
    shQuote(log_path)
  )
  status <- system(full_cmd)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  if (!identical(status, 0L)) {
    stop(sprintf("Task failed: %s. See %s", label, log_path))
  }
  message(sprintf("[task %d/%d] completed in %.1f seconds: %s", task_index, n_tasks, elapsed, label))
  invisible(list(log_path = log_path, elapsed_seconds = elapsed))
}

run_all_fast_multiplier_sequential <- function(n_cores = 12L,
                                               derivative_mc_size = 1000L,
                                               derivative_mc_seed = 20260613L,
                                               output_root = file.path("output", "fast_multiplier"),
                                               validation_B = 1000L,
                                               validation_M_outer = 5L,
                                               kill_orphans_first = TRUE) {
  scripts_dir <- dirname(resolve_run_all_fast_path("scripts", "run_all_fast_multiplier_sequential.R"))
  if (isTRUE(kill_orphans_first)) {
    kill_orphan_rsock_workers(verbose = TRUE)
  }

  log_dir <- file.path(output_root, "logs")
  validation_output <- file.path(
    output_root,
    sprintf("validation_all_models_12cores_B%d_M%d", as.integer(validation_B), as.integer(validation_M_outer))
  )
  paper_output <- file.path(output_root, "paper_results")

  tasks <- list(
    list(
      label = "Validation old_vs_fast all models",
      args = c(
        file.path(scripts_dir, "run_fast_multiplier_validation_all_models.R"),
        sprintf("--output_root=%s", validation_output),
        sprintf("--B=%d", as.integer(validation_B)),
        sprintf("--M_outer=%d", as.integer(validation_M_outer)),
        sprintf("--derivative_mc_size=%d", as.integer(derivative_mc_size)),
        sprintf("--base_seed=%d", 20260613L),
        sprintf("--n_cores=%d", as.integer(n_cores))
      )
    ),
    list(
      label = "Paper composite calibration and real data fast",
      args = c(
        file.path(scripts_dir, "run_fast_multiplier_paper_results.R"),
        "--mode=paper_only",
        sprintf("--output_root=%s", paper_output),
        sprintf("--n_cores=%d", as.integer(n_cores)),
        sprintf("--derivative_mc_size=%d", as.integer(derivative_mc_size)),
        sprintf("--derivative_mc_seed=%d", as.integer(derivative_mc_seed)),
        sprintf("--validation_M_outer=%d", as.integer(validation_M_outer))
      )
    )
  )

  progress_bar <- utils::txtProgressBar(min = 0, max = length(tasks), style = 3)
  on.exit(close(progress_bar), add = TRUE)

  results <- vector("list", length(tasks))
  for (i in seq_along(tasks)) {
    results[[i]] <- run_task(
      task_index = i,
      n_tasks = length(tasks),
      label = tasks[[i]]$label,
      command_args = tasks[[i]]$args,
      log_dir = log_dir
    )
    utils::setTxtProgressBar(progress_bar, i)
  }

  message("All sequential fast-multiplier tasks completed.")
  invisible(results)
}

if (sys.nframe() == 0L) {
  args <- parse_run_all_args(commandArgs(trailingOnly = TRUE))
  run_all_fast_multiplier_sequential(
    n_cores = as.integer(args$n_cores %||% 12L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L),
    derivative_mc_seed = as.integer(args$derivative_mc_seed %||% 20260613L),
    output_root = args$output_root %||% file.path("output", "fast_multiplier"),
    validation_B = as.integer(args$validation_B %||% 1000L),
    validation_M_outer = as.integer(args$validation_M_outer %||% 5L),
    kill_orphans_first = !identical(tolower(as.character(args$kill_orphans_first %||% "true")), "false")
  )
}


%RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript scripts/run_fast_multiplier_tests.R
%
%RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript scripts/run_fast_multiplier_validation_all_models.R \
%  --output_root=output/fast_multiplier/validation_all_models_12cores_B1000_M5 \
%  --B=1000 \
%  --M_outer=5 \
%  --derivative_mc_size=1000 \
%  --base_seed=20260613 \
%  --n_cores=12
%
%RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript scripts/run_fast_multiplier_paper_results.R \
%  --mode=calibration \
%  --output_root=output/fast_multiplier/paper_results \
%  --n_cores=12 \
%  --derivative_mc_size=1000 \
%  --derivative_mc_seed=20260613
%
%RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript scripts/run_logistic_gaussian_dataset_screening.R \
%  --B=1000 \
%  --n_cores=12 \
%  --bootstrap_mode=composite_multiplier \
%  --output_dir=output/fast_multiplier/paper_results/real_data/logistic_gaussian
%
%RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript -e 'source("real_data/wind/run_risoe_125m_screening_ks_cvm.R"); run_risoe_125m_screening_ks_cvm(output_dir="output/fast_multiplier/paper_results/real_data/risoe_125m", B=1000L, n_cores=12L, bootstrap_method="fast_multiplier")'
%
%RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript -e 'source("real_data/sunspots/run_sunspots_weighted_mixture_rolling_windows_gof.R"); run_sunspots_weighted_mixture_rolling_windows_gof(output_dir="output/fast_multiplier/paper_results/real_data/sunspots_weighted_mixture", statistics="ks", B=1000L, n_cores=12L, bootstrap_method="fast_multiplier")'