resolve_rotmix_progressive_path <- function(...) {
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

source(resolve_rotmix_progressive_path("bootstrap", "calibration_study.R"))

parse_named_args_rotmix_progressive <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
  output <- vector("list", length(args))
  names(output) <- rep("", length(args))

  for (i in seq_along(args)) {
    arg <- substring(args[[i]], 3L)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[i]] <- value
    names(output)[[i]] <- key
  }

  output
}

default_rotmix_composite_progressive_scenarios <- function() {
  list(
    default_rotational_beta_mixture2_composite_calibration_scenarios()[[1L]],
    default_rotational_logitnormal_mixture2_composite_calibration_scenarios()[[1L]]
  )
}

run_rotmix_composite_progressive_B <- function(output_root = file.path(
                                                 "output",
                                                 "bootstrap_calibration",
                                                 "rotational_mixtures_composite_progressive_B"
                                               ),
                                               B_values = c(25L, 50L, 100L, 200L, 300L, 400L, 500L),
                                               B_target = 500L,
                                               n_values = 50L,
                                               M_outer = 500L,
                                               n_cores_outer = 12L,
                                               statistics = c("ks", "cvm"),
                                               alphas = c(0.01, 0.05, 0.10),
                                               seed = 20260601L,
                                               max_runtime_seconds = 5400,
                                               show_progress = TRUE,
                                               verbose = TRUE) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  scenarios <- default_rotmix_composite_progressive_scenarios()
  B_values <- sort(unique(as.integer(B_values)))
  B_values <- B_values[B_values > 0L]
  if (!(as.integer(B_target) %in% B_values)) {
    B_values <- sort(unique(c(B_values, as.integer(B_target))))
  }

  manifest_rows <- list()
  estimate_history <- list()

  for (i in seq_along(B_values)) {
    B_i <- as.integer(B_values[[i]])
    run_dir <- file.path(output_root, sprintf("B_%03d", B_i))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

    start_time <- proc.time()[["elapsed"]]
    result_i <- run_bootstrap_calibration_study(
      scenarios = scenarios,
      n_values = as.integer(n_values),
      M_outer = as.integer(M_outer),
      B = B_i,
      alpha_nominal = 0.05,
      alphas = as.numeric(alphas),
      statistics = statistics,
      n_cores_outer = as.integer(n_cores_outer),
      seed = as.integer(seed) + i,
      output_dir = run_dir,
      show_progress = isTRUE(show_progress),
      verbose = isTRUE(verbose)
    )
    elapsed_i <- proc.time()[["elapsed"]] - start_time
    projected_target <- elapsed_i * as.numeric(B_target) / as.numeric(B_i)
    estimate_row <- data.frame(
      B = B_i,
      elapsed_seconds = elapsed_i,
      projected_seconds_at_target = projected_target,
      projected_minutes_at_target = projected_target / 60,
      stop_threshold_seconds = as.numeric(max_runtime_seconds),
      stringsAsFactors = FALSE
    )
    estimate_history[[length(estimate_history) + 1L]] <- estimate_row
    utils::write.csv(
      do.call(rbind, estimate_history),
      file = file.path(output_root, "progressive_runtime_estimates.csv"),
      row.names = FALSE
    )

    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      B = B_i,
      output_dir = run_dir,
      raw_csv = result_i$raw_csv,
      summary_csv = result_i$summary_csv,
      elapsed_seconds = elapsed_i,
      projected_seconds_at_target = projected_target,
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      do.call(rbind, manifest_rows),
      file = file.path(output_root, "progressive_manifest.csv"),
      row.names = FALSE
    )

    if (B_i < as.integer(B_target) && projected_target > as.numeric(max_runtime_seconds)) {
      return(list(
        status = "stopped_early",
        stopping_B = B_i,
        elapsed_seconds = elapsed_i,
        projected_seconds_at_target = projected_target,
        output_root = output_root,
        manifest = do.call(rbind, manifest_rows),
        estimates = do.call(rbind, estimate_history)
      ))
    }
  }

  list(
    status = "completed",
    output_root = output_root,
    manifest = do.call(rbind, manifest_rows),
    estimates = do.call(rbind, estimate_history)
  )
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_rotmix_progressive(commandArgs(trailingOnly = TRUE))

  result <- run_rotmix_composite_progressive_B(
    output_root = args$output_root %||% file.path(
      "output",
      "bootstrap_calibration",
      "rotational_mixtures_composite_progressive_B"
    ),
    B_values = if (!is.null(args$B_values)) {
      as.integer(strsplit(args$B_values, ",", fixed = TRUE)[[1L]])
    } else {
      c(25L, 50L, 100L, 200L, 300L, 400L, 500L)
    },
    B_target = if (!is.null(args$B_target)) as.integer(args$B_target) else 500L,
    n_values = if (!is.null(args$n_values)) as.integer(strsplit(args$n_values, ",", fixed = TRUE)[[1L]]) else 50L,
    M_outer = if (!is.null(args$M_outer)) as.integer(args$M_outer) else 500L,
    n_cores_outer = if (!is.null(args$n_cores_outer)) as.integer(args$n_cores_outer) else 12L,
    statistics = if (!is.null(args$statistics)) {
      strsplit(tolower(args$statistics), ",", fixed = TRUE)[[1L]]
    } else {
      c("ks", "cvm")
    },
    alphas = if (!is.null(args$alphas)) {
      as.numeric(strsplit(args$alphas, ",", fixed = TRUE)[[1L]])
    } else {
      c(0.01, 0.05, 0.10)
    },
    seed = if (!is.null(args$seed)) as.integer(args$seed) else 20260601L,
    max_runtime_seconds = if (!is.null(args$max_runtime_seconds)) as.numeric(args$max_runtime_seconds) else 5400,
    show_progress = if (!is.null(args$show_progress)) isTRUE(as.logical(args$show_progress)) else TRUE,
    verbose = if (!is.null(args$verbose)) isTRUE(as.logical(args$verbose)) else TRUE
  )

  saveRDS(result, file = file.path(result$output_root, "progressive_result.rds"))
  message(sprintf("Progressive status: %s", result$status))
  if (identical(result$status, "stopped_early")) {
    message(sprintf(
      "Stopped at B=%d because projected runtime at B=%d is %.1f minutes.",
      result$stopping_B,
      if (!is.null(args$B_target)) as.integer(args$B_target) else 500L,
      result$projected_seconds_at_target / 60
    ))
  }
}
