resolve_rotational_mixture_calibration_path <- function(...) {
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

calibration_study_script_rotational_mixtures <- resolve_rotational_mixture_calibration_path(
  "bootstrap",
  "calibration_study.R"
)
source(calibration_study_script_rotational_mixtures)

rotmix_expected_raw_rows <- function(M_outer,
                                     n_values,
                                     statistics) {
  as.integer(M_outer) * length(as.integer(n_values)) * length(statistics)
}

rotmix_expected_summary_rows <- function(n_values,
                                         statistics,
                                         alphas) {
  length(as.integer(n_values)) * length(statistics) * length(as.numeric(alphas))
}

rotmix_is_scenario_output_complete <- function(output_dir,
                                               expected_raw_rows,
                                               expected_summary_rows) {
  raw_csv <- file.path(output_dir, "bootstrap_calibration_raw.csv")
  summary_csv <- file.path(output_dir, "bootstrap_calibration_summary.csv")

  if (!file.exists(raw_csv) || !file.exists(summary_csv)) {
    return(FALSE)
  }

  raw_df <- try(utils::read.csv(raw_csv, stringsAsFactors = FALSE), silent = TRUE)
  summary_df <- try(utils::read.csv(summary_csv, stringsAsFactors = FALSE), silent = TRUE)

  if (inherits(raw_df, "try-error") || inherits(summary_df, "try-error")) {
    return(FALSE)
  }

  nrow(raw_df) >= expected_raw_rows && nrow(summary_df) >= expected_summary_rows
}

rotmix_run_single_scenario_with_checkpoint <- function(scenario,
                                                       output_dir,
                                                       n_values,
                                                       M_outer,
                                                       B,
                                                       statistics,
                                                       alphas,
                                                       n_cores_outer,
                                                       seed,
                                                       show_progress = TRUE,
                                                       verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  expected_raw <- rotmix_expected_raw_rows(
    M_outer = M_outer,
    n_values = n_values,
    statistics = statistics
  )
  expected_summary <- rotmix_expected_summary_rows(
    n_values = n_values,
    statistics = statistics,
    alphas = alphas
  )

  if (rotmix_is_scenario_output_complete(
    output_dir = output_dir,
    expected_raw_rows = expected_raw,
    expected_summary_rows = expected_summary
  )) {
    message(sprintf("[SKIP] Scenario '%s' already complete at %s", scenario$id, output_dir))
    return(list(
      scenario_id = scenario$id,
      output_dir = output_dir,
      skipped = TRUE,
      raw_csv = file.path(output_dir, "bootstrap_calibration_raw.csv"),
      summary_csv = file.path(output_dir, "bootstrap_calibration_summary.csv")
    ))
  }

  message(sprintf("[RUN ] Scenario '%s' -> %s", scenario$id, output_dir))
  result <- run_bootstrap_calibration_study(
    scenarios = list(scenario),
    n_values = as.integer(n_values),
    M_outer = as.integer(M_outer),
    B = as.integer(B),
    alpha_nominal = 0.05,
    alphas = as.numeric(alphas),
    statistics = statistics,
    n_cores_outer = as.integer(n_cores_outer),
    seed = as.integer(seed),
    output_dir = output_dir,
    show_progress = isTRUE(show_progress),
    verbose = isTRUE(verbose)
  )

  list(
    scenario_id = scenario$id,
    output_dir = output_dir,
    skipped = FALSE,
    raw_csv = result$raw_csv,
    summary_csv = result$summary_csv
  )
}

rotmix_consolidate_sequential_outputs <- function(run_rows,
                                                  consolidated_dir) {
  dir.create(consolidated_dir, recursive = TRUE, showWarnings = FALSE)

  summary_all <- do.call(rbind, lapply(run_rows$summary_csv, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$source_summary_csv <- path
    x
  }))
  raw_all <- do.call(rbind, lapply(run_rows$raw_csv, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$source_raw_csv <- path
    x
  }))

  summary_out <- file.path(consolidated_dir, "bootstrap_calibration_summary_all_scenarios.csv")
  raw_out <- file.path(consolidated_dir, "bootstrap_calibration_raw_all_scenarios.csv")
  manifest_out <- file.path(consolidated_dir, "sequential_manifest.csv")

  utils::write.csv(summary_all, summary_out, row.names = FALSE)
  utils::write.csv(raw_all, raw_out, row.names = FALSE)
  utils::write.csv(run_rows, manifest_out, row.names = FALSE)

  list(
    summary_csv = summary_out,
    raw_csv = raw_out,
    manifest_csv = manifest_out
  )
}

default_rotational_mixture_selected_scenarios <- function() {
  list(
    default_rotational_beta_mixture2_simple_calibration_scenarios()[[1L]],
    default_rotational_logitnormal_mixture2_simple_calibration_scenarios()[[1L]],
    default_rotational_beta_mixture2_composite_calibration_scenarios()[[1L]],
    default_rotational_logitnormal_mixture2_composite_calibration_scenarios()[[1L]]
  )
}

run_rotational_mixture_calibration_b500_sequential <- function(
  output_root = file.path(
    "output",
    "bootstrap_calibration",
    "rotational_mixtures_simple_composite_M500_B500_bootlocal"
  ),
  n_values = c(50L, 100L, 200L),
  M_outer = 500L,
  B = 500L,
  n_cores_outer = 12L,
  statistics = c("ks", "cvm"),
  alphas = c(0.01, 0.05, 0.10),
  seed = 20260601L,
  show_progress = TRUE,
  verbose = TRUE
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  scenarios <- default_rotational_mixture_selected_scenarios()
  run_rows <- list()

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    scenario_seed <- as.integer(seed) + i
    scenario_dir <- file.path(output_root, sprintf("%02d_%s", i, scenario$id))

    run_info <- rotmix_run_single_scenario_with_checkpoint(
      scenario = scenario,
      output_dir = scenario_dir,
      n_values = as.integer(n_values),
      M_outer = as.integer(M_outer),
      B = as.integer(B),
      statistics = statistics,
      alphas = as.numeric(alphas),
      n_cores_outer = as.integer(n_cores_outer),
      seed = scenario_seed,
      show_progress = isTRUE(show_progress),
      verbose = isTRUE(verbose)
    )

    run_rows[[length(run_rows) + 1L]] <- data.frame(
      order_id = i,
      scenario_id = run_info$scenario_id,
      output_dir = run_info$output_dir,
      skipped = run_info$skipped,
      raw_csv = run_info$raw_csv,
      summary_csv = run_info$summary_csv,
      stringsAsFactors = FALSE
    )
  }

  manifest_df <- do.call(rbind, run_rows)
  consolidated <- rotmix_consolidate_sequential_outputs(
    run_rows = manifest_df,
    consolidated_dir = file.path(output_root, "_consolidated")
  )

  result <- list(
    output_root = output_root,
    manifest = manifest_df,
    consolidated = consolidated,
    config = list(
      n_values = as.integer(n_values),
      M_outer = as.integer(M_outer),
      B = as.integer(B),
      n_cores_outer = as.integer(n_cores_outer),
      statistics = statistics,
      alphas = as.numeric(alphas),
      seed = as.integer(seed)
    )
  )

  saveRDS(result, file = file.path(output_root, "run_result.rds"))
  result
}

parse_named_args_rotmix_calibration <- function(args) {
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

if (sys.nframe() == 0L) {
  args <- parse_named_args_rotmix_calibration(commandArgs(trailingOnly = TRUE))

  result <- run_rotational_mixture_calibration_b500_sequential(
    output_root = args$output_root %||% file.path(
      "output",
      "bootstrap_calibration",
      "rotational_mixtures_simple_composite_M500_B500_bootlocal"
    ),
    n_values = if (!is.null(args$n_values)) {
      as.integer(strsplit(args$n_values, ",", fixed = TRUE)[[1L]])
    } else {
      c(50L, 100L, 200L)
    },
    M_outer = if (!is.null(args$M_outer)) as.integer(args$M_outer) else 500L,
    B = if (!is.null(args$B)) as.integer(args$B) else 500L,
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
    show_progress = if (!is.null(args$show_progress)) isTRUE(as.logical(args$show_progress)) else TRUE,
    verbose = if (!is.null(args$verbose)) isTRUE(as.logical(args$verbose)) else TRUE
  )

  message(sprintf("Consolidated summary: %s", result$consolidated$summary_csv))
}
