resolve_comets_jp_long_path <- function(...) {
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

short_benchmark_script <- resolve_comets_jp_long_path(
  "scripts",
  "run_comets_distance_profile_jp_short_benchmark.R"
)
source(short_benchmark_script)

load_comets_distance_profile_data_jp_long <- function() {
  if (!requireNamespace("sphunif", quietly = TRUE)) {
    stop("Package `sphunif` is required for the comet analyses.")
  }

  data("comets", package = "sphunif")
  comets$normal <- cbind(
    sin(comets$i) * sin(comets$om),
    -sin(comets$i) * cos(comets$om),
    cos(comets$i)
  )

  valid_rows <-
    !is.na(comets$class) &
    is.finite(comets$per_y) &
    is.finite(comets$i) &
    is.finite(comets$om) &
    !is.na(comets$frag)

  dropped_incomplete <- sum(!valid_rows)
  if (dropped_incomplete > 0L) {
    warning(
      sprintf(
        "Dropping %d comet rows with incomplete orbital elements before long-period filtering.",
        dropped_incomplete
      ),
      call. = FALSE
    )
  }

  long_selector <-
    valid_rows &
    !(comets$class %in% c("HYP", "PAR")) &
    comets$per_y >= 200 &
    !comets$frag
  comets_long <- comets[long_selector, , drop = FALSE]

  normal_matrix <- as.matrix(comets_long$normal)
  finite_rows <- apply(normal_matrix, 1L, function(r) all(is.finite(r)))
  dropped_nonfinite <- sum(!finite_rows)
  if (dropped_nonfinite > 0L) {
    warning(
      sprintf(
        "Dropping %d long-period comet rows with non-finite normal coordinates after filtering.",
        dropped_nonfinite
      ),
      call. = FALSE
    )
    comets_long <- comets_long[finite_rows, , drop = FALSE]
  }

  list(
    raw = comets,
    long = comets_long
  )
}

run_comets_distance_profile_jp_long_benchmark <- function(output_root = NULL,
                                                          B_values = c(10L, 30L, 50L, 100L, 200L, 500L, 1000L),
                                                          statistic = "ks",
                                                          n_cores = 12L,
                                                          ks_t_points = 250L,
                                                          base_seed = 20260529L,
                                                          distance_type = "geodesic",
                                                          control = list(
                                                            jp_profile_method = "tabulated",
                                                            jp_profile_n_u = 1025L,
                                                            jp_profile_n_delta = 257L,
                                                            jp_vmf_switch_abs_kappa_psi = 1e-3,
                                                            jp_mle_max_abs_kappa_psi = 6
                                                          )) {
  if (is.null(output_root)) {
    output_root <- file.path(
      "output",
      "comets_distance_profile_jp",
      paste0("run_", timestamp_tag_jp_comets(), "_long_", statistic, "_benchmark")
    )
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  B_values <- unique(as.integer(B_values))
  B_values <- B_values[is.finite(B_values) & B_values > 0L]
  if (length(B_values) == 0L) {
    stop("`B_values` must contain at least one positive integer.")
  }

  comets_data <- load_comets_distance_profile_data_jp_long()
  dataset_summary <- data.frame(
    dataset = "long_period",
    n = nrow(comets_data$long),
    ambient_dim = ncol(comets_data$long$normal),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    dataset_summary,
    file = file.path(output_root, "dataset_summary.csv"),
    row.names = FALSE
  )

  write_lines_if_possible_jp_comets(
    capture.output(sessionInfo()),
    file.path(output_root, "sessionInfo.txt")
  )

  saveRDS(
    list(
      output_root = output_root,
      benchmark_B_values = B_values,
      benchmark_M_values = B_values,
      statistic = statistic,
      n_cores = as.integer(n_cores),
      ks_t_points = as.integer(ks_t_points),
      base_seed = as.integer(base_seed),
      distance_type = distance_type,
      dataset = "long_period",
      model = "jp_composite"
    ),
    file = file.path(output_root, "run_config.rds")
  )

  manifest_rows <- list()
  stage_results <- list()

  for (i in seq_along(B_values)) {
    B_value <- B_values[[i]]
    stage_id <- sprintf("%02d", i)
    stage_name <- sprintf("long_%s_M%d_B%d", statistic, B_value, B_value)
    stage_dir <- file.path(output_root, paste0(stage_id, "_", stage_name))
    stage_start <- Sys.time()

    message(sprintf(
      "[Long-period JP %s] %d/%d: M = B = %d with %d cores",
      toupper(statistic),
      i,
      length(B_values),
      B_value,
      as.integer(n_cores)
    ))

    stage_results[[stage_name]] <- run_jp_short_comet_stage(
      data_matrix = comets_data$long$normal,
      dataset_label = sprintf("Long-period JP %s", toupper(statistic)),
      statistic = statistic,
      stage_dir = stage_dir,
      B = B_value,
      M_value = B_value,
      n_cores = n_cores,
      seed = base_seed + i,
      distance_type = distance_type,
      ks_t_points = ks_t_points,
      control = control
    )

    manifest_rows <- append_manifest_row_jp_comets(
      manifest_rows = manifest_rows,
      stage_id = stage_id,
      stage_label = sprintf("Long-period JP %s M=B=%d", toupper(statistic), B_value),
      status = "completed",
      stage_dir = stage_dir,
      started_at = stage_start,
      finished_at = Sys.time()
    )

    utils::write.csv(
      do.call(rbind, manifest_rows),
      file = file.path(output_root, "pipeline_manifest.csv"),
      row.names = FALSE
    )

    benchmark_summary <- do.call(
      rbind,
      lapply(stage_results, function(bundle) bundle$summary)
    )
    utils::write.csv(
      benchmark_summary,
      file = file.path(output_root, "benchmark_summary.csv"),
      row.names = FALSE
    )
  }

  pipeline_result <- list(
    output_root = output_root,
    dataset_summary = dataset_summary,
    stages = stage_results,
    manifest = do.call(rbind, manifest_rows),
    engine = "multiplier_bootstrap_gof",
    method = "distance_profiles"
  )

  saveRDS(pipeline_result, file = file.path(output_root, "pipeline_result.rds"))
  pipeline_result
}

parse_named_args_jp_long_comets <- function(args) {
  parse_named_args_jp_comets(args)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_jp_long_comets(commandArgs(trailingOnly = TRUE))

  B_values <- if (!is.null(args$B_values)) {
    as.integer(strsplit(args$B_values, ",", fixed = TRUE)[[1L]])
  } else {
    c(10L, 30L, 50L, 100L, 200L, 500L, 1000L)
  }

  output_root <- if (!is.null(args$output_root)) args$output_root else NULL
  statistic <- if (!is.null(args$statistic)) tolower(args$statistic) else "ks"
  n_cores <- if (!is.null(args$n_cores)) as.integer(args$n_cores) else 12L
  ks_t_points <- if (!is.null(args$ks_t_points)) as.integer(args$ks_t_points) else 250L
  base_seed <- if (!is.null(args$seed)) as.integer(args$seed) else 20260529L
  distance_type <- if (!is.null(args$distance_type)) args$distance_type else "geodesic"

  run_comets_distance_profile_jp_long_benchmark(
    output_root = output_root,
    B_values = B_values,
    statistic = statistic,
    n_cores = n_cores,
    ks_t_points = ks_t_points,
    base_seed = base_seed,
    distance_type = distance_type
  )
}
