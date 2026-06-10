resolve_comets_small_circle_path <- function(...) {
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

model_specs_script_path_sc_circle <- resolve_comets_small_circle_path("bootstrap", "model_specs.R")
small_circle_model_spec_script_path <- resolve_comets_small_circle_path("bootstrap", "small_circle_model_spec.R")
multiplier_bootstrap_script_path_sc_circle <- resolve_comets_small_circle_path("bootstrap", "multiplier_bootstrap.R")
utils_script_path_sc_circle <- resolve_comets_small_circle_path("utils.R")
comets_utils_script_path_sc_circle <- resolve_comets_small_circle_path(
  "real_data",
  "comets",
  "utils_comets_data.R"
)

source(model_specs_script_path_sc_circle)
source(small_circle_model_spec_script_path)
source(multiplier_bootstrap_script_path_sc_circle)
source(utils_script_path_sc_circle)
source(comets_utils_script_path_sc_circle)

make_small_circle_spec <- get("make_small_circle_spec", mode = "function")
multiplier_bootstrap_gof <- get("multiplier_bootstrap_gof", mode = "function")
generate_canonical_lattice <- get("generate_canonical_lattice", mode = "function")

timestamp_tag_small_circle_comets <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

write_lines_if_possible_small_circle_comets <- function(lines, path) {
  writeLines(as.character(lines), con = path)
  invisible(path)
}

load_existing_small_circle_benchmark_summary <- function(output_root) {
  summary_path <- file.path(output_root, "benchmark_summary.csv")
  if (!file.exists(summary_path)) {
    return(NULL)
  }

  utils::read.csv(summary_path, stringsAsFactors = FALSE)
}

run_small_circle_benchmark_stage <- function(data_matrix,
                                             statistic,
                                             B_value,
                                             n_cores,
                                             seed,
                                             distance_type,
                                             ks_t_points,
                                             control) {
  requested_cores <- min(as.integer(n_cores), 12L)

  run_once <- function(stage_cores) {
    multiplier_bootstrap_gof(
      data = data_matrix,
      spec = make_small_circle_spec(distance_type = distance_type),
      null = list(type = "composite"),
      statistics = statistic,
      ks_grid = if (identical(statistic, "ks")) make_small_circle_ks_grid_comets(B_value, ks_t_points) else NULL,
      B = B_value,
      alpha = 0.05,
      n_cores = as.integer(stage_cores),
      seed = as.integer(seed),
      keep = list(observed_process = TRUE, bootstrap_statistics = TRUE, bootstrap_thetas = TRUE),
      control = control
    )
  }

  stage_try <- tryCatch(
    run_once(requested_cores),
    error = identity
  )

  if (!inherits(stage_try, "error")) {
    return(stage_try)
  }

  stage_message <- conditionMessage(stage_try)
  socket_failure <- grepl(
    "unserialize\\(node\\$con\\)|error reading from connection|Cluster setup failed|failed to connect",
    stage_message,
    ignore.case = TRUE
  )

  if (!socket_failure || requested_cores <= 1L) {
    stop(stage_try)
  }

  message(sprintf(
    "Stage B=M=%d failed with %d cores due to cluster/socket error; retrying with 1 core.",
    as.integer(B_value),
    requested_cores
  ))
  run_once(1L)
}

load_comets_distance_profile_data_small_circle <- function() {
  load_comets_real_data(finite_normals = "none")
}

make_small_circle_ks_grid_comets <- function(M_value,
                                             ks_t_points = 250L) {
  list(
    omega_grid = generate_canonical_lattice(as.integer(M_value), dim = 3),
    t_grid = seq(1e-8, pi - 1e-8, length.out = as.integer(ks_t_points))
  )
}

summarize_small_circle_model_result <- function(result,
                                                dataset_label,
                                                statistic,
                                                M_value) {
  theta_hat <- result$observed$theta_hat
  inference <- result$inference[[statistic]]
  theta_star <- result$bootstrap$theta_star

  kappa_star_range <- if (is.null(theta_star) || length(theta_star) == 0L) {
    c(NA_real_, NA_real_)
  } else {
    range(vapply(theta_star, `[[`, numeric(1), "kappa"))
  }
  nu_star_range <- if (is.null(theta_star) || length(theta_star) == 0L) {
    c(NA_real_, NA_real_)
  } else {
    range(vapply(theta_star, `[[`, numeric(1), "nu"))
  }

  data.frame(
    dataset = dataset_label,
    model = "small_circle_composite",
    null_type = "composite",
    statistic = statistic,
    spec_name = result$diagnostics$spec_name,
    engine = result$diagnostics$engine,
    method = result$diagnostics$method,
    weighted_mle = isTRUE(result$diagnostics$weighted_mle),
    observed_statistic = inference$observed,
    critical_value = inference$critical_value,
    p_value = inference$p_value,
    reject = inference$reject,
    mu_hat_1 = theta_hat$mu[[1L]],
    mu_hat_2 = theta_hat$mu[[2L]],
    mu_hat_3 = theta_hat$mu[[3L]],
    kappa_hat = theta_hat$kappa,
    nu_hat = theta_hat$nu,
    kappa_star_min = kappa_star_range[[1L]],
    kappa_star_max = kappa_star_range[[2L]],
    nu_star_min = nu_star_range[[1L]],
    nu_star_max = nu_star_range[[2L]],
    B = result$bootstrap$B,
    M = as.integer(M_value),
    n = result$diagnostics$n,
    n_cores = result$diagnostics$n_cores,
    elapsed_seconds = result$diagnostics$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

run_comets_distance_profile_small_circle_benchmark <- function(output_root = NULL,
                                                               dataset = c("short", "long"),
                                                               B_values = c(10L, 30L, 50L, 100L, 200L, 500L, 1000L),
                                                               statistic = "ks",
                                                               n_cores = 12L,
                                                               ks_t_points = 250L,
                                                               base_seed = 20260531L,
                                                               distance_type = "geodesic",
                                                               control = list(
                                                                 small_circle_profile_method = "legendre",
                                                                 small_circle_L_max = 200L,
                                                                 small_circle_quad_n = 400L,
                                                                 small_circle_tol = 1e-10,
                                                                 small_circle_optim_control = list(maxit = 300L, reltol = 1e-9)
                                                               )) {
  dataset <- match.arg(dataset)
  if (is.null(output_root)) {
    output_root <- file.path(
      "output",
      "comets",
      "small_circle",
      paste0("run_", timestamp_tag_small_circle_comets(), "_", dataset, "_", statistic, "_benchmark")
    )
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  B_values <- unique(as.integer(B_values))
  B_values <- B_values[is.finite(B_values) & B_values > 0L]
  comets_data <- load_comets_distance_profile_data_small_circle()
  data_matrix <- if (identical(dataset, "short")) comets_data$short$normal else comets_data$long$normal
  dataset_label <- if (identical(dataset, "short")) "short_period" else "long_period"

  write_comets_dataset_summary(
    path = file.path(output_root, "dataset_summary.csv"),
    dataset_label = dataset_label,
    data_matrix = data_matrix
  )
  write_lines_if_possible_small_circle_comets(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))

  existing_summary <- load_existing_small_circle_benchmark_summary(output_root)
  summary_rows <- list()
  for (i in seq_along(B_values)) {
    B_value <- B_values[[i]]
    stage_file <- file.path(output_root, sprintf("stage_%02d_M%d_B%d.rds", i, B_value, B_value))

    if (file.exists(stage_file)) {
      existing_row <- NULL
      if (!is.null(existing_summary)) {
        row_idx <- which(
          existing_summary$dataset == dataset_label &
            existing_summary$statistic == statistic &
            existing_summary$B == B_value &
            existing_summary$M == B_value
        )
        if (length(row_idx) >= 1L) {
          existing_row <- existing_summary[row_idx[[1L]], , drop = FALSE]
        }
      }

      if (is.null(existing_row)) {
        result_existing <- readRDS(stage_file)
        existing_row <- summarize_small_circle_model_result(
          result = result_existing,
          dataset_label = dataset_label,
          statistic = statistic,
          M_value = B_value
        )
      }

      summary_rows[[length(summary_rows) + 1L]] <- existing_row
      next
    }

    result <- run_small_circle_benchmark_stage(
      data_matrix = data_matrix,
      statistic = statistic,
      B_value = B_value,
      n_cores = n_cores,
      seed = as.integer(base_seed + i),
      distance_type = distance_type,
      ks_t_points = ks_t_points,
      control = control
    )

    summary_rows[[i]] <- summarize_small_circle_model_result(
      result = result,
      dataset_label = dataset_label,
      statistic = statistic,
      M_value = B_value
    )
    utils::write.csv(do.call(rbind, summary_rows), file = file.path(output_root, "benchmark_summary.csv"), row.names = FALSE)
    saveRDS(result, file = file.path(output_root, sprintf("stage_%02d_M%d_B%d.rds", i, B_value, B_value)))
  }

  invisible(output_root)
}

parse_named_args_small_circle_comets <- function(args) {
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
  args <- parse_named_args_small_circle_comets(commandArgs(trailingOnly = TRUE))
  run_comets_distance_profile_small_circle_benchmark(
    output_root = args$output_root %||% NULL,
    dataset = args$dataset %||% "short",
    B_values = if (!is.null(args$B_values)) as.integer(strsplit(args$B_values, ",", fixed = TRUE)[[1L]]) else c(10L, 30L, 50L, 100L, 200L, 500L, 1000L),
    statistic = args$statistic %||% "ks",
    n_cores = if (!is.null(args$n_cores)) as.integer(args$n_cores) else 12L,
    ks_t_points = if (!is.null(args$ks_t_points)) as.integer(args$ks_t_points) else 250L,
    base_seed = if (!is.null(args$seed)) as.integer(args$seed) else 20260531L,
    distance_type = args$distance_type %||% "geodesic"
  )
}
