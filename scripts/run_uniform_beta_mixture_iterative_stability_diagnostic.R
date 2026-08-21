# Iterative endpoint-removal stability diagnostic for the long-period
# uniform-beta mixture.  This script does not alter the bootstrap engine.

resolve_iterative_diagnostic_path <- function(...) {
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

source(resolve_iterative_diagnostic_path("utils.R"))
source(resolve_iterative_diagnostic_path("bootstrap", "model_specs.R"))
source(resolve_iterative_diagnostic_path("bootstrap", "multiplier_bootstrap.R"))
source(resolve_iterative_diagnostic_path(
  "real_data", "comets", "utils_comets_data.R"
))
source(resolve_iterative_diagnostic_path("scripts", "path_helpers.R"))
source(resolve_iterative_diagnostic_path(
  "bootstrap", "uniform_beta_mixture_model_spec.R"
))

parse_named_args_iterative_diagnostic <- function(args) {
  output <- list()

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (!startsWith(arg, "--")) {
      i <- i + 1L
      next
    }

    token <- substring(arg, 3L)
    pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]

    if (length(pieces) > 1L) {
      value <- paste(pieces[-1L], collapse = "=")
    } else if (
      i < length(args) &&
        !startsWith(args[[i + 1L]], "--")
    ) {
      value <- args[[i + 1L]]
      i <- i + 1L
    } else {
      value <- "TRUE"
    }

    output[[key]] <- value
    i <- i + 1L
  }

  output
}

args <- parse_named_args_iterative_diagnostic(
  commandArgs(trailingOnly = TRUE)
)

get_arg_iterative <- function(name, default) {
  if (is.null(args[[name]])) default else args[[name]]
}

output_root_iterative <- get_arg_iterative(
  "output_root",
  file.path(
    "real_data", "comets", "mixture",
    "uniform_beta_mixture_iterative_stability_long_period"
  )
)

max_iterations_iterative <- as.integer(get_arg_iterative(
  "max_iterations",
  20L
))

if (!is.finite(max_iterations_iterative) || max_iterations_iterative < 1L) {
  stop("`max_iterations` must be a positive integer.")
}

dir.create(
  output_root_iterative,
  recursive = TRUE,
  showWarnings = FALSE
)

# Same model-fitting controls as the original Uniform-Beta experiments.
control_uniform_iterative <- list(
  uniform_beta_mixture_profile_method = "legendre",
  uniform_beta_mixture_quad_n = 100L,
  uniform_beta_mixture_optim_control = list(
    maxit = 350L,
    reltol = 1e-9
  ),
  uniform_beta_mixture_shape_lower = NULL,
  uniform_beta_mixture_shape_upper = NULL,

  # This is used only to force the final diagnostic bootstrap through
  # the existing fast branch. It does not alter the model fit.
  uniform_beta_mixture_fast_shape_regular_eps = -1,
  fast_bootstrap_chunk_size = NULL,
  progress_bar = TRUE
)

# Same final fast-bootstrap configuration as the forced-fast run.
bootstrap_B_iterative <- 5000L
bootstrap_n_cores_iterative <- 3L
bootstrap_seed_ks_iterative <- 20261603L
bootstrap_seed_cvm_iterative <- 20262603L
bootstrap_method_iterative <- "fast_multiplier"
distance_type_iterative <- "geodesic"

comets_data_iterative <- load_comets_real_data(
  finite_normals = "long"
)

x_current <- as.matrix(comets_data_iterative$long$normal)
current_original_indices <- seq_len(nrow(x_current))

spec_iterative <- make_uniform_beta_mixture_spec(
  distance_type = distance_type_iterative
)

compute_y_iterative <- function(x, theta) {
  z <- pmin(
    pmax(as.numeric(x %*% theta$mu), -1),
    1
  )
  (1 + z) / 2
}

prepare_fast_diagnostics_iterative <- function(x, theta) {
  ks_grid <- make_sample_unique_distance_ks_grid()

  ks_prep <- prepare_ks_observed_data(
    data = x,
    spec = spec_iterative,
    theta_hat = theta,
    ks_grid = ks_grid,
    control = control_uniform_iterative,
    light = TRUE,
    share_cvm_statistic = TRUE
  )

  cvm_prep <- prepare_cvm_observed_data_from_sample_ks(
    data = x,
    spec = spec_iterative,
    theta_hat = theta,
    ks_prep = ks_prep,
    control = control_uniform_iterative
  )

  fast_prep <- spec_fast_multiplier_prepare(
    spec = spec_iterative,
    data = x,
    theta_hat = theta,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control_uniform_iterative
  )

  if (isTRUE(fast_prep$fallback_to_reestimated)) {
    stop("Fast diagnostic preparation unexpectedly requested fallback.")
  }

  fast_prep
}

extract_vhat_diagnostics_iterative <- function(fast_prep) {
  V_raw <- as.matrix(fast_prep$Vhat)
  V_sym <- (V_raw + t(V_raw)) / 2

  eigenvalues <- eigen(
    V_sym,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  V_inverse <- tryCatch(
    fast_multiplier_solve_vhat(
      V_raw,
      diag(nrow(V_raw)),
      label = "iterative Vhat diagnostic"
    ),
    error = function(e) NULL
  )

  spectral_norm_V <- norm(V_raw, type = "2")
  spectral_norm_V_inverse <- if (is.null(V_inverse)) {
    NA_real_
  } else {
    norm(V_inverse, type = "2")
  }

  condition_number <- if (is.null(V_inverse)) {
    Inf
  } else {
    spectral_norm_V * spectral_norm_V_inverse
  }

  vhat_diagnostics <- fast_prep$vhat_diagnostics
  if (is.null(vhat_diagnostics)) {
    vhat_diagnostics <- list()
  }

  list(
    V_raw = V_raw,
    V_inverse = V_inverse,
    eigenvalues = eigenvalues,
    min_abs_eigenvalue = min(abs(eigenvalues)),
    max_abs_eigenvalue = max(abs(eigenvalues)),
    spectral_norm_V = spectral_norm_V,
    spectral_norm_V_inverse = spectral_norm_V_inverse,
    condition_number = condition_number,
    rcond = as.numeric(
      vhat_diagnostics$Vhat_rcond %||% NA_real_
    ),
    fast_condition_number = as.numeric(
      vhat_diagnostics$Vhat_condition_number %||% NA_real_
    )
  )
}

save_iteration_vhat_iterative <- function(
    iteration,
    vhat_diagnostics) {
  iteration_tag <- sprintf("iteration_%03d", iteration)

  utils::write.csv(
    vhat_diagnostics$V_raw,
    file = file.path(
      output_root_iterative,
      paste0("Vhat_", iteration_tag, ".csv")
    ),
    row.names = TRUE
  )

  utils::write.csv(
    data.frame(
      iteration = iteration,
      eigen_index = seq_along(vhat_diagnostics$eigenvalues),
      eigenvalue = vhat_diagnostics$eigenvalues,
      abs_eigenvalue = abs(vhat_diagnostics$eigenvalues)
    ),
    file = file.path(
      output_root_iterative,
      paste0("eigenvalues_", iteration_tag, ".csv")
    ),
    row.names = FALSE
  )
}

iteration_rows_iterative <- list()
removed_history_iterative <- list()
converged_iterative <- FALSE
stop_reason_iterative <- NA_character_
final_theta_iterative <- NULL
final_x_iterative <- NULL
final_original_indices_iterative <- NULL
final_y_iterative <- NULL

for (iteration in 0:max_iterations_iterative) {
  theta_current <- fit_uniform_beta_mixture_theta(
    data = x_current,
    null = list(type = "composite"),
    control = control_uniform_iterative
  )

  y_current <- compute_y_iterative(
    x = x_current,
    theta = theta_current
  )

  near_zero_current <- y_current < 1e-5
  near_one_current <- y_current > 1 - 1e-5
  violating_current <- near_zero_current | near_one_current

  violating_positions <- which(violating_current)
  removed_original_indices <- current_original_indices[violating_positions]
  removed_y_values <- y_current[violating_positions]
  removed_endpoints <- ifelse(
    removed_y_values < 1e-5,
    "zero",
    "one"
  )

  fast_prep_current <- prepare_fast_diagnostics_iterative(
    x = x_current,
    theta = theta_current
  )

  vhat_current <- extract_vhat_diagnostics_iterative(
    fast_prep_current
  )

  iteration_rows_iterative[[length(iteration_rows_iterative) + 1L]] <- data.frame(
    iteration = iteration,
    n = nrow(x_current),
    removed_original_indices = paste(
      removed_original_indices,
      collapse = ","
    ),
    n_removed = length(removed_original_indices),
    n_near_zero = sum(near_zero_current),
    prop_near_zero = mean(near_zero_current),
    n_near_one = sum(near_one_current),
    prop_near_one = mean(near_one_current),
    min_y = min(y_current),
    max_y = max(y_current),
    mu_hat_1 = theta_current$mu[1],
    mu_hat_2 = theta_current$mu[2],
    mu_hat_3 = theta_current$mu[3],
    w_hat = theta_current$weight_uniform,
    alpha_hat = theta_current$alpha,
    beta_hat = theta_current$beta,
    min_abs_eigenvalue = vhat_current$min_abs_eigenvalue,
    max_abs_eigenvalue = vhat_current$max_abs_eigenvalue,
    spectral_norm_V = vhat_current$spectral_norm_V,
    spectral_norm_V_inverse = vhat_current$spectral_norm_V_inverse,
    condition_number_V = vhat_current$condition_number,
    rcond = vhat_current$rcond,
    fast_reported_condition_number = vhat_current$fast_condition_number,
    stringsAsFactors = FALSE
  )

  removed_history_iterative[[length(removed_history_iterative) + 1L]] <- if (
    length(removed_original_indices) == 0L
  ) {
    data.frame(
      iteration = integer(0),
      original_index = integer(0),
      y_at_removal = numeric(0),
      endpoint = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      iteration = rep.int(iteration, length(removed_original_indices)),
      original_index = removed_original_indices,
      y_at_removal = removed_y_values,
      endpoint = removed_endpoints,
      stringsAsFactors = FALSE
    )
  }

  save_iteration_vhat_iterative(
    iteration = iteration,
    vhat_diagnostics = vhat_current
  )

  utils::write.csv(
    data.frame(
      iteration = iteration,
      original_index = current_original_indices,
      y_i = y_current,
      near_zero = near_zero_current,
      near_one = near_one_current,
      stringsAsFactors = FALSE
    ),
    file = file.path(
      output_root_iterative,
      sprintf("y_values_iteration_%03d.csv", iteration)
    ),
    row.names = FALSE
  )

  cat(
    sprintf(
      "Iteration %d: n=%d, removed=%d, min(y)=%.17g, max(y)=%.17g\n",
      iteration,
      nrow(x_current),
      length(removed_original_indices),
      min(y_current),
      max(y_current)
    )
  )

  final_theta_iterative <- theta_current
  final_x_iterative <- x_current
  final_original_indices_iterative <- current_original_indices
  final_y_iterative <- y_current

  if (length(removed_original_indices) == 0L) {
    converged_iterative <- TRUE
    stop_reason_iterative <- "no_endpoint_violations"
    break
  }

  if (iteration >= max_iterations_iterative) {
    converged_iterative <- FALSE
    stop_reason_iterative <- "maximum_iterations_reached"
    break
  }

  x_current <- x_current[!violating_current, , drop = FALSE]
  current_original_indices <- current_original_indices[!violating_current]
}

iteration_summary_iterative <- do.call(
  rbind,
  iteration_rows_iterative
)

removed_history_nonempty <- removed_history_iterative[
  vapply(
    removed_history_iterative,
    nrow,
    integer(1)
  ) > 0L
]

removed_history_df_iterative <- if (length(removed_history_nonempty) == 0L) {
  data.frame(
    iteration = integer(0),
    original_index = integer(0),
    y_at_removal = numeric(0),
    endpoint = character(0),
    stringsAsFactors = FALSE
  )
} else {
  do.call(rbind, removed_history_nonempty)
}

utils::write.csv(
  iteration_summary_iterative,
  file = file.path(
    output_root_iterative,
    "iterative_stability_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  removed_history_df_iterative,
  file = file.path(
    output_root_iterative,
    "removed_observations_history.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  data.frame(
    original_index = final_original_indices_iterative,
    y_i = final_y_iterative,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_root_iterative,
    "final_dataset_y_values.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  data.frame(
    convergence = converged_iterative,
    stop_reason = stop_reason_iterative,
    final_iteration = tail(
      iteration_summary_iterative$iteration,
      1L
    ),
    final_n = nrow(final_x_iterative),
    final_min_y = min(final_y_iterative),
    final_max_y = max(final_y_iterative),
    final_mu_1 = final_theta_iterative$mu[1],
    final_mu_2 = final_theta_iterative$mu[2],
    final_mu_3 = final_theta_iterative$mu[3],
    final_w_hat = final_theta_iterative$weight_uniform,
    final_alpha_hat = final_theta_iterative$alpha,
    final_beta_hat = final_theta_iterative$beta,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_root_iterative,
    "convergence_summary.csv"
  ),
  row.names = FALSE
)

# ================================================================
# Final fast multiplier bootstrap only after the iterative diagnostic.
# ================================================================

run_final_fast_statistic_iterative <- function(
    statistic_name,
    statistic_seed) {
  ks_grid <- if (identical(statistic_name, "ks")) {
    make_sample_unique_distance_ks_grid()
  } else {
    NULL
  }

  result <- multiplier_bootstrap_gof(
    data = final_x_iterative,
    spec = spec_iterative,
    null = list(type = "composite"),
    statistics = statistic_name,
    ks_grid = ks_grid,
    B = bootstrap_B_iterative,
    alpha = 0.05,
    n_cores = bootstrap_n_cores_iterative,
    seed = statistic_seed,
    observed_theta_hat = final_theta_iterative,
    bootstrap_method = bootstrap_method_iterative,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = utils::modifyList(
      control_uniform_iterative,
      list(
        progress_bar = TRUE,
        progress_label = paste(
          "final",
          statistic_name,
          sep = " / "
        )
      )
    )
  )

  if (!identical(
    result$diagnostics$effective_bootstrap_method,
    "fast_multiplier"
  )) {
    stop(sprintf(
      "Final %s bootstrap did not use effective fast multiplier.",
      statistic_name
    ))
  }

  result
}

final_ks_result_iterative <- run_final_fast_statistic_iterative(
  statistic_name = "ks",
  statistic_seed = bootstrap_seed_ks_iterative
)

final_cvm_result_iterative <- run_final_fast_statistic_iterative(
  statistic_name = "cvm",
  statistic_seed = bootstrap_seed_cvm_iterative
)

saveRDS(
  final_ks_result_iterative,
  file = file.path(output_root_iterative, "final_gof_ks_fast.rds")
)

saveRDS(
  final_cvm_result_iterative,
  file = file.path(output_root_iterative, "final_gof_cvm_fast.rds")
)

utils::write.csv(
  data.frame(
    statistic = c("ks", "cvm"),
    observed = c(
      final_ks_result_iterative$inference$ks$observed,
      final_cvm_result_iterative$inference$cvm$observed
    ),
    p_value = c(
      final_ks_result_iterative$inference$ks$p_value,
      final_cvm_result_iterative$inference$cvm$p_value
    ),
    critical_value = c(
      final_ks_result_iterative$inference$ks$critical_value,
      final_cvm_result_iterative$inference$cvm$critical_value
    ),
    bootstrap_mean = c(
      mean(final_ks_result_iterative$bootstrap$statistics$ks),
      mean(final_cvm_result_iterative$bootstrap$statistics$cvm)
    ),
    bootstrap_sd = c(
      sd(final_ks_result_iterative$bootstrap$statistics$ks),
      sd(final_cvm_result_iterative$bootstrap$statistics$cvm)
    ),
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_root_iterative,
    "final_fast_bootstrap_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  data.frame(
    B = bootstrap_B_iterative,
    n_cores = bootstrap_n_cores_iterative,
    seed_ks = bootstrap_seed_ks_iterative,
    seed_cvm = bootstrap_seed_cvm_iterative,
    bootstrap_method = bootstrap_method_iterative,
    effective_bootstrap_method_ks = final_ks_result_iterative$diagnostics$effective_bootstrap_method,
    effective_bootstrap_method_cvm = final_cvm_result_iterative$diagnostics$effective_bootstrap_method,
    ks_grid_mode = "sample_points_unique_distances",
    distance_type = distance_type_iterative,
    final_n = nrow(final_x_iterative),
    final_alpha_hat = final_theta_iterative$alpha,
    final_beta_hat = final_theta_iterative$beta,
    final_w_hat = final_theta_iterative$weight_uniform,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_root_iterative,
    "final_fast_bootstrap_configuration.csv"
  ),
  row.names = FALSE
)

cat("Iterative diagnostic and final fast bootstrap written to:\n")
cat(normalizePath(output_root_iterative), "\n")
