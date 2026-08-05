#!/usr/bin/env Rscript

# Standalone diagnostics for the temporal two-beta law used in the joint
# cycle-23 model. This runner intentionally does not execute GOF statistics,
# joint profiles, or any parametric bootstrap.

resolve_sunspots_temporal_diagnostics_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_sunspots_temporal_diagnostics_path(
  "real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"
))

sunspots_joint_temporal_sample_definitions <- function(input_csv) {
  if (!file.exists(input_csv)) stop(sprintf("Input CSV not found: %s", input_csv))
  input <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  if (!all(c("cycle", "date") %in% names(input)) || any(input$cycle != 23L)) {
    stop("`input_csv` must contain only cycle-23 observations with a `date` column.")
  }
  calendar_day <- as.Date(as.POSIXct(input$date, tz = "UTC"), tz = "UTC")
  if (anyNA(calendar_day)) stop("The input `date` column contains invalid timestamps.")
  data.frame(
    sample = c("full", "central"),
    sample_id = c("full_cycle", "central_interval"),
    label = c("Full cycle 23", "Central interval"),
    start_date = c(as.character(min(calendar_day)), "1997-06-01"),
    end_date = c(as.character(max(calendar_day) + 1), "2006-01-01"),
    stringsAsFactors = FALSE
  )
}

normalize_sunspots_joint_temporal_samples <- function(samples) {
  samples <- unique(tolower(trimws(as.character(samples))))
  if (length(samples) == 0L || any(!samples %in% c("full", "central"))) {
    stop("`samples` must be one or both of 'full' and 'central'.")
  }
  samples
}

sunspots_joint_temporal_sensitivity_seeds <- function(dequantization_seed,
                                                       dequantization_seeds = NULL,
                                                       n_default = 20L) {
  dequantization_seed <- as.integer(dequantization_seed)
  if (length(dequantization_seed) != 1L || !is.finite(dequantization_seed)) {
    stop("`dequantization_seed` must be one finite integer.")
  }
  if (is.null(dequantization_seeds)) {
    n_default <- as.integer(n_default)
    if (!is.finite(n_default) || n_default < 1L) stop("`n_default` must be a positive integer.")
    dequantization_seeds <- dequantization_seed + seq.int(0L, n_default - 1L)
  }
  dequantization_seeds <- as.integer(dequantization_seeds)
  if (length(dequantization_seeds) == 0L || any(!is.finite(dequantization_seeds))) {
    stop("`dequantization_seeds` must contain finite integers.")
  }
  unique(c(dequantization_seed, dequantization_seeds))
}

sunspots_joint_temporal_fit_row <- function(fit, n, sample_definition, seed) {
  criteria <- sunspots_joint_time_information_criteria(fit$loglik, n = n, n_parameters = 5L)
  boundary <- fit$boundary_diagnostics
  out <- data.frame(
    sample = sample_definition$sample,
    sample_id = sample_definition$sample_id,
    sample_label = sample_definition$label,
    start_date = sample_definition$start_date,
    end_date_exclusive = sample_definition$end_date,
    dequantization_seed = as.integer(seed),
    n = as.integer(n),
    weight1 = fit$weight1,
    alpha1 = fit$alpha1,
    beta1 = fit$beta1,
    alpha2 = fit$alpha2,
    beta2 = fit$beta2,
    mean1 = fit$mean1,
    mean2 = fit$mean2,
    loglik = fit$loglik,
    aic = criteria$aic,
    bic = criteria$bic,
    convergence = as.integer(fit$opt$convergence %||% NA_integer_),
    convergence_message = as.character(fit$opt$message %||% ""),
    n_starts = as.integer(fit$n_starts),
    n_finite_fits = as.integer(fit$n_successful_starts),
    n_converged_fits = as.integer(fit$n_converged_starts),
    selected_converged = isTRUE(fit$selected_converged),
    canonical_component_order = fit$mean1 <= fit$mean2,
    raw_component_order_swapped = isTRUE(fit$component_swapped),
    component_mean_gap = fit$mean2 - fit$mean1,
    any_near_lower_bound = any(boundary$near_lower_bound),
    any_near_upper_bound = any(boundary$near_upper_bound),
    any_near_bound = any(boundary$near_any_bound),
    stringsAsFactors = FALSE
  )
  for (row_index in seq_len(nrow(boundary))) {
    parameter <- boundary$parameter[[row_index]]
    out[[paste0(parameter, "_lower_bound")]] <- boundary$lower_bound[[row_index]]
    out[[paste0(parameter, "_upper_bound")]] <- boundary$upper_bound[[row_index]]
    out[[paste0(parameter, "_distance_to_lower")]] <- boundary$distance_to_lower[[row_index]]
    out[[paste0(parameter, "_distance_to_upper")]] <- boundary$distance_to_upper[[row_index]]
    out[[paste0(parameter, "_near_lower_bound")]] <- boundary$near_lower_bound[[row_index]]
    out[[paste0(parameter, "_near_upper_bound")]] <- boundary$near_upper_bound[[row_index]]
    out[[paste0(parameter, "_near_any_bound")]] <- boundary$near_any_bound[[row_index]]
  }
  out
}

sunspots_joint_temporal_sensitivity_summary <- function(fits) {
  parameters <- c("weight1", "alpha1", "beta1", "alpha2", "beta2")
  parameter_summary <- do.call(rbind, lapply(parameters, function(parameter) {
    values <- fits[[parameter]]
    data.frame(
      parameter = parameter,
      n_seeds = nrow(fits),
      mean = mean(values), sd = stats::sd(values), min = min(values), max = max(values),
      lower_bound_frequency = mean(fits[[paste0(parameter, "_near_lower_bound")]]),
      upper_bound_frequency = mean(fits[[paste0(parameter, "_near_upper_bound")]]),
      any_bound_frequency = mean(fits[[paste0(parameter, "_near_any_bound")]]),
      stringsAsFactors = FALSE
    )
  }))
  overview <- data.frame(
    n_seeds = nrow(fits),
    selected_converged_frequency = mean(fits$selected_converged),
    canonical_component_order_frequency = mean(fits$canonical_component_order),
    raw_component_order_swapped_frequency = mean(fits$raw_component_order_swapped),
    mean_gap_min = min(fits$component_mean_gap),
    mean_gap_max = max(fits$component_mean_gap),
    stringsAsFactors = FALSE
  )
  list(parameter_summary = parameter_summary, overview = overview)
}

sunspots_joint_temporal_rank_windows <- function(data) {
  required <- c("s", "calendar_day", "calendar_day_start", "dequantized_timestamp")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf("`data` is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  n <- nrow(data)
  if (n < 1L) stop("`data` must contain at least one observation.")
  lower_levels <- c(0, 0.10, 0.20, 0.40, 0.60)
  upper_levels <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  center_levels <- c(0.05, 0.15, 0.30, 0.50, 0.70)
  tie_breaker <- if ("NOAA" %in% names(data)) data$NOAA else seq_len(n)
  order_index <- order(data$s, tie_breaker, seq_len(n), na.last = TRUE)
  rows <- lapply(seq_along(lower_levels), function(index) {
    first_rank <- floor(lower_levels[[index]] * n) + 1L
    last_rank <- floor(upper_levels[[index]] * n)
    center_rank <- min(max(as.integer(ceiling(center_levels[[index]] * n)), 1L), n)
    center_index <- order_index[[center_rank]]
    member_indices <- if (last_rank >= first_rank) order_index[first_rank:last_rank] else integer(0L)
    data.frame(
      lower_rank_level = lower_levels[[index]],
      upper_rank_level = upper_levels[[index]],
      center_rank_level = center_levels[[index]],
      center_rank = center_rank,
      center_empirical_quantile = center_rank / n,
      n = length(member_indices),
      min_s = if (length(member_indices)) min(data$s[member_indices]) else NA_real_,
      max_s = if (length(member_indices)) max(data$s[member_indices]) else NA_real_,
      min_date = if (length(member_indices)) as.character(min(data$calendar_day[member_indices])) else NA_character_,
      max_date = if (length(member_indices)) as.character(max(data$calendar_day[member_indices])) else NA_character_,
      center_s = data$s[[center_index]],
      center_date = as.character(data$calendar_day[[center_index]]),
      center_dequantized_timestamp = as.character(data$dequantized_timestamp[[center_index]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

sunspots_joint_temporal_plot_default_fit <- function(data, fit, output_dir, sample_label) {
  s_grid <- seq(1e-5, 1 - 1e-5, length.out = 1001L)
  total_density <- sunspots_joint_time_density(s_grid, fit)
  component1 <- fit$weight1 * stats::dbeta(s_grid, fit$alpha1, fit$beta1)
  component2 <- (1 - fit$weight1) * stats::dbeta(s_grid, fit$alpha2, fit$beta2)
  pit <- sunspots_joint_time_cdf(data$s, fit)
  parameter_label <- sprintf(
    "w=%.6f; (a1,b1)=(%.6f,%.6f); (a2,b2)=(%.6f,%.6f)",
    fit$weight1, fit$alpha1, fit$beta1, fit$alpha2, fit$beta2
  )

  density_path <- file.path(output_dir, "temporal_density_components.png")
  grDevices::png(density_path, width = 1500, height = 950, res = 150)
  hist(data$s, breaks = 40L, freq = FALSE, col = "#d9e8f5", border = "white",
       xlim = c(0, 1), ylim = c(0, max(total_density, component1, component2) * 1.08),
       xlab = "Dequantized first-record day (S)", ylab = "Density",
       main = sprintf("%s: temporal two-beta fit", sample_label))
  lines(s_grid, total_density, col = "#111111", lwd = 3)
  lines(s_grid, component1, col = "#0072B2", lwd = 2, lty = 2)
  lines(s_grid, component2, col = "#D55E00", lwd = 2, lty = 3)
  legend("topright", legend = c("Empirical histogram", "Fitted total", "Weighted beta 1", "Weighted beta 2", parameter_label),
         col = c("#d9e8f5", "#111111", "#0072B2", "#D55E00", NA),
         lwd = c(8, 3, 2, 2, NA), lty = c(1, 1, 2, 3, NA), bty = "n", cex = 0.85)
  grDevices::dev.off()

  cdf_path <- file.path(output_dir, "temporal_empirical_fitted_cdf.png")
  grDevices::png(cdf_path, width = 1500, height = 950, res = 150)
  ordered_s <- sort(data$s)
  plot(ordered_s, seq_along(ordered_s) / length(ordered_s), type = "s", xlim = c(0, 1), ylim = c(0, 1),
       col = "#111111", lwd = 2, xlab = "Dequantized first-record day (S)", ylab = "CDF",
       main = sprintf("%s: empirical and fitted temporal CDF", sample_label))
  lines(s_grid, sunspots_joint_time_cdf(s_grid, fit), col = "#8B0000", lwd = 3)
  abline(0, 1, col = "#999999", lty = 3)
  legend("topleft", legend = c("Empirical CDF", "Fitted CDF"), col = c("#111111", "#8B0000"), lwd = c(2, 3), bty = "n")
  grDevices::dev.off()

  qq_path <- file.path(output_dir, "temporal_pit_qq.png")
  grDevices::png(qq_path, width = 1500, height = 950, res = 150)
  qqplot(stats::ppoints(length(pit)), sort(pit), xlim = c(0, 1), ylim = c(0, 1),
         pch = 16, cex = 0.55, col = "#0072B2", xlab = "Uniform quantiles",
         ylab = "Temporal PIT quantiles", main = sprintf("%s: temporal PIT Q-Q plot", sample_label))
  abline(0, 1, col = "#8B0000", lwd = 2)
  grDevices::dev.off()

  pit_path <- file.path(output_dir, "temporal_pit_histogram.png")
  grDevices::png(pit_path, width = 1500, height = 950, res = 150)
  hist(pit, breaks = 20L, freq = FALSE, xlim = c(0, 1), col = "#d9e8f5", border = "white",
       xlab = "Temporal PIT", ylab = "Density", main = sprintf("%s: temporal PIT histogram", sample_label))
  abline(h = 1, col = "#8B0000", lwd = 2)
  grDevices::dev.off()

  list(density = density_path, cdf = cdf_path, qq = qq_path, pit_histogram = pit_path)
}

run_sunspots_cycle23_temporal_beta_diagnostics <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
    output_root = file.path("real_data", "sunspots", "output", "cycle23_temporal_beta_diagnostics"),
    samples = c("full", "central"),
    dequantization_seed = 20260712L,
    dequantization_seeds = NULL,
    control = list()) {
  samples <- normalize_sunspots_joint_temporal_samples(samples)
  seed_values <- sunspots_joint_temporal_sensitivity_seeds(
    dequantization_seed = dequantization_seed, dequantization_seeds = dequantization_seeds
  )
  definitions <- sunspots_joint_temporal_sample_definitions(input_csv)
  definitions <- definitions[definitions$sample %in% samples, , drop = FALSE]
  results <- vector("list", nrow(definitions))

  for (definition_index in seq_len(nrow(definitions))) {
    definition <- definitions[definition_index, , drop = FALSE]
    sample_dir <- file.path(output_root, definition$sample_id)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    fit_rows <- vector("list", length(seed_values))
    default_data <- NULL
    default_fit <- NULL
    for (seed_index in seq_along(seed_values)) {
      seed <- seed_values[[seed_index]]
      data <- prepare_sunspots_joint_time_data(
        input_csv = input_csv, start_date = definition$start_date,
        end_date = definition$end_date, dequantization_seed = seed
      )
      fit <- fit_sunspots_joint_time_beta_mixture2(data$s, control = control)
      fit_rows[[seed_index]] <- sunspots_joint_temporal_fit_row(fit, nrow(data), definition, seed)
      if (identical(seed, as.integer(dequantization_seed))) {
        default_data <- data
        default_fit <- fit
      }
    }
    fits <- do.call(rbind, fit_rows)
    summary <- sunspots_joint_temporal_sensitivity_summary(fits)
    default_row <- fits[fits$dequantization_seed == as.integer(dequantization_seed), , drop = FALSE]
    windows <- sunspots_joint_temporal_rank_windows(default_data)
    default_data$temporal_fitted_density <- sunspots_joint_time_density(default_data$s, default_fit)
    default_data$temporal_pit <- sunspots_joint_time_cdf(default_data$s, default_fit)
    plot_paths <- sunspots_joint_temporal_plot_default_fit(default_data, default_fit, sample_dir, definition$label)

    utils::write.csv(default_row, file.path(sample_dir, "temporal_fit_default_seed.csv"), row.names = FALSE)
    utils::write.csv(default_fit$boundary_diagnostics, file.path(sample_dir, "temporal_boundary_diagnostics_default_seed.csv"), row.names = FALSE)
    utils::write.csv(fits, file.path(sample_dir, "temporal_beta_sensitivity_fits.csv"), row.names = FALSE)
    utils::write.csv(summary$parameter_summary, file.path(sample_dir, "temporal_beta_sensitivity_parameter_summary.csv"), row.names = FALSE)
    utils::write.csv(summary$overview, file.path(sample_dir, "temporal_beta_sensitivity_overview.csv"), row.names = FALSE)
    utils::write.csv(windows, file.path(sample_dir, "future_time_rank_windows.csv"), row.names = FALSE)
    utils::write.csv(default_data, file.path(sample_dir, "temporal_retained_data_default_seed.csv"), row.names = FALSE)
    utils::write.csv(data.frame(
      sample = definition$sample, sample_id = definition$sample_id,
      start_date = definition$start_date, end_date_exclusive = definition$end_date,
      default_dequantization_seed = as.integer(dequantization_seed),
      sensitivity_seeds = paste(seed_values, collapse = ","), stringsAsFactors = FALSE
    ), file.path(sample_dir, "temporal_diagnostic_metadata.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(sample_dir, "sessionInfo.txt"))
    results[[definition$sample]] <- list(
      fits = fits, summary = summary, windows = windows, plots = plot_paths, output_dir = sample_dir
    )
  }
  invisible(results)
}

parse_sunspots_joint_temporal_diagnostic_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --samples=full,central --input_csv=PATH --output_root=PATH ",
        "--dequantization_seed=INTEGER --dequantization_seeds=INTEGER,...\n"
      ))
      quit(save = "no", status = 0L)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop(sprintf("Invalid option: %s", arg))
    key <- parts[[1L]]
    value <- parts[[2L]]
    if (identical(key, "samples")) out$samples <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (key %in% c("input_csv", "output_root")) out[[key]] <- value
    if (identical(key, "dequantization_seed")) out[[key]] <- as.integer(value)
    if (identical(key, "dequantization_seeds")) {
      out[[key]] <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
    }
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(run_sunspots_cycle23_temporal_beta_diagnostics, parse_sunspots_joint_temporal_diagnostic_args())
}