#!/usr/bin/env Rscript

# Build the final real-data logistic-Gaussian screening without rerunning
# compatible B = 5000 results from the existing paper run.

resolve_consolidated_screening_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_consolidated_screening_path(
  "real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"
))
source(resolve_consolidated_screening_path("scripts", "path_helpers.R"))

canonical_result_aliases <- c(Aar_oxides = "AarMajorOxides")

canonical_screening_reference_dir <- function() {
  canonical_logistic_gaussian_screening_dir("fast", "paper_results_B5000_sampleks")
}

canonical_screening_output_dir <- function() {
  canonical_logistic_gaussian_screening_dir("fast", "paper_results_B5000_consolidated")
}

result_path_for_dataset <- function(base_dir, dataset_name) {
  file.path(base_dir, sprintf("%s_results.rds", slugify_dataset_name(dataset_name)))
}

reference_result_path <- function(dataset_name,
                                  reference_dir = canonical_screening_reference_dir()) {
  stored_name <- unname(canonical_result_aliases[dataset_name])
  if (is.na(stored_name)) {
    stored_name <- dataset_name
  }
  result_path_for_dataset(reference_dir, stored_name)
}

is_compatible_b5000_result <- function(result, dataset_name) {
  if (!identical(result$data_prep$status %||% "failed", "ok") ||
      !identical(as.integer(result$settings$B), 5000L) ||
      !identical(result$bootstrap$mode, "composite_multiplier") ||
      !identical(result$settings$omega_grid_type, "sample_points") ||
      !identical(result$settings$t_grid_type, "sample_distances")) {
    return(FALSE)
  }

  if (!is.null(composition_registry[[dataset_name]])) {
    entry <- composition_registry[[dataset_name]]
    if (!identical(result$data_prep$component_names, entry$parts) ||
        isTRUE(result$data_prep$has_zeros) || isTRUE(result$data_prep$has_missing) ||
        sum(duplicated(result$data_prep$X_comp)) > 0L ||
        max(abs(rowSums(result$data_prep$X_closed) - 1)) > 1e-12) {
      return(FALSE)
    }
  }

  TRUE
}

classify_consolidated_dataset <- function(dataset_name,
                                          reference_dir = canonical_screening_reference_dir(),
                                          output_dir = NULL) {
  candidate_paths <- c(
    if (!is.null(output_dir)) result_path_for_dataset(output_dir, dataset_name),
    reference_result_path(dataset_name, reference_dir)
  )
  for (path in unique(candidate_paths[file.exists(candidate_paths)])) {
    result <- tryCatch(readRDS(path), error = function(e) e)
    if (!inherits(result, "error") && is_compatible_b5000_result(result, dataset_name)) {
      is_current_output <- !is.null(output_dir) && identical(
        normalizePath(path, mustWork = FALSE),
        normalizePath(result_path_for_dataset(output_dir, dataset_name), mustWork = FALSE)
      )
      return(list(
        dataset = dataset_name,
        initial_status = if (is_current_output) "completed" else "existing_valid_result",
        result_path = path,
        reason = if (is_current_output) "Completed in the consolidated B = 5000 run." else "Compatible composite-multiplier B = 5000 result."
      ))
    }
  }

  prepared <- tryCatch(prepare_composition_dataset(dataset_name), error = function(e) e)
  if (inherits(prepared, "error")) {
    return(list(
      dataset = dataset_name,
      initial_status = "failed",
      result_path = NA_character_,
      reason = conditionMessage(prepared)
    ))
  }
  if (!identical(prepared$status, "ok")) {
    return(list(
      dataset = dataset_name,
      initial_status = "failed",
      result_path = NA_character_,
      reason = paste(prepared$notes %||% "Data source was not found.", collapse = " | ")
    ))
  }
  if (isTRUE(prepared$has_zeros) || isTRUE(prepared$has_missing) ||
      (prepared$n_duplicate_rows %||% 0L) > 0L) {
    return(list(
      dataset = dataset_name,
      initial_status = "excluded",
      result_path = NA_character_,
      reason = "The source data fail the no-zeros, no-missing-values, no-duplicate-rows requirement for the Logistic Gaussian fit."
    ))
  }

  list(
    dataset = dataset_name,
    initial_status = "needs_run",
    result_path = NA_character_,
    reason = "No compatible B = 5000 result is available."
  )
}

consolidated_manifest <- function(reference_dir = canonical_screening_reference_dir(),
                                  output_dir = NULL) {
  datasets <- default_logistic_gaussian_screening_datasets()
  manifest <- lapply(
    datasets,
    classify_consolidated_dataset,
    reference_dir = reference_dir,
    output_dir = output_dir
  )
  data.frame(
    dataset = vapply(manifest, `[[`, character(1), "dataset"),
    initial_status = vapply(manifest, `[[`, character(1), "initial_status"),
    existing_result_path = vapply(manifest, `[[`, character(1), "result_path"),
    reason = vapply(manifest, `[[`, character(1), "reason"),
    stringsAsFactors = FALSE
  )
}

make_consolidated_overview_plot <- function(summary_df, output_dir) {
  keep <- summary_df[
    summary_df$result_status %in% c("existing_valid_result", "completed") &
      is.finite(summary_df$ks_pvalue) & is.finite(summary_df$cvm_pvalue),
    , drop = FALSE
  ]
  if (nrow(keep) == 0L) {
    return(NULL)
  }

  plot_df <- rbind(
    data.frame(dataset = keep$dataset, statistic = "KS", p_value = keep$ks_pvalue, data_source = keep$data_source),
    data.frame(dataset = keep$dataset, statistic = "CvM", p_value = keep$cvm_pvalue, data_source = keep$data_source)
  )
  plot_df$dataset <- factor(plot_df$dataset, levels = rev(unique(keep$dataset)))
  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dataset, y = p_value, colour = data_source, shape = statistic)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Bootstrap p-value", colour = "Source", shape = "Statistic") +
    ggplot2::theme_minimal(base_size = 11)
  path <- file.path(output_dir, "plots", "consolidated_ks_cvm_pvalues.png")
  ggplot2::ggsave(path, plot, width = 10, height = 12, dpi = 220)
  path
}

validate_consolidated_summary <- function(summary_df) {
  canonical_names <- names(composition_registry)
  external_names <- external_logistic_gaussian_screening_datasets()
  excluded_compositions <- c("Blood23", "Glacial", "Skulls", "jura259")

  stopifnot(length(composition_registry) == 28L)
  stopifnot(!anyDuplicated(names(composition_registry)))
  stopifnot(setequal(summary_df$dataset[summary_df$data_source == "compositions"], canonical_names))
  stopifnot(all(external_names %in% summary_df$dataset))
  stopifnot(!anyDuplicated(summary_df$dataset))
  stopifnot(!any(excluded_compositions %in% summary_df$dataset))
  stopifnot(!any(startsWith(summary_df$dataset, "sa.")))
  stopifnot(!anyNA(summary_df$result_status))
  stopifnot(!any(summary_df$result_status == "manual_review"))
  completed <- summary_df$result_status %in% c("existing_valid_result", "completed")
  stopifnot(all(summary_df$B[completed] == 5000L))
  invisible(TRUE)
}

run_consolidated_logistic_gaussian_screening <- function(
    output_dir = canonical_screening_output_dir(),
    reference_dir = canonical_screening_reference_dir(),
    B = 5000L,
    max_centers = 100L,
    n_t = 60L,
    seed = 123L,
    n_cores = 8L,
    make_plots = TRUE,
    verbose = TRUE) {
  stopifnot(identical(as.integer(B), 5000L))
  dirs <- ensure_logistic_gaussian_screening_directories(output_dir)
  manifest <- consolidated_manifest(reference_dir, output_dir)
  utils::write.csv(manifest, file.path(dirs$metadata, "consolidated_manifest.csv"), row.names = FALSE)

  to_run <- manifest$dataset[manifest$initial_status == "needs_run"]
  fresh_batch <- NULL
  if (length(to_run) > 0L) {
    fresh_batch <- run_logistic_gaussian_screening_batch(
      dataset_names = to_run,
      B = B,
      max_centers = max_centers,
      n_t = n_t,
      bootstrap_mode = "composite_multiplier",
      seed = seed,
      n_cores = n_cores,
      control = list(logistic_gaussian_quadform_method = "hbe"),
      omega_grid_type = "sample_points",
      t_grid_type = "sample_distances",
      make_plots = FALSE,
      output_dir = output_dir,
      verbose = verbose
    )
  }

  summary_rows <- lapply(seq_len(nrow(manifest)), function(i) {
    dataset_name <- manifest$dataset[[i]]
    initial_status <- manifest$initial_status[[i]]
    if (initial_status %in% c("existing_valid_result", "completed")) {
      result <- readRDS(manifest$existing_result_path[[i]])
      row <- make_logistic_gaussian_screening_summary_row(result)
      metadata <- screening_dataset_metadata(dataset_name)
      row$dataset <- dataset_name
      row$data_source <- metadata$data_source
      row$source_object <- metadata$source_object
      row$selected_parts <- metadata$selected_parts
      if (!nzchar(row$selected_parts)) {
        row$selected_parts <- paste(result$data_prep$component_names %||% character(0), collapse = ",")
      }
      row$result_status <- initial_status
      return(row)
    }

    if (identical(initial_status, "needs_run")) {
      path <- result_path_for_dataset(output_dir, dataset_name)
      if (file.exists(path)) {
        result <- readRDS(path)
        row <- make_logistic_gaussian_screening_summary_row(result)
        row$result_status <- if (identical(row$status, "ok")) "completed" else "failed"
        return(row)
      }
      row <- fresh_batch$summary[fresh_batch$summary$dataset == dataset_name, , drop = FALSE]
      row$result_status <- "failed"
      return(row)
    }

    metadata <- screening_dataset_metadata(dataset_name)
    final_status <- if (identical(initial_status, "excluded")) "excluded" else "failed"
    data.frame(
      dataset = dataset_name,
      data_source = metadata$data_source,
      source_object = metadata$source_object,
      selected_parts = metadata$selected_parts,
      status = final_status,
      result_status = final_status,
      source_package = NA_character_,
      source_dataset_name = NA_character_,
      n = NA_integer_, D = NA_integer_, zeros = NA_character_, missing = NA_character_,
      duplicate_rows = NA_integer_, missing_rows_removed = NA_integer_, zero_replacement = NA_real_,
      n_zeros_replaced = NA_integer_, min_component = NA_real_, ridge_added = NA_real_,
      min_eigenvalue = NA_real_, max_eigenvalue = NA_real_, condition_number = NA_real_,
      boundary_epsilon = NA_real_, omega_grid_construction = NA_character_, n_centers = NA_integer_,
      t_grid_max = NA_real_, ks_statistic = NA_real_, ks_pvalue = NA_real_, cvm_statistic = NA_real_,
      cvm_pvalue = NA_real_, mardia_skew_pvalue = NA_real_, mardia_kurtosis_pvalue = NA_real_,
      shapiro_min_pvalue = NA_real_, diagnosis = final_status, use_in_paper = "no", why = manifest$reason[[i]],
      bootstrap_mode = "composite_multiplier", bootstrap_engine = NA_character_, B = NA_integer_,
      seed = NA_integer_, elapsed_seconds = NA_real_, notes = manifest$reason[[i]], stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)
  validate_consolidated_summary(summary_df)
  summary_path <- file.path(output_dir, "summary_logistic_gaussian_screening.csv")
  utils::write.csv(summary_df, summary_path, row.names = FALSE)
  overview_plot <- if (isTRUE(make_plots)) make_consolidated_overview_plot(summary_df, output_dir) else NULL

  list(
    manifest = manifest,
    summary = summary_df,
    summary_csv = summary_path,
    overview_plot = overview_plot,
    fresh_batch = fresh_batch
  )
}

if (sys.nframe() == 0L) {
  invisible(run_consolidated_logistic_gaussian_screening())
}
