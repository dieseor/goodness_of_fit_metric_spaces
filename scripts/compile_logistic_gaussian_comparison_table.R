#!/usr/bin/env Rscript

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

resolve_comparison_table_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_comparison_table_path(
  "real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"
))

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

coerce_shapiro_min <- function(p_values) {
  value <- suppressWarnings(min(as.numeric(p_values), na.rm = TRUE))
  if (!is.finite(value)) NA_real_ else value
}

extract_row <- function(result_path) {
  result <- readRDS(result_path)
  data.frame(
    dataset = result$dataset_name %||% NA_character_,
    n = result$data_prep$n %||% NA_integer_,
    D = result$data_prep$D %||% NA_integer_,
    ks_p = result$inference$ks$p_value %||% NA_real_,
    cvm_p = result$inference$cvm$p_value %||% NA_real_,
    mardia_skew_p = result$diagnostics$mardia$skewness_p_value %||% NA_real_,
    mardia_kurt_p = result$diagnostics$mardia$kurtosis_p_value %||% NA_real_,
    shapiro_min_p = coerce_shapiro_min(result$diagnostics$shapiro_p_values),
    diagnosis = result$classification$diagnosis %||% NA_character_,
    result_path = result_path,
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

dataset_names <- if (is.null(args$datasets)) {
  default_logistic_gaussian_screening_datasets()
} else {
  trimws(strsplit(args$datasets, ",", fixed = TRUE)[[1]])
}

exclude_datasets <- if (is.null(args$exclude)) character(0) else trimws(strsplit(args$exclude, ",", fixed = TRUE)[[1]])
dataset_names <- setdiff(dataset_names, exclude_datasets)

search_dirs <- if (is.null(args$search_dirs)) {
  c(
    file.path("output", "calibration", "bootstrap", "logistic_gaussian", "definitive15_B1000_12cores")
  )
} else {
  trimws(strsplit(args$search_dirs, ",", fixed = TRUE)[[1]])
}

output_csv <- args$output_csv %||% file.path("output", "logistic_gaussian_comparison_table.csv")

dataset_to_slug <- function(name) {
  gsub("[^a-z0-9]+", "_", tolower(name))
}

rows <- lapply(dataset_names, function(dataset_name) {
  slug <- dataset_to_slug(dataset_name)
  candidates <- file.path(search_dirs, sprintf("%s_results.rds", slug))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    return(data.frame(
      dataset = dataset_name,
      n = NA_integer_,
      D = NA_integer_,
      ks_p = NA_real_,
      cvm_p = NA_real_,
      mardia_skew_p = NA_real_,
      mardia_kurt_p = NA_real_,
      shapiro_min_p = NA_real_,
      diagnosis = "pending",
      result_path = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  extract_row(existing[[1L]])
})

comparison_df <- do.call(rbind, rows)
comparison_df <- comparison_df[order(comparison_df$dataset), , drop = FALSE]

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(comparison_df, output_csv, row.names = FALSE)

cat("Wrote comparison CSV to ", output_csv, "\n", sep = "")
print(comparison_df, row.names = FALSE)
