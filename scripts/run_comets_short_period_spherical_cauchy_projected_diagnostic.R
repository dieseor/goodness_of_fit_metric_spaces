suppressPackageStartupMessages({
  source(file.path("utils.R"))
  source(file.path("distance_profiles", "spherical_cauchy_projected_poor_mans_diagnostic.R"))
  source(file.path("real_data", "comets", "utils_comets_data.R"))
  source(file.path("scripts", "path_helpers.R"))
})

run_sc_projected_poor_mans_diagnostic <- get("run_sc_projected_poor_mans_diagnostic", mode = "function")

load_comets_short_period_for_projection_sc <- function() {
  comets_data <- load_comets_real_data(finite_normals = "short")
  as.matrix(comets_data$short$normal)
}

parse_named_args_sc_projected <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(pieces) >= 2L) {
      key <- trimws(pieces[[1L]])
      val <- trimws(paste(pieces[-1L], collapse = "="))
      if (nzchar(key)) {
        out[[key]] <- val
      }
    }
  }
  out
}

parse_numeric_csv <- function(text, default) {
  if (is.null(text) || !nzchar(text)) {
    return(default)
  }
  parts <- strsplit(text, ",", fixed = TRUE)[[1L]]
  vals <- suppressWarnings(as.numeric(trimws(parts)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    default
  } else {
    vals
  }
}

run_comets_short_period_spherical_cauchy_projected_diagnostic <- function(
    output_dir = canonical_comets_spherical_cauchy_dir("projected_diagnostic"),
    sample_label = "short_period_comets_sc",
    rho_grid = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95),
    n_cores = 1L,
    add_rho_hat = TRUE,
    save_plot = TRUE,
    verbose = TRUE) {
  data_matrix <- load_comets_short_period_for_projection_sc()

  run_sc_projected_poor_mans_diagnostic(
    data = data_matrix,
    rho_grid = rho_grid,
    output_dir = output_dir,
    sample_label = sample_label,
    n_cores = as.integer(n_cores),
    add_rho_hat = isTRUE(add_rho_hat),
    save_plot = isTRUE(save_plot),
    verbose = isTRUE(verbose)
  )
}

if (identical(environment(), globalenv())) {
  args <- parse_named_args_sc_projected(commandArgs(trailingOnly = TRUE))

  output_dir <- if (!is.null(args$output_dir) && nzchar(args$output_dir)) {
    args$output_dir
  } else {
    canonical_comets_spherical_cauchy_dir("projected_diagnostic")
  }

  sample_label <- if (!is.null(args$sample_label) && nzchar(args$sample_label)) {
    args$sample_label
  } else {
    "short_period_comets_sc"
  }

  rho_grid <- parse_numeric_csv(
    args$rho_grid,
    default = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95)
  )

  n_cores <- if (!is.null(args$n_cores) && nzchar(args$n_cores)) {
    as.integer(args$n_cores)
  } else {
    1L
  }
  if (!is.finite(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }

  add_rho_hat <- if (!is.null(args$add_rho_hat) && nzchar(args$add_rho_hat)) {
    tolower(args$add_rho_hat) %in% c("1", "true", "t", "yes", "y")
  } else {
    TRUE
  }

  invisible(run_comets_short_period_spherical_cauchy_projected_diagnostic(
    output_dir = output_dir,
    sample_label = sample_label,
    rho_grid = rho_grid,
    n_cores = n_cores,
    add_rho_hat = add_rho_hat,
    save_plot = TRUE,
    verbose = TRUE
  ))
}
