source(file.path("distance_profiles", "rotational_mixtures_distance_profile_analysis.R"))

parse_named_args_rotmix_dp <- function(args) {
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
  args <- parse_named_args_rotmix_dp(commandArgs(trailingOnly = TRUE))

  model_name <- if (!is.null(args$model_name)) args$model_name else "rotational_beta_mixture2"
  distance_type <- if (!is.null(args$distance_type)) args$distance_type else "geodesic"
  sample_sizes <- if (!is.null(args$sample_sizes)) {
    as.integer(strsplit(args$sample_sizes, ",", fixed = TRUE)[[1L]])
  } else {
    c(50L, 200L)
  }

  result <- run_rotational_mixture_distance_profile_analysis(
    model_name = model_name,
    distance_type = distance_type,
    output_dir = args$output_dir %||% NULL,
    sample_sizes = sample_sizes,
    n_simulations = if (!is.null(args$n_simulations)) as.integer(args$n_simulations) else 10L,
    validation_n = if (!is.null(args$validation_n)) as.integer(args$validation_n) else 5000L,
    n_points = if (!is.null(args$n_points)) as.integer(args$n_points) else 200L,
    seed = if (!is.null(args$seed)) as.integer(args$seed) else 123L,
    save_plots = if (!is.null(args$save_plots)) isTRUE(as.logical(args$save_plots)) else TRUE
  )

  message(sprintf("Rotational mixture distance-profile output: %s", result$output_dir))
}
