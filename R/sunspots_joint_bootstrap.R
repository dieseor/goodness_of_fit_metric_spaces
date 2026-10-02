multiplier_bootstrap_sunspots_joint_time_space <- function(
    data,
    null,
    statistics = c("ks", "cvm"),
    ks_grid = NULL,
    B = 5000,
    alpha = 0.05,
    multipliers = NULL,
    n_cores = 1,
    seed = NULL,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(),
    distance_profile_backend = "cpp") {
  if (!identical(bootstrap_method, "fast_multiplier")) {
    stop("The joint sunspots model supports only `bootstrap_method = 'fast_multiplier'`.",
         call. = FALSE)
  }
  if (!is.list(null) || !identical(null$type, "composite")) {
    stop("The fast joint sunspots test currently requires `null = list(type = 'composite')`.",
         call. = FALSE)
  }
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  if (!is.null(control$hemisphere_regression) &&
      !identical(control$hemisphere_regression, "shared")) {
    stop("The published joint sunspots model requires `hemisphere_regression = 'shared'`.",
         call. = FALSE)
  }
  control$hemisphere_regression <- "shared"
  if (is.null(ks_grid)) ks_grid <- make_sample_unique_distance_ks_grid()
  multiplier_bootstrap_gof(
    data = data,
    spec = make_sunspots_joint_time_space_spec("shared"),
    null = null,
    statistics = statistics,
    ks_grid = ks_grid,
    B = B,
    alpha = alpha,
    multipliers = multipliers,
    n_cores = n_cores,
    seed = seed,
    bootstrap_method = bootstrap_method,
    keep = keep,
    control = control,
    distance_profile_backend = distance_profile_backend
  )
}
