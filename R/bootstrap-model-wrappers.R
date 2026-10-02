multiplier_bootstrap_normal <- function(data,
                                        null,
                                        statistics = c("ks", "cvm"),
                                        ks_grid = NULL,
                                        B = 5000,
                                        alpha = 0.05,
                                        multipliers = NULL,
                                        n_cores = 1,
                                        seed = NULL,
                                        bootstrap_method = c("reestimated", "fast_multiplier"),
                                        keep = list(
                                          observed_process = TRUE,
                                          bootstrap_statistics = TRUE,
                                          bootstrap_thetas = FALSE
                                        ),
                                        control = list(),
                                        unknown_param = NULL,
                                        distance_profile_backend = c("r", "cpp")) {
  spec <- make_normal_spec(unknown_param = unknown_param)
  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
multiplier_bootstrap_mvnormal <- function(data,
                                          null,
                                          statistics = c("ks", "cvm"),
                                          ks_grid = NULL,
                                          B = 5000,
                                          alpha = 0.05,
                                          multipliers = NULL,
                                          n_cores = 1,
                                          seed = NULL,
                                          bootstrap_method = c("reestimated", "fast_multiplier"),
                                          keep = list(
                                            observed_process = TRUE,
                                            bootstrap_statistics = TRUE,
                                            bootstrap_thetas = FALSE
                                          ),
                                          control = list(),
                                          unknown_param = "both",
                                          distance_profile_backend = c("r", "cpp"),
                                          fast_multiplier_backend = c("cpp", "r"),
                                          fuse_ks_cvm = TRUE,
                                          cache_block_corrections = c(
                                            "auto", "true", "false"
                                          )) {
  fast_multiplier_backend <- normalize_fast_multiplier_backend(
    fast_multiplier_backend
  )
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(fuse_ks_cvm)
  cache_block_corrections <- normalize_fast_multiplier_cache(
    cache_block_corrections
  )
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  requested_method <- tolower(as.character(control$derivative_method %||% "auto"))
  effective_method <- if (requested_method %in% c("auto", "deterministic")) {
    "score_mc"
  } else {
    requested_method
  }
  selection_source <- if (!is.null(control$derivative_method)) {
    if (identical(requested_method, "auto")) "explicit_auto" else "explicit"
  } else {
    "model_default"
  }
  control$derivative_method <- effective_method
  control$derivative_mc_size <- as.integer(control$derivative_mc_size %||% 10000L)
  spec <- make_mvnormal_spec(unknown_param = unknown_param)
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
  result$diagnostics$fast_multiplier_backend_requested <-
    fast_multiplier_backend
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <-
    fuse_ks_cvm
  result$diagnostics$fast_multiplier_cache_corrections_requested <-
    cache_block_corrections
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$derivative_method_selection_source <- selection_source
  if (!identical(
      result$diagnostics$effective_bootstrap_method,
      "fast_multiplier"
  )) {
    result$diagnostics$fast_multiplier_backend_effective <- "r"
    result$diagnostics$fast_multiplier_fuse_ks_cvm_effective <- FALSE
    result$diagnostics$fast_multiplier_cache_corrections_effective <- FALSE
  }
  result
}
multiplier_bootstrap_vmf <- function(data,
                                     null,
                                     statistics = c("ks", "cvm"),
                                     ks_grid = NULL,
                                     B = 5000,
                                     alpha = 0.05,
                                     multipliers = NULL,
                                     n_cores = 1,
                                     seed = NULL,
                                     bootstrap_method = c("reestimated", "fast_multiplier"),
                                     keep = list(
                                       observed_process = TRUE,
                                       bootstrap_statistics = TRUE,
                                       bootstrap_thetas = FALSE
                                     ),
                                     control = list(),
                                     distance_type = c("chordal", "geodesic"),
                                     unknown_param = "xi",
                                     distance_profile_backend = c("r", "cpp"),
                                     fast_multiplier_backend = c("cpp", "r"),
                                     fast_multiplier_cpp_kernel = c(
                                       "contiguous_double", "legacy"
                                     ),
                                     fuse_ks_cvm = TRUE,
                                     cache_block_corrections = c(
                                       "auto", "true", "false"
                                     )) {
  distance_type <- match.arg(distance_type)
  fast_multiplier_backend <- normalize_fast_multiplier_backend(
    if (missing(fast_multiplier_backend)) {
      control$fast_multiplier_backend %||% "cpp"
    } else {
      fast_multiplier_backend
    }
  )
  fast_multiplier_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(
    if (missing(fast_multiplier_cpp_kernel)) {
      control$fast_multiplier_cpp_kernel %||% "contiguous_double"
    } else {
      fast_multiplier_cpp_kernel
    }
  )
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(
    if (missing(fuse_ks_cvm)) {
      control$fast_multiplier_fuse_ks_cvm %||% TRUE
    } else {
      fuse_ks_cvm
    }
  )
  cache_block_corrections <- normalize_fast_multiplier_cache(
    if (missing(cache_block_corrections)) {
      control$fast_multiplier_cache_corrections %||% "auto"
    } else {
      cache_block_corrections
    }
  )
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_cpp_kernel <- fast_multiplier_cpp_kernel
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  spec <- make_vmf_spec(
    distance_type = distance_type,
    unknown_param = unknown_param
  )

  legacy_mc_control <- is.null(control$derivative_method) &&
    (!is.null(control$derivative_mc_size) ||
       !is.null(control$derivative_mc_seed))
  requested_method <- control$derivative_method %||% "quadrature"
  selection_source <- if (!is.null(control$derivative_method)) {
    "explicit"
  } else if (legacy_mc_control) {
    "model_default_legacy_mc_controls_ignored"
  } else {
    "model_default"
  }
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_selection_source <- selection_source
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$fast_multiplier_backend_requested <-
    fast_multiplier_backend
  result$diagnostics$fast_multiplier_cpp_kernel_requested <-
    fast_multiplier_cpp_kernel
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <-
    fuse_ks_cvm
  result$diagnostics$fast_multiplier_cache_corrections_requested <-
    cache_block_corrections
  if (!identical(
      result$diagnostics$effective_bootstrap_method,
      "fast_multiplier"
  )) {
    result$diagnostics$fast_multiplier_backend_effective <- "r"
    result$diagnostics$fast_multiplier_cpp_kernel_effective <- "not_used"
    result$diagnostics$fast_multiplier_fuse_ks_cvm_effective <- FALSE
    result$diagnostics$fast_multiplier_cache_corrections_effective <- FALSE
  }
  result
}
multiplier_bootstrap_hvmf <- function(data,
                                      null,
                                      statistics = c("ks", "cvm"),
                                      ks_grid = NULL,
                                      B = 5000,
                                      alpha = 0.05,
                                      multipliers = NULL,
                                      n_cores = 1,
                                      seed = NULL,
                                      bootstrap_method = c("reestimated", "fast_multiplier"),
                                      keep = list(
                                        observed_process = TRUE,
                                        bootstrap_statistics = TRUE,
                                        bootstrap_thetas = FALSE
                                      ),
                                      control = list(),
                                      unknown_param = "both",
                                      distance_profile_backend = c("r", "cpp"),
                                      fast_multiplier_backend = c("cpp", "r"),
                                      fast_multiplier_cpp_kernel = c(
                                        "contiguous_double", "legacy"
                                      ),
                                      fuse_ks_cvm = TRUE,
                                      cache_block_corrections = c(
                                        "auto", "true", "false"
                                      )) {
  fast_multiplier_backend <- normalize_fast_multiplier_backend(
    if (missing(fast_multiplier_backend)) {
      control$fast_multiplier_backend %||% "cpp"
    } else {
      fast_multiplier_backend
    }
  )
  fast_multiplier_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(
    if (missing(fast_multiplier_cpp_kernel)) {
      control$fast_multiplier_cpp_kernel %||% "contiguous_double"
    } else {
      fast_multiplier_cpp_kernel
    }
  )
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(
    if (missing(fuse_ks_cvm)) {
      control$fast_multiplier_fuse_ks_cvm %||% TRUE
    } else {
      fuse_ks_cvm
    }
  )
  cache_block_corrections <- normalize_fast_multiplier_cache(
    if (missing(cache_block_corrections)) {
      control$fast_multiplier_cache_corrections %||% "auto"
    } else {
      cache_block_corrections
    }
  )
  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_cpp_kernel <- fast_multiplier_cpp_kernel
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections
  spec <- make_hvmf_spec(unknown_param = unknown_param)

  legacy_mc_control <- is.null(control$derivative_method) &&
    (!is.null(control$derivative_mc_size) ||
       !is.null(control$derivative_mc_seed))
  requested_method <- control$derivative_method %||% "quadrature"
  selection_source <- if (!is.null(control$derivative_method)) {
    "explicit"
  } else if (legacy_mc_control) {
    "model_default_legacy_mc_controls_ignored"
  } else {
    "model_default"
  }
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
  result$diagnostics$fast_multiplier_backend_requested <-
    fast_multiplier_backend
  result$diagnostics$fast_multiplier_cpp_kernel_requested <-
    fast_multiplier_cpp_kernel
  result$diagnostics$fast_multiplier_fuse_ks_cvm_requested <-
    fuse_ks_cvm
  result$diagnostics$fast_multiplier_cache_corrections_requested <-
    cache_block_corrections
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_selection_source <- selection_source
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  if (!identical(
      result$diagnostics$effective_bootstrap_method,
      "fast_multiplier"
  )) {
    result$diagnostics$fast_multiplier_backend_effective <- "r"
    result$diagnostics$fast_multiplier_fuse_ks_cvm_effective <- FALSE
    result$diagnostics$fast_multiplier_cache_corrections_effective <- FALSE
  }
  result
}
multiplier_bootstrap_logistic_gaussian <- function(data,
                                                   null,
                                                   statistics = c("ks", "cvm"),
                                                   ks_grid = NULL,
                                                   B = 5000,
                                                   alpha = 0.05,
                                                   multipliers = NULL,
                                                   n_cores = 1,
                                                   seed = NULL,
                                                   bootstrap_method = c("reestimated", "fast_multiplier"),
                                                   keep = list(
                                                     observed_process = TRUE,
                                                     bootstrap_statistics = TRUE,
                                                     bootstrap_thetas = FALSE
                                                   ),
                                                   control = list(),
                                                   unknown_param = "both",
                                                   distance_profile_backend = c("r", "cpp")) {
  requested_method <- tolower(as.character(control$derivative_method %||% "auto"))
  effective_method <- if (requested_method %in% c("auto", "deterministic")) {
    "score_mc"
  } else {
    requested_method
  }
  selection_source <- if (!is.null(control$derivative_method)) {
    if (identical(requested_method, "auto")) "explicit_auto" else "explicit"
  } else {
    "model_default"
  }
  control$derivative_method <- effective_method
  control$derivative_mc_size <- as.integer(control$derivative_mc_size %||% 10000L)
  spec <- make_logistic_gaussian_spec(unknown_param = unknown_param)

  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
  result$diagnostics$derivative_method_requested <- requested_method
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$derivative_method_selection_source <- selection_source
  result
}
multiplier_bootstrap_uniform_beta_mixture <- function(data,
                                                      null,
                                                      statistics = c("ks", "cvm"),
                                                      ks_grid = NULL,
                                                      B = 5000,
                                                      alpha = 0.05,
                                                      multipliers = NULL,
                                                      n_cores = 1,
                                                      seed = NULL,
                                                      bootstrap_method = c("reestimated", "fast_multiplier"),
                                                      keep = list(
                                                        observed_process = TRUE,
                                                        bootstrap_statistics = TRUE,
                                                        bootstrap_thetas = FALSE
                                                      ),
                                                      control = list(),
                                                      distance_type = c("chordal", "geodesic"),
                                                      distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_uniform_beta_mixture_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
multiplier_bootstrap_cardioid <- function(data,
                                          null,
                                          k,
                                          statistics = c("ks", "cvm"),
                                          ks_grid = NULL,
                                          B = 5000,
                                          alpha = 0.05,
                                          multipliers = NULL,
                                          n_cores = 1,
                                          seed = NULL,
                                          bootstrap_method = c("reestimated", "fast_multiplier"),
                                          keep = list(
                                            observed_process = TRUE,
                                            bootstrap_statistics = TRUE,
                                            bootstrap_thetas = FALSE
                                          ),
                                          control = list(),
                                          distance_type = c("chordal", "geodesic"),
                                          unknown_param = "both",
                                          distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_cardioid_spec(
    k = as.integer(k),
    distance_type = distance_type,
    unknown_param = unknown_param
  )

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
multiplier_bootstrap_small_circle <- function(data,
                                              null,
                                              statistics = c("ks", "cvm"),
                                              ks_grid = NULL,
                                              B = 5000,
                                              alpha = 0.05,
                                              multipliers = NULL,
                                              n_cores = 1,
                                              seed = NULL,
                                              bootstrap_method = c("reestimated", "fast_multiplier"),
                                              keep = list(
                                                observed_process = TRUE,
                                                bootstrap_statistics = TRUE,
                                                bootstrap_thetas = FALSE
                                                ),
                                                control = list(),
                                                distance_type = c("chordal", "geodesic"),
                                                distance_profile_backend = c("r", "cpp")) {
  distance_type <- match.arg(distance_type)
  spec <- make_small_circle_spec(distance_type = distance_type)

  multiplier_bootstrap_gof(
    data = data,
    spec = spec,
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
