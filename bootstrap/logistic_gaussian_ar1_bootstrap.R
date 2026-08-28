# Public bootstrap wrapper for the restricted Logistic-Gaussian AR(1) family
#   ilr(X) ~ N_d(mu, R_d(rho)).

if (!exists("multiplier_bootstrap_gof", mode = "function")) {
  candidates <- c(
    file.path("bootstrap", "multiplier_bootstrap.R"),
    "multiplier_bootstrap.R",
    file.path("..", "bootstrap", "multiplier_bootstrap.R"),
    file.path("..", "..", "bootstrap", "multiplier_bootstrap.R")
  )
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("Could not locate bootstrap/multiplier_bootstrap.R.")
  source(path)
}

if (!exists("make_logistic_gaussian_ar1_spec", mode = "function")) {
  candidates <- c(
    file.path("bootstrap", "logistic_gaussian_ar1_model_spec.R"),
    "logistic_gaussian_ar1_model_spec.R",
    file.path("..", "bootstrap", "logistic_gaussian_ar1_model_spec.R"),
    file.path("..", "..", "bootstrap", "logistic_gaussian_ar1_model_spec.R")
  )
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) {
    stop("Could not locate bootstrap/logistic_gaussian_ar1_model_spec.R.")
  }
  source(path)
}

multiplier_bootstrap_logistic_gaussian_ar1 <- function(
    data,
    null,
    statistics = c("ks", "cvm"),
    ks_grid = NULL,
    B = 5000L,
    alpha = 0.05,
    multipliers = NULL,
    n_cores = 1L,
    seed = NULL,
    bootstrap_method = c("reestimated", "fast_multiplier"),
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(),
    distance_profile_backend = c("r", "cpp"),
    fast_multiplier_backend = c("cpp", "r"),
    fast_multiplier_cpp_kernel = c("contiguous_double", "legacy"),
    fuse_ks_cvm = TRUE,
    cache_block_corrections = c("auto", "true", "false")) {

  fast_multiplier_backend <-
    normalize_fast_multiplier_backend(fast_multiplier_backend)
  fast_multiplier_cpp_kernel <-
    normalize_fast_multiplier_cpp_kernel(fast_multiplier_cpp_kernel)
  fuse_ks_cvm <- normalize_fast_multiplier_fusion(fuse_ks_cvm)
  cache_block_corrections <-
    normalize_fast_multiplier_cache(cache_block_corrections)

  control$fast_multiplier_backend <- fast_multiplier_backend
  control$fast_multiplier_cpp_kernel <- fast_multiplier_cpp_kernel
  control$fast_multiplier_fuse_ks_cvm <- fuse_ks_cvm
  control$fast_multiplier_cache_corrections <- cache_block_corrections

  derivative_method_was_supplied <- !is.null(control$derivative_method)
  requested_derivative_method <-
    tolower(as.character(control$derivative_method %||% "auto"))

  if (length(requested_derivative_method) != 1L ||
      !requested_derivative_method %in% c("auto", "score_mc")) {
    stop(paste(
      "The Logistic-Gaussian AR(1) fast multiplier supports only",
      "`derivative_method = 'score_mc'` (or 'auto')."
    ))
  }

  control$derivative_method <- "score_mc"
  control$derivative_mc_size <-
    as.integer(control$derivative_mc_size %||% 10000L)

  result <- multiplier_bootstrap_gof(
    data = data,
    spec = make_logistic_gaussian_ar1_spec(),
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
  result$diagnostics$derivative_method_requested <-
    requested_derivative_method
  result$diagnostics$derivative_method_effective <-
    result$diagnostics$derivative_method %||% NA_character_
  result$diagnostics$derivative_method_selection_source <-
    if (!derivative_method_was_supplied) {
      "model_default"
    } else if (identical(requested_derivative_method, "auto")) {
      "explicit_auto"
    } else {
      "explicit"
    }

  result
}
