gof_test_families <- function() {
  c("normal", "mvnormal", "restricted_spiked_normal", "normal_sigma_Id",
    "vmf", "vmf_fixed_kappa", "hvmf", "logistic_gaussian",
    "logistic_gaussian_ar1", "cardioid", "small_circle",
    "uniform_beta_mixture", "sunspots_joint_time_space")
}

gof_test <- function(data,
                     h0,
                     theta = NULL,
                     bootstrap = c("fast", "reestimated"),
                     B = 1000,
                     statistics = c("ks", "cvm"),
                     seed = NULL,
                     n_cores = 1,
                     ...) {
  families <- gof_test_families()
  if (missing(h0) || length(h0) != 1L || !h0 %in% families) {
    stop("`h0` must be one of: ", paste(sort(families), collapse = ", "), ".", call. = FALSE)
  }
  bootstrap <- match.arg(bootstrap)
  wrapper <- get(paste0("multiplier_bootstrap_", h0), mode = "function")
  extra <- list(...)

  args <- list(
    data = data,
    null = if (is.null(theta)) list(type = "composite") else list(type = "simple", theta = theta),
    statistics = statistics,
    B = B,
    seed = seed,
    n_cores = n_cores
  )
  if ("bootstrap_method" %in% names(formals(wrapper))) {
    args$bootstrap_method <- if (bootstrap == "fast") "fast_multiplier" else "reestimated"
  } else if (bootstrap == "fast") {
    stop("The fast bootstrap is not available for h0 = \"", h0,
         "\"; use bootstrap = \"reestimated\".", call. = FALSE)
  }
  # The fast bootstrap evaluates KS and CvM together and does not keep the
  # observed process.
  if (bootstrap == "fast" && is.null(extra$keep)) {
    args$keep <- list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE)
  }
  # Under a composite normal null, estimate both the mean and the variance.
  if (h0 == "normal" && is.null(theta) && is.null(extra$unknown_param)) {
    args$unknown_param <- "both"
  }
  args[names(extra)] <- extra
  result <- do.call(wrapper, args)

  used <- result$diagnostics$effective_bootstrap_method %||% args$bootstrap_method %||% "reestimated"
  used <- if (identical(used, "fast_multiplier")) "fast" else "reestimated"
  fallback <- bootstrap == "fast" && used == "reestimated"
  reason <- result$diagnostics$fallback_reason %||% NA_character_
  result$bootstrap <- list(requested = bootstrap, used = used, fallback = fallback,
                           fallback_reason = if (fallback) reason else NA_character_)
  if (fallback) {
    warning("The fast bootstrap could not be used",
            if (!is.na(reason)) paste0(" (", reason, ")"),
            "; the re-estimated bootstrap was used instead.", call. = FALSE)
  }
  result
}
