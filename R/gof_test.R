gof_test_families <- function() {
  sub("^multiplier_bootstrap_", "", setdiff(
    grep("^multiplier_bootstrap_", ls(topenv(environment(gof_test_families))), value = TRUE),
    "multiplier_bootstrap_gof"
  ))
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
  do.call(wrapper, args)
}
