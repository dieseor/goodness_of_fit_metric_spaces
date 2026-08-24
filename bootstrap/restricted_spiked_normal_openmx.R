# Independent OpenMx reference fit for the restricted-spiked Gaussian model.
# This file is intentionally not sourced by model_specs.R: OpenMx is a test
# dependency, whereas production simulations use the profiled custom MLE.

restricted_spiked_normal_openmx_available <- function() {
  requireNamespace("OpenMx", quietly = TRUE)
}

fit_restricted_spiked_normal_openmx <- function(data,
                                                start_theta = NULL,
                                                start_lambda = NULL,
                                                silent = TRUE,
                                                control = list()) {
  if (!restricted_spiked_normal_openmx_available()) {
    stop("OpenMx is required for the independent restricted-spiked MLE reference fit.")
  }
  x <- normalize_restricted_spiked_normal_data(data, control)
  d <- ncol(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(d))

  custom_start <- fit_restricted_spiked_normal_theta(
    x, null = list(type = "composite"), control = control
  )
  start_theta <- as.numeric(start_theta %||% custom_start$theta)
  start_lambda <- as.numeric(start_lambda %||% custom_start$lambda)[1L]
  if (length(start_theta) != d || any(!is.finite(start_theta)) ||
      !is.finite(start_lambda) || start_lambda <= 0 || sum(start_theta^2) == 0) {
    stop("OpenMx restricted-spiked starts must have a non-zero finite theta and strictly positive lambda.")
  }

  theta_labels <- paste0("theta_", seq_len(d))
  model <- OpenMx::mxModel(
    "restricted_spiked_normal_reference",
    OpenMx::mxMatrix(
      type = "Full", nrow = d, ncol = 1L, free = TRUE,
      values = matrix(start_theta, ncol = 1L), labels = theta_labels,
      name = "theta"
    ),
    OpenMx::mxMatrix(
      type = "Full", nrow = 1L, ncol = 1L, free = TRUE,
      # OpenMx supports a closed lower bound only.  A fitted value numerically
      # at zero is rejected below by the shared strict model validator.
      values = matrix(start_lambda, nrow = 1L), lbound = matrix(0, nrow = 1L),
      labels = "lambda_parameter", name = "lambda"
    ),
    OpenMx::mxMatrix(type = "Iden", nrow = d, ncol = d, name = "I"),
    OpenMx::mxAlgebra(theta %*% t(theta) / sum(theta * theta), name = "uu"),
    OpenMx::mxAlgebra(I + lambda[1, 1] * uu, name = "Sigma"),
    OpenMx::mxAlgebra(t(theta), name = "mean_row"),
    OpenMx::mxData(observed = x, type = "raw"),
    OpenMx::mxExpectationNormal(
      covariance = "Sigma", means = "mean_row", dimnames = colnames(x)
    ),
    OpenMx::mxFitFunctionML()
  )
  fit <- OpenMx::mxRun(model, silent = isTRUE(silent))
  theta_hat <- as.numeric(fit$theta@values)
  lambda_hat <- as.numeric(fit$lambda@values)[1L]
  normalized <- normalize_restricted_spiked_normal_theta(
    list(theta = theta_hat, lambda = lambda_hat), ambient_dim = d, control = control
  )
  loglik_manual <- nrow(x) * restricted_spiked_normal_loglik(x, normalized, control = control)
  loglik_openmx <- -0.5 * as.numeric(fit$output$fit)
  c(
    normalized,
    list(
      loglik = loglik_manual,
      loglik_openmx = loglik_openmx,
      fit_status_code = as.integer(fit$output$status$code %||% NA_integer_),
      fit_status_message = as.character(fit$output$status$status %||% NA_character_),
      openmx_fit = fit,
      start_theta = start_theta,
      start_lambda = start_lambda
    )
  )
}

fit_restricted_spiked_normal_openmx_multistart <- function(data,
                                                           starts,
                                                           silent = TRUE,
                                                           control = list()) {
  if (!is.list(starts) || !length(starts)) {
    stop("`starts` must be a non-empty list of theta/lambda start lists.")
  }
  fits <- lapply(starts, function(start) {
    tryCatch(
      fit_restricted_spiked_normal_openmx(
        data = data,
        start_theta = start$theta,
        start_lambda = start$lambda,
        silent = silent,
        control = control
      ),
      error = function(error) list(error = conditionMessage(error))
    )
  })
  successful <- which(!vapply(fits, function(fit) !is.null(fit$error), logical(1L)))
  if (!length(successful)) {
    stop("Every OpenMx restricted-spiked reference fit failed.")
  }
  loglik <- vapply(fits[successful], `[[`, numeric(1L), "loglik")
  list(
    fits = fits,
    successful = successful,
    best = fits[[successful[[which.max(loglik)]]]],
    loglik_range = range(loglik),
    all_converged = all(vapply(fits[successful], function(fit) {
      identical(fit$fit_status_code, 0L)
    }, logical(1L)))
  )
}
