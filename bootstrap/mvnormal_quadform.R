# EN DUDA (2026-07-26): shared numerical infrastructure for multivariate
# normal distance profiles.  This replaces the logistic-Gaussian-only
# dispatcher in production paths, pending the project owner's final review.
#
# For Q = sum_j lambda_j chi^2_{h_j}(delta_j), CompQuadForm returns Qq=P(Q>q)
# whereas distance profiles require F_Q(q).  All exact calls below therefore
# convert with 1-Qq.  HBE is intentionally an explicit approximation only.

mvnormal_quadform_control <- function(control, name, default, aliases = character()) {
  keys <- c(paste0("mvnormal_quadform_", name), aliases)
  for (key in keys) {
    value <- control[[key]]
    if (!is.null(value)) return(value)
  }
  default
}

mvnormal_quadform_with_label <- function(control, label) {
  if (is.null(control$mvnormal_quadform_mc_label) && is.null(control$logistic_gaussian_quadform_mc_label)) {
    control$mvnormal_quadform_mc_label <- label
  }
  control
}

mvnormal_quadform_validate_scalar <- function(value, default, lower = 0, integer = FALSE) {
  value <- suppressWarnings(as.numeric(value)[1L])
  if (!is.finite(value) || value <= lower) value <- default
  if (integer) value <- as.integer(value)
  value
}

mvnormal_quadform_settings <- function(control = list()) {
  method <- tolower(as.character(mvnormal_quadform_control(
    control, "method", "auto", aliases = "logistic_gaussian_quadform_method"
  ))[1L])
  allowed <- c("auto", "farebrother", "imhof", "davies", "hbe")
  if (!is.finite(match(method, allowed))) {
    stop("`mvnormal_quadform_method` must be one of 'auto', 'farebrother', 'imhof', 'davies', or 'hbe'.")
  }

  settings <- list(
    method = method,
    condition_threshold = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(
        control, "condition_threshold", 1e4,
        aliases = c("mvnormal_quadform_cond_threshold", "logistic_gaussian_ill_conditioned_cond_threshold")
      ), 1e4
    ),
    log_a0_threshold = mvnormal_quadform_validate_scalar(mvnormal_quadform_control(
      control, "log_a0_threshold", log(.Machine$double.xmin),
      aliases = "logistic_gaussian_quadform_log_a0_threshold"
    ), log(.Machine$double.xmin), lower = -Inf),
    farebrother_eps = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "eps", 1e-8, aliases = "logistic_gaussian_quadform_eps"),
      1e-8
    ),
    farebrother_maxit = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "maxit", 100000L, aliases = "logistic_gaussian_quadform_maxit"),
      100000L, lower = 999, integer = TRUE
    ),
    imhof_epsabs = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "imhof_epsabs", 1e-8, aliases = "logistic_gaussian_quadform_imhof_epsabs"),
      1e-8
    ),
    imhof_epsrel = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "imhof_epsrel", 1e-8, aliases = "logistic_gaussian_quadform_imhof_epsrel"),
      1e-8
    ),
    imhof_limit = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "imhof_limit", 20000L, aliases = "logistic_gaussian_quadform_imhof_limit"),
      20000L, lower = 999, integer = TRUE
    ),
    rescue_abserr = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(
        control, "rescue_abserr", 1e-4,
        aliases = c("logistic_gaussian_quadform_rescue_abserr", "logistic_gaussian_quadform_imhof_max_abserr")
      ), 1e-4
    ),
    davies_acc = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "davies_acc", 1e-4, aliases = c("logistic_gaussian_quadform_acc", "logistic_gaussian_quadform_davies_acc")),
      1e-4
    ),
    davies_lim = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "davies_lim", 200000L, aliases = c("logistic_gaussian_quadform_lim", "logistic_gaussian_quadform_davies_lim")),
      200000L, lower = 999, integer = TRUE
    ),
    mc_conf_level = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "mc_conf_level", 0.95, aliases = "logistic_gaussian_quadform_mc_conf_level"),
      0.95
    ),
    mc_abs_error = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "mc_abs_error", 1e-4, aliases = "logistic_gaussian_quadform_mc_abs_error"),
      1e-4
    ),
    mc_batch_size = mvnormal_quadform_validate_scalar(
      mvnormal_quadform_control(control, "mc_batch_size", 1000000L, aliases = "logistic_gaussian_quadform_mc_batch_size"),
      1000000L, lower = 0, integer = TRUE
    ),
    mc_seed = {
      seed <- suppressWarnings(as.integer(mvnormal_quadform_control(
        control, "mc_seed", 20260725L, aliases = "logistic_gaussian_quadform_mc_seed"
      ))[1L])
      if (is.na(seed)) 20260725L else seed
    },
    mc_label = as.character(mvnormal_quadform_control(
      control, "mc_label", "unspecified MVN/LG profile", aliases = "logistic_gaussian_quadform_mc_label"
    ))[1L]
  )
  if (settings$mc_conf_level >= 1) settings$mc_conf_level <- 0.95
  settings
}

mvnormal_quadform_canonicalize <- function(q, lambda, h, delta) {
  q <- as.numeric(q)
  lambda <- as.numeric(lambda)
  h <- as.numeric(h)
  delta <- as.numeric(delta)
  if (!length(q) || any(!is.finite(q))) stop("`q` must be a non-empty finite numeric vector.")
  if (!length(lambda) || length(h) != length(lambda) || length(delta) != length(lambda) ||
      any(!is.finite(lambda)) || any(!is.finite(h)) || any(!is.finite(delta)) ||
      any(lambda <= 0) || any(h <= 0) || any(delta < 0)) {
    stop("`lambda`, `h`, and `delta` must have equal positive finite length, with non-negative `delta`.")
  }
  scale <- max(lambda)
  lambda <- lambda / scale
  q <- q / scale
  ordering <- order(lambda)
  lambda <- lambda[ordering]
  h <- h[ordering]
  delta <- delta[ordering]
  starts <- c(TRUE, lambda[-1L] != lambda[-length(lambda)])
  groups <- cumsum(starts)
  lambda <- lambda[starts]
  h <- as.numeric(rowsum(h, group = groups, reorder = FALSE)[, 1L])
  delta <- as.numeric(rowsum(delta, group = groups, reorder = FALSE)[, 1L])
  lambda_min <- min(lambda)
  list(
    q = q, lambda = lambda, h = h, delta = delta, scale = scale,
    condition_number = max(lambda) / lambda_min,
    log_a0 = 0.5 * sum(h * log(lambda_min / lambda)) - 0.5 * sum(delta)
  )
}

mvnormal_quadform_capture <- function(fun) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  list(value = value, warnings = unique(warnings))
}

# EN DUDA (2026-07-26): one narrow indirection keeps the numerical backends
# testable without changing their production arguments or control flow.
mvnormal_quadform_backend_call <- function(method, ...) {
  do.call(getExportedValue("CompQuadForm", method), list(...))
}

mvnormal_quadform_raw_cdf <- function(captured) {
  if (inherits(captured$value, "error") || is.null(captured$value$Qq)) return(NA_real_)
  raw <- suppressWarnings(as.numeric(1 - captured$value$Qq)[1L])
  if (is.finite(raw)) raw else NA_real_
}

# This projection is valid only when the reported absolute numerical error is
# large enough to cover the excursion outside [0,1].
mvnormal_quadform_safe_clip <- function(raw_probability, abserr) {
  if (length(raw_probability) != 1L || length(abserr) != 1L ||
      !is.finite(raw_probability) || !is.finite(abserr) || abserr < 0 ||
      raw_probability < -abserr || raw_probability > 1 + abserr) {
    return(NA_real_)
  }
  pmin(pmax(raw_probability, 0), 1)
}

mvnormal_quadform_run_farebrother <- function(q, canonical, settings) {
  captured <- mvnormal_quadform_capture(function() mvnormal_quadform_backend_call("farebrother",
    q = q, lambda = canonical$lambda, h = canonical$h, delta = canonical$delta,
    eps = settings$farebrother_eps, maxit = settings$farebrother_maxit, mode = 1L
  ))
  result <- captured$value
  ifault <- if (!inherits(result, "error")) as.integer(result$ifault %||% NA_integer_) else NA_integer_
  raw <- mvnormal_quadform_raw_cdf(captured)
  list(
    probability = if (is.finite(raw) && raw >= 0 && raw <= 1 && identical(ifault, 0L)) raw else NA_real_,
    raw_probability = raw, ifault = ifault, warnings = captured$warnings,
    error = if (inherits(result, "error")) conditionMessage(result) else NA_character_,
    settings = list(eps = settings$farebrother_eps, maxit = settings$farebrother_maxit, mode = 1L)
  )
}

mvnormal_quadform_run_imhof <- function(q, canonical, epsabs, epsrel, limit, accepted_error) {
  captured <- mvnormal_quadform_capture(function() mvnormal_quadform_backend_call("imhof",
    q = q, lambda = canonical$lambda, h = canonical$h, delta = canonical$delta,
    epsabs = epsabs, epsrel = epsrel, limit = limit
  ))
  result <- captured$value
  abserr <- if (!inherits(result, "error")) suppressWarnings(as.numeric(result$abserr)[1L]) else NA_real_
  raw <- mvnormal_quadform_raw_cdf(captured)
  clipped <- mvnormal_quadform_safe_clip(raw, abserr)
  list(
    probability = if (is.finite(abserr) && abserr <= accepted_error) clipped else NA_real_,
    raw_probability = raw, abserr = abserr, warnings = captured$warnings,
    error = if (inherits(result, "error")) conditionMessage(result) else NA_character_,
    settings = list(epsabs = epsabs, epsrel = epsrel, limit = limit, accepted_error = accepted_error)
  )
}

mvnormal_quadform_run_davies <- function(q, canonical, settings) {
  captured <- mvnormal_quadform_capture(function() mvnormal_quadform_backend_call("davies",
    q = q, lambda = canonical$lambda, h = canonical$h, delta = canonical$delta,
    acc = settings$davies_acc, lim = settings$davies_lim
  ))
  result <- captured$value
  ifault <- if (!inherits(result, "error")) as.integer(result$ifault %||% NA_integer_) else NA_integer_
  raw <- mvnormal_quadform_raw_cdf(captured)
  list(
    probability = if (identical(ifault, 0L) && is.finite(raw) && raw >= 0 && raw <= 1) raw else NA_real_,
    raw_probability = raw, ifault = ifault, warnings = captured$warnings,
    error = if (inherits(result, "error")) conditionMessage(result) else NA_character_,
    settings = list(acc = settings$davies_acc, lim = settings$davies_lim)
  )
}

mvnormal_quadform_cp_interval <- function(successes, n, conf_level) {
  alpha <- 1 - conf_level
  lower <- ifelse(successes == 0L, 0, stats::qbeta(alpha / 2, successes, n - successes + 1))
  upper <- ifelse(successes == n, 1, stats::qbeta(1 - alpha / 2, successes + 1, n - successes))
  cbind(lower = lower, upper = upper)
}

# Serial, shared simulation for every unresolved radius of one centre.  The
# returned probability is the empirical mean, not the midpoint of this CI.
mvnormal_quadform_mc_cdf <- function(q, canonical, settings) {
  q <- as.numeric(q)
  if (!length(q) || any(!is.finite(q))) stop("MC requires finite thresholds.")
  start_time <- proc.time()[["elapsed"]]
  message(sprintf(
    "[mvnormal_quadform MC] start | %s | tolerance=%g | thresholds=%d",
    settings$mc_label, settings$mc_abs_error, length(q)
  ))
  completed <- FALSE
  n_total <- 0L
  successes <- integer(length(q))
  interval <- matrix(NA_real_, nrow = length(q), ncol = 2L,
                     dimnames = list(NULL, c("lower", "upper")))
  on.exit({
    if (!completed) {
      elapsed <- proc.time()[["elapsed"]] - start_time
      estimate_text <- if (n_total > 0L) {
        paste(formatC(successes / n_total, digits = 7L, format = "fg"), collapse = ", ")
      } else {
        "unavailable"
      }
      interval_text <- if (n_total > 0L && all(is.finite(interval))) {
        paste(apply(interval, 1L, function(x) paste(formatC(x, digits = 7L, format = "fg"), collapse = ",")), collapse = "; ")
      } else {
        "unavailable"
      }
      message(sprintf(
        "[mvnormal_quadform MC] interrupted | %s | tolerance=%g | B=%d | estimate=[%s] | CI95=[%s] | elapsed=%.2fs",
        settings$mc_label, settings$mc_abs_error, n_total, estimate_text, interval_text, elapsed
      ))
    }
  }, add = TRUE)

  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_rng) {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (is.na(settings$mc_seed)) stop("`mvnormal_quadform_mc_seed` must be a finite integer.")
  set.seed(settings$mc_seed)

  repeat {
    batch <- settings$mc_batch_size
    draws <- numeric(batch)
    for (j in seq_along(canonical$lambda)) {
      draws <- draws + canonical$lambda[[j]] * stats::rchisq(
        batch, df = canonical$h[[j]], ncp = canonical$delta[[j]]
      )
    }
    sorted_draws <- sort(draws)
    successes <- successes + findInterval(q, sorted_draws)
    n_total <- n_total + batch
    interval <- mvnormal_quadform_cp_interval(successes, n_total, settings$mc_conf_level)
    if (all((interval[, "upper"] - interval[, "lower"]) / 2 <= settings$mc_abs_error)) break
  }

  probability <- successes / n_total
  elapsed <- proc.time()[["elapsed"]] - start_time
  completed <- TRUE
  message(sprintf(
    "[mvnormal_quadform MC] complete | %s | B=%d | estimate=[%s] | CI95=[%s] | elapsed=%.2fs",
    settings$mc_label, n_total,
    paste(formatC(probability, digits = 7L, format = "fg"), collapse = ", "),
    paste(apply(interval, 1L, function(x) paste(formatC(x, digits = 7L, format = "fg"), collapse = ",")), collapse = "; "),
    elapsed
  ))
  list(probability = probability, interval = interval, n = n_total, elapsed = elapsed)
}

mvnormal_quadform_one_auto <- function(q, canonical, settings) {
  pre_screen <- canonical$condition_number > settings$condition_threshold ||
    canonical$log_a0 < settings$log_a0_threshold
  attempts <- list()
  if (!pre_screen) {
    attempts$farebrother <- mvnormal_quadform_run_farebrother(q, canonical, settings)
    if (is.finite(attempts$farebrother$probability)) {
      return(list(probability = attempts$farebrother$probability, status = "farebrother", attempts = attempts,
                  pre_screen = FALSE))
    }
  }

  attempts$imhof_usual <- mvnormal_quadform_run_imhof(
    q, canonical, settings$imhof_epsabs, settings$imhof_epsrel,
    settings$imhof_limit, settings$rescue_abserr
  )
  if (is.finite(attempts$imhof_usual$probability)) {
    return(list(probability = attempts$imhof_usual$probability, status = "imhof_usual", attempts = attempts,
                pre_screen = pre_screen))
  }

  attempts$imhof_rescue <- mvnormal_quadform_run_imhof(
    q, canonical, settings$rescue_abserr, settings$rescue_abserr,
    settings$imhof_limit, settings$rescue_abserr
  )
  if (is.finite(attempts$imhof_rescue$probability)) {
    return(list(probability = attempts$imhof_rescue$probability, status = "imhof_rescue", attempts = attempts,
                pre_screen = pre_screen))
  }

  attempts$davies <- mvnormal_quadform_run_davies(q, canonical, settings)
  if (is.finite(attempts$davies$probability)) {
    return(list(probability = attempts$davies$probability, status = "davies", attempts = attempts,
                pre_screen = pre_screen))
  }
  list(probability = NA_real_, status = "needs_mc", attempts = attempts, pre_screen = pre_screen)
}

mvnormal_quadform_one_explicit <- function(q, canonical, settings) {
  method <- settings$method
  if (identical(method, "hbe")) {
    probability <- sphunif::p_wschisq(
      x = q, weights = canonical$lambda, dfs = canonical$h, ncps = canonical$delta, method = "HBE"
    )
    if (!is.finite(probability) || probability < 0 || probability > 1) {
      stop("Explicit HBE returned an invalid probability.")
    }
    return(list(probability = probability, status = "hbe", attempts = list(), pre_screen = NA))
  }
  if (identical(method, "farebrother")) {
    attempt <- mvnormal_quadform_run_farebrother(q, canonical, settings)
    if (!is.finite(attempt$probability)) stop("Explicit Farebrother evaluation failed.")
    return(list(probability = attempt$probability, status = "farebrother", attempts = list(farebrother = attempt), pre_screen = NA))
  }
  if (identical(method, "imhof")) {
    attempt <- mvnormal_quadform_run_imhof(q, canonical, settings$imhof_epsabs, settings$imhof_epsrel,
                                            settings$imhof_limit, settings$rescue_abserr)
    if (!is.finite(attempt$probability)) stop("Explicit Imhof evaluation did not reach the requested accuracy.")
    return(list(probability = attempt$probability, status = "imhof", attempts = list(imhof = attempt), pre_screen = NA))
  }
  attempt <- mvnormal_quadform_run_davies(q, canonical, settings)
  if (!is.finite(attempt$probability)) stop("Explicit Davies evaluation failed (ifault must be zero).")
  list(probability = attempt$probability, status = "davies", attempts = list(davies = attempt), pre_screen = NA)
}

# Vector thresholds share the same canonical Q.  If exact methods fail at more
# than one radius, the terminal MC route is called once for all such radii.
mvnormal_quadform_cdf_diagnostics <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)), control = list()) {
  canonical <- mvnormal_quadform_canonicalize(q, lambda, h, delta)
  settings <- mvnormal_quadform_settings(control)
  results <- vector("list", length(canonical$q))
  for (i in seq_along(canonical$q)) {
    threshold <- canonical$q[[i]]
    if (threshold <= 0) {
      results[[i]] <- list(probability = 0, status = "analytic_below_support", attempts = list(), pre_screen = NA)
    } else if (length(canonical$lambda) == 1L && !identical(settings$method, "hbe")) {
      results[[i]] <- list(
        probability = stats::pchisq(threshold / canonical$lambda, df = canonical$h, ncp = canonical$delta),
        status = "analytic_noncentral_chisquare", attempts = list(), pre_screen = NA
      )
    } else if (identical(settings$method, "auto")) {
      results[[i]] <- mvnormal_quadform_one_auto(threshold, canonical, settings)
    } else {
      results[[i]] <- mvnormal_quadform_one_explicit(threshold, canonical, settings)
    }
  }
  need_mc <- which(vapply(results, function(x) identical(x$status, "needs_mc"), logical(1)))
  if (length(need_mc)) {
    mc <- mvnormal_quadform_mc_cdf(canonical$q[need_mc], canonical, settings)
    for (j in seq_along(need_mc)) {
      i <- need_mc[[j]]
      results[[i]]$probability <- mc$probability[[j]]
      results[[i]]$status <- "mc"
      results[[i]]$mc <- mc[c("interval", "n", "elapsed")]
    }
  }
  list(
    probability = vapply(results, `[[`, numeric(1), "probability"),
    details = results,
    canonical = canonical,
    settings = settings
  )
}

mvnormal_quadform_cdf <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)), control = list()) {
  mvnormal_quadform_cdf_diagnostics(q, lambda, h, delta, control)$probability
}

# Backwards-compatible historical name; it returns the CDF despite the older
# name saying 'tail probability'.
logistic_gaussian_quadform_tail_probability <- function(q, lambda, h, delta, control = list()) {
  probabilities <- mvnormal_quadform_cdf(q, lambda, h, delta, control)
  if (length(probabilities) == 1L) probabilities[[1L]] else probabilities
}

evaluate_mvnorm_distance_profile <- function(shift, t_values, eigenvalues_full,
                                             eigenvectors_full, positive_idx,
                                             tol = 1e-12, control = list()) {
  shift <- as.numeric(shift)
  t_values <- as.numeric(t_values)
  if (any(!is.finite(shift)) || any(!is.finite(t_values))) stop("`shift` and `t_values` must be finite.")
  nu <- as.vector(crossprod(eigenvectors_full, shift))
  null_const <- sum(nu[!positive_idx]^2)
  threshold_sq <- t_values^2 - null_const
  output <- numeric(length(t_values))
  positive <- threshold_sq >= -tol
  if (!any(positive_idx)) {
    output[positive] <- 1
    return(output)
  }
  threshold_sq[threshold_sq < 0] <- 0
  if (any(positive)) {
    lambda <- as.numeric(eigenvalues_full[positive_idx])
    delta <- nu[positive_idx]^2 / lambda
    output[positive] <- mvnormal_quadform_cdf(threshold_sq[positive], lambda, rep.int(1, length(lambda)), delta, control)
  }
  pmin(pmax(output, 0), 1)
}

evaluate_mvnorm_distance_profile_matrix <- function(shift_matrix, t_matrix,
                                                    eigenvalues_full, eigenvectors_full,
                                                    positive_idx, tol = 1e-12,
                                                    control = list()) {
  shift_matrix <- as.matrix(shift_matrix)
  t_matrix <- as.matrix(t_matrix)
  if (nrow(shift_matrix) == 0L || ncol(shift_matrix) == 0L || nrow(t_matrix) != nrow(shift_matrix) ||
      any(!is.finite(shift_matrix)) || any(!is.finite(t_matrix))) {
    stop("`shift_matrix` and `t_matrix` must be finite, non-empty and have equal row counts.")
  }
  nu <- shift_matrix %*% eigenvectors_full
  null_const <- if (any(!positive_idx)) rowSums(nu[, !positive_idx, drop = FALSE]^2) else rep.int(0, nrow(nu))
  threshold_sq <- t_matrix^2 - matrix(null_const, nrow(t_matrix), ncol(t_matrix))
  positive <- threshold_sq >= -tol
  output <- matrix(0, nrow(t_matrix), ncol(t_matrix))
  if (!any(positive_idx)) {
    output[positive] <- 1
    return(output)
  }
  threshold_sq[threshold_sq < 0] <- 0
  lambda <- as.numeric(eigenvalues_full[positive_idx])
  for (i in seq_len(nrow(t_matrix))) {
    idx <- which(positive[i, ])
    if (!length(idx)) next
    delta <- nu[i, positive_idx]^2 / lambda
    output[i, idx] <- mvnormal_quadform_cdf(threshold_sq[i, idx], lambda, rep.int(1, length(lambda)), delta, control)
  }
  pmin(pmax(output, 0), 1)
}
