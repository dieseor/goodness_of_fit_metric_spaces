source("bootstrap/calibration_study.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

normalize_mu <- function(mu) {
  mu <- as.numeric(mu)
  mu / sqrt(sum(mu^2))
}

angle_error_deg <- function(mu_hat, mu_true) {
  mu_hat <- normalize_mu(mu_hat)
  mu_true <- normalize_mu(mu_true)
  acos(pmax(pmin(sum(mu_hat * mu_true), 1), -1)) * 180 / pi
}

format_mu_hat <- function(mu) {
  paste(sprintf("%.10f", as.numeric(mu)), collapse = ",")
}

safe_loglik <- function(theta, x, prob_weights) {
  if (is.null(theta) || any(!is.finite(unlist(theta[c("mu", "kappa", "psi")]))) ) {
    return(NA_real_)
  }
  jp_weighted_loglik_s2_prepared(
    mu = theta$mu,
    kappa = theta$kappa,
    psi = theta$psi,
    x = x,
    prob_weights = prob_weights
  )
}

make_mu_starts <- function(x,
                           resultant_tol = 0.08,
                           n_extra_when_uniform = 12L,
                           n_extra_when_concentrated = 4L) {
  starts <- list(c(0, 0, 1), c(0, 0, -1))

  resultant <- colMeans(x)
  resultant_norm <- sqrt(sum(resultant^2))
  if (is.finite(resultant_norm) && resultant_norm > resultant_tol) {
    starts[[length(starts) + 1L]] <- resultant / resultant_norm
  }

  n_extra <- if (is.finite(resultant_norm) && resultant_norm <= resultant_tol) {
    as.integer(n_extra_when_uniform)
  } else {
    as.integer(n_extra_when_concentrated)
  }

  lattice <- generate_canonical_lattice(n = n_extra, dim = 3)
  for (i in seq_len(nrow(lattice))) {
    starts[[length(starts) + 1L]] <- lattice[i, ]
  }

  start_matrix <- do.call(rbind, starts)
  start_matrix <- start_matrix / sqrt(rowSums(start_matrix^2))

  # Remove near-duplicates to keep optimization costs controlled.
  keep <- !duplicated(round(start_matrix, 6))
  start_matrix[keep, , drop = FALSE]
}

optimize_mu_given_kappa_psi <- function(x,
                                        prob_weights,
                                        kappa,
                                        psi,
                                        mu_starts,
                                        maxit_mu = 200L,
                                        reltol_mu = 1e-8,
                                        method_mu = "Nelder-Mead") {
  best <- list(loglik = -Inf, mu = NULL, status = "no_finite_start")

  objective <- function(raw_mu) {
    mu <- jp_mu_s2_from_raw(raw_mu[[1L]], raw_mu[[2L]])
    value <- jp_weighted_loglik_s2_prepared(
      mu = mu,
      kappa = kappa,
      psi = psi,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      return(Inf)
    }
    -value
  }

  for (i in seq_len(nrow(mu_starts))) {
    raw_start <- jp_mu_s2_to_raw(mu_starts[i, ])
    fit <- try(
      stats::optim(
        par = raw_start,
        fn = objective,
        method = method_mu,
        control = list(maxit = as.integer(maxit_mu), reltol = as.numeric(reltol_mu))
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || !is.list(fit) || !is.finite(fit$value)) {
      next
    }

    mu_hat <- jp_mu_s2_from_raw(fit$par[[1L]], fit$par[[2L]])
    loglik <- -fit$value
    if (is.finite(loglik) && loglik > best$loglik) {
      best <- list(loglik = loglik, mu = mu_hat, status = ifelse(fit$convergence == 0L, "ok", "optim_nonzero_conv"))
    }
  }

  best
}

profile_branch_candidates <- function(x,
                                      prob_weights,
                                      psi_sign,
                                      kappa_grid,
                                      psi_abs_grid,
                                      max_abs_kappa_psi,
                                      mu_starts,
                                      maxit_mu,
                                      reltol_mu,
                                      method_mu = "Nelder-Mead") {
  candidates <- list()
  idx <- 1L

  for (kappa in kappa_grid) {
    if (!is.finite(kappa) || kappa < 0) {
      next
    }

    for (psi_abs in psi_abs_grid) {
      if (!is.finite(psi_abs) || psi_abs <= 0) {
        next
      }

      psi <- psi_sign * psi_abs
      if (abs(kappa * psi) > max_abs_kappa_psi) {
        next
      }

      mu_fit <- optimize_mu_given_kappa_psi(
        x = x,
        prob_weights = prob_weights,
        kappa = kappa,
        psi = psi,
        mu_starts = mu_starts,
        maxit_mu = maxit_mu,
        reltol_mu = reltol_mu,
        method_mu = method_mu
      )

      if (!is.finite(mu_fit$loglik)) {
        next
      }

      candidates[[idx]] <- list(
        mu = mu_fit$mu,
        kappa = kappa,
        psi = psi,
        loglik = mu_fit$loglik,
        branch = ifelse(psi_sign > 0, "positive", "negative"),
        source = "profile_kappa_psi_grid"
      )
      idx <- idx + 1L
    }
  }

  if (length(candidates) == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind, lapply(candidates, function(candidate) {
    data.frame(
      mu1 = candidate$mu[[1L]],
      mu2 = candidate$mu[[2L]],
      mu3 = candidate$mu[[3L]],
      kappa = candidate$kappa,
      psi = candidate$psi,
      loglik = candidate$loglik,
      branch = candidate$branch,
      source = candidate$source,
      stringsAsFactors = FALSE
    )
  }))

  out[order(-out$loglik), , drop = FALSE]
}

refine_profile_candidate <- function(candidate,
                                     x,
                                     prob_weights,
                                     psi_min,
                                     max_abs_kappa_psi,
                                     maxit_refine = 500L,
                                     reltol_refine = 1e-10) {
  sign_branch <- ifelse(candidate$psi >= 0, 1L, -1L)
  psi_abs <- abs(candidate$psi)
  if (!is.finite(psi_abs) || psi_abs <= psi_min) {
    return(list(theta = jp_weighted_vmf_mle_s2(x, prob_weights), status = "refine_near_zero_to_vmf"))
  }

  a_start <- abs(candidate$kappa * candidate$psi)
  raw_start <- c(
    jp_mu_s2_to_raw(c(candidate$mu1, candidate$mu2, candidate$mu3)),
    log(max(a_start, 1e-6)),
    log(max(psi_abs - psi_min, 1e-8))
  )

  objective <- function(raw) {
    if (!all(is.finite(raw)) || raw[[3L]] > 700 || raw[[4L]] > 700) {
      return(Inf)
    }
    a_value <- exp(raw[[3L]])
    if (!is.finite(a_value) || a_value > max_abs_kappa_psi) {
      return(Inf)
    }
    jp_neg_weighted_loglik_s2_raw(
      raw = raw,
      sign_branch = sign_branch,
      x = x,
      prob_weights = prob_weights,
      psi_min = psi_min,
      max_abs_kappa_psi = max_abs_kappa_psi
    )
  }

  fit_nm <- try(
    stats::optim(
      par = raw_start,
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = as.integer(maxit_refine), reltol = as.numeric(reltol_refine))
    ),
    silent = TRUE
  )

  if (inherits(fit_nm, "try-error") || !is.list(fit_nm) || !is.finite(fit_nm$value)) {
    return(list(theta = NULL, status = "refine_nm_failed"))
  }

  fit_bfgs <- try(
    stats::optim(
      par = fit_nm$par,
      fn = objective,
      method = "BFGS",
      control = list(maxit = as.integer(maxit_refine), reltol = as.numeric(reltol_refine))
    ),
    silent = TRUE
  )

  best_fit <- if (!inherits(fit_bfgs, "try-error") && is.list(fit_bfgs) && is.finite(fit_bfgs$value)) fit_bfgs else fit_nm
  params <- jp_params_s2_from_raw(best_fit$par, sign_branch = sign_branch, psi_min = psi_min)
  if (is.null(params)) {
    return(list(theta = NULL, status = "refine_params_invalid"))
  }

  if (!is.finite(abs(params$kappa * params$psi)) || abs(params$kappa * params$psi) > max_abs_kappa_psi) {
    return(list(theta = NULL, status = "refine_cap_violation"))
  }

  theta <- list(mu = params$mu, kappa = params$kappa, psi = params$psi)
  list(theta = theta, status = ifelse(best_fit$convergence == 0L, "ok", "refine_nonzero_conv"))
}

compute_profile_context <- function(x,
                                    prob_weights,
                                    theta_true,
                                    robust_control) {
  vmf_candidate <- jp_weighted_vmf_mle_s2(x, prob_weights)
  loglik_vmf_candidate <- safe_loglik(vmf_candidate, x, prob_weights)
  loglik_true_theta <- jp_weighted_loglik_s2_prepared(
    mu = theta_true$mu,
    kappa = theta_true$kappa,
    psi = theta_true$psi,
    x = x,
    prob_weights = prob_weights
  )

  branch_candidates_pos <- profile_branch_candidates(
    x = x,
    prob_weights = prob_weights,
    psi_sign = 1L,
    kappa_grid = robust_control$kappa_grid,
    psi_abs_grid = robust_control$psi_abs_grid,
    max_abs_kappa_psi = robust_control$max_abs_kappa_psi,
    mu_starts = robust_control$mu_starts,
    maxit_mu = robust_control$maxit_mu,
    reltol_mu = robust_control$reltol_mu,
    method_mu = robust_control$method_mu
  )

  branch_candidates_neg <- profile_branch_candidates(
    x = x,
    prob_weights = prob_weights,
    psi_sign = -1L,
    kappa_grid = robust_control$kappa_grid,
    psi_abs_grid = robust_control$psi_abs_grid,
    max_abs_kappa_psi = robust_control$max_abs_kappa_psi,
    mu_starts = robust_control$mu_starts,
    maxit_mu = robust_control$maxit_mu,
    reltol_mu = robust_control$reltol_mu,
    method_mu = robust_control$method_mu
  )

  best_positive <- if (nrow(branch_candidates_pos) > 0L) {
    branch_candidates_pos[which.max(branch_candidates_pos$loglik), , drop = FALSE]
  } else {
    data.frame()
  }

  best_negative <- if (nrow(branch_candidates_neg) > 0L) {
    branch_candidates_neg[which.max(branch_candidates_neg$loglik), , drop = FALSE]
  } else {
    data.frame()
  }

  list(
    vmf_candidate = vmf_candidate,
    loglik_vmf_candidate = loglik_vmf_candidate,
    loglik_true_theta = loglik_true_theta,
    branch_candidates_pos = branch_candidates_pos,
    branch_candidates_neg = branch_candidates_neg,
    loglik_best_positive_branch = if (nrow(best_positive) > 0L) best_positive$loglik[[1L]] else NA_real_,
    loglik_best_negative_branch = if (nrow(best_negative) > 0L) best_negative$loglik[[1L]] else NA_real_
  )
}

run_observed_method <- function(method,
                                x,
                                theta_true,
                                scenario_id,
                                n,
                                outer_id,
                                prob_weights,
                                fast_control,
                                precise_control,
                                robust_control,
                                profile_context) {
  t0 <- proc.time()[["elapsed"]]

  vmf_candidate <- profile_context$vmf_candidate
  loglik_vmf_candidate <- profile_context$loglik_vmf_candidate
  loglik_true_theta <- profile_context$loglik_true_theta
  branch_candidates_pos <- profile_context$branch_candidates_pos
  branch_candidates_neg <- profile_context$branch_candidates_neg
  loglik_best_positive_branch <- profile_context$loglik_best_positive_branch
  loglik_best_negative_branch <- profile_context$loglik_best_negative_branch

  theta_hat <- NULL
  convergence_status <- "ok"
  loglik_profile_refine_final <- NA_real_
  profile_refine_selected_branch <- NA_character_
  profile_refine_selected_candidate <- NA_character_
  loglik_refined_positive_best <- NA_real_
  loglik_refined_negative_best <- NA_real_
  kappa_refined_positive_best <- NA_real_
  psi_refined_positive_best <- NA_real_
  kappa_refined_negative_best <- NA_real_
  psi_refined_negative_best <- NA_real_
  delta_loglik_refined_neg_minus_pos <- NA_real_

  if (identical(method, "current_fast")) {
    fit <- try(jp_mle_s2_weighted(data = x, weights = prob_weights, control = fast_control), silent = TRUE)
    if (inherits(fit, "try-error")) {
      convergence_status <- paste("error", as.character(fit))
    } else {
      theta_hat <- fit
    }
  } else if (identical(method, "current_precise")) {
    fit <- try(jp_mle_s2_weighted(data = x, weights = prob_weights, control = precise_control), silent = TRUE)
    if (inherits(fit, "try-error")) {
      convergence_status <- paste("error", as.character(fit))
    } else {
      theta_hat <- fit
    }
  } else if (identical(method, "profile_refine")) {
    candidates_to_refine <- data.frame()
    if (nrow(branch_candidates_pos) > 0L) {
      candidates_to_refine <- rbind(candidates_to_refine, head(branch_candidates_pos[order(-branch_candidates_pos$loglik), , drop = FALSE], robust_control$top_k_per_branch))
    }
    if (nrow(branch_candidates_neg) > 0L) {
      candidates_to_refine <- rbind(candidates_to_refine, head(branch_candidates_neg[order(-branch_candidates_neg$loglik), , drop = FALSE], robust_control$top_k_per_branch))
    }

    refined <- list()
    refined_status <- character(0)
    if (nrow(candidates_to_refine) > 0L) {
      for (i in seq_len(nrow(candidates_to_refine))) {
        candidate <- candidates_to_refine[i, , drop = FALSE]
        candidate_branch <- as.character(candidate$branch[[1L]])
        refined_fit <- refine_profile_candidate(
          candidate = candidate,
          x = x,
          prob_weights = prob_weights,
          psi_min = robust_control$psi_min,
          max_abs_kappa_psi = robust_control$max_abs_kappa_psi,
          maxit_refine = robust_control$maxit_refine,
          reltol_refine = robust_control$reltol_refine
        )
        refined_status <- c(refined_status, refined_fit$status)
        if (!is.null(refined_fit$theta)) {
          refined[[length(refined) + 1L]] <- list(
            theta = refined_fit$theta,
            source = paste0("refined_", candidate_branch),
            selected_branch = candidate_branch
          )
        }
      }
    }

    refined[[length(refined) + 1L]] <- list(
      theta = vmf_candidate,
      source = "vMF",
      selected_branch = "vMF"
    )

    if (length(refined) == 0L) {
      convergence_status <- "profile_no_candidates"
    } else {
      refined_loglik <- vapply(refined, function(entry) safe_loglik(entry$theta, x, prob_weights), numeric(1))

      positive_idx <- which(vapply(refined, function(entry) identical(entry$selected_branch, "positive"), logical(1)))
      negative_idx <- which(vapply(refined, function(entry) identical(entry$selected_branch, "negative"), logical(1)))

      if (length(positive_idx) > 0L && any(is.finite(refined_loglik[positive_idx]))) {
        best_pos_idx <- positive_idx[which.max(refined_loglik[positive_idx])]
        loglik_refined_positive_best <- refined_loglik[[best_pos_idx]]
        kappa_refined_positive_best <- refined[[best_pos_idx]]$theta$kappa
        psi_refined_positive_best <- refined[[best_pos_idx]]$theta$psi
      }

      if (length(negative_idx) > 0L && any(is.finite(refined_loglik[negative_idx]))) {
        best_neg_idx <- negative_idx[which.max(refined_loglik[negative_idx])]
        loglik_refined_negative_best <- refined_loglik[[best_neg_idx]]
        kappa_refined_negative_best <- refined[[best_neg_idx]]$theta$kappa
        psi_refined_negative_best <- refined[[best_neg_idx]]$theta$psi
      }

      if (is.finite(loglik_refined_positive_best) && is.finite(loglik_refined_negative_best)) {
        delta_loglik_refined_neg_minus_pos <- loglik_refined_negative_best - loglik_refined_positive_best
      }

      if (all(!is.finite(refined_loglik))) {
        convergence_status <- "profile_no_finite_refined"
      } else {
        best_idx <- which.max(refined_loglik)
        theta_hat <- refined[[best_idx]]$theta
        loglik_profile_refine_final <- refined_loglik[[best_idx]]
        profile_refine_selected_branch <- refined[[best_idx]]$selected_branch
        profile_refine_selected_candidate <- refined[[best_idx]]$source
        convergence_status <- paste(c("ok", unique(refined_status)), collapse = ";")
      }
    }
  } else {
    stop(sprintf("Unsupported method: %s", method))
  }

  elapsed_seconds <- proc.time()[["elapsed"]] - t0
  loglik_hat <- safe_loglik(theta_hat, x, prob_weights)

  candidate_names <- c("vMF", "best_positive_branch", "best_negative_branch")
  candidate_logliks <- c(loglik_vmf_candidate, loglik_best_positive_branch, loglik_best_negative_branch)
  if (identical(method, "profile_refine")) {
    candidate_names <- c(candidate_names, "final_profile_refine")
    candidate_logliks <- c(candidate_logliks, loglik_profile_refine_final)
  } else {
    candidate_names <- c(candidate_names, "current_theta_hat")
    candidate_logliks <- c(candidate_logliks, loglik_hat)
  }

  method_vs_candidates <- data.frame(
    candidate = candidate_names,
    loglik = candidate_logliks,
    stringsAsFactors = FALSE
  )
  method_vs_candidates <- method_vs_candidates[is.finite(method_vs_candidates$loglik), , drop = FALSE]
  winner <- if (nrow(method_vs_candidates) > 0L) method_vs_candidates$candidate[[which.max(method_vs_candidates$loglik)]] else NA_character_

  if (is.null(theta_hat)) {
    return(data.frame(
      scenario = scenario_id,
      n = n,
      outer_id = outer_id,
      method = method,
      mu_hat = NA_character_,
      mu1 = NA_real_,
      mu2 = NA_real_,
      mu3 = NA_real_,
      kappa_hat = NA_real_,
      psi_hat = NA_real_,
      loglik_hat = NA_real_,
      loglik_true_theta = loglik_true_theta,
      loglik_vMF_candidate = loglik_vmf_candidate,
      loglik_best_positive_branch = loglik_best_positive_branch,
      loglik_best_negative_branch = loglik_best_negative_branch,
      loglik_profile_refine_final = if (identical(method, "profile_refine")) loglik_profile_refine_final else NA_real_,
      loglik_refined_positive_best = if (identical(method, "profile_refine")) loglik_refined_positive_best else NA_real_,
      loglik_refined_negative_best = if (identical(method, "profile_refine")) loglik_refined_negative_best else NA_real_,
      kappa_refined_positive_best = if (identical(method, "profile_refine")) kappa_refined_positive_best else NA_real_,
      psi_refined_positive_best = if (identical(method, "profile_refine")) psi_refined_positive_best else NA_real_,
      kappa_refined_negative_best = if (identical(method, "profile_refine")) kappa_refined_negative_best else NA_real_,
      psi_refined_negative_best = if (identical(method, "profile_refine")) psi_refined_negative_best else NA_real_,
      delta_loglik_refined_neg_minus_pos = if (identical(method, "profile_refine")) delta_loglik_refined_neg_minus_pos else NA_real_,
      mu_angle_error_deg = NA_real_,
      abs_kappa_error = NA_real_,
      abs_psi_error = NA_real_,
      psi_sign_correct = NA,
      elapsed_seconds = elapsed_seconds,
      convergence_status = convergence_status,
      winner_among_candidates = winner,
      profile_refine_selected_branch = if (identical(method, "profile_refine")) profile_refine_selected_branch else NA_character_,
      profile_refine_selected_candidate = if (identical(method, "profile_refine")) profile_refine_selected_candidate else NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    scenario = scenario_id,
    n = n,
    outer_id = outer_id,
    method = method,
    mu_hat = format_mu_hat(theta_hat$mu),
    mu1 = as.numeric(theta_hat$mu[[1L]]),
    mu2 = as.numeric(theta_hat$mu[[2L]]),
    mu3 = as.numeric(theta_hat$mu[[3L]]),
    kappa_hat = as.numeric(theta_hat$kappa),
    psi_hat = as.numeric(theta_hat$psi),
    loglik_hat = loglik_hat,
    loglik_true_theta = loglik_true_theta,
    loglik_vMF_candidate = loglik_vmf_candidate,
    loglik_best_positive_branch = loglik_best_positive_branch,
    loglik_best_negative_branch = loglik_best_negative_branch,
    loglik_profile_refine_final = if (identical(method, "profile_refine")) loglik_profile_refine_final else NA_real_,
    loglik_refined_positive_best = if (identical(method, "profile_refine")) loglik_refined_positive_best else NA_real_,
    loglik_refined_negative_best = if (identical(method, "profile_refine")) loglik_refined_negative_best else NA_real_,
    kappa_refined_positive_best = if (identical(method, "profile_refine")) kappa_refined_positive_best else NA_real_,
    psi_refined_positive_best = if (identical(method, "profile_refine")) psi_refined_positive_best else NA_real_,
    kappa_refined_negative_best = if (identical(method, "profile_refine")) kappa_refined_negative_best else NA_real_,
    psi_refined_negative_best = if (identical(method, "profile_refine")) psi_refined_negative_best else NA_real_,
    delta_loglik_refined_neg_minus_pos = if (identical(method, "profile_refine")) delta_loglik_refined_neg_minus_pos else NA_real_,
    mu_angle_error_deg = angle_error_deg(theta_hat$mu, theta_true$mu),
    abs_kappa_error = abs(as.numeric(theta_hat$kappa) - as.numeric(theta_true$kappa)),
    abs_psi_error = abs(as.numeric(theta_hat$psi) - as.numeric(theta_true$psi)),
    psi_sign_correct = sign(as.numeric(theta_hat$psi)) == sign(as.numeric(theta_true$psi)),
    elapsed_seconds = elapsed_seconds,
    convergence_status = convergence_status,
    winner_among_candidates = winner,
    profile_refine_selected_branch = if (identical(method, "profile_refine")) profile_refine_selected_branch else NA_character_,
    profile_refine_selected_candidate = if (identical(method, "profile_refine")) profile_refine_selected_candidate else NA_character_,
    stringsAsFactors = FALSE
  )
}

summarize_method_group <- function(df_group) {
  metric_summary <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q90 = as.numeric(stats::quantile(x, probs = 0.9, na.rm = TRUE)),
      max = max(x, na.rm = TRUE)
    )
  }

  angle <- metric_summary(df_group$mu_angle_error_deg)
  kappa <- metric_summary(df_group$abs_kappa_error)
  psi <- metric_summary(df_group$abs_psi_error)
  elapsed <- metric_summary(df_group$elapsed_seconds)

  data.frame(
    scenario = df_group$scenario[[1L]],
    n = df_group$n[[1L]],
    method = df_group$method[[1L]],
    mu_angle_error_deg_mean = angle[["mean"]],
    mu_angle_error_deg_median = angle[["median"]],
    mu_angle_error_deg_q90 = angle[["q90"]],
    mu_angle_error_deg_max = angle[["max"]],
    abs_kappa_error_mean = kappa[["mean"]],
    abs_kappa_error_median = kappa[["median"]],
    abs_kappa_error_q90 = kappa[["q90"]],
    abs_kappa_error_max = kappa[["max"]],
    abs_psi_error_mean = psi[["mean"]],
    abs_psi_error_median = psi[["median"]],
    abs_psi_error_q90 = psi[["q90"]],
    abs_psi_error_max = psi[["max"]],
    psi_sign_correct_mean = mean(as.numeric(df_group$psi_sign_correct), na.rm = TRUE),
    elapsed_seconds_mean = elapsed[["mean"]],
    elapsed_seconds_median = elapsed[["median"]],
    failure_rate = mean(!grepl("^ok", df_group$convergence_status), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

run_jp_observed_mle_method_comparison <- function(n_values = c(50L, 100L, 200L),
                                                  M_outer = 100L,
                                                  seed = 123L,
                                                  n_cores = 10L,
                                                  output_dir = file.path("output", "bootstrap_calibration", "jp_observed_mle_method_comparison"),
                                                  methods = c("current_fast", "current_precise", "profile_refine"),
                                                  fast_control = list(
                                                    jp_mle_max_abs_kappa_psi = 6,
                                                    jp_mle_maxit = 20L,
                                                    jp_mle_reltol = 1e-4,
                                                    jp_profile_n_u = 1025L,
                                                    jp_profile_n_delta = 257L
                                                  ),
                                                  precise_control = list(
                                                    jp_mle_max_abs_kappa_psi = 6,
                                                    jp_mle_maxit = 200L,
                                                    jp_mle_reltol = 1e-10,
                                                    jp_profile_n_u = 1025L,
                                                    jp_profile_n_delta = 257L
                                                  ),
                                                  robust_control = list(
                                                    max_abs_kappa_psi = 6,
                                                    kappa_grid = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0),
                                                    psi_abs_grid = c(0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0),
                                                    method_mu = "Nelder-Mead",
                                                    maxit_mu = 150L,
                                                    reltol_mu = 1e-8,
                                                    top_k_per_branch = 3L,
                                                    psi_min = 1e-3,
                                                    maxit_refine = 500L,
                                                    reltol_refine = 1e-10
                                                  )) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  scenarios <- default_jp_composite_calibration_scenarios()
  all_rows <- list()
  idx <- 1L

  for (scenario in scenarios) {
    scenario_id <- scenario$id
    theta_true <- scenario$sample_params

    for (n in as.integer(n_values)) {
      cat(sprintf("\n[Observed MLE comparison] Scenario=%s n=%d M_outer=%d\n", scenario_id, n, as.integer(M_outer)))

      worker <- function(m) {
        set.seed(as.integer(seed) + 100000L * match(n, n_values) + 1000L * m)
        x <- simulate_h0_sample(scenario = scenario, n = n, replicate_id = m)
        x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
        prob_weights <- rep.int(1 / nrow(x), nrow(x))

        robust_control_local <- robust_control
        robust_control_local$mu_starts <- make_mu_starts(x)

        profile_context <- compute_profile_context(
          x = x,
          prob_weights = prob_weights,
          theta_true = theta_true,
          robust_control = robust_control_local
        )

        do.call(rbind, lapply(methods, function(method) {
          run_observed_method(
            method = method,
            x = x,
            theta_true = theta_true,
            scenario_id = scenario_id,
            n = n,
            outer_id = m,
            prob_weights = prob_weights,
            fast_control = fast_control,
            precise_control = precise_control,
            robust_control = robust_control_local,
            profile_context = profile_context
          )
        }))
      }

      rows_list <- if (as.integer(n_cores) > 1L) {
        parallel::mclapply(seq_len(as.integer(M_outer)), worker, mc.cores = as.integer(n_cores), mc.preschedule = TRUE)
      } else {
        lapply(seq_len(as.integer(M_outer)), worker)
      }

      block_rows <- do.call(rbind, rows_list)
      all_rows[[idx]] <- block_rows
      idx <- idx + 1L
    }
  }

  raw_df <- do.call(rbind, all_rows)
  rownames(raw_df) <- NULL

  groups <- split(raw_df, interaction(raw_df$scenario, raw_df$n, raw_df$method, drop = TRUE))
  summary_df <- do.call(rbind, lapply(groups, summarize_method_group))
  rownames(summary_df) <- NULL

  branch_winners <- subset(raw_df, grepl("kappa_1_psi_0p5", scenario))
  winner_summary <- aggregate(
    outer_id ~ scenario + n + method + winner_among_candidates,
    data = transform(branch_winners, outer_id = 1L),
    FUN = sum
  )
  names(winner_summary)[names(winner_summary) == "outer_id"] <- "count"

  raw_csv <- file.path(output_dir, "jp_observed_mle_method_comparison_raw.csv")
  summary_csv <- file.path(output_dir, "jp_observed_mle_method_comparison_summary.csv")
  winner_csv <- file.path(output_dir, "jp_observed_mle_branch_winners_truepsi_pos.csv")

  utils::write.csv(raw_df, raw_csv, row.names = FALSE)
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)
  utils::write.csv(winner_summary, winner_csv, row.names = FALSE)

  list(
    raw = raw_df,
    summary = summary_df,
    winner_summary = winner_summary,
    raw_csv = raw_csv,
    summary_csv = summary_csv,
    winner_csv = winner_csv
  )
}

if (sys.nframe() == 0L) {
  result <- run_jp_observed_mle_method_comparison(
    n_values = c(50L, 100L, 200L),
    M_outer = 100L,
    seed = 123L,
    n_cores = 10L
  )

  cat("raw_csv=", result$raw_csv, "\n", sep = "")
  cat("summary_csv=", result$summary_csv, "\n", sep = "")
  cat("winner_csv=", result$winner_csv, "\n", sep = "")
}

run_jp_jones_pewsey_population_sanity_check <- function(n = 50000L,
                                                        seed = 123L,
                                                        n_cores = 1L,
                                                        output_dir = file.path("output", "bootstrap_calibration", "jp_jones_pewsey_population_sanity"),
                                                        fixed_mu_kappa_grid = seq(0, 2.5, by = 0.1),
                                                        fixed_mu_psi_grid = seq(-1.5, 1.5, by = 0.1),
                                                        profile_surface_kappa_grid = seq(0, 2.0, by = 0.2),
                                                        profile_surface_psi_grid = seq(-1.2, 1.2, by = 0.2),
                                                        robust_control_override = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  scenario <- default_jp_composite_calibration_scenarios()[[1L]]
  theta_true <- scenario$sample_params

  set.seed(seed)
  x <- simulate_h0_sample(scenario = scenario, n = as.integer(n), replicate_id = 1L)
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- rep.int(1 / nrow(x), nrow(x))

  robust_control <- list(
    max_abs_kappa_psi = 6,
    kappa_grid = seq(0, 3, length.out = 61L),
    psi_abs_grid = c(seq(0.05, 1.0, length.out = 40L), seq(1.1, 2.0, length.out = 10L)),
    method_mu = "Nelder-Mead",
    maxit_mu = 200L,
    reltol_mu = 1e-9,
    top_k_per_branch = 5L,
    psi_min = 1e-3,
    maxit_refine = 800L,
    reltol_refine = 1e-11
  )
  if (!is.null(robust_control_override)) {
    robust_control <- utils::modifyList(robust_control, robust_control_override)
  }
  robust_control$mu_starts <- make_mu_starts(x, n_extra_when_concentrated = 8L)

  profile_context <- compute_profile_context(
    x = x,
    prob_weights = prob_weights,
    theta_true = theta_true,
    robust_control = robust_control
  )

  profile_row <- run_observed_method(
    method = "profile_refine",
    x = x,
    theta_true = theta_true,
    scenario_id = paste0(scenario$id, "_population_sanity"),
    n = nrow(x),
    outer_id = 1L,
    prob_weights = prob_weights,
    fast_control = list(jp_mle_max_abs_kappa_psi = 6, jp_mle_maxit = 20L, jp_mle_reltol = 1e-4),
    precise_control = list(jp_mle_max_abs_kappa_psi = 6, jp_mle_maxit = 200L, jp_mle_reltol = 1e-10),
    robust_control = robust_control,
    profile_context = profile_context
  )

  current_precise_row <- run_observed_method(
    method = "current_precise",
    x = x,
    theta_true = theta_true,
    scenario_id = paste0(scenario$id, "_population_sanity"),
    n = nrow(x),
    outer_id = 1L,
    prob_weights = prob_weights,
    fast_control = list(jp_mle_max_abs_kappa_psi = 6, jp_mle_maxit = 20L, jp_mle_reltol = 1e-4),
    precise_control = list(jp_mle_max_abs_kappa_psi = 6, jp_mle_maxit = 500L, jp_mle_reltol = 1e-11),
    robust_control = robust_control,
    profile_context = profile_context
  )

  best_pos_row <- if (nrow(profile_context$branch_candidates_pos) > 0L) {
    profile_context$branch_candidates_pos[which.max(profile_context$branch_candidates_pos$loglik), , drop = FALSE]
  } else {
    data.frame()
  }
  best_neg_row <- if (nrow(profile_context$branch_candidates_neg) > 0L) {
    profile_context$branch_candidates_neg[which.max(profile_context$branch_candidates_neg$loglik), , drop = FALSE]
  } else {
    data.frame()
  }

  build_population_row <- function(method_name,
                                   theta_hat,
                                   elapsed_seconds,
                                   loglik_hat,
                                   profile_row_ref = profile_row,
                                   current_precise_ref = current_precise_row) {
    data.frame(
      scenario = paste0(scenario$id, "_population_sanity"),
      n = nrow(x),
      outer_id = 1L,
      method = method_name,
      mu_hat = if (!is.null(theta_hat)) format_mu_hat(theta_hat$mu) else NA_character_,
      mu1 = if (!is.null(theta_hat)) as.numeric(theta_hat$mu[[1L]]) else NA_real_,
      mu2 = if (!is.null(theta_hat)) as.numeric(theta_hat$mu[[2L]]) else NA_real_,
      mu3 = if (!is.null(theta_hat)) as.numeric(theta_hat$mu[[3L]]) else NA_real_,
      kappa_hat = if (!is.null(theta_hat)) as.numeric(theta_hat$kappa) else NA_real_,
      psi_hat = if (!is.null(theta_hat)) as.numeric(theta_hat$psi) else NA_real_,
      loglik_hat = as.numeric(loglik_hat),
      loglik_true_theta = profile_context$loglik_true_theta,
      loglik_vMF_candidate = profile_context$loglik_vmf_candidate,
      loglik_best_positive_branch = profile_context$loglik_best_positive_branch,
      loglik_best_negative_branch = profile_context$loglik_best_negative_branch,
      loglik_profile_refine_final = profile_row_ref$loglik_profile_refine_final[[1L]],
      loglik_refined_positive_best = profile_row_ref$loglik_refined_positive_best[[1L]],
      loglik_refined_negative_best = profile_row_ref$loglik_refined_negative_best[[1L]],
      kappa_refined_positive_best = profile_row_ref$kappa_refined_positive_best[[1L]],
      psi_refined_positive_best = profile_row_ref$psi_refined_positive_best[[1L]],
      kappa_refined_negative_best = profile_row_ref$kappa_refined_negative_best[[1L]],
      psi_refined_negative_best = profile_row_ref$psi_refined_negative_best[[1L]],
      delta_loglik_refined_neg_minus_pos = profile_row_ref$delta_loglik_refined_neg_minus_pos[[1L]],
      mu_angle_error_deg = if (!is.null(theta_hat)) angle_error_deg(theta_hat$mu, theta_true$mu) else NA_real_,
      abs_kappa_error = if (!is.null(theta_hat)) abs(as.numeric(theta_hat$kappa) - as.numeric(theta_true$kappa)) else NA_real_,
      abs_psi_error = if (!is.null(theta_hat)) abs(as.numeric(theta_hat$psi) - as.numeric(theta_true$psi)) else NA_real_,
      psi_sign_correct = if (!is.null(theta_hat)) sign(as.numeric(theta_hat$psi)) == sign(as.numeric(theta_true$psi)) else NA,
      elapsed_seconds = as.numeric(elapsed_seconds),
      convergence_status = if (identical(method_name, "current_precise")) current_precise_ref$convergence_status[[1L]] else "ok",
      winner_among_candidates = if (identical(method_name, "current_precise")) current_precise_ref$winner_among_candidates[[1L]] else NA_character_,
      profile_refine_selected_branch = if (identical(method_name, "profile_refine_final")) profile_row_ref$profile_refine_selected_branch[[1L]] else NA_character_,
      profile_refine_selected_candidate = if (identical(method_name, "profile_refine_final")) profile_row_ref$profile_refine_selected_candidate[[1L]] else NA_character_,
      stringsAsFactors = FALSE
    )
  }

  t_true <- proc.time()[["elapsed"]]
  theta_true_loglik <- safe_loglik(theta_true, x, prob_weights)
  elapsed_true <- proc.time()[["elapsed"]] - t_true

  t_vmf <- proc.time()[["elapsed"]]
  vmf_theta <- profile_context$vmf_candidate
  loglik_vmf <- safe_loglik(vmf_theta, x, prob_weights)
  elapsed_vmf <- proc.time()[["elapsed"]] - t_vmf

  t_pos <- proc.time()[["elapsed"]]
  pos_theta <- if (nrow(best_pos_row) > 0L) {
    list(mu = as.numeric(best_pos_row[1L, c("mu1", "mu2", "mu3")]),
         kappa = as.numeric(best_pos_row$kappa[[1L]]),
         psi = as.numeric(best_pos_row$psi[[1L]]))
  } else {
    NULL
  }
  loglik_pos <- if (!is.null(pos_theta)) safe_loglik(pos_theta, x, prob_weights) else NA_real_
  elapsed_pos <- proc.time()[["elapsed"]] - t_pos

  t_neg <- proc.time()[["elapsed"]]
  neg_theta <- if (nrow(best_neg_row) > 0L) {
    list(mu = as.numeric(best_neg_row[1L, c("mu1", "mu2", "mu3")]),
         kappa = as.numeric(best_neg_row$kappa[[1L]]),
         psi = as.numeric(best_neg_row$psi[[1L]]))
  } else {
    NULL
  }
  loglik_neg <- if (!is.null(neg_theta)) safe_loglik(neg_theta, x, prob_weights) else NA_real_
  elapsed_neg <- proc.time()[["elapsed"]] - t_neg

  t_rpos <- proc.time()[["elapsed"]]
  rpos_theta <- if (is.finite(profile_row$loglik_refined_positive_best[[1L]])) {
    list(mu = pos_theta$mu,
         kappa = profile_row$kappa_refined_positive_best[[1L]],
         psi = profile_row$psi_refined_positive_best[[1L]])
  } else {
    NULL
  }
  loglik_rpos <- if (!is.null(rpos_theta)) safe_loglik(rpos_theta, x, prob_weights) else NA_real_
  elapsed_rpos <- proc.time()[["elapsed"]] - t_rpos

  t_rneg <- proc.time()[["elapsed"]]
  rneg_theta <- if (is.finite(profile_row$loglik_refined_negative_best[[1L]])) {
    list(mu = neg_theta$mu,
         kappa = profile_row$kappa_refined_negative_best[[1L]],
         psi = profile_row$psi_refined_negative_best[[1L]])
  } else {
    NULL
  }
  loglik_rneg <- if (!is.null(rneg_theta)) safe_loglik(rneg_theta, x, prob_weights) else NA_real_
  elapsed_rneg <- proc.time()[["elapsed"]] - t_rneg

  profile_final_theta <- list(
    mu = c(profile_row$mu1[[1L]], profile_row$mu2[[1L]], profile_row$mu3[[1L]]),
    kappa = profile_row$kappa_hat[[1L]],
    psi = profile_row$psi_hat[[1L]]
  )

  current_precise_theta <- list(
    mu = c(current_precise_row$mu1[[1L]], current_precise_row$mu2[[1L]], current_precise_row$mu3[[1L]]),
    kappa = current_precise_row$kappa_hat[[1L]],
    psi = current_precise_row$psi_hat[[1L]]
  )

  rows <- do.call(rbind, list(
    build_population_row("true_theta", theta_true, elapsed_true, theta_true_loglik),
    build_population_row("vMF_candidate", vmf_theta, elapsed_vmf, loglik_vmf),
    build_population_row("best_positive_branch", pos_theta, elapsed_pos, loglik_pos),
    build_population_row("best_negative_branch", neg_theta, elapsed_neg, loglik_neg),
    build_population_row("best_refined_positive", rpos_theta, elapsed_rpos, loglik_rpos),
    build_population_row("best_refined_negative", rneg_theta, elapsed_rneg, loglik_rneg),
    build_population_row("profile_refine_final", profile_final_theta, profile_row$elapsed_seconds[[1L]], profile_row$loglik_hat[[1L]]),
    build_population_row("current_precise", current_precise_theta, current_precise_row$elapsed_seconds[[1L]], current_precise_row$loglik_hat[[1L]])
  ))

  mu_true <- normalize_mu(theta_true$mu)
  fixed_mu_rows <- list()
  profile_surface_rows <- list()
  idx_fixed <- 1L
  idx_surface <- 1L

  for (kappa in as.numeric(fixed_mu_kappa_grid)) {
    for (psi in as.numeric(fixed_mu_psi_grid)) {
      if (!is.finite(kappa) || !is.finite(psi) || kappa < 0) {
        next
      }
      if (psi != 0 && abs(kappa * psi) > robust_control$max_abs_kappa_psi) {
        next
      }

      theta_fixed <- list(mu = mu_true, kappa = kappa, psi = psi)
      loglik_fixed <- safe_loglik(theta_fixed, x, prob_weights)

      if (abs(psi) <= robust_control$psi_min) {
        mu_profile_fit <- optimize_mu_given_kappa_psi(
          x = x,
          prob_weights = prob_weights,
          kappa = kappa,
          psi = 0,
          mu_starts = robust_control$mu_starts,
          maxit_mu = robust_control$maxit_mu,
          reltol_mu = robust_control$reltol_mu,
          method_mu = robust_control$method_mu
        )
      } else {
        mu_profile_fit <- optimize_mu_given_kappa_psi(
          x = x,
          prob_weights = prob_weights,
          kappa = kappa,
          psi = psi,
          mu_starts = robust_control$mu_starts,
          maxit_mu = robust_control$maxit_mu,
          reltol_mu = robust_control$reltol_mu,
          method_mu = robust_control$method_mu
        )
      }

      loglik_profile_mu <- mu_profile_fit$loglik
      fixed_mu_rows[[idx_fixed]] <- data.frame(
        kappa = kappa,
        psi = psi,
        loglik_fixed_mu_true = loglik_fixed,
        loglik_profile_mu_optimized = loglik_profile_mu,
        stringsAsFactors = FALSE
      )
      idx_fixed <- idx_fixed + 1L
    }
  }

  for (kappa in as.numeric(profile_surface_kappa_grid)) {
    for (psi in as.numeric(profile_surface_psi_grid)) {
      if (!is.finite(kappa) || !is.finite(psi) || kappa < 0) {
        next
      }
      if (psi != 0 && abs(kappa * psi) > robust_control$max_abs_kappa_psi) {
        next
      }

      theta_fixed <- list(mu = mu_true, kappa = kappa, psi = psi)
      loglik_fixed <- safe_loglik(theta_fixed, x, prob_weights)
      mu_profile_fit <- optimize_mu_given_kappa_psi(
        x = x,
        prob_weights = prob_weights,
        kappa = kappa,
        psi = psi,
        mu_starts = robust_control$mu_starts,
        maxit_mu = robust_control$maxit_mu,
        reltol_mu = robust_control$reltol_mu,
        method_mu = robust_control$method_mu
      )

      profile_surface_rows[[idx_surface]] <- data.frame(
        kappa = kappa,
        psi = psi,
        loglik_fixed_mu_true = loglik_fixed,
        loglik_profile_mu_optimized = mu_profile_fit$loglik,
        stringsAsFactors = FALSE
      )
      idx_surface <- idx_surface + 1L
    }
  }

  fixed_mu_grid <- do.call(rbind, fixed_mu_rows)
  fixed_mu_grid <- fixed_mu_grid[order(-fixed_mu_grid$loglik_fixed_mu_true), , drop = FALSE]
  profile_surface_df <- do.call(rbind, profile_surface_rows)
  profile_surface_df <- profile_surface_df[order(-profile_surface_df$loglik_profile_mu_optimized), , drop = FALSE]

  row_csv <- file.path(output_dir, "jp_population_sanity_rows.csv")
  grid_csv <- file.path(output_dir, "jp_population_sanity_fixed_mu_grid.csv")
  profile_surface_csv <- file.path(output_dir, "jp_population_profile_surface.csv")
  utils::write.csv(rows, row_csv, row.names = FALSE)
  utils::write.csv(fixed_mu_grid, grid_csv, row.names = FALSE)
  utils::write.csv(profile_surface_df, profile_surface_csv, row.names = FALSE)

  print(rows[, c(
    "method", "kappa_hat", "psi_hat", "loglik_hat", "loglik_true_theta",
    "loglik_vMF_candidate", "loglik_best_positive_branch", "loglik_best_negative_branch",
    "loglik_refined_positive_best", "loglik_refined_negative_best",
    "delta_loglik_refined_neg_minus_pos", "profile_refine_selected_branch",
    "psi_sign_correct", "elapsed_seconds"
  )])
  cat("Top fixed-mu grid rows:\n")
  print(head(fixed_mu_grid, 20L))

  cap_violations <- sum(abs(rows$kappa_hat * rows$psi_hat) > robust_control$max_abs_kappa_psi, na.rm = TRUE)
  cat(sprintf("Cap violations in population sanity rows (|kappa*psi| > %.2f): %d\n", robust_control$max_abs_kappa_psi, cap_violations))

  list(
    rows = rows,
    fixed_mu_grid = fixed_mu_grid,
    profile_surface = profile_surface_df,
    row_csv = row_csv,
    grid_csv = grid_csv,
    profile_surface_csv = profile_surface_csv
  )
}