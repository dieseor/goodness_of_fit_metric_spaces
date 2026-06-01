resolve_jp_profile_refine_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {    run_comets_distance_profile_jp_short_observed_method_benchmark(
      output_root = "output/comets_distance_profile_jp/run_20260530_short_ks_bfgs_obs_bootlocal_B1000_12cores",
      observed_methods = "standard_bfgs",
      B_values = 1000L,
      statistic = "ks",
      n_cores = 12L,
      base_seed = 20260540L
    )
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

utils_path_jp_profile_refine <- resolve_jp_profile_refine_path("utils.R")
if (!exists("jp_mle_s2_weighted", mode = "function")) {
  source(utils_path_jp_profile_refine)
}

generate_canonical_lattice <- get("generate_canonical_lattice", mode = "function")
jp_mu_s2_from_raw <- get("jp_mu_s2_from_raw", mode = "function")
jp_weighted_loglik_s2_prepared <- get("jp_weighted_loglik_s2_prepared", mode = "function")
jp_mu_s2_to_raw <- get("jp_mu_s2_to_raw", mode = "function")
jp_weighted_vmf_mle_s2 <- get("jp_weighted_vmf_mle_s2", mode = "function")
jp_neg_weighted_loglik_s2_raw <- get("jp_neg_weighted_loglik_s2_raw", mode = "function")
jp_params_s2_from_raw <- get("jp_params_s2_from_raw", mode = "function")
jp_normalize_unit_matrix <- get("jp_normalize_unit_matrix", mode = "function")
jp_normalize_probability_weights <- get("jp_normalize_probability_weights", mode = "function")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

make_jp_profile_refine_mu_starts <- function(x,
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
  keep <- !duplicated(round(start_matrix, 6))
  start_matrix[keep, , drop = FALSE]
}

optimize_jp_mu_given_kappa_psi <- function(x,
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

profile_jp_branch_candidates <- function(x,
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

      mu_fit <- optimize_jp_mu_given_kappa_psi(
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

      candidates[[idx]] <- data.frame(
        mu1 = mu_fit$mu[[1L]],
        mu2 = mu_fit$mu[[2L]],
        mu3 = mu_fit$mu[[3L]],
        kappa = kappa,
        psi = psi,
        loglik = mu_fit$loglik,
        branch = ifelse(psi_sign > 0, "positive", "negative"),
        source = "profile_kappa_psi_grid",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  if (length(candidates) == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind, candidates)
  out[order(-out$loglik), , drop = FALSE]
}

refine_jp_profile_candidate <- function(candidate,
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

  list(
    theta = list(mu = params$mu, kappa = params$kappa, psi = params$psi),
    status = ifelse(best_fit$convergence == 0L, "ok", "refine_nonzero_conv")
  )
}

fit_jp_theta_profile_refine_observed_details <- function(data,
                                                         weights = NULL,
                                                         control = list()) {
  t0 <- proc.time()[["elapsed"]]
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("`fit_jp_theta_profile_refine_observed_details()` currently supports only S^2 data.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  local_control <- list(
    max_abs_kappa_psi = 6,
    kappa_grid = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0),
    psi_abs_grid = c(0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0),
    method_mu = "Nelder-Mead",
    maxit_mu = 150L,
    reltol_mu = 1e-8,
    top_k_per_branch = 3L,
    psi_min = 1e-3,
    maxit_refine = 500L,
    reltol_refine = 1e-10,
    mu_starts = NULL
  )
  local_control[names(control)] <- control

  if (is.null(local_control$mu_starts)) {
    local_control$mu_starts <- make_jp_profile_refine_mu_starts(x)
  }

  vmf_candidate <- jp_weighted_vmf_mle_s2(x, prob_weights)
  branch_candidates_pos <- profile_jp_branch_candidates(
    x = x,
    prob_weights = prob_weights,
    psi_sign = 1L,
    kappa_grid = local_control$kappa_grid,
    psi_abs_grid = local_control$psi_abs_grid,
    max_abs_kappa_psi = local_control$max_abs_kappa_psi,
    mu_starts = local_control$mu_starts,
    maxit_mu = local_control$maxit_mu,
    reltol_mu = local_control$reltol_mu,
    method_mu = local_control$method_mu
  )
  branch_candidates_neg <- profile_jp_branch_candidates(
    x = x,
    prob_weights = prob_weights,
    psi_sign = -1L,
    kappa_grid = local_control$kappa_grid,
    psi_abs_grid = local_control$psi_abs_grid,
    max_abs_kappa_psi = local_control$max_abs_kappa_psi,
    mu_starts = local_control$mu_starts,
    maxit_mu = local_control$maxit_mu,
    reltol_mu = local_control$reltol_mu,
    method_mu = local_control$method_mu
  )

  candidates_to_refine <- data.frame()
  if (nrow(branch_candidates_pos) > 0L) {
    candidates_to_refine <- rbind(candidates_to_refine, head(branch_candidates_pos[order(-branch_candidates_pos$loglik), , drop = FALSE], local_control$top_k_per_branch))
  }
  if (nrow(branch_candidates_neg) > 0L) {
    candidates_to_refine <- rbind(candidates_to_refine, head(branch_candidates_neg[order(-branch_candidates_neg$loglik), , drop = FALSE], local_control$top_k_per_branch))
  }

  refined <- list()
  refined_status <- character(0)
  if (nrow(candidates_to_refine) > 0L) {
    for (i in seq_len(nrow(candidates_to_refine))) {
      candidate <- candidates_to_refine[i, , drop = FALSE]
      candidate_branch <- as.character(candidate$branch[[1L]])
      refined_fit <- refine_jp_profile_candidate(
        candidate = candidate,
        x = x,
        prob_weights = prob_weights,
        psi_min = local_control$psi_min,
        max_abs_kappa_psi = local_control$max_abs_kappa_psi,
        maxit_refine = local_control$maxit_refine,
        reltol_refine = local_control$reltol_refine
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

  refined_loglik <- vapply(refined, function(entry) {
    jp_weighted_loglik_s2_prepared(
      mu = entry$theta$mu,
      kappa = entry$theta$kappa,
      psi = entry$theta$psi,
      x = x,
      prob_weights = prob_weights
    )
  }, numeric(1))

  if (all(!is.finite(refined_loglik))) {
    stop("Profile-refine observed JP fit failed: no finite refined candidate was found.")
  }

  best_idx <- which.max(refined_loglik)
  best_theta <- refined[[best_idx]]$theta

  list(
    theta = best_theta,
    elapsed_seconds = proc.time()[["elapsed"]] - t0,
    convergence_status = paste(c("ok", unique(refined_status)), collapse = ";"),
    selected_branch = refined[[best_idx]]$selected_branch,
    selected_candidate = refined[[best_idx]]$source,
    loglik = refined_loglik[[best_idx]],
    branch_candidates_pos = branch_candidates_pos,
    branch_candidates_neg = branch_candidates_neg
  )
}

fit_jp_theta_profile_refine_observed <- function(data,
                                                 weights = NULL,
                                                 control = list()) {
  fit_jp_theta_profile_refine_observed_details(
    data = data,
    weights = weights,
    control = control
  )$theta
}