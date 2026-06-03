source("utils.R")

rotmix_logit_boundary_flags <- function(theta, control = list()) {
  weight_eps <- as.numeric(control$logitnormal_mixture2_weight_eps %||% 0.01)
  sd_lower <- as.numeric(control$logitnormal_mixture2_sd_lower %||% 0.05)
  sd_upper <- as.numeric(control$logitnormal_mixture2_sd_upper %||% 5)
  mean_lower <- as.numeric(control$logitnormal_mixture2_mean_lower %||% -8)
  mean_upper <- as.numeric(control$logitnormal_mixture2_mean_upper %||% 8)
  clip_means <- isTRUE(control$logitnormal_mixture2_clip_means %||% FALSE)
  tol <- 1e-10

  data.frame(
    est_weight1 = theta$weight1,
    est_mean1 = theta$mean1,
    est_sd1 = theta$sd1,
    est_mean2 = theta$mean2,
    est_sd2 = theta$sd2,
    at_weight_lower = abs(theta$weight1 - weight_eps) <= tol,
    at_weight_upper = abs(theta$weight1 - (1 - weight_eps)) <= tol,
    at_sd1_lower = abs(theta$sd1 - sd_lower) <= tol,
    at_sd1_upper = abs(theta$sd1 - sd_upper) <= tol,
    at_sd2_lower = abs(theta$sd2 - sd_lower) <= tol,
    at_sd2_upper = abs(theta$sd2 - sd_upper) <= tol,
    at_mean1_lower = clip_means && abs(theta$mean1 - mean_lower) <= tol,
    at_mean1_upper = clip_means && abs(theta$mean1 - mean_upper) <= tol,
    at_mean2_lower = clip_means && abs(theta$mean2 - mean_lower) <= tol,
    at_mean2_upper = clip_means && abs(theta$mean2 - mean_upper) <= tol,
    mean_clip_active = clip_means,
    mean_lower = mean_lower,
    mean_upper = mean_upper,
    sd_lower = sd_lower,
    sd_upper = sd_upper,
    weight_eps = weight_eps
  )
}

run_logit_diagnostic_config <- function(config_name,
                                        control,
                                        n_values = c(50L, 200L),
                                        reps = 8L,
                                        seed = 20260601L) {
  theta_true <- logitnormal_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(-0.25, 0.3, 0.92), arg_name = "mu", min_length = 3L),
    weight1 = 0.45,
    mean1 = -1.0,
    sd1 = 0.45,
    mean2 = 1.15,
    sd2 = 0.55
  ))

  rows <- list()
  idx <- 1L
  set.seed(seed)
  for (n in as.integer(n_values)) {
    for (rep_id in seq_len(reps)) {
      x <- r_sph_logitnormal_mixture2(
        n = n,
        mu = theta_true$mu,
        weight1 = theta_true$weight1,
        mean1 = theta_true$mean1,
        sd1 = theta_true$sd1,
        mean2 = theta_true$mean2,
        sd2 = theta_true$sd2
      )

      elapsed <- system.time(
        fit <- logitnormal_mixture2_mle_s2_weighted(
          x = x,
          control = control
        )
      )[["elapsed"]]

      ll_true <- logitnormal_mixture2_weighted_loglik_s2(
        mu = theta_true$mu,
        weight1 = theta_true$weight1,
        mean1 = theta_true$mean1,
        sd1 = theta_true$sd1,
        mean2 = theta_true$mean2,
        sd2 = theta_true$sd2,
        x = x
      )
      ll_hat <- logitnormal_mixture2_weighted_loglik_s2(
        mu = fit$mu,
        weight1 = fit$weight1,
        mean1 = fit$mean1,
        sd1 = fit$sd1,
        mean2 = fit$mean2,
        sd2 = fit$sd2,
        x = x
      )

      rows[[idx]] <- cbind(
        data.frame(
          config = config_name,
          model = "logit",
          n = n,
          rep = rep_id,
          elapsed = elapsed,
          conv = fit$opt$convergence,
          loglik_true = ll_true,
          loglik_hat = ll_hat,
          loglik_gain = ll_hat - ll_true,
          mu_alignment = sum(theta_true$mu * fit$mu),
          weight1_abs_err = abs(fit$weight1 - theta_true$weight1),
          mean1_abs_err = abs(fit$mean1 - theta_true$mean1),
          mean2_abs_err = abs(fit$mean2 - theta_true$mean2),
          sd1_abs_err = abs(fit$sd1 - theta_true$sd1),
          sd2_abs_err = abs(fit$sd2 - theta_true$sd2),
          stringsAsFactors = FALSE
        ),
        rotmix_logit_boundary_flags(fit, control = control)
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

run_beta_reference_diagnostic <- function(n_values = c(50L, 200L),
                                          reps = 8L,
                                          seed = 20260611L) {
  theta_true <- beta_mixture2_normalize_theta(list(
    mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
    weight1 = 0.4,
    alpha1 = 2.5,
    beta1 = 8,
    alpha2 = 9,
    beta2 = 2.5
  ))

  rows <- list()
  idx <- 1L
  set.seed(seed)
  for (n in as.integer(n_values)) {
    for (rep_id in seq_len(reps)) {
      x <- r_sph_beta_mixture2(
        n = n,
        mu = theta_true$mu,
        weight1 = theta_true$weight1,
        alpha1 = theta_true$alpha1,
        beta1 = theta_true$beta1,
        alpha2 = theta_true$alpha2,
        beta2 = theta_true$beta2
      )

      elapsed <- system.time(
        fit <- beta_mixture2_mle_s2_weighted(
          x = x,
          control = list(beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
        )
      )[["elapsed"]]

      ll_true <- beta_mixture2_weighted_loglik_s2(
        mu = theta_true$mu,
        weight1 = theta_true$weight1,
        alpha1 = theta_true$alpha1,
        beta1 = theta_true$beta1,
        alpha2 = theta_true$alpha2,
        beta2 = theta_true$beta2,
        x = x
      )
      ll_hat <- beta_mixture2_weighted_loglik_s2(
        mu = fit$mu,
        weight1 = fit$weight1,
        alpha1 = fit$alpha1,
        beta1 = fit$beta1,
        alpha2 = fit$alpha2,
        beta2 = fit$beta2,
        x = x
      )

      rows[[idx]] <- data.frame(
        config = "beta_reference",
        model = "beta",
        n = n,
        rep = rep_id,
        elapsed = elapsed,
        conv = fit$opt$convergence,
        loglik_true = ll_true,
        loglik_hat = ll_hat,
        loglik_gain = ll_hat - ll_true,
        mu_alignment = sum(theta_true$mu * fit$mu),
        weight1_abs_err = abs(fit$weight1 - theta_true$weight1),
        p1_abs_err = abs(fit$alpha1 - theta_true$alpha1),
        p2_abs_err = abs(fit$beta1 - theta_true$beta1),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

main <- function() {
  output_dir <- file.path("output", "rotational_logitnormal_mle_diagnostic")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  configs <- list(
    baseline_clipped = list(
      logitnormal_mixture2_clip_means = TRUE,
      logitnormal_mixture2_mean_lower = -8,
      logitnormal_mixture2_mean_upper = 8,
      logitnormal_mixture2_sd_lower = 0.05,
      logitnormal_mixture2_sd_upper = 5,
      logitnormal_mixture2_weight_eps = 0.01,
      logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    ),
    clipped_mean15 = list(
      logitnormal_mixture2_clip_means = TRUE,
      logitnormal_mixture2_mean_lower = -15,
      logitnormal_mixture2_mean_upper = 15,
      logitnormal_mixture2_sd_lower = 0.05,
      logitnormal_mixture2_sd_upper = 5,
      logitnormal_mixture2_weight_eps = 0.01,
      logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    ),
    clipped_sd015_weps005 = list(
      logitnormal_mixture2_clip_means = TRUE,
      logitnormal_mixture2_mean_lower = -15,
      logitnormal_mixture2_mean_upper = 15,
      logitnormal_mixture2_sd_lower = 0.15,
      logitnormal_mixture2_sd_upper = 5,
      logitnormal_mixture2_weight_eps = 0.05,
      logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    ),
    unrestricted_sd015_weps005 = list(
      logitnormal_mixture2_clip_means = FALSE,
      logitnormal_mixture2_mean_lower = -15,
      logitnormal_mixture2_mean_upper = 15,
      logitnormal_mixture2_sd_lower = 0.15,
      logitnormal_mixture2_sd_upper = 5,
      logitnormal_mixture2_weight_eps = 0.05,
      logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
    )
  )

  logit_rows <- do.call(rbind, lapply(names(configs), function(name) {
    message("Running config: ", name)
    run_logit_diagnostic_config(name, configs[[name]])
  }))
  beta_rows <- run_beta_reference_diagnostic()

  utils::write.csv(logit_rows, file.path(output_dir, "logit_diagnostic_raw.csv"), row.names = FALSE)
  utils::write.csv(beta_rows, file.path(output_dir, "beta_reference_raw.csv"), row.names = FALSE)

  logit_summary <- aggregate(
    cbind(elapsed, loglik_gain, mu_alignment, weight1_abs_err, mean1_abs_err, mean2_abs_err, sd1_abs_err, sd2_abs_err,
          at_weight_lower, at_weight_upper, at_sd1_lower, at_sd1_upper, at_sd2_lower, at_sd2_upper,
          at_mean1_lower, at_mean1_upper, at_mean2_lower, at_mean2_upper) ~ config + n,
    data = logit_rows,
    FUN = mean
  )
  beta_summary <- aggregate(
    cbind(elapsed, loglik_gain, mu_alignment, weight1_abs_err, p1_abs_err, p2_abs_err) ~ config + n,
    data = beta_rows,
    FUN = mean
  )

  utils::write.csv(logit_summary, file.path(output_dir, "logit_diagnostic_summary.csv"), row.names = FALSE)
  utils::write.csv(beta_summary, file.path(output_dir, "beta_reference_summary.csv"), row.names = FALSE)
  message("Diagnostic output: ", output_dir)
}

if (sys.nframe() == 0L) {
  main()
}
