source(file.path("bootstrap", "calibration_study.R"))

run_mle_validation_rotmix <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_validation <- function(model_name, n, reps = 8L) {
    rows <- vector("list", reps)
    for (b in seq_len(reps)) {
      if (identical(model_name, "beta")) {
        theta <- rotational_beta_mixture2_normalize_theta(list(
          mu = jp_normalize_unit_vector(c(0.2, -0.35, 0.915), arg_name = "mu", min_length = 3L),
          weight1 = 0.4,
          alpha1 = 2.5,
          beta1 = 8,
          alpha2 = 9,
          beta2 = 2.5
        ))
        x <- r_sph_rotational_beta_mixture2(
          n = n,
          mu = theta$mu,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2
        )
        t_obs <- system.time(
          fit <- rotational_beta_mixture2_mle_s2_weighted(
            x = x,
            control = list(rotational_beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
          )
        )[["elapsed"]]
        w <- abs(stats::rnorm(n))
        w <- w / sum(w)
        t_boot <- system.time(
          fit_w <- rotational_beta_mixture2_mle_s2_weighted(
            x = x,
            weights = w,
            control = list(
              theta_start = fit,
              rotational_beta_mixture2_start_theta = fit,
              rotational_beta_mixture2_warm_start_only = TRUE,
              rotational_beta_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
            )
          )
        )[["elapsed"]]
        rows[[b]] <- data.frame(
          model = "beta",
          n = n,
          rep = b,
          obs_elapsed = t_obs,
          boot_refit_elapsed = t_boot,
          mu_alignment = sum(theta$mu * fit$mu),
          weight1_abs_err = abs(fit$weight1 - theta$weight1),
          p1_abs_err = abs(fit$alpha1 - theta$alpha1),
          p2_abs_err = abs(fit$beta1 - theta$beta1),
          conv_obs = fit$opt$convergence,
          conv_boot = fit_w$opt$convergence
        )
      } else {
        theta <- rotational_logitnormal_mixture2_normalize_theta(list(
          mu = jp_normalize_unit_vector(c(-0.25, 0.3, 0.92), arg_name = "mu", min_length = 3L),
          weight1 = 0.45,
          mean1 = -1.0,
          sd1 = 0.45,
          mean2 = 1.15,
          sd2 = 0.55
        ))
        x <- r_sph_rotational_logitnormal_mixture2(
          n = n,
          mu = theta$mu,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2
        )
        t_obs <- system.time(
          fit <- rotational_logitnormal_mixture2_mle_s2_weighted(
            x = x,
            control = list(rotational_logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9))
          )
        )[["elapsed"]]
        w <- abs(stats::rnorm(n))
        w <- w / sum(w)
        t_boot <- system.time(
          fit_w <- rotational_logitnormal_mixture2_mle_s2_weighted(
            x = x,
            weights = w,
            control = list(
              theta_start = fit,
              rotational_logitnormal_mixture2_start_theta = fit,
              rotational_logitnormal_mixture2_warm_start_only = TRUE,
              rotational_logitnormal_mixture2_optim_control = list(maxit = 250L, reltol = 1e-9)
            )
          )
        )[["elapsed"]]
        rows[[b]] <- data.frame(
          model = "logit",
          n = n,
          rep = b,
          obs_elapsed = t_obs,
          boot_refit_elapsed = t_boot,
          mu_alignment = sum(theta$mu * fit$mu),
          weight1_abs_err = abs(fit$weight1 - theta$weight1),
          p1_abs_err = abs(fit$mean1 - theta$mean1),
          p2_abs_err = abs(fit$mean2 - theta$mean2),
          conv_obs = fit$opt$convergence,
          conv_boot = fit_w$opt$convergence
        )
      }
    }
    do.call(rbind, rows)
  }

  set.seed(20260601L)
  raw <- do.call(rbind, c(
    lapply(c(50L, 200L), function(n) run_validation("beta", n, reps = 8L)),
    lapply(c(50L, 200L), function(n) run_validation("logit", n, reps = 8L))
  ))
  utils::write.csv(raw, file.path(output_dir, "mle_validation_raw.csv"), row.names = FALSE)
  invisible(raw)
}

run_simple_calibrations_rotmix <- function(output_dir) {
  run_bootstrap_calibration_study(
    scenarios = list(
      default_rotational_beta_mixture2_simple_calibration_scenarios()[[1L]],
      default_rotational_logitnormal_mixture2_simple_calibration_scenarios()[[1L]]
    ),
    n_values = c(50L, 200L),
    M_outer = 10L,
    B = 100L,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1L,
    seed = 20260601L,
    output_dir = output_dir,
    show_progress = FALSE,
    verbose = TRUE
  )
}

run_beta_composite_short_rotmix <- function(output_dir) {
  run_bootstrap_calibration_study(
    scenarios = list(default_rotational_beta_mixture2_composite_calibration_scenarios()[[1L]]),
    n_values = c(50L, 200L),
    M_outer = 10L,
    B = 100L,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = 1L,
    seed = 20260602L,
    output_dir = output_dir,
    show_progress = FALSE,
    verbose = TRUE
  )
}

run_comets_rotmix_stage <- function(dataset_name, output_root) {
  source(file.path("scripts", "run_comets_rotational_mixtures_short_long.R"))
  run_comets_rotational_mixtures_short_long(
    output_root = output_root,
    datasets = dataset_name,
    models = c("rotational_beta_mixture2", "rotational_logitnormal_mixture2"),
    B = 100L,
    statistics = c("ks", "cvm"),
    n_cores = 12L,
    seed = if (identical(dataset_name, "short")) 20260611L else 20260612L,
    M_value = 30L,
    ks_t_points = 120L,
    distance_type = "geodesic"
  )
}

main <- function() {
  root <- file.path("output", "rotmix_unattended_sequence_20260601")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  message("[1/5] MLE validation")
  run_mle_validation_rotmix(file.path(root, "01_mle_validation"))

  message("[2/5] Simple calibrations")
  run_simple_calibrations_rotmix(file.path(root, "02_simple_calibrations"))

  message("[3/5] Beta composite short calibration")
  run_beta_composite_short_rotmix(file.path(root, "03_beta_composite_short"))

  message("[4/5] Short-period comets")
  run_comets_rotmix_stage("short", file.path(root, "04_comets_short"))

  message("[5/5] Long-period comets")
  run_comets_rotmix_stage("long", file.path(root, "05_comets_long"))

  message("DONE: ", root)
}

if (sys.nframe() == 0L) {
  main()
}
