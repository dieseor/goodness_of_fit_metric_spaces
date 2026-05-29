source("bootstrap/calibration_study.R")

# ============================================================
# JP theta diagnostic: observed MLE and bootstrap weighted MLEs
# Same calibration sizes as composite experiment:
#   M_outer = 100, B = 99, n = 50, 100, 200
# Same JP fast controls:
#   maxit = 20, reltol = 1e-4, profile grid irrelevant here
# ============================================================

jp_bootstrap_fast_control <- list(
  jp_mle_maxit = 20L,
  jp_mle_reltol = 1e-4,
  jp_profile_n_u = 1025L,
  jp_profile_n_delta = 257L
)

jp_observed_high_precision_control <- list(
  jp_mle_maxit = 200L,
  jp_mle_reltol = 1e-10,
  jp_profile_n_u = 1025L,
  jp_profile_n_delta = 257L
)

M_outer_diag <- 100L
B_diag <- 99L
n_values_diag <- c(50L, 100L, 200L)
seed_diag <- 123L
n_cores_diag <- 12L
run_theta_diagnostic_outer_replication <- function(m,
                                                   scenario,
                                                   theta_true,
                                                   scenario_id,
                                                   n,
                                                   n_values,
                                                   B,
                                                   seed,
                                                   observed_control,
                                                   bootstrap_control) {
  sample_seed <- seed + 100000L * match(n, n_values) + 1000L * m
  x <- simulate_one_sample(
    scenario = scenario,
    n = n,
    outer_id = m,
    seed_base = sample_seed
  )

  theta_hat <- fit_theta_same_as_bootstrap(
    scenario = scenario,
    data = x,
    weights = NULL,
    theta_start = NULL,
    control_override = observed_control
  )

  observed_df <- extract_theta_row(
    theta = theta_hat,
    theta_true = theta_true,
    scenario_id = scenario_id,
    n = n,
    outer_id = m,
    boot_id = NA_integer_
  )

  set.seed(seed + 900000L + 100000L * match(n, n_values) + 1000L * m)
  raw_multiplier_matrix <- matrix(stats::rexp(B * n, rate = 1), nrow = B, ncol = n)
  normalized_multiplier_matrix <- raw_multiplier_matrix / rowMeans(raw_multiplier_matrix)

  bootstrap_rows <- vector("list", B)

  for (b in seq_len(B)) {
    theta_star <- try(
      fit_theta_same_as_bootstrap(
        scenario = scenario,
        data = x,
        weights = normalized_multiplier_matrix[b, ],
        theta_start = theta_hat,
        control_override = bootstrap_control
      ),
      silent = TRUE
    )

    if (inherits(theta_star, "try-error")) {
      bootstrap_rows[[b]] <- data.frame(
        scenario = scenario_id,
        n = n,
        outer_id = m,
        boot_id = b,
        is_bootstrap = TRUE,
        mu1 = NA_real_, mu2 = NA_real_, mu3 = NA_real_,
        kappa = NA_real_,
        psi = NA_real_,
        true_mu1 = theta_true$mu[[1L]],
        true_mu2 = theta_true$mu[[2L]],
        true_mu3 = theta_true$mu[[3L]],
        true_kappa = theta_true$kappa,
        true_psi = theta_true$psi,
        mu_angle_error_deg = NA_real_,
        kappa_error = NA_real_,
        psi_error = NA_real_,
        abs_kappa_error = NA_real_,
        abs_psi_error = NA_real_,
        psi_sign_correct = NA,
        stringsAsFactors = FALSE
      )
    } else {
      bootstrap_rows[[b]] <- extract_theta_row(
        theta = theta_star,
        theta_true = theta_true,
        scenario_id = scenario_id,
        n = n,
        outer_id = m,
        boot_id = b
      )
    }
  }

  list(
    observed = observed_df,
    bootstrap = do.call(rbind, bootstrap_rows)
  )
}

scenarios_all <- default_jp_composite_calibration_scenarios()

scenario_pos <- scenarios_all[[1]]
scenario_pos$control <- modifyList(scenario_pos$control %||% list(), jp_bootstrap_fast_control)

scenario_neg <- scenarios_all[[2]]
scenario_neg$control <- modifyList(scenario_neg$control %||% list(), jp_bootstrap_fast_control)

# True theta values matching the two default JP composite scenarios.
# If default_jp_composite_calibration_scenarios() uses a different true mu,
# change it here.
true_theta_pos <- list(mu = c(0, 0, 1), kappa = 1, psi = 0.5)
true_theta_neg <- list(mu = c(0, 0, 1), kappa = 2, psi = -0.5)

normalize_mu <- function(mu) {
  mu <- as.numeric(mu)
  mu / sqrt(sum(mu^2))
}

angle_error_deg <- function(mu_hat, mu_true) {
  mu_hat <- normalize_mu(mu_hat)
  mu_true <- normalize_mu(mu_true)
  acos(pmax(pmin(sum(mu_hat * mu_true), 1), -1)) * 180 / pi
}

extract_theta_row <- function(theta, theta_true, scenario_id, n, outer_id, boot_id = NA_integer_) {
  mu_hat <- normalize_mu(theta$mu)
  mu_true <- normalize_mu(theta_true$mu)

  data.frame(
    scenario = scenario_id,
    n = n,
    outer_id = outer_id,
    boot_id = boot_id,
    is_bootstrap = !is.na(boot_id),
    mu1 = mu_hat[[1L]],
    mu2 = mu_hat[[2L]],
    mu3 = mu_hat[[3L]],
    kappa = as.numeric(theta$kappa),
    psi = as.numeric(theta$psi),
    true_mu1 = mu_true[[1L]],
    true_mu2 = mu_true[[2L]],
    true_mu3 = mu_true[[3L]],
    true_kappa = as.numeric(theta_true$kappa),
    true_psi = as.numeric(theta_true$psi),
    mu_angle_error_deg = angle_error_deg(mu_hat, mu_true),
    kappa_error = as.numeric(theta$kappa) - as.numeric(theta_true$kappa),
    psi_error = as.numeric(theta$psi) - as.numeric(theta_true$psi),
    abs_kappa_error = abs(as.numeric(theta$kappa) - as.numeric(theta_true$kappa)),
    abs_psi_error = abs(as.numeric(theta$psi) - as.numeric(theta_true$psi)),
    psi_sign_correct = sign(as.numeric(theta$psi)) == sign(as.numeric(theta_true$psi)),
    stringsAsFactors = FALSE
  )
}

simulate_one_sample <- function(scenario, n, outer_id, seed_base) {
  # Prefer the calibration helper if it exists, because that is closest
  # to what run_bootstrap_calibration_study() uses.
  if (exists("simulate_h0_sample", mode = "function")) {
    out <- try(
      simulate_h0_sample(
        scenario = scenario,
        n = n,
        replicate_id = outer_id,
        seed = seed_base + outer_id
      ),
      silent = TRUE
    )
    if (!inherits(out, "try-error")) {
      return(out)
    }

    out <- try(
      simulate_h0_sample(
        scenario = scenario,
        n = n,
        replicate_id = outer_id
      ),
      silent = TRUE
    )
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }

  set.seed(seed_base + outer_id)

  if (is.function(scenario$simulate)) {
    return(scenario$simulate(n))
  }
  if (is.function(scenario$generator)) {
    return(scenario$generator(n))
  }

  stop("Could not simulate from scenario. Inspect scenario names/functions.")
}

fit_theta_same_as_bootstrap <- function(scenario,
                                        data,
                                        weights = NULL,
                                        theta_start = NULL,
                                        control_override = NULL) {
  spec <- scenario$spec %||% scenario$model_spec

  if (is.null(spec) && is.character(scenario$model) && identical(scenario$model, "jp")) {
    spec <- make_jp_spec(distance_type = scenario$distance_type %||% "geodesic")
  }

  if (is.null(spec) && is.list(scenario$model) && is.function(scenario$model$fit_theta)) {
    spec <- scenario$model
  }

  if (!is.list(spec)) {
    stop(
      "Could not find the model spec in the scenario. Tried `scenario$spec`, ",
      "`scenario$model_spec`, reconstructed JP from `scenario$model == \"jp\"`, ",
      "and list-valued `scenario$model`. Available scenario fields: ",
      paste(names(scenario), collapse = ", ")
    )
  }

  if (!is.function(spec$fit_theta)) {
    stop(
      "The model spec does not contain a callable `fit_theta` function. ",
      "Available spec fields: ", paste(names(spec), collapse = ", ")
    )
  }

  null <- scenario$null
  control <- scenario$control %||% list()
  if (!is.null(control_override)) {
    control <- modifyList(control, control_override)
  }

  if (!is.null(theta_start) && grepl("^jp_", as.character(spec$name))) {
    control$jp_mle_start_theta <- theta_start
    control$jp_mle_warm_start_only <- TRUE
    control$jp_mle_bootstrap_refit <- TRUE
  } else if (!is.null(theta_start) && is.null(control$jp_mle_start_theta)) {
    control$jp_mle_start_theta <- theta_start
  }

  if (exists("spec_normalize_data", mode = "function")) {
    data_normalized <- spec_normalize_data(spec, data, control)
  } else if (is.function(spec$normalize_data)) {
    data_normalized <- spec$normalize_data(data, control)
  } else {
    data_normalized <- data
  }

  spec$fit_theta(
    data = data_normalized,
    weights = weights,
    null = null,
    control = control
  )
}

run_theta_diagnostic_one_scenario <- function(scenario,
                                              theta_true,
                                              scenario_id,
                                              n_values = n_values_diag,
                                              M_outer = M_outer_diag,
                                              B = B_diag,
                                              seed = seed_diag,
                                              n_cores = n_cores_diag,
                                              observed_control = jp_observed_high_precision_control,
                                              bootstrap_control = jp_bootstrap_fast_control) {
  observed_by_n <- list()
  bootstrap_by_n <- list()

  for (n in n_values) {
    cat("\n============================================================\n")
    cat("Theta diagnostic:", scenario_id, "n =", n, "\n")
    cat("M_outer =", M_outer, "B =", B, "n_cores =", n_cores, "\n")
    cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("============================================================\n")

    t0 <- proc.time()[["elapsed"]]

    if (n_cores > 1L) {
      results_m <- parallel::mclapply(
        X = seq_len(M_outer),
        FUN = run_theta_diagnostic_outer_replication,
        scenario = scenario,
        theta_true = theta_true,
        scenario_id = scenario_id,
        n = n,
        n_values = n_values,
        B = B,
        seed = seed,
        observed_control = observed_control,
        bootstrap_control = bootstrap_control,
        mc.cores = n_cores,
        mc.preschedule = TRUE
      )
    } else {
      results_m <- lapply(
        seq_len(M_outer),
        run_theta_diagnostic_outer_replication,
        scenario = scenario,
        theta_true = theta_true,
        scenario_id = scenario_id,
        n = n,
        n_values = n_values,
        B = B,
        seed = seed,
        observed_control = observed_control,
        bootstrap_control = bootstrap_control
      )
    }

    observed_by_n[[as.character(n)]] <- do.call(rbind, lapply(results_m, `[[`, "observed"))
    bootstrap_by_n[[as.character(n)]] <- do.call(rbind, lapply(results_m, `[[`, "bootstrap"))

    t1 <- proc.time()[["elapsed"]]
    cat("End:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Elapsed seconds:", round(t1 - t0, 2), "\n")
  }

  list(
    observed = do.call(rbind, observed_by_n),
    bootstrap = do.call(rbind, bootstrap_by_n)
  )
}

summarise_theta_errors <- function(df) {
  aggregate(
    cbind(
      mu_angle_error_deg,
      abs_kappa_error,
      abs_psi_error,
      psi_sign_correct
    ) ~ scenario + n + is_bootstrap,
    data = transform(df, psi_sign_correct = as.numeric(psi_sign_correct)),
    FUN = function(x) {
      c(
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        q90 = as.numeric(stats::quantile(x, 0.9, na.rm = TRUE)),
        max = max(x, na.rm = TRUE)
      )
    }
  )
}

summarise_bootstrap_relative_to_hat <- function(observed_df, bootstrap_df) {
  key_obs <- paste(observed_df$scenario, observed_df$n, observed_df$outer_id, sep = "__")
  obs_lookup <- split(observed_df, key_obs)

  rows <- lapply(seq_len(nrow(bootstrap_df)), function(i) {
    row <- bootstrap_df[i, ]
    key <- paste(row$scenario, row$n, row$outer_id, sep = "__")
    obs <- obs_lookup[[key]][1, ]

    data.frame(
      scenario = row$scenario,
      n = row$n,
      outer_id = row$outer_id,
      boot_id = row$boot_id,
      boot_mu_angle_from_hat_deg = angle_error_deg(
        c(row$mu1, row$mu2, row$mu3),
        c(obs$mu1, obs$mu2, obs$mu3)
      ),
      boot_kappa_minus_hat = row$kappa - obs$kappa,
      boot_psi_minus_hat = row$psi - obs$psi,
      abs_boot_kappa_minus_hat = abs(row$kappa - obs$kappa),
      abs_boot_psi_minus_hat = abs(row$psi - obs$psi),
      boot_psi_sign_equals_hat = sign(row$psi) == sign(obs$psi),
      stringsAsFactors = FALSE
    )
  })

  rel_df <- do.call(rbind, rows)

  aggregate(
    cbind(
      boot_mu_angle_from_hat_deg,
      abs_boot_kappa_minus_hat,
      abs_boot_psi_minus_hat,
      boot_psi_sign_equals_hat
    ) ~ scenario + n,
    data = transform(rel_df, boot_psi_sign_equals_hat = as.numeric(boot_psi_sign_equals_hat)),
    FUN = function(x) {
      c(
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        q90 = as.numeric(stats::quantile(x, 0.9, na.rm = TRUE)),
        max = max(x, na.rm = TRUE)
      )
    }
  )
}

output_dir <- file.path("output", "bootstrap_calibration", "jp_theta_diagnostics_obsHigh_bootFast_M100_B99")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

diag_pos <- run_theta_diagnostic_one_scenario(
  scenario = scenario_pos,
  theta_true = true_theta_pos,
  scenario_id = "jp_composite_kappa_1_psi_0p5",
  n_values = n_values_diag,
  M_outer = M_outer_diag,
  B = B_diag,
  seed = seed_diag,
  n_cores = n_cores_diag,
  observed_control = jp_observed_high_precision_control,
  bootstrap_control = jp_bootstrap_fast_control
)

diag_neg <- run_theta_diagnostic_one_scenario(
  scenario = scenario_neg,
  theta_true = true_theta_neg,
  scenario_id = "jp_composite_kappa_2_psi_neg0p5",
  n_values = n_values_diag,
  M_outer = M_outer_diag,
  B = B_diag,
  seed = seed_diag,
  n_cores = n_cores_diag,
  observed_control = jp_observed_high_precision_control,
  bootstrap_control = jp_bootstrap_fast_control
)

theta_observed <- rbind(diag_pos$observed, diag_neg$observed)
theta_bootstrap <- rbind(diag_pos$bootstrap, diag_neg$bootstrap)

theta_error_summary_observed <- summarise_theta_errors(theta_observed)
theta_error_summary_bootstrap <- summarise_theta_errors(theta_bootstrap)
theta_bootstrap_relative_to_hat_summary <- rbind(
  summarise_bootstrap_relative_to_hat(diag_pos$observed, diag_pos$bootstrap),
  summarise_bootstrap_relative_to_hat(diag_neg$observed, diag_neg$bootstrap)
)

write.csv(theta_observed,
          file.path(output_dir, "theta_observed_all.csv"),
          row.names = FALSE)
write.csv(theta_bootstrap,
          file.path(output_dir, "theta_bootstrap_all.csv"),
          row.names = FALSE)
write.csv(theta_error_summary_observed,
          file.path(output_dir, "theta_error_summary_observed.csv"),
          row.names = FALSE)
write.csv(theta_error_summary_bootstrap,
          file.path(output_dir, "theta_error_summary_bootstrap.csv"),
          row.names = FALSE)
write.csv(theta_bootstrap_relative_to_hat_summary,
          file.path(output_dir, "theta_bootstrap_relative_to_hat_summary.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat("Theta diagnostics finished.\n")
cat("Saved in:", output_dir, "\n")
cat("Files:\n")
cat(" - theta_observed_all.csv\n")
cat(" - theta_bootstrap_all.csv\n")
cat(" - theta_error_summary_observed.csv\n")
cat(" - theta_error_summary_bootstrap.csv\n")
cat(" - theta_bootstrap_relative_to_hat_summary.csv\n")
cat("============================================================\n\n")

cat("\nObserved theta error summary:\n")
print(theta_error_summary_observed)

cat("\nBootstrap theta error summary relative to true theta:\n")
print(theta_error_summary_bootstrap)

cat("\nBootstrap theta variation relative to theta_hat:\n")
print(theta_bootstrap_relative_to_hat_summary)

cat("\nControl used for theta_hat (observed):\n")
print(jp_observed_high_precision_control)
cat("\nControl used for theta_star (bootstrap refit):\n")
print(jp_bootstrap_fast_control)