resolve_fast_validation_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}


`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

# Helper function to rbind data.frames with differing columns
rbind_fill_frames <- function(frames) {
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(frames, names)))
  frames_aligned <- lapply(frames, function(df) {
    missing_names <- setdiff(all_names, names(df))
    for (nm in missing_names) {
      df[[nm]] <- NA
    }
    df[, all_names, drop = FALSE]
  })
  do.call(rbind, frames_aligned)
}

source(resolve_fast_validation_path("bootstrap", "multiplier_bootstrap.R"))
source(resolve_fast_validation_path("scripts", "path_helpers.R"))

make_sample_ks_grid <- function(spec,
                                data,
                                control = list(),
                                n_omega = 6L,
                                n_t = 6L) {
  data_normalized <- spec_normalize_data(spec, data, control)
  n <- spec_n_obs(spec, data_normalized, control)
  omega_idx <- unique(round(seq.int(1L, n, length.out = min(n_omega, n))))

  if (is.matrix(data_normalized) || is.data.frame(data_normalized)) {
    omega_grid <- data_normalized[omega_idx, , drop = FALSE]
  } else if (inherits(data_normalized, "logistic_gaussian_simplex_data")) {
    omega_grid <- data_normalized$simplex[omega_idx, , drop = FALSE]
  } else {
    omega_grid <- data_normalized[omega_idx]
  }

  distance_matrix <- spec$distance_matrix(data_normalized, omega_grid, control)
  distance_values <- sort(unique(as.numeric(distance_matrix)))
  distance_values <- distance_values[is.finite(distance_values) & distance_values > 0]
  if (length(distance_values) == 0L) {
    distance_values <- c(0.25, 0.5, 0.75)
  }
  probs <- unique(seq(0.15, 0.90, length.out = min(n_t, length(distance_values))))
  t_grid <- as.numeric(stats::quantile(distance_values, probs = probs, names = FALSE, type = 8))
  t_grid <- sort(unique(t_grid[t_grid > 0]))
  if (length(t_grid) == 0L) {
    t_grid <- distance_values[seq_len(min(length(distance_values), 3L))]
  }

  list(omega_grid = omega_grid, t_grid = t_grid)
}

validation_extract_theta_diagnostics <- function(case,
                                                 theta_hat,
                                                 control = list()) {
  if (is.null(theta_hat) || !is.list(theta_hat)) {
    return(list())
  }

  out <- list()
  if (!is.null(theta_hat$mu)) {
    mu_hat <- as.numeric(theta_hat$mu)
    for (j in seq_along(mu_hat)) {
      out[[sprintf("mu_hat_%d", j)]] <- mu_hat[[j]]
    }
  }

  if (identical(case$model, "cardioid")) {
    out$k <- as.integer(case$wrapper_args$k %||% theta_hat$k %||% NA_integer_)
    out$rho_hat <- as.numeric(theta_hat$rho %||% NA_real_)
  } else if (identical(case$model, "jp")) {
    kappa_hat <- as.numeric(theta_hat$kappa %||% NA_real_)
    psi_hat <- as.numeric(theta_hat$psi %||% NA_real_)
    out$kappa_hat <- kappa_hat
    out$psi_hat <- psi_hat
    out$alpha_hat <- if (is.finite(kappa_hat) && is.finite(psi_hat)) tanh(kappa_hat * psi_hat) else NA_real_
    out$beta_hat <- if (is.finite(psi_hat) && psi_hat != 0) 1 / psi_hat else NA_real_
    out$jp_vmf_limit_branch <- as.integer(jp_is_near_zero_vmf_s2(
      ambient_dim = length(theta_hat$mu),
      kappa = kappa_hat,
      psi = psi_hat,
      abs_kappa_psi_tol = as.numeric(control$jp_vmf_switch_abs_kappa_psi %||% jp_vmf_near_zero_abs_kappa_psi_default)
    ))
  } else if (identical(case$model, "beta_mixture2")) {
    shape_lower <- as.numeric(control$beta_mixture2_shape_lower %||% 0.05)
    shape_upper <- as.numeric(control$beta_mixture2_shape_upper %||% 1e3)
    weight_eps <- as.numeric(control$beta_mixture2_weight_eps %||% 0.01)
    alpha1_hat <- as.numeric(theta_hat$alpha1 %||% NA_real_)
    beta1_hat <- as.numeric(theta_hat$beta1 %||% NA_real_)
    alpha2_hat <- as.numeric(theta_hat$alpha2 %||% NA_real_)
    beta2_hat <- as.numeric(theta_hat$beta2 %||% NA_real_)
    shape_values <- c(alpha1_hat, beta1_hat, alpha2_hat, beta2_hat)
    out$weight_hat <- as.numeric(theta_hat$weight1 %||% NA_real_)
    out$alpha1_hat <- alpha1_hat
    out$beta1_hat <- beta1_hat
    out$alpha2_hat <- alpha2_hat
    out$beta2_hat <- beta2_hat
    out$min_shape_hat <- min(shape_values, na.rm = TRUE)
    out$any_shape_le_1 <- as.integer(any(shape_values <= 1, na.rm = TRUE))
    out$any_shape_at_lower_bound <- as.integer(any(abs(shape_values - shape_lower) <= 1e-10, na.rm = TRUE))
    out$any_shape_at_upper_bound <- as.integer(any(abs(shape_values - shape_upper) <= 1e-10, na.rm = TRUE))
    out$weight_at_boundary <- as.integer(
      is.finite(out$weight_hat) &&
        (abs(out$weight_hat - weight_eps) <= 1e-10 || abs(out$weight_hat - (1 - weight_eps)) <= 1e-10)
    )
  } else if (identical(case$model, "small_circle_symmetric_mixture2")) {
    out$kappa_hat <- as.numeric(theta_hat$kappa %||% NA_real_)
    out$nu_hat <- as.numeric(theta_hat$nu %||% NA_real_)
  } else if (identical(case$model, "watson")) {
    out$kappa_hat <- as.numeric(theta_hat$kappa %||% NA_real_)
  } else if (identical(case$model, "small_circle_weighted_mixture2")) {
    out$weight_hat <- as.numeric(theta_hat$pi %||% NA_real_)
    out$kappa1_hat <- as.numeric(theta_hat$kappa1 %||% NA_real_)
    out$nu1_hat <- as.numeric(theta_hat$nu1 %||% NA_real_)
    out$kappa2_hat <- as.numeric(theta_hat$kappa2 %||% NA_real_)
    out$nu2_hat <- as.numeric(theta_hat$nu2 %||% NA_real_)
  }

  out
}

make_validation_cases <- function(n_override = NULL) {
  n_value <- if (is.null(n_override)) 100L else as.integer(n_override)
  if (!is.finite(n_value) || n_value <= 0L) {
    stop("`n_override` must be a strictly positive integer when supplied.")
  }

  list(
    list(
      model = "normal",
      scenario = "normal_both",
      n = n_value,
      spec_fn = function() make_normal_spec(unknown_param = "both"),
      wrapper = multiplier_bootstrap_normal,
      wrapper_args = list(null = list(type = "composite"), unknown_param = "both"),
      sample_fn = function(n) stats::rnorm(n, mean = 0.4, sd = 1.1),
      control = list()
    ),
    list(
      model = "vmf",
      scenario = "vmf_kappa2",
      n = n_value,
      spec_fn = function() make_vmf_spec(distance_type = "geodesic", unknown_param = "xi"),
      wrapper = multiplier_bootstrap_vmf,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic", unknown_param = "xi"),
      sample_fn = function(n) rotasym::r_vMF(n, mu = c(1, 0, 0), kappa = 2),
      control = list()
    ),
    list(
      model = "jp",
      scenario = "jp_regular",
      n = n_value,
      spec_fn = function() make_jp_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_jp,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_jp(n, mu = c(0, 0, 1), kappa = 2, psi = 0.7),
      control = list()
    ),
    list(
      model = "hvmf",
      scenario = "hvmf_kappa3",
      n = n_value,
      spec_fn = function() make_hvmf_spec(unknown_param = "both"),
      wrapper = multiplier_bootstrap_hvmf,
      wrapper_args = list(null = list(type = "composite"), unknown_param = "both"),
      sample_fn = function(n) rhvmf_h2_polar(n, mu = c(cosh(0.35), sinh(0.35), 0), kappa = 3),
      control = list()
    ),
    list(
      model = "logistic_gaussian",
      scenario = "logistic_gaussian_both",
      n = n_value,
      spec_fn = function() make_logistic_gaussian_spec(unknown_param = "both"),
      wrapper = multiplier_bootstrap_logistic_gaussian,
      wrapper_args = list(null = list(type = "composite"), unknown_param = "both"),
      sample_fn = function(n) rlogistic_gaussian_simplex(
        n = n,
        mu_ilr = c(0.3, -0.2),
        Sigma_ilr = matrix(c(0.45, 0.10, 0.10, 0.30), nrow = 2)
      ),
      control = list(logistic_gaussian_quadform_method = "hbe")
    ),
    list(
      model = "spherical_cauchy",
      scenario = "spherical_cauchy_rho04",
      n = n_value,
      spec_fn = function() make_spherical_cauchy_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_spherical_cauchy,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_spherical_cauchy(n, mu = c(0, 0, 1), rho = 0.4),
      control = list()
    ),
    list(
      model = "small_circle",
      scenario = "small_circle_kappa8_nu03",
      n = n_value,
      spec_fn = function() make_small_circle_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_small_circle,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_small_circle(n, mu = c(0, 0, 1), kappa = 8, nu = 0.3),
      control = list()
    ),
    list(
      model = "watson",
      scenario = "watson_kappa8",
      n = n_value,
      spec_fn = function() make_watson_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_watson,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_watson(n, mu = c(0, 0, 1), kappa = 8),
      control = list()
    ),
    list(
      model = "cardioid",
      scenario = "cardioid_k2_rho04",
      n = n_value,
      spec_fn = function() make_cardioid_spec(k = 2L, distance_type = "geodesic", unknown_param = "both"),
      wrapper = multiplier_bootstrap_cardioid,
      wrapper_args = list(null = list(type = "composite"), k = 2L, distance_type = "geodesic", unknown_param = "both"),
      sample_fn = function(n) r_sph_car(n, mu = c(0, 0, 1), rho = 0.4, k = 2),
      control = list()
    ),
    list(
      model = "beta_mixture2",
      scenario = "beta_mixture2_asymmetric",
      n = n_value,
      spec_fn = function() make_beta_mixture2_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_beta_mixture2,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_beta_mixture2(n, mu = c(0, 0, 1), weight1 = 0.4, alpha1 = 2, beta1 = 7, alpha2 = 7, beta2 = 2),
      control = list()
    ),
    list(
      model = "small_circle_symmetric_mixture2",
      scenario = "small_circle_symmetric_mixture2",
      n = n_value,
      spec_fn = function() make_small_circle_symmetric_mixture2_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_small_circle_symmetric_mixture2,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_small_circle_symmetric_mixture2(n, mu = c(0, 0, 1), kappa = 10, nu = 0.4),
      control = list()
    ),
    list(
      model = "small_circle_weighted_mixture2",
      scenario = "small_circle_weighted_mixture2",
      n = n_value,
      spec_fn = function() make_small_circle_weighted_mixture2_spec(distance_type = "geodesic"),
      wrapper = multiplier_bootstrap_small_circle_weighted_mixture2,
      wrapper_args = list(null = list(type = "composite"), distance_type = "geodesic"),
      sample_fn = function(n) r_sph_small_circle_weighted_mixture2(n, mu = c(0, 0, 1), pi = 0.55, kappa1 = 12, nu1 = 0.45, kappa2 = 8, nu2 = 0.35),
      control = list()
    )
  )
}

run_validation_case <- function(case,
                                replicate_id,
                                B,
                                derivative_mc_size,
                                sample_seed,
                                bootstrap_seed,
                                derivative_mc_seed,
                                n_cores = 12L,
                                statistics = c("ks", "cvm")) {
  set.seed(sample_seed)
  data <- case$sample_fn(case$n)
  spec <- case$spec_fn()
  theta_hat_fit <- spec$fit_theta(
    data = data,
    weights = NULL,
    null = case$wrapper_args$null,
    control = case$control
  )
  theta_diag <- validation_extract_theta_diagnostics(
    case = case,
    theta_hat = theta_hat_fit,
    control = case$control
  )
  ks_grid <- make_sample_ks_grid(spec, data, control = case$control)

  base_args <- c(
    list(
      data = data,
      statistics = statistics,
      ks_grid = ks_grid,
      B = B,
      alpha = 0.05,
      seed = bootstrap_seed,
      n_cores = as.integer(n_cores),
      keep = list(observed_process = FALSE, bootstrap_statistics = TRUE)
    ),
    case$wrapper_args
  )

  old_result <- tryCatch(
    do.call(case$wrapper, c(
      base_args,
      list(
        bootstrap_method = "reestimated",
        control = case$control
      )
    )),
    error = identity
  )

  fast_result <- tryCatch(
    do.call(case$wrapper, c(
      base_args,
      list(
        bootstrap_method = "fast_multiplier",
        control = utils::modifyList(
          case$control,
          list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_mc_seed),
            fast_multiplier_cvm_block_size = case$n
          )
        )
      )
    )),
    error = identity
  )

  if (inherits(old_result, "error") || inherits(fast_result, "error")) {
    return(do.call(rbind, lapply(statistics, function(stat_name) {
      data.frame(
        model = case$model,
        scenario = case$scenario,
        n = case$n,
        B = B,
        statistic = stat_name,
        replicate_id = replicate_id,
        observed_statistic = NA_real_,
        derivative_method = if (inherits(fast_result, "error")) NA_character_ else fast_result$diagnostics$derivative_method,
        derivative_mc_size = if (inherits(fast_result, "error")) NA_integer_ else fast_result$diagnostics$derivative_mc_size,
        derivative_mc_seed = if (inherits(fast_result, "error")) NA_integer_ else fast_result$diagnostics$derivative_mc_seed,
        bootstrap_method = if (inherits(fast_result, "error")) "fast_multiplier" else fast_result$diagnostics$bootstrap_method,
        effective_bootstrap_method = if (inherits(fast_result, "error")) NA_character_ else fast_result$diagnostics$effective_bootstrap_method,
        fallback_reason = if (inherits(fast_result, "error")) NA_character_ else fast_result$diagnostics$fallback_reason,
        common_observed_seconds = if (inherits(fast_result, "error")) NA_real_ else fast_result$diagnostics$common_observed_seconds,
        old_p_value = NA_real_,
        fast_p_value = NA_real_,
        abs_p_value_diff = NA_real_,
        old_prep_seconds = if (inherits(old_result, "error")) NA_real_ else old_result$diagnostics$old_prep_seconds,
        old_loop_seconds = if (inherits(old_result, "error")) NA_real_ else old_result$diagnostics$old_loop_seconds,
        old_total_seconds = if (inherits(old_result, "error")) NA_real_ else old_result$diagnostics$old_total_seconds,
        fast_prep_seconds = if (inherits(fast_result, "error")) NA_real_ else fast_result$diagnostics$fast_prep_seconds,
        fast_loop_seconds = if (inherits(fast_result, "error")) NA_real_ else fast_result$diagnostics$fast_loop_seconds,
        fast_total_seconds = if (inherits(fast_result, "error")) NA_real_ else fast_result$diagnostics$fast_total_seconds,
        speedup_factor = NA_real_,
        old_reject_0_01 = NA_integer_,
        old_reject_0_05 = NA_integer_,
        old_reject_0_10 = NA_integer_,
        fast_reject_0_01 = NA_integer_,
        fast_reject_0_05 = NA_integer_,
        fast_reject_0_10 = NA_integer_,
        old_error = if (inherits(old_result, "error")) conditionMessage(old_result) else NA_character_,
        fast_error = if (inherits(fast_result, "error")) conditionMessage(fast_result) else NA_character_,
        stringsAsFactors = FALSE
      ) |> cbind(as.data.frame(theta_diag, stringsAsFactors = FALSE))
    })))
  }

  do.call(rbind, lapply(statistics, function(stat_name) {
    old_inf <- old_result$inference[[stat_name]]
    fast_inf <- fast_result$inference[[stat_name]]
    data.frame(
      model = case$model,
      scenario = case$scenario,
      n = case$n,
      B = B,
      statistic = stat_name,
      replicate_id = replicate_id,
      observed_statistic = old_inf$observed,
      derivative_method = fast_result$diagnostics$derivative_method,
      derivative_mc_size = fast_result$diagnostics$derivative_mc_size,
      derivative_mc_seed = fast_result$diagnostics$derivative_mc_seed,
      bootstrap_method = fast_result$diagnostics$bootstrap_method,
      effective_bootstrap_method = fast_result$diagnostics$effective_bootstrap_method,
      fallback_reason = fast_result$diagnostics$fallback_reason,
      common_observed_seconds = fast_result$diagnostics$common_observed_seconds,
      old_p_value = old_inf$p_value,
      fast_p_value = fast_inf$p_value,
      abs_p_value_diff = abs(old_inf$p_value - fast_inf$p_value),
      old_prep_seconds = old_result$diagnostics$old_prep_seconds,
      old_loop_seconds = old_result$diagnostics$old_loop_seconds,
      old_total_seconds = old_result$diagnostics$old_total_seconds,
      fast_prep_seconds = fast_result$diagnostics$fast_prep_seconds,
      fast_loop_seconds = fast_result$diagnostics$fast_loop_seconds,
      fast_total_seconds = fast_result$diagnostics$fast_total_seconds,
      speedup_factor = old_result$diagnostics$old_total_seconds / fast_result$diagnostics$fast_total_seconds,
      old_reject_0_01 = as.integer(old_inf$p_value <= 0.01),
      old_reject_0_05 = as.integer(old_inf$p_value <= 0.05),
      old_reject_0_10 = as.integer(old_inf$p_value <= 0.10),
      fast_reject_0_01 = as.integer(fast_inf$p_value <= 0.01),
      fast_reject_0_05 = as.integer(fast_inf$p_value <= 0.05),
      fast_reject_0_10 = as.integer(fast_inf$p_value <= 0.10),
      old_error = NA_character_,
      fast_error = NA_character_,
      stringsAsFactors = FALSE
    ) |> cbind(as.data.frame(theta_diag, stringsAsFactors = FALSE))
  }))
}

build_validation_summary <- function(raw_df) {
  safe_mean <- function(x) {
    if (!any(is.finite(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }
  safe_median <- function(x) {
    if (!any(is.finite(x))) {
      return(NA_real_)
    }
    stats::median(x, na.rm = TRUE)
  }
  groups <- split(raw_df, list(raw_df$model, raw_df$scenario, raw_df$statistic), drop = TRUE)
  rows <- lapply(groups, function(df) {
    valid_idx <- is.finite(df$old_p_value) & is.finite(df$fast_p_value)
    corr_value <- if (sum(valid_idx) > 1L) stats::cor(df$old_p_value[valid_idx], df$fast_p_value[valid_idx]) else NA_real_
    data.frame(
      model = df$model[[1L]],
      scenario = df$scenario[[1L]],
      statistic = df$statistic[[1L]],
      M_outer = nrow(df),
      successful_pairs = sum(valid_idx),
      old_p_value_mean = safe_mean(df$old_p_value),
      fast_p_value_mean = safe_mean(df$fast_p_value),
      p_value_correlation = corr_value,
      abs_diff_mean = safe_mean(df$abs_p_value_diff),
      abs_diff_median = safe_median(df$abs_p_value_diff),
      old_rejection_rate_0_05 = safe_mean(df$old_reject_0_05),
      fast_rejection_rate_0_05 = safe_mean(df$fast_reject_0_05),
      old_runtime_mean = safe_mean(df$old_total_seconds),
      fast_runtime_mean = safe_mean(df$fast_total_seconds),
      median_speedup = safe_median(df$speedup_factor),
      fallback_count = sum(!is.na(df$fallback_reason) & nzchar(df$fallback_reason)),
      old_error_count = sum(!is.na(df$old_error) & nzchar(df$old_error)),
      fast_error_count = sum(!is.na(df$fast_error) & nzchar(df$fast_error)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_fast_multiplier_validation_all_models <- function(output_root = canonical_fast_multiplier_validation_dir("validation_all_models"),
                                                      B = 5000L,
                                                      M_outer = 5L,
                                                      derivative_mc_size = 1000L,
                                                      base_seed = 20260613L,
                                                      n_cores = 12L,
                                                      models = NULL,
                                                      n_override = NULL) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  cases <- make_validation_cases(n_override = n_override)
  if (!is.null(models)) {
    keep_models <- as.character(models)
    cases <- Filter(function(case) case$model %in% keep_models, cases)
    if (length(cases) == 0L) {
      stop("`models` did not match any validation case.")
    }
  }
  rows <- list()
  idx <- 1L
  raw_csv <- file.path(output_root, "validation_raw.csv")
  summary_csv <- file.path(output_root, "validation_summary.csv")

  for (case in cases) {
    for (replicate_id in seq_len(M_outer)) {
      message(sprintf("[fast validation] %s replicate %d/%d", case$model, replicate_id, M_outer))
      rows[[idx]] <- run_validation_case(
        case = case,
        replicate_id = replicate_id,
        B = as.integer(B),
        derivative_mc_size = as.integer(derivative_mc_size),
        sample_seed = as.integer(base_seed + 1000L * idx),
        bootstrap_seed = as.integer(base_seed + 2000L * idx),
        derivative_mc_seed = as.integer(base_seed + 3000L * idx),
        n_cores = as.integer(n_cores),
        statistics = c("ks", "cvm")
      )
      partial_raw <- rbind_fill_frames(rows[seq_len(idx)])
      utils::write.csv(partial_raw, raw_csv, row.names = FALSE)
      utils::write.csv(build_validation_summary(partial_raw), summary_csv, row.names = FALSE)
      idx <- idx + 1L
    }
  }

  raw_df <- rbind_fill_frames(rows)
  summary_df <- build_validation_summary(raw_df)
  utils::write.csv(raw_df, raw_csv, row.names = FALSE)
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  list(output_root = output_root, raw_csv = raw_csv, summary_csv = summary_csv, raw = raw_df, summary = summary_df)
}

parse_fast_validation_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

fast_validation_arg <- function(args, key, default = NULL) {
  if (key %in% names(args)) {
    return(args[[key]])
  }
  default
}

if (sys.nframe() == 0L) {
  args <- parse_fast_validation_args(commandArgs(trailingOnly = TRUE))
  result <- run_fast_multiplier_validation_all_models(
    output_root = fast_validation_arg(args, "output_root", canonical_fast_multiplier_validation_dir("validation_all_models")),
    B = as.integer(fast_validation_arg(args, "B", 1000L)),
    M_outer = as.integer(fast_validation_arg(args, "M_outer", 5L)),
    derivative_mc_size = as.integer(fast_validation_arg(args, "derivative_mc_size", 1000L)),
    base_seed = as.integer(fast_validation_arg(args, "base_seed", 20260613L)),
    n_cores = as.integer(fast_validation_arg(args, "n_cores", 12L)),
    models = {
      models_arg <- fast_validation_arg(args, "models", NULL)
      if (is.null(models_arg)) NULL else strsplit(models_arg, ",", fixed = TRUE)[[1L]]
    },
    n_override = {
      n_arg <- fast_validation_arg(args, "n", NULL)
      if (is.null(n_arg)) NULL else as.integer(n_arg)
    }
  )
  cat("Raw CSV:", result$raw_csv, "\n")
  cat("Summary CSV:", result$summary_csv, "\n")
}
