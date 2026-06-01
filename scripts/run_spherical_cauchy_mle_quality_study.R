resolve_spherical_cauchy_mle_study_path <- function(...) {
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

utils_script_path_spherical_cauchy_mle <- resolve_spherical_cauchy_mle_study_path("utils.R")
source(utils_script_path_spherical_cauchy_mle)

angle_error_degrees <- function(mu_hat, mu_true) {
  cosine <- sum(mu_hat * mu_true)
  cosine <- min(max(cosine, -1), 1)
  acos(cosine) * 180 / pi
}

finite_difference_gradient <- function(fn, par, eps = 1e-6) {
  grad <- numeric(length(par))
  for (j in seq_along(par)) {
    step <- rep(0, length(par))
    step[[j]] <- eps
    grad[[j]] <- (fn(par + step) - fn(par - step)) / (2 * eps)
  }
  grad
}

run_consistency_checks_spherical_cauchy <- function(seed = 20260601L) {
  set.seed(as.integer(seed))

  x_grad <- jp_normalize_unit_matrix(rbind(
    c(0.0, 0.1, 1.0),
    c(0.1, 0.0, 1.0),
    c(-0.1, 0.05, 1.0),
    c(0.05, -0.1, 1.0)
  ))
  prob_weights <- rep(1 / nrow(x_grad), nrow(x_grad))
  phi <- c(0.1, -0.05, 0.2)

  grad_analytic <- spherical_cauchy_weighted_loglik_phi_grad(phi, x_grad, prob_weights)
  grad_fd <- finite_difference_gradient(
    fn = function(phi_arg) spherical_cauchy_weighted_loglik_phi(phi_arg, x_grad, prob_weights),
    par = phi,
    eps = 1e-6
  )

  x <- r_sph_spherical_cauchy(24, mu = c(0, 0, 1), rho = 0.3)
  weights <- c(2, 1, 3, 2, rep(1, nrow(x) - 4L))
  expanded_x <- x[rep(seq_len(nrow(x)), times = weights), , drop = FALSE]

  fit_unweighted <- spherical_cauchy_mle_s2_weighted(
    x,
    weights = NULL,
    control = list(spherical_cauchy_maxit = 300L)
  )
  fit_equal <- spherical_cauchy_mle_s2_weighted(
    x,
    weights = rep(1, nrow(x)),
    control = list(spherical_cauchy_maxit = 300L, theta_start = fit_unweighted)
  )
  fit_weighted <- spherical_cauchy_mle_s2_weighted(
    x,
    weights = weights,
    control = list(spherical_cauchy_maxit = 300L)
  )
  fit_replicated <- spherical_cauchy_mle_s2_weighted(
    expanded_x,
    weights = NULL,
    control = list(spherical_cauchy_maxit = 300L, theta_start = fit_weighted)
  )

  data.frame(
    check_name = c(
      "gradient_max_abs_difference",
      "equal_weights_vs_unweighted_delta_rho",
      "equal_weights_vs_unweighted_mu_angle_deg",
      "weighted_vs_replicated_delta_rho",
      "weighted_vs_replicated_mu_angle_deg"
    ),
    value = c(
      max(abs(grad_analytic - grad_fd)),
      abs(fit_equal$rho - fit_unweighted$rho),
      angle_error_degrees(fit_equal$mu, fit_unweighted$mu),
      abs(fit_weighted$rho - fit_replicated$rho),
      angle_error_degrees(fit_weighted$mu, fit_replicated$mu)
    ),
    stringsAsFactors = FALSE
  )
}

run_one_mle_fit_spherical_cauchy <- function(rho_true,
                                             n,
                                             mu_true,
                                             replicate_id,
                                             maxit,
                                             reltol,
                                             optim_method,
                                             use_gradient,
                                             seed_base) {
  seed_i <- as.integer(seed_base + 100000L * as.integer(rho_true * 1000) + 1000L * n + replicate_id)
  set.seed(seed_i)

  x <- r_sph_spherical_cauchy(n = n, mu = mu_true, rho = rho_true)

  fit <- try(
    spherical_cauchy_mle_s2_weighted(
      data = x,
      weights = NULL,
      control = list(
        spherical_cauchy_maxit = as.integer(maxit),
        spherical_cauchy_reltol = as.numeric(reltol),
        spherical_cauchy_optim_method = as.character(optim_method),
        spherical_cauchy_use_gradient = isTRUE(use_gradient)
      )
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(data.frame(
      rho_true = rho_true,
      n = as.integer(n),
      replicate_id = as.integer(replicate_id),
      converged = FALSE,
      objective_convergence_code = NA_integer_,
      rho_hat = NA_real_,
      rho_error = NA_real_,
      abs_rho_error = NA_real_,
      mu_angle_deg = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  conv_code <- as.integer(fit$opt$convergence)
  rho_hat <- as.numeric(fit$rho)
  mu_hat <- as.numeric(fit$mu)

  data.frame(
    rho_true = rho_true,
    n = as.integer(n),
    replicate_id = as.integer(replicate_id),
    converged = isTRUE(conv_code == 0L),
    objective_convergence_code = conv_code,
    rho_hat = rho_hat,
    rho_error = rho_hat - rho_true,
    abs_rho_error = abs(rho_hat - rho_true),
    mu_angle_deg = angle_error_degrees(mu_hat, mu_true),
    stringsAsFactors = FALSE
  )
}

run_spherical_cauchy_mle_quality_study <- function(
  output_root = file.path("output", "spherical_cauchy_mle_quality"),
  n_values = c(50L, 100L, 200L, 500L),
  rho_values = c(0.1, 0.3, 0.7),
  n_replicates = 100L,
  mu_true = c(0, 0, 1),
  maxit = 500L,
  reltol = 1e-10,
  optim_method = "BFGS",
  use_gradient = TRUE,
  seed = 20260601L
) {
  n_values <- as.integer(n_values)
  rho_values <- as.numeric(rho_values)
  n_replicates <- as.integer(n_replicates)
  seed <- as.integer(seed)

  if (length(n_values) == 0L || any(!is.finite(n_values)) || any(n_values < 5L)) {
    stop("`n_values` must contain integers >= 5.")
  }
  if (length(rho_values) == 0L || any(!is.finite(rho_values)) || any(rho_values < 0) || any(rho_values >= 1)) {
    stop("`rho_values` must contain finite values in [0, 1).")
  }
  if (length(n_replicates) != 1L || !is.finite(n_replicates) || n_replicates < 10L) {
    stop("`n_replicates` must be an integer >= 10.")
  }

  mu_true <- jp_normalize_unit_vector(mu_true, arg_name = "`mu_true`", min_length = 3L)
  if (length(mu_true) != 3L) {
    stop("`mu_true` must have length 3.")
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(output_root, paste0("run_", timestamp))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  grid <- expand.grid(
    rho_true = rho_values,
    n = n_values,
    replicate_id = seq_len(n_replicates),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "Running spherical Cauchy MLE quality study with %d fits (%d rho values x %d n values x %d replicates).",
    nrow(grid),
    length(rho_values),
    length(n_values),
    n_replicates
  ))

  t0 <- Sys.time()
  raw_rows <- lapply(seq_len(nrow(grid)), function(i) {
    row <- grid[i, , drop = FALSE]
    run_one_mle_fit_spherical_cauchy(
      rho_true = row$rho_true,
      n = row$n,
      mu_true = mu_true,
      replicate_id = row$replicate_id,
      maxit = maxit,
      reltol = reltol,
      optim_method = optim_method,
      use_gradient = use_gradient,
      seed_base = seed
    )
  })
  elapsed_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  raw_df <- do.call(rbind, raw_rows)

  summary_split <- split(raw_df, list(raw_df$rho_true, raw_df$n), drop = TRUE)
  summary_rows <- lapply(summary_split, function(df) {
    converged <- df$converged & is.finite(df$rho_hat) & is.finite(df$mu_angle_deg)
    df_ok <- df[converged, , drop = FALSE]

    if (nrow(df_ok) == 0L) {
      return(data.frame(
        rho_true = unique(df$rho_true),
        n = unique(df$n),
        n_replicates = nrow(df),
        n_converged = 0L,
        convergence_rate = 0,
        rho_bias = NA_real_,
        rho_rmse = NA_real_,
        rho_mae = NA_real_,
        mu_angle_mean_deg = NA_real_,
        mu_angle_median_deg = NA_real_,
        mu_angle_q90_deg = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      rho_true = unique(df$rho_true),
      n = unique(df$n),
      n_replicates = nrow(df),
      n_converged = nrow(df_ok),
      convergence_rate = nrow(df_ok) / nrow(df),
      rho_bias = mean(df_ok$rho_error),
      rho_rmse = sqrt(mean(df_ok$rho_error^2)),
      rho_mae = mean(df_ok$abs_rho_error),
      mu_angle_mean_deg = mean(df_ok$mu_angle_deg),
      mu_angle_median_deg = stats::median(df_ok$mu_angle_deg),
      mu_angle_q90_deg = stats::quantile(df_ok$mu_angle_deg, probs = 0.9, names = FALSE),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(summary_df$rho_true, summary_df$n), , drop = FALSE]

  consistency_df <- run_consistency_checks_spherical_cauchy(seed = seed + 999L)

  utils::write.csv(raw_df, file.path(run_dir, "mle_quality_raw.csv"), row.names = FALSE)
  utils::write.csv(summary_df, file.path(run_dir, "mle_quality_summary.csv"), row.names = FALSE)
  utils::write.csv(consistency_df, file.path(run_dir, "mle_consistency_checks.csv"), row.names = FALSE)

  saveRDS(
    list(
      run_dir = run_dir,
      n_values = n_values,
      rho_values = rho_values,
      n_replicates = n_replicates,
      mu_true = mu_true,
      maxit = as.integer(maxit),
      reltol = as.numeric(reltol),
      optim_method = as.character(optim_method),
      use_gradient = isTRUE(use_gradient),
      seed = seed,
      elapsed_seconds = elapsed_seconds,
      n_total_fits = nrow(grid)
    ),
    file = file.path(run_dir, "run_config.rds")
  )

  writeLines(
    c(
      sprintf("run_dir: %s", run_dir),
      sprintf("elapsed_seconds: %.3f", elapsed_seconds),
      sprintf("n_total_fits: %d", nrow(grid)),
      sprintf("n_replicates: %d", n_replicates),
      sprintf("n_values: %s", paste(n_values, collapse = ",")),
      sprintf("rho_values: %s", paste(rho_values, collapse = ","))
    ),
    con = file.path(run_dir, "study_metadata.txt")
  )

  list(
    run_dir = run_dir,
    raw = raw_df,
    summary = summary_df,
    consistency_checks = consistency_df,
    elapsed_seconds = elapsed_seconds
  )
}

parse_named_args_spherical_cauchy_mle <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
  if (length(args) == 0L) {
    return(list())
  }

  output <- vector("list", length(args))
  names(output) <- rep("", length(args))

  for (i in seq_along(args)) {
    arg <- substring(args[[i]], 3L)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[i]] <- value
    names(output)[[i]] <- key
  }

  output
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_spherical_cauchy_mle(commandArgs(trailingOnly = TRUE))

  output_root <- if (!is.null(args$output_root)) args$output_root else file.path("output", "spherical_cauchy_mle_quality")
  n_values <- if (!is.null(args$n_values)) as.integer(strsplit(args$n_values, ",", fixed = TRUE)[[1L]]) else c(50L, 100L, 200L, 500L)
  rho_values <- if (!is.null(args$rho_values)) as.numeric(strsplit(args$rho_values, ",", fixed = TRUE)[[1L]]) else c(0.1, 0.3, 0.7)
  n_replicates <- if (!is.null(args$n_replicates)) as.integer(args$n_replicates) else 100L
  maxit <- if (!is.null(args$maxit)) as.integer(args$maxit) else 500L
  reltol <- if (!is.null(args$reltol)) as.numeric(args$reltol) else 1e-10
  optim_method <- if (!is.null(args$optim_method)) args$optim_method else "BFGS"
  use_gradient <- if (!is.null(args$use_gradient)) {
    tolower(args$use_gradient) %in% c("1", "true", "t", "yes", "y")
  } else {
    TRUE
  }
  seed <- if (!is.null(args$seed)) as.integer(args$seed) else 20260601L

  result <- run_spherical_cauchy_mle_quality_study(
    output_root = output_root,
    n_values = n_values,
    rho_values = rho_values,
    n_replicates = n_replicates,
    maxit = maxit,
    reltol = reltol,
    optim_method = optim_method,
    use_gradient = use_gradient,
    seed = seed
  )

  message(sprintf("Study completed. Outputs at: %s", result$run_dir))
}
