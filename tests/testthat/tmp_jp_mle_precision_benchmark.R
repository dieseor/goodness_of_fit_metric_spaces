source("bootstrap/calibration_study.R")

vmf_mle_s2_closed_form <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  resultant <- colSums(x)
  resultant_norm <- sqrt(sum(resultant^2))

  if (!is.finite(resultant_norm) || resultant_norm <= .Machine$double.eps) {
    return(list(mu = c(0, 0, 1), kappa = 0, rbar = 0))
  }

  mu_hat <- resultant / resultant_norm
  rbar <- resultant_norm / n

  if (rbar < 1e-10) {
    kappa_hat <- 0
  } else if (rbar >= 1 - 1e-10) {
    kappa_hat <- 1 / max(1 - rbar, .Machine$double.eps)
  } else {
    A3 <- function(kappa) {
      if (kappa < 1e-6) {
        kappa / 3
      } else {
        1 / tanh(kappa) - 1 / kappa
      }
    }
    upper <- max(10, 2 / (1 - rbar))
    while (A3(upper) < rbar && upper < 1e6) {
      upper <- 2 * upper
    }
    kappa_hat <- uniroot(
      function(kappa) A3(kappa) - rbar,
      lower = 0,
      upper = upper,
      tol = 1e-10
    )$root
  }

  list(mu = mu_hat, kappa = kappa_hat, rbar = rbar)
}

time_per_call <- function(fun, inner_reps = 1L) {
  inner_reps <- as.integer(inner_reps)
  elapsed <- system.time({
    for (j in seq_len(inner_reps)) {
      fun()
    }
  })[["elapsed"]]
  elapsed / inner_reps
}

summarize_times <- function(times) {
  c(
    median = median(times),
    mean = mean(times),
    min = min(times),
    max = max(times)
  )
}

benchmark_jp_mle_precision <- function(n_values = c(50, 100, 200),
                                       reps = 10L,
                                       seed = 123L,
                                       n_cores = 10L,
                                       configs = list(
                                         precise_80_1e_6 = list(maxit = 80L, reltol = 1e-6),
                                         mid_30_1e_5 = list(maxit = 30L, reltol = 1e-5),
                                         fast_20_1e_4 = list(maxit = 20L, reltol = 1e-4),
                                         very_fast_10_1e_3 = list(maxit = 10L, reltol = 1e-3)
                                       )) {
  base_control <- list(
    jp_mle_maxit = 250L,
    jp_mle_reltol = 1e-9,
    jp_mle_psi_min = 1e-3,
    jp_mle_psi_abs_starts = c(0.25, 0.5, 1, 2)
  )

  fit_levels <- c("vMF_closed_form", "JP_fast_warm_marginal", "JP_fast_warm_total")

  task_grid <- expand.grid(
    n = n_values,
    config = names(configs),
    fit = fit_levels,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  task_grid$seed <- as.integer(seed +
    100000L * task_grid$n +
    10000L * match(task_grid$config, names(configs)) +
    1000L * match(task_grid$fit, fit_levels) +
    task_grid$rep)

  run_one <- function(task) {
    n <- task$n
    config_name <- task$config
    fit_name <- task$fit
    set.seed(task$seed)

    x_jp <- r_sph_jp(n, mu = c(0, 0, 1), kappa = 1.5, psi = 0.5)
    w_uniform <- rep.int(1, n)
    w_random <- runif(n, min = 0.75, max = 1.25)

    cfg <- configs[[config_name]]
    fast_control <- modifyList(base_control, list(
      jp_mle_maxit = as.integer(cfg$maxit),
      jp_mle_reltol = as.numeric(cfg$reltol),
      jp_mle_psi_abs_starts = c(0.5),
      jp_mle_sign_branches = c(1L),
      jp_mle_warm_start_only = TRUE,
      jp_mle_bootstrap_refit = TRUE
    ))

    inner_reps <- switch(
      fit_name,
      vMF_closed_form = 1000L,
      JP_fast_warm_marginal = 100L,
      JP_fast_warm_total = 1L,
      1L
    )

    elapsed <- switch(
      fit_name,
      vMF_closed_form = time_per_call(
        function() vmf_mle_s2_closed_form(x_jp),
        inner_reps = inner_reps
      ),
      JP_fast_warm_marginal = {
        theta_start <- jp_mle_s2_weighted(
          x_jp,
          weights = w_uniform,
          control = modifyList(base_control, list(
            jp_mle_maxit = as.integer(cfg$maxit),
            jp_mle_reltol = as.numeric(cfg$reltol)
          ))
        )
        time_per_call(
          function() {
            jp_mle_s2_weighted(
              x_jp,
              weights = w_random,
              control = modifyList(fast_control, list(jp_mle_start_theta = theta_start))
            )
          },
          inner_reps = inner_reps
        )
      },
      JP_fast_warm_total = time_per_call(
        function() {
          theta_start <- jp_mle_s2_weighted(
            x_jp,
            weights = w_uniform,
            control = modifyList(base_control, list(
              jp_mle_maxit = as.integer(cfg$maxit),
              jp_mle_reltol = as.numeric(cfg$reltol)
            ))
          )
          jp_mle_s2_weighted(
            x_jp,
            weights = w_random,
            control = modifyList(fast_control, list(jp_mle_start_theta = theta_start))
          )
        },
        inner_reps = inner_reps
      ),
      stop("Unknown fit type: ", fit_name)
    )

    data.frame(
      model = ifelse(grepl("^JP", fit_name), "JP", "vMF"),
      fit = fit_name,
      config = config_name,
      n = n,
      rep = task$rep,
      inner_reps = inner_reps,
      maxit = ifelse(fit_name == "vMF_closed_form", NA_integer_, as.integer(cfg$maxit)),
      reltol = ifelse(fit_name == "vMF_closed_form", NA_real_, as.numeric(cfg$reltol)),
      elapsed_sec = elapsed,
      stringsAsFactors = FALSE
    )
  }

  task_list <- split(task_grid, seq_len(nrow(task_grid)))

  if (n_cores <= 1L) {
    raw <- lapply(task_list, run_one)
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    parallel::clusterEvalQ(cl, {
      source("bootstrap/calibration_study.R")
      NULL
    })
    parallel::clusterExport(
      cl,
      c("configs", "base_control", "vmf_mle_s2_closed_form", "time_per_call"),
      envir = environment()
    )
    raw <- parallel::parLapply(cl, task_list, run_one)
  }

  raw_df <- do.call(rbind, raw)
  rownames(raw_df) <- NULL

  groups <- split(raw_df, list(raw_df$model, raw_df$fit, raw_df$config, raw_df$n), drop = TRUE)
  summary_df <- do.call(rbind, lapply(groups, function(df) {
    times <- summarize_times(df$elapsed_sec)
    baseline <- median(raw_df$elapsed_sec[
      raw_df$fit == "vMF_closed_form" & raw_df$n == unique(df$n)
    ])
    data.frame(
      model = unique(df$model),
      fit = unique(df$fit),
      config = unique(df$config),
      n = unique(df$n),
      reps = nrow(df),
      inner_reps = unique(df$inner_reps),
      maxit = unique(df$maxit),
      reltol = unique(df$reltol),
      median_sec = unname(times[["median"]]),
      mean_sec = unname(times[["mean"]]),
      min_sec = unname(times[["min"]]),
      max_sec = unname(times[["max"]]),
      ratio_vs_vmf_median = unname(times[["median"]]) / baseline,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_df) <- NULL

  summary_df <- summary_df[order(summary_df$n, summary_df$fit, summary_df$config), , drop = FALSE]

  list(summary = summary_df, raw = raw_df)
}

bench_jp_precision <- benchmark_jp_mle_precision(
  n_values = c(50, 100, 200),
  reps = 10L,
  seed = 123L,
  n_cores = 10L
)

bench_jp_precision$summary
