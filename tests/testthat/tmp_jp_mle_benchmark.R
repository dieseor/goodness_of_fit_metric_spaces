source("bootstrap/calibration_study.R")

jp_test_control <- list(
  jp_mle_maxit = 250L,
  jp_mle_reltol = 1e-9,
  jp_mle_psi_min = 1e-3,
  jp_mle_psi_abs_starts = c(0.25, 0.5, 1, 2)
)

jp_fast_control <- modifyList(jp_test_control, list(
  jp_mle_maxit = 80L,
  jp_mle_reltol = 1e-6,
  jp_mle_psi_abs_starts = c(0.5),
  jp_mle_sign_branches = c(1L)
))

summarize_times <- function(times) {
  c(
    median = median(times),
    mean = mean(times),
    min = min(times),
    max = max(times)
  )
}

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
  if (inner_reps < 1L) {
    stop("inner_reps must be positive")
  }
  elapsed <- system.time({
    for (j in seq_len(inner_reps)) {
      force(fun())
    }
  })[["elapsed"]]
  elapsed / inner_reps
}

make_jp_mle_benchmark_tasks <- function(n_values = c(50, 100, 200),
                                        reps = 10L,
                                        seed = 123L) {
  task_grid <- expand.grid(
    n = n_values,
    fit = c(
      "vMF_closed_form",
      "full_unweighted",
      "full_weighted",
      "fast_weighted_warm_marginal",
      "fast_weighted_warm_total"
    ),
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  fit_levels <- c(
    "vMF_closed_form",
    "full_unweighted",
    "full_weighted",
    "fast_weighted_warm_marginal",
    "fast_weighted_warm_total"
  )
  task_grid$seed <- as.integer(seed + 100000L * task_grid$n +
    1000L * match(task_grid$fit, fit_levels) +
    task_grid$rep)
  task_grid
}

run_one_jp_mle_benchmark_task <- function(task) {
  n <- task$n
  fit_name <- task$fit
  set.seed(task$seed)

  x_jp <- r_sph_jp(n, mu = c(0, 0, 1), kappa = 1.5, psi = 0.5)
  w_uniform <- rep.int(1, n)
  w_random <- runif(n, min = 0.75, max = 1.25)

  inner_reps <- switch(
    fit_name,
    vMF_closed_form = 1000L,
    fast_weighted_warm_marginal = 100L,
    full_unweighted = 1L,
    full_weighted = 1L,
    fast_weighted_warm_total = 1L,
    1L
  )

  elapsed <- switch(
    fit_name,
    vMF_closed_form = {
      time_per_call(
        function() vmf_mle_s2_closed_form(x_jp),
        inner_reps = inner_reps
      )
    },
    full_unweighted = {
      time_per_call(
        function() {
          jp_mle_s2_weighted(
            x_jp,
            weights = w_uniform,
            control = jp_test_control
          )
        },
        inner_reps = inner_reps
      )
    },
    full_weighted = {
      time_per_call(
        function() {
          jp_mle_s2_weighted(
            x_jp,
            weights = w_random,
            control = jp_test_control
          )
        },
        inner_reps = inner_reps
      )
    },
    fast_weighted_warm_marginal = {
      theta_start <- jp_mle_s2_weighted(
        x_jp,
        weights = w_uniform,
        control = jp_fast_control
      )
      time_per_call(
        function() {
          jp_mle_s2_weighted(
            x_jp,
            weights = w_random,
            control = modifyList(
              jp_fast_control,
              list(jp_mle_start_theta = theta_start)
            )
          )
        },
        inner_reps = inner_reps
      )
    },
    fast_weighted_warm_total = {
      time_per_call(
        function() {
          theta_start <- jp_mle_s2_weighted(
            x_jp,
            weights = w_uniform,
            control = jp_fast_control
          )
          jp_mle_s2_weighted(
            x_jp,
            weights = w_random,
            control = modifyList(
              jp_fast_control,
              list(jp_mle_start_theta = theta_start)
            )
          )
        },
        inner_reps = inner_reps
      )
    },
    stop("Unknown fit type: ", fit_name)
  )

  data.frame(
    model = "JP",
    fit = fit_name,
    n = n,
    rep = task$rep,
    seed = task$seed,
    inner_reps = inner_reps,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

benchmark_jp_mle <- function(n_values = c(50, 100, 200),
                             reps = 10L,
                             seed = 123L,
                             n_cores = 1L) {
  tasks <- make_jp_mle_benchmark_tasks(
    n_values = n_values,
    reps = reps,
    seed = seed
  )

  task_list <- split(tasks, seq_len(nrow(tasks)))

  if (n_cores <= 1L) {
    raw_results <- lapply(task_list, run_one_jp_mle_benchmark_task)
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
      c(
        "jp_test_control",
        "jp_fast_control",
        "vmf_mle_s2_closed_form",
        "time_per_call",
        "run_one_jp_mle_benchmark_task"
      ),
      envir = environment()
    )

    raw_results <- parallel::parLapply(
      cl,
      task_list,
      run_one_jp_mle_benchmark_task
    )
  }

  raw_df <- do.call(rbind, raw_results)
  rownames(raw_df) <- NULL

  split_results <- split(raw_df, list(raw_df$n, raw_df$fit), drop = TRUE)
  summary_df <- do.call(rbind, lapply(split_results, function(df) {
    times <- summarize_times(df$elapsed_sec)
    data.frame(
      model = "JP",
      fit = unique(df$fit),
      n = unique(df$n),
      reps = nrow(df),
      inner_reps = unique(df$inner_reps),
      median_sec = unname(times[["median"]]),
      mean_sec = unname(times[["mean"]]),
      min_sec = unname(times[["min"]]),
      max_sec = unname(times[["max"]]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_df) <- NULL
  summary_df <- summary_df[order(summary_df$n, summary_df$fit), , drop = FALSE]

  list(
    summary = summary_df,
    raw = raw_df
  )
}

bench_jp <- benchmark_jp_mle(
  n_values = c(50, 100, 200),
  reps = 10L,
  seed = 123L,
  n_cores = 10L
)

bench_jp$summary

# Optional: save results
# dir.create("output/mle_benchmarks", recursive = TRUE, showWarnings = FALSE)
# write.csv(bench_jp$summary, "output/mle_benchmarks/jp_mle_benchmark_summary.csv", row.names = FALSE)
# write.csv(bench_jp$raw, "output/mle_benchmarks/jp_mle_benchmark_raw.csv", row.names = FALSE)