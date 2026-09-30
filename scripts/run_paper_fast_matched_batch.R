#!/usr/bin/env Rscript

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OMP_THREAD_LIMIT = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

source("utils.R")
source("utils.R")
source("utils.R")
source("utils.R")
source("utils.R")
ensure_distance_profile_cpp_loaded()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript run_paper_fast_matched_batch.R SCENARIO D REP_START REP_END")
}

scenario <- as.character(args[[1L]])
d <- as.integer(args[[2L]])
rep_start <- as.integer(args[[3L]])
rep_end <- as.integer(args[[4L]])

if (!scenario %in% paste0("S", 1:8) ||
    !d %in% c(2L, 5L) ||
    !is.finite(rep_start) || !is.finite(rep_end) ||
    rep_start < 1L || rep_end > 1000L || rep_start > rep_end) {
  stop("Invalid batch specification.")
}

n <- 100L
B <- 1000L
alpha <- 0.05

root <- file.path(
  "simulation_results",
  "paper_fast_vs_reestimated_n100_beta0"
)
seed_path <- file.path(root, "fast_seed_manifest.csv")
output_root <- Sys.getenv("PAPER_FAST_MATCHED_OUTPUT_DIR", unset = "")
if (!grepl("^simulation_results/paper_fast_matched_c3_[0-9]{8}_[0-9]{6}$", output_root) ||
    !dir.exists(output_root)) {
  stop("PAPER_FAST_MATCHED_OUTPUT_DIR must be an existing new timestamped directory.")
}
output_dir <- file.path(output_root, "fast_batches")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(seed_path)) {
  stop("Missing fast_seed_manifest.csv")
}

seeds <- utils::read.csv(seed_path, stringsAsFactors = FALSE)

jobs <- seeds[
  seeds$scenario == scenario &
    seeds$d == d &
    seeds$n == 100L &
    abs(seeds$beta) < 1e-14 &
    seeds$rep >= rep_start &
    seeds$rep <= rep_end,
  ,
  drop = FALSE
]

jobs <- jobs[order(jobs$rep), , drop = FALSE]

expected_reps <- seq.int(rep_start, rep_end)

if (nrow(jobs) != length(expected_reps) ||
    !identical(as.integer(jobs$rep), expected_reps)) {
  stop("Seed manifest does not contain exactly the requested rep range.")
}

output_path <- file.path(
  output_dir,
  sprintf(
    "%s_d%d_rep%04d_%04d.csv",
    scenario, d, rep_start, rep_end
  )
)

lock_path <- paste0(output_path, ".lock")
if (!dir.create(lock_path, showWarnings = FALSE)) {
  stop("This batch is already locked: ", output_path)
}
on.exit(unlink(lock_path, recursive = TRUE, force = TRUE), add = TRUE)

empty_results <- function() {
  data.frame(
    scenario = character(),
    d = integer(),
    n = integer(),
    beta = numeric(),
    rep = integer(),
    seed_data = integer(),
    seed_bootstrap = integer(),
    status = character(),
    error_message = character(),
    warning_message = character(),
    ks_observed = numeric(),
    cvm_observed = numeric(),
    ks_pvalue = numeric(),
    cvm_pvalue = numeric(),
    ks_reject = logical(),
    cvm_reject = logical(),
    effective_bootstrap_method = character(),
    elapsed_seconds = numeric(),
    slurm_job_id = character(),
    stringsAsFactors = FALSE
  )
}

atomic_write <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) {
    unlink(tmp, force = TRUE)
    stop("Could not atomically write ", path)
  }
}

existing <- if (file.exists(output_path)) {
  utils::read.csv(output_path, stringsAsFactors = FALSE)
} else {
  empty_results()
}

if (!all(names(empty_results()) %in% names(existing))) {
  stop("Existing batch output has incompatible schema.")
}

done <- if (nrow(existing)) {
  as.integer(existing$rep[existing$status == "ok"])
} else {
  integer()
}

jobs <- jobs[!jobs$rep %in% done, , drop = FALSE]

unit_vector <- function(p, j = 1L) {
  z <- numeric(p)
  z[[j]] <- 1
  z
}

generate_exact_null_sample <- function(job) {
  set.seed(as.integer(job$seed_data))

  ## S1 production RNG order:
  ## labels first, then the Gaussian sample.
  if (scenario == "S1") {
    theta0 <- rep(1 / sqrt(d), d)
    u0 <- theta0
    sigma0 <- diag(d) + 2 * tcrossprod(u0)

    labels <- stats::rbinom(n, 1L, 0) == 1L
    x <- mvtnorm::rmvnorm(
      n,
      mean = theta0,
      sigma = sigma0
    )

    if (any(labels)) {
      stop("Impossible non-null S1 label at beta=0.")
    }
    return(x)
  }

  ## S2 production: Gaussian first, mixture uniforms afterwards.
  if (scenario == "S2") {
    x <- matrix(
      stats::rnorm(n * d),
      nrow = n,
      ncol = d
    )
    invisible(stats::runif(n) < 0)
    return(x)
  }

  ## S3 production: LG null sample first, mixture uniforms afterwards.
  if (scenario == "S3") {
    x <- rlogistic_gaussian_simplex(
      n = n,
      mu_ilr = numeric(d),
      Sigma_ilr = diag(d)
    )
    invisible(stats::runif(n) < 0)
    return(x)
  }

  ## S4 production: LG AR(1) null sample first, mixture uniforms afterwards.
  if (scenario == "S4") {
    x <- rlogistic_gaussian_ar1(
      n = n,
      mu_ilr = numeric(d),
      rho = 0.5
    )
    invisible(stats::runif(n) < 0)
    return(x)
  }

  ## S5--S8 production consumes the mixture uniforms before the null sample.
  if (scenario == "S5") {
    invisible(stats::runif(n) < 0)
    return(rotasym::r_vMF(
      n,
      mu = unit_vector(d + 1L),
      kappa = 2
    ))
  }

  if (scenario == "S6") {
    invisible(stats::runif(n) < 0)
    return(rotasym::r_vMF(
      n,
      mu = unit_vector(d + 1L),
      kappa = 2 * d
    ))
  }

  if (scenario == "S7") {
    invisible(stats::runif(n) < 0)
    return(rhvmf_polar(
      n,
      mu = c(sqrt(2), unit_vector(d)),
      kappa = d
    ))
  }

  if (scenario == "S8") {
    invisible(stats::runif(n) < 0)
    return(rhvmf_polar(
      n,
      mu = c(sqrt(2), unit_vector(d)),
      kappa = 2 * sqrt(d)
    ))
  }

  stop("Unsupported scenario.")
}

run_fast <- function(x, bootstrap_seed) {
  derivative_seed <- as.integer(
    (as.numeric(bootstrap_seed) + 10000019) %% 2147483647
  )
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = alpha,
    n_cores = 1L,
    seed = as.integer(bootstrap_seed),
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = FALSE
    ),
    distance_profile_backend = "r"
  )

  control_common <- list(
    derivative_method = if (scenario %in% paste0("S", 1:4))
      "score_mc" else "quadrature",
    derivative_mc_size = 10000L,
    derivative_mc_seed = derivative_seed,
    fast_multiplier_cvm_block_size = 50L,
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fast_multiplier_fuse_ks_cvm = TRUE,
    fast_multiplier_cache_corrections = if (scenario == "S3") "true" else "auto",
    fast_multiplier_stream_chunk_size = 100L,
    progress_bar = FALSE
  )
  fast_args <- list(
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fuse_ks_cvm = TRUE,
    cache_block_corrections = if (scenario == "S3") "true" else "auto"
  )

  if (scenario == "S1") {
    return(do.call(
      multiplier_bootstrap_restricted_spiked_normal,
      c(common, list(control = control_common), fast_args)
    ))
  }

  if (scenario == "S2") {
    return(do.call(
      multiplier_bootstrap_normal_sigma_Id,
      c(common, list(control = control_common), fast_args)
    ))
  }

  if (scenario == "S3") {
    common$null <- list(
      type = "composite",
      fixed = list(Sigma_ilr = diag(d))
    )

    return(do.call(
      multiplier_bootstrap_logistic_gaussian,
      c(
        common,
        list(
          unknown_param = "mu",
          control = control_common
        )
      )
    ))
  }

  if (scenario == "S4") {
    return(do.call(
      multiplier_bootstrap_logistic_gaussian_ar1,
      c(common, list(control = control_common), fast_args)
    ))
  }

  if (scenario == "S5") {
    return(do.call(
      multiplier_bootstrap_vmf_fixed_kappa,
      c(
        common,
        list(
          kappa = 2,
          distance_type = "geodesic",
          control = c(
            control_common,
            list(vmf_profile_method = "integral")
          )
        ),
        fast_args
      )
    ))
  }

  if (scenario == "S6") {
    return(do.call(
      multiplier_bootstrap_vmf,
      c(
        common,
        list(
          unknown_param = "xi",
          distance_type = "geodesic",
          control = c(
            control_common,
            list(vmf_profile_method = "integral")
          )
        ),
        fast_args
      )
    ))
  }

  if (scenario %in% c("S7", "S8")) {
    return(do.call(
      multiplier_bootstrap_hvmf,
      c(
        common,
        list(
          unknown_param = "both",
          control = c(
            control_common,
            list(hvmf_profile_method = "integral")
          )
        ),
        fast_args
      )
    ))
  }

  stop("Unsupported scenario.")
}

run_one <- function(job) {
  started <- proc.time()[["elapsed"]]
  warnings <- character()

  base <- data.frame(
    scenario = scenario,
    d = d,
    n = n,
    beta = 0,
    rep = as.integer(job$rep),
    seed_data = as.integer(job$seed_data),
    seed_bootstrap = as.integer(job$seed_bootstrap),
    status = "ok",
    error_message = NA_character_,
    warning_message = NA_character_,
    ks_observed = NA_real_,
    cvm_observed = NA_real_,
    ks_pvalue = NA_real_,
    cvm_pvalue = NA_real_,
    ks_reject = NA,
    cvm_reject = NA,
    effective_bootstrap_method = NA_character_,
    elapsed_seconds = NA_real_,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    stringsAsFactors = FALSE
  )

  tryCatch({
    x <- generate_exact_null_sample(job)

    fit <- withCallingHandlers(
      run_fast(x, job$seed_bootstrap),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    base$ks_observed <- fit$inference$ks$observed
    base$cvm_observed <- fit$inference$cvm$observed
    base$ks_pvalue <- fit$inference$ks$p_value
    base$cvm_pvalue <- fit$inference$cvm$p_value
    base$ks_reject <- fit$inference$ks$reject
    base$cvm_reject <- fit$inference$cvm$reject
    base$effective_bootstrap_method <-
      fit$diagnostics$effective_bootstrap_method %||% NA_character_

    if (!identical(
      base$effective_bootstrap_method,
      "fast_multiplier"
    )) {
      stop("Requested fast multiplier bootstrap was not effective.")
    }

    base$warning_message <- if (length(warnings)) {
      paste(unique(warnings), collapse = " | ")
    } else {
      NA_character_
    }

    base$elapsed_seconds <- proc.time()[["elapsed"]] - started
    base
  }, error = function(e) {
    base$status <- "error"
    base$error_message <- conditionMessage(e)
    base$warning_message <- if (length(warnings)) {
      paste(unique(warnings), collapse = " | ")
    } else {
      NA_character_
    }
    base$elapsed_seconds <- proc.time()[["elapsed"]] - started
    base
  })
}

cores <- as.integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")
)
cores <- max(1L, cores)

cat(sprintf(
  "%s d=%d reps=%d:%d pending=%d cores=%d B=%d\n",
  scenario, d, rep_start, rep_end, nrow(jobs), cores, B
))

if (nrow(jobs)) {
  for (first in seq.int(1L, nrow(jobs), by = cores)) {
    idx <- seq.int(
      first,
      min(first + cores - 1L, nrow(jobs))
    )

    rows <- split(jobs[idx, , drop = FALSE], seq_along(idx))

    wave <- if (length(rows) == 1L) {
      list(run_one(rows[[1L]]))
    } else {
      parallel::mclapply(
        rows,
        run_one,
        mc.cores = length(rows),
        mc.preschedule = FALSE,
        mc.set.seed = FALSE,
        mc.cleanup = TRUE
      )
    }

    existing <- rbind(existing, do.call(rbind, wave))

    existing <- existing[
      order(existing$rep, existing$status != "ok"),
      ,
      drop = FALSE
    ]

    ## If a previous failed attempt and a successful retry coexist,
    ## retain the successful row.
    existing <- existing[
      !duplicated(existing$rep),
      ,
      drop = FALSE
    ]

    atomic_write(existing, output_path)

    cat(sprintf(
      "completed through rep %d; ok=%d/%d\n",
      max(jobs$rep[idx]),
      sum(existing$status == "ok"),
      length(expected_reps)
    ))
    flush.console()
  }
}

if (any(existing$status != "ok")) {
  stop("Batch completed with one or more failed replications; see CSV.")
}

if (!setequal(
  as.integer(existing$rep[existing$status == "ok"]),
  expected_reps
)) {
  stop("Batch ended without all requested replications.")
}

cat(sprintf(
  "BATCH COMPLETE: %s d=%d reps=%d:%d\n",
  scenario, d, rep_start, rep_end
))
