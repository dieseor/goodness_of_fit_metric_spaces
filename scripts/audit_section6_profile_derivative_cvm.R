#!/usr/bin/env Rscript

# Independent numerical audit of three ingredients of the Section 6 Normal
# implementation: the distance profile, its score derivative, and the CvM
# functional.  This is deliberately a diagnostic script: it does not alter
# production code or any simulation result.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    bits <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[bits[[1L]]]] <- if (length(bits) == 1L) "TRUE" else paste(bits[-1L], collapse = "=")
  }
  out
}

csv_integer <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(as.integer(default))
  as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dimensions <- csv_integer(args$dimensions, c(2L, 5L, 10L))
n <- as.integer(args$n %||% 100L)
n_fixed_samples <- as.integer(args$n_fixed_samples %||% 2L)
n_reference <- as.integer(args$n_reference %||% 500000L)
n_deriv <- csv_integer(args$n_deriv, c(1000L, 10000L, 50000L))
cores <- as.integer(args$cores %||% 10L)
seed <- as.integer(args$seed %||% 20260806L)
output_dir <- args$output_dir %||% file.path(
  "simulation_results", "section6_new_scenarios",
  "audit_normal_profile_derivative_cvm"
)

if (any(!is.finite(dimensions)) || any(dimensions < 2L) || n < 10L ||
    n_fixed_samples < 1L || n_reference < max(n_deriv) ||
    any(n_deriv < 100L) || cores < 1L) {
  stop("Invalid dimensions or Monte Carlo sizes.")
}

source(file.path("bootstrap", "multiplier_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

e1 <- function(d) c(1, rep.int(0, d - 1L))
section6_sigma <- function(d) {
  Sigma <- diag(d)
  Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.75
  Sigma
}

gaussian_score <- function(x, theta) {
  x <- as.matrix(x)
  d <- ncol(x)
  Sigma_inv <- solve(theta$Sigma)
  centered <- sweep(x, 2L, theta$mu, "-")
  score_mu <- centered %*% t(Sigma_inv)
  score_sigma <- t(vapply(seq_len(nrow(x)), function(i) {
    rr <- centered[i, , drop = FALSE]
    fast_multiplier_sym_score_to_vech(
      0.5 * (Sigma_inv %*% crossprod(rr) %*% Sigma_inv - Sigma_inv)
    )
  }, numeric(d * (d + 1L) / 2L)))
  cbind(score_mu, score_sigma)
}

theta_vector <- function(theta) c(theta$mu, fast_multiplier_vech(theta$Sigma))
theta_from_vector <- function(value, d) {
  list(mu = value[seq_len(d)], Sigma = fast_multiplier_ivech(value[-seq_len(d)], d))
}

profile_value <- function(spec, omega, radius, theta) {
  as.numeric(spec$profile_eval(
    omega, radius, theta,
    control = list(mvnormal_quadform_method = "auto")
  ))
}

numeric_profile_derivative <- function(spec, omega, radius, theta, step = 1e-4) {
  d <- length(theta$mu)
  value <- theta_vector(theta)
  output <- numeric(length(value))
  for (j in seq_along(value)) {
    h <- step * max(1, abs(value[[j]]))
    plus <- value
    minus <- value
    plus[[j]] <- plus[[j]] + h
    minus[[j]] <- minus[[j]] - h
    theta_plus <- theta_from_vector(plus, d)
    theta_minus <- theta_from_vector(minus, d)
    # The selected Section 6 covariance has smallest eigenvalue 0.25, so the
    # symmetric perturbations used here remain positive definite by a margin.
    if (min(eigen(theta_plus$Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0 ||
        min(eigen(theta_minus$Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop("Numerical derivative left the positive-definite covariance cone.")
    }
    output[[j]] <- (profile_value(spec, omega, radius, theta_plus) -
      profile_value(spec, omega, radius, theta_minus)) / (2 * h)
  }
  output
}

relative_norm <- function(x, y) sqrt(sum((x - y)^2)) / max(sqrt(sum(y^2)), 1e-14)

profile_oracle_rows <- function(d, seed_value) {
  set.seed(seed_value)
  mu <- 0.5 * e1(d)
  theta <- normalize_mvnormal_theta(list(mu = mu, Sigma = diag(d)))
  spec <- make_mvnormal_spec("both")
  x <- mvtnorm::rmvnorm(80L, mean = mu, sigma = diag(d))
  centers <- rbind(mu, x[c(1L, 20L, 40L, 80L), , drop = FALSE])
  radii <- unique(as.numeric(quantile(
    as.vector(spec$distance_matrix(x, centers, list())),
    probs = c(.01, .10, .50, .90, .99), names = FALSE, type = 8L
  )))
  rows <- vector("list", nrow(centers) * length(radii))
  k <- 1L
  for (i in seq_len(nrow(centers))) for (radius in radii) {
    got <- profile_value(spec, centers[i, ], radius, theta)
    truth <- stats::pchisq(radius^2, df = d, ncp = sum((centers[i, ] - mu)^2))
    rows[[k]] <- data.frame(
      d = d, center = i, radius = radius, evaluated = got, oracle = truth,
      absolute_error = abs(got - truth), stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  do.call(rbind, rows)
}

derivative_rows_one_sample <- function(d, sample_id, seed_value) {
  set.seed(seed_value)
  spec <- make_mvnormal_spec("both")
  mu <- 0.5 * e1(d)
  Sigma <- section6_sigma(d)
  x <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
  theta <- spec$fit_theta(x, NULL, list(type = "composite"), list())
  influence_obs <- cbind(
    sweep(x, 2L, theta$mu, "-"),
    t(vapply(seq_len(nrow(x)), function(i) {
      rr <- sweep(x[i, , drop = FALSE], 2L, theta$mu, "-")
      fast_multiplier_vech(crossprod(rr) - theta$Sigma)
    }, numeric(d * (d + 1L) / 2L)))
  )

  distance <- spec$distance_matrix(x, x, list())
  selected <- cbind(
    center = c(1L, max(1L, floor(n / 2L)), n),
    column = c(max(1L, floor(.25 * n)), max(1L, floor(.50 * n)), max(1L, floor(.75 * n)))
  )
  # Use three non-degenerate cells of precisely the sample-centred grid used
  # by the KS/CvM statistics.  The high-precision reference is accumulated
  # over independent chunks.  This uses the requested cores without retaining
  # a n_reference by p score matrix in every child process.
  centers <- x[selected[, "center"], , drop = FALSE]
  radii <- vapply(seq_len(nrow(selected)), function(g) {
    sort(distance[selected[g, "center"], ])[selected[g, "column"]]
  }, numeric(1))
  n_chunks <- min(cores, n_reference)
  chunk_sizes <- rep.int(floor(n_reference / n_chunks), n_chunks)
  chunk_sizes[seq_len(n_reference %% n_chunks)] <-
    chunk_sizes[seq_len(n_reference %% n_chunks)] + 1L
  reference_chunk <- function(chunk_id) {
    set.seed(seed_value + 500000L + chunk_id)
    aux <- mvtnorm::rmvnorm(chunk_sizes[[chunk_id]], theta$mu, theta$Sigma)
    score <- gaussian_score(aux, theta)
    indicator <- vapply(seq_len(nrow(centers)), function(g) {
      as.numeric(spec$distance_matrix(aux, centers[g, , drop = FALSE], list()) <= radii[[g]])
    }, numeric(nrow(aux)))
    list(
      sum = crossprod(indicator, score),
      sumsq = crossprod(indicator, score^2)
    )
  }
  reference_chunks <- if (.Platform$OS.type == "unix" && n_chunks > 1L) {
    parallel::mclapply(seq_len(n_chunks), reference_chunk,
      mc.cores = n_chunks, mc.set.seed = FALSE
    )
  } else lapply(seq_len(n_chunks), reference_chunk)
  reference_sum <- Reduce(`+`, lapply(reference_chunks, `[[`, "sum"))
  reference_sumsq <- Reduce(`+`, lapply(reference_chunks, `[[`, "sumsq"))
  derivative_reference <- reference_sum / n_reference
  reference_se <- sqrt(pmax(reference_sumsq / n_reference - derivative_reference^2, 0) / n_reference)
  rows <- list()
  row_id <- 1L
  for (g in seq_len(nrow(selected))) {
    i <- selected[g, "center"]
    radius <- radii[[g]]
    derivative_ref <- derivative_reference[g, ]
    # An independent numerical derivative addresses the identity
    # d F / d theta = E[1{d(X,omega)<=t} psi_theta(X)], rather than merely
    # comparing two MC calculations of the same expression.
    derivative_numeric <- numeric_profile_derivative(spec, x[i, ], radius, theta)
    correction_ref <- drop(influence_obs %*% derivative_ref)
    correction_numeric <- drop(influence_obs %*% derivative_numeric)
    rows[[row_id]] <- data.frame(
      audit = "numeric_vs_reference", d = d, sample_id = sample_id,
      grid_center = i, radius = radius, n_aux = n_reference,
      derivative_relative_error = relative_norm(derivative_numeric, derivative_ref),
      derivative_max_abs_error = max(abs(derivative_numeric - derivative_ref)),
      correction_relative_error = relative_norm(correction_numeric, correction_ref),
      correction_max_abs_error = max(abs(correction_numeric - correction_ref)),
      reference_max_coordinate_mcse = max(reference_se[g, ]),
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
    for (n_aux in n_deriv) {
      aux <- mvtnorm::rmvnorm(n_aux, theta$mu, theta$Sigma)
      derivative_mc <- drop(crossprod(
        as.numeric(spec$distance_matrix(aux, x[i, , drop = FALSE], list()) <= radius),
        gaussian_score(aux, theta)
      )) / n_aux
      correction_mc <- drop(influence_obs %*% derivative_mc)
      rows[[row_id]] <- data.frame(
        audit = "mc_vs_reference", d = d, sample_id = sample_id,
        grid_center = i, radius = radius, n_aux = n_aux,
        derivative_relative_error = relative_norm(derivative_mc, derivative_ref),
        derivative_max_abs_error = max(abs(derivative_mc - derivative_ref)),
        correction_relative_error = relative_norm(correction_mc, correction_ref),
        correction_max_abs_error = max(abs(correction_mc - correction_ref)),
        reference_max_coordinate_mcse = max(reference_se[g, ]),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
  do.call(rbind, rows)
}

cvm_rows <- function(seed_value) {
  set.seed(seed_value)
  out <- list()
  cases <- list(
    normal = list(
      x = mvtnorm::rmvnorm(13L, c(.5, -.2), matrix(c(1, .25, .25, .8), 2L)),
      spec = make_mvnormal_spec("both")
    ),
    lg = list(
      x = rlogistic_gaussian_simplex(13L, c(.2, -.1), matrix(c(1, .25, .25, .8), 2L)),
      spec = make_logistic_gaussian_spec("both")
    ),
    vmf = list(
      x = normalize_vmf_data(rotasym::r_vMF(13L, c(1, 0, 0), 2)),
      spec = make_vmf_spec("geodesic", "xi")
    )
  )
  for (family in names(cases)) {
    item <- cases[[family]]
    theta <- item$spec$fit_theta(item$x, NULL, list(type = "composite"), list())
    prep_dense <- prepare_cvm_observed_data(item$x, item$spec, theta, list(), light = FALSE)
    prep_light <- prepare_cvm_observed_data(item$x, item$spec, theta, list(cvm_block_size = 4L), light = TRUE)
    dist <- item$spec$distance_matrix(item$x, item$x, list())
    empirical <- vapply(seq_len(nrow(dist)), function(i) {
      vapply(dist[i, ], function(t) mean(dist[i, ] <= t), numeric(1))
    }, numeric(nrow(dist)))
    empirical <- t(empirical)
    theoretical <- compute_theoretical_sample_profile_matrix(
      item$spec, item$x, dist, theta, list()
    )
    manual <- nrow(dist) * mean((empirical - theoretical)^2)
    multiplier <- normalize_multiplier_weights(
      generate_multiplier_matrix(
        B = 1L, n = nrow(dist), multiplier_spec = resolve_multiplier_spec(NULL),
        seed = seed_value + match(family, names(cases))
      )[1L, ]
    )
    theta_star <- item$spec$fit_theta(
      item$x, multiplier, list(type = "composite"), list()
    )
    weighted_empirical <- outer(seq_len(nrow(dist)), seq_len(ncol(dist)),
      Vectorize(function(i, j) sum(multiplier[dist[i, ] <= dist[i, j]]) / nrow(dist)
    ))
    theoretical_star <- compute_theoretical_sample_profile_matrix(
      item$spec, item$x, dist, theta_star, list()
    )
    weighted_manual <- nrow(dist) * mean(
      ((weighted_empirical - theoretical_star) - (empirical - theoretical))^2
    )
    engine_weighted <- run_bootstrap_chunk(
      weight_chunk = matrix(multiplier, nrow = 1L), spec = item$spec,
      data = item$x, null = list(type = "composite"), control = list(),
      scale_factor = 1, cvm_prep = prep_dense, want_cvm = TRUE,
      theta_start = theta
    )$cvm[[1L]]
    out[[family]] <- data.frame(
      family = family,
      manual_discrete_cvm = manual,
      dense_cvm = prep_dense$statistic,
      light_cvm = prep_light$statistic,
      dense_absolute_error = abs(manual - prep_dense$statistic),
      light_absolute_error = abs(manual - prep_light$statistic),
      manual_weighted_bootstrap_cvm = weighted_manual,
      engine_weighted_bootstrap_cvm = engine_weighted,
      weighted_bootstrap_absolute_error = abs(weighted_manual - engine_weighted),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

run_parallel <- function(x, fun) {
  if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(x, fun, mc.cores = min(cores, length(x)), mc.set.seed = FALSE)
  } else lapply(x, fun)
}

message("Profile oracle audit (Normal spherical) ...")
profile <- do.call(rbind, run_parallel(seq_along(dimensions), function(i) {
  profile_oracle_rows(dimensions[[i]], seed + 1000L * i)
}))
utils::write.csv(profile, file.path(output_dir, "normal_spherical_profile_oracle.csv"), row.names = FALSE)

message("CvM discrete-functional audit ...")
cvm <- cvm_rows(seed + 9000L)
utils::write.csv(cvm, file.path(output_dir, "cvm_discrete_functional_audit.csv"), row.names = FALSE)

tasks <- expand.grid(d = dimensions, sample_id = seq_len(n_fixed_samples), KEEP.OUT.ATTRS = FALSE)
message(sprintf("Derivative audit: %d fixed fitted samples; reference N=%d ...", nrow(tasks), n_reference))
# Each task itself parallelizes its high-precision reference over `cores`.
# Running these tasks concurrently would multiply the memory requirement.
derivative <- do.call(rbind, lapply(seq_len(nrow(tasks)), function(i) {
  derivative_rows_one_sample(
    d = tasks$d[[i]], sample_id = tasks$sample_id[[i]],
    seed_value = seed + 100000L * tasks$d[[i]] + 1000L * tasks$sample_id[[i]]
  )
}))
utils::write.csv(derivative, file.path(output_dir, "normal_derivative_audit.csv"), row.names = FALSE)

summary <- data.frame(
  component = c("normal_spherical_profile", "cvm_observed_discrete_functional",
                "cvm_weighted_bootstrap_functional"),
  max_absolute_error = c(
    max(profile$absolute_error),
    max(c(cvm$dense_absolute_error, cvm$light_absolute_error)),
    max(cvm$weighted_bootstrap_absolute_error)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
derivative_summary <- aggregate(
  cbind(derivative_max_abs_error, derivative_relative_error,
        correction_max_abs_error, correction_relative_error) ~ audit + n_aux,
  derivative, function(z) median(z)
)
utils::write.csv(derivative_summary, file.path(output_dir, "normal_derivative_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  dimensions = paste(dimensions, collapse = ","), n = n,
  n_fixed_samples = n_fixed_samples, n_reference = n_reference,
  n_deriv = paste(n_deriv, collapse = ","), cores = cores, seed = seed
), file.path(output_dir, "config.csv"), row.names = FALSE)
message("Audit written to: ", normalizePath(output_dir))
