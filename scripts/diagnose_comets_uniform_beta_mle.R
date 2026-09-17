## Standalone empirical audit of the comet Uniform--Beta likelihood.
## This script only reads production code/data and writes diagnostic artifacts.

source(file.path("utils.R"))
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "uniform_beta_mixture_model_spec.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path("real_data", "comets", "utils_comets_data.R"))

parse_args <- function(x) {
  out <- list()
  for (arg in x[startsWith(x, "--")]) {
    z <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[z[[1L]]]] <- if (length(z) == 1L) TRUE else paste(z[-1L], collapse = "=")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
output_root <- args$output_root %||% file.path(
  "real_data", "comets", "diagnostics", "uniform_beta_mle_audit_pilot"
)
n_multistart <- as.integer(args$n_multistart %||% 40L)
n_local_perturb <- as.integer(args$n_local_perturb %||% 20L)
n_bootstrap <- as.integer(args$n_bootstrap %||% 10L)
n_cores <- as.integer(args$n_cores %||% 3L)
audit_seed <- as.integer(args$seed %||% 20260915L)
if (any(!is.finite(c(n_multistart, n_local_perturb, n_bootstrap, n_cores))) ||
    min(n_multistart, n_local_perturb, n_bootstrap, n_cores) < 1L) {
  stop("All run sizes must be positive integers.")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

production_control <- list(
  uniform_beta_mixture_profile_method = "legendre",
  uniform_beta_mixture_quad_n = 100L,
  uniform_beta_mixture_optim_control = list(maxit = 350L, reltol = 1e-9)
)
bootstrap_control <- utils::modifyList(production_control, list(
  uniform_beta_mixture_bootstrap_n_starts = 1L,
  uniform_beta_mixture_bootstrap_optim_control = list(maxit = 80L, reltol = 1e-6)
))

reference_root <- file.path(
  "real_data", "reruns", "paper_main_realdata_B1000_3cores_20260831_113532",
  "comets", "uniform_beta"
)
reference_dirs <- c(
  short = file.path(reference_root, "01_short_period_uniform_beta_mixture"),
  long = file.path(reference_root, "02_long_period_uniform_beta_mixture")
)

read_theta <- function(path) {
  z <- utils::read.csv(path, stringsAsFactors = FALSE)[1L, ]
  list(
    mu = as.numeric(z[c("mu_1", "mu_2", "mu_3")]),
    weight_uniform = z$weight_uniform,
    alpha = z$alpha,
    beta = z$beta
  )
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(message = conditionMessage(e)), class = "audit_error")
  )
  list(value = value, warnings = unique(warnings))
}

surface_loglik <- function(theta, x, eps = 1e-12) {
  sum(d_sph_uniform_beta_mixture_s2(
    x = x, mu = theta$mu, weight_uniform = theta$weight_uniform,
    alpha = theta$alpha, beta = theta$beta, log = TRUE, eps = eps
  ))
}

theta_vector <- function(theta) {
  c(theta$mu, theta$weight_uniform, theta$alpha, theta$beta)
}

angle <- function(a, b) acos(pmin(pmax(sum(a * b), -1), 1))

nearest_distances <- function(mu, x) {
  dots <- pmin(pmax(as.numeric(x %*% mu), -1), 1)
  c(to_mu = min(acos(dots)), to_antipode = min(acos(-dots)))
}

tangent_basis <- function(mu) {
  axis <- diag(3L)[, which.min(abs(mu))]
  e1 <- axis - sum(axis * mu) * mu
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- c(
    mu[[2L]] * e1[[3L]] - mu[[3L]] * e1[[2L]],
    mu[[3L]] * e1[[1L]] - mu[[1L]] * e1[[3L]],
    mu[[1L]] * e1[[2L]] - mu[[2L]] * e1[[1L]]
  )
  cbind(e1, e2)
}

exp_map <- function(mu0, basis, v) {
  r <- sqrt(sum(v^2))
  if (r < 1e-14) return(mu0)
  cos(r) * mu0 + sin(r) * drop(basis %*% (v / r))
}

chart_theta <- function(q, theta0, basis) {
  list(
    mu = exp_map(theta0$mu, basis, q[1:2]),
    weight_uniform = rotational_bounded_weight(q[[3L]], weight_eps = 0.01),
    alpha = rotational_positive_parameter(q[[4L]], lower = 0.05, upper = 1e3),
    beta = rotational_positive_parameter(q[[5L]], lower = 0.05, upper = 1e3)
  )
}

chart_origin <- function(theta) {
  c(0, 0, stats::qlogis(theta$weight_uniform), log(theta$alpha), log(theta$beta))
}

local_diagnostics <- function(theta, x, n_perturb, seed) {
  basis <- tangent_basis(theta$mu)
  q0 <- chart_origin(theta)
  fn <- function(q) surface_loglik(chart_theta(q, theta, basis), x)
  central_gradient <- function(q, h) vapply(seq_along(q), function(j) {
    step <- rep(0, length(q)); step[[j]] <- h
    (fn(q + step) - fn(q - step)) / (2 * h)
  }, numeric(1L))
  central_hessian <- function(q, h) {
    p <- length(q); out <- matrix(0, p, p); f0 <- fn(q)
    for (j in seq_len(p)) {
      ej <- rep(0, p); ej[[j]] <- h
      out[j, j] <- (fn(q + ej) - 2 * f0 + fn(q - ej)) / h^2
      if (j < p) for (k in (j + 1L):p) {
        ek <- rep(0, p); ek[[k]] <- h
        out[j, k] <- out[k, j] <-
          (fn(q + ej + ek) - fn(q + ej - ek) - fn(q - ej + ek) + fn(q - ej - ek)) / (4 * h^2)
      }
    }
    out
  }
  endpoint_distance <- min(nearest_distances(theta$mu, x))
  hessian_step <- min(1e-5, endpoint_distance / 50)
  gradient_steps <- 10^seq(-3, -8, by = -1)
  gradient_by_step <- do.call(rbind, lapply(gradient_steps, function(h) data.frame(
    step = h, parameter = c("mu_tangent_1", "mu_tangent_2", "logit_w", "log_alpha", "log_beta"),
    gradient = central_gradient(q0, h)
  )))
  gradient <- central_gradient(q0, min(1e-7, endpoint_distance / 100))
  hessian <- central_hessian(q0, hessian_step)
  eig <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values

  set.seed(seed)
  scales <- c(2e-3, 2e-3, 2e-2, 2e-2, 2e-2)
  perturbations <- lapply(seq_len(n_perturb), function(i) {
    start <- q0 + stats::rnorm(5L, sd = scales)
    cw <- capture_warnings(stats::optim(
      start, function(q) -fn(q), method = "BFGS",
      control = production_control$uniform_beta_mixture_optim_control
    ))
    if (inherits(cw$value, "audit_error")) {
      return(data.frame(run = i, convergence = NA_integer_, loglik = NA_real_,
                        angle_to_reported = NA_real_, parameter_distance = NA_real_,
                        warnings = paste(c(cw$warnings, cw$value$message), collapse = " | ")))
    }
    fit <- chart_theta(cw$value$par, theta, basis)
    data.frame(
      run = i, convergence = cw$value$convergence, loglik = -cw$value$value,
      angle_to_reported = angle(fit$mu, theta$mu),
      parameter_distance = sqrt(sum((cw$value$par[3:5] - q0[3:5])^2)),
      warnings = paste(cw$warnings, collapse = " | ")
    )
  })

  list(
    gradient = gradient,
    gradient_by_step = gradient_by_step,
    hessian = hessian,
    hessian_step = hessian_step,
    eigenvalues = eig,
    negative_hessian_condition = if (all(eig < 0)) max(-eig) / min(-eig) else NA_real_,
    perturbations = do.call(rbind, perturbations)
  )
}

single_start_fit <- function(theta0, x, weights = NULL, optim_control = NULL) {
  control <- production_control
  if (!is.null(optim_control)) control$uniform_beta_mixture_optim_control <- optim_control
  par0 <- uniform_beta_mixture_pack_par(theta0, control = control)$par
  prob_weights <- if (is.null(weights)) rep(1 / nrow(x), nrow(x)) else
    jp_normalize_probability_weights(weights, nrow(x))
  objective <- function(par) {
    theta <- uniform_beta_mixture_unpack_par(par, control = control)
    value <- -uniform_beta_mixture_weighted_loglik_s2(
      mu = theta$mu, weight_uniform = theta$weight_uniform,
      alpha = theta$alpha, beta = theta$beta, x = x,
      prob_weights = prob_weights,
      eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
    )
    if (is.finite(value)) value else .Machine$double.xmax / 100
  }
  cw <- capture_warnings(stats::optim(
    par0, objective, method = "BFGS",
    control = control$uniform_beta_mixture_optim_control
  ))
  if (inherits(cw$value, "audit_error")) return(c(cw, list(theta = NULL)))
  c(cw, list(theta = uniform_beta_mixture_unpack_par(cw$value$par, control = control),
             objective = objective))
}

make_multistarts <- function(x, reported, n, seed) {
  set.seed(seed)
  production <- uniform_beta_mixture_start_thetas_s2(x, control = production_control)
  dots <- as.numeric(x %*% reported$mu)
  targeted <- list(
    reported,
    list(mu = -reported$mu, weight_uniform = reported$weight_uniform,
         alpha = reported$beta, beta = reported$alpha),
    list(mu = x[which.max(dots), ], weight_uniform = reported$weight_uniform,
         alpha = reported$alpha, beta = reported$beta),
    list(mu = -x[which.min(dots), ], weight_uniform = reported$weight_uniform,
         alpha = reported$alpha, beta = reported$beta),
    list(mu = x[which.max(dots), ], weight_uniform = 0.1, alpha = 2, beta = 0.2),
    list(mu = -x[which.min(dots), ], weight_uniform = 0.1, alpha = 0.2, beta = 2)
  )
  directions <- lapply(production, `[[`, "mu")
  chosen <- unique(round(seq(1, nrow(x), length.out = min(16L, nrow(x)))))
  directions <- c(directions, lapply(chosen, function(i) x[i, ]),
                  lapply(chosen, function(i) -x[i, ]))
  while (length(targeted) + length(directions) < n) {
    z <- stats::rnorm(3L)
    directions[[length(directions) + 1L]] <- z / sqrt(sum(z^2))
  }
  shape_pairs <- rbind(
    c(0.20, 0.20), c(0.70, 0.70), c(0.70, 2), c(2, 0.70),
    c(1.20, 1.20), c(5, 2), c(20, 0.70), c(40, 1.50)
  )
  weights <- c(0.03, 0.12, 0.30, 0.55, 0.80)
  generated <- lapply(seq_len(max(0L, n - length(targeted))), function(i) list(
    mu = directions[[i]],
    weight_uniform = weights[1L + ((i - 1L) %% length(weights))],
    alpha = shape_pairs[1L + ((i - 1L) %% nrow(shape_pairs)), 1L],
    beta = shape_pairs[1L + ((i - 1L) %% nrow(shape_pairs)), 2L]
  ))
  c(targeted, generated)[seq_len(n)]
}

reported_distances <- function(theta, reported) {
  direct_parameter <- sqrt(sum((
    c(stats::qlogis(theta$weight_uniform), log(theta$alpha), log(theta$beta)) -
      c(stats::qlogis(reported$weight_uniform), log(reported$alpha), log(reported$beta))
  )^2))
  swapped_parameter <- sqrt(sum((
    c(stats::qlogis(theta$weight_uniform), log(theta$alpha), log(theta$beta)) -
      c(stats::qlogis(reported$weight_uniform), log(reported$beta), log(reported$alpha))
  )^2))
  c(
    signed_angle = angle(theta$mu, reported$mu),
    antipodal_angle = angle(theta$mu, -reported$mu),
    equivalent_angle = min(angle(theta$mu, reported$mu), angle(theta$mu, -reported$mu)),
    direct_parameter = direct_parameter,
    swapped_parameter = swapped_parameter,
    equivalent_parameter = min(direct_parameter, swapped_parameter)
  )
}

fit_row <- function(fit, start, run, x, reported) {
  if (is.null(fit$theta)) {
    return(data.frame(
      run = run, convergence = NA_integer_, message = fit$value$message,
      loglik = NA_real_, AIC = NA_real_, BIC = NA_real_,
      mu_1 = NA_real_, mu_2 = NA_real_, mu_3 = NA_real_,
      weight_uniform = NA_real_, alpha = NA_real_, beta = NA_real_,
      distance_to_nearest = NA_real_, antipode_distance_to_nearest = NA_real_,
      angle_to_reported = NA_real_, angle_to_reported_equivalent = NA_real_,
      transformed_parameter_distance = NA_real_,
      same_reported_basin = FALSE, warnings = paste(fit$warnings, collapse = " | ")))
  }
  theta <- fit$theta
  ll <- surface_loglik(theta, x)
  near <- nearest_distances(theta$mu, x)
  rd <- reported_distances(theta, reported)
  data.frame(
    run = run, convergence = fit$value$convergence,
    message = fit$value$message %||% "", loglik = ll,
    AIC = -2 * ll + 10, BIC = -2 * ll + 5 * log(nrow(x)),
    mu_1 = theta$mu[[1L]], mu_2 = theta$mu[[2L]], mu_3 = theta$mu[[3L]],
    weight_uniform = theta$weight_uniform, alpha = theta$alpha, beta = theta$beta,
    distance_to_nearest = near[[1L]], antipode_distance_to_nearest = near[[2L]],
    angle_to_reported = rd[["signed_angle"]],
    angle_to_reported_equivalent = rd[["equivalent_angle"]],
    transformed_parameter_distance = rd[["equivalent_parameter"]],
    same_reported_basin = (rd[["signed_angle"]] < 1e-3 && rd[["direct_parameter"]] < 0.05) ||
      (rd[["antipodal_angle"]] < 1e-3 && rd[["swapped_parameter"]] < 0.05),
    warnings = paste(fit$warnings, collapse = " | ")
  )
}

assign_basins <- function(z) {
  z$basin <- NA_integer_
  good <- which(is.finite(z$loglik))
  centers <- integer()
  for (i in good[order(z$loglik[good], decreasing = TRUE)]) {
    assigned <- FALSE
    for (b in seq_along(centers)) {
      j <- centers[[b]]
      ti <- list(mu = as.numeric(z[i, c("mu_1", "mu_2", "mu_3")]),
                 weight_uniform = z$weight_uniform[[i]], alpha = z$alpha[[i]], beta = z$beta[[i]])
      tj <- list(mu = as.numeric(z[j, c("mu_1", "mu_2", "mu_3")]),
                 weight_uniform = z$weight_uniform[[j]], alpha = z$alpha[[j]], beta = z$beta[[j]])
      rd <- reported_distances(ti, tj)
      equivalent <- (rd[["signed_angle"]] < 0.01 && rd[["direct_parameter"]] < 0.15) ||
        (rd[["antipodal_angle"]] < 0.01 && rd[["swapped_parameter"]] < 0.15)
      if (equivalent) {
        z$basin[[i]] <- b
        assigned <- TRUE
        break
      }
    }
    if (!assigned) {
      centers <- c(centers, i)
      z$basin[[i]] <- length(centers)
    }
  }
  z
}

contribution_tables <- function(theta, x, label) {
  values <- d_sph_uniform_beta_mixture_s2(
    x, theta$mu, theta$weight_uniform, theta$alpha, theta$beta, log = TRUE
  )
  dots <- pmin(pmax(as.numeric(x %*% theta$mu), -1), 1)
  all <- data.frame(
    observation = seq_len(nrow(x)), log_density = values,
    y = (dots + 1) / 2, distance_to_mu = acos(dots),
    distance_to_antipode = acos(-dots), fit = label
  )
  q <- stats::quantile(values, probs = c(0, .5, .9, .95, .99, .995, .999, 1), names = FALSE)
  excess <- values + log(4 * pi)
  ordered_excess <- sort(excess, decreasing = TRUE)
  total_excess <- sum(excess)
  summary <- data.frame(
    fit = label, probability = c(0, .5, .9, .95, .99, .995, .999, 1),
    log_density = q, total_loglik = sum(values), total_excess_over_uniform = total_excess,
    top1_excess_fraction = ordered_excess[[1L]] / total_excess,
    top5_excess_fraction = sum(ordered_excess[seq_len(5L)]) / total_excess,
    top10_excess_fraction = sum(ordered_excess[seq_len(10L)]) / total_excess
  )
  list(all = all, top = all[order(all$log_density, decreasing = TRUE)[seq_len(min(20L, nrow(all)))], ],
       summary = summary)
}

loglik_target_stable <- function(mu, weight, alpha, beta, x, target_index,
                                 endpoint, delta) {
  y <- (pmin(pmax(as.numeric(x %*% mu), -1), 1) + 1) / 2
  log_y <- log(y)
  log_one_minus_y <- log1p(-y)
  if (endpoint == "beta") {
    log_y[[target_index]] <- 2 * log(cos(delta / 2))
    log_one_minus_y[[target_index]] <- 2 * log(sin(delta / 2))
  } else {
    log_y[[target_index]] <- 2 * log(sin(delta / 2))
    log_one_minus_y[[target_index]] <- 2 * log(cos(delta / 2))
  }
  log_beta <- (alpha - 1) * log_y + (beta - 1) * log_one_minus_y - lbeta(alpha, beta)
  sum(-log(4 * pi) + rotational_logsumexp2(log(weight), log1p(-weight) + log_beta))
}

profile_singularity <- function(theta, x, target_index, endpoint) {
  target <- if (endpoint == "beta") x[target_index, ] else -x[target_index, ]
  initial <- angle(theta$mu, target)
  tangent <- theta$mu - sum(theta$mu * target) * target
  if (sqrt(sum(tangent^2)) < 1e-12) tangent <- tangent_basis(target)[, 1L]
  tangent <- tangent / sqrt(sum(tangent^2))
  distances <- unique(exp(seq(log(max(initial, 1e-7)), log(1e-14), length.out = 70L)))
  q_start <- c(stats::qlogis(theta$weight_uniform), log(theta$alpha), log(theta$beta))
  rows <- lapply(distances, function(delta) {
    mu <- cos(delta) * target + sin(delta) * tangent
    fixed <- loglik_target_stable(mu, theta$weight_uniform, theta$alpha, theta$beta,
                                  x, target_index, endpoint, delta)
    opt <- stats::optim(q_start, function(q) -loglik_target_stable(
      mu, rotational_bounded_weight(q[[1L]], 0.01),
      rotational_positive_parameter(q[[2L]], 0.05, 1e3),
      rotational_positive_parameter(q[[3L]], 0.05, 1e3),
      x, target_index, endpoint, delta
    ), method = "BFGS", control = list(maxit = 200L, reltol = 1e-9))
    q_start <<- opt$par
    data.frame(
      observation = target_index, endpoint = endpoint, distance = delta,
      fixed_loglik = fixed, conditional_loglik = -opt$value,
      conditional_weight = rotational_bounded_weight(opt$par[[1L]], 0.01),
      conditional_alpha = rotational_positive_parameter(opt$par[[2L]], 0.05, 1e3),
      conditional_beta = rotational_positive_parameter(opt$par[[3L]], 0.05, 1e3),
      conditional_convergence = opt$convergence
    )
  })
  do.call(rbind, rows)
}

bootstrap_refits <- function(x, reported, B, seed, cores) {
  weights <- generate_multiplier_matrix(
    B, nrow(x), resolve_multiplier_spec(NULL), seed = seed
  )
  weights <- weights / rowMeans(weights)
  worker <- function(b) {
    control <- bootstrap_control
    control$uniform_beta_mixture_start_theta <- reported
    control$uniform_beta_mixture_warm_start_only <- TRUE
    control$uniform_beta_mixture_n_starts <- 1L
    control$uniform_beta_mixture_optim_control <-
      control$uniform_beta_mixture_bootstrap_optim_control %||% list(maxit = 80L, reltol = 1e-6)
    cw <- capture_warnings(uniform_beta_mixture_mle_s2_weighted(
      x, weights = weights[b, ], control = control
    ))
    if (inherits(cw$value, "audit_error")) {
      return(data.frame(replicate = b, convergence = NA_integer_, weighted_loglik = NA_real_,
                        function_evaluations = NA_integer_, gradient_evaluations = NA_integer_,
                        optim_message = cw$value$message,
                        mu_1 = NA_real_, mu_2 = NA_real_, mu_3 = NA_real_, weight_uniform = NA_real_,
                        alpha = NA_real_, beta = NA_real_, distance_to_nearest = NA_real_,
                        antipode_distance_to_nearest = NA_real_, angle_to_reported = NA_real_,
                        angle_to_reported_equivalent = NA_real_,
                        transformed_parameter_distance = NA_real_, same_reported_basin = FALSE,
                        exactly_unchanged_from_start = FALSE,
                        warnings = paste(c(cw$warnings, cw$value$message), collapse = " | ")))
    }
    z <- cw$value
    near <- nearest_distances(z$mu, x)
    rd <- reported_distances(z, reported)
    data.frame(
      replicate = b, convergence = z$opt$convergence,
      weighted_loglik = nrow(x) * (z$loglik - log(4 * pi)),
      function_evaluations = as.integer(z$opt$counts[["function"]]),
      gradient_evaluations = as.integer(z$opt$counts[["gradient"]]),
      optim_message = z$opt$message %||% "",
      mu_1 = z$mu[[1L]], mu_2 = z$mu[[2L]], mu_3 = z$mu[[3L]],
      weight_uniform = z$weight_uniform, alpha = z$alpha, beta = z$beta,
      distance_to_nearest = near[[1L]], antipode_distance_to_nearest = near[[2L]],
      angle_to_reported = rd[["signed_angle"]],
      angle_to_reported_equivalent = rd[["equivalent_angle"]],
      transformed_parameter_distance = rd[["equivalent_parameter"]],
      same_reported_basin = (rd[["signed_angle"]] < 0.15 && rd[["direct_parameter"]] < 1) ||
        (rd[["antipodal_angle"]] < 0.15 && rd[["swapped_parameter"]] < 1),
      exactly_unchanged_from_start = rd[["signed_angle"]] < 1e-12 && rd[["direct_parameter"]] < 1e-10,
      warnings = paste(cw$warnings, collapse = " | ")
    )
  }
  if (.Platform$OS.type == "unix" && cores > 1L) {
    do.call(rbind, parallel::mclapply(seq_len(B), worker, mc.cores = min(cores, B), mc.preschedule = TRUE))
  } else do.call(rbind, lapply(seq_len(B), worker))
}

write_plot <- function(path, expr) {
  grDevices::png(path, width = 1100, height = 800, res = 130)
  on.exit(grDevices::dev.off())
  force(expr)
}

data <- load_comets_real_data(finite_normals = "both")
samples <- list(short = as.matrix(data$short$normal), long = as.matrix(data$long$normal))
reported <- lapply(reference_dirs, function(d) read_theta(file.path(d, "theta_hat.csv")))

all_reproduction <- list()
for (dataset in names(samples)) {
  message("[audit] ", dataset)
  x <- samples[[dataset]]
  theta_ref <- reported[[dataset]]
  dataset_dir <- file.path(output_root, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  reproduction <- capture_warnings(uniform_beta_mixture_mle_s2_weighted(
    x, control = production_control
  ))
  if (inherits(reproduction$value, "audit_error")) stop(reproduction$value$message)
  theta_reproduced <- reproduction$value
  ll <- surface_loglik(theta_reproduced, x)
  production_objective <- function(p) {
    z <- uniform_beta_mixture_unpack_par(p, production_control)
    -uniform_beta_mixture_weighted_loglik_s2(
      mu = z$mu, weight_uniform = z$weight_uniform,
      alpha = z$alpha, beta = z$beta, x = x
    )
  }
  prod_grad <- vapply(seq_along(reproduction$value$opt$par), function(j) {
    h <- 1e-7; step <- rep(0, length(reproduction$value$opt$par)); step[[j]] <- h
    (production_objective(reproduction$value$opt$par + step) -
       production_objective(reproduction$value$opt$par - step)) / (2 * h)
  }, numeric(1L))
  rep_row <- data.frame(
    dataset = dataset, n = nrow(x), loglik = ll, AIC = -2 * ll + 10,
    BIC = -2 * ll + 5 * log(nrow(x)), convergence = theta_reproduced$opt$convergence,
    production_parameter_gradient_norm = sqrt(sum(prod_grad^2)),
    mu_1 = theta_reproduced$mu[[1L]], mu_2 = theta_reproduced$mu[[2L]],
    mu_3 = theta_reproduced$mu[[3L]], weight_uniform = theta_reproduced$weight_uniform,
    alpha = theta_reproduced$alpha, beta = theta_reproduced$beta,
    max_abs_parameter_difference_from_saved = max(abs(theta_vector(theta_reproduced) - theta_vector(theta_ref))),
    warnings = paste(reproduction$warnings, collapse = " | ")
  )
  all_reproduction[[dataset]] <- rep_row
  utils::write.csv(rep_row, file.path(dataset_dir, "reproduction.csv"), row.names = FALSE)

  local <- local_diagnostics(theta_ref, x, n_local_perturb, audit_seed + if (dataset == "short") 1L else 2L)
  utils::write.csv(data.frame(parameter = c("mu_tangent_1", "mu_tangent_2", "logit_w", "log_alpha", "log_beta"),
                              gradient = local$gradient),
                   file.path(dataset_dir, "local_gradient.csv"), row.names = FALSE)
  utils::write.csv(local$gradient_by_step,
                   file.path(dataset_dir, "local_gradient_by_step.csv"), row.names = FALSE)
  utils::write.csv(local$hessian, file.path(dataset_dir, "local_hessian.csv"), row.names = FALSE)
  utils::write.csv(data.frame(eigenvalue = local$eigenvalues,
                              negative_hessian_condition = local$negative_hessian_condition,
                              finite_difference_step = local$hessian_step),
                   file.path(dataset_dir, "local_hessian_eigenvalues.csv"), row.names = FALSE)
  utils::write.csv(local$perturbations, file.path(dataset_dir, "local_perturbations.csv"), row.names = FALSE)

  starts <- make_multistarts(x, theta_ref, n_multistart, audit_seed + if (dataset == "short") 11L else 12L)
  fits <- lapply(seq_along(starts), function(i) single_start_fit(starts[[i]], x))
  multistart <- do.call(rbind, Map(function(fit, start, i) fit_row(fit, start, i, x, theta_ref),
                                  fits, starts, seq_along(starts)))
  multistart <- assign_basins(multistart)
  utils::write.csv(multistart, file.path(dataset_dir, "multistart_runs.csv"), row.names = FALSE)
  basin_summary <- aggregate(run ~ basin, multistart[is.finite(multistart$loglik), ], length)
  names(basin_summary)[[2L]] <- "frequency"
  best_by_basin <- aggregate(loglik ~ basin, multistart[is.finite(multistart$loglik), ], max)
  basin_summary <- merge(basin_summary, best_by_basin, by = "basin")
  utils::write.csv(basin_summary[order(basin_summary$loglik, decreasing = TRUE), ],
                   file.path(dataset_dir, "basin_summary.csv"), row.names = FALSE)

  best_i <- which.max(multistart$loglik)
  best_theta <- list(mu = as.numeric(multistart[best_i, c("mu_1", "mu_2", "mu_3")]),
                     weight_uniform = multistart$weight_uniform[[best_i]],
                     alpha = multistart$alpha[[best_i]], beta = multistart$beta[[best_i]])
  contrib <- lapply(list(reported = theta_ref, multistart_best = best_theta), function(th)
    contribution_tables(th, x, if (identical(th, theta_ref)) "reported" else "multistart_best"))
  utils::write.csv(do.call(rbind, lapply(contrib, `[[`, "summary")),
                   file.path(dataset_dir, "contribution_summary.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, lapply(contrib, `[[`, "top")),
                   file.path(dataset_dir, "top_contributions.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, lapply(contrib, `[[`, "all")),
                   file.path(dataset_dir, "all_contributions.csv"), row.names = FALSE)

  dots <- as.numeric(x %*% theta_ref$mu)
  targets <- list(beta = which.max(dots))
  if (dataset == "long") targets$alpha <- which.min(dots)
  profiles <- do.call(rbind, lapply(names(targets), function(endpoint)
    profile_singularity(theta_ref, x, targets[[endpoint]], endpoint)))
  profiles$reported_loglik <- surface_loglik(theta_ref, x)
  utils::write.csv(profiles, file.path(dataset_dir, "singularity_profiles.csv"), row.names = FALSE)

  bootstrap_seed <- if (dataset == "short") 20261602L else 20261603L
  boot <- bootstrap_refits(x, theta_ref, n_bootstrap, bootstrap_seed, n_cores)
  boot <- assign_basins(transform(boot, run = replicate, loglik = weighted_loglik,
                                 AIC = NA_real_, BIC = NA_real_, message = ""))
  boot$run <- NULL; boot$loglik <- NULL; boot$AIC <- NULL; boot$BIC <- NULL; boot$message <- NULL
  utils::write.csv(boot, file.path(dataset_dir, "bootstrap_refits.csv"), row.names = FALSE)

  write_plot(file.path(dataset_dir, "multistart_loglik.png"), {
    graphics::plot(multistart$distance_to_nearest, multistart$loglik, log = "x", pch = 19,
                   col = ifelse(multistart$same_reported_basin, "#2166ac", "#b2182b"),
                   xlab = "Distance from fitted mu to nearest observation (radians)",
                   ylab = "Log-likelihood", main = paste(dataset, "multistart"))
    graphics::abline(h = surface_loglik(theta_ref, x), lty = 2)
  })
  write_plot(file.path(dataset_dir, "singularity_profiles.png"), {
    graphics::plot(range(profiles$distance), range(c(profiles$fixed_loglik, profiles$conditional_loglik)),
                   log = "x", type = "n", xlab = "Distance to target (radians)",
                   ylab = "Log-likelihood", main = paste(dataset, "singularity profiles"))
    for (endpoint in unique(profiles$endpoint)) {
      z <- profiles[profiles$endpoint == endpoint, ]
      graphics::lines(z$distance, z$fixed_loglik, col = if (endpoint == "beta") 2 else 4)
      graphics::lines(z$distance, z$conditional_loglik, col = if (endpoint == "beta") 2 else 4, lty = 2)
    }
    graphics::abline(h = surface_loglik(theta_ref, x), lty = 3)
  })
}

utils::write.csv(do.call(rbind, all_reproduction), file.path(output_root, "reproduction_all.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
writeLines(c(
  paste("created_at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("n_multistart:", n_multistart), paste("n_local_perturb:", n_local_perturb),
  paste("n_bootstrap:", n_bootstrap), paste("n_cores:", n_cores),
  "production objective: mean projected log-density, BFGS, 12 starts, maxit=350, reltol=1e-9",
  "surface AIC/BIC: full S2 log-density with k=5",
  "bootstrap refits: Exp(1) row-normalized weights; observed-fit warm start only; maxit=80; reltol=1e-6",
  "production endpoint clamp: eps=1e-12; stable singularity profiles replace the target contribution analytically"
), file.path(output_root, "manifest.txt"))
message("Audit artifacts: ", output_root)
