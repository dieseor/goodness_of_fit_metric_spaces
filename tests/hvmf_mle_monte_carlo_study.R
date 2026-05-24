source("utils.R")

hyperbolic_geodesic_distance_h2 <- function(x, y, tol = 1e-12) {
  inner <- minkowski_inner_product(x, y)
  value <- max(-inner, 1 + tol)
  acosh(value)
}

true_typeiv_hvmf_center <- function() {
  t0 <- qnorm(0.25, mean = 0, sd = 1 / 4)
  spatial_direction <- c(1 / sqrt(2), 1 / sqrt(2))

  c(
    cosh(abs(t0)),
    sinh(abs(t0)) * sign(t0) * spatial_direction[[1]],
    sinh(abs(t0)) * sign(t0) * spatial_direction[[2]]
  )
}

summarize_hvmf_mle_folder <- function(folder, kappa_true) {
  files <- sort(list.files(folder, pattern = "\\.csv$", full.names = TRUE))
  if (length(files) == 0L) {
    stop(sprintf("No CSV files found in `%s`.", folder))
  }

  xi_true <- true_typeiv_hvmf_center()
  estimates <- lapply(files, function(path) {
    dat <- read.csv(path)
    xyz <- as.matrix(dat[, c("V1", "V2", "V3")])
    fit <- hvmf_mle_h2(xyz)

    c(
      xi1 = fit$xi[[1]],
      xi2 = fit$xi[[2]],
      xi3 = fit$xi[[3]],
      kappa = fit$kappa,
      angular_error = hyperbolic_geodesic_distance_h2(fit$xi, xi_true)
    )
  })

  estimates <- do.call(rbind, estimates)
  estimates <- as.data.frame(estimates)

  data.frame(
    folder = folder,
    n_replicates = nrow(estimates),
    n = nrow(read.csv(files[[1]])),
    kappa_true = kappa_true,
    xi_true_1 = xi_true[[1]],
    xi_true_2 = xi_true[[2]],
    xi_true_3 = xi_true[[3]],
    mean_kappa_hat = mean(estimates$kappa),
    sd_kappa_hat = sd(estimates$kappa),
    bias_kappa_hat = mean(estimates$kappa) - kappa_true,
    rmse_kappa_hat = sqrt(mean((estimates$kappa - kappa_true)^2)),
    mean_xi1_hat = mean(estimates$xi1),
    mean_xi2_hat = mean(estimates$xi2),
    mean_xi3_hat = mean(estimates$xi3),
    rmse_xi1_hat = sqrt(mean((estimates$xi1 - xi_true[[1]])^2)),
    rmse_xi2_hat = sqrt(mean((estimates$xi2 - xi_true[[2]])^2)),
    rmse_xi3_hat = sqrt(mean((estimates$xi3 - xi_true[[3]])^2)),
    mean_angular_error = mean(estimates$angular_error),
    median_angular_error = median(estimates$angular_error),
    q90_angular_error = unname(stats::quantile(estimates$angular_error, 0.9)),
    max_angular_error = max(estimates$angular_error)
  )
}

run_hvmf_mle_monte_carlo_study <- function(base_dir = file.path("data", "hvmf_typeiv_calibration"),
                                            kappas = c(50, 200),
                                            n_values = c(50, 100, 200),
                                            output_csv = file.path("output", "hvmf_mle_monte_carlo_summary.csv")) {
  summaries <- list()

  for (kappa in kappas) {
    for (n in n_values) {
      folder <- file.path(base_dir, paste0("kappa", kappa), paste0("n", n))
      summaries[[paste(kappa, n, sep = "_")]] <- summarize_hvmf_mle_folder(folder, kappa_true = kappa)
    }
  }

  summary_df <- do.call(rbind, summaries)
  rownames(summary_df) <- NULL

  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_df, output_csv, row.names = FALSE)

  summary_df
}

if (sys.nframe() == 0L) {
  summary_df <- run_hvmf_mle_monte_carlo_study()
  print(summary_df)
}
