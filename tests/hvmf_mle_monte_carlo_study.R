source("utils.R")

hyperbolic_geodesic_distance_h2 <- function(x, y, tol = 1e-12) {
  acosh(max(-minkowski_inner_product(x, y), 1 + tol))
}

run_hvmf_mle_monte_carlo_study <- function(kappas = c(5, 25, 50),
                                            n_values = c(50, 100, 200),
                                            M = 100L,
                                            mu = c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2)),
                                            output_csv = file.path("output", "hvmf_mle_monte_carlo_summary.csv")) {
  rows <- list()
  index <- 1L
  for (kappa in kappas) {
    for (n in n_values) {
      estimates <- t(vapply(seq_len(M), function(rep_id) {
        fit <- hvmf_mle_h2(rhvmf_h2_polar(n = n, mu = mu, kappa = kappa))
        c(kappa = fit$kappa, distance = hyperbolic_geodesic_distance_h2(fit$mu, mu))
      }, numeric(2)))
      rows[[index]] <- data.frame(
        kappa_true = kappa,
        n = n,
        M = M,
        mean_kappa_hat = mean(estimates[, "kappa"]),
        sd_kappa_hat = stats::sd(estimates[, "kappa"]),
        mean_center_distance = mean(estimates[, "distance"]),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  summary_df <- do.call(rbind, rows)
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary_df, output_csv, row.names = FALSE)
  summary_df
}

if (sys.nframe() == 0L) {
  print(run_hvmf_mle_monte_carlo_study())
}
