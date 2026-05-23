source("utils.R")
suppressPackageStartupMessages(library(ggplot2))

dir.create("output/gaussian_process_normal", recursive = TRUE, showWarnings = FALSE)

simulate_empirical_supremum <- function(n, M, t_grid, r_sampler, p_cdf, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- numeric(M)
  F_t <- p_cdf(t_grid)

  for (m in seq_len(M)) {
    x <- sort(r_sampler(n))
    F_n_t <- findInterval(t_grid, x, rightmost.closed = TRUE) / n
    out[m] <- max(abs(sqrt(n) * (F_n_t - F_t)))
  }

  out
}

simulate_limit_supremum_ecdf <- function(M, F_t, seed = NULL) {
  cov_matrix <- outer(F_t, F_t, function(a, b) pmin(a, b) - a * b)
  simulate_limit_gaussian(cov_matrix, M = M, seed = seed)
}

make_hist_plot <- function(empirical_values, limit_values, title_txt, subtitle_txt, out_path) {
  df_emp <- data.frame(values = as.numeric(empirical_values), process = "Empirical (n=10)")
  df_lim <- data.frame(values = as.numeric(limit_values), process = "Limit (G)")
  df <- rbind(df_emp, df_lim)

  p <- ggplot(df, aes(x = values, fill = process)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 70,
      alpha = 0.45,
      position = "identity",
      color = NA
    ) +
    scale_fill_manual(values = c("Empirical (n=10)" = "#F8766D", "Limit (G)" = "#000000")) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Supremum of sqrt(n) (F_n(t) - F(t))",
      y = "Density",
      fill = "Process"
    ) +
    theme_minimal(base_size = 14)

  ggsave(out_path, plot = p, width = 12, height = 8, dpi = 300)
}

run_case <- function(case_name, n, M, t_grid, r_sampler, p_cdf, out_file, empirical_seed, limit_seed) {
  cat("===", case_name, "===\n")
  cat("n =", n, "| M =", M, "| grid size =", length(t_grid), "\n")
  cat("empirical_seed =", empirical_seed, "| limit_seed =", limit_seed, "\n")

  empirical_values <- simulate_empirical_supremum(
    n = n,
    M = M,
    t_grid = t_grid,
    r_sampler = r_sampler,
    p_cdf = p_cdf,
    seed = empirical_seed
  )

  limit_values <- simulate_limit_supremum_ecdf(
    M = M,
    F_t = p_cdf(t_grid),
    seed = limit_seed
  )

  make_hist_plot(
    empirical_values = empirical_values,
    limit_values = limit_values,
    title_txt = paste0(case_name, ": Histograms only (grid=5, seeded)"),
    subtitle_txt = paste0("n=", n, ", M=", M, ", |t-grid|=", length(t_grid)),
    out_path = out_file
  )

  cat("SAVED:", out_file, "\n\n")
}

# Shared settings
n <- 10
M <- 10000

# 1) Uniform(0,1) with t-grid of 5 points
# Use interior points to avoid exact 0/1 boundaries.
t_grid_uniform <- seq(0.1, 0.9, length.out = 5)
out_uniform <- "output/gaussian_process_normal/prueba_ecdf_uniform_0_1_n10_M10000_hist_grid5_seeded.png"

run_case(
  case_name = "Uniform(0,1)",
  n = n,
  M = M,
  t_grid = t_grid_uniform,
  r_sampler = function(k) runif(k, min = 0, max = 1),
  p_cdf = function(t) punif(t, min = 0, max = 1),
  out_file = out_uniform,
  empirical_seed = 123,
  limit_seed = 123
)

# 2) Normal(mu=3, var=0.01 => sd=0.1) with t-grid of 5 points
mu <- 3
sd <- 0.1
u_grid <- seq(0.1, 0.9, length.out = 5)
t_grid_normal <- qnorm(u_grid, mean = mu, sd = sd)
out_normal <- "output/gaussian_process_normal/prueba_ecdf_normal_mu3_var001_n10_M10000_hist_grid5_seeded.png"

run_case(
  case_name = "Normal(mu=3, var=0.01)",
  n = n,
  M = M,
  t_grid = t_grid_normal,
  r_sampler = function(k) rnorm(k, mean = mu, sd = sd),
  p_cdf = function(t) pnorm(t, mean = mu, sd = sd),
  out_file = out_normal,
  empirical_seed = 456,
  limit_seed = 456
)

cat("DONE\n")
