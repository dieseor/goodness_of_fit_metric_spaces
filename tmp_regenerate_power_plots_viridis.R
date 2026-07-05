safe_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

plot_title_suffix <- function(df) {
  first <- df[1L, , drop = FALSE]

  if (identical(first$scenario, "vmf_s2_antipodal")) {
    return("vMF on S^2: antipodal mixture")
  }
  if (identical(first$scenario, "vmf_s2_non_antipodal")) {
    return(sprintf("vMF on S^2: non-antipodal mixture, gamma = %s deg", first$gamma_deg))
  }
  if (identical(first$scenario, "logistic_gaussian_simplex_d3")) {
    return(sprintf(
      "Simplex D = 3: LG vs Dirichlet (%s, %s, %s)",
      format(first$dirichlet_alpha1, trim = TRUE),
      format(first$dirichlet_alpha2, trim = TRUE),
      format(first$dirichlet_alpha3, trim = TRUE)
    ))
  }

  sprintf("%s / %s", first$scenario, first$alternative)
}

plot_group_slug <- function(df) {
  first <- df[1L, , drop = FALSE]
  pieces <- c(first$scenario, first$alternative)

  if (!is.na(first$gamma_deg)) {
    pieces <- c(pieces, sprintf("gamma_%s", first$gamma_deg))
  }
  if (!is.na(first$dirichlet_a)) {
    pieces <- c(pieces, sprintf("a_%s", first$dirichlet_a))
  } else if (!is.na(first$dirichlet_alpha1) && !is.na(first$dirichlet_alpha2) && !is.na(first$dirichlet_alpha3)) {
    pieces <- c(
      pieces,
      sprintf(
        "alpha_%s_%s_%s",
        format(first$dirichlet_alpha1, trim = TRUE),
        format(first$dirichlet_alpha2, trim = TRUE),
        format(first$dirichlet_alpha3, trim = TRUE)
      )
    )
  }

  safe_slug(paste(pieces, collapse = "_"))
}

save_power_plot <- function(summary_subset, stat_name, file_path) {
  n_values <- sort(unique(summary_subset$n))
  beta_values <- sort(unique(summary_subset$beta))
  colors <- viridisLite::viridis(length(n_values), option = "D", begin = 0.15, end = 0.9)
  y_column <- if (identical(stat_name, "ks")) "power_ks_005" else "power_cvm_005"
  calibration_band <- c(0.0365, 0.0635)

  grDevices::png(filename = file_path, width = 1800, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(6.2, 8.2, 2.2, 1.2))

  plot(
    NA,
    xlim = range(beta_values),
    ylim = c(0, 1),
    xlab = "",
    ylab = expression(P[M](p <= 0.05)),
    cex.axis = 1.5,
    cex.lab = 1.5
  )
  usr <- graphics::par("usr")
  x_pad <- 0.01 * diff(usr[1:2])
  grid(col = "grey78")
  graphics::rect(
    xleft = usr[[1L]] - x_pad,
    ybottom = calibration_band[[1L]],
    xright = usr[[2L]] + x_pad,
    ytop = calibration_band[[2L]],
    col = grDevices::adjustcolor("grey75", alpha.f = 0.35),
    border = NA
  )
  graphics::abline(h = 0.05, lty = 2, lwd = 2, col = "black")
  graphics::mtext(expression(beta), side = 1, line = 2.5, cex = 2)

  for (i in seq_along(n_values)) {
    n_i <- n_values[[i]]
    df_i <- summary_subset[summary_subset$n == n_i, , drop = FALSE]
    df_i <- df_i[order(df_i$beta), , drop = FALSE]
    lines(df_i$beta, df_i[[y_column]], col = colors[[i]], lwd = 3)
    points(df_i$beta, df_i[[y_column]], col = colors[[i]], pch = 19, cex = 1.1)
  }

  legend(
    "topleft",
    legend = sprintf("n = %d", n_values),
    col = colors,
    lwd = 3,
    pch = 19,
    bty = "n",
    cex = 2,
    pt.cex = 2
  )
}

save_alpha_curve_plot <- function(alpha_subset, stat_name, file_path) {
  n_values <- sort(unique(alpha_subset$n))
  colors <- viridisLite::viridis(length(n_values), option = "D", begin = 0.15, end = 0.9)
  y_column <- if (identical(stat_name, "ks")) "rejection_prob_ks" else "rejection_prob_cvm"
  beta_value <- unique(alpha_subset$beta)

  grDevices::png(filename = file_path, width = 1800, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(6.2, 8.2, 2.2, 1.2))

  plot(
    c(0, 1),
    c(0, 1),
    type = "n",
    xlab = "",
    ylab = expression(P[M](p <= alpha)),
    cex.axis = 1.5,
    cex.lab = 1.5
  )
  grid(col = "grey85")
  abline(0, 1, lty = 2, lwd = 2, col = "black")
  graphics::mtext(expression(alpha), side = 1, line = 2.5, cex = 2)

  for (i in seq_along(n_values)) {
    n_i <- n_values[[i]]
    df_i <- alpha_subset[alpha_subset$n == n_i, , drop = FALSE]
    df_i <- df_i[order(df_i$alpha), , drop = FALSE]
    lines(df_i$alpha, df_i[[y_column]], col = colors[[i]], lwd = 3)
  }

  legend(
    "bottomright",
    legend = sprintf("n = %d", n_values),
    col = colors,
    lwd = 3,
    bty = "n",
    cex = 2
  )
}

save_all_plots <- function(summary_df, alpha_curve_df, output_dir) {
  plots_dir <- file.path(output_dir, "plots")
  power_dir <- file.path(plots_dir, "power_vs_beta")
  alpha_dir <- file.path(plots_dir, "alpha_curves")

  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(summary_df) > 0L) {
    group_key <- paste(
      summary_df$scenario,
      summary_df$alternative,
      ifelse(is.na(summary_df$gamma_deg), "NA", format(summary_df$gamma_deg, trim = TRUE)),
      ifelse(is.na(summary_df$dirichlet_a), "NA", format(summary_df$dirichlet_a, trim = TRUE)),
      sep = "|"
    )
    summary_groups <- split(summary_df, group_key, drop = TRUE)

    for (df in summary_groups) {
      group_slug <- plot_group_slug(df)
      save_power_plot(df, stat_name = "ks", file_path = file.path(power_dir, sprintf("%s_power_ks.png", group_slug)))
      save_power_plot(df, stat_name = "cvm", file_path = file.path(power_dir, sprintf("%s_power_cvm.png", group_slug)))
    }
  }

  if (nrow(alpha_curve_df) > 0L) {
    group_key <- paste(
      alpha_curve_df$scenario,
      alpha_curve_df$alternative,
      ifelse(is.na(alpha_curve_df$gamma_deg), "NA", format(alpha_curve_df$gamma_deg, trim = TRUE)),
      ifelse(is.na(alpha_curve_df$dirichlet_a), "NA", format(alpha_curve_df$dirichlet_a, trim = TRUE)),
      alpha_curve_df$beta,
      sep = "|"
    )
    alpha_groups <- split(alpha_curve_df, group_key, drop = TRUE)

    for (df in alpha_groups) {
      group_slug <- sprintf("%s_beta_%s", plot_group_slug(df), safe_slug(format(df$beta[[1L]], trim = TRUE)))
      save_alpha_curve_plot(df, stat_name = "ks", file_path = file.path(alpha_dir, sprintf("%s_alpha_ks.png", group_slug)))
      save_alpha_curve_plot(df, stat_name = "cvm", file_path = file.path(alpha_dir, sprintf("%s_alpha_cvm.png", group_slug)))
    }
  }
}

summary_df <- utils::read.csv("simulation_results/power_mixtures_pilot/power_summary_005.csv", stringsAsFactors = FALSE)
alpha_curve_df <- utils::read.csv("simulation_results/power_mixtures_pilot/alpha_curve_summary.csv", stringsAsFactors = FALSE)
save_all_plots(summary_df, alpha_curve_df, "simulation_results/power_mixtures_pilot")
