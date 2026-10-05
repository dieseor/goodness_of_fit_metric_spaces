#!/usr/bin/env Rscript

# Rebuild the 12 active manuscript PNGs without rerunning any simulations.
# Bandwidth adjustment multipliers are global by sample size, so every figure
# necessarily uses the same multiplier for a given n.

script_path <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
if (length(script_path) != 1L || !nzchar(script_path)) stop("Run this script with Rscript.")
project_root <- dirname(dirname(normalizePath(script_path)))
setwd(project_root)

suppressPackageStartupMessages(library(ggplot2))
`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else TRUE
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- args$input_dir %||%
  "output/convergence/gaussian_process/vmf_paper_grid10x10/simulation_results"
output_dir <- args$output_dir %||%
  "output/convergence/gaussian_process/vmf_paper_grid10x10/plots_bcv"
M <- as.integer(args$M %||% 10000L)
h0 <- args$h0 %||% "both"
if (!h0 %in% c("simple", "composite", "both")) {
  stop("--h0 must be one of: simple, composite, both.")
}
adjustments <- c(
  `50` = as.numeric(args$adjust_50 %||% 1),
  `100` = as.numeric(args$adjust_100 %||% 1),
  `500` = as.numeric(args$adjust_500 %||% 1)
)
if (any(!is.finite(adjustments)) || any(adjustments <= 0)) {
  stop("All bandwidth adjustment multipliers must be positive finite numbers.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rds_paths <- sort(Sys.glob(file.path(input_dir, sprintf("*M%d_grid10x10.rds", M))))
if (length(rds_paths) != 6L) {
  stop(sprintf("Expected six M=%d simulation RDS files in %s; found %d.", M, input_dir, length(rds_paths)))
}
if (h0 == "simple") {
  rds_paths <- rds_paths[startsWith(basename(rds_paths), "simple_")]
} else if (h0 == "composite") {
  rds_paths <- rds_paths[startsWith(basename(rds_paths), "comp_")]
}
expected_selected <- if (h0 == "both") 6L else 3L
if (length(rds_paths) != expected_selected) {
  stop(sprintf("Expected %d RDS files for h0=%s; found %d.", expected_selected, h0, length(rds_paths)))
}

kde_axis_size <- 19 * 1.5 * 1.2
kde_legend_size <- 19 * 2
qq_scale <- 8 / 12  # Both PNG types have the same final width in the paper.

make_density_plot <- function(result, show_legend, is_composite) {
  n_values <- as.integer(result$n_values)
  empirical_colors <- scales::hue_pal()(length(n_values))
  labels <- c(paste0("n=", n_values), "G")
  limit_label <- if (is_composite) {
    quote("‖" * widetilde("𝔾")[0] * "‖"[infinity])
  } else {
    quote("‖" * "𝔾"[0] * "‖"[infinity])
  }
  legend_labels <- as.expression(c(lapply(n_values, function(n) bquote(n == .(n))), list(limit_label)))
  colors <- setNames(c(empirical_colors, "#000000"), labels)
  linetypes <- setNames(c(rep("solid", length(n_values)), "dashed"), labels)

  density_data <- do.call(rbind, lapply(n_values, function(n) {
    data.frame(values = result$empirical_data[[as.character(n)]], label = paste0("n=", n))
  }))
  density_data <- rbind(
    density_data,
    data.frame(values = result$limit_values_vec[is.finite(result$limit_values_vec)], label = "G")
  )
  density_data$label <- factor(density_data$label, levels = labels)

  plot <- ggplot()
  for (n in n_values) {
    label <- paste0("n=", n)
    plot <- plot + geom_density(
      data = density_data[density_data$label == label, , drop = FALSE],
      aes(x = values, color = label, linetype = label),
      bw = "bcv",
      adjust = unname(adjustments[[as.character(n)]]),
      fill = NA,
      linewidth = 1,
      trim = TRUE,
      show.legend = TRUE,
      key_glyph = draw_key_path
    )
  }
  plot +
    geom_density(
      data = density_data[density_data$label == "G", , drop = FALSE],
      aes(x = values, color = label, linetype = label),
      bw = "bcv",
      adjust = 1,
      fill = NA,
      linewidth = 1,
      trim = TRUE,
      show.legend = TRUE,
      key_glyph = draw_key_path
    ) +
    scale_color_manual(values = colors, breaks = labels, labels = legend_labels, drop = FALSE) +
    scale_linetype_manual(values = linetypes, breaks = labels, drop = FALSE) +
    labs(x = "Supremum of the process", y = "Density", color = "Process", linetype = "Process") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = kde_axis_size),
      axis.title = element_text(size = kde_axis_size),
      legend.position = if (show_legend) c(0.98, 0.98) else "none",
      legend.justification = c("right", "top"),
      legend.text = element_text(size = kde_legend_size, family = "STIX Two Math",
                                 margin = margin(t = 4, b = 4)),
      legend.title = element_text(size = kde_legend_size, margin = margin(b = 8)),
      legend.key.width = grid::unit(0.85, "in"),
      legend.key.height = grid::unit(0.6, "in")
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", "21"),
                                               fill = NA, alpha = 1,
                                               linewidth = c(2.5, 2.5, 2.5, 2))),
      linetype = "none"
    )
}

for (path in rds_paths) {
  result <- readRDS(path)
  prefix <- tools::file_path_sans_ext(basename(path))
  show_legend <- grepl("kappa0.5_", prefix, fixed = TRUE)
  is_composite <- startsWith(prefix, "comp_")
  ggsave(file.path(output_dir, paste0(prefix, ".png")), make_density_plot(result, show_legend, is_composite),
         width = 12, height = 8, dpi = 300)
  qq_plot <- result$qq_plot + labs(color = expression(n)) + theme(
    axis.text = element_text(size = kde_axis_size * qq_scale),
    axis.title = element_text(size = kde_axis_size * qq_scale),
    legend.position = if (show_legend) c(0.02, 0.98) else "none",
    legend.justification = c("left", "top"),
    legend.text = element_text(size = kde_legend_size * qq_scale, family = "STIX Two Math",
                               margin = margin(t = 4 * qq_scale, b = 4 * qq_scale)),
    legend.title = element_text(size = kde_legend_size * qq_scale, family = "STIX Two Math",
                                margin = margin(b = 8 * qq_scale)),
    legend.key.height = grid::unit(0.6 * qq_scale, "in")
  ) + guides(color = guide_legend(override.aes = list(size = 3)))
  ggsave(file.path(output_dir, paste0("qq_", prefix, ".png")), qq_plot,
         width = 8, height = 6, dpi = 300)
}

utils::write.csv(
  data.frame(n = as.integer(names(adjustments)), bandwidth_adjust = unname(adjustments)),
  file.path(output_dir, "bandwidth_adjustments.csv"),
  row.names = FALSE
)
message("Created 12 PNG files in: ", normalizePath(output_dir))
