#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(compositions)
  library(ggplot2)
})

output_file <- "/Users/Diego/Library/CloudStorage/Dropbox/Apps/Overleaf/Goodness-of-fit for distributions on metric spaces/AoS/img/clam_east_west_ilr.png"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

data(ClamEast, package = "compositions")
data(ClamWest, package = "compositions")

x_east <- as.matrix(ClamEast)
x_west <- as.matrix(ClamWest)
x_all <- rbind(x_east, x_west)
group <- c(rep("ClamEast", nrow(x_east)), rep("ClamWest", nrow(x_west)))

z_all <- as.matrix(ilr(acomp(x_all)))
if (ncol(z_all) < 2L) {
  stop("The clam datasets must have ilr dimension at least 2.")
}

plot_df <- data.frame(
  ilr1 = z_all[, 1],
  ilr2 = z_all[, 2],
  group = factor(group, levels = c("ClamEast", "ClamWest"))
)

plt <- ggplot(plot_df, aes(x = ilr1, y = ilr2, color = group)) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_color_manual(values = c("ClamEast" = "#2166ac", "ClamWest" = "#b2182b")) +
  labs(
    x = "First ilr coordinate",
    y = "Second ilr coordinate",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )

ggsave(
  filename = output_file,
  plot = plt,
  width = 5.6,
  height = 4.3,
  dpi = 200,
  bg = "transparent"
)

cat(output_file, "\n")
