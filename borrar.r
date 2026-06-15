## Plot the three-part compositional datasets used in the logistic Gaussian
## real-data screening on the simplex.
##
## D = 3 datasets: ArcticLake, Metabolites, Sediments, SkyeAFM.

suppressPackageStartupMessages({
  library(compositions)
  library(ggplot2)
})

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------

as_composition_df <- function(x, parts = NULL, dataset_name) {
  x <- as.data.frame(x)
  if (is.null(parts)) {
    numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
    parts <- numeric_cols[seq_len(3)]
  }
  missing_parts <- setdiff(parts, names(x))
  if (length(missing_parts) > 0) {
    stop(
      "Dataset '", dataset_name, "' is missing columns: ",
      paste(missing_parts, collapse = ", "),
      call. = FALSE
    )
  }

  comp <- x[, parts, drop = FALSE]
  comp <- comp[stats::complete.cases(comp), , drop = FALSE]
  comp <- as.data.frame(apply(comp, 2, as.numeric))
  names(comp) <- c("x1", "x2", "x3")

  row_sums <- rowSums(comp)
  comp <- comp / row_sums

  ## Barycentric coordinates in an equilateral triangle:
  ## vertex 1 = (0,0), vertex 2 = (1,0), vertex 3 = (1/2, sqrt(3)/2).
  data.frame(
    dataset = dataset_name,
    x = comp$x2 + 0.5 * comp$x3,
    y = (sqrt(3) / 2) * comp$x3,
    part1 = comp$x1,
    part2 = comp$x2,
    part3 = comp$x3
  )
}

## -----------------------------------------------------------------------------
## Load datasets from the compositions package
## -----------------------------------------------------------------------------

data("ArcticLake", package = "compositions")
data("Metabolites", package = "compositions")
data("Sediments", package = "compositions")
data("SkyeAFM", package = "compositions")

plot_data <- rbind(
  as_composition_df(ArcticLake, c("sand", "silt", "clay"), "ArcticLake"),
  as_composition_df(Metabolites, dataset_name = "Metabolites"),
  as_composition_df(Sediments, c("sand", "silt", "clay"), "Sediments"),
  as_composition_df(SkyeAFM, c("A", "F", "M"), "SkyeAFM")
)

simplex_boundary <- data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, sqrt(3) / 2, 0)
)

## -----------------------------------------------------------------------------
## Plot
## -----------------------------------------------------------------------------

p <- ggplot() +
  geom_path(
    data = simplex_boundary,
    aes(x = x, y = y),
    linewidth = 0.6
  ) +
  geom_point(
    data = plot_data,
    aes(x = x, y = y, colour = dataset),
    size = 2.2,
    alpha = 0.85
  ) +
  coord_equal(xlim = c(-0.04, 1.04), ylim = c(-0.03, sqrt(3) / 2 + 0.04), expand = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    colour = "Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.93, 0.91),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA),
    legend.margin = margin(2, 2, 2, 2)
  )

print(p)

ggsave(
  filename = "simplex_D3_datasets.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)