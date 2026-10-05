args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript --vanilla scripts/update_paper_fig1_radius.R <positive-radius> [nonnegative-amplitude-factor]"
if (!(length(args) %in% 1:2)) stop(usage)
radius <- suppressWarnings(as.numeric(args[1]))
factor <- if (length(args) == 2L) suppressWarnings(as.numeric(args[2])) else 1
if (anyNA(c(radius, factor)) || !is.finite(radius) || !is.finite(factor) ||
    radius <= 0 || factor < 0) {
  stop(usage)
}

source("scripts/gaussian_process_s1_visualization.R")
saved <- readRDS("output/convergence/gaussian_process/s1/limit_gaussian_s1_exact_kappa2_angles16_t200_tdomain_result.rds")
if (radius >= 1 + saved$ray_extension) {
  stop("The circle radius must be smaller than the fixed outer ray radius (", 1 + saved$ray_extension, ").")
}
if (factor != 1) {
  max_g <- max(abs(saved$curve_data$g))
  extra_scale <- if (max_g > 0) (factor - 1) * saved$lateral_scale / max_g else 0
  saved$curve_data$x_circle <- saved$curve_data$x_circle -
    extra_scale * saved$curve_data$g * sin(saved$curve_data$theta)
  saved$curve_data$y_circle <- saved$curve_data$y_circle +
    extra_scale * saved$curve_data$g * cos(saved$curve_data$theta)
}

paper_img <- Sys.getenv("PAPER_IMG_DIR")
if (!nzchar(paper_img)) stop("Set PAPER_IMG_DIR to the manuscript image directory.")
git_active <- system2("git", c("-C", shQuote(paper_img), "rev-parse", "--is-inside-work-tree"), stdout = TRUE, stderr = FALSE)
if (!identical(git_active, "true")) stop("The paper image directory must be inside an active Git repository.")
image <- file.path(paper_img, "limit_gaussian_s1_exact_kappa2_angles16_t200_rainbow_left.png")
if (!file.exists(image)) stop("Figure 1 image not found: ", image)

panel <- replot_limit_gaussian_s1_vmf_from_result(
  saved, color_scheme = "rainbow", circle_radius = radius,
  ray_linewidth = 0.3, ray_linetype = "dashed"
)$left_panel
ggplot2::ggsave(image, panel, width = 6, height = 6, dpi = 300)
cat("Updated ", image, " with circle radius ", radius,
    " and amplitude factor ", factor, "\n", sep = "")
