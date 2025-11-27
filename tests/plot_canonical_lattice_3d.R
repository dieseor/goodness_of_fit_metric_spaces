# Interactive 3D visualization of canonical lattice on sphere using rgl
# This will work in RStudio and can be rotated with mouse!

library(rotasym)
library(rgl)

source("R/utils.R")
source("R/gaussian_process_vmf.R")

cat("=== INTERACTIVE 3D VISUALIZATION OF CANONICAL LATTICE ===\n\n")

# Generate lattice
n_points <- 30
cat("Generating canonical lattice with", n_points, "points...\n")
omega_grid <- generate_canonical_lattice(n_points)

# Check points are on sphere
norms <- sqrt(rowSums(omega_grid^2))
cat("Norms range: [", min(norms), ",", max(norms), "]\n")
if (max(abs(norms - 1)) > 1e-10) {
  stop("ERROR: Points are not on unit sphere!")
} else {
  cat("✓ All points are on unit sphere\n\n")
}

# Create 3D plot
cat("Creating interactive 3D plot...\n")
cat("In RStudio: The plot will appear in the 'Plots' pane or a new window\n")
cat("Use mouse to rotate, zoom, and pan the view!\n\n")

# Open 3D device
open3d(windowRect = c(50, 50, 850, 850))

# Set background
bg3d("white")

# Draw the sphere surface (semi-transparent)
theta <- seq(0, pi, length.out = 50)
phi <- seq(0, 2*pi, length.out = 50)

# Create mesh for sphere
x_sphere <- outer(sin(theta), cos(phi))
y_sphere <- outer(sin(theta), sin(phi))
z_sphere <- outer(cos(theta), rep(1, length(phi)))

# Plot sphere surface
persp3d(x_sphere, y_sphere, z_sphere, 
        col = "lightblue", alpha = 0.2, 
        lit = TRUE, smooth = TRUE,
        add = TRUE)

# Add coordinate axes
axes3d(c('x', 'y', 'z'), color = "gray")
title3d(main = paste("Canonical Lattice on S² (n =", n_points, ")"), 
        xlab = "X", ylab = "Y", zlab = "Z",
        cex = 1.5)

# Plot lattice points
# Color by Z coordinate (blue = -1, red = +1)
z_values <- omega_grid[, 3]
colors <- colorRampPalette(c("blue", "red"))(100)
color_indices <- cut(z_values, breaks = 100, labels = FALSE)
point_colors <- colors[color_indices]

# Add points
spheres3d(omega_grid[, 1], omega_grid[, 2], omega_grid[, 3],
          radius = 0.03, color = point_colors)

# Add a few reference circles for orientation
# Equator
circle_theta <- seq(0, 2*pi, length.out = 100)
lines3d(cos(circle_theta), sin(circle_theta), rep(0, 100),
        col = "gray", lwd = 2)

# Prime meridian
lines3d(rep(0, 100), cos(seq(0, 2*pi, length.out = 100)), 
        sin(seq(0, 2*pi, length.out = 100)),
        col = "gray", lwd = 1)

# Another meridian (perpendicular)
lines3d(cos(seq(0, 2*pi, length.out = 100)), rep(0, 100),
        sin(seq(0, 2*pi, length.out = 100)),
        col = "gray", lwd = 1)

# Set nice viewing angle
view3d(theta = 45, phi = 30, zoom = 0.8)

cat("✓ 3D plot created!\n\n")
cat("Controls:\n")
cat("  - Click and drag to rotate\n")
cat("  - Scroll to zoom\n")
cat("  - Right-click and drag to pan\n")
cat("  - Press ESC or close window to exit\n\n")

# Save snapshot
snapshot_file <- file.path("output", "canonical_lattice_3d_snapshot.png")
rgl.snapshot(snapshot_file)
cat("Snapshot saved to:", snapshot_file, "\n\n")

# Also save as interactive WebGL (can open in browser!)
webgl_file <- file.path("output", "canonical_lattice_3d_interactive.html")

# Create the WebGL export
cat("Exporting interactive WebGL version...\n")
tryCatch({
  htmlwidgets::saveWidget(
    rglwidget(width = 800, height = 800),
    file = normalizePath(webgl_file),
    selfcontained = TRUE
  )
  cat("✓ Interactive HTML saved to:", webgl_file, "\n")
  cat("  Open this file in a web browser to interact with the 3D plot!\n\n")
}, error = function(e) {
  cat("Note: Could not save interactive HTML:", e$message, "\n")
  cat("The rgl window is still open for interaction.\n\n")
})

# Print statistics
cat("=== LATTICE STATISTICS ===\n")
cat("Total points:", nrow(omega_grid), "\n")
cat("Upper hemisphere (z > 0):", sum(omega_grid[, 3] > 0), "\n")
cat("Lower hemisphere (z < 0):", sum(omega_grid[, 3] < 0), "\n\n")

# Compute pairwise distances
n <- nrow(omega_grid)
chordal_dists <- numeric(n * (n - 1) / 2)
idx <- 1
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    chordal_dists[idx] <- sqrt(sum((omega_grid[i, ] - omega_grid[j, ])^2))
    idx <- idx + 1
  }
}

cat("Pairwise chordal distances:\n")
cat("  Min:", round(min(chordal_dists), 4), "\n")
cat("  Max:", round(max(chordal_dists), 4), "\n")
cat("  Mean:", round(mean(chordal_dists), 4), "\n")
cat("  Median:", round(median(chordal_dists), 4), "\n\n")

cat("=== VISUALIZATION COMPLETE ===\n")
cat("\nThe 3D plot window should be open for interaction.\n")
cat("Press ENTER to close and exit, or manually close the window...\n")

# In RStudio, keep the plot open
if (interactive()) {
  cat("\n(Plot will stay open in RStudio - close window when done)\n")
} else {
  readline()
  rgl.close()
}
