#!/usr/bin/env Rscript

# Post-processes the repeated-MLE diagnostic only.  It never calls or changes
# the production estimator.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1L]] else stop("Usage: Rscript ... <diagnostic-directory>")
source("utils.R")

fibonacci_sphere <- function(n) {
  i <- seq_len(n) - .5; z <- 1 - 2 * i / n; phi <- pi * (3 - sqrt(5)) * i
  cbind(sqrt(pmax(0, 1 - z^2)) * cos(phi), sqrt(pmax(0, 1 - z^2)) * sin(phi), z)
}
grid <- fibonacci_sphere(512L)
sqrt_density <- function(z) sqrt((4 * pi / nrow(grid)) / 2 * d_sph_uniform_beta_mixture_s2(
  grid, c(z$mu1, z$mu2, z$mu3), z$weight_uniform, z$alpha, z$beta, log = FALSE))

summaries <- list()
for (label in c("short", "long")) {
  fits <- utils::read.csv(file.path(input, sprintf("%s_all_fits.csv", label)), check.names = FALSE)
  reference <- utils::read.csv(file.path(input, sprintf("%s_direct_production_reproduction.csv", label)))
  ok <- is.finite(fits$loglik) & (is.na(fits$error) | !nzchar(fits$error))
  near <- pmin(fits$nearest_mu, fits$nearest_antimu) <= 1e-5
  interior <- fits[ok & !fits$clipping_active, , drop = FALSE]
  S <- t(vapply(seq_len(nrow(interior)), function(i) sqrt_density(interior[i, ]), numeric(nrow(grid))))
  refS <- sqrt_density(reference[1L, ])
  Href <- sqrt(rowSums((S - rep(refS, each = nrow(S)))^2))
  H <- matrix(sqrt(pmax(0, outer(rowSums(S^2), rowSums(S^2), "+") - 2 * tcrossprod(S))),
    nrow = nrow(interior), ncol = nrow(interior))

  # Complete linkage prevents the chaining that a connected-component rule can
  # introduce.  Every pair within a reported cluster is within 0.025.
  cluster0 <- stats::cutree(stats::hclust(stats::as.dist(H), method = "complete"), h = .025)
  ord <- order(-tabulate(cluster0), -vapply(split(interior$loglik, cluster0), max, numeric(1)))
  cluster <- match(cluster0, ord)
  interior$density_basin <- cluster
  interior$hellinger_to_reported <- Href
  interior$near_endpoint_1e5 <- pmin(interior$nearest_mu, interior$nearest_antimu) <= 1e-5

  # Orient the exact UB symmetry relative to the reported direction before
  # describing parameter variability; density clustering itself needs no such convention.
  dotref <- as.matrix(interior[, c("mu1", "mu2", "mu3")]) %*% as.numeric(reference[1L, c("mu1", "mu2", "mu3")])
  flip <- as.numeric(dotref) < 0
  interior$mu1_canonical <- ifelse(flip, -interior$mu1, interior$mu1)
  interior$mu2_canonical <- ifelse(flip, -interior$mu2, interior$mu2)
  interior$mu3_canonical <- ifelse(flip, -interior$mu3, interior$mu3)
  interior$alpha_canonical <- ifelse(flip, interior$beta, interior$alpha)
  interior$beta_canonical <- ifelse(flip, interior$alpha, interior$beta)
  utils::write.csv(interior, file.path(input, sprintf("%s_complete_linkage_assignments.csv", label)), row.names = FALSE)

  by_basin <- split(interior, interior$density_basin)
  basin <- do.call(rbind, lapply(by_basin, function(z) data.frame(
    basin = z$density_basin[[1L]], n = nrow(z), fraction_all_starts = nrow(z) / nrow(fits),
    loglik_min = min(z$loglik), loglik_median = median(z$loglik), loglik_max = max(z$loglik),
    H_to_reported_min = min(z$hellinger_to_reported), H_to_reported_median = median(z$hellinger_to_reported), H_to_reported_max = max(z$hellinger_to_reported),
    near_endpoint_fraction = mean(z$near_endpoint_1e5),
    mu1_median = median(z$mu1_canonical), mu2_median = median(z$mu2_canonical), mu3_median = median(z$mu3_canonical),
    weight_median = median(z$weight_uniform), alpha_median = median(z$alpha_canonical), beta_median = median(z$beta_canonical),
    weight_sd = sd(z$weight_uniform), alpha_sd = sd(z$alpha_canonical), beta_sd = sd(z$beta_canonical)
  )))
  basin <- basin[order(-basin$n, -basin$loglik_max), ]
  utils::write.csv(basin, file.path(input, sprintf("%s_complete_linkage_density_basins.csv", label)), row.names = FALSE)

  relevant <- basin$basin[basin$n >= 5L]
  representatives <- do.call(rbind, lapply(relevant, function(k) {
    ii <- which(interior$density_basin == k)
    medoid <- ii[which.min(rowMeans(H[ii, ii, drop = FALSE]))]
    best <- ii[which.max(interior$loglik[ii])]
    rbind(data.frame(role = "medoid", interior[medoid, c("density_basin", "start_id", "mu1_canonical", "mu2_canonical", "mu3_canonical", "weight_uniform", "alpha_canonical", "beta_canonical", "loglik", "AIC", "nearest_mu", "nearest_antimu", "near_endpoint_1e5", "hellinger_to_reported")]),
          data.frame(role = "best_loglik", interior[best, c("density_basin", "start_id", "mu1_canonical", "mu2_canonical", "mu3_canonical", "weight_uniform", "alpha_canonical", "beta_canonical", "loglik", "AIC", "nearest_mu", "nearest_antimu", "near_endpoint_1e5", "hellinger_to_reported")]))
  }))
  utils::write.csv(representatives, file.path(input, sprintf("%s_representative_density_basins.csv", label)), row.names = FALSE)
  summaries[[label]] <- data.frame(dataset = label, starts = nrow(fits), convergence_zero = sum(fits$convergence == 0L, na.rm = TRUE),
    clipping_active = sum(fits$clipping_active, na.rm = TRUE), near_endpoint_1e5 = sum(near, na.rm = TRUE),
    nonclipped = nrow(interior), reported_like_H_le_0025 = sum(Href <= .025), reported_start_basin = interior$density_basin[match(1L, interior$start_id)])
}
utils::write.csv(do.call(rbind, summaries), file.path(input, "repeated_mle_overall_summary.csv"), row.names = FALSE)
writeLines(c(
  "Read this directory in order:",
  "1. repeated_mle_overall_summary.csv: convergence, clipping, and close-to-data rates.",
  "2. short_complete_linkage_density_basins.csv and long_complete_linkage_density_basins.csv: stable fitted-density basins.",
  "3. *_representative_density_basins.csv: medoid and highest-likelihood fit for each basin with at least five starts.",
  "4. *_complete_linkage_assignments.csv: every nonclipped fit, its basin, Hellinger distance to the reported fit, and canonicalized parameters.",
  "The reported fit is the start_id=1 row.  'Near endpoint' means min geodesic distance to a data point or its antipode <= 1e-5 radians; it is intentionally broader than the implementation's literal clipping flag.",
  "Clusters are complete-linkage Hellinger clusters (threshold .025) on a 512-point equal-area sphere grid, so no cluster is created solely from parameter distance."
), file.path(input, "READING_GUIDE.txt"))
message("Wrote post-processing to: ", normalizePath(input))
