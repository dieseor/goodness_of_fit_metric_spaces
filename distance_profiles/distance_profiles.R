# Distance Profile Analysis Framework
# Functions for validating theoretical distance profiles against empirical data

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sphunif)
  library(MASS)
  library(mvtnorm)
  library(viridis)
  library(gridExtra)
})

#' Euclidean distance between two points
#' @param x1 First point (vector)
#' @param x2 Second point (vector)
#' @return Euclidean distance
euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}

#' Compute distances from omega to all points in a dataset
#' @param omega Reference point (vector)
#' @param data Matrix with observations in rows
#' @return Vector of distances
compute_distances <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  apply(data, 1, function(x) euclidean_distance(omega, x))
}
