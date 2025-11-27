source("R/utils.R")
source("R/gaussian_process_vmf.R")
library(rotasym)

MU <- c(0,0,1)
KAPPA <- 2.0
omega_grid <- generate_canonical_lattice(10)
t_grid <- seq(0, 2, length.out=10)

cat("F_theoretical for omega[1,]:\n")
probs <- compute_distance_profile_vmf(omega_grid[1,], MU, KAPPA, t_grid, "chordal")
print(round(probs, 4))
cat("\nLast value should be ~1, not 0!\n")
