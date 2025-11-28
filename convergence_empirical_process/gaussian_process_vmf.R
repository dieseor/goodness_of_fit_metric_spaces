# ============================================================================
# GAUSSIAN PROCESS FUNCTIONS - VON MISES-FISHER DISTRIBUTION
# ============================================================================
# Functions specific to the von Mises-Fisher distribution on the unit sphere
# for Gaussian process limit theory in goodness-of-fit testing.
# ============================================================================

# Load generic utilities
source(file.path("utils.R"))


# Load required libraries
library(rotasym)  # For vMF functions (c_vMF, r_vMF)
library(parallel)



# ============================================================================
# DISTANCE PROFILE AND PROBABILITIES
# ============================================================================

#' Compute one row of the covariance matrix (VECTORIZED over all columns)
#' @param idx Row index in flattened grid
#' @param omega_grid Matrix of omega points (n_omega × q)
#' @param t_grid Vector of t values (length n_t)
#' @param mu Mean direction
#' @param kappa Concentration parameter
#' @param distance_type Distance type
#' @param mc_samples Monte Carlo samples (n_mc × q)
#' @param A_q_kappa Pre-computed A_q(κ)
#' @param var_X Pre-computed Var(X)
#' @return Vector of covariances for the entire row (length n_omega × n_t)
row_cov_vmf <- function(idx, omega_grid, t_grid, mu, kappa, distance_type = "chordal", mc_samples, A_q_kappa, var_X, precomp = NULL, col_idxs = NULL, h0 = c("simple","composite"), unknown_param = NULL) {
  # Prefer explicit passing of h0/unknown_param as arguments; fall back to defaults
  h0 <- match.arg(h0)
  if (!is.null(unknown_param)) unknown_param <- match.arg(unknown_param, choices = c('xi'))
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  q <- length(mu) - 1  # Dimension of sphere S^q in R^{q+1}
  n_mc <- nrow(mc_samples)
  
  # Decode row index to (omega1, t1)
  omega1_idx <- ((idx - 1) %% n_omega) + 1
  t1_idx <- ((idx - 1) %/% n_omega) + 1
  omega1 <- omega_grid[omega1_idx, ]
  t1 <- t_grid[t1_idx]
  
  # Pre-compute for omega1, t1 (scalars/vectors)
  F1 <- theoretical_distance_profile_vmf(omega1, mu, kappa, t1, distance_type)
  E1 <- compute_conditional_expectation_vmf(omega1, t1, mu, kappa, distance_type, mc_samples)
  m1 <- E1 - A_q_kappa * mu
  # inverse(Var(X)) for vectorized quadratic computations
  inv_var_X <- tryCatch(solve(var_X), error = function(e) {
    warning('Var(X) is singular or near-singular; falling back to MASS::ginv pseudo-inverse')
    return(MASS::ginv(var_X))
  })
  
  # =========================================================================
  # VECTORIZED: Pre-compute distances from mc_samples to omega1
  # =========================================================================
  dots_1 <- mc_samples %*% omega1  # (n_mc × 1)
  if (distance_type == "chordal") {
    dists_1 <- sqrt(2 * (1 - dots_1))
  } else {
    # Ensure numerical stability before calling acos
    dots_1 <- check_dot_products(dots_1)
    dists_1 <- acos(dots_1)
  }
  in_ball_1 <- dists_1 <= t1  # (n_mc × 1) logical
  
  # =========================================================================
  # VECTORIZED: Compute distances between omega1 and ALL omega2
  # =========================================================================
  dot_prods_omega <- omega_grid %*% omega1  # (n_omega × 1)
  if (distance_type == "chordal") {
    d_omega_vec <- sqrt(2 * (1 - dot_prods_omega))  # (n_omega × 1)
  } else {
    dot_prods_omega <- check_dot_products(dot_prods_omega)
    d_omega_vec <- acos(dot_prods_omega)
  }
  
  # =========================================================================
  # VECTORIZED: Compute F2, E2, m2 for ALL (omega2, t2) pairs
  # =========================================================================
  # Use precomputed F2/E2/dists if provided, else compute them using centralized helper
  if (!is.null(precomp)) {
    F2_matrix <- precomp$F2_matrix
    E2_array <- precomp$E2_array
    E2_mat <- precomp$E2_mat
    F2_vec <- precomp$F2_vec
    m2_mat <- precomp$m2_mat
    dists_all <- precomp$dists_all
    in_ball_list <- precomp$in_ball_list
    n_mc <- nrow(dists_all)
  } else {
    precomp_loc <- compute_precomp_vmf(mc_samples, omega_grid, t_grid, mu, kappa, distance_type, A_q_kappa)
    F2_matrix <- precomp_loc$F2_matrix
    E2_array <- precomp_loc$E2_array
    E2_mat <- precomp_loc$E2_mat
    F2_vec <- precomp_loc$F2_vec
    m2_mat <- precomp_loc$m2_mat
    dists_all <- precomp_loc$dists_all
    in_ball_list <- precomp_loc$in_ball_list
    n_mc <- nrow(dists_all)
  }
  
  # Note: E2_array/E2_mat have been computed above or via precomp
  
  # =========================================================================
  # VECTORIZED: Compute joint probabilities for ALL (omega2, t2) pairs
  # =========================================================================
  # Vectorized computation of joint probabilities for all (omega2, t2) pairs
  # Create logical matrix: (n_omega × n_t), TRUE if omega2 can overlap for (t1, t2)
  t1_vec <- rep(t1, n_t)
  t2_vec <- t_grid
  # For each t2, get can_overlap vector (n_omega × n_t)
  can_overlap_mat <- outer(d_omega_vec, t1_vec + t2_vec + 1e-10, "<=")
  # For each omega2 and t2, compute joint indicator and mean
  # Precompute in_ball_2 matrix: (n_mc × n_omega × n_t)
  # dists_all: (n_mc × n_omega), t_grid: (n_t)
  # For each omega2, create (n_mc × n_t) logical matrix
  # Fully vectorized: in_ball_2_array (n_mc × n_omega × n_t)
  # Create in_ball_2_array: logical array (n_mc x n_omega x n_t)
  # Compute P_joint_matrix by iterating over t (n_t usually small) and computing counts
  P_joint_matrix <- matrix(0, nrow = n_omega, ncol = n_t)
  for (k in seq_len(n_t)) {
    in_ball_t_k <- in_ball_list[[k]]
    # in_ball_t_k has dimensions (n_mc x n_omega); compute counts via matrix multiply
    # joint_counts = t(in_ball_t_k) %*% in_ball_1  (vectorized/matrix method)
    joint_counts <- as.numeric(t(in_ball_t_k) %*% as.numeric(in_ball_1))
    P_joint_matrix[, k] <- joint_counts / n_mc
  }
  P_joint_matrix[!can_overlap_mat] <- 0
  
  # =========================================================================
  # VECTORIZED: Compute covariance for ALL (omega2, t2) pairs
  # =========================================================================
  # Fully vectorized: Compute covariance for ALL (omega2, t2) pairs
  # Prepare m2 matrix (n_total × q+1)
  # ...
  # Use the precomputed E2_mat/F2_vec from earlier (or from precomp)
  F2_vec <- as.vector(F2_matrix)  # flatten with omega as fastest index
  P_joint_vec <- as.vector(P_joint_matrix)
  
  # m2: center E2
  m2_mat <- sweep(E2_mat, 2, A_q_kappa * mu, "-")
  
  # quadratic form for each row (vectorized)
  # m1' (inv_var_X) is a (1 x q+1) row vector; multiplying by each m2 gives a scalar
  L_vec <- as.numeric(t(m1) %*% inv_var_X) # length q+1
  # If computing only a subset of columns (e.g., upper-triangle), restrict the expensive multiplier
  cols_to_consider <- if (is.null(col_idxs)) seq_len(n_total) else col_idxs
  F2_vec_sub <- F2_vec[cols_to_consider]
  m2_mat_sub <- m2_mat[cols_to_consider, , drop = FALSE]
  P_joint_vec_sub <- P_joint_vec[cols_to_consider]
  # Compute quadratic forms deterministically using matrix multiplication
  # Use matrix multiply: (m2 %*% L_vec)
  quadratic_forms_sub <- as.numeric(m2_mat_sub %*% L_vec)
  # Determine whether to include composite correction based on h0
  C_base_sub <- P_joint_vec_sub - F1 * F2_vec_sub
  if (!is.null(h0) && h0 == 'composite') {
    composite_corr_sub <- - F1 * F2_vec_sub * quadratic_forms_sub
    cov_row_sub <- C_base_sub + composite_corr_sub
  } else {
    # simple null: return base covariance (joint - F1*F2)
    cov_row_sub <- C_base_sub
  }

  cov_row <- rep(NA_real_, n_total)
  cov_row[cols_to_consider] <- cov_row_sub
  
  return(cov_row)
}

 
cov_vmf <- function(omega_grid, t_grid, mu, kappa, distance_type = "chordal", n_mc_samples = 1000, n_cores = 10, seed = NULL, upper_triangle = FALSE, mc_samples = NULL, h0 = c("simple","composite"), unknown_param = NULL) {
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  q <- length(mu) - 1  # dimension of sphere S^q in R^{q+1}
  
  # Normalize h0 and unknown_param defaults
  h0 <- match.arg(h0)
  cat("Creating vMF covariance matrix of size", n_total, "x", n_total, "\n")
  cat("Using", n_cores, "cores\n")
  
  # Pre-compute A_q(κ)
  A_q_kappa <- besselI(kappa, nu = (q + 1) / 2, expon.scaled = TRUE) / 
               besselI(kappa, nu = (q - 1) / 2, expon.scaled = TRUE)
  
  
  # Generate Monte Carlo samples from vMF(μ, κ)
  if (is.null(mc_samples)) {
    if (!is.null(seed)) set.seed(seed)
    cat("Generating", n_mc_samples, "samples from vMF(μ, κ)...\n")
    t_mc_start <- Sys.time()
    mc_samples <- rotasym::r_vMF(n = n_mc_samples, mu = mu, kappa = kappa)
    t_mc_end <- Sys.time()
    cat("Samples generated in", round(as.numeric(difftime(t_mc_end, t_mc_start, units = "secs")), 2), "seconds.\n")
  } else {
    # Use supplied MC samples (assumes they have correct dimensions)
    n_mc_samples <- nrow(mc_samples)
    cat("Using provided", n_mc_samples, "samples from vMF(μ, κ)...\n")
  }
  
  # Compute Var(X) using theoretical formula
  # Var(X) = (A_q/κ)I + [1 - A_q² - (q+1)A_q/κ]μμ'
  scalar_coef <- 1 - A_q_kappa^2 - ((q + 1) * A_q_kappa / kappa)
  var_X <- (A_q_kappa / kappa) * diag(q + 1) + scalar_coef * outer(mu, mu)
  
  # Setup parallel cluster
  # Force single-threaded BLAS/OpenMP to reduce non-deterministic numeric operations across processes
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)  # Ensure cluster is stopped even if interrupted
  # Set seed for workers if provided, else default to 42
  clusterSetRNGStream(cl, iseed = ifelse(is.null(seed), 42, seed))
  
  # Export functions and variables to workers
  clusterEvalQ(cl, {
    library(rotasym)
    source(file.path("utils.R"))
    # Do not source vmf_cov_functions.R on workers to avoid function definition collisions
    # Exporting the master functions via clusterExport below ensures workers execute
  })
  # Reduce BLAS/OpenMP threads on workers for deterministic arithmetic if possible
  clusterEvalQ(cl, {
    try({
      
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
      }
    }, silent = TRUE)
  })

  # Ensure omega_grid is a matrix and t_grid is a vector before exporting
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  # Index mapping for flattened grid columns
  omega_idx_vec <- ((0:(n_total - 1)) %% n_omega) + 1
  t_idx_vec <- ((0:(n_total - 1)) %/% n_omega) + 1

  # --- Precompute shared intermediates for vectorized rows (avoid per-row recomputation) ---
  cat("Precomputing intermediates for vMF covariance (shared across rows)...\n")
  precomp_loc <- compute_precomp_vmf(mc_samples, omega_grid, t_grid, mu, kappa, distance_type, A_q_kappa)
  # add index vectors used for flattening to the precomp list
  precomp <- c(precomp_loc, list(omega_idx_vec = omega_idx_vec, t_idx_vec = t_idx_vec))

  # Ensure the h0 value is exported as a plain variable to workers (avoid quoting or name shadowing)
  h0_val <- h0

  # Export all needed functions and variables from the main environment
  clusterExport(cl, c(
    "check_dot_products", "sphere_distance", 
    "compute_conditional_expectation_vmf",
    "row_cov_vmf", "compute_precomp_vmf",
    "omega_grid", "t_grid", "mu", "kappa", "distance_type",
    "mc_samples", "A_q_kappa", "var_X", "precomp", "omega_idx_vec", "t_idx_vec",
    "h0", "h0_val", "unknown_param"
  ), envir = environment())
  # Export control flags
  clusterExport(cl, c("upper_triangle", "n_total"), envir = environment())
  
  start_time <- Sys.time()

  
  # Parallel computation: create chunked tasks for dynamic load balancing and reduced RPC overhead
  # Choose a chunk_size so each core has multiple chunks for better load balancing
  default_chunks_per_core <- 4
  chunk_size <- max(1, floor(n_total / (n_cores * default_chunks_per_core)))
  # Split into contiguous chunks for better cache locality and fewer RPCs
  row_indices <- seq_len(n_total)
  row_chunks <- split(row_indices, ceiling(seq_along(row_indices) / chunk_size))
  
  cat("Computing", sum(sapply(row_chunks, length)), "rows in parallel\n")
  
  # Compute rows in parallel
  results <- parLapply(cl, row_chunks, function(row_indices) {
    local_rows <- list()
    for (i in row_indices) {
      # Execute the row computation on the worker; errors will propagate

        row_vec <- row_cov_vmf(i, omega_grid, t_grid, mu, kappa,
                    distance_type, mc_samples,
                    A_q_kappa, var_X, precomp = precomp,
                  col_idxs = if (upper_triangle) seq(i, n_total) else NULL,
                   h0 = h0_val, unknown_param = unknown_param)
        diag_output <- NULL
      local_rows[[length(local_rows) + 1]] <- list(i = i, row = row_vec, diag = diag_output)
    }
    return(local_rows)
  })


  # Print diagnostics and errors from all workers (only non-empty diagnostics)
  for (worker_results in results) {
    for (result in worker_results) {
      if (!is.null(result$diag) && length(result$diag) > 0 && any(nzchar(result$diag))) {
        cat(paste(result$diag, collapse = "\n"), "\n")
      }
      if (!is.null(result$error)) cat(sprintf("[ERROR] idx=%d: %s\n", result$i, result$error))
    }
  }
  
  # Assemble matrix
  cov_matrix <- matrix(0, n_total, n_total)
  for (worker_results in results) {
    for (result in worker_results) {
      if (!is.null(result$row) && length(result$row) == n_total) {
        non_na_cols <- which(!is.na(result$row))
        if (length(non_na_cols) > 0) {
          # Only assign upper triangle part (i <= j) so each symmetric pair is written exactly once
          upper_cols <- non_na_cols[non_na_cols >= result$i]
          if (length(upper_cols) > 0) {
            cov_matrix[result$i, upper_cols] <- result$row[upper_cols]
            # Mirror for symmetry: set symmetric entries exactly to the same values
            cov_matrix[upper_cols, result$i] <- result$row[upper_cols]
          }
        }
      }
    }
  }

  # No automatic PSD enforcement: we do not alter covariance matrices here.

  end_time <- Sys.time()
  time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("Covariance matrix created in", round(time_elapsed, 1), "seconds!\n\n")
  
  # Check positive definiteness (for diagnostic purposes only)
  eigs <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eigs)
  max_eig <- max(eigs)
  n_negative <- sum(eigs < -1e-10)
  
  cat("Covariance matrix diagnostics:\n")
  cat("  Eigenvalue range: [", format(min_eig, scientific = TRUE), ", ", 
      format(max_eig, scientific = TRUE), "]\n", sep = "")
  cat("  Negative eigenvalues:", n_negative, "\n")
  
  if (n_negative > 0) {
    # Matrix has negative eigenvalues; do not auto-correct the covariance matrix here
    cat("  WARNING: Matrix has", n_negative, "negative eigenvalues!\n")
    cat("  Note: No automatic covariance correction is applied; returning matrix as-is.\n")
  } else {
    cat("  ✓ Matrix is positive semi-definite\n")
  }
  cat("\n")
  # Parity checks removed.
  
  return(cov_matrix)
}

# ============================================================================
# PROCESS SIMULATION
# ============================================================================

#' Simulate the limiting Gaussian process for vMF
#' @param omega_grid Matrix of omega points (m × q)
#' @param t_grid Vector of t values
#' @param mu Mean direction of vMF
#' @param kappa Concentration parameter
#' @param distance_type Either "chordal" or "geodesic"
#' @param M Number of Monte Carlo simulations
#' @param n_mc_samples Number of MC samples for covariance computation
#' @param n_cores Number of cores
#' @return Vector of M supremum values from the Gaussian process
simulate_limit_gaussian_vmf <- function(omega_grid, t_grid, mu, kappa,
                                       distance_type = "chordal",
                                       M = 10000,
                                       n_mc_samples = 100000,
                                       n_cores = 10,
                                       seed = NULL,
                                       h0 = c("simple","composite"),
                                       unknown_param = NULL) {
  h0 <- match.arg(h0)
  cat("=== Simulating Gaussian Process Limit for vMF Distribution ===\n")
  
  # Create covariance matrix; for composite null, create_covariance_matrix_vmf
  # will accept the h0/unknown_param arguments.
  cov_matrix <- cov_vmf(omega_grid, t_grid, mu, kappa,
                                            distance_type, n_mc_samples, n_cores, seed = seed, h0 = h0, unknown_param = unknown_param)
  
  # Use generic function to simulate
  supremum_values <- simulate_limit_gaussian(cov_matrix, M)
  
  return(supremum_values)
}

#' Simulate the empirical process for vMF (parallelized)
#' @param omega_grid Matrix of omega points
#' @param t_grid Vector of t values
#' @param n Sample size
#' @param mu Mean direction of vMF
#' @param kappa Concentration parameter
#' @param distance_type Either "chordal" or "geodesic"
#' @param M Number of Monte Carlo simulations
#' @param n_cores Number of cores
#' @return Vector of M supremum values from the empirical process
simulate_empirical_process_vmf <- function(omega_grid, t_grid, n, mu, kappa,
                                          distance_type = "chordal",
                                          M = 10000,
                                          n_cores = 10,
                                          seed = NULL,
                                          h0 = c("simple","composite"),
                                          unknown_param = NULL) {
  h0 <- match.arg(h0)
  cat("=== Simulating Empirical Process for vMF Distribution ===\n")
  
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  q <- length(mu) - 1  # dimension of sphere S^q in R^{q+1}
  
  # Pre-compute theoretical values for all (omega, t) pairs if simple null
  cat("Computing theoretical distance profiles...\n")
  if (h0 == "simple") {
    F_theoretical_matrix <- t(apply(omega_grid, 1, function(omega) {
      theoretical_distance_profile_vmf(omega, mu, kappa, t_grid, distance_type)
    }))
  } else {
    # Placeholder for composite null: the master process does not compute a single global
    # theoretical matrix; workers will compute F_theoretical using the sample-estimated
    # xi (MLE) per simulation and use that in the scaled differences.
    F_theoretical_matrix <- NULL
  }
  
  # Setup parallel cluster
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)  # Ensure cluster is stopped even if interrupted
  # Allow caller to pass seed for reproducibility; fallback to 123
  clusterSetRNGStream(cl, iseed = ifelse(is.null(seed), 123, seed))
  
  # Setup workers with correct directory and sources
  clusterEvalQ(cl, source("utils.R"))
  
  # Export function to workers
  clusterEvalQ(cl, {
    library(rotasym)
    # Ensure utils and required packages are available on workers
    library(movMF)
    source(file.path('utils.R'))
  })
  
    clusterExport(cl, c("check_dot_products", "omega_grid", "t_grid", "n", "mu", "kappa"
  , "distance_type",
                      "F_theoretical_matrix", "n_omega", "n_t", "q", "h0", "unknown_param", "theoretical_distance_profile_vmf", "compute_mle_xi"),
                envir = environment())
  # Export compute_mle_xi from utils so workers can call it for composite nulls
  clusterExport(cl, c("compute_mle_xi"), envir = environment())
  
  start_time <- Sys.time()
  
  # Parallel computation
  supremum_values <- {
    # Distribute simulations across cores
    sim_chunks <- lapply(1:n_cores, function(i) seq(from = i, to = M, by = n_cores))
    
    cat("Computing", M, "simulations in parallel across", n_cores, "cores...\n")
    cat("Sample size n =", n, "\n")
    
    results <- parLapply(cl, sim_chunks, function(sim_indices) {
      local_supremums <- numeric(length(sim_indices))
      
      # Force omega_grid and t_grid to correct types
      omega_grid_matrix <- as.matrix(omega_grid)
      t_grid_vec <- as.numeric(t_grid)
      
      for (idx in seq_along(sim_indices)) {
  # Step 1: Generate sample
  sample_data <- rotasym::r_vMF(n = n, mu = mu, kappa = kappa)
        
        # Step 2: Compute distances
        dot_products <- sample_data %*% t(omega_grid_matrix)
        #dot_products <- check_dot_products(dot_products)
        
        if (distance_type == "chordal") {
          distance_matrix <- sqrt(2 * (1 - dot_products))
        } else {  # geodesic
          distance_matrix <- acos(dot_products)
        }
        
        # Step 3: Compute empirical CDF (VECTORIZED like Normal case!)
        # sapply over t_grid, each iteration computes colMeans
        F_hat_matrix <- sapply(t_grid_vec, function(t_val) {
          colMeans(distance_matrix <= t_val)
        })
        # Step 4: Compute theoretical matrix (simple vs composite handling)
        if (h0 == "simple") {
          F_theoretical_local <- F_theoretical_matrix
        } else {
          # Composite null: estimate xi (vector) from sample; we use the 
          # MLE function defined in R/utils.R.
          xi_hat <- compute_mle_xi(sample_data)
          # Convert to mu_hat and kappa_hat placeholders
          kappa_hat <- sqrt(sum(xi_hat^2))
          mu_hat <- xi_hat / kappa_hat
          F_theoretical_local <- t(apply(omega_grid, 1, function(omega) {
            theoretical_distance_profile_vmf(omega, mu_hat, kappa_hat, t_grid_vec, distance_type)
          }))
        }
        
        # Step 5: Compute supremum
        scaled_diff_matrix <- sqrt(n) * abs(F_hat_matrix - F_theoretical_local)
        local_supremums[idx] <- max(scaled_diff_matrix)
      }
      
      return(local_supremums)
    })
    
    # Combine results
    all_supremums <- numeric(M)
    for (worker_idx in 1:n_cores) {
      indices <- sim_chunks[[worker_idx]]
      all_supremums[indices] <- results[[worker_idx]]
    }
    
    all_supremums
  }
  
  end_time <- Sys.time()
  time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("Empirical process simulated in", round(time_elapsed, 1), "seconds\n\n")
  
  # cat("Supremum statistics, empirical process:\n")
  # cat("  Mean:", round(mean(supremum_values), 4), "\n")
  # cat("  Median:", round(median(supremum_values), 4), "\n")
  # cat("  Max:", round(max(supremum_values), 4), "\n")
  # cat("  Min:", round(min(supremum_values), 4), "\n\n")
  
  return(supremum_values)
}

# ============================================================================
# VISUALIZATION AND ANALYSIS
# ============================================================================

#' Visualize convergence of empirical process to limit for vMF
#' @param n_values Vector of sample sizes to compare
#' @param mu Mean direction of vMF distribution
#' @param kappa Concentration parameter
#' @param distance_type Either "chordal" or "geodesic"
#' @param omega_points Number of omega grid points (canonical lattice)
#' @param t_points Number of t grid points
#' @param M Number of Monte Carlo simulations
#' @param n_mc_samples Number of MC samples for covariance
#' @param n_cores Number of cores for parallelization
#' @return List with all simulated values and convergence plot
visualize_convergence_to_limit_vmf <- function(n_values = c(50, 100, 500),
                                              mu = c(0, 0, 1),
                                              kappa = 2.0,
                                              distance_type = "chordal",
                                              omega_points = 5,
                                              t_points = 5,
                                              omega_grid = NULL,
                                              t_grid = NULL,
                                              M = 10000,
                                              n_mc_samples = 100000,
                                              n_cores = 10,
                                              seed = NULL,
                                              density_adjust = 1.0,
                                              h0 = c("simple","composite"),
                                              unknown_param = NULL,
                                              xlim = NULL) {
  h0 <- match.arg(h0)
  
  # Use provided omega_grid/t_grid if present, otherwise create automatically
  if (is.null(omega_grid)) {
    omega_grid <- generate_canonical_lattice(omega_points, dim = length(mu))
  } else {
    omega_points <- nrow(omega_grid)
  }

    if (is.null(t_grid)) {
      t_max <- if (distance_type == "chordal") 2 - 1e-8 else pi - 1e-8
      t_grid <- seq(0 + 1e-8, t_max, length.out = t_points)
    } else {
      t_points <- length(t_grid)
    }

    n_omega <- nrow(omega_grid)
    n_t <- length(t_grid)
    n_total <- n_omega * n_t
    q <- length(mu) - 1

    # MC samples
    mc_samples <- rotasym::r_vMF(n = n_mc_samples, mu = mu, kappa = kappa)
    A_q_kappa <- besselI(kappa, nu = (q + 1) / 2, expon.scaled = TRUE) /
                 besselI(kappa, nu = (q - 1) / 2, expon.scaled = TRUE)
    scalar_coef <- 1 - A_q_kappa^2 - ((q + 1) * A_q_kappa / kappa)
    var_X <- (A_q_kappa / kappa) * diag(q + 1) + scalar_coef * outer(mu, mu)

    # Generate t grid with a small epsilon to avoid numerical issues

    # Vectorized version
    t_vectorized_start <- Sys.time()
    # Use the same MC samples for the vectorized assembly so both methods are comparable
    cov_matrix_vectorized <- cov_vmf(
      omega_grid, t_grid, mu, kappa, distance_type, n_mc_samples, n_cores,
      mc_samples = mc_samples, seed = seed, h0 = h0, unknown_param = unknown_param
    )
    t_vectorized_end <- Sys.time()
    vectorized_time <- as.numeric(difftime(t_vectorized_end, t_vectorized_start, units = "secs"))
    cat("  Vectorized (parallelized):", round(vectorized_time, 2), "seconds\n\n")


    # Use the same seed for vectorized limit simulation to compare identical Gaussian draws
    limit_values_vec <- simulate_limit_gaussian(cov_matrix_vectorized, M, seed = seed)

    # Store runtime information to annotate the plot and return to caller (captured later)
  empirical_data <- list()
  for (n in n_values) {
    cat("Simulating empirical process with n =", n, "...\n")
    empirical_values <- simulate_empirical_process_vmf(
      omega_grid, t_grid, n, mu, kappa, distance_type, M, n_cores, seed = seed, h0 = h0, unknown_param = unknown_param
    )
    empirical_data[[as.character(n)]] <- empirical_values
  }
  
  # Prepare data for plotting (include both limit simulations)
  n_max <- max(n_values)
  
  # Density curves: All n values
  density_data <- data.frame()
  for (n in n_values) {
    n_char <- as.character(n)
    density_data <- rbind(
      density_data,
      data.frame(
        values = empirical_data[[n_char]],
        n = n,
        label = paste0("n=", n)
      )
    )
  }
  
  # Add limit process to density data
  density_data_full <- density_data
  if (!is.null(limit_values_vec)) {
    density_data_full <- rbind(density_data_full, data.frame(values = limit_values_vec[!is.na(limit_values_vec)], n = Inf, label = "G"))
  }
  # Create color palette (limit is black now)
  limit_color <- "#000000"  # vectorized
  n_colors <- length(n_values)
  empirical_colors <- rev(rainbow(n_colors))
  empirical_labels <- paste0("n=", n_values)
  all_labels <- c(empirical_labels, "G")
  all_colors_vec <- c(empirical_colors, limit_color)

  all_colors <- setNames(all_colors_vec, all_labels)
  linetype_vals <- rep("solid", length(all_labels))
  names(linetype_vals) <- all_labels
  # No histogram fills (we no longer draw histograms)

  all_colors <- setNames(all_colors_vec, all_labels)
  linetype_vals <- rep("solid", length(all_labels))
  names(linetype_vals) <- all_labels

  fill_labels <- c("G", paste0("Empirical (n=", n_max, ")"))
  fill_colors_vec <- c(limit_color, "#FF3333")

  # Create plot (moved here so density_data_full exists)
  library(ggplot2)
  p_convergence <- ggplot() +
    geom_density(
      data = density_data_full,
      aes(x = values, color = label, linetype = label),
      linewidth = 1.0,
      adjust = density_adjust,
      trim = TRUE,
      show.legend = TRUE,
      alpha = 0.4
    ) +
    scale_color_manual(values = all_colors) +
    scale_linetype_manual(values = linetype_vals) +
    labs(
      x = "Supremum of the process",
      y = "Density",
      color = "Process",
      linetype = "Process"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19),
      axis.title.y = element_text(size = 19)
    )
  if (!is.null(xlim) && length(xlim) == 2 && !any(is.na(xlim))) {
    p_convergence <- p_convergence + coord_cartesian(xlim = xlim)
  }
  # Rely on ggplot2 automatic legend in the top-right
  p_convergence <- p_convergence +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text = element_text(size = 19),
      legend.title = element_text(size = 19)
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = "solid", shape = NA, size = 2)),
      linetype = guide_legend(override.aes = list(linetype = "solid", shape = NA, size = 2))
    )


  # Statistical comparison
  cat("=== Kolmogorov-Smirnov Tests vs Limit ===\n")
  ks_results <- list()
  for (n in n_values) {
    n_char <- as.character(n)
    if (!is.null(limit_values_vec)) {
      ks_test_vec <- safe_ks_test(limit_values_vec, empirical_data[[n_char]])
      ks_results[[n_char]] <- list(vectorized = ks_test_vec)
      pval_vec <- if (!is.null(ks_test_vec$p.value)) ks_test_vec$p.value else NA_real_
      if (is.na(pval_vec)) {
        cat(sprintf("n=%5d vs Vectorized limit: p-value = NA (test not run / invalid input)\n", n))
      } else {
        cat(sprintf("n=%5d vs Vectorized limit: p-value = %.4f", n, pval_vec))
        if (pval_vec > 0.05) {
          cat(" -> Not significantly different ✓\n")
        } else {
          cat(" -> Significantly different\n")
        }
      }
    } else {
      ks_results[[n_char]] <- list(vectorized = NA)
      cat(sprintf("n=%5d vs Vectorized limit: SKIPPED (limit simulation not available)\n", n))
    }
    
  }
  cat("\n")
  
  return(list(
    limit_values_vec = limit_values_vec,
    cov_matrix_vectorized = cov_matrix_vectorized,
    # NOTE: PSD corrections are not applied internally, so no 'corrected' flags are returned
    empirical_data = empirical_data,
    n_values = n_values,
    omega_grid = omega_grid,
    t_grid = t_grid,
    plot = p_convergence,
    ks_results = ks_results
  ))
}

cat("vMF distribution Gaussian process functions loaded successfully!\n")
