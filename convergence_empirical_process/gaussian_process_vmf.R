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
# SPHERE GRID GENERATION
# ============================================================================

#' Generate canonical lattice on unit sphere (Fibonacci sphere)
#' Based on golden ratio spiral for quasi-uniform distribution
#' @param n Number of points to generate
#' @return Matrix of size (n, 3) with unit vectors on sphere
generate_canonical_lattice <- function(n, dim = 3) {
  #' Generate canonical lattice on unit sphere
  #' If dim == 3, use Fibonacci sphere (quasi-uniform)
  #' If dim > 3, use normalized Gaussian vectors
  #' @param n Number of points
  #' @param dim Ambient dimension (default 3)
  #' @return Matrix of size (n, dim) with unit vectors on sphere
  if (dim == 3) {
    golden_ratio <- (1 + sqrt(5)) / 2
    i <- seq(0, n - 1)
    theta <- 2 * pi * i / golden_ratio
    phi <- acos(1 - 2 * (i + 0.5) / n)
    x <- cos(theta) * sin(phi)
    y <- sin(theta) * sin(phi)
    z <- cos(phi)
    mat <- cbind(x, y, z)
    return(mat)
  } else if (dim > 1) {
    # Simple: normalizar vectores gaussianos
    mat <- matrix(rnorm(n * dim), nrow = n, ncol = dim)
    mat <- t(apply(mat, 1, function(v) v / sqrt(sum(v^2))))
    return(mat)
  } else {
    stop("Dimension must be >= 2")
  }
}

# ============================================================================
# DISTANCE FUNCTIONS
# ============================================================================

#' Check dot products are valid (internal helper)
#' Throws error if values are outside [-1, 1]
#' @param dot_products Scalar or vector/matrix of dot products
#' @param tolerance Tolerance for checking (default 1e-10)
#' @return Original dot products (no modification!)
check_dot_products <- function(dot_products, tolerance = 1e-10) {
  if (any(dot_products < -1 - tolerance) || any(dot_products > 1 + tolerance)) {
    bad_values <- dot_products[dot_products < -1 - tolerance | dot_products > 1 + tolerance]
    stop(paste0("Invalid dot products detected: range [", 
                round(min(bad_values), 10), ", ", round(max(bad_values), 10), 
                "]. Check that inputs are unit vectors."))
  }
  # Return unchanged - NO CLAMPING!
  return(dot_products)
}

#' Compute distance between points on sphere
#' @param omega1 Point on sphere (vector of length q)
#' @param omega2 Point on sphere (vector of length q)
#' @param distance_type Either "chordal" or "geodesic"
#' @return Distance value
sphere_distance <- function(omega1, omega2, distance_type = "chordal") {
  omega1 <- omega1 / sqrt(sum(omega1^2))
  omega2 <- omega2 / sqrt(sum(omega2^2))
  dot_product <- sum(omega1 * omega2)
  dot_product <- pmax(pmin(dot_product, 1), -1)  # Clamp to [-1, 1]
  
  if (distance_type == "chordal") {
    return(sqrt(2 * (1 - dot_product)))
  } else if (distance_type == "geodesic") {
    return(acos(dot_product))
  } else {
    stop("distance_type must be 'chordal' or 'geodesic'")
  }
}


#' Precompute common intermediates for vMF vectorized covariance
#' This helper centralizes computation of dists_all, F2_matrix, E2_array, E2_mat, F2_vec, m2_mat and in_ball_list.
#' @param mc_samples Monte Carlo sample (n_mc × q)
#' @param omega_grid Matrix (n_omega × q)
#' @param t_grid Vector (n_t)
#' @param mu Mean direction (q+1)
#' @param kappa Concentration parameter
#' @param distance_type Either 'chordal' or 'geodesic'
#' @param A_q_kappa Precomputed A_q(kappa)
#' @return A named list with dists_all, F2_matrix, E2_array, E2_mat, F2_vec, m2_mat, in_ball_list
compute_precomp_vmf <- function(mc_samples, omega_grid, t_grid, mu, kappa, distance_type, A_q_kappa) {
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  q <- length(mu) - 1
  n_mc <- nrow(mc_samples)

  # Dot products and distances from all MC samples to all omegas
  dots_all <- mc_samples %*% t(omega_grid)
  if (distance_type != "chordal") dots_all <- check_dot_products(dots_all)
  dists_all <- if (distance_type == "chordal") sqrt(2 * (1 - dots_all)) else acos(dots_all)

  # F2_matrix: marginal probabilities for each omega2 and t
  F2_matrix <- t(apply(omega_grid, 1, function(omega2) {
    theoretical_distance_profile_vmf(omega2, mu, kappa, t_grid, distance_type)
  }))

  # Build E2_mat (n_total x q+1) using a simpler and explicit rbind-by-t approach
  in_ball_list <- vector("list", n_t)
  E2_rows_by_t <- lapply(seq_len(n_t), function(k) {
    in_ball_k <- dists_all <= t_grid[k]
    in_ball_list[[k]] <<- in_ball_k
    counts_k <- colSums(in_ball_k)
    if (all(counts_k == 0)) return(matrix(0, nrow = n_omega, ncol = q + 1))
    sums_mat <- t(mc_samples) %*% (in_ball_k * 1) # (q+1) x n_omega
    counts_k_safe <- counts_k
    counts_k_safe[counts_k_safe == 0] <- 1
    means_mat <- sweep(sums_mat, 2, counts_k_safe, "/")
    if (any(counts_k == 0)) means_mat[, counts_k == 0] <- 0
    # return n_omega x q+1 block with rows ordered by omega index
    return(t(means_mat))
  })
  E2_mat <- do.call(rbind, E2_rows_by_t)
  # Build E2_array (n_omega x n_t x q+1) for compatibility with older code paths
  E2_array <- array(NA, dim = c(n_omega, n_t, q + 1))
  for (k in seq_len(n_t)) {
    E2_array[, k, ] <- E2_rows_by_t[[k]]
  }

  F2_vec <- as.vector(F2_matrix)
  m2_mat <- sweep(E2_mat, 2, A_q_kappa * mu, "-")

  return(list(
    dists_all = dists_all,
    F2_matrix = F2_matrix,
    E2_array = E2_array,
    E2_mat = E2_mat,
    F2_vec = F2_vec,
    m2_mat = m2_mat,
    in_ball_list = in_ball_list
  ))
}

#' Compute distance matrix between sample and grid (vectorized)
#' @param sample_matrix Matrix (n × q) of sample points
#' @param omega_grid Matrix (m × q) of grid points
#' @param distance_type Either "chordal" or "geodesic"
#' @return Distance matrix (n × m)
compute_distance_matrix_vmf <- function(sample_matrix, omega_grid, distance_type = "chordal") {
  # Compute dot products: (n × q) %*% (q × m) = (n × m)
  dot_products <- sample_matrix %*% t(omega_grid)
  dot_products <- check_dot_products(dot_products)
  
  if (distance_type == "chordal") {
    return(sqrt(2 * (1 - dot_products)))
  } else if (distance_type == "geodesic") {
    # Validate dot products before invoking acos to ensure valid domain
    dot_products <- check_dot_products(dot_products)
    return(acos(dot_products))
  } else {
    stop("distance_type must be 'chordal' or 'geodesic'")
  }
}

# ============================================================================
# DISTANCE PROFILE AND PROBABILITIES
# ============================================================================

# NOTE: theoretical_distance_profile_vmf is defined in utils.R (original from vmf_distance_profile_analysis.R)

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
  # Debug prints removed for production
  
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
  # Debug prints removed for production
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
    # the same code as the master; avoid sourcing local debug variants on workers.
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
  # No propagation of debug options to workers; debug code removed

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
    "compute_joint_probability_vmf", "compute_conditional_expectation_vmf",
    "row_cov_vmf", "compute_precomp_vmf",
    "omega_grid", "t_grid", "mu", "kappa", "distance_type",
    "mc_samples", "A_q_kappa", "var_X", "precomp", "omega_idx_vec", "t_idx_vec",
    "h0", "h0_val", "unknown_param"
  ), envir = environment())
  # Export control flags
  clusterExport(cl, c("upper_triangle", "n_total"), envir = environment())
  
  start_time <- Sys.time()

  # Non-parallel test removed after debugging
  
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
visualize_convergence_to_limit_vmf <- function(n_values = c(10, 50, 100, 500),
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
                                              verbose = FALSE,
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
  if (verbose) cat('DEBUG: After generating omega grid; n_omega =', nrow(omega_grid), '\n')

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
  if (verbose) cat('DEBUG: After empirical process simulation. n_values length =', length(n_values), '\n')
  
  # Prepare data for plotting (include both limit simulations)
  n_max <- max(n_values)
  hist_values <- empirical_data[[as.character(n_max)]]
  hist_process <- rep(paste0("Empirical (n=", n_max, ")"), each = M)
  # add vectorized limit
  if (!is.null(limit_values_vec)) {
    hist_values <- c(hist_values, limit_values_vec[!is.na(limit_values_vec)])
    hist_process <- c(hist_process, rep("Limit Gaussian (Vectorized)", each = sum(!is.na(limit_values_vec))))
  }

  histogram_data <- data.frame(values = hist_values, process = hist_process)
  if (verbose) cat('DEBUG: histogram_data rows =', nrow(histogram_data), '\n')
  
  # Density curves: All n values
  density_data <- data.frame()
  if (verbose) cat('DEBUG: Created empty density_data\n')
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
  if (verbose) cat('DEBUG: Built density_data_full rows =', nrow(density_data_full), '\n')
  if (verbose) cat('DEBUG: After assign: density_data exists? ', exists('density_data'), ' rows=', nrow(density_data), '\n')
  if (verbose) cat('DEBUG: After assign: density_data_full exists? ', exists('density_data_full'), ' rows=', nrow(density_data_full), '\n')
  if (!is.null(limit_values_vec)) {
    density_data_full <- rbind(density_data_full, data.frame(values = limit_values_vec[!is.na(limit_values_vec)], n = Inf, label = "Limit Gaussian (Vectorized)"))
  }
  # Create color palette (must be before building the plot)
  limit_color <- "#0066CC"  # vectorized
  n_colors <- length(n_values)
  empirical_colors <- rev(rainbow(n_colors))
  empirical_labels <- paste0("n=", n_values)
  all_labels <- c(empirical_labels, "Limit Gaussian (Vectorized)")
  all_colors_vec <- c(empirical_colors, limit_color)

  all_colors <- setNames(all_colors_vec, all_labels)
  linetype_vals <- rep("solid", length(all_labels))
  names(linetype_vals) <- all_labels
  # Fill palette for histogram (limit and empirical)
  fill_labels <- c("Limit Gaussian (Vectorized)", paste0("Empirical (n=", n_max, ")"))
  fill_colors_vec <- c(limit_color, "#FF3333")

  # Reassign palette mappings after possible h0 transformations
  all_colors <- setNames(all_colors_vec, all_labels)
  linetype_vals <- rep("solid", length(all_labels))
  names(linetype_vals) <- all_labels

  fill_labels <- c("Limit Gaussian (Vectorized)", paste0("Empirical (n=", n_max, ")"))
  fill_colors_vec <- c(limit_color, "#FF3333")

  # Create plot (moved here so density_data_full exists)
  library(ggplot2)
  p_convergence <- ggplot() +
    geom_histogram(
      data = histogram_data,
      aes(x = values, y = after_stat(density), fill = process),
      alpha = 0.3,
      position = "identity",
      bins = 50,
      show.legend = FALSE
    ) +
    geom_density(
      data = density_data_full,
      aes(x = values, color = label, linetype = label),
      linewidth = 1.0,
      adjust = density_adjust,
      trim = TRUE,
      show.legend = FALSE,
      alpha = 0.4
    ) +
    scale_color_manual(values = all_colors) +
    scale_linetype_manual(values = linetype_vals) +
    scale_fill_manual(values = setNames(fill_colors_vec, fill_labels)) +
    labs(
      x = "Supremum of the process",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19),
      axis.title.y = element_text(size = 19)
    ) +
    # Apply x-axis limits if provided, otherwise keep the default auto-scale
    if (!is.null(xlim) && length(xlim) == 2 && !any(is.na(xlim))) {
      p_convergence <- p_convergence + coord_cartesian(xlim = xlim)
    }
  
  # Manual legend positioning: compute positions relative to x-axis data limits
  # Determine max x value that should be visible in the plot
  if (!is.null(xlim) && length(xlim) == 2 && !any(is.na(xlim))) {
    x_upper <- xlim[2]
  } else {
    # compute from observed histogram and limit draws
    x_candidates <- c(histogram_data$values)
    if (!is.null(limit_values_vec)) x_candidates <- c(x_candidates, limit_values_vec[!is.na(limit_values_vec)])
    # fallback to at least 2.0 as a reasonable default (chordal's upper bound)
    x_data_max <- if (length(x_candidates) > 0) max(x_candidates, na.rm = TRUE) else 2.0
    x_upper <- max(2.0, x_data_max * 1.05) # small buffer so contents don't overlap legend
  }
  # Legend position is now a fraction of x_upper to keep it within plot bounds
  legend_x <- x_upper * 0.76  
  legend_y_start <- 1.55
  legend_spacing <- 0.08  
  # For composite null, add 0.02 more vertical spacing
  if (!is.null(h0) && h0 == 'composite') {
    legend_spacing <- legend_spacing + 0.02
  }
  legend_line_length <- x_upper * 0.03
  legend_rect_height <- 0.023
  # rectangle half-width for legend boxes scaled to data range, with a minimum
  legend_rect_half_width <- max(0.2, x_upper * 0.12)
  legend_line_x_start <- legend_x - legend_rect_half_width * 0.35
  # Increase horizontal separation between color line/rect and text by 0.3 units
  legend_text_x <- legend_line_x_start + legend_line_length + 0.05
  
  # Add legend elements (build a single sequence to ensure homogeneous spacing)
  p_convergence <- p_convergence +
    annotate("rect",
             xmin = legend_x - legend_rect_half_width, xmax = legend_x + legend_rect_half_width,
             ymin = legend_y_start - (n_colors + 2) * legend_spacing - 0.03,
             ymax = legend_y_start + 0.02,
             fill = "transparent")

  # Build list of legend entries in order: empirical densities, limits, histograms
  legend_entries <- list()
  for (i in seq_along(n_values)) {
    legend_entries[[length(legend_entries) + 1]] <- list(type = 'density', label = paste0('italic(n)==', n_values[i]), color = all_colors[i], parse = TRUE)
  }
  # Add limit entries
  # Use plotmath label "'𝔾'[mu[theta[0]]]" for vectorized limit in all cases (keeps original look)
  # Use plotmath rendering for the vectorized limit label in all cases (including simple)
  limit_text_label <- "'𝔾'[mu[theta[0]]]"
  limit_text_parse <- TRUE
  # Add limit entries into the legend_entries vector
  legend_entries[[length(legend_entries) + 1]] <- list(type = 'limit', label = limit_text_label, color = limit_color, parse = limit_text_parse)

  # Add histogram fills (e.g., empirical) as legend entries (red box) - put them after limits
  # Keep the label matching the italic(n) style used for density entries
  legend_entries[[length(legend_entries) + 1]] <- list(type = 'hist', label = paste0('italic(n)==', n_max), color = '#FF3333')

  # Now compute y positions for all entries uniformly
  num_entries <- length(legend_entries)
  y_positions <- legend_y_start - (seq_len(num_entries) - 1) * legend_spacing

  # Draw each legend entry by type
  for (idx in seq_along(legend_entries)) {
    entry <- legend_entries[[idx]]
    y_pos <- y_positions[idx]
    if (entry$type == 'density') {
      p_convergence <- p_convergence +
        annotate('segment', x = legend_line_x_start, xend = legend_line_x_start + legend_line_length, y = y_pos, yend = y_pos, color = entry$color, linetype = 'solid', linewidth = 1.0) +
        annotate('text', x = legend_text_x, y = y_pos, label = entry$label, parse = isTRUE(entry$parse), hjust = 0, size = 5.5, color = 'black')
    } else if (entry$type == 'limit') {
      if (isTRUE(entry$parse)) {
        p_convergence <- p_convergence +
          annotate('segment', x = legend_line_x_start, xend = legend_line_x_start + legend_line_length, y = y_pos, yend = y_pos, color = entry$color, linetype = 'solid', linewidth = 1.0) +
          annotate('text', x = legend_text_x, y = y_pos, label = entry$label, parse = TRUE, hjust = 0, size = 5.5, color = 'black')
      } else {
        p_convergence <- p_convergence +
          annotate('segment', x = legend_line_x_start, xend = legend_line_x_start + legend_line_length, y = y_pos, yend = y_pos, color = entry$color, linetype = 'solid', linewidth = 1.0) +
          annotate('text', x = legend_text_x, y = y_pos, label = entry$label, parse = FALSE, hjust = 0, size = 5.5, color = 'black')
      }
    } else if (entry$type == 'hist') {
      p_convergence <- p_convergence +
        annotate('rect', xmin = legend_line_x_start, xmax = legend_line_x_start + legend_line_length, ymin = y_pos - legend_rect_height, ymax = y_pos + legend_rect_height, fill = entry$color, alpha = 0.3, color = NA) +
        annotate('text', x = legend_text_x, y = y_pos, label = entry$label, parse = TRUE, hjust = 0, size = 5.5, color = 'black')
    }
  }


  # Limit Gaussian histogram (it didn't appear automatically, so I add it manually)
  # Histogram boxes
  y_pos_hist1 <- legend_y_start - (n_colors + 2) * legend_spacing
  p_convergence <- p_convergence +
    annotate("rect",
         xmin = legend_line_x_start, xmax = legend_line_x_start + legend_line_length,
             ymin = y_pos_hist1 - legend_rect_height, ymax = y_pos_hist1 + legend_rect_height,
             fill = limit_color, alpha = 0.3, color = NA) +
    annotate("text",
         x = legend_text_x,
             y = y_pos_hist1,
             label = "'𝔾'[mu[theta[0]]]",
             parse = TRUE,
             hjust = 0,
             size = 5.5,
             color = "black")
  if (verbose) cat('DEBUG: legend y_positions =', paste(round(y_positions, 4), collapse = ','), '\n')

  # Note: any PSD corrections are not applied internally; when a covariance is non-PSD the simulation will error during sampling and you must fix the covariance externally.
  # (Removed automatic correction annotation) No annotation position computation needed
  # No correction annotations are applied; when a covariance is non-PSD sampling will error and you must fix the covariance upstream.
  
    # Removed older histogram coordinate computations; the histogram legend entry is drawn above in the unified legend loop.
  
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
