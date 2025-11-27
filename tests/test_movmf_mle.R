# Test of movMF MLE estimation for different kappa values
# Investigating why supervisor insists on kappa=1

library(movMF)
library(sphunif) 
library(pracma)

# Function to test MLE estimation
test_mle_estimation <- function(mu_true, kappa_true, n_samples = 1000, n_tests = 10) {
  
  cat("=== Testing MLE for mu =", paste(round(mu_true, 4), collapse = ", "), 
      ", kappa =", kappa_true, "===\n")
  
  results <- data.frame(
    test = 1:n_tests,
    kappa_estimated = numeric(n_tests),
    mu_norm_error = numeric(n_tests),
    convergence = logical(n_tests)
  )
  
  for (i in 1:n_tests) {
    # Generate sample from vMF
    set.seed(123 + i)
    X <- r_vMF(n_samples, mu_true, kappa_true)
    
    # Estimate parameters using movMF
    tryCatch({
      mle_result <- movMF(X, k = 1)  # k=1 means single vMF component
      
      # Extract estimates
      theta_hat <- mle_result$theta[[1]]  # Estimated theta (mu * kappa)
      kappa_hat <- norm(theta_hat, type = "2")  # ||theta|| = kappa
      mu_hat <- theta_hat / kappa_hat  # mu = theta / ||theta||
      
      # Store results
      results$kappa_estimated[i] <- kappa_hat
      results$mu_norm_error[i] <- norm(mu_hat - mu_true, type = "2")
      results$convergence[i] <- TRUE
      
    }, error = function(e) {
      cat("Error in test", i, ":", e$message, "\n")
      results$kappa_estimated[i] <- NA
      results$mu_norm_error[i] <- NA
      results$convergence[i] <- FALSE
    })
  }
  
  # Summary statistics
  converged_results <- results[results$convergence, ]
  if (nrow(converged_results) > 0) {
    cat("Converged tests:", nrow(converged_results), "/", n_tests, "\n")
    cat("Kappa - True:", kappa_true, ", Mean estimated:", 
        round(mean(converged_results$kappa_estimated, na.rm = TRUE), 4), "\n")
    cat("Kappa - Bias:", 
        round(mean(converged_results$kappa_estimated, na.rm = TRUE) - kappa_true, 4), "\n")
    cat("Kappa - SD:", 
        round(sd(converged_results$kappa_estimated, na.rm = TRUE), 4), "\n")
    cat("Mu error - Mean:", 
        round(mean(converged_results$mu_norm_error, na.rm = TRUE), 4), "\n")
    cat("Mu error - SD:", 
        round(sd(converged_results$mu_norm_error, na.rm = TRUE), 4), "\n")
  } else {
    cat("NO CONVERGENCE in any test!\n")
  }
  cat("\n")
  
  return(results)
}

# Test different kappa values
mu_test <- c(1, 0, 0)  # Simple mu on S^2

cat("TESTING DIFFERENT KAPPA VALUES WITH movMF\n")
cat("==========================================\n\n")

# Test various kappa values
kappa_values <- c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0)

all_results <- list()
for (k in kappa_values) {
  all_results[[paste0("kappa_", k)]] <- test_mle_estimation(mu_test, k, n_samples = 1000)
}

# Let's also check what movMF actually optimizes
cat("=== CHECKING movMF DOCUMENTATION AND METHOD ===\n")

# Generate a sample with kappa = 2
set.seed(123)
X_test <- r_vMF(1000, mu_test, 2.0)

# Fit with movMF
fit_result <- movMF(X_test, k = 1)

cat("movMF fitted object structure:\n")
str(fit_result)

cat("\nmovMF theta (estimated):\n")
print(fit_result$theta[[1]])

cat("\nEstimated kappa (||theta||):", norm(fit_result$theta[[1]], type = "2"), "\n")
cat("True kappa:", 2.0, "\n")

# Let's check if movMF is actually doing MLE or something else
cat("\n=== CHECKING IF movMF USES TRUE MLE ===\n")

# The theoretical MLE for vMF can be computed analytically:
# mu_mle = (sum of X_i) / ||sum of X_i||
# kappa_mle is the solution to: A_d(kappa) = ||sum of X_i|| / n
# where A_d(kappa) = I_{d/2}(kappa) / I_{d/2-1}(kappa) is the ratio of modified Bessel functions

sample_mean <- colMeans(X_test)
R_bar <- norm(sample_mean, type = "2")  # ||sample mean||
mu_mle_theoretical <- sample_mean / R_bar

cat("Theoretical MLE mu:", paste(round(mu_mle_theoretical, 4), collapse = ", "), "\n")
cat("movMF estimated mu:", paste(round(fit_result$theta[[1]] / norm(fit_result$theta[[1]], type = "2"), 4), collapse = ", "), "\n")
cat("Difference in mu:", norm(mu_mle_theoretical - fit_result$theta[[1]] / norm(fit_result$theta[[1]], type = "2"), type = "2"), "\n")

# For kappa, we need to solve the equation A_d(kappa) = R_bar
# This is typically done numerically
d <- length(mu_test)  # dimension
A_d_target <- R_bar  # target value for A_d(kappa)

cat("R_bar (||sample mean||):", R_bar, "\n")
cat("Target A_d(kappa) value:", A_d_target, "\n")

# Let's see what A_d function looks like for the estimated kappa
kappa_estimated <- norm(fit_result$theta[[1]], type = "2")
A_d_estimated <- besselI(kappa_estimated, d/2) / besselI(kappa_estimated, d/2 - 1)

cat("Estimated kappa:", kappa_estimated, "\n")
cat("A_d(estimated kappa):", A_d_estimated, "\n")
cat("Should equal R_bar if true MLE:", R_bar, "\n")
cat("Difference:", abs(A_d_estimated - R_bar), "\n")