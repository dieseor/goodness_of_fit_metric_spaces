`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/Diego/Desktop/Codigo/goodness_of_fit_metric_spaces/tests/stress_launcher_smoke/cache/R/goodness_of_fit_metric_spaces/sourceCpp/sourceCpp-x86_64-apple-darwin20-1.0.14/sourcecpp_17e6f2963b538/sourceCpp_2.so')

cpp_dp_legendre_matrix <- Rcpp:::sourceCppFunction(function(x, l_max) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_dp_legendre_matrix')
cpp_dp_projection_cdf_legendre_matrix <- Rcpp:::sourceCppFunction(function(x_matrix, r, coefficients, enforce_bounds = TRUE) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_dp_projection_cdf_legendre_matrix')
cpp_dp_normal_profile <- Rcpp:::sourceCppFunction(function(omega, t_values, mu, sigma) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_dp_normal_profile')
cpp_dp_normal_profile_matrix <- Rcpp:::sourceCppFunction(function(omega, t_matrix, mu, sigma) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_dp_normal_profile_matrix')
cpp_dp_weighted_sample_profile_rows <- Rcpp:::sourceCppFunction(function(order_matrix, rank_matrix, normalized_weights) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_dp_weighted_sample_profile_rows')
cpp_dp_weighted_sample_profile_linear <- Rcpp:::sourceCppFunction(function(order_matrix, rank_linear_index, normalized_weights) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_dp_weighted_sample_profile_linear')
cpp_profile_lookup_tensor_local_polynomial <- Rcpp:::sourceCppFunction(function(values, t_grid, geometry_grid, kappa_grid, t_query, geometry_query, kappa_query, stencil_size = 6L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_profile_lookup_tensor_local_polynomial')
cpp_fast_sample_ks_cvm_stats <- Rcpp:::sourceCppFunction(function(centered_weights, score_block, obs_order_matrix, tie_end_matrix, correction_matrix, scale_factor, compute_ks = TRUE, compute_cvm = TRUE) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_fast_sample_ks_cvm_stats')
cpp_fast_sample_ks_cvm_stats_contiguous_double <- Rcpp:::sourceCppFunction(function(centered_weights, score_block, obs_order_matrix, tie_end_matrix, correction_matrix, scale_factor, compute_ks = TRUE, compute_cvm = TRUE) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_fast_sample_ks_cvm_stats_contiguous_double')
cpp_sunspots_joint_profile_block <- Rcpp:::sourceCppFunction(function(radii, rho, center_s, time_nodes, time_weights, coefficients) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_sunspots_joint_profile_block')

rm(`.sourceCpp_1_DLLInfo`)
