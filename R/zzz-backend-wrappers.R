# Original distance-profile wrapper groups.
install_distance_profile_backend_wrappers(
  c(
    "theoretical_distance_profile_normal",
    "distance_profile_small_circle_weighted_mixture2",
    "distance_profile_small_circle_weighted_mixture2_grid"
  ),
  envir = environment(),
  cpp_supported = TRUE
)

install_distance_profile_backend_wrappers(
  c(
    "theoretical_distance_profile_hvmf",
    "theoretical_distance_profile_vmf_s1_chordal",
    "distance_profile_vmf_s2_integral",
    "theoretical_distance_profile_vmf",
    "theoretical_distance_profile_vmf_s2_fast",
    "distance_profile_vmf_s2_grid",
    "distance_profile_vmf_s2_legendre",
    "distance_profile_vmf_s2_legendre_grid",
    "distance_profile_vmf_s2_legendre_cvm_grid",
    "distance_profile_vmf_s2_cvm_grid",
    "distance_profile_spherical_cauchy",
    "distance_profile_spherical_cauchy_grid",
    "distance_profile_spherical_cauchy_cvm_grid",
    "distance_profile_small_circle",
    "distance_profile_small_circle_grid",
    "distance_profile_small_circle_cvm_grid",
    "small_circle_distance_profile_integral",
    "distance_profile_watson",
    "distance_profile_watson_grid",
    "distance_profile_watson_cvm_grid",
    "distance_profile_small_circle_symmetric_mixture2",
    "distance_profile_small_circle_symmetric_mixture2_grid",
    "distance_profile_small_circle_symmetric_mixture2_cvm_grid",
    "distance_profile_beta_mixture2",
    "distance_profile_beta_mixture2_grid",
    "distance_profile_beta_mixture2_cvm_grid",
    "distance_profile_uniform_beta_mixture",
    "distance_profile_uniform_beta_mixture_grid",
    "distance_profile_uniform_beta_mixture_cvm_grid",
    "distance_profile_logitnormal_mixture2",
    "distance_profile_logitnormal_mixture2_grid",
    "distance_profile_logitnormal_mixture2_cvm_grid",
    "rotational_distance_profile_integral"
  ),
  envir = environment(),
  cpp_supported = FALSE
)


install_distance_profile_backend_wrappers(
  c(
    "evaluate_mvnorm_distance_profile",
    "evaluate_mvnorm_distance_profile_matrix",
    "axial_truncnorm_mixture_distance_profile",
    "axial_truncnorm_mixture_distance_profile_grid"
  ),
  envir = environment(),
  cpp_supported = FALSE
)

install_distance_profile_backend_wrappers(
  "theoretical_distance_profile_cardioid",
  envir = environment(),
  cpp_supported = FALSE
)
