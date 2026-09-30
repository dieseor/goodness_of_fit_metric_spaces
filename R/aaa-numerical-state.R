# Mutable caches and numerical constants from utils.R.
vmf_s2_legendre_cache <- new.env(parent = emptyenv())
jp_cache_env <- new.env(parent = emptyenv())
jp_vmf_near_zero_abs_kappa_psi_default <- 0.001
small_circle_gauss_legendre_cache <- new.env(parent = emptyenv())
rotational_gauss_hermite_cache <- new.env(parent = emptyenv())
