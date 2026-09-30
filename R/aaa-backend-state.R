# State for the optional compiled distance-profile backend.
.distance_profile_cpp_state <- new.env(parent = emptyenv())
.distance_profile_cpp_state$loaded <- FALSE
.distance_profile_cpp_state$exports <- new.env(parent = emptyenv())
.distance_profile_cpp_state$active_backend <- "r"
