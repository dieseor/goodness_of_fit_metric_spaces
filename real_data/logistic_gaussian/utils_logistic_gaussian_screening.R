Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

suppressPackageStartupMessages({
  library(compositions)
  library(MASS)
  library(ggplot2)
})

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) {
    rhs
  } else {
    lhs
  }
}

logistic_gaussian_screening_warning <- paste(
  "This workflow performs the logistic Gaussian goodness-of-fit analysis through the global bootstrap pipeline used by the calibration code.",
  "The default mode is the composite multiplier bootstrap for the logistic Gaussian model, with weighted re-estimation of the unknown parameters in every bootstrap replicate.",
  "For the logistic Gaussian model with Aitchison distance, theoretical distance profiles are evaluated through the global logistic Gaussian model specification used by the calibration code, not by local screening-specific profile code.",
  "The legacy labels composite_parametric and plugin_simple_null are accepted only for backwards compatibility; the screening analysis should be run with bootstrap_mode = 'composite_multiplier'."
)

default_logistic_gaussian_screening_datasets <- function(include_external = TRUE) {
  datasets <- c(
    "SkyeAFM",
    "SkyeLavasComplete",
    "AarMajorOxides",
    "Sediments",
    "HouseholdExp",
    "ClamEast",
    "ClamWest",
    "ClamCombined"
  )
  if (isTRUE(include_external)) {
    datasets <- c(
      datasets,
      "FerrettiGut",
      "FerrettiOral",
      "Shi2015",
      "HongKongBudgetsA",
      "HongKongBudgetsB",
      "HongKongBudgetsCombined"
    )
  }
  datasets
}

logistic_gaussian_screening_dataset_registry <- function() {
  list(
    SkyeAFM = list(
      source_type = "compositions_data",
      candidate_names = "SkyeAFM",
      candidate_packages = "compositions",
      compositional_columns = c("A", "F", "M"),
      notes = c(
        "AFM percentages for 23 aphyric Skye lavas.",
        "Rows already sum to 100 before closure."
      )
    ),
    SkyeLavasComplete = list(
      source_type = "search_r_packages",
      candidate_names = c("Skye", "skye", "SkyeLavas", "skyeLavas", "SkyeLava", "skyeLava", "lavas", "Lava"),
      candidate_packages = c("compositions", "robCompositions", "compositional", "Hotelling", "bayesm"),
      compositional_columns = NULL,
      expected_n = 32L,
      expected_D = 10L,
      notes = c(
        "Complete Skye lavas dataset from Aitchison/Thompson--Esson--Duncan if reproducibly available.",
        "If the dataset cannot be found in installed packages, the result is recorded as not_found."
      )
    ),
    AarMajorOxides = list(
      source_type = "compositions_data",
      candidate_names = "Aar",
      candidate_packages = "compositions",
      compositional_columns = c("SiO2", "TiO2", "Al2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "Fe2O3t"),
      notes = c(
        "Major-oxide geochemical compositions from the Aar massif.",
        "Selected oxide parts are treated compositionally after closure."
      )
    ),
    Sediments = list(
      source_type = "compositions_data",
      candidate_names = "Sediments",
      candidate_packages = "compositions",
      compositional_columns = c("sand", "silt", "clay"),
      notes = c(
        "Classical ternary sediment composition dataset.",
        "The sample type metadata column is excluded from the compositional analysis."
      )
    ),
    HouseholdExp = list(
      source_type = "compositions_data",
      candidate_names = "HouseholdExp",
      candidate_packages = "compositions",
      compositional_columns = c("Housing", "Food", "Other", "Services"),
      notes = c(
        "Canonical household budget shares dataset used in compositional PCA examples.",
        "The `Sex` metadata column is excluded from the compositional analysis."
      )
    ),
    Aar = list(
      source_type = "compositions_data",
      candidate_names = "Aar",
      candidate_packages = "compositions",
      compositional_columns = c("SiO2", "TiO2", "Al2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "Fe2O3t"),
      notes = c(
        "Alias entry for the Aar major-oxide geochemical compositions.",
        "Selected oxide parts are treated compositionally after closure."
      )
    ),
    Metabolites = list(
      source_type = "compositions_data",
      candidate_names = "Metabolites",
      candidate_packages = "compositions",
      compositional_columns = c("met1", "met2", "met3"),
      notes = c(
        "Steroid metabolite compositions.",
        "The class label column `Type` is excluded from the compositional analysis."
      )
    ),
    SerumProtein = list(
      source_type = "compositions_data",
      candidate_names = "SerumProtein",
      candidate_packages = "compositions",
      compositional_columns = c("a", "b", "c", "d"),
      notes = c(
        "Serum protein compositions.",
        "The class label column `Type` is excluded from the compositional analysis."
      )
    ),
    WhiteCells = list(
      source_type = "compositions_data",
      candidate_names = "WhiteCells",
      candidate_packages = "compositions",
      compositional_columns = c("mG", "mL", "mM", "iG", "iL", "iM"),
      notes = c(
        "White-cell composition profiles.",
        "All six columns are treated compositionally after closure."
      )
    ),
    Boxite = list(
      source_type = "compositions_data",
      candidate_names = "Boxite",
      candidate_packages = "compositions",
      compositional_columns = c("A", "B", "C", "D", "E"),
      notes = c(
        "Boxite compositions.",
        "The depth covariate is excluded from the compositional analysis."
      )
    ),
    arcticLake = list(
      source_type = "search_r_packages",
      candidate_names = c("arcticLake", "ArcticLake"),
      candidate_packages = c("robCompositions", "compositions"),
      compositional_columns = c("sand", "silt", "clay"),
      notes = c(
        "Arctic lake sediment ternary compositions (robCompositions/compositions).",
        "Depth metadata is excluded when present."
      )
    ),
    expenditures = list(
      source_type = "search_r_packages",
      candidate_names = "expenditures",
      candidate_packages = c("robCompositions", "Compositional"),
      compositional_columns = c("housing", "foodstuffs", "alcohol", "other", "services"),
      notes = c(
        "Synthetic household expenditure composition toy dataset.",
        "All five expenditure categories are treated compositionally."
      )
    ),
    expendituresEU = list(
      source_type = "search_r_packages",
      candidate_names = "expendituresEU",
      candidate_packages = c("robCompositions", "Compositional"),
      compositional_columns = c("Food", "Alcohol", "Clothing", "Housing", "Furnishings", "Health", "Transport", "Communications", "Recreation", "Education", "Restaurants", "Other"),
      notes = c(
        "EU expenditure compositions.",
        "All expenditure categories are treated compositionally."
      )
    ),
    coffee = list(
      source_type = "search_r_packages",
      candidate_names = "coffee",
      candidate_packages = c("robCompositions", "Compositional"),
      compositional_columns = c("sort", "acit", "metpyr", "furfu", "furfualc", "dimeth", "met5"),
      notes = c(
        "Coffee composition dataset from robCompositions.",
        "All seven columns are treated compositionally."
      )
    ),
    alcohol = list(
      source_type = "search_r_packages",
      candidate_names = "alcohol",
      candidate_packages = c("robCompositions", "Compositional"),
      compositional_columns = c("beer", "wine", "spirits", "other"),
      notes = c(
        "Alcohol consumption composition by beverage type.",
        "Metadata columns (`country`, `year`) are excluded from the compositional analysis."
      )
    ),
    ArcticLake = list(
      source_type = "compositions_data",
      candidate_names = "ArcticLake",
      candidate_packages = "compositions",
      compositional_columns = c("sand", "silt", "clay"),
      notes = c(
        "Sediment composition data; depth is excluded from the compositional analysis.",
        "Closure is applied because rows are only approximately summing to 100.",
        "Aitchison indicates that marginal logistic-normality is not the natural model here; regression on depth is more appropriate."
      )
    ),
    ClamEast = list(
      source_type = "compositions_data",
      candidate_names = "ClamEast",
      candidate_packages = "compositions",
      compositional_columns = c("dl", "dm", "ds", "ll", "lm", "ls"),
      notes = c(
        "East Bay clam color-size compositions.",
        "No automatic zero replacement is performed in the main analysis."
      )
    ),
    ClamWest = list(
      source_type = "compositions_data",
      candidate_names = "ClamWest",
      candidate_packages = "compositions",
      compositional_columns = c("dl", "dm", "ds", "ll", "lm", "ls"),
      notes = c(
        "West Bay clam color-size compositions.",
        "No automatic zero replacement is performed in the main analysis."
      )
    ),
    ClamCombined = list(
      source_type = "combined_prepared",
      components = c("ClamEast", "ClamWest"),
      group_variable = "site",
      notes = c(
        "Combined East and West clam color-size compositions.",
        "A site label is retained for ilr diagnostics because the combined sample may be a mixture."
      )
    ),
    FerrettiGut = list(
      source_type = "curated_metagenomic_data",
      study = "FerrettiP_2018",
      body_site = "stool",
      infant_timepoint = "Day 1",
      top_taxa = 4L,
      notes = c(
        "Ferretti 2018 gut/stool microbiome dataset from curatedMetagenomicData, if available.",
        "The target subset follows Fang--Subedi: adults plus infant Day 1 samples, genus level, top 4 taxa plus Other."
      )
    ),
    FerrettiOral = list(
      source_type = "curated_metagenomic_data",
      study = "FerrettiP_2018",
      body_site = "oral",
      infant_timepoint = "Day 3",
      top_taxa = 4L,
      notes = c(
        "Ferretti 2018 oral microbiome dataset from curatedMetagenomicData, if available.",
        "The target subset follows Fang--Subedi: adults plus infant Day 3 samples, genus level, top 4 taxa plus Other."
      )
    ),
    Shi2015 = list(
      source_type = "curated_metagenomic_data",
      study = "ShiB_2015",
      top_taxa = 4L,
      notes = c(
        "Shi 2015 subgingival microbiome dataset from curatedMetagenomicData, if available.",
        "The target subset follows Fang--Subedi: periodontitis and recovered samples, genus level, top 4 taxa plus Other."
      )
    ),
    HongKongBudgetsA = list(
      source_type = "search_r_packages",
      candidate_names = c("HongKong", "HongKongBudgets", "Household", "HouseholdBudgets", "Budgets", "Expenditure", "Expenditures", "HK", "HKBudgets"),
      candidate_packages = c("compositions", "robCompositions", "compositional"),
      subgroup = "A",
      notes = c(
        "Hong Kong household expenditure budgets, housing category A, if reproducibly available.",
        "Aitchison reports 41 households in category A."
      )
    ),
    HongKongBudgetsB = list(
      source_type = "search_r_packages",
      candidate_names = c("HongKong", "HongKongBudgets", "Household", "HouseholdBudgets", "Budgets", "Expenditure", "Expenditures", "HK", "HKBudgets"),
      candidate_packages = c("compositions", "robCompositions", "compositional"),
      subgroup = "B",
      notes = c(
        "Hong Kong household expenditure budgets, housing category B, if reproducibly available.",
        "Aitchison reports 42 households in category B."
      )
    ),
    HongKongBudgetsCombined = list(
      source_type = "combined_prepared",
      components = c("HongKongBudgetsA", "HongKongBudgetsB"),
      group_variable = "housing_category",
      notes = c(
        "Combined Hong Kong household expenditure budgets, categories A and B, if both are reproducibly available.",
        "The grouped analyses remain the primary ones because categories A and B may have different distributions."
      )
    )
  )
}

slugify_dataset_name <- function(name) {
  gsub("[^a-z0-9]+", "_", tolower(name))
}

logistic_gaussian_screening_directories <- function(base_dir = file.path("output", "real_data", "logistic_gaussian", "screening")) {
  list(
    base = base_dir,
    plots = file.path(base_dir, "plots"),
    metadata = file.path(base_dir, "metadata")
  )
}

ensure_logistic_gaussian_screening_directories <- function(base_dir = file.path("output", "real_data", "logistic_gaussian", "screening")) {
  dirs <- logistic_gaussian_screening_directories(base_dir = base_dir)
  for (path in unname(unlist(dirs))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}

# Additional helper functions for data loading and preparation
make_not_found_data_prep <- function(name, notes = character(0), searched = character(0)) {
  list(
    dataset_name = name,
    X_raw = NULL,
    X_comp = matrix(numeric(0), nrow = 0L, ncol = 0L),
    X_closed = matrix(numeric(0), nrow = 0L, ncol = 0L),
    n = NA_integer_,
    D = NA_integer_,
    has_zeros = NA,
    has_missing = NA,
    n_missing_rows_removed = NA_integer_,
    min_entry = NA_real_,
    row_sums_before_closure = numeric(0),
    component_names = character(0),
    notes = c(notes, sprintf("Dataset was not found reproducibly. Searched: %s", paste(searched, collapse = ", "))),
    status = "not_found"
  )
}

load_named_dataset_from_package <- function(dataset_name, package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(NULL)
  }
  env <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    {
      utils::data(list = dataset_name, package = package, envir = env)
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  if (!isTRUE(loaded) || !exists(dataset_name, envir = env, inherits = FALSE)) {
    return(NULL)
  }
  get(dataset_name, envir = env, inherits = FALSE)
}

load_first_available_dataset <- function(candidate_names, candidate_packages) {
  searched <- character(0)
  for (pkg in candidate_packages) {
    for (nm in candidate_names) {
      searched <- c(searched, sprintf("%s::%s", pkg, nm))
      obj <- load_named_dataset_from_package(nm, pkg)
      if (!is.null(obj)) {
        return(list(object = obj, name = nm, package = pkg, searched = searched))
      }
    }
  }
  list(object = NULL, name = NA_character_, package = NA_character_, searched = searched)
}

infer_numeric_compositional_columns <- function(x_raw_df,
                                                expected_D = NULL,
                                                excluded_name_patterns = c("depth", "total", "count", "id", "site", "group", "class", "type", "category", "sex", "age", "bmi")) {
  numeric_columns <- names(x_raw_df)[vapply(x_raw_df, is.numeric, logical(1))]
  if (length(numeric_columns) == 0L) {
    return(character(0))
  }

  keep <- rep(TRUE, length(numeric_columns))
  lower_names <- tolower(numeric_columns)
  for (pattern in excluded_name_patterns) {
    keep <- keep & !grepl(pattern, lower_names, fixed = TRUE)
  }
  candidate_columns <- numeric_columns[keep]

  if (!is.null(expected_D) && length(candidate_columns) >= expected_D) {
    return(candidate_columns[seq_len(expected_D)])
  }

  if (length(candidate_columns) >= 2L) {
    return(candidate_columns)
  }

  numeric_columns
}

apply_zero_replacement <- function(x, zero_replacement = NULL) {
  x <- as.matrix(x)
  if (is.null(zero_replacement) || !any(x == 0, na.rm = TRUE)) {
    return(list(x = x, n_zeros_replaced = 0L, zero_replacement = NA_real_))
  }
  if (!is.finite(zero_replacement) || zero_replacement <= 0) {
    stop("`zero_replacement` must be a positive finite number when zeros are present.")
  }
  n_zeros <- sum(x == 0, na.rm = TRUE)
  x[x == 0] <- zero_replacement
  list(x = x, n_zeros_replaced = n_zeros, zero_replacement = zero_replacement)
}

collapse_to_top_taxa_plus_other <- function(count_or_abundance_matrix, top_taxa = 4L) {
  x <- as.matrix(count_or_abundance_matrix)
  storage.mode(x) <- "double"
  if (ncol(x) <= top_taxa + 1L) {
    return(x)
  }
  taxa_totals <- colSums(x, na.rm = TRUE)
  top_idx <- order(taxa_totals, decreasing = TRUE)[seq_len(as.integer(top_taxa))]
  other_idx <- setdiff(seq_len(ncol(x)), top_idx)
  out <- cbind(x[, top_idx, drop = FALSE], Other = rowSums(x[, other_idx, drop = FALSE], na.rm = TRUE))
  colnames(out)[seq_along(top_idx)] <- colnames(x)[top_idx]
  out
}

prepare_curated_metagenomic_placeholder <- function(name, spec) {
  notes <- c(
    spec$notes,
    "Automatic retrieval from curatedMetagenomicData is intentionally conservative in this utility file.",
    "If this row is marked not_found, implement or run the dataset-specific retrieval in the runner once the local Bioconductor objects are available."
  )
  make_not_found_data_prep(
    name = name,
    notes = notes,
    searched = c("curatedMetagenomicData", spec$study %||% NA_character_)
  )
}

close_composition_rows <- function(x) {
  row_sums <- rowSums(x)
  if (any(!is.finite(row_sums)) || any(row_sums <= 0)) {
    stop("Every compositional row must have strictly positive finite sum before closure.")
  }
  x / row_sums
}

coerce_ilr_matrix <- function(z) {
  z <- as.matrix(z)
  if (nrow(z) == 0L) {
    stop("The ilr matrix cannot be empty.")
  }
  if (is.null(dim(z)) || ncol(z) == 0L) {
    z <- matrix(as.numeric(z), ncol = 1L)
  }
  storage.mode(z) <- "double"
  z
}

ilr_basis_for_dimension <- function(D) {
  compositions::ilrBase(D = D)
}

ilr_transform_closed <- function(x_closed, V = NULL) {
  z <- compositions::ilr(compositions::acomp(x_closed), V = V)
  coerce_ilr_matrix(z)
}

inverse_ilr_to_closed <- function(z, V = NULL) {
  z <- coerce_ilr_matrix(z)
  x <- compositions::ilrInv(z, V = V)
  x <- as.matrix(compositions::acomp(x))
  x / rowSums(x)
}

distance_matrix_from_ilr <- function(z_data, z_centers) {
  z_data <- coerce_ilr_matrix(z_data)
  z_centers <- coerce_ilr_matrix(z_centers)
  if (ncol(z_data) != ncol(z_centers)) {
    stop("The ilr dimensions of `z_data` and `z_centers` do not match.")
  }

  data_sq <- rowSums(z_data^2)
  center_sq <- rowSums(z_centers^2)
  sq_distances <- outer(data_sq, center_sq, FUN = "+") - 2 * (z_data %*% t(z_centers))
  sqrt(pmax(sq_distances, 0))
}

compute_profile_matrix_from_distances <- function(distance_matrix, t_grid) {
  distance_matrix <- as.matrix(distance_matrix)
  t_grid <- as.numeric(t_grid)
  n_centers <- ncol(distance_matrix)
  output <- vapply(t_grid, function(t_value) {
    colMeans(distance_matrix <= t_value)
  }, numeric(n_centers))
  matrix(output, nrow = n_centers, ncol = length(t_grid))
}

ensure_logistic_gaussian_model_spec_available <- function() {
  if (!exists("make_logistic_gaussian_spec", mode = "function")) {
    candidate_paths <- c(
      file.path("bootstrap", "model_specs.R"),
      file.path("..", "bootstrap", "model_specs.R"),
      file.path("model_specs.R")
    )

    sourced <- FALSE
    for (path in candidate_paths) {
      if (file.exists(path)) {
        source(path)
        sourced <- TRUE
        break
      }
    }

    if (!isTRUE(sourced) || !exists("make_logistic_gaussian_spec", mode = "function")) {
      stop("Could not find make_logistic_gaussian_spec(). Source bootstrap/model_specs.R before running the logistic Gaussian analysis.")
    }
  }

  invisible(TRUE)
}

ensure_logistic_gaussian_bootstrap_available <- function() {
  ensure_logistic_gaussian_model_spec_available()
  if (!exists("multiplier_bootstrap_logistic_gaussian", mode = "function")) {
    candidate_paths <- c(
      file.path("bootstrap", "multiplier_bootstrap.R"),
      file.path("..", "bootstrap", "multiplier_bootstrap.R"),
      file.path("multiplier_bootstrap.R")
    )

    sourced <- FALSE
    for (path in candidate_paths) {
      if (file.exists(path)) {
        source(path)
        sourced <- TRUE
        break
      }
    }

    if (!isTRUE(sourced) || !exists("multiplier_bootstrap_logistic_gaussian", mode = "function")) {
      stop("Could not find multiplier_bootstrap_logistic_gaussian(). Source bootstrap/multiplier_bootstrap.R before running the logistic Gaussian analysis.")
    }
  }

  invisible(TRUE)
}

make_fitted_logistic_gaussian_theta <- function(fit, ambient_dim) {
  list(
    mu_ilr = fit$mu_hat,
    Sigma_ilr = fit$Sigma_hat,
    ambient_dim = ambient_dim
  )
}

evaluate_fitted_logistic_gaussian_profile <- function(omega_grid,
                                                      t_grid,
                                                      fit,
                                                      control = list()) {
  ensure_logistic_gaussian_model_spec_available()
  spec <- make_logistic_gaussian_spec(unknown_param = "both")

  if (is.null(spec$profile_matrix_eval) || !is.function(spec$profile_matrix_eval)) {
    stop("make_logistic_gaussian_spec() must return a callable `profile_matrix_eval` component.")
  }

  theta <- make_fitted_logistic_gaussian_theta(
    fit = fit,
    ambient_dim = ncol(as.matrix(omega_grid))
  )

  spec$profile_matrix_eval(
    omega_grid = omega_grid,
    t_grid = t_grid,
    theta = theta,
    control = control
  )
}

compute_screening_statistics_from_profiles <- function(empirical_profile, theoretical_profile, n_obs) {
  process_matrix <- sqrt(n_obs) * (empirical_profile - theoretical_profile)
  list(
    process_matrix = process_matrix,
    ks = max(abs(process_matrix)),
    cvm = mean(process_matrix^2)
  )
}

standardize_screening_bootstrap_mode <- function(bootstrap_mode) {
  bootstrap_mode <- as.character(bootstrap_mode)
  if (bootstrap_mode %in% c("composite_multiplier", "composite_parametric")) {
    return("composite_multiplier")
  }
  if (identical(bootstrap_mode, "plugin_simple_null")) {
    stop("bootstrap_mode = 'plugin_simple_null' is no longer supported by the screening workflow. Use bootstrap_mode = 'composite_multiplier'.")
  }
  stop(sprintf("Unsupported bootstrap_mode: %s", bootstrap_mode))
}

find_named_numeric_value <- function(x, candidate_names) {
  if (is.null(x)) {
    return(NA_real_)
  }
  if (is.list(x)) {
    for (candidate_name in candidate_names) {
      if (!is.null(x[[candidate_name]]) && length(x[[candidate_name]]) == 1L && is.numeric(x[[candidate_name]])) {
        return(as.numeric(x[[candidate_name]]))
      }
    }
    for (element in x) {
      value <- find_named_numeric_value(element, candidate_names)
      if (is.finite(value)) {
        return(value)
      }
    }
  }
  NA_real_
}

find_named_numeric_vector <- function(x, candidate_names) {
  if (is.null(x)) {
    return(numeric(0))
  }
  if (is.list(x)) {
    for (candidate_name in candidate_names) {
      if (!is.null(x[[candidate_name]]) && is.numeric(x[[candidate_name]]) && length(x[[candidate_name]]) > 1L) {
        return(as.numeric(x[[candidate_name]]))
      }
    }
    for (element in x) {
      value <- find_named_numeric_vector(element, candidate_names)
      if (length(value) > 0L) {
        return(value)
      }
    }
  }
  numeric(0)
}

extract_metric_result <- function(bootstrap_result, statistic_name) {
  statistic_name <- as.character(statistic_name)
  statistic_node <- NULL

  for (container_name in c("inference", "statistics", "statistic", "results", "observed")) {
    if (!is.null(bootstrap_result[[container_name]][[statistic_name]])) {
      statistic_node <- bootstrap_result[[container_name]][[statistic_name]]
      break
    }
  }
  if (is.null(statistic_node) && !is.null(bootstrap_result[[statistic_name]])) {
    statistic_node <- bootstrap_result[[statistic_name]]
  }
  if (is.null(statistic_node)) {
    statistic_node <- bootstrap_result
  }

  observed_statistic <- find_named_numeric_value(
    statistic_node,
    c("statistic", "observed_statistic", "observed", paste0(statistic_name, "_statistic"), paste0("observed_", statistic_name))
  )
  p_value <- find_named_numeric_value(
    statistic_node,
    c("p_value", "pvalue", "p.value", "p", paste0(statistic_name, "_p_value"), paste0(statistic_name, "_pvalue"))
  )
  bootstrap_statistics <- find_named_numeric_vector(
    statistic_node,
    c("bootstrap_statistics", "bootstrap_statistic", "statistics_star", "statistic_star", "star_statistics", "values", paste0(statistic_name, "_bootstrap_statistics"))
  )

  if (!is.finite(observed_statistic)) {
    observed_statistic <- find_named_numeric_value(
      bootstrap_result,
      c(paste0(statistic_name, "_statistic"), paste0("observed_", statistic_name))
    )
  }
  if (!is.finite(p_value)) {
    p_value <- find_named_numeric_value(
      bootstrap_result,
      c(paste0(statistic_name, "_p_value"), paste0(statistic_name, "_pvalue"))
    )
  }
  if (length(bootstrap_statistics) == 0L) {
    bootstrap_statistics <- find_named_numeric_vector(
      bootstrap_result,
      c(paste0(statistic_name, "_bootstrap_statistics"), paste0(statistic_name, "_statistics"))
    )
  }

  list(
    statistic = observed_statistic,
    p_value = p_value,
    bootstrap_statistics = bootstrap_statistics
  )
}

safe_matrix_cov <- function(z) {
  sigma_hat <- stats::cov(z)
  sigma_hat <- as.matrix(sigma_hat)
  if (is.null(dim(sigma_hat))) {
    sigma_hat <- matrix(as.numeric(sigma_hat), nrow = 1L, ncol = 1L)
  }
  sigma_hat
}

fit_logistic_gaussian_plugin <- function(x_closed,
                                         ridge = 1e-8,
                                         eigen_tol = 1e-10) {
  x_closed <- as.matrix(x_closed)
  if (nrow(x_closed) < 2L) {
    stop("At least two observations are required to fit a Logistic Gaussian model.")
  }
  if (ncol(x_closed) < 2L) {
    stop("The simplex dimension must be at least 2.")
  }

  ilr_basis <- ilr_basis_for_dimension(ncol(x_closed))
  z <- ilr_transform_closed(x_closed, V = ilr_basis)
  mu_hat <- colMeans(z)
  sigma_raw <- safe_matrix_cov(z)
  sigma_raw <- 0.5 * (sigma_raw + t(sigma_raw))

  eigen_raw <- eigen(sigma_raw, symmetric = TRUE, only.values = TRUE)$values
  min_eigen_raw <- min(eigen_raw)
  ridge_added <- 0
  sigma_used <- sigma_raw

  if (!is.finite(min_eigen_raw)) {
    stop("The sample covariance has non-finite eigenvalues.")
  }

  if (min_eigen_raw <= eigen_tol) {
    ridge_added <- max(as.numeric(ridge), eigen_tol - min_eigen_raw + as.numeric(ridge))
    sigma_used <- sigma_raw + ridge_added * diag(ncol(sigma_raw))
  }

  eigen_used <- eigen(sigma_used, symmetric = TRUE, only.values = TRUE)$values
  positive_eigen <- eigen_used[eigen_used > eigen_tol]
  condition_number <- if (length(positive_eigen) == 0L) {
    Inf
  } else {
    max(positive_eigen) / min(positive_eigen)
  }

  list(
    mu_hat = mu_hat,
    Sigma_hat = sigma_used,
    Sigma_hat_raw = sigma_raw,
    ridge_added = ridge_added,
    eigenvalues_raw = eigen_raw,
    eigenvalues = eigen_used,
    condition_number = condition_number,
    ilr_basis = ilr_basis,
    Z = z
  )
}

simulate_fitted_logistic_gaussian <- function(n, fit, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  z_sim <- MASS::mvrnorm(n = as.integer(n), mu = fit$mu_hat, Sigma = fit$Sigma_hat)
  z_sim <- coerce_ilr_matrix(z_sim)
  x_sim <- inverse_ilr_to_closed(z_sim, V = fit$ilr_basis)
  if (any(!is.finite(x_sim))) {
    stop("The simulated Logistic Gaussian sample contains non-finite entries.")
  }
  if (any(x_sim <= 0)) {
    stop("The simulated Logistic Gaussian sample left the simplex interior numerically.")
  }
  x_sim / rowSums(x_sim)
}

prepare_composition_dataset <- function(name) {
  registry <- logistic_gaussian_screening_dataset_registry()
  spec <- registry[[name]]
  if (is.null(spec)) {
    stop(sprintf("Unsupported dataset: %s", name))
  }

  if (identical(spec$source_type, "combined_prepared")) {
    prepared_components <- lapply(spec$components, prepare_composition_dataset)
    if (any(vapply(prepared_components, function(x) identical(x$status %||% "ok", "not_found"), logical(1)))) {
      return(make_not_found_data_prep(
        name = name,
        notes = spec$notes,
        searched = spec$components
      ))
    }
    component_names <- lapply(prepared_components, function(x) x$component_names)
    if (length(unique(vapply(component_names, paste, character(1), collapse = "||"))) != 1L) {
      stop(sprintf("Cannot combine dataset %s because component names do not match.", name))
    }
    x_comp <- do.call(rbind, lapply(prepared_components, function(x) x$X_comp))
    group <- rep(spec$components, times = vapply(prepared_components, function(x) nrow(x$X_comp), integer(1)))
    x_raw_df <- as.data.frame(x_comp, stringsAsFactors = FALSE)
    x_raw_df[[spec$group_variable %||% "group"]] <- group
    x_closed <- close_composition_rows(x_comp)
    return(list(
      dataset_name = name,
      X_raw = x_raw_df,
      X_comp = x_comp,
      X_closed = x_closed,
      n = nrow(x_closed),
      D = ncol(x_closed),
      has_zeros = any(x_comp == 0),
      has_missing = FALSE,
      n_missing_rows_removed = 0L,
      min_entry = min(x_comp),
      row_sums_before_closure = rowSums(x_comp),
      component_names = colnames(x_comp),
      group = group,
      group_variable = spec$group_variable %||% "group",
      notes = spec$notes,
      status = "ok"
    ))
  }

  if (identical(spec$source_type, "curated_metagenomic_data")) {
    return(prepare_curated_metagenomic_placeholder(name, spec))
  }

  if (identical(spec$source_type, "compositions_data")) {
    loaded <- load_first_available_dataset(
      candidate_names = spec$candidate_names,
      candidate_packages = spec$candidate_packages
    )
  } else if (identical(spec$source_type, "search_r_packages")) {
    loaded <- load_first_available_dataset(
      candidate_names = spec$candidate_names,
      candidate_packages = spec$candidate_packages
    )
  } else {
    stop(sprintf("Unsupported source_type for dataset %s: %s", name, spec$source_type))
  }

  if (is.null(loaded$object)) {
    return(make_not_found_data_prep(
      name = name,
      notes = spec$notes,
      searched = loaded$searched
    ))
  }

  x_raw <- loaded$object
  x_raw_df <- as.data.frame(x_raw, stringsAsFactors = FALSE)

  comp_columns <- spec$compositional_columns
  if (is.null(comp_columns)) {
    comp_columns <- infer_numeric_compositional_columns(
      x_raw_df,
      expected_D = spec$expected_D %||% NULL
    )
  }
  if (length(comp_columns) < 2L) {
    stop(sprintf("Dataset %s did not yield at least two compositional columns.", name))
  }

  if (is.numeric(comp_columns)) {
    x_comp_df <- x_raw_df[, comp_columns, drop = FALSE]
  } else {
    missing_columns <- setdiff(comp_columns, names(x_raw_df))
    if (length(missing_columns) > 0L) {
      stop(sprintf("Dataset %s is missing expected compositional columns: %s", name, paste(missing_columns, collapse = ", ")))
    }
    x_comp_df <- x_raw_df[, comp_columns, drop = FALSE]
  }

  x_comp_numeric <- data.matrix(x_comp_df)
  if (!is.matrix(x_comp_numeric) || nrow(x_comp_numeric) == 0L || ncol(x_comp_numeric) < 2L) {
    stop(sprintf("Dataset %s did not yield a valid compositional matrix.", name))
  }

  has_missing <- anyNA(x_comp_numeric)
  n_missing_rows <- sum(!stats::complete.cases(x_comp_numeric))
  x_comp_complete <- x_comp_numeric[stats::complete.cases(x_comp_numeric), , drop = FALSE]

  if (nrow(x_comp_complete) == 0L) {
    stop(sprintf("Dataset %s has no complete compositional rows after removing missing values.", name))
  }

  has_negative <- any(x_comp_complete < 0)
  if (isTRUE(has_negative)) {
    stop(sprintf("Dataset %s contains negative compositional entries.", name))
  }

  zero_info <- apply_zero_replacement(
    x = x_comp_complete,
    zero_replacement = spec$zero_replacement %||% NULL
  )
  x_comp_complete <- zero_info$x

  has_zeros <- any(x_comp_complete == 0)
  row_sums_before_closure <- rowSums(x_comp_complete)
  x_closed <- close_composition_rows(x_comp_complete)

  closure_error <- max(abs(rowSums(x_closed) - 1))
  if (closure_error > 1e-10) {
    stop(sprintf("Closure failed numerically for dataset %s.", name))
  }

  notes <- c(
    spec$notes,
    sprintf("Loaded as %s::%s.", loaded$package, loaded$name)
  )
  if (has_missing) {
    notes <- c(
      notes,
      sprintf("Removed %d rows with missing compositional entries before closure.", n_missing_rows)
    )
  }
  if (zero_info$n_zeros_replaced > 0L) {
    notes <- c(
      notes,
      sprintf("Replaced %d zero entries by %g before closure.", zero_info$n_zeros_replaced, zero_info$zero_replacement)
    )
  }
  if (has_zeros) {
    notes <- c(
      notes,
      "Dataset still contains zeros after preprocessing, so the standard Logistic Gaussian fit is problematic."
    )
  }

  list(
    dataset_name = name,
    X_raw = x_raw_df,
    X_comp = x_comp_complete,
    X_closed = x_closed,
    n = nrow(x_closed),
    D = ncol(x_closed),
    has_zeros = has_zeros,
    has_missing = has_missing,
    n_missing_rows_removed = n_missing_rows,
    min_entry = min(x_comp_complete),
    row_sums_before_closure = row_sums_before_closure,
    component_names = colnames(x_comp_complete),
    source_package = loaded$package,
    source_dataset_name = loaded$name,
    zero_replacement = zero_info$zero_replacement,
    n_zeros_replaced = zero_info$n_zeros_replaced,
    notes = notes,
    status = "ok"
  )
}

choose_screening_centers <- function(x_closed, max_centers = 100L, seed = 123L) {
  x_closed <- as.matrix(x_closed)
  n <- nrow(x_closed)
  max_centers <- as.integer(max_centers)
  if (!is.finite(max_centers) || max_centers <= 0L) {
    stop("`max_centers` must be a strictly positive integer.")
  }

  if (n <= max_centers) {
    center_indices <- seq_len(n)
  } else {
    set.seed(seed)
    center_indices <- sort(sample.int(n, size = max_centers, replace = FALSE))
  }

  list(
    center_indices = center_indices,
    omega = x_closed[center_indices, , drop = FALSE]
  )
}

choose_screening_t_grid <- function(distance_observed,
                                    distance_null = NULL,
                                    n_t = 60L,
                                    probs = NULL) {
  if (is.null(probs)) {
    probs <- seq(0.01, 0.99, length.out = as.integer(n_t))
  }
  probs <- as.numeric(probs)
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    stop("`probs` must lie strictly inside (0, 1).")
  }

  pooled_distances <- as.numeric(distance_observed)
  if (!is.null(distance_null)) {
    pooled_distances <- c(pooled_distances, as.numeric(distance_null))
  }
  pooled_distances <- pooled_distances[is.finite(pooled_distances)]
  if (length(pooled_distances) == 0L) {
    stop("Could not build a threshold grid because no finite distances were available.")
  }

  t_grid <- as.numeric(stats::quantile(pooled_distances, probs = probs, names = FALSE, type = 8))
  t_grid <- unique(sort(t_grid))

  if (length(t_grid) < 2L) {
    max_distance <- max(pooled_distances)
    if (!is.finite(max_distance) || max_distance <= 0) {
      max_distance <- 1
    }
    t_grid <- seq(0, max_distance, length.out = max(2L, length(probs)))
  }

  t_grid
}

mardia_multivariate_normality <- function(z, tol = 1e-10) {
  z <- coerce_ilr_matrix(z)
  n <- nrow(z)
  p <- ncol(z)

  centered <- sweep(z, 2L, colMeans(z), FUN = "-")
  sigma_hat <- safe_matrix_cov(z)
  sigma_hat <- 0.5 * (sigma_hat + t(sigma_hat))
  eigen_sigma <- eigen(sigma_hat, symmetric = TRUE)
  adjusted_values <- pmax(eigen_sigma$values, tol)
  sigma_inv <- eigen_sigma$vectors %*% diag(1 / adjusted_values, nrow = p) %*% t(eigen_sigma$vectors)
  mahal_matrix <- centered %*% sigma_inv %*% t(centered)
  mahal_diag <- diag(mahal_matrix)

  skewness <- sum(mahal_matrix^3) / (n^2)
  kurtosis <- mean(mahal_diag^2)
  skew_stat <- n * skewness / 6
  skew_df <- p * (p + 1) * (p + 2) / 6
  skew_p_value <- stats::pchisq(skew_stat, df = skew_df, lower.tail = FALSE)

  kurt_mean <- p * (p + 2)
  kurt_sd <- sqrt(8 * p * (p + 2) / n)
  kurt_z <- (kurtosis - kurt_mean) / kurt_sd
  kurt_p_value <- 2 * stats::pnorm(-abs(kurt_z))

  list(
    skewness = skewness,
    skewness_statistic = skew_stat,
    skewness_df = skew_df,
    skewness_p_value = skew_p_value,
    kurtosis = kurtosis,
    kurtosis_z = kurt_z,
    kurtosis_p_value = kurt_p_value
  )
}

compute_marginal_shapiro_pvalues <- function(z) {
  z <- coerce_ilr_matrix(z)
  if (nrow(z) < 3L || nrow(z) > 5000L) {
    return(rep(NA_real_, ncol(z)))
  }

  vapply(seq_len(ncol(z)), function(j) {
    stats::shapiro.test(z[, j])$p.value
  }, numeric(1))
}

build_qq_plot_data <- function(z, mu_hat, Sigma_hat) {
  z <- coerce_ilr_matrix(z)
  n <- nrow(z)
  p <- ncol(z)
  probs <- (seq_len(n) - 0.5) / n

  do.call(rbind, lapply(seq_len(p), function(j) {
    sigma_j <- sqrt(max(Sigma_hat[j, j], 0))
    data.frame(
      coordinate = factor(sprintf("ilr_%d", j), levels = sprintf("ilr_%d", seq_len(p))),
      theoretical = stats::qnorm(probs, mean = mu_hat[j], sd = sigma_j),
      sample = sort(z[, j]),
      stringsAsFactors = FALSE
    )
  }))
}

save_ilr_qq_plot <- function(result, file) {
  qq_df <- build_qq_plot_data(
    z = result$fit$Z,
    mu_hat = result$fit$mu_hat,
    Sigma_hat = result$fit$Sigma_hat
  )

  plt <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
    geom_abline(slope = 1, intercept = 0, color = "#b2182b", linewidth = 0.5) +
    geom_point(color = "#2166ac", size = 1.4, alpha = 0.8) +
    facet_wrap(~ coordinate, scales = "free", ncol = 2) +
    labs(
      title = sprintf("%s: marginal ilr QQ plots", result$dataset_name),
      x = "Theoretical normal quantiles",
      y = "Sample ilr quantiles"
    ) +
    theme_minimal(base_size = 11)

  ggplot2::ggsave(filename = file, plot = plt, width = 7, height = 4 + 1.5 * ceiling(ncol(result$fit$Z) / 2), dpi = 180)
  invisible(file)
}

save_ilr_scatter_plot <- function(result, file) {
  z <- coerce_ilr_matrix(result$fit$Z)
  if (ncol(z) != 2L) {
    return(invisible(NULL))
  }

  df <- data.frame(z1 = z[, 1], z2 = z[, 2])
  plt <- ggplot(df, aes(x = z1, y = z2)) +
    geom_point(color = "#1b7837", size = 2, alpha = 0.8) +
    labs(
      title = sprintf("%s: ilr scatter", result$dataset_name),
      x = "ilr_1",
      y = "ilr_2"
    ) +
    theme_minimal(base_size = 11)

  ggplot2::ggsave(filename = file, plot = plt, width = 5.5, height = 4.5, dpi = 180)
  invisible(file)
}

project_simplex_to_ternary <- function(x_closed) {
  x_closed <- as.matrix(x_closed)
  if (ncol(x_closed) != 3L) {
    stop("A ternary projection is only available for 3-part compositions.")
  }

  data.frame(
    x = x_closed[, 2] + 0.5 * x_closed[, 3],
    y = sqrt(3) * 0.5 * x_closed[, 3]
  )
}

save_ternary_plot <- function(result, file) {
  x_closed <- as.matrix(result$data_prep$X_closed)
  if (ncol(x_closed) != 3L) {
    return(invisible(NULL))
  }

  proj <- project_simplex_to_ternary(x_closed)
  component_names <- result$data_prep$component_names
  triangle <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3) / 2, 0)
  )
  labels <- data.frame(
    x = c(-0.03, 1.03, 0.5),
    y = c(-0.03, -0.03, sqrt(3) / 2 + 0.03),
    label = component_names
  )

  plt <- ggplot() +
    geom_path(data = triangle, aes(x = x, y = y), color = "grey30", linewidth = 0.6) +
    geom_point(data = proj, aes(x = x, y = y), color = "#762a83", size = 2, alpha = 0.85) +
    geom_text(data = labels, aes(x = x, y = y, label = label), size = 3.8) +
    coord_equal() +
    labs(
      title = sprintf("%s: ternary plot", result$dataset_name),
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 11)

  ggplot2::ggsave(filename = file, plot = plt, width = 5.5, height = 4.8, dpi = 180)
  invisible(file)
}

save_screening_plots <- function(result, base_dir) {
  dirs <- ensure_logistic_gaussian_screening_directories(base_dir)
  slug <- slugify_dataset_name(result$dataset_name)
  plot_paths <- list(
    ilr_qq = file.path(dirs$plots, sprintf("%s_ilr_qq.png", slug)),
    ilr_scatter = file.path(dirs$plots, sprintf("%s_ilr_scatter.png", slug)),
    ternary = file.path(dirs$plots, sprintf("%s_ternary.png", slug))
  )

  save_ilr_qq_plot(result, plot_paths$ilr_qq)
  save_ilr_scatter_plot(result, plot_paths$ilr_scatter)
  save_ternary_plot(result, plot_paths$ternary)

  plot_paths
}

screening_problem_reasons <- function(data_prep, fit, ilr_dim) {
  reasons <- character(0)

  if (isTRUE(data_prep$has_zeros)) {
    reasons <- c(reasons, "zeros in the compositional sample")
  }
  if (isTRUE(data_prep$has_missing) && data_prep$n_missing_rows_removed > 0.1 * nrow(data_prep$X_raw)) {
    reasons <- c(reasons, "more than 10% of rows removed because of missing values")
  }
  if (data_prep$n <= ilr_dim) {
    reasons <- c(reasons, "sample size is not larger than the ilr dimension")
  }
  if (!is.finite(fit$condition_number) || fit$condition_number > 1e12) {
    reasons <- c(reasons, "ill-conditioned fitted covariance")
  }
  if (fit$ridge_added > 1e-4) {
    reasons <- c(reasons, "substantial ridge regularization was required")
  }

  reasons
}

classify_logistic_gaussian_screening <- function(result) {
  fit <- result$fit
  data_prep <- result$data_prep
  diagnostics <- result$diagnostics
  p_min <- min(result$inference$ks$p_value, result$inference$cvm$p_value)
  problem_reasons <- screening_problem_reasons(data_prep, fit, ilr_dim = ncol(fit$Z))

  mardia_p_min <- min(
    diagnostics$mardia$skewness_p_value,
    diagnostics$mardia$kurtosis_p_value,
    na.rm = TRUE
  )
  shapiro_p_min <- min(diagnostics$shapiro_p_values, na.rm = TRUE)
  if (!is.finite(shapiro_p_min)) {
    shapiro_p_min <- NA_real_
  }

  diagnostics_reasonable <- isTRUE(mardia_p_min > 0.01) &&
    (is.na(shapiro_p_min) || shapiro_p_min > 0.01)

  if (length(problem_reasons) > 0L) {
    return(list(
      diagnosis = "problematic",
      use_in_paper = "no",
      why = paste(problem_reasons, collapse = "; ")
    ))
  }

  if (p_min <= 0.01) {
    return(list(
      diagnosis = "reject",
      use_in_paper = "no",
      why = "At least one logistic Gaussian GOF p-value is at most 0.01."
    ))
  }

  if (p_min > 0.10 && diagnostics_reasonable) {
    why <- "Logistic Gaussian GOF p-values exceed 0.10 and the ilr diagnostics are not clearly contradictory."
    if (data_prep$n < 30L) {
      why <- paste(why, "The sample is small, so lack of power remains a serious caveat.")
    }
    return(list(
      diagnosis = "promising",
      use_in_paper = "yes",
      why = why
    ))
  }

  list(
    diagnosis = "borderline",
    use_in_paper = "no",
    why = "Evidence is mixed: the GOF p-values and/or ilr diagnostics are not fully convincing."
  )
}

compute_seed_sensitivity_summary <- function(reference_result,
                                             alt_result) {
  list(
    assessed = TRUE,
    reference_seed = reference_result$settings$seed,
    alternative_seed = alt_result$settings$seed,
    reference_B = reference_result$settings$B,
    alternative_B = alt_result$settings$B,
    delta_ks_p_value = alt_result$inference$ks$p_value - reference_result$inference$ks$p_value,
    delta_cvm_p_value = alt_result$inference$cvm$p_value - reference_result$inference$cvm$p_value
  )
}

make_problematic_screening_result <- function(dataset_name,
                                              data_prep,
                                              reason,
                                              started_at,
                                              dataset_started_proc,
                                              B,
                                              max_centers,
                                              n_t,
                                              seed,
                                              alpha,
                                              ridge,
                                              bootstrap_mode,
                                              n_cores,
                                              output_dir,
                                              make_plots,
                                              save_outputs) {
  result <- list(
    dataset_name = dataset_name,
    screening_type = "logistic_gaussian_composite_or_screening",
    warning = logistic_gaussian_screening_warning,
    data_prep = data_prep,
    fit = NULL,
    grid = list(
      center_indices = integer(0),
      omega = NULL,
      t_grid = numeric(0),
      probs_t = seq(0.01, 0.99, length.out = n_t),
      max_centers = max_centers
    ),
    theoretical_profile_info = list(
      method = "not_evaluated",
      note = "No theoretical profile was evaluated because the dataset was problematic or not found.",
      theoretical_profile = NULL
    ),
    observed = list(
      empirical_profile = NULL,
      theoretical_profile = NULL,
      process_matrix = NULL
    ),
    bootstrap = list(
      mode = bootstrap_mode,
      ks_statistics = numeric(0),
      cvm_statistics = numeric(0),
      bootstrap_fit_summaries = NULL
    ),
    inference = list(
      ks = list(statistic = NA_real_, p_value = NA_real_, bootstrap_replicates = B),
      cvm = list(statistic = NA_real_, p_value = NA_real_, bootstrap_replicates = B)
    ),
    diagnostics = list(
      mardia = list(
        skewness = NA_real_,
        skewness_statistic = NA_real_,
        skewness_df = NA_real_,
        skewness_p_value = NA_real_,
        kurtosis = NA_real_,
        kurtosis_z = NA_real_,
        kurtosis_p_value = NA_real_
      ),
      shapiro_p_values = NA_real_
    ),
    classification = list(
      diagnosis = if (identical(data_prep$status %||% "ok", "not_found")) "not_found" else "problematic",
      use_in_paper = "no",
      why = reason
    ),
    settings = list(
      B = B,
      max_centers = max_centers,
      n_t = n_t,
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      bootstrap_mode = bootstrap_mode,
      n_cores = n_cores
    ),
    runtime = list(
      started_at = started_at,
      ended_at = Sys.time(),
      elapsed_seconds = unname(proc.time()[["elapsed"]] - dataset_started_proc)
    ),
    sensitivity = list(assessed = FALSE)
  )

  if (isTRUE(save_outputs)) {
    dirs <- ensure_logistic_gaussian_screening_directories(output_dir)
    slug <- slugify_dataset_name(dataset_name)
    result_path <- file.path(dirs$base, sprintf("%s_results.rds", slug))
    saveRDS(result, file = result_path)
    result$output_paths <- list(result_rds = result_path)
  }

  result
}

run_logistic_gaussian_screening <- function(dataset_name,
                                            B = 1000L,
                                            max_centers = 100L,
                                            n_t = 60L,
                                            probs_t = NULL,
                                            bootstrap_mode = c("composite_multiplier", "composite_parametric", "plugin_simple_null"),
                                            seed = 123L,
                                            alpha = 0.05,
                                            ridge = 1e-8,
                                            n_cores = 1L,
                                            control = list(),
                                            make_plots = TRUE,
                                            save_outputs = TRUE,
                                            output_dir = file.path("output", "calibration", "bootstrap", "logistic_gaussian", "composite"),
                                            run_seed_sensitivity = FALSE,
                                            sensitivity_B = NULL,
                                            verbose = TRUE) {
  bootstrap_mode <- standardize_screening_bootstrap_mode(match.arg(bootstrap_mode))
  started_at <- Sys.time()
  dataset_started_proc <- proc.time()[["elapsed"]]

  if (verbose) {
    message(sprintf("[%s] screening %s", format(started_at, "%Y-%m-%d %H:%M:%S"), dataset_name))
  }

  data_prep <- prepare_composition_dataset(dataset_name)
  if (identical(data_prep$status %||% "ok", "not_found")) {
    return(make_problematic_screening_result(
      dataset_name = dataset_name,
      data_prep = data_prep,
      reason = "Dataset was not found reproducibly and was therefore not analysed.",
      started_at = started_at,
      dataset_started_proc = dataset_started_proc,
      B = B,
      max_centers = max_centers,
      n_t = n_t,
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      bootstrap_mode = bootstrap_mode,
      n_cores = n_cores,
      output_dir = output_dir,
      make_plots = make_plots,
      save_outputs = save_outputs
    ))
  }
  notes <- data_prep$notes
  if (isTRUE(data_prep$has_zeros)) {
    notes <- c(notes, "Main screening fit skipped because zero replacement is intentionally not automatic.")
    data_prep$notes <- notes
    return(make_problematic_screening_result(
      dataset_name = dataset_name,
      data_prep = data_prep,
      reason = "The dataset contains zeros, so the standard Logistic Gaussian screening fit was not attempted.",
      started_at = started_at,
      dataset_started_proc = dataset_started_proc,
      B = B,
      max_centers = max_centers,
      n_t = n_t,
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      bootstrap_mode = bootstrap_mode,
      n_cores = n_cores,
      output_dir = output_dir,
      make_plots = make_plots,
      save_outputs = save_outputs
    ))
  }

  fit <- fit_logistic_gaussian_plugin(
    x_closed = data_prep$X_closed,
    ridge = ridge
  )

  centers <- choose_screening_centers(
    x_closed = data_prep$X_closed,
    max_centers = max_centers,
    seed = seed
  )

  z_centers <- ilr_transform_closed(centers$omega, V = fit$ilr_basis)
  distance_observed <- distance_matrix_from_ilr(fit$Z, z_centers)

  t_grid <- choose_screening_t_grid(
    distance_observed = distance_observed,
    distance_null = NULL,
    n_t = n_t,
    probs = probs_t
  )

  theoretical_profile <- evaluate_fitted_logistic_gaussian_profile(
    omega_grid = centers$omega,
    t_grid = t_grid,
    fit = fit
  )
  empirical_profile <- compute_profile_matrix_from_distances(distance_observed, t_grid)
  empirical_profile <- compute_profile_matrix_from_distances(distance_observed, t_grid)
  observed_stats_grid <- compute_screening_statistics_from_profiles(
    empirical_profile = empirical_profile,
    theoretical_profile = theoretical_profile,
    n_obs = data_prep$n
  )

  ensure_logistic_gaussian_bootstrap_available()
  bootstrap_result <- multiplier_bootstrap_logistic_gaussian(
    data = data_prep$X_closed,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = list(
      omega_grid = centers$omega,
      t_grid = t_grid
    ),
    B = B,
    alpha = alpha,
    n_cores = n_cores,
    seed = seed,
    control = control,
    unknown_param = "both"
  )

  ks_extracted <- extract_metric_result(bootstrap_result, "ks")
  cvm_extracted <- extract_metric_result(bootstrap_result, "cvm")

  inference <- list(
    ks = list(
      statistic = ks_extracted$statistic,
      p_value = ks_extracted$p_value,
      bootstrap_replicates = B
    ),
    cvm = list(
      statistic = cvm_extracted$statistic,
      p_value = cvm_extracted$p_value,
      bootstrap_replicates = B
    )
  )

  diagnostics <- list(
    mardia = mardia_multivariate_normality(fit$Z),
    shapiro_p_values = compute_marginal_shapiro_pvalues(fit$Z)
  )

  diagnostics <- list(
    mardia = mardia_multivariate_normality(fit$Z),
    shapiro_p_values = compute_marginal_shapiro_pvalues(fit$Z)
  )

  result <- list(
    dataset_name = dataset_name,
    screening_type = "logistic_gaussian_composite_or_screening",
    warning = logistic_gaussian_screening_warning,
    data_prep = data_prep,
    fit = fit,
    grid = list(
      center_indices = centers$center_indices,
      omega = centers$omega,
      t_grid = t_grid,
      probs_t = probs_t %||% seq(0.01, 0.99, length.out = n_t),
      max_centers = max_centers
    ),
    theoretical_profile_info = list(
      method = "global_logistic_gaussian_model_spec",
      note = "Theoretical distance profiles are evaluated by make_logistic_gaussian_spec()$profile_matrix_eval, the same model-spec path used by the calibration code.",
      theoretical_profile = theoretical_profile
    ),
    observed = list(
      empirical_profile = empirical_profile,
      theoretical_profile = theoretical_profile,
      process_matrix = observed_stats_grid$process_matrix,
      grid_ks = observed_stats_grid$ks,
      grid_cvm = observed_stats_grid$cvm
    ),
    bootstrap = list(
      mode = bootstrap_mode,
      engine = "multiplier_bootstrap_logistic_gaussian",
      raw_result = bootstrap_result,
      ks_statistics = ks_extracted$bootstrap_statistics,
      cvm_statistics = cvm_extracted$bootstrap_statistics,
      bootstrap_fit_summaries = NULL
    ),
    inference = inference,
    diagnostics = diagnostics,
    settings = list(
      B = B,
      max_centers = max_centers,
      n_t = length(t_grid),
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      bootstrap_mode = bootstrap_mode,
      n_cores = n_cores,
      control = control
    ),
    runtime = list(
      started_at = started_at
    )
  )

  classification <- classify_logistic_gaussian_screening(result)
  result$classification <- classification

  if (isTRUE(run_seed_sensitivity)) {
    alt_seed <- as.integer(seed) + 1000L
    alt_result <- run_logistic_gaussian_screening(
      dataset_name = dataset_name,
      B = sensitivity_B %||% max(99L, ceiling(B / 4)),
      max_centers = max_centers,
      n_t = n_t,
      probs_t = probs_t,
      bootstrap_mode = bootstrap_mode,
      seed = alt_seed,
      alpha = alpha,
      ridge = ridge,
      n_cores = n_cores,
      control = control,
      make_plots = FALSE,
      save_outputs = FALSE,
      output_dir = output_dir,
      run_seed_sensitivity = FALSE,
      verbose = FALSE
    )
    result$sensitivity <- compute_seed_sensitivity_summary(result, alt_result)
  } else {
    result$sensitivity <- list(assessed = FALSE)
  }

  result$runtime$ended_at <- Sys.time()
  result$runtime$elapsed_seconds <- unname(proc.time()[["elapsed"]] - dataset_started_proc)

  if (isTRUE(save_outputs)) {
    dirs <- ensure_logistic_gaussian_screening_directories(output_dir)
    slug <- slugify_dataset_name(dataset_name)
    result_path <- file.path(dirs$base, sprintf("%s_results.rds", slug))
    saveRDS(result, file = result_path)
    result$output_paths <- list(result_rds = result_path)

    if (isTRUE(make_plots)) {
      result$output_paths$plots <- save_screening_plots(result, base_dir = output_dir)
      saveRDS(result, file = result_path)
    }
  }

  result
}

make_logistic_gaussian_screening_summary_row <- function(result) {
  fit <- result$fit
  mardia <- result$diagnostics$mardia
  shapiro_min_pvalue <- suppressWarnings(min(result$diagnostics$shapiro_p_values, na.rm = TRUE))
  if (!is.finite(shapiro_min_pvalue)) {
    shapiro_min_pvalue <- NA_real_
  }

  data.frame(
    dataset = result$dataset_name,
    status = result$data_prep$status %||% "ok",
    source_package = result$data_prep$source_package %||% NA_character_,
    source_dataset_name = result$data_prep$source_dataset_name %||% NA_character_,
    n = result$data_prep$n,
    D = result$data_prep$D,
    zeros = if (isTRUE(result$data_prep$has_zeros)) "yes" else "no",
    missing = if (isTRUE(result$data_prep$has_missing)) "yes" else "no",
    missing_rows_removed = result$data_prep$n_missing_rows_removed,
    zero_replacement = result$data_prep$zero_replacement %||% NA_real_,
    n_zeros_replaced = result$data_prep$n_zeros_replaced %||% NA_integer_,
    min_component = result$data_prep$min_entry,
    ridge_added = if (is.null(fit)) NA_real_ else fit$ridge_added,
    min_eigenvalue = if (is.null(fit)) NA_real_ else min(fit$eigenvalues),
    max_eigenvalue = if (is.null(fit)) NA_real_ else max(fit$eigenvalues),
    condition_number = if (is.null(fit)) NA_real_ else fit$condition_number,
    ks_statistic = result$inference$ks$statistic,
    ks_pvalue = result$inference$ks$p_value,
    cvm_statistic = result$inference$cvm$statistic,
    cvm_pvalue = result$inference$cvm$p_value,
    mardia_skew_pvalue = mardia$skewness_p_value,
    mardia_kurtosis_pvalue = mardia$kurtosis_p_value,
    shapiro_min_pvalue = shapiro_min_pvalue,
    diagnosis = result$classification$diagnosis,
    use_in_paper = result$classification$use_in_paper,
    why = result$classification$why,
    bootstrap_mode = result$bootstrap$mode,
    bootstrap_engine = result$bootstrap$engine %||% NA_character_,
    B = result$settings$B,
    seed = result$settings$seed,
    elapsed_seconds = result$runtime$elapsed_seconds,
    notes = paste(result$data_prep$notes %||% character(0), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

run_logistic_gaussian_screening_batch <- function(dataset_names = default_logistic_gaussian_screening_datasets(),
                                                  B = 1000L,
                                                  max_centers = 100L,
                                                  n_t = 60L,
                                                  probs_t = NULL,
                                                  bootstrap_mode = "composite_multiplier",
                                                  seed = 123L,
                                                  alpha = 0.05,
                                                  ridge = 1e-8,
                                                  n_cores = 1L,
                                                  control = list(),
                                                  make_plots = TRUE,
                                                  output_dir = file.path("output", "calibration", "bootstrap", "logistic_gaussian", "composite"),
                                                  run_seed_sensitivity = FALSE,
                                                  verbose = TRUE) {
  dirs <- ensure_logistic_gaussian_screening_directories(output_dir)
  package_versions <- vapply(
    c("compositions", "MASS", "ggplot2"),
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  )

  summary_rows <- vector("list", length(dataset_names))
  names(summary_rows) <- dataset_names
  result_paths <- character(length(dataset_names))
  names(result_paths) <- dataset_names

  for (i in seq_along(dataset_names)) {
    dataset_name <- dataset_names[[i]]
    dataset_seed <- as.integer(seed) + 100L * (i - 1L)
    result <- tryCatch(
      run_logistic_gaussian_screening(
        dataset_name = dataset_name,
        B = B,
        max_centers = max_centers,
        n_t = n_t,
        probs_t = probs_t,
        bootstrap_mode = bootstrap_mode,
        seed = dataset_seed,
        alpha = alpha,
        ridge = ridge,
        n_cores = n_cores,
        control = control,
        make_plots = make_plots,
        save_outputs = TRUE,
        output_dir = output_dir,
        run_seed_sensitivity = run_seed_sensitivity,
        verbose = verbose
      ),
      error = function(e) {
        list(
          dataset_name = dataset_name,
          screening_type = "logistic_gaussian_composite_or_screening",
          warning = logistic_gaussian_screening_warning,
          data_prep = list(
            n = NA_integer_,
            D = NA_integer_,
            has_zeros = NA,
            has_missing = NA,
            n_missing_rows_removed = NA_integer_,
            min_entry = NA_real_
          ),
          fit = NULL,
          inference = list(
            ks = list(statistic = NA_real_, p_value = NA_real_),
            cvm = list(statistic = NA_real_, p_value = NA_real_)
          ),
          diagnostics = list(
            mardia = list(skewness_p_value = NA_real_, kurtosis_p_value = NA_real_),
            shapiro_p_values = NA_real_
          ),
          classification = list(
            diagnosis = "problematic",
            use_in_paper = "no",
            why = conditionMessage(e)
          ),
          bootstrap = list(mode = bootstrap_mode, engine = NA_character_),
          settings = list(B = B, seed = dataset_seed, n_cores = n_cores, control = control),
          runtime = list(elapsed_seconds = NA_real_)
        )
      }
    )
    summary_rows[[i]] <- make_logistic_gaussian_screening_summary_row(result)
    result_paths[[i]] <- result$output_paths$result_rds %||% NA_character_
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_csv <- file.path(dirs$base, "summary_logistic_gaussian_screening.csv")
  utils::write.csv(summary_df, file = summary_csv, row.names = FALSE)

  warning_file <- file.path(dirs$metadata, "screening_warning.txt")
  writeLines(logistic_gaussian_screening_warning, con = warning_file)

  session_info_file <- file.path(dirs$metadata, "session_info.txt")
  writeLines(capture.output(sessionInfo()), con = session_info_file)

  package_versions_file <- file.path(dirs$metadata, "package_versions.csv")
  utils::write.csv(
    data.frame(package = names(package_versions), version = unname(package_versions), stringsAsFactors = FALSE),
    file = package_versions_file,
    row.names = FALSE
  )

  config_file <- file.path(dirs$metadata, "run_config.rds")
  saveRDS(
    list(
      dataset_names = dataset_names,
      B = B,
      max_centers = max_centers,
      n_t = n_t,
      probs_t = probs_t %||% seq(0.01, 0.99, length.out = n_t),
      bootstrap_mode = bootstrap_mode,
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      n_cores = n_cores,
      control = control,
      result_paths = result_paths
    ),
    file = config_file
  )

  list(
    summary = summary_df,
    summary_csv = summary_csv,
    session_info_file = session_info_file,
    package_versions_file = package_versions_file,
    warning_file = warning_file,
    config_file = config_file
  )
}
