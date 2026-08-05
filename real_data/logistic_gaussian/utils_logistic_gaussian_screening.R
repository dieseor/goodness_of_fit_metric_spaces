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
  "The default quadratic-form evaluator is the shared MVN/logistic-Gaussian auto dispatcher: it uses exact CompQuadForm methods with diagnostics and reserves Monte Carlo for a terminal numerical rescue.",
  "The HBE approximation from sphunif::p_wschisq() remains available only when explicitly requested for historical reproduction or auditing.",
  "For the logistic Gaussian model with Aitchison distance, theoretical distance profiles are evaluated through the global logistic Gaussian model specification used by the calibration code, not by local screening-specific profile code.",
  "The exploratory label plugin_simple_null runs a simple-null plug-in screening: the logistic Gaussian parameters are estimated once from the data, then treated as fixed, the theoretical profile is approximated by Monte Carlo, and the bootstrap samples are generated from that fitted null without re-estimation.",
  "Therefore, plugin_simple_null p-values are exploratory only and are not valid p-values for the composite null."
)

normalize_logistic_gaussian_screening_control <- function(control = list()) {
  control <- control %||% list()
  # EN DUDA (2026-07-26): generic shared dispatcher control.  The legacy
  # logistic-Gaussian name remains accepted by mvnormal_quadform.R.
  if (is.null(control$mvnormal_quadform_method) && is.null(control$logistic_gaussian_quadform_method)) {
    control$mvnormal_quadform_method <- "auto"
  }
  control
}

composition_registry <- list(
  Aar_oxides = list(
    data_entry = "Aar",
    object_name = "Aar",
    parts = c("SiO2", "TiO2", "Al2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "Fe2O3t")
  ),
  Activity10 = list(
    data_entry = "Activity10",
    object_name = "Activity10",
    parts = c("teac", "cons", "admi", "rese", "wake", "slee")
  ),
  Activity31 = list(
    data_entry = "Activity31",
    object_name = "Activity31",
    parts = c("teac", "cons", "admi", "rese", "wake", "slee")
  ),
  AnimalVegetation = list(
    data_entry = "AnimalVegetation",
    object_name = "AnimalVegetation",
    parts = c("disc", "spick", "din", "spin")
  ),
  ArcticLake = list(
    data_entry = "ArcticLake",
    object_name = "ArcticLake",
    parts = c("sand", "silt", "clay")
  ),
  Bayesite = list(
    data_entry = "Bayesite",
    object_name = "Bayesite",
    parts = c("A", "B", "C", "D")
  ),
  Boxite = list(
    data_entry = "Boxite",
    object_name = "Boxite",
    parts = c("A", "B", "C", "D", "E")
  ),
  ClamEast = list(
    data_entry = "ClamEast",
    object_name = "ClamEast",
    parts = c("dl", "dm", "ds", "ll", "lm", "ls")
  ),
  ClamWest = list(
    data_entry = "ClamWest",
    object_name = "ClamWest",
    parts = c("dl", "dm", "ds", "ll", "lm", "ls")
  ),
  Coxite = list(
    data_entry = "Coxite",
    object_name = "Coxite",
    parts = c("A", "B", "C", "D", "E")
  ),
  DiagnosticProb = list(
    data_entry = "DiagnosticProb",
    object_name = "DiagnosticProb",
    parts = c("A", "B", "C")
  ),
  Firework = list(
    data_entry = "Firework",
    object_name = "Firework",
    parts = c("a", "b", "c", "d", "e")
  ),
  Hongite = list(
    data_entry = "Hongite",
    object_name = "Hongite",
    parts = c("A", "B", "C", "D", "E")
  ),
  HouseholdExp = list(
    data_entry = "HouseholdExp",
    object_name = "HouseholdExp",
    parts = c("Housing", "Food", "Other", "Services")
  ),
  Hydrochem = list(
    data_entry = "Hydrochem",
    object_name = "Hydrochem",
    parts = c("H", "Na", "K", "Mg", "Ca", "Sr", "Ba", "NH4", "Cl", "NO3", "PO4", "SO4", "HCO3", "TOC")
  ),
  juraset = list(
    data_entry = "juraset",
    object_name = "juraset",
    parts = c("Cd", "Cu", "Pb", "Co", "Cr", "Ni", "Zn")
  ),
  Kongite = list(
    data_entry = "Kongite",
    object_name = "Kongite",
    parts = c("A", "B", "C", "D", "E")
  ),
  Metabolites = list(
    data_entry = "Metabolites",
    object_name = "Metabolites",
    parts = c("met1", "met2", "met3")
  ),
  PogoJump = list(
    data_entry = "PogoJump",
    object_name = "PogoJump",
    parts = c("yat", "yee", "sam")
  ),
  Sediments = list(
    data_entry = "Sediments",
    object_name = "Sediments",
    parts = c("sand", "silt", "clay")
  ),
  SerumProtein = list(
    data_entry = "SerumProtein",
    object_name = "SerumProtein",
    parts = c("a", "b", "c", "d")
  ),
  ShiftOperators = list(
    data_entry = "ShiftOperators",
    object_name = "ShiftOperators",
    parts = c("A", "B", "C", "D")
  ),
  SkyeAFM = list(
    data_entry = "SkyeAFM",
    object_name = "SkyeAFM",
    parts = c("A", "F", "M")
  ),
  Supervisor = list(
    data_entry = "Supervisor",
    object_name = "Supervisor",
    parts = c("C", "D", "E", "F")
  ),
  WhiteCells_microscopic = list(
    data_entry = "WhiteCells",
    object_name = "WhiteCells",
    parts = c("mG", "mL", "mM")
  ),
  WhiteCells_image = list(
    data_entry = "WhiteCells",
    object_name = "WhiteCells",
    parts = c("iG", "iL", "iM")
  ),
  Yatquat_preference = list(
    data_entry = "Yatquat",
    object_name = "Yatquat",
    parts = c("prFL", "prSK", "prST")
  ),
  Yatquat_panel = list(
    data_entry = "Yatquat",
    object_name = "Yatquat",
    parts = c("paFL", "paSK", "paST")
  )
)

stopifnot(length(composition_registry) == 28L)
stopifnot(!anyDuplicated(names(composition_registry)))

external_logistic_gaussian_screening_datasets <- function() {
  c(
    "SkyeLavas",
    "SkyeLavasAitchison32",
    "ClamCombined",
    "coffee",
    "expenditures",
    "expendituresEU",
    "alcohol"
  )
}

default_logistic_gaussian_screening_datasets <- function() {
  c(names(composition_registry), external_logistic_gaussian_screening_datasets())
}

make_skye_lavas_blocks <- function() {
  make_block <- function(sample_id, rock_type,
                         SiO2, Al2O3, Fe2O3, MgO, CaO,
                         Na2O, K2O, TiO2, P2O5, MnO) {
    data.frame(
      sample_id = as.character(sample_id),
      rock_type = as.character(rock_type),
      SiO2 = SiO2,
      Al2O3 = Al2O3,
      Fe2O3 = Fe2O3,
      MgO = MgO,
      CaO = CaO,
      Na2O = Na2O,
      K2O = K2O,
      TiO2 = TiO2,
      P2O5 = P2O5,
      MnO = MnO,
      stringsAsFactors = FALSE
    )
  }

  b1 <- make_block(
    sample_id = c(937, 976, 925, 974, 924, 950, 891, 928, 929, 290),
    rock_type = rep("Basalt", 10),
    SiO2 = c(46.31, 45.51, 45.19, 47.71, 45.68, 45.60, 46.80, 46.73, 46.83, 46.10),
    Al2O3 = c(14.18, 14.13, 13.51, 14.46, 14.02, 14.16, 14.57, 14.59, 14.90, 13.92),
    Fe2O3 = c(12.32, 12.84, 13.32, 11.49, 13.20, 13.58, 12.98, 12.00, 11.67, 12.73),
    MgO = c(12.74, 12.65, 12.27, 10.50, 11.52, 11.59, 10.99, 10.16, 9.86, 10.65),
    CaO = c(9.62, 9.41, 8.30, 9.86, 8.52, 8.88, 9.03, 9.06, 9.63, 9.87),
    Na2O = c(2.51, 2.47, 2.89, 2.48, 2.66, 2.70, 2.78, 2.53, 2.92, 2.66),
    K2O = c(0.34, 0.47, 0.37, 0.60, 0.47, 0.43, 0.40, 0.38, 0.44, 0.54),
    TiO2 = c(1.53, 1.66, 2.03, 1.30, 1.95, 1.93, 1.65, 1.60, 1.44, 1.82),
    P2O5 = c(0.16, 0.18, 0.22, 0.16, 0.22, 0.21, 0.19, 0.16, 0.17, 0.20),
    MnO = c(0.18, 0.19, 0.16, 0.18, 0.19, 0.19, 0.18, 0.17, 0.16, 0.19)
  )

  b2 <- make_block(
    sample_id = c(927, 921, 968, 977, 953, 952, 923, 289, 285, 922, 280, 926),
    rock_type = rep("Basalt", 12),
    SiO2 = c(45.39, 45.05, 45.38, 46.88, 47.09, 48.55, 46.14, 46.60, 46.67, 45.88, 45.28, 46.34),
    Al2O3 = c(14.48, 14.01, 15.11, 15.17, 14.80, 14.44, 14.77, 15.65, 15.79, 13.80, 15.64, 14.55),
    Fe2O3 = c(12.75, 13.94, 13.01, 12.33, 11.55, 10.94, 12.76, 11.79, 13.07, 13.10, 14.13, 13.65),
    MgO = c(10.55, 11.47, 10.28, 9.70, 9.09, 8.60, 9.85, 8.92, 9.57, 9.58, 10.13, 9.78),
    CaO = c(8.81, 9.04, 9.58, 9.90, 10.41, 9.77, 8.81, 10.78, 9.73, 8.73, 10.12, 9.11),
    Na2O = c(2.29, 2.70, 2.48, 2.61, 2.43, 2.52, 2.72, 2.30, 2.96, 2.96, 2.74, 2.84),
    K2O = c(0.48, 0.35, 0.45, 0.51, 0.55, 0.62, 0.53, 0.46, 0.48, 0.53, 0.25, 0.57),
    TiO2 = c(1.90, 2.15, 1.79, 1.42, 1.40, 1.24, 1.83, 1.42, 1.82, 1.99, 1.82, 1.90),
    P2O5 = c(0.21, 0.23, 0.20, 0.17, 0.18, 0.15, 0.20, 0.16, 0.24, 0.23, 0.19, 0.25),
    MnO = c(0.19, 0.20, 0.20, 0.19, 0.19, 0.18, 0.18, 0.22, 0.22, 0.17, 0.22, 0.22)
  )

  b3 <- make_block(
    sample_id = c(932, 279, 892, 949, 216, 282, 896, 276, 284, 967),
    rock_type = rep("Basalt", 10),
    SiO2 = c(44.89, 45.44, 48.01, 44.12, 44.72, 44.60, 47.05, 44.37, 46.54, 47.84),
    Al2O3 = c(16.04, 15.18, 15.59, 14.95, 16.15, 15.81, 15.28, 15.40, 16.05, 17.47),
    Fe2O3 = c(14.76, 14.54, 11.70, 14.81, 14.62, 14.41, 13.71, 14.98, 13.06, 12.97),
    MgO = c(9.97, 9.58, 7.63, 9.28, 8.42, 8.21, 7.79, 8.43, 6.98, 6.66),
    CaO = c(9.29, 8.88, 10.29, 9.31, 9.93, 9.70, 8.83, 9.61, 10.23, 10.70),
    Na2O = c(2.81, 3.35, 2.59, 2.78, 2.79, 2.80, 3.29, 2.83, 2.82, 2.83),
    K2O = c(0.22, 0.47, 0.62, 0.41, 0.35, 0.48, 0.70, 0.39, 0.61, 0.51),
    TiO2 = c(2.05, 2.14, 1.27, 2.63, 2.25, 2.30, 1.60, 2.58, 1.99, 1.25),
    P2O5 = c(0.20, 0.27, 0.17, 0.28, 0.26, 0.27, 0.26, 0.28, 0.29, 0.17),
    MnO = c(0.21, 0.20, 0.17, 0.23, 0.23, 0.23, 0.20, 0.24, 0.26, 0.26)
  )

  b_extra <- make_block(
    sample_id = 982,
    rock_type = "Basalt",
    SiO2 = 46.84,
    Al2O3 = 15.73,
    Fe2O3 = 11.68,
    MgO = 9.72,
    CaO = 12.63,
    Na2O = 1.78,
    K2O = 0.04,
    TiO2 = 1.11,
    P2O5 = 0.09,
    MnO = 0.23
  )

  list(
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b_extra = b_extra
  )
}

make_skye_lavas_33_basalt_raw <- function() {
  blocks <- make_skye_lavas_blocks()
  output <- rbind(blocks$b1, blocks$b2, blocks$b3, blocks$b_extra)
  rownames(output) <- output$sample_id
  output
}

make_skye_lavas_aitchison_32_raw <- function() {
  blocks <- make_skye_lavas_blocks()
  output <- rbind(blocks$b1, blocks$b2, blocks$b3)
  rownames(output) <- output$sample_id
  output
}

logistic_gaussian_screening_dataset_registry <- function() {
  canonical_entries <- lapply(composition_registry, function(entry) {
    list(
      source_type = "compositions_data",
      candidate_names = entry$data_entry,
      candidate_packages = "compositions",
      compositional_columns = entry$parts,
      source_object = entry$object_name,
      strict_compositional_validation = TRUE,
      notes = "Canonical compositional-column selection from the compositions package."
    )
  })
  names(canonical_entries) <- names(composition_registry)

  c(
    canonical_entries,
    list(
      SkyeLavas = list(
        source_type = "local_constructed",
        constructor = make_skye_lavas_33_basalt_raw,
        compositional_columns = c("SiO2", "Al2O3", "Fe2O3", "MgO", "CaO", "Na2O", "K2O", "TiO2", "P2O5", "MnO"),
        notes = c(
          "Skye lavas dataset reconstructed from Thompson, Esson and Duncan (1972), Table 2.",
          "This version contains the 33 basalt specimens from the first three basalt blocks plus sample 982, and uses the 10 major oxides as compositional parts."
        )
      ),
      SkyeLavasAitchison32 = list(
        source_type = "local_constructed",
        constructor = make_skye_lavas_aitchison_32_raw,
        compositional_columns = c("SiO2", "Al2O3", "Fe2O3", "MgO", "CaO", "Na2O", "K2O", "TiO2", "P2O5", "MnO"),
        notes = c(
          "Aitchison's 32 basalt Skye lavas reconstructed from Thompson, Esson and Duncan (1972), Table 2.",
          "This version keeps only the first three basalt blocks and uses the 10 major oxides as compositional parts."
        )
      ),
      SkyeLavasComplete = list(
        source_type = "local_constructed",
        constructor = make_skye_lavas_33_basalt_raw,
        compositional_columns = c("SiO2", "Al2O3", "Fe2O3", "MgO", "CaO", "Na2O", "K2O", "TiO2", "P2O5", "MnO"),
        notes = c(
          "Alias of SkyeLavas: Skye lavas reconstructed from Thompson, Esson and Duncan (1972), Table 2.",
          "This alias currently points to the 33-basalt version, i.e. the first three basalt blocks plus sample 982."
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
        compositional_columns = c("acit", "metpyr", "furfu", "furfualc", "dimeth", "met5"),
        notes = c(
          "Coffee composition dataset from robCompositions.",
          "The factor `sort` is a coffee-type label and is excluded from the compositional analysis.",
          "The six volatile compounds are treated compositionally."
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
      ClamCombined = list(
        source_type = "combined_prepared",
        components = c("ClamEast", "ClamWest"),
        group_variable = "site",
        source_object = "ClamCombined",
        notes = c(
          "Combined East and West clam color-size compositions.",
          "A site label is retained for ilr diagnostics because the combined sample may be a mixture."
        )
      )
      )
    )
}

slugify_dataset_name <- function(name) {
  gsub("[^a-z0-9]+", "_", tolower(name))
}

logistic_gaussian_screening_directories <- function(base_dir = canonical_logistic_gaussian_screening_dir("slow")) {
  list(
    base = base_dir,
    plots = file.path(base_dir, "plots"),
    metadata = file.path(base_dir, "metadata")
  )
}

ensure_logistic_gaussian_screening_directories <- function(base_dir = canonical_logistic_gaussian_screening_dir("slow")) {
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
    return("plugin_simple_null")
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
  statistic_node_source <- NA_character_

  for (container_name in c("inference", "statistics", "statistic", "results", "observed")) {
    if (!is.null(bootstrap_result[[container_name]][[statistic_name]])) {
      statistic_node <- bootstrap_result[[container_name]][[statistic_name]]
      statistic_node_source <- sprintf("%s$%s", container_name, statistic_name)
      break
    }
  }
  if (is.null(statistic_node) && !is.null(bootstrap_result[[statistic_name]])) {
    statistic_node <- bootstrap_result[[statistic_name]]
    statistic_node_source <- statistic_name
  }
  if (is.null(statistic_node)) {
    statistic_node <- bootstrap_result
    statistic_node_source <- "bootstrap_result"
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
    bootstrap_statistics = bootstrap_statistics,
    statistic_node_source = statistic_node_source,
    statistic_node_names = if (is.list(statistic_node)) paste(names(statistic_node), collapse = ";") else NA_character_
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

estimate_theoretical_profile_from_null_sample <- function(omega_grid,
                                                          t_grid,
                                                          x_null_sample,
                                                          fit) {
  z_centers <- ilr_transform_closed(as.matrix(omega_grid), V = fit$ilr_basis)
  z_null <- ilr_transform_closed(as.matrix(x_null_sample), V = fit$ilr_basis)
  distance_null <- distance_matrix_from_ilr(z_null, z_centers)
  compute_profile_matrix_from_distances(distance_null, t_grid)
}

plugin_parametric_bootstrap_logistic_gaussian_screening <- function(data_closed,
                                                                    fit,
                                                                    omega_grid,
                                                                    t_grid,
                                                                    theoretical_profile,
                                                                    observed_ks,
                                                                    observed_cvm,
                                                                    B,
                                                                    seed,
                                                                    n_cores = 1L) {
  data_closed <- as.matrix(data_closed)
  n <- nrow(data_closed)
  z_centers <- ilr_transform_closed(as.matrix(omega_grid), V = fit$ilr_basis)
  B <- as.integer(B)
  n_cores <- max(1L, as.integer(n_cores))

  run_one <- function(b) {
    x_boot <- simulate_fitted_logistic_gaussian(
      n = n,
      fit = fit,
      seed = as.integer(seed) + b
    )
    z_boot <- ilr_transform_closed(x_boot, V = fit$ilr_basis)
    distance_boot <- distance_matrix_from_ilr(z_boot, z_centers)
    empirical_boot <- compute_profile_matrix_from_distances(distance_boot, t_grid)
    compute_screening_statistics_from_profiles(
      empirical_profile = empirical_boot,
      theoretical_profile = theoretical_profile,
      n_obs = n
    )
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L && B > 1L) {
    bootstrap_stats <- parallel::mclapply(
      X = seq_len(B),
      FUN = run_one,
      mc.cores = min(n_cores, B)
    )
  } else {
    bootstrap_stats <- lapply(seq_len(B), run_one)
  }

  ks_statistics <- vapply(bootstrap_stats, function(x) x$ks, numeric(1))
  cvm_statistics <- vapply(bootstrap_stats, function(x) x$cvm, numeric(1))

  list(
    ks_statistics = ks_statistics,
    cvm_statistics = cvm_statistics,
    p_value_ks = mean(ks_statistics >= observed_ks),
    p_value_cvm = mean(cvm_statistics >= observed_cvm)
  )
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
      n_duplicate_rows = sum(duplicated(x_comp)),
      n_missing_rows_removed = 0L,
      min_entry = min(x_comp),
      row_sums_before_closure = rowSums(x_comp),
      component_names = colnames(x_comp),
      source_object = spec$source_object %||% name,
      group = group,
      group_variable = spec$group_variable %||% "group",
      notes = spec$notes,
      status = "ok"
    ))
  }

  if (identical(spec$source_type, "curated_metagenomic_data")) {
    return(prepare_curated_metagenomic_placeholder(name, spec))
  }

  if (identical(spec$source_type, "local_constructed")) {
    if (!is.function(spec$constructor)) {
      stop(sprintf("Dataset %s has source_type 'local_constructed' but no valid constructor.", name))
    }
    loaded <- list(
      object = spec$constructor(),
      package = "local",
      name = name,
      searched = name
    )
  } else if (identical(spec$source_type, "compositions_data")) {
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

  strict_validation <- isTRUE(spec$strict_compositional_validation)
  if (strict_validation && anyNA(x_comp_numeric)) {
    stop(sprintf("Canonical compositions dataset %s contains missing selected parts.", name))
  }
  if (strict_validation && any(is.infinite(x_comp_numeric))) {
    stop(sprintf("Canonical compositions dataset %s contains infinite selected parts.", name))
  }
  if (strict_validation && any(x_comp_numeric == 0)) {
    stop(sprintf("Canonical compositions dataset %s contains zero selected parts.", name))
  }
  if (strict_validation && any(x_comp_numeric < 0)) {
    stop(sprintf("Canonical compositions dataset %s contains negative selected parts.", name))
  }
  if (strict_validation && anyDuplicated(x_comp_numeric)) {
    stop(sprintf("Canonical compositions dataset %s contains duplicate rows after selecting its registered parts.", name))
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

  n_duplicate_rows <- sum(duplicated(x_comp_complete))

  zero_info <- apply_zero_replacement(
    x = x_comp_complete,
    zero_replacement = spec$zero_replacement %||% NULL
  )
  x_comp_complete <- zero_info$x

  has_zeros <- any(x_comp_complete == 0)
  row_sums_before_closure <- rowSums(x_comp_complete)
  x_closed <- close_composition_rows(x_comp_complete)

  closure_error <- max(abs(rowSums(x_closed) - 1))
  closure_tolerance <- if (strict_validation) 1e-12 else 1e-10
  if (closure_error > closure_tolerance) {
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
    n_duplicate_rows = n_duplicate_rows,
    n_missing_rows_removed = n_missing_rows,
    min_entry = min(x_comp_complete),
    row_sums_before_closure = row_sums_before_closure,
    component_names = colnames(x_comp_complete),
    source_object = spec$source_object %||% loaded$name,
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

simplex_lattice_count <- function(ambient_dim, level) {
  choose(level + ambient_dim - 1L, ambient_dim - 1L)
}

simplex_lattice_indices_recursive <- function(ambient_dim, level) {
  if (ambient_dim == 1L) {
    return(matrix(level, nrow = 1L, ncol = 1L))
  }

  blocks <- lapply(0:level, function(first_part) {
    tail_block <- simplex_lattice_indices_recursive(ambient_dim - 1L, level - first_part)
    cbind(first_part, tail_block)
  })
  do.call(rbind, blocks)
}

choose_simplex_lattice_level <- function(ambient_dim, target_points, overshoot_factor = 1.25) {
  target_points <- as.integer(target_points)
  if (target_points <= 0L) {
    stop("`target_points` must be strictly positive.")
  }

  level <- 1L
  counts <- integer(0)
  repeat {
    count_level <- simplex_lattice_count(ambient_dim, level)
    counts <- c(counts, count_level)
    if (count_level >= target_points * overshoot_factor || level >= 25L) {
      break
    }
    level <- level + 1L
  }

  candidate_levels <- seq_along(counts)
  best <- candidate_levels[[which.min(abs(counts - target_points))]]
  as.integer(best)
}

build_fixed_simplex_omega_grid <- function(ambient_dim,
                                           max_centers = 100L,
                                           boundary_epsilon = NULL) {
  ambient_dim <- as.integer(ambient_dim)
  max_centers <- as.integer(max_centers)
  if (ambient_dim < 2L) {
    stop("`ambient_dim` must be at least 2.")
  }
  if (max_centers <= 0L) {
    stop("`max_centers` must be strictly positive.")
  }
  if (is.null(boundary_epsilon)) {
    boundary_epsilon <- boundary_epsilon <- 0.015 / ambient_dim
  }
  if (!is.finite(boundary_epsilon) || boundary_epsilon <= 0 || boundary_epsilon >= 1 / ambient_dim) {
    stop("`boundary_epsilon` must belong to (0, 1 / ambient_dim).")
  }

  level <- choose_simplex_lattice_level(
    ambient_dim = ambient_dim,
    target_points = max_centers
  )
  lattice_idx <- simplex_lattice_indices_recursive(
    ambient_dim = ambient_dim,
    level = level
  )
  omega_grid <- lattice_idx / level
  omega_grid <- (1 - ambient_dim * boundary_epsilon) * omega_grid + boundary_epsilon
  omega_grid <- omega_grid / rowSums(omega_grid)

  list(
    omega = omega_grid,
    lattice_level = level,
    n_centers = nrow(omega_grid),
    boundary_epsilon = boundary_epsilon,
    construction = "fixed_simplex_lattice"
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

build_fixed_t_grid_logistic_gaussian <- function(fit,
                                                 omega_grid,
                                                 n_t = 60L,
                                                 tail_prob = 1e-8) {
  omega_grid <- as.matrix(omega_grid)
  omega_ilr <- ilr_transform_closed(omega_grid, V = fit$ilr_basis)
  shift_norms <- sqrt(rowSums((omega_ilr - matrix(
    rep(fit$mu_hat, each = nrow(omega_ilr)),
    nrow = nrow(omega_ilr),
    ncol = length(fit$mu_hat)
  ))^2))

  lambda_max <- max(fit$eigenvalues, 0)
  ilr_dim <- ncol(fit$Z)

  if (!is.finite(tail_prob) || tail_prob <= 0 || tail_prob >= 1) {
    stop("`tail_prob` must belong to (0, 1).")
  }

  # The KS grid in t should cover the distances that are relevant under the
  # fitted logistic Gaussian model. In ilr coordinates, the fitted model is
  # approximately N(mu_hat, Sigma_hat). If Y ~ N(0, I_q), q = ilr_dim, then
  # ||Z - mu_hat||_2 <= sqrt(lambda_max(Sigma_hat)) * ||Y||_2 in the
  # conservative worst direction, and ||Y||_2^2 has a chi-square_q law.
  # Therefore this radius is a high-probability conservative bound for the
  # fitted Gaussian spread, avoiding an arbitrary fixed multiplier such as
  # eight effective standard deviations.
  gaussian_radius <- sqrt(max(lambda_max, 0) * stats::qchisq(
    p = 1 - tail_prob,
    df = ilr_dim
  ))

  t_max <- max(shift_norms) + gaussian_radius

  if (!is.finite(t_max) || t_max <= 0) {
    t_max <- 1
  }

  list(
    t_grid = seq(0, t_max, length.out = as.integer(n_t)),
    t_max = t_max,
    tail_prob = tail_prob,
    construction = "fixed_interval_from_fitted_high_probability_radius"
  )
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

build_combined_simplex_dataset_plot_data <- function(dataset_names) {
  ternary_datasets <- character(0)
  plot_blocks <- list()

  for (dataset_name in dataset_names) {
    data_prep <- tryCatch(
      prepare_composition_dataset(dataset_name),
      error = function(e) NULL
    )
    if (is.null(data_prep) || !identical(data_prep$status %||% "ok", "ok")) {
      next
    }
    if (!is.matrix(data_prep$X_closed) || ncol(data_prep$X_closed) != 3L) {
      next
    }

    proj <- project_simplex_to_ternary(data_prep$X_closed)
    plot_blocks[[length(plot_blocks) + 1L]] <- data.frame(
      dataset = dataset_name,
      x = proj$x,
      y = proj$y,
      stringsAsFactors = FALSE
    )
    ternary_datasets <- c(ternary_datasets, dataset_name)
  }

  if (length(plot_blocks) == 0L) {
    return(NULL)
  }

  list(
    plot_data = do.call(rbind, plot_blocks),
    dataset_names = ternary_datasets
  )
}

save_combined_simplex_dataset_plot <- function(dataset_names, file) {
  plot_info <- build_combined_simplex_dataset_plot_data(dataset_names)
  if (is.null(plot_info)) {
    return(invisible(NULL))
  }

  simplex_boundary <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3) / 2, 0)
  )

  plt <- ggplot() +
    geom_path(
      data = simplex_boundary,
      aes(x = x, y = y),
      linewidth = 0.6
    ) +
    geom_point(
      data = plot_info$plot_data,
      aes(x = x, y = y, colour = dataset),
      size = 2.2,
      alpha = 0.85
    ) +
    coord_equal(
      xlim = c(-0.04, 1.04),
      ylim = c(-0.03, sqrt(3) / 2 + 0.04),
      expand = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL,
      colour = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = c(0.93, 0.91),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA),
      legend.margin = margin(2, 2, 2, 2)
    )

  ggplot2::ggsave(filename = file, plot = plt, width = 7, height = 5, dpi = 300)
  invisible(file)
}

save_fixed_simplex_omega_grid_plot <- function(max_centers, boundary_epsilon, file) {
  centers <- build_fixed_simplex_omega_grid(
    ambient_dim = 3L,
    max_centers = max_centers,
    boundary_epsilon = boundary_epsilon
  )
  proj <- project_simplex_to_ternary(centers$omega)
  simplex_boundary <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3) / 2, 0)
  )

  plt <- ggplot() +
    geom_path(
      data = simplex_boundary,
      aes(x = x, y = y),
      linewidth = 0.6,
      colour = "grey25"
    ) +
    geom_point(
      data = proj,
      aes(x = x, y = y),
      colour = "#2166ac",
      size = 2.1,
      alpha = 0.95
    ) +
    coord_equal(
      xlim = c(-0.04, 1.04),
      ylim = c(-0.03, sqrt(3) / 2 + 0.04),
      expand = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf(
        "Fixed simplex lattice for KS (D = 3, %d centers, boundary epsilon = %.6g)",
        nrow(centers$omega),
        centers$boundary_epsilon
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  ggplot2::ggsave(filename = file, plot = plt, width = 6, height = 5, dpi = 300)
  invisible(file)
}

save_batch_screening_plots <- function(dataset_names,
                                       base_dir,
                                       max_centers,
                                       boundary_epsilon) {
  dirs <- ensure_logistic_gaussian_screening_directories(base_dir)
  unlink(file.path(dirs$plots, "*.png"))

  combined_plot <- file.path(dirs$plots, "simplex_d3_datasets.png")
  ks_grid_plot <- file.path(dirs$plots, "simplex_fixed_ks_lattice.png")

  output_paths <- list(
    simplex_datasets = save_combined_simplex_dataset_plot(
      dataset_names = dataset_names,
      file = combined_plot
    ),
    simplex_fixed_ks_lattice = save_fixed_simplex_omega_grid_plot(
      max_centers = max_centers,
      boundary_epsilon = boundary_epsilon,
      file = ks_grid_plot
    )
  )

  output_paths[!vapply(output_paths, is.null, logical(1))]
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

  if (length(problem_reasons) > 0L) {
    return(list(
      diagnosis = "problematic",
      use_in_paper = "no",
      why = paste(problem_reasons, collapse = "; ")
    ))
  }

  if (!isTRUE(diagnostics$computed %||% TRUE)) {
    return(list(
      diagnosis = "not_assessed",
      use_in_paper = NA_character_,
      why = "Auxiliary screening diagnostics were disabled; only the KS/CvM bootstrap inference was computed."
    ))
  }

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

strip_auxiliary_screening_artifacts <- function(result) {
  result$grid <- utils::modifyList(
    result$grid %||% list(),
    list(
      center_indices = integer(0),
      omega = NULL,
      t_grid = numeric(0),
      omega_grid_construction = "not_computed",
      omega_grid_lattice_level = NA_integer_,
      omega_grid_size = NA_integer_,
      omega_grid_boundary_epsilon = NA_real_,
      t_grid_construction = "not_computed",
      t_grid_t_max = NA_real_,
      t_grid_tail_prob = NA_real_
    )
  )
  result$theoretical_profile_info <- list(
    method = "not_computed",
    note = "The auxiliary theoretical profile was not evaluated.",
    theoretical_profile = NULL
  )
  result$observed <- list(
    empirical_profile = NULL,
    theoretical_profile = NULL,
    process_matrix = NULL,
    grid_ks = NA_real_,
    grid_cvm = NA_real_
  )
  result$diagnostics <- list(
    computed = FALSE,
    note = "Auxiliary screening profile and normality diagnostics were not computed."
  )
  result$settings$compute_auxiliary_diagnostics <- FALSE
  result$classification <- classify_logistic_gaussian_screening(result)
  result
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
                                              bootstrap_method,
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
      bootstrap_method = bootstrap_method,
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
                                            B = 5000L,
                                            max_centers = 100L,
                                            n_t = 60L,
                                            t_grid_tail_prob = 1e-8,
                                            boundary_epsilon = NULL,
                                            probs_t = NULL,
                                            bootstrap_mode = c("composite_multiplier", "composite_parametric", "plugin_simple_null"),
                                            seed = 123L,
                                            alpha = 0.05,
                                            ridge = 1e-8,
                                            n_cores = 1L,
                                            bootstrap_method = "reestimated",
                                            bootstrap_keep = list(
                                              observed_process = TRUE,
                                              bootstrap_statistics = TRUE,
                                              bootstrap_thetas = FALSE
                                            ),
                                            compute_auxiliary_diagnostics = TRUE,
                                            control = list(),
                                            omega_grid_type = c("sample_points", "fixed_simplex_lattice"),
                                            t_grid_type = c("sample_distances", "fixed_fitted_scale", "sample_quantiles"),
                                            null_mc_size = 20000L,
                                            make_plots = TRUE,
                                            save_outputs = TRUE,
                                            output_dir = canonical_logistic_gaussian_screening_dir("slow"),
                                            run_seed_sensitivity = FALSE,
                                            sensitivity_B = NULL,
                                            verbose = TRUE) {
  bootstrap_mode <- standardize_screening_bootstrap_mode(match.arg(bootstrap_mode))
  omega_grid_type <- match.arg(omega_grid_type)
  t_grid_type <- match.arg(t_grid_type)
  control <- normalize_logistic_gaussian_screening_control(control)
  if (!isTRUE(compute_auxiliary_diagnostics) && identical(bootstrap_mode, "plugin_simple_null")) {
    stop("`plugin_simple_null` requires the auxiliary screening profile; use a composite bootstrap when it is disabled.")
  }
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
      bootstrap_method = bootstrap_method,
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
      bootstrap_method = bootstrap_method,
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

  centers <- NULL
  t_grid_info <- list(
    construction = "not_computed",
    t_max = NA_real_,
    tail_prob = NA_real_
  )
  t_grid <- numeric(0)
  theoretical_profile <- NULL
  empirical_profile <- NULL
  observed_stats_grid <- list(process_matrix = NULL, ks = NA_real_, cvm = NA_real_)

  if (isTRUE(compute_auxiliary_diagnostics)) {
    centers <- if (identical(omega_grid_type, "fixed_simplex_lattice")) {
      build_fixed_simplex_omega_grid(
        ambient_dim = data_prep$D,
        max_centers = max_centers,
        boundary_epsilon = boundary_epsilon
      )
    } else {
      sample_centers <- choose_screening_centers(
        x_closed = data_prep$X_closed,
        max_centers = max_centers,
        seed = seed
      )
      c(
        sample_centers,
        list(
          lattice_level = NA_integer_,
          n_centers = nrow(sample_centers$omega),
          construction = "sample_points"
        )
      )
    }

    z_centers <- ilr_transform_closed(centers$omega, V = fit$ilr_basis)
    distance_observed <- distance_matrix_from_ilr(fit$Z, z_centers)

    t_grid_info <- if (identical(t_grid_type, "fixed_fitted_scale")) {
      build_fixed_t_grid_logistic_gaussian(
        fit = fit,
        omega_grid = centers$omega,
        n_t = n_t,
        tail_prob = t_grid_tail_prob
      )
    } else if (identical(t_grid_type, "sample_distances")) {
      list(
        t_grid = sort(unique(as.numeric(distance_observed))),
        t_max = NA_real_,
        tail_prob = NA_real_,
        construction = "sample_distances"
      )
    } else {
      list(
        t_grid = choose_screening_t_grid(
          distance_observed = distance_observed,
          distance_null = NULL,
          n_t = n_t,
          probs = probs_t
        ),
        t_max = NA_real_,
        tail_prob = NA_real_,
        construction = "sample_quantiles"
      )
    }
    t_grid <- t_grid_info$t_grid

    theoretical_profile <- if (identical(bootstrap_mode, "plugin_simple_null")) {
      x_null_mc <- simulate_fitted_logistic_gaussian(
        n = null_mc_size,
        fit = fit,
        seed = as.integer(seed) + 500000L
      )
      estimate_theoretical_profile_from_null_sample(
        omega_grid = centers$omega,
        t_grid = t_grid,
        x_null_sample = x_null_mc,
        fit = fit
      )
    } else {
      evaluate_fitted_logistic_gaussian_profile(
        omega_grid = centers$omega,
        t_grid = t_grid,
        fit = fit,
        control = control
      )
    }
    empirical_profile <- compute_profile_matrix_from_distances(distance_observed, t_grid)
    observed_stats_grid <- compute_screening_statistics_from_profiles(
      empirical_profile = empirical_profile,
      theoretical_profile = theoretical_profile,
      n_obs = data_prep$n
    )
  }

  if (identical(bootstrap_mode, "plugin_simple_null")) {
    bootstrap_result <- plugin_parametric_bootstrap_logistic_gaussian_screening(
      data_closed = data_prep$X_closed,
      fit = fit,
      omega_grid = centers$omega,
      t_grid = t_grid,
      theoretical_profile = theoretical_profile,
      observed_ks = observed_stats_grid$ks,
      observed_cvm = observed_stats_grid$cvm,
      B = B,
      seed = seed,
      n_cores = n_cores
    )
    ks_extracted <- list(
      statistic = observed_stats_grid$ks,
      p_value = bootstrap_result$p_value_ks,
      bootstrap_statistics = bootstrap_result$ks_statistics
    )
    cvm_extracted <- list(
      statistic = observed_stats_grid$cvm,
      p_value = bootstrap_result$p_value_cvm,
      bootstrap_statistics = bootstrap_result$cvm_statistics
    )
  } else {
    ensure_logistic_gaussian_bootstrap_available()
    bootstrap_result <- multiplier_bootstrap_logistic_gaussian(
      data = data_prep$X_closed,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = B,
      alpha = alpha,
      n_cores = n_cores,
      seed = seed,
      bootstrap_method = bootstrap_method,
      keep = bootstrap_keep,
      control = control,
      unknown_param = "both"
    )

    ks_extracted <- extract_metric_result(bootstrap_result, "ks")
    cvm_extracted <- extract_metric_result(bootstrap_result, "cvm")
  }

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

  diagnostics <- if (isTRUE(compute_auxiliary_diagnostics)) {
    list(
      computed = TRUE,
      mardia = mardia_multivariate_normality(fit$Z),
      shapiro_p_values = compute_marginal_shapiro_pvalues(fit$Z)
    )
  } else {
    list(
      computed = FALSE,
      note = "Auxiliary screening profile and normality diagnostics were not computed."
    )
  }

  result <- list(
    dataset_name = dataset_name,
    screening_type = "logistic_gaussian_composite_or_screening",
    warning = logistic_gaussian_screening_warning,
    data_prep = data_prep,
    fit = fit,
    grid = list(
      center_indices = if (isTRUE(compute_auxiliary_diagnostics)) centers$center_indices %||% NA_integer_ else integer(0),
      omega = if (isTRUE(compute_auxiliary_diagnostics)) centers$omega else NULL,
      t_grid = t_grid,
      probs_t = probs_t %||% seq(0.01, 0.99, length.out = n_t),
      max_centers = max_centers,
      omega_grid_type = omega_grid_type,
      omega_grid_construction = if (isTRUE(compute_auxiliary_diagnostics)) centers$construction %||% omega_grid_type else "not_computed",
      omega_grid_lattice_level = if (isTRUE(compute_auxiliary_diagnostics)) centers$lattice_level %||% NA_integer_ else NA_integer_,
      omega_grid_size = if (isTRUE(compute_auxiliary_diagnostics)) centers$n_centers %||% nrow(centers$omega) else NA_integer_,
      omega_grid_boundary_epsilon = if (isTRUE(compute_auxiliary_diagnostics)) centers$boundary_epsilon %||% NA_real_ else NA_real_,
      t_grid_type = t_grid_type,
      t_grid_construction = t_grid_info$construction %||% t_grid_type,
      t_grid_t_max = t_grid_info$t_max %||% NA_real_,
      t_grid_tail_prob = t_grid_info$tail_prob %||% NA_real_
    ),
    theoretical_profile_info = list(
      method = if (!isTRUE(compute_auxiliary_diagnostics)) "not_computed" else if (identical(bootstrap_mode, "plugin_simple_null")) "monte_carlo_from_fitted_null" else "global_logistic_gaussian_model_spec",
      note = if (!isTRUE(compute_auxiliary_diagnostics)) {
        "The auxiliary theoretical profile was not evaluated."
      } else if (identical(bootstrap_mode, "plugin_simple_null")) {
        "Theoretical distance profiles are approximated by a large Monte Carlo sample from the fitted logistic Gaussian null with parameters treated as fixed."
      } else {
        "Theoretical distance profiles are evaluated by make_logistic_gaussian_spec()$profile_matrix_eval, the same model-spec path used by the calibration code."
      },
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
      bootstrap_method = bootstrap_method,
      engine = if (identical(bootstrap_mode, "plugin_simple_null")) "plugin_parametric_bootstrap_logistic_gaussian_screening" else "multiplier_bootstrap_logistic_gaussian",
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
      n_t = n_t,
      t_grid_tail_prob = t_grid_tail_prob,
      boundary_epsilon = centers$boundary_epsilon %||% NA_real_,
      probs_t = probs_t %||% seq(0.01, 0.99, length.out = n_t),
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      bootstrap_mode = bootstrap_mode,
      bootstrap_keep = bootstrap_keep,
      compute_auxiliary_diagnostics = isTRUE(compute_auxiliary_diagnostics),
      n_cores = n_cores,
      omega_grid_type = omega_grid_type,
      t_grid_type = t_grid_type,
      null_mc_size = null_mc_size,
      quadform_method = control$mvnormal_quadform_method %||% control$logistic_gaussian_quadform_method %||% NA_character_,
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
      t_grid_tail_prob = t_grid_tail_prob,
      boundary_epsilon = boundary_epsilon,
      probs_t = probs_t,
      bootstrap_mode = bootstrap_mode,
      seed = alt_seed,
      alpha = alpha,
      ridge = ridge,
      n_cores = n_cores,
      bootstrap_keep = bootstrap_keep,
      compute_auxiliary_diagnostics = compute_auxiliary_diagnostics,
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

screening_dataset_metadata <- function(dataset_name) {
  canonical_entry <- composition_registry[[dataset_name]]
  if (!is.null(canonical_entry)) {
    return(list(
      data_source = "compositions",
      source_object = canonical_entry$object_name,
      selected_parts = paste(canonical_entry$parts, collapse = ",")
    ))
  }

  spec <- logistic_gaussian_screening_dataset_registry()[[dataset_name]]
  list(
    data_source = "external",
    source_object = spec$source_object %||% spec$study %||% dataset_name,
    selected_parts = paste(spec$compositional_columns %||% character(0), collapse = ",")
  )
}

make_logistic_gaussian_screening_summary_row <- function(result) {
  fit <- result$fit
  mardia <- result$diagnostics$mardia %||% list(
    skewness_p_value = NA_real_,
    kurtosis_p_value = NA_real_
  )
  metadata <- screening_dataset_metadata(result$dataset_name)
  selected_parts <- metadata$selected_parts
  if (!nzchar(selected_parts)) {
    selected_parts <- paste(result$data_prep$component_names %||% character(0), collapse = ",")
  }
  shapiro_min_pvalue <- suppressWarnings(min(result$diagnostics$shapiro_p_values %||% NA_real_, na.rm = TRUE))
  if (!is.finite(shapiro_min_pvalue)) {
    shapiro_min_pvalue <- NA_real_
  }

  data.frame(
    dataset = result$dataset_name,
    data_source = metadata$data_source,
    source_object = metadata$source_object,
    selected_parts = selected_parts,
    status = result$data_prep$status %||% "ok",
    result_status = if (identical(result$data_prep$status %||% "ok", "ok")) "completed" else "failed",
    source_package = result$data_prep$source_package %||% NA_character_,
    source_dataset_name = result$data_prep$source_dataset_name %||% NA_character_,
    n = result$data_prep$n,
    D = result$data_prep$D,
    zeros = if (isTRUE(result$data_prep$has_zeros)) "yes" else "no",
    missing = if (isTRUE(result$data_prep$has_missing)) "yes" else "no",
    duplicate_rows = result$data_prep$n_duplicate_rows %||% NA_integer_,
    missing_rows_removed = result$data_prep$n_missing_rows_removed,
    zero_replacement = result$data_prep$zero_replacement %||% NA_real_,
    n_zeros_replaced = result$data_prep$n_zeros_replaced %||% NA_integer_,
    min_component = result$data_prep$min_entry,
    ridge_added = if (is.null(fit)) NA_real_ else fit$ridge_added,
    min_eigenvalue = if (is.null(fit)) NA_real_ else min(fit$eigenvalues),
    max_eigenvalue = if (is.null(fit)) NA_real_ else max(fit$eigenvalues),
    condition_number = if (is.null(fit)) NA_real_ else fit$condition_number,
    boundary_epsilon = result$settings$boundary_epsilon %||% result$grid$omega_grid_boundary_epsilon %||% NA_real_,
    omega_grid_construction = result$grid$omega_grid_construction %||% NA_character_,
    n_centers = result$grid$omega_grid_size %||% NA_integer_,
    t_grid_max = result$grid$t_grid_t_max %||% NA_real_,
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
                                                  B = 5000L,
                                                  max_centers = 100L,
                                                  n_t = 60L,
                                                  t_grid_tail_prob = 1e-8,
                                                  boundary_epsilon = NULL,
                                                  probs_t = NULL,
                                                  bootstrap_mode = "composite_multiplier",
                                                  seed = 123L,
                                                  alpha = 0.05,
                                                  ridge = 1e-8,
                                                  n_cores = 1L,
                                                  bootstrap_method = "reestimated",
                                                  control = list(),
                                                  omega_grid_type = "sample_points",
                                                  t_grid_type = "sample_distances",
                                                  make_plots = TRUE,
                                                  output_dir = canonical_logistic_gaussian_screening_dir("slow"),
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
        t_grid_tail_prob = t_grid_tail_prob,
        boundary_epsilon = boundary_epsilon,
        probs_t = probs_t,
        bootstrap_mode = bootstrap_mode,
        seed = dataset_seed,
        alpha = alpha,
        ridge = ridge,
        n_cores = n_cores,
        bootstrap_method = bootstrap_method,
        control = control,
        omega_grid_type = omega_grid_type,
        t_grid_type = t_grid_type,
        make_plots = FALSE,
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
            status = "failed",
            n = NA_integer_,
            D = NA_integer_,
            has_zeros = NA,
            has_missing = NA,
            n_duplicate_rows = NA_integer_,
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
          bootstrap = list(mode = bootstrap_mode, engine = NA_character_, bootstrap_method = bootstrap_method),
          settings = list(B = B, seed = dataset_seed, n_cores = n_cores, bootstrap_method = bootstrap_method, control = control),
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
      t_grid_tail_prob = t_grid_tail_prob,
      boundary_epsilon = boundary_epsilon,
      probs_t = probs_t %||% seq(0.01, 0.99, length.out = n_t),
      bootstrap_mode = bootstrap_mode,
      seed = seed,
      alpha = alpha,
      ridge = ridge,
      n_cores = n_cores,
      bootstrap_method = bootstrap_method,
      control = control,
      result_paths = result_paths
    ),
    file = config_file
  )

  plot_paths <- NULL
  if (isTRUE(make_plots)) {
    plot_paths <- save_batch_screening_plots(
      dataset_names = dataset_names,
      base_dir = output_dir,
      max_centers = max_centers,
      boundary_epsilon = boundary_epsilon
    )
  }

  list(
    summary = summary_df,
    summary_csv = summary_csv,
    session_info_file = session_info_file,
    package_versions_file = package_versions_file,
    warning_file = warning_file,
    config_file = config_file,
    plot_paths = plot_paths
  )
}
source(file.path("scripts", "path_helpers.R"))
