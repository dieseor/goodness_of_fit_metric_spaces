source(file.path("real_data", "comets", "utils_comets_data.R"))

build_comets_audit_summary <- function(output_path = file.path("real_data", "comets", "comets_audit_summary.csv")) {
  comets_data <- load_comets_real_data(finite_normals = "both")

  small_circle_short <- utils::read.csv(
    file.path("real_data", "comets", "small_circle", "short_comets", "benchmark_summary.csv"),
    stringsAsFactors = FALSE
  )
  small_circle_long <- utils::read.csv(
    file.path("real_data", "comets", "small_circle", "long_comets", "benchmark_summary.csv"),
    stringsAsFactors = FALSE
  )
  jp_ks <- utils::read.csv(
    file.path("real_data", "comets", "jp", "ks_benchmark_main", "benchmark_summary.csv"),
    stringsAsFactors = FALSE
  )
  jp_cvm <- utils::read.csv(
    file.path("real_data", "comets", "jp", "short_cvm_nm_full_boot_vs_local_fast_B50_100_200_500_1000_comparison.csv"),
    stringsAsFactors = FALSE
  )
  beta_mix <- utils::read.csv(
    file.path("real_data", "comets", "mixture", "beta_mixture2_short_long_B1000", "comets_rotational_mixtures_summary.csv"),
    stringsAsFactors = FALSE
  )

  latest_small_circle_short <- small_circle_short[which.max(small_circle_short$B), ]
  latest_small_circle_long <- small_circle_long[which.max(small_circle_long$B), ]
  latest_jp_ks <- jp_ks[which.max(jp_ks$B), ]
  latest_jp_cvm <- jp_cvm[which.max(jp_cvm$B), ]

  add_row <- function(dataset,
                      model,
                      n,
                      tests_available,
                      ks_p_value,
                      cvm_p_value,
                      source_path,
                      notes) {
    data.frame(
      dataset = dataset,
      model = model,
      n = n,
      tests_available = tests_available,
      ks_p_value = ks_p_value,
      cvm_p_value = cvm_p_value,
      source_path = source_path,
      notes = notes,
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    add_row(
      dataset = "short_period",
      model = "small_circle",
      n = nrow(comets_data$short),
      tests_available = "CvM",
      ks_p_value = NA_real_,
      cvm_p_value = latest_small_circle_short$p_value,
      source_path = "real_data/comets/small_circle/short_comets/benchmark_summary.csv",
      notes = "Only CvM benchmark stored; latest row uses B=M=1000."
    ),
    add_row(
      dataset = "long_period",
      model = "small_circle",
      n = nrow(comets_data$long),
      tests_available = "CvM",
      ks_p_value = NA_real_,
      cvm_p_value = latest_small_circle_long$p_value,
      source_path = "real_data/comets/small_circle/long_comets/benchmark_summary.csv",
      notes = "Only CvM benchmark stored; latest row uses B=M=1000."
    ),
    add_row(
      dataset = "short_period",
      model = "JP",
      n = nrow(comets_data$short),
      tests_available = "KS and CvM",
      ks_p_value = latest_jp_ks$p_value,
      cvm_p_value = latest_jp_cvm$p_value_full_boot_nm,
      source_path = "real_data/comets/jp/{ks_benchmark_main,short_cvm_nm_full_boot_vs_local_fast_B50_100_200_500_1000_comparison.csv}",
      notes = "KS main benchmark exists only for short-period data; CvM comparison file also exists for short-period data."
    ),
    add_row(
      dataset = "long_period",
      model = "JP",
      n = nrow(comets_data$long),
      tests_available = "none stored",
      ks_p_value = NA_real_,
      cvm_p_value = NA_real_,
      source_path = "",
      notes = "No long-period JP benchmark artifacts were found in real_data/comets."
    ),
    add_row(
      dataset = "short_period",
      model = "beta_mixture2",
      n = nrow(comets_data$short),
      tests_available = "KS and CvM",
      ks_p_value = beta_mix$gof_ks_p_value[beta_mix$dataset == "short_period"][1L],
      cvm_p_value = beta_mix$gof_cvm_p_value[beta_mix$dataset == "short_period"][1L],
      source_path = "real_data/comets/mixture/beta_mixture2_short_long_B1000/comets_rotational_mixtures_summary.csv",
      notes = "Authoritative real_data copy uses B=1000 and M=60."
    ),
    add_row(
      dataset = "long_period",
      model = "beta_mixture2",
      n = nrow(comets_data$long),
      tests_available = "KS and CvM",
      ks_p_value = beta_mix$gof_ks_p_value[beta_mix$dataset == "long_period"][1L],
      cvm_p_value = beta_mix$gof_cvm_p_value[beta_mix$dataset == "long_period"][1L],
      source_path = "real_data/comets/mixture/beta_mixture2_short_long_B1000/comets_rotational_mixtures_summary.csv",
      notes = "A stale duplicate exists under output/real_data/comets with B=10 and M=50; do not use it for the AoS table."
    )
  )

  summary_df <- do.call(rbind, rows)
  utils::write.csv(summary_df, file = output_path, row.names = FALSE)
  summary_df
}

if (sys.nframe() == 0L) {
  build_comets_audit_summary()
}
