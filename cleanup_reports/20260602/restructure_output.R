base <- 'output'
structured_root <- file.path(base, 'structured')
dir.create(structured_root, recursive = TRUE, showWarnings = FALSE)

areas <- c('calibration','real_data','distance_profiles','convergence','diagnostics','comets','other')
for (a in areas) dir.create(file.path(structured_root, a), recursive = TRUE, showWarnings = FALSE)

classify_area <- function(x) {
  lx <- tolower(x)
  if (grepl('bootstrap|calibration', lx)) return('calibration')
  if (grepl('screening|realdata|real_data|dataset', lx)) return('real_data')
  if (grepl('distance|geodesic|chordal|poor_mans|plots', lx)) return('distance_profiles')
  if (grepl('gaussian_process|convergence', lx)) return('convergence')
  if (grepl('diagnostic|validation|mle', lx)) return('diagnostics')
  if (grepl('comets', lx)) return('comets')
  'other'
}

skip <- c('_catalog', 'structured')
all_top <- list.dirs(base, full.names = FALSE, recursive = FALSE)
move_candidates <- setdiff(all_top, skip)

log_rows <- list()
for (nm in sort(move_candidates)) {
  src <- file.path(base, nm)
  if (!dir.exists(src)) next
  area <- classify_area(nm)
  dst <- file.path(structured_root, area, nm)

  status <- 'skipped'
  note <- ''
  if (normalizePath(src, winslash='/', mustWork=FALSE) == normalizePath(dst, winslash='/', mustWork=FALSE)) {
    status <- 'noop'
  } else if (dir.exists(dst)) {
    status <- 'exists'
    note <- 'target_exists'
  } else {
    ok <- file.rename(src, dst)
    if (isTRUE(ok)) {
      status <- 'moved'
    } else {
      status <- 'failed'
      note <- 'rename_failed'
    }
  }

  log_rows[[length(log_rows) + 1L]] <- data.frame(
    name = nm,
    area = area,
    old_path = src,
    new_path = dst,
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )
}

log_df <- do.call(rbind, log_rows)
write.csv(log_df, file = 'cleanup_reports/20260602/output_restructure_map_20260602.csv', row.names = FALSE)

cat('Moved:', sum(log_df$status == 'moved'), '\n')
cat('Exists:', sum(log_df$status == 'exists'), '\n')
cat('Failed:', sum(log_df$status == 'failed'), '\n')
