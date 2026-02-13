# terralink raster runner (QGIS-free)
#
# Usage in R:
#   source("/Users/benbishop/projects/terralink/inst/scripts/run_raster.R")

clear_environment <- function() {
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  invisible(gc())
}
# Uncomment to clear workspace on run.
# clear_environment()

pkg_path <- "/Users/benbishop/projects/terralink"
if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(pkg_path, "DESCRIPTION"))) {
  suppressMessages(devtools::load_all(pkg_path, quiet = TRUE))
} else if (requireNamespace("terralink", quietly = TRUE)) {
  # Package installed; namespace will be used below.
} else {
  stop("terralink is not available. Install it or update pkg_path.")
}

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("terra is required. Please install.packages('terra').")
}

choose_option <- function(options, title) {
  if (requireNamespace("tcltk", quietly = TRUE)) {
    pick <- tcltk::tk_select.list(options, title = title, multiple = FALSE)
    if (length(pick) && nzchar(pick)) {
      return(pick)
    }
  }
  choice <- utils::menu(options, title = title)
  if (choice > 0) {
    return(options[[choice]])
  }
  options[[1]]
}

print_value_counts <- function(raster, max_rows = 10) {
  vals <- terra::freq(raster)
  if (is.null(vals) || nrow(vals) == 0) {
    cat("Value counts: none\n")
    return(invisible(NULL))
  }
  vals <- vals[order(vals$count, decreasing = TRUE), , drop = FALSE]
  if (nrow(vals) > max_rows) {
    vals <- vals[seq_len(max_rows), , drop = FALSE]
  }
  cat("Top raster values by count:\n")
  print(vals)
  invisible(NULL)
}

diagnose_patch_counts <- function(raster, patch_values, min_patch_size) {
  mask <- terralink::raster_mask_from_values(raster, patch_values)
  stats <- list()
  for (conn in c(4, 8)) {
    labels <- terralink::label_patches(mask, connectivity = conn)
    freq <- terra::freq(labels)
    if (is.null(freq) || nrow(freq) == 0) {
      cat("Connectivity", conn, ": 0 raw patches\n")
      stats[[as.character(conn)]] <- list(raw = 0L, kept = 0L)
      next
    }
    freq <- freq[!is.na(freq$value) & freq$value != 0, , drop = FALSE]
    raw_count <- nrow(freq)
    kept <- freq[freq$count >= min_patch_size, , drop = FALSE]
    cat(
      "Connectivity", conn,
      ": raw patches =", raw_count,
      ", patches >= min_patch_size =", nrow(kept),
      "\n"
    )
    stats[[as.character(conn)]] <- list(raw = raw_count, kept = nrow(kept))
  }
  stats
}

# ---------------------------
# USER SETTINGS (edit here)
# ---------------------------
input_path <- ""   # Leave blank to choose a file.
output_dir <- ""   # Leave blank to create <input>_output folder.
strategy <- "most_connectivity"  # "most_connectivity" or "largest_network".

units <- "pixels"  # "pixels", "metric", or "imperial".

patch_values <- c(12)       # Habitat cell values in the raster.
patch_ranges <- list()      # Optional list of ranges, e.g. list(c(5, 7), c(10, 12)).
obstacle_values <- c()      # Optional impassable values.
obstacle_ranges <- list()   # Optional impassable ranges.

budget <- 100               # Total corridor budget (units above).
min_patch_size <- 20        # Minimum patch size (units above).
min_corridor_width <- 3     # Minimum corridor width (units above).
max_search_distance <- 300  # Maximum search distance (units above).
allow_bottlenecks <- FALSE
patch_connectivity <- 8     # 4 or 8. Leave NA to choose on run.
auto_patch_connectivity <- TRUE  # If TRUE, choose 4/8 based on patch counts.
# ---------------------------

if (!nzchar(input_path)) {
  input_path <- file.choose()
}
input_path <- normalizePath(input_path, mustWork = TRUE)

if (!nzchar(output_dir)) {
  stem <- tools::file_path_sans_ext(basename(input_path))
  output_dir <- file.path(dirname(input_path), paste0(stem, "_output"))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (is.na(patch_connectivity)) {
  pick <- choose_option(c("4", "8"), "Pixel neighborhood (4 or 8)")
  patch_connectivity <- as.integer(pick)
}

r <- terra::rast(input_path)

cat("Raster input:", input_path, "\n")
cat("Output dir:", output_dir, "\n")
cat("Strategy:", strategy, "\n")
cat("Patch connectivity:", patch_connectivity, "\n")
print_value_counts(r)
counts <- diagnose_patch_counts(r, patch_values, min_patch_size)
if (isTRUE(auto_patch_connectivity) && !is.na(patch_connectivity)) {
  raw4 <- counts[["4"]]$raw
  raw8 <- counts[["8"]]$raw
  if (raw4 > raw8 && raw4 > 1 && patch_connectivity != 4) {
    cat("Auto-switching patch_connectivity to 4 based on patch counts.\n")
    patch_connectivity <- 4
  } else if (raw8 > raw4 && raw8 > 1 && patch_connectivity != 8) {
    cat("Auto-switching patch_connectivity to 8 based on patch counts.\n")
    patch_connectivity <- 8
  }
}

result <- terralink::terralink_raster(
  raster = r,
  patch_values = patch_values,
  patch_ranges = patch_ranges,
  budget = budget,
  strategy = strategy,
  min_patch_size = min_patch_size,
  min_corridor_width = min_corridor_width,
  max_search_distance = max_search_distance,
  obstacle_values = obstacle_values,
  obstacle_ranges = obstacle_ranges,
  allow_bottlenecks = allow_bottlenecks,
  patch_connectivity = patch_connectivity,
  units = units
)

saveRDS(result, file.path(output_dir, "terralink_result.rds"))
if (!is.null(result$corridor_raster)) {
  terra::writeRaster(result$corridor_raster, file.path(output_dir, "corridors.tif"), overwrite = TRUE)
}
if (!is.null(result$patches)) {
  terra::writeRaster(result$patches, file.path(output_dir, "patches.tif"), overwrite = TRUE)
}
if (!is.null(result$patch_table)) {
  utils::write.csv(result$patch_table, file.path(output_dir, "patches.csv"), row.names = FALSE)
}
if (!is.null(result$corridors)) {
  utils::write.csv(result$corridors, file.path(output_dir, "corridors.csv"), row.names = FALSE)
}
if (!is.null(result$contiguous_raster)) {
  terra::writeRaster(result$contiguous_raster, file.path(output_dir, "contiguous.tif"), overwrite = TRUE)
}
if (!is.null(result$metrics_report)) {
  writeLines(result$metrics_report, con = file.path(output_dir, "landscape_metrics.txt"))
}

cat("Saved result to:", file.path(output_dir, "terralink_result.rds"), "\n")
cat("Reload with: result <- readRDS(\"", file.path(output_dir, "terralink_result.rds"), "\")\n", sep = "")
print(result$summary)
if (!is.null(result$summary$patches) && result$summary$patches < 2) {
  cat("Warning: only one patch found; no corridors can be created.\n")
}
if (!is.null(result$summary$corridors_used) && result$summary$corridors_used == 0) {
  cat("Warning: no corridors selected. Check patch_values/min_patch_size/max_search_distance/budget.\n")
}
