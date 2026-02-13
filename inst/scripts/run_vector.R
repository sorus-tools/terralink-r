# terralink vector runner (QGIS-free)
#
# Usage in R:
#   source("/Users/benbishop/projects/terralink/inst/scripts/run_vector.R")

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

# ---------------------------
# USER SETTINGS (edit here)
# ---------------------------
input_path <- ""   # Leave blank to choose a file.
output_dir <- ""   # Leave blank to create <input>_output folder.
strategy <- "most_connectivity"  # "most_connectivity" or "largest_network".
units <- "metric"  # "metric" or "imperial".

budget <- 50        # Corridor budget (ha/ac).
min_patch_size <- 5 # Minimum patch size (ha/ac).
min_corridor_width <- 200  # Corridor width (m/ft).
max_search_distance <- 5000 # Max search distance (m/ft).

obstacles_path <- ""  # Optional obstacles layer path (roads, urban, etc).
obstacle_resolution <- NULL # Optional raster resolution for obstacle routing.
# NOTE: obstacle-aware routing uses optional packages: gdistance, raster, sp.
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

obstacles <- NULL
if (nzchar(obstacles_path)) {
  obstacles <- obstacles_path
}

patches <- sf::st_read(input_path, quiet = TRUE)

result <- terralink::terralink_vector(
  patches = patches,
  budget = budget,
  strategy = strategy,
  min_patch_size = min_patch_size,
  min_corridor_width = min_corridor_width,
  max_search_distance = max_search_distance,
  obstacle_layers = obstacles,
  obstacle_resolution = obstacle_resolution,
  units = units
)

saveRDS(result, file.path(output_dir, "terralink_result.rds"))
if (!is.null(result$corridors)) {
  sf::st_write(result$corridors, file.path(output_dir, "corridors.gpkg"), delete_layer = TRUE, quiet = TRUE)
}
if (!is.null(result$networks)) {
  sf::st_write(result$networks, file.path(output_dir, "networks.gpkg"), delete_layer = TRUE, quiet = TRUE)
}
if (!is.null(result$metrics_report)) {
  writeLines(result$metrics_report, con = file.path(output_dir, "landscape_metrics.txt"))
}

cat("Saved result to:", file.path(output_dir, "terralink_result.rds"), "\n")
print(result$summary)
