# TerraLink RStudio Tryout Script
# Purpose: quick, clean manual testing of TerraLink in RStudio.
#
# Usage pattern:
# 1) Run the "Setup" section once per R session.
# 2) Run either the Raster section OR the Vector section.
# 3) Edit only the paths and parameters in the section you are running.
#
# This script is intentionally section-based (not meant to be run top-to-bottom every time).


# Setup ----

# Keep this TRUE while testing package code in this repo.
use_local_source <- TRUE
local_repo_path <- "/Users/benbishop/projects/terralink"

if (isTRUE(use_local_source)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload first: install.packages('pkgload')", call. = FALSE)
  }
  pkgload::load_all(local_repo_path, quiet = TRUE)
}

suppressPackageStartupMessages({
  library(terralink)
  library(terra)
  library(sf)
})

cat("Using terralink from:\n")
cat(normalizePath(getNamespaceInfo("terralink", "path"), winslash = "/"), "\n")
cat("terralink version:", as.character(utils::packageVersion("terralink")), "\n\n")

tl_assert_file <- function(path, label) {
  if (!is.character(path) || length(path) != 1 || !nzchar(path) || !file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

tl_assert_dir <- function(path, label) {
  if (!is.character(path) || length(path) != 1 || !nzchar(path) || !dir.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}


# Raster Tryout ----

# Inputs
raster_input_path <- "/Users/benbishop/Downloads/smallest_barbados.tif"
raster_output_dir <- "/Users/benbishop/R projects/LS test/outputs"

# Parameters
raster_args <- list(
  patch_values = 12,
  patch_ranges = NULL,
  obstacle_values = NULL,
  obstacle_ranges = NULL,
  budget = 50,
  min_patch_size = 2,
  min_corridor_width = 3,
  max_search_distance = 600,
  units = "pixels",
  allow_bottlenecks = FALSE,
  strategy = "circuit_utility",
  verbose = 1,
  progress = FALSE
)

tl_assert_file(raster_input_path, "Raster input")
tl_assert_dir(raster_output_dir, "Raster output directory")
raster_run_dir <- file.path(raster_output_dir, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(raster_run_dir, recursive = TRUE, showWarnings = FALSE)

raster_result <- do.call(
  terralink_raster,
  c(list(raster = raster_input_path), raster_args)
)

print(raster_result)
plot(raster_result)

raster_prefix <- paste0(tools::file_path_sans_ext(basename(raster_input_path)), "_r_tryout")
raster_outputs <- write_terralink_raster_outputs(
  result = raster_result,
  output_dir = raster_run_dir,
  prefix = raster_prefix,
  overwrite = TRUE,
  output_paths = list()
)

cat("Raster summary:\n")
print(raster_result$summary)
cat("\nRaster output files:\n")
print(raster_outputs)

raster_run <- list(result = raster_result, outputs = raster_outputs)


# Vector Tryout ----

# Inputs
vector_patches_path <- "/Users/benbishop/Downloads/small barbados vector.shp"
vector_output_dir <- "/Users/benbishop/R projects/LS test/outputs"

# Optional: obstacle layer(s). Leave NULL to disable.
vector_obstacle_layers <- NULL
# Example:
# vector_obstacle_layers <- c("/Users/benbishop/Downloads/obstacles.gpkg")

# Parameters
vector_args <- list(
  budget = 1,
  min_patch_size = 1,
  min_corridor_width = 20,
  max_search_distance = 5000,
  units = "metric",
  return_crs = "input",
  strategy = "circuit_utility",
  verbose = 1,
  progress = FALSE
)

tl_assert_file(vector_patches_path, "Vector patches input")
tl_assert_dir(vector_output_dir, "Vector output directory")
vector_run_dir <- file.path(vector_output_dir, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(vector_run_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(vector_obstacle_layers)) {
  for (p in vector_obstacle_layers) tl_assert_file(p, "Vector obstacle layer")
}

vector_call <- c(list(patches = vector_patches_path), vector_args)
if (!is.null(vector_obstacle_layers)) {
  vector_call$obstacle_layers <- vector_obstacle_layers
}

vector_result <- do.call(terralink_vector, vector_call)

print(vector_result)
plot(vector_result)

vector_prefix <- paste0(tools::file_path_sans_ext(basename(vector_patches_path)), "_v_tryout")
vector_outputs <- write_terralink_vector_outputs(
  result = vector_result,
  output_dir = vector_run_dir,
  prefix = vector_prefix,
  overwrite = TRUE,
  output_paths = list()
)

cat("Vector summary:\n")
print(vector_result$summary)
cat("\nVector output files:\n")
print(vector_outputs)

vector_run <- list(result = vector_result, outputs = vector_outputs)


# How to run in RStudio ----

# 1) Run Setup once.
# 2) Highlight and run the entire Raster block (top line to bottom line of that block).
# 3) Highlight and run the entire Vector block (top line to bottom line of that block).
