# TerraLink raster example (Barbados-style settings)
#
# Usage:
#   source(system.file("scripts", "example_raster_barbados.R", package = "terralink"))
#
# Edit `raster_path` and `output_dir` below, then run this whole script.

suppressPackageStartupMessages({
  library(terralink)
  library(terra)
})

raster_path <- "/path/to/smallest_barbados.tif"
output_dir <- file.path(tempdir(), "terralink_raster_example")

if (!file.exists(raster_path)) {
  stop("Set `raster_path` to an existing .tif file before running this example.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result <- terralink_raster(
  raster = raster_path,
  patch_values = 12,
  budget = 50,
  min_patch_size = 2,
  min_corridor_width = 3,
  max_search_distance = 600,
  units = "pixels",
  allow_bottlenecks = FALSE,
  strategy = "circuit_utility",
  progress = TRUE
)

print(result$summary)
plot(result)

outputs <- write_terralink_raster_outputs(
  result = result,
  output_dir = output_dir,
  prefix = "barbados_raster_example",
  overwrite = TRUE
)
print(outputs)
