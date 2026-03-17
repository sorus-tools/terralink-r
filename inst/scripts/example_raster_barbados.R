# TerraLink raster example (Barbados-style settings)
#
# Usage:
#   source(system.file("scripts", "example_raster_barbados.R", package = "terralink"))
#
# By default this uses packaged sample data. Replace `raster_path` if needed.

suppressPackageStartupMessages({
  library(terralink)
  library(terra)
})

raster_path <- terralink_sample_data("raster")
output_dir <- file.path(tempdir(), "terralink_raster_example")

if (!file.exists(raster_path)) {
  stop("No sample raster found. Set `raster_path` to an existing .tif file before running this example.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result <- terralink_raster(
  raster = raster_path,
  patch_values = 1,
  obstacle_values = 2,
  budget = 220,
  min_patch_size = 10,
  min_corridor_width = 3,
  max_search_distance = 30,
  units = "pixels",
  allow_bottlenecks = FALSE,
  strategy = "most_connected_habitat",
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
