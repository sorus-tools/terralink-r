# TerraLink vector example (Barbados-style settings)
#
# Usage:
#   source(system.file("scripts", "example_vector_barbados.R", package = "terralink"))
#
# Edit `patches_path` and `output_dir` below, then run this whole script.

suppressPackageStartupMessages({
  library(terralink)
  library(sf)
})

patches_path <- "/path/to/small_barbados_vector.shp"
output_dir <- file.path(tempdir(), "terralink_vector_example")

if (!file.exists(patches_path)) {
  stop("Set `patches_path` to an existing vector file before running this example.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result <- terralink_vector(
  patches = patches_path,
  budget = 1,
  min_patch_size = 1,
  min_corridor_width = 20,
  max_search_distance = 5000,
  units = "metric",
  return_crs = "input",
  strategy = "circuit_utility",
  progress = TRUE
)

print(result$summary)
plot(result)

outputs <- write_terralink_vector_outputs(
  result = result,
  output_dir = output_dir,
  prefix = "barbados_vector_example",
  overwrite = TRUE
)
print(outputs)
