# TerraLink vector example (Barbados-style settings)
#
# Usage:
#   source(system.file("scripts", "example_vector_barbados.R", package = "terralink"))
#
# By default this uses packaged sample data. Replace `patches_path` if needed.

suppressPackageStartupMessages({
  library(terralink)
  library(sf)
})

patches_path <- terralink_sample_data("vector")
output_dir <- file.path(tempdir(), "terralink_vector_example")

if (!file.exists(patches_path)) {
  stop("No sample vector found. Set `patches_path` to an existing vector file before running this example.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result <- terralink_vector(
  patches = patches_path,
  budget = 2.5,
  min_patch_size = 0.1,
  min_corridor_width = 20,
  max_search_distance = 500,
  units = "metric",
  return_crs = "input",
  strategy = "most_connected_habitat",
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
