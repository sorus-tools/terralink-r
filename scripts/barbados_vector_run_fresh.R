if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("/Users/benbishop/projects/terralink")) {
  pkgload::load_all("/Users/benbishop/projects/terralink", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(terralink)
  library(sf)
})

input_shp <- "/Users/benbishop/Downloads/small barbados vector.shp"
if (!file.exists(input_shp)) stop("Input vector not found: ", input_shp, call. = FALSE)

result <- terralink_vector(
  patches = input_shp,
  budget = 1,
  min_patch_size = 1,
  min_corridor_width = 20,
  max_search_distance = 5000,
  units = "metric",
  return_crs = "input",
  strategy = "circuit_utility",
  verbose = 0,
  progress = FALSE
)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(tempdir(), paste0("barbados_vector_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- write_terralink_vector_outputs(
  result = result,
  output_dir = out_dir,
  prefix = "barbados_vector_r_output",
  overwrite = TRUE,
  output_paths = list()
)

cat("Summary:\n")
cat("  corridors_used = ", result$summary$corridors_used, "\n", sep = "")
cat("  budget_used    = ", result$summary$budget_used, "\n", sep = "")
cat("  budget_total   = ", result$summary$budget_total, "\n", sep = "")
cat("  candidate_edges= ", result$summary$candidate_edges, "\n", sep = "")
cat("  candidate_pairs= ", result$summary$candidate_pairs, "\n", sep = "")
cat("  possible_pairs = ", result$summary$possible_pairs, "\n", sep = "")
cat("  patches        = ", result$summary$patches, "\n", sep = "")

cat("Outputs:\n")
print(outputs)
