if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("/Users/benbishop/projects/terralink")) {
  pkgload::load_all("/Users/benbishop/projects/terralink", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(terralink)
  library(sf)
})

input_shp <- "/Users/benbishop/Downloads/small barbados vector.shp"
reference_gpkg <- "/Users/benbishop/Downloads/correctoutput_vector.gpkg"
reference_layer <- "corridors_most_connectivity"

if (!file.exists(input_shp)) stop("Input vector not found: ", input_shp, call. = FALSE)
if (!file.exists(reference_gpkg)) stop("Reference GPKG not found: ", reference_gpkg, call. = FALSE)

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
out_dir <- file.path(tempdir(), paste0("barbados_vector_parity_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- write_terralink_vector_outputs(
  result = result,
  output_dir = out_dir,
  prefix = "r_vector_parity",
  overwrite = TRUE,
  output_paths = list()
)

out_layers <- st_layers(outputs$corridors)$name
out_layer <- if ("corridors" %in% out_layers) "corridors" else out_layers[[1]]
ref <- st_read(reference_gpkg, layer = reference_layer, quiet = TRUE)
out <- st_read(outputs$corridors, layer = out_layer, quiet = TRUE)

ref_u <- st_union(st_make_valid(ref))
out_u <- st_union(st_make_valid(out))
area_ref <- as.numeric(st_area(ref_u))
area_out <- as.numeric(st_area(out_u))
area_inter <- as.numeric(st_area(st_intersection(ref_u, out_u)))
area_union <- as.numeric(st_area(st_union(ref_u, out_u)))
area_symdiff <- as.numeric(st_area(st_sym_difference(ref_u, out_u)))
iou <- if (area_union > 0) area_inter / area_union else NA_real_

cat("Summary:\n")
cat("  corridors_used = ", result$summary$corridors_used, "\n", sep = "")
cat("  budget_used    = ", result$summary$budget_used, "\n", sep = "")
cat("  budget_total   = ", result$summary$budget_total, "\n", sep = "")
cat("  candidate_edges= ", result$summary$candidate_edges, "\n", sep = "")

cat("Geometry comparison:\n")
cat("  ref_features   = ", nrow(ref), "\n", sep = "")
cat("  out_features   = ", nrow(out), "\n", sep = "")
cat("  area_ref_m2    = ", area_ref, "\n", sep = "")
cat("  area_out_m2    = ", area_out, "\n", sep = "")
cat("  intersection_m2= ", area_inter, "\n", sep = "")
cat("  union_m2       = ", area_union, "\n", sep = "")
cat("  symdiff_m2     = ", area_symdiff, "\n", sep = "")
cat("  IoU            = ", sprintf("%.8f", iou), "\n", sep = "")

cat("Outputs:\n")
cat("  output_gpkg = ", outputs$corridors, "\n", sep = "")
cat("  reference   = ", reference_gpkg, " (layer ", reference_layer, ")\n", sep = "")
