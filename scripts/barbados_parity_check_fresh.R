if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("/Users/benbishop/projects/terralink")) {
  pkgload::load_all("/Users/benbishop/projects/terralink", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(terralink)
  library(terra)
})

input_tif <- "/Users/benbishop/Downloads/smallest_barbados.tif"
reference_tif <- "/Users/benbishop/Downloads/correctoutput.tif"

if (!file.exists(input_tif)) stop("Input raster not found: ", input_tif, call. = FALSE)
if (!file.exists(reference_tif)) stop("Reference raster not found: ", reference_tif, call. = FALSE)

result <- terralink_raster(
  raster = input_tif,
  patch_values = 12,
  budget = 50,
  min_patch_size = 2,
  min_corridor_width = 3,
  max_search_distance = 600,
  units = "pixels",
  allow_bottlenecks = FALSE,
  strategy = "circuit_utility",
  verbose = 0,
  progress = FALSE
)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(tempdir(), paste0("barbados_parity_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_tif <- file.path(out_dir, "r_output_corridor.tif")
writeRaster(result$corridor_raster, out_tif, overwrite = TRUE)

ref <- rast(reference_tif)
out <- rast(out_tif)

rv <- values(ref, mat = FALSE)
ov <- values(out, mat = FALSE)

same <- sum((is.na(rv) & is.na(ov)) | (!is.na(rv) & !is.na(ov) & rv == ov))
total <- length(rv)
diff <- total - same
max_abs_diff <- max(abs(ifelse(is.na(rv), 0, rv) - ifelse(is.na(ov), 0, ov)))

ref_on <- which(!is.na(rv) & rv > 0)
out_on <- which(!is.na(ov) & ov > 0)
inter <- length(intersect(ref_on, out_on))
uni <- length(union(ref_on, out_on))
jaccard <- if (uni == 0) NA_real_ else inter / uni

cat("Summary:\n")
cat("  corridors_used =", result$summary$corridors_used, "\n")
cat("  budget_used    =", result$summary$budget_used, "\n")
cat("  candidate_edges=", result$summary$candidate_edges, "\n")
cat("Comparison:\n")
cat("  exact_equal_cells =", same, "of", total, "\n")
cat("  diff_cells        =", diff, "\n")
cat("  max_abs_diff      =", max_abs_diff, "\n")
cat("  jaccard           =", sprintf("%.6f", jaccard), "\n")
cat("Output:\n")
cat("  output_tif =", out_tif, "\n")
cat("  reference  =", reference_tif, "\n")

if (diff == 0) {
  cat("RESULT: PASS (output matches reference exactly)\n")
} else {
  cat("RESULT: FAIL (output does not match reference)\n")
}
