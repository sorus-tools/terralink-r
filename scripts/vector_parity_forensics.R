suppressPackageStartupMessages({
  library(sf)
})

if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("/Users/benbishop/projects/terralink")) {
  pkgload::load_all("/Users/benbishop/projects/terralink", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(terralink)
})

parse_pair_key <- function(x) {
  vals <- suppressWarnings(as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2) return(NA_character_)
  paste(min(vals), max(vals), sep = "-")
}

pair_key <- function(a, b) {
  paste(pmin(as.integer(a), as.integer(b)), pmax(as.integer(a), as.integer(b)), sep = "-")
}

args <- commandArgs(trailingOnly = TRUE)
input_shp <- if (length(args) >= 1) args[[1]] else "/Users/benbishop/Downloads/small barbados vector.shp"
reference_gpkg <- if (length(args) >= 2) args[[2]] else "/Users/benbishop/Downloads/correctoutput_vector.gpkg"
reference_layer <- if (length(args) >= 3) args[[3]] else "corridors_most_connectivity"

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

out <- result$corridors
ref <- st_read(reference_gpkg, layer = reference_layer, quiet = TRUE)
out$key <- pair_key(out$patch1, out$patch2)
ref$key <- vapply(ref$patch_ids, parse_pair_key, character(1))

if (!isTRUE(st_crs(out) == st_crs(ref))) {
  ref <- st_transform(ref, st_crs(out))
}

out_u <- st_union(st_make_valid(out))
ref_u <- st_union(st_make_valid(ref))
area_out <- as.numeric(st_area(out_u))
area_ref <- as.numeric(st_area(ref_u))
area_inter <- as.numeric(st_area(st_intersection(out_u, ref_u)))
area_union <- as.numeric(st_area(st_union(out_u, ref_u)))
iou <- if (area_union > 0) area_inter / area_union else NA_real_

cat("Global summary\n")
cat("  corridors_used = ", result$summary$corridors_used, "\n", sep = "")
cat("  budget_used    = ", result$summary$budget_used, "\n", sep = "")
cat("  area_out_m2    = ", area_out, "\n", sep = "")
cat("  area_ref_m2    = ", area_ref, "\n", sep = "")
cat("  global_iou     = ", sprintf("%.8f", iou), "\n\n", sep = "")

all_keys <- sort(unique(c(out$key, ref$key)))
rows <- vector("list", length(all_keys))
for (idx in seq_along(all_keys)) {
  k <- all_keys[[idx]]
  out_k <- out[out$key == k, , drop = FALSE]
  ref_k <- ref[ref$key == k, , drop = FALSE]
  if (nrow(out_k) == 0 || nrow(ref_k) == 0) {
    rows[[idx]] <- data.frame(
      key = k,
      area_out_m2 = if (nrow(out_k) > 0) as.numeric(st_area(st_union(st_make_valid(out_k)))) else 0,
      area_ref_m2 = if (nrow(ref_k) > 0) as.numeric(st_area(st_union(st_make_valid(ref_k)))) else 0,
      iou = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }
  go <- st_union(st_make_valid(out_k))
  gr <- st_union(st_make_valid(ref_k))
  ao <- as.numeric(st_area(go))
  ar <- as.numeric(st_area(gr))
  inter <- as.numeric(st_area(st_intersection(go, gr)))
  uni <- as.numeric(st_area(st_union(go, gr)))
  rows[[idx]] <- data.frame(
    key = k,
    area_out_m2 = ao,
    area_ref_m2 = ar,
    iou = if (uni > 0) inter / uni else NA_real_,
    stringsAsFactors = FALSE
  )
}

pair_df <- do.call(rbind, rows)
pair_df <- pair_df[order(pair_df$iou, na.last = TRUE), , drop = FALSE]

cat("Per-pair comparison (worst first)\n")
print(pair_df, row.names = FALSE)
