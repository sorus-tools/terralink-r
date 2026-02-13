suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("/Users/benbishop/projects/terralink")) {
    pkgload::load_all("/Users/benbishop/projects/terralink", quiet = TRUE)
  }
  library(terralink)
  library(sf)
})

make_rect <- function(xmin, ymin, width, height) {
  st_polygon(list(rbind(
    c(xmin, ymin),
    c(xmin + width, ymin),
    c(xmin + width, ymin + height),
    c(xmin, ymin + height),
    c(xmin, ymin)
  )))
}

build_synthetic_patches <- function(crs = 32618) {
  polys <- list(
    make_rect(0, 0, 100, 100),
    make_rect(220, 10, 100, 100),
    make_rect(440, 0, 100, 100),
    make_rect(660, 25, 100, 100),
    make_rect(1000, 0, 100, 100),
    make_rect(1200, 10, 100, 100),
    make_rect(280, 320, 100, 100),
    make_rect(520, 320, 100, 100),
    make_rect(760, 320, 100, 100),
    make_rect(1320, 300, 100, 100)
  )
  st_sf(
    patch_id = seq_along(polys),
    geometry = st_sfc(polys, crs = crs)
  )
}

run_synthetic_validation <- function(output_root = tempdir()) {
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(output_root, paste0("terralink_vector_synth_", run_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  patches <- build_synthetic_patches()
  input_gpkg <- file.path(out_dir, "synthetic_patches.gpkg")
  st_write(patches, input_gpkg, layer = "patches", delete_dsn = TRUE, quiet = TRUE)

  cfg <- list(
    budget = 1.2,
    min_patch_size = 0.5,
    min_corridor_width = 30,
    max_search_distance = 300,
    units = "metric",
    strategy = "circuit_utility"
  )

  result <- terralink_vector(
    patches = input_gpkg,
    budget = cfg$budget,
    min_patch_size = cfg$min_patch_size,
    min_corridor_width = cfg$min_corridor_width,
    max_search_distance = cfg$max_search_distance,
    units = cfg$units,
    strategy = cfg$strategy,
    return_crs = "input",
    verbose = 1,
    progress = FALSE
  )

  outputs <- write_terralink_vector_outputs(
    result = result,
    output_dir = out_dir,
    prefix = "synthetic_vector",
    overwrite = TRUE,
    output_paths = list()
  )

  corridors <- result$corridors
  corridor_count <- if (inherits(corridors, "sf")) nrow(corridors) else 0L
  pre_components <- nrow(patches)
  post_components <- if (!is.null(result$networks) && inherits(result$networks, "sf")) nrow(result$networks) else pre_components

  overlap_m2 <- 0
  if (corridor_count > 0) {
    patch_union <- st_union(st_geometry(result$patches))
    overlap_geom <- suppressWarnings(tryCatch(st_intersection(st_union(st_geometry(corridors)), patch_union), error = function(e) NULL))
    if (!is.null(overlap_geom) && length(overlap_geom) > 0 && !isTRUE(any(st_is_empty(overlap_geom)))) {
      overlap_m2 <- as.numeric(st_area(overlap_geom))
    }
  }

  checks <- c(
    corridor_count_positive = corridor_count > 0,
    budget_not_exceeded = isTRUE(result$summary$budget_used <= result$summary$budget_total + 1e-9),
    corridor_geometry_valid = if (corridor_count > 0) all(st_is_valid(corridors)) else TRUE,
    corridor_patch_overlap_small = overlap_m2 < 1e-6,
    components_not_increasing = post_components <= pre_components
  )

  png(file.path(out_dir, "synthetic_result_plot.png"), width = 1200, height = 800)
  plot(st_geometry(result$patches), border = "grey30", col = "grey95", main = "Synthetic TerraLink Vector Validation")
  if (corridor_count > 0) {
    plot(st_geometry(corridors), add = TRUE, border = "red3", col = grDevices::adjustcolor("orange", alpha.f = 0.45))
  }
  dev.off()

  cat("Synthetic vector validation summary\n")
  cat("  Input:", input_gpkg, "\n")
  cat("  Output dir:", out_dir, "\n")
  cat("  Corridors used:", result$summary$corridors_used, "\n")
  cat("  Budget used:", sprintf("%.6f", result$summary$budget_used), "of", sprintf("%.6f", result$summary$budget_total), "\n")
  cat("  Components pre/post:", pre_components, "/", post_components, "\n")
  cat("  Corridor overlap with patches (m2):", sprintf("%.8f", overlap_m2), "\n")
  cat("Checks:\n")
  for (nm in names(checks)) {
    cat("  ", nm, "=", checks[[nm]], "\n", sep = "")
  }

  if (!all(checks)) {
    stop("Synthetic validation failed one or more checks.", call. = FALSE)
  }

  invisible(list(
    out_dir = out_dir,
    input_gpkg = input_gpkg,
    outputs = outputs,
    summary = result$summary,
    checks = checks
  ))
}

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else tempdir()
run_synthetic_validation(output_root = output_root)
