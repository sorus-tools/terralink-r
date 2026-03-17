#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(sf)
  library(terra)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Unable to determine script path.")
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

parse_args <- function(args) {
  out <- list()
  idx <- 1L
  while (idx <= length(args)) {
    key <- args[[idx]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    value <- TRUE
    if (idx < length(args) && !startsWith(args[[idx + 1L]], "--")) {
      value <- args[[idx + 1L]]
      idx <- idx + 1L
    }
    out[[sub("^--", "", key)]] <- value
    idx <- idx + 1L
  }
  out
}

safe_name <- function(x) {
  gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
}

copy_input <- function(src, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(src, dest, overwrite = TRUE)) {
    stop(sprintf("Failed to copy %s to %s", src, dest))
  }
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

write_vector_input <- function(x, dest, layer) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest)) {
    file.remove(dest)
  }
  sf::st_write(x, dest, layer = layer, quiet = TRUE)
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

build_review_obstacle <- function(patches) {
  bbox <- sf::st_bbox(patches)
  width <- unname(bbox[["xmax"]] - bbox[["xmin"]])
  height <- unname(bbox[["ymax"]] - bbox[["ymin"]])

  xmin <- unname(bbox[["xmin"]] + width * 0.48)
  xmax <- unname(bbox[["xmin"]] + width * 0.56)
  ymin <- unname(bbox[["ymin"]] + height * 0.28)
  ymax <- unname(bbox[["ymin"]] + height * 0.73)

  geom <- sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(xmin, ymin),
      c(xmax, ymin),
      c(xmax, ymax),
      c(xmin, ymax),
      c(xmin, ymin)
    ))),
    crs = sf::st_crs(patches)
  )

  sf::st_sf(
    obstacle = "review_barrier",
    geometry = geom
  )
}

build_review_raster <- function(patches, obstacles, dest, resolution = 25) {
  combined_bbox <- sf::st_bbox(sf::st_union(
    sf::st_geometry(patches),
    sf::st_geometry(obstacles)
  ))
  crs_wkt <- sf::st_crs(patches)$wkt %||% sf::st_crs(patches)$input
  base <- terra::rast(
    xmin = unname(combined_bbox[["xmin"]] - resolution),
    xmax = unname(combined_bbox[["xmax"]] + resolution),
    ymin = unname(combined_bbox[["ymin"]] - resolution),
    ymax = unname(combined_bbox[["ymax"]] + resolution),
    resolution = resolution,
    crs = crs_wkt
  )
  values(base) <- 0

  habitat <- terra::rasterize(terra::vect(patches), base, field = 1, background = 0)
  barrier <- terra::rasterize(terra::vect(obstacles), base, field = 9, background = NA)
  output <- habitat
  output[!is.na(barrier)] <- 9

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(output, dest, overwrite = TRUE, datatype = "INT1U")
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

write_summary_csv <- function(path, result) {
  summary_df <- data.frame(
    item = names(result$summary),
    value = vapply(result$summary, function(x) paste(as.character(x), collapse = ", "), character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary_df, path, row.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_metrics_txt <- function(path, result) {
  writeLines(as.character(result$metrics_report %||% character(0)), con = path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

record_stats <- function(result) {
  metrics <- result$metrics %||% list()
  summary <- result$summary %||% list()
  list(
    corridors_used = as.integer(summary$corridors_used %||% 0L),
    budget_used = as.numeric(summary$budget_used %||% 0),
    total_connected = as.numeric(metrics$total_connected_habitat_area_post %||% 0),
    largest_network = as.numeric(metrics$largest_network_area_post %||% 0),
    habitat_availability_post = as.numeric(metrics$habitat_availability_post %||% NA_real_),
    mean_reachable_area_post = as.numeric(metrics$mean_reachable_area_post %||% NA_real_),
    largest_reachable_habitat_cluster_post = as.numeric(metrics$largest_reachable_habitat_cluster_post %||% NA_real_),
    mean_effective_resistance_post = as.numeric(metrics$mean_effective_resistance_post %||% NA_real_),
    landscape_fluidity_post = as.numeric(metrics$landscape_fluidity_post %||% NA_real_)
  )
}

raster_area_scale_ha <- function(raster_path) {
  if (is.null(raster_path) || !file.exists(raster_path)) {
    return(1)
  }
  rast <- terra::rast(raster_path)
  res_xy <- terra::res(rast)
  if (length(res_xy) < 2 || any(!is.finite(res_xy[1:2]))) {
    return(1)
  }
  abs(as.numeric(res_xy[[1]]) * as.numeric(res_xy[[2]])) / 10000
}

normalize_review_stats <- function(stats, input_type, raster_scale_ha = 1) {
  if (!identical(input_type, "raster") || !is.finite(raster_scale_ha) || raster_scale_ha <= 0) {
    return(stats)
  }
  area_keys <- c(
    "budget_used",
    "total_connected",
    "largest_network",
    "habitat_availability_post",
    "mean_reachable_area_post",
    "largest_reachable_habitat_cluster_post"
  )
  for (key in area_keys) {
    value <- stats[[key]] %||% NA_real_
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      stats[[key]] <- as.numeric(value) * raster_scale_ha
    }
  }
  stats
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(file.path(dirname(script_path()), ".."), winslash = "/", mustWork = TRUE)
  devtools::load_all(repo_root, quiet = TRUE)

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  root_dir <- normalizePath(
    args$root %||% file.path(repo_root, "tmp", paste0("dummy_visual_review_", timestamp)),
    winslash = "/",
    mustWork = FALSE
  )
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
  inputs_dir <- file.path(root_dir, "inputs")
  r_runs_dir <- file.path(root_dir, "r_runs")
  dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(r_runs_dir, recursive = TRUE, showWarnings = FALSE)

  input_vector_src <- terralink_sample_data("synthetic_vector")

  input_vector_path <- copy_input(input_vector_src, file.path(inputs_dir, "dummy_vector_patches.gpkg"))
  input_patches <- sf::st_read(input_vector_path, quiet = TRUE)
  input_obstacles <- build_review_obstacle(input_patches)
  input_obstacle_path <- write_vector_input(
    input_obstacles,
    file.path(inputs_dir, "dummy_vector_impassable.gpkg"),
    "dummy_vector_impassable"
  )
  input_raster_path <- build_review_raster(
    input_patches,
    input_obstacles,
    file.path(inputs_dir, "dummy_raster.tif")
  )
  raster_scale_ha <- raster_area_scale_ha(input_raster_path)

  vector_layer_name <- sf::st_layers(input_vector_path)$name[[1]]
  obstacle_layer_name <- sf::st_layers(input_obstacle_path)$name[[1]]

  strategies <- c(
    "most_connected_habitat",
    "largest_single_network",
    "reachable_habitat_advanced",
    "landscape_fluidity"
  )
  obstacle_modes <- c("no_obstacles", "with_obstacles")

  vector_params <- list(
    budget = 18,
    min_patch_size = 0.01,
    min_corridor_width = 60,
    max_search_distance = 900,
    units = "metric",
    species_dispersal_distance = 800,
    species_dispersal_kernel = "exponential",
    patch_quality_field = "quality",
    patch_area_scaling = "sqrt"
  )
  raster_params <- list(
    patch_values = c(1),
    obstacle_values = c(9),
    budget = 120,
    min_patch_size = 4,
    min_corridor_width = 1,
    max_search_distance = 30,
    units = "pixels",
    species_dispersal_distance = 12,
    species_dispersal_kernel = "exponential",
    patch_area_scaling = "sqrt"
  )

  run_records <- list()
  summary_rows <- list()

  for (input_type in c("vector", "raster")) {
    for (obstacle_mode in obstacle_modes) {
      for (strategy in strategies) {
        run_key <- safe_name(paste("r", input_type, obstacle_mode, strategy, sep = "__"))
        run_dir <- file.path(r_runs_dir, run_key)
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

        if (identical(input_type, "vector")) {
          result <- terralink_vector(
            patches = input_vector_path,
            obstacle_layers = if (identical(obstacle_mode, "with_obstacles")) input_obstacle_path else NULL,
            budget = vector_params$budget,
            min_patch_size = vector_params$min_patch_size,
            min_corridor_width = vector_params$min_corridor_width,
            max_search_distance = vector_params$max_search_distance,
            units = vector_params$units,
            strategy = strategy,
            species_dispersal_distance = if (identical(strategy, "reachable_habitat_advanced")) vector_params$species_dispersal_distance else NULL,
            species_dispersal_kernel = vector_params$species_dispersal_kernel,
            patch_quality_field = if (identical(strategy, "reachable_habitat_advanced")) vector_params$patch_quality_field else NULL,
            patch_area_scaling = vector_params$patch_area_scaling
          )
          gpkg_path <- file.path(run_dir, paste0(run_key, ".gpkg"))
          if (inherits(result$corridors, "sf") && nrow(result$corridors) > 0) {
            sf::st_write(result$corridors, gpkg_path, layer = "corridors", delete_layer = TRUE, quiet = TRUE)
          }
          if (inherits(result$networks, "sf") && nrow(result$networks) > 0) {
            sf::st_write(result$networks, gpkg_path, layer = "networks", delete_layer = TRUE, quiet = TRUE)
          }
          summary_path <- write_summary_csv(file.path(run_dir, paste0(run_key, "_summary.csv")), result)
          metrics_path <- write_metrics_txt(file.path(run_dir, paste0(run_key, "_metrics.txt")), result)
          rds_path <- file.path(run_dir, paste0(run_key, ".rds"))
          saveRDS(result, rds_path)
          extra_paths <- list(
            gpkg_path = normalizePath(gpkg_path, winslash = "/", mustWork = TRUE),
            corridor_layer = "corridors",
            network_layer = "networks",
            summary_path = summary_path,
            metrics_path = metrics_path,
            result_rds_path = normalizePath(rds_path, winslash = "/", mustWork = TRUE)
          )
        } else {
          result <- terralink_raster(
            raster = input_raster_path,
            patch_values = raster_params$patch_values,
            obstacle_values = if (identical(obstacle_mode, "with_obstacles")) raster_params$obstacle_values else NULL,
            budget = raster_params$budget,
            min_patch_size = raster_params$min_patch_size,
            min_corridor_width = raster_params$min_corridor_width,
            max_search_distance = raster_params$max_search_distance,
            units = raster_params$units,
            strategy = strategy,
            species_dispersal_distance = if (identical(strategy, "reachable_habitat_advanced")) raster_params$species_dispersal_distance else NULL,
            species_dispersal_kernel = raster_params$species_dispersal_kernel,
            patch_area_scaling = raster_params$patch_area_scaling
          )
          gpkg_path <- file.path(run_dir, paste0(run_key, ".gpkg"))
          if (inherits(result$corridors, "sf") && nrow(result$corridors) > 0) {
            sf::st_write(result$corridors, gpkg_path, layer = "corridors", delete_layer = TRUE, quiet = TRUE)
          }
          if (inherits(result$networks, "sf") && nrow(result$networks) > 0) {
            sf::st_write(result$networks, gpkg_path, layer = "networks", delete_layer = TRUE, quiet = TRUE)
          }
          corridor_raster_path <- file.path(run_dir, paste0(run_key, "_corridors.tif"))
          contiguous_raster_path <- file.path(run_dir, paste0(run_key, "_contiguous.tif"))
          patches_raster_path <- file.path(run_dir, paste0(run_key, "_patches.tif"))
          if (!is.null(result$corridor_raster)) {
            terra::writeRaster(result$corridor_raster, corridor_raster_path, overwrite = TRUE)
          }
          if (!is.null(result$contiguous_raster)) {
            terra::writeRaster(result$contiguous_raster, contiguous_raster_path, overwrite = TRUE)
          }
          if (!is.null(result$patches)) {
            terra::writeRaster(result$patches, patches_raster_path, overwrite = TRUE)
          }
          summary_path <- write_summary_csv(file.path(run_dir, paste0(run_key, "_summary.csv")), result)
          metrics_path <- write_metrics_txt(file.path(run_dir, paste0(run_key, "_metrics.txt")), result)
          rds_path <- file.path(run_dir, paste0(run_key, ".rds"))
          saveRDS(result, rds_path)
          extra_paths <- list(
            gpkg_path = normalizePath(gpkg_path, winslash = "/", mustWork = TRUE),
            corridor_layer = "corridors",
            network_layer = "networks",
            corridor_raster_path = if (file.exists(corridor_raster_path)) normalizePath(corridor_raster_path, winslash = "/", mustWork = TRUE) else NULL,
            contiguous_raster_path = if (file.exists(contiguous_raster_path)) normalizePath(contiguous_raster_path, winslash = "/", mustWork = TRUE) else NULL,
            patches_raster_path = if (file.exists(patches_raster_path)) normalizePath(patches_raster_path, winslash = "/", mustWork = TRUE) else NULL,
            summary_path = summary_path,
            metrics_path = metrics_path,
            result_rds_path = normalizePath(rds_path, winslash = "/", mustWork = TRUE)
          )
        }

        stats <- normalize_review_stats(
          record_stats(result),
          input_type = input_type,
          raster_scale_ha = raster_scale_ha
        )
        run_records[[length(run_records) + 1L]] <- c(
          list(
            engine = "R",
            input_type = input_type,
            obstacle_mode = obstacle_mode,
            strategy = strategy,
            run_key = run_key,
            run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
            display_label = paste("R", input_type, obstacle_mode, strategy, sep = " | "),
            stats = stats
          ),
          extra_paths
        )
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          engine = "R",
          input_type = input_type,
          obstacle_mode = obstacle_mode,
          strategy = strategy,
          corridors_used = stats$corridors_used,
          budget_used = stats$budget_used,
          total_connected = stats$total_connected,
          largest_network = stats$largest_network,
          habitat_availability_post = stats$habitat_availability_post,
          mean_reachable_area_post = stats$mean_reachable_area_post,
          largest_reachable_habitat_cluster_post = stats$largest_reachable_habitat_cluster_post,
          mean_effective_resistance_post = stats$mean_effective_resistance_post,
          landscape_fluidity_post = stats$landscape_fluidity_post,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary_csv_path <- file.path(root_dir, "r_run_summary.csv")
  utils::write.csv(do.call(rbind, summary_rows), summary_csv_path, row.names = FALSE)

  manifest <- list(
    root_dir = normalizePath(root_dir, winslash = "/", mustWork = TRUE),
    group_name = paste("TerraLink Dummy Visual Review", timestamp),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    inputs = list(
      vector_patches_path = input_vector_path,
      vector_patches_layer = vector_layer_name,
      vector_obstacles_path = input_obstacle_path,
      vector_obstacles_layer = obstacle_layer_name,
      raster_path = input_raster_path
    ),
    strategies = as.list(strategies),
    params = list(vector = vector_params, raster = raster_params),
    r_runs = run_records,
    r_summary_csv = normalizePath(summary_csv_path, winslash = "/", mustWork = TRUE)
  )
  manifest_path <- file.path(root_dir, "review_manifest.json")
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

  cat(
    jsonlite::toJSON(
      list(
        root_dir = normalizePath(root_dir, winslash = "/", mustWork = TRUE),
        manifest_path = normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
      ),
      auto_unbox = TRUE
    ),
    "\n",
    sep = ""
  )
}

main()
