# Output writers for TerraLink

#' Write raster outputs to disk
#'
#' @param result Result list from terralink_raster.
#' @param output_dir Directory to write outputs.
#' @param prefix Optional name prefix for outputs.
#' @param output_paths Named list of explicit file paths to override defaults.
#' @param overwrite Whether to overwrite existing files.
#' @return Named list of written file paths.
#' @export
write_terralink_raster_outputs <- function(result, output_dir, prefix = NULL, overwrite = TRUE, output_paths = list()) {
  if (is.null(output_dir) || !nzchar(output_dir)) {
    terralink_abort("output_dir must be provided.", class = "terralink_error_input")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(prefix) || !nzchar(prefix)) {
    prefix <- terralink_safe_name(result$summary$input_name %||% "terralink")
  } else {
    prefix <- terralink_safe_name(prefix)
  }

  paths <- list()

  if (!is.null(result$corridor_raster)) {
    path <- output_paths$corridor_raster %||% file.path(output_dir, paste0("corridors_", prefix, ".tif"))
    terra::writeRaster(result$corridor_raster, path, overwrite = overwrite)
    paths$corridor_raster <- path
  }
  if (!is.null(result$contiguous_raster)) {
    path <- output_paths$contiguous_raster %||% file.path(output_dir, paste0("terralink_contiguous_", prefix, ".tif"))
    terra::writeRaster(result$contiguous_raster, path, overwrite = overwrite)
    paths$contiguous_raster <- path
  }
  if (!is.null(result$patches)) {
    path <- output_paths$patches %||% file.path(output_dir, paste0("patches_", prefix, ".tif"))
    terra::writeRaster(result$patches, path, overwrite = overwrite)
    paths$patches <- path
  }

  if (!is.null(result$patch_table)) {
    path <- output_paths$patch_table %||% file.path(output_dir, paste0("patches_", prefix, ".csv"))
    utils::write.csv(result$patch_table, path, row.names = FALSE)
    paths$patch_table <- path
  }
  if (!is.null(result$corridors)) {
    path <- output_paths$corridors %||% file.path(output_dir, paste0("corridors_", prefix, ".csv"))
    utils::write.csv(result$corridors, path, row.names = FALSE)
    paths$corridors <- path
  }
  if (!is.null(result$summary)) {
    path <- output_paths$summary %||% file.path(output_dir, paste0("terralink_raster_summary_", prefix, ".csv"))
    summary_df <- data.frame(
      item = names(result$summary),
      value = vapply(result$summary, function(x) paste(as.character(x), collapse = ", "), character(1))
    )
    utils::write.csv(summary_df, path, row.names = FALSE)
    paths$summary <- path
  }
  if (!is.null(result$metrics_report)) {
    path <- output_paths$metrics %||% file.path(output_dir, paste0("landscape_metrics_", prefix, ".txt"))
    writeLines(result$metrics_report, con = path)
    paths$metrics <- path
  }

  rds_path <- output_paths$result %||% file.path(output_dir, paste0("terralink_result_", prefix, ".rds"))
  saveRDS(result, rds_path)
  paths$result <- rds_path

  paths
}

#' Write vector outputs to disk
#'
#' @param result Result list from terralink_vector.
#' @param output_dir Directory to write outputs.
#' @param prefix Optional name prefix for outputs.
#' @param output_paths Named list of explicit file paths to override defaults.
#' @param overwrite Whether to overwrite existing files.
#' @return Named list of written file paths.
#' @export
write_terralink_vector_outputs <- function(result, output_dir, prefix = NULL, overwrite = TRUE, output_paths = list()) {
  if (is.null(output_dir) || !nzchar(output_dir)) {
    terralink_abort("output_dir must be provided.", class = "terralink_error_input")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(prefix) || !nzchar(prefix)) {
    prefix <- terralink_safe_name(result$summary$input_name %||% "terralink")
  } else {
    prefix <- terralink_safe_name(prefix)
  }

  paths <- list()
  gpkg_path <- output_paths$gpkg %||% file.path(output_dir, paste0("terralink_vector_", prefix, ".gpkg"))
  if (!is.null(result$corridors)) {
    sf::st_write(result$corridors, gpkg_path, layer = "corridors", delete_layer = TRUE, quiet = TRUE)
    paths$corridors <- gpkg_path
  }
  if (!is.null(result$networks)) {
    sf::st_write(result$networks, gpkg_path, layer = "networks", delete_layer = TRUE, quiet = TRUE)
    paths$networks <- gpkg_path
  }

  if (!is.null(result$summary)) {
    path <- output_paths$summary %||% file.path(output_dir, paste0("terralink_vector_summary_", prefix, ".csv"))
    summary_df <- data.frame(
      item = names(result$summary),
      value = vapply(result$summary, function(x) paste(as.character(x), collapse = ", "), character(1))
    )
    utils::write.csv(summary_df, path, row.names = FALSE)
    paths$summary <- path
  }
  if (!is.null(result$metrics_report)) {
    path <- output_paths$metrics %||% file.path(output_dir, paste0("landscape_metrics_", prefix, ".txt"))
    writeLines(result$metrics_report, con = path)
    paths$metrics <- path
  }

  rds_path <- output_paths$result %||% file.path(output_dir, paste0("terralink_result_", prefix, ".rds"))
  saveRDS(result, rds_path)
  paths$result <- rds_path

  paths
}
