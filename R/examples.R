#' Locate packaged TerraLink example scripts
#'
#' @param type Which example scripts to return: `"all"`, `"raster"`, or `"vector"`.
#' @return Character vector of absolute file paths.
#' @export
#' @examples
#' terralink_examples()
#' terralink_examples("raster")
terralink_examples <- function(type = c("all", "raster", "vector")) {
  type <- match.arg(type)
  scripts_dir <- system.file("scripts", package = "terralink")
  if (!nzchar(scripts_dir) || !dir.exists(scripts_dir)) {
    return(character(0))
  }

  all_scripts <- c(
    file.path(scripts_dir, "example_raster_barbados.R"),
    file.path(scripts_dir, "example_vector_barbados.R"),
    file.path(scripts_dir, "run_raster.R"),
    file.path(scripts_dir, "run_vector.R")
  )
  all_scripts <- unique(all_scripts[file.exists(all_scripts)])

  if (type == "all") {
    return(all_scripts)
  }
  if (type == "raster") {
    return(all_scripts[grepl("raster", basename(all_scripts), ignore.case = TRUE)])
  }
  all_scripts[grepl("vector", basename(all_scripts), ignore.case = TRUE)]
}

#' Locate packaged TerraLink sample data files
#'
#' @param type Which sample data path to return: `"all"`, `"raster"`, `"vector"`,
#'   `"obstacle"`, `"synthetic_raster"`, `"synthetic_vector"`, or `"synthetic_obstacle"`.
#' @return For `"all"`, a named character vector of absolute file paths.
#'   For other values, a single absolute file path (character scalar) or `character(0)`
#'   when unavailable.
#' @export
#' @examples
#' terralink_sample_data()
#' terralink_sample_data("raster")
terralink_sample_data <- function(type = c("all", "raster", "vector", "obstacle", "synthetic_raster", "synthetic_vector", "synthetic_obstacle")) {
  type <- match.arg(type)

  files <- c(
    raster = system.file("extdata", "habitat.tif", package = "terralink"),
    vector = system.file("extdata", "patches.gpkg", package = "terralink"),
    obstacle = system.file("extdata", "impassable.gpkg", package = "terralink"),
    synthetic_raster = system.file("extdata", "synthetic_habitat.tif", package = "terralink"),
    synthetic_vector = system.file("extdata", "synthetic_patches.gpkg", package = "terralink"),
    synthetic_obstacle = system.file("extdata", "synthetic_impassable.gpkg", package = "terralink")
  )
  files <- files[nzchar(files) & file.exists(files)]

  if (type == "all") {
    return(files)
  }

  key <- switch(
    type,
    raster = "raster",
    vector = "vector",
    obstacle = "obstacle",
    synthetic_raster = "synthetic_raster",
    synthetic_vector = "synthetic_vector",
    synthetic_obstacle = "synthetic_obstacle"
  )
  out <- files[[key]]
  if (is.null(out)) character(0) else unname(out)
}
