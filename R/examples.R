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
