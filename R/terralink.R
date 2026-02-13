#' Run TerraLink optimization on abstract nodes and edges
#'
#' @param nodes Named numeric vector of patch weights.
#' @param edges Data frame with columns u, v, id, cost.
#' @param budget Numeric budget for corridor costs.
#' @param loop_fraction Fraction of budget available for loops.
#' @param max_redundancy Maximum redundant edges per component.
#' @return List with selected edges and component summaries.
#' @export
terralink_engine <- function(nodes, edges, budget, loop_fraction = 0.05, max_redundancy = 2) {
  optimize_network(nodes = nodes, edges = edges, budget = budget, loop_fraction = loop_fraction, max_redundancy = max_redundancy)
}

#' Run TerraLink corridor analysis on raster data
#'
#' @param raster SpatRaster or raster path.
#' @param patch_values Vector of habitat cell values.
#' @param patch_ranges Optional list of habitat ranges.
#' @param budget Total corridor budget (units defined by `units`).
#' @param budget_pixels Back-compat alias for budget (pixels).
#' @param strategy Either "largest_network", "circuit_utility", or "most_connectivity" (alias).
#' @param min_patch_size Minimum patch size (units defined by `units`).
#' @param min_corridor_width Minimum corridor width (units defined by `units`).
#' @param max_search_distance Maximum search distance (units defined by `units`).
#' @param obstacle_values Optional raster values that block corridors.
#' @param obstacle_ranges Optional list of impassable ranges.
#' @param allow_bottlenecks Whether to allow bottlenecks.
#' @param patch_connectivity Connectivity for patch labeling (4 or 8).
#' @param units "pixels", "metric", or "imperial".
#' @param allow_large Allow processing very large rasters.
#' @param max_pair_checks Limit for candidate pair checks.
#' @param max_candidates Limit for candidate corridors.
#' @param verbose Verbosity level (0-2).
#' @param progress Show progress bars.
#' @param obstacle_strategy Behavior when gdistance is unavailable and obstacles are provided.
#' @param output_dir Optional output directory.
#' @param output_prefix Optional output name prefix.
#' @param output_paths Optional named list of explicit output file paths.
#' @param write_outputs Whether to write outputs to disk.
#' @param keep_candidates Keep candidate list in output.
#' @return List containing corridors, metrics, and output layers.
#' @examples
#' r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
#' terra::values(r) <- 0
#' terra::values(r)[1:10] <- 1
#' terra::values(r)[300:320] <- 1
#'
#' result <- terralink_raster(
#'   raster = r,
#'   patch_values = 1,
#'   budget = 200,
#'   min_patch_size = 3,
#'   min_corridor_width = 3,
#'   max_search_distance = 50,
#'   units = "pixels"
#' )
#' result$summary
#' @export
terralink_raster <- function(
  raster,
  patch_values = NULL,
  patch_ranges = NULL,
  budget = NULL,
  budget_pixels = NULL,
  strategy = "circuit_utility",
  min_patch_size = 10,
  min_corridor_width = 3,
  max_search_distance = 100,
  obstacle_values = NULL,
  obstacle_ranges = NULL,
  allow_bottlenecks = FALSE,
  patch_connectivity = 8,
  units = "pixels",
  allow_large = FALSE,
  max_pair_checks = 2000000,
  max_candidates = 200000,
  verbose = 0,
  progress = FALSE,
  obstacle_strategy = c("error", "straight_line", "disable_obstacles"),
  output_dir = NULL,
  output_prefix = NULL,
  output_paths = NULL,
  write_outputs = FALSE,
  keep_candidates = FALSE
) {
  if (is.null(budget)) {
    if (is.null(budget_pixels)) {
      terralink_abort("budget must be provided.", class = "terralink_error_input")
    }
    if (!missing(units) && units != "pixels") {
      terralink_abort(
        "budget_pixels provided with non-pixel units.",
        class = "terralink_error_scale",
        fix = c("Set units = 'pixels'", "Provide budget in metric/imperial units")
      )
    }
    units <- "pixels"
    budget <- budget_pixels
  }
  result <- run_raster_analysis(
    raster = raster,
    patch_values = patch_values,
    patch_ranges = patch_ranges,
    budget = budget,
    strategy = strategy,
    min_patch_size = min_patch_size,
    min_corridor_width = min_corridor_width,
    max_search_distance = max_search_distance,
    obstacle_values = obstacle_values,
    obstacle_ranges = obstacle_ranges,
    allow_bottlenecks = allow_bottlenecks,
    patch_connectivity = patch_connectivity,
    units = units,
    allow_large = allow_large,
    max_pair_checks = max_pair_checks,
    max_candidates = max_candidates,
    verbose = verbose,
    progress = progress,
    obstacle_strategy = obstacle_strategy,
    keep_candidates = keep_candidates
  )
  if (isTRUE(write_outputs)) {
    if (is.null(output_dir) || !nzchar(output_dir)) {
      output_dir <- terralink_default_output_dir(raster)
    }
    result$output_paths <- write_terralink_raster_outputs(
      result,
      output_dir,
      output_prefix,
      output_paths = output_paths %||% list()
    )
  }
  result
}

#' Run TerraLink corridor analysis on vector patches
#'
#' @param patches sf object with polygon geometry (one feature per patch) or path.
#' @param budget Total corridor area budget (ha/ac).
#' @param strategy Either "largest_network", "circuit_utility", or "most_connectivity" (alias).
#' @param min_patch_size Minimum patch area (ha/ac).
#' @param min_corridor_width Minimum corridor width (m/ft).
#' @param max_search_distance Maximum search distance (m/ft).
#' @param obstacle_layers Optional obstacle layers (sf or path).
#' @param obstacle_resolution Raster resolution for obstacle routing.
#' @param units Unit system (metric or imperial).
#' @param max_pair_checks Limit for candidate pair checks.
#' @param max_candidates Limit for candidate corridors.
#' @param verbose Verbosity level (0-2).
#' @param progress Show progress bars.
#' @param obstacle_strategy Behavior when gdistance is unavailable and obstacles are provided.
#' @param return_crs CRS for outputs ("input" or "utm").
#' @param output_dir Optional output directory.
#' @param output_prefix Optional output name prefix.
#' @param output_paths Optional named list of explicit output file paths.
#' @param write_outputs Whether to write outputs to disk.
#' @param keep_candidates Keep candidate list in output.
#' @return List containing corridors sf, networks sf, and metrics.
#' @examples
#' p1 <- sf::st_polygon(list(rbind(c(0, 0), c(0, 10), c(10, 10), c(10, 0), c(0, 0))))
#' p2 <- sf::st_polygon(list(rbind(c(30, 0), c(30, 10), c(40, 10), c(40, 0), c(30, 0))))
#' patches <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(p1, p2), crs = 32618)
#'
#' result <- terralink_vector(
#'   patches = patches,
#'   budget = 1,
#'   min_patch_size = 0.001,
#'   min_corridor_width = 5,
#'   max_search_distance = 200,
#'   units = "metric"
#' )
#' result$summary
#' @export
terralink_vector <- function(
  patches,
  budget,
  strategy = "circuit_utility",
  min_patch_size = NULL,
  min_corridor_width = 100,
  max_search_distance = 5000,
  obstacle_layers = NULL,
  obstacle_resolution = NULL,
  units = "metric",
  max_pair_checks = 2000000,
  max_candidates = 200000,
  verbose = 0,
  progress = FALSE,
  obstacle_strategy = c("error", "straight_line", "disable_obstacles"),
  return_crs = c("input", "utm"),
  output_dir = NULL,
  output_prefix = NULL,
  output_paths = NULL,
  write_outputs = FALSE,
  keep_candidates = FALSE
) {
  result <- run_vector_analysis(
    patches = patches,
    budget = budget,
    strategy = strategy,
    min_patch_size = min_patch_size,
    min_corridor_width = min_corridor_width,
    max_search_distance = max_search_distance,
    obstacle_layers = obstacle_layers,
    obstacle_resolution = obstacle_resolution,
    units = units,
    max_pair_checks = max_pair_checks,
    max_candidates = max_candidates,
    verbose = verbose,
    progress = progress,
    obstacle_strategy = obstacle_strategy,
    return_crs = return_crs,
    keep_candidates = keep_candidates
  )
  if (isTRUE(write_outputs)) {
    if (is.null(output_dir) || !nzchar(output_dir)) {
      output_dir <- terralink_default_output_dir(patches)
    }
    result$output_paths <- write_terralink_vector_outputs(
      result,
      output_dir,
      output_prefix,
      output_paths = output_paths %||% list()
    )
  }
  result
}

#' Run TerraLink with a single entry point
#'
#' @param mode "raster" or "vector".
#' @param input Raster path/SpatRaster or sf/path.
#' @param ... Parameters forwarded to terralink_raster or terralink_vector.
#' @return Result list.
#' @export
terralink_run <- function(mode = c("raster", "vector"), input, ...) {
  mode <- match.arg(mode)
  if (mode == "raster") {
    terralink_raster(raster = input, ...)
  } else {
    terralink_vector(patches = input, ...)
  }
}
