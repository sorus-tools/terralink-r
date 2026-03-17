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
#' Identifies habitat patches from a classified raster, builds candidate
#' corridors between nearby patches, and selects an optimal corridor network
#' under a budget constraint. Raster inputs are polygonized and routed through
#' TerraLink's vector engine, matching the current QGIS plugin workflow.
#'
#' @param raster SpatRaster or file path to a single-band raster.
#' @param patch_values Integer vector of cell values that represent habitat
#'   (e.g., \code{c(1, 3)}). At least one of \code{patch_values} or
#'   \code{patch_ranges} must be provided.
#' @param patch_ranges Optional list of length-2 numeric vectors giving
#'   inclusive value ranges that define habitat (e.g., \code{list(c(1, 3))}).
#' @param budget Total corridor budget. When \code{units = "pixels"} this is
#'   the number of corridor cells allowed. When \code{units = "metric"} or
#'   \code{"imperial"}, this is the total corridor \strong{area} in hectares or
#'   acres. A reasonable starting point is 5--20 percent of total habitat area.
#' @param budget_pixels Back-compat alias for \code{budget} in pixel units.
#'   Use \code{budget} instead.
#' @param strategy Character string selecting the optimization objective.
#'   One of \code{"most_connected_habitat"} (default; maximizes total
#'   structurally connected habitat area), \code{"largest_single_network"}
#'   (maximizes the single largest connected component),
#'   \code{"landscape_fluidity"} (maximizes ease of movement and route
#'   redundancy), or \code{"reachable_habitat_advanced"} (maximizes
#'   dispersal-weighted reachable habitat; set
#'   \code{species_dispersal_distance}).
#' @param min_patch_size Numeric. Minimum patch size to include. In pixel
#'   units when \code{units = "pixels"} (number of cells), in hectares when
#'   \code{units = "metric"}, or acres when \code{units = "imperial"}.
#'   Patches smaller than this are dropped before analysis. Default: 10.
#' @param min_corridor_width Numeric. Minimum corridor width. In pixel units
#'   when \code{units = "pixels"} (cell widths), in meters when
#'   \code{units = "metric"}, or feet when \code{units = "imperial"}.
#'   Controls the buffer applied to corridor center-lines. Default: 3.
#' @param corridor_cell_assignment Character string controlling how corridor
#'   cells are valued in the output raster. One of
#'   \code{"sum_total_network_area"} (default; total network area),
#'   \code{"sum_direct_connected_patches"} (area of the two directly linked
#'   patches), or \code{"efficiency"} (cost-efficiency score).
#' @param max_search_distance Numeric. Maximum distance between patch edges
#'   to consider a candidate corridor. Same unit system as
#'   \code{min_corridor_width} (pixels / meters / feet). Increase this if few
#'   or no corridors are found. Default: 100.
#' @param obstacle_values Optional integer vector of raster cell values that
#'   represent impassable barriers (e.g., roads, water bodies).
#' @param obstacle_ranges Optional list of length-2 numeric vectors giving
#'   inclusive value ranges for obstacles.
#' @param allow_bottlenecks Logical. If \code{TRUE}, corridors narrower than
#'   \code{min_corridor_width} are still allowed when no wider path exists.
#'   Default: \code{FALSE}.
#' @param patch_connectivity Integer, 4 or 8. Pixel connectivity rule for
#'   grouping habitat cells into patches. 8 (default) includes diagonal
#'   neighbors; 4 uses only cardinal neighbors.
#' @param units Character string specifying the unit system:
#'   \code{"pixels"} (raster cell units), \code{"metric"} (hectares for area,
#'   meters for distance), or \code{"imperial"} (acres for area, feet for
#'   distance). Default: \code{"pixels"}.
#' @param allow_large Logical. Set to \code{TRUE} to allow processing rasters
#'   with more than 10 million cells. Default: \code{FALSE}.
#' @param max_pair_checks Integer. Upper limit on the number of patch pairs
#'   evaluated during candidate generation. Increase for landscapes with many
#'   patches; decrease if running out of memory. Default: 2,000,000.
#' @param max_candidates Integer. Upper limit on total candidate corridors
#'   retained. Default: 200,000.
#' @param verbose Integer verbosity level: 0 = silent, 1 = progress messages,
#'   2 = detailed diagnostics. Default: 0.
#' @param progress Logical. Show progress bars during long operations.
#'   Default: \code{FALSE}.
#' @param obstacle_strategy Character string controlling behavior when obstacle
#'   values are provided but the \pkg{gdistance} package is not installed.
#'   One of \code{"error"} (default; stop with an informative error),
#'   \code{"straight_line"} (fall back to straight-line corridors, ignoring
#'   obstacles), or \code{"disable_obstacles"} (silently drop obstacle data).
#' @param species_dispersal_distance Numeric. Typical movement distance for
#'   the focal species, in the same distance units as the analysis (pixels /
#'   meters / feet). Used by \code{"reachable_habitat_advanced"} and by
#'   habitat availability metrics. If \code{NULL} (default), the
#'   \code{max_search_distance} value is used as a proxy.
#' @param species_dispersal_kernel Character string. Dispersal probability
#'   kernel. Currently only \code{"exponential"} is supported (default).
#' @param min_patch_area_for_species Numeric. Minimum patch area (in analysis
#'   area units) for a patch to be included in species-level habitat
#'   availability calculations. Default: 0.
#' @param patch_area_scaling Character string controlling how patch area is
#'   transformed before weighting in habitat availability calculations.
#'   \code{"sqrt"} (default) applies square-root scaling, giving moderate
#'   weight to large patches. \code{"log"} applies logarithmic scaling
#'   (log(1 + area)), reducing the influence of very large patches.
#' @param mobility_detour_cap Numeric. Maximum detour ratio used by
#'   graph-based fluidity metrics. Controls how much longer an indirect
#'   route can be relative to the straight-line distance before it is
#'   considered non-functional. Default: 8.
#' @param redundancy_method Character string selecting the flow redundancy
#'   calculation method. \code{"ime"} (default) uses Inverse Mean
#'   Effective-resistance, measuring how many alternative routes exist.
#'   \code{"fri"} uses the Flow Redundancy Index, an alternative based on
#'   current-flow theory.
#' @param metric_weights Optional named numeric vector with elements
#'   \code{"mesh"}, \code{"lcc"}, \code{"pc"}, and \code{"flow"} controlling
#'   the composite connectivity score blend. Values are normalized to sum to
#'   1. If \code{NULL}, equal weights are used.
#' @param weight_m Optional numeric override for the mesh component weight
#'   in the composite score.
#' @param weight_lcc Optional numeric override for the LCC component weight.
#' @param weight_pc Optional numeric override for the PC component weight.
#' @param weight_f Optional numeric override for the flow component weight.
#' @param output_dir Optional character path. Directory for writing output
#'   files when \code{write_outputs = TRUE}.
#' @param output_prefix Optional character string prepended to output
#'   file names.
#' @param output_paths Optional named list of explicit output file paths,
#'   overriding the default naming convention.
#' @param write_outputs Logical. If \code{TRUE}, write output rasters and
#'   CSV files to disk. Default: \code{FALSE}.
#' @param keep_candidates Logical. If \code{TRUE}, include the full candidate
#'   corridor table in the result (useful for debugging). Default: \code{FALSE}.
#'
#' @return An object of class \code{"terralink_result"} (a list) with the
#'   following elements:
#'   \itemize{
#'     \item \code{corridors}: Data frame or sf object of selected corridor
#'       geometries with columns \code{patch1}, \code{patch2},
#'       \code{corridor_area}, \code{corridor_length}, \code{connected_area},
#'       and \code{network_area}.
#'     \item \code{patches}: SpatRaster of labeled patch cells (raster mode).
#'     \item \code{patch_table}: Data frame of patch attributes (id, area,
#'       centroid coordinates).
#'     \item \code{networks}: sf object of connected network polygons (one
#'       feature per component of patches + corridors).
#'     \item \code{corridor_raster}: SpatRaster where corridor cells are
#'       assigned values according to \code{corridor_cell_assignment}.
#'     \item \code{contiguous_raster}: SpatRaster labeling each contiguous
#'       patch-corridor network.
#'     \item \code{strategy}: The strategy key that was used.
#'     \item \code{summary}: Named list with run overview including
#'       \code{budget_total}, \code{budget_used}, \code{corridors_used},
#'       \code{candidate_edges}, \code{patches}, \code{strategy},
#'       \code{units}.
#'     \item \code{metrics}: Named list of PRE/POST landscape metrics. Each
#'       metric has a \code{_pre} (before corridors) and \code{_post} (after
#'       corridors) value. Key metrics: \code{total_connected_habitat_area},
#'       \code{largest_network_area}, \code{habitat_availability},
#'       \code{mean_effective_resistance} (lower is better),
#'       \code{mesh_norm}, \code{lcc}, \code{pc}, \code{flow_redundancy},
#'       \code{strategic_mobility}, \code{landscape_fluidity},
#'       \code{composite_connectivity}.
#'     \item \code{metrics_report}: Character vector with a human-readable
#'       PRE/POST metrics table. Print with
#'       \code{cat(result$metrics_report, sep = "\\n")}.
#'     \item \code{strategy_stats}: Named list of strategy-specific
#'       optimization statistics (e.g., primary vs. redundant links).
#'     \item \code{mode}: Character string \code{"raster"}.
#'     \item \code{inputs}: Named list echoing key input parameters.
#'     \item \code{run_stats}: Named list with \code{elapsed_s},
#'       \code{candidate_edges}, \code{candidate_pairs}.
#'     \item \code{warnings}: Character vector of any warnings raised during
#'       the run.
#'     \item \code{diagnostics}: List of diagnostic messages (e.g., why no
#'       corridors were selected).
#'   }
#'   The object has \code{print()}, \code{summary()}, and \code{plot()}
#'   methods.
#'
#' @section Parameter guidance:
#' \itemize{
#'   \item \strong{budget}: A practical starting point is often around 5--20
#'     percent of total habitat area. Run several budget levels and compare
#'     PRE/POST metrics.
#'   \item \strong{min_patch_size}: Use this to exclude patches too small to
#'     function as habitat in your planning context. For raster mode, 5--20
#'     pixels is a common starting range; for real landscapes, 1--10 ha can be
#'     a reasonable first pass.
#'   \item \strong{min_corridor_width}: Should reflect the minimum width for
#'     species movement. Depending on species and landscape context, 30--100 m
#'     is a common starting range for terrestrial mammals and 10--30 m for
#'     some birds.
#'   \item \strong{max_search_distance}: Should be at or above the maximum
#'     distance the focal species can cross non-habitat. 500--5000 m is a
#'     common starting range; increase if 0 corridors are generated.
#'   \item \strong{species_dispersal_distance}: Set to the focal species'
#'     typical natal or daily movement range. This directly affects the
#'     \code{"reachable_habitat_advanced"} strategy and habitat availability
#'     metrics.
#' }
#'
#' @examples
#' r <- terra::rast(
#'   nrows = 6, ncols = 6,
#'   xmin = 0, xmax = 600,
#'   ymin = 0, ymax = 600,
#'   crs = "EPSG:3857"
#' )
#' vals <- rep(0, terra::ncell(r))
#' vals[c(1, 2, 7, 8, 29, 30, 35, 36)] <- 1
#' terra::values(r) <- vals
#'
#' result <- terralink_raster(
#'   raster = r,
#'   patch_values = 1,
#'   budget = 15,
#'   min_patch_size = 2,
#'   min_corridor_width = 1,
#'   max_search_distance = 12,
#'   units = "pixels"
#' )
#' result$summary
#'
#' # Access PRE/POST metrics
#' result$metrics$largest_network_area_pre
#' result$metrics$largest_network_area_post
#'
#' # Print the full metrics report
#' cat(result$metrics_report, sep = "\n")
#' @export
terralink_raster <- function(
  raster,
  patch_values = NULL,
  patch_ranges = NULL,
  budget = NULL,
  budget_pixels = NULL,
  strategy = "most_connected_habitat",
  min_patch_size = 10,
  min_corridor_width = 3,
  corridor_cell_assignment = "sum_total_network_area",
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
  species_dispersal_distance = NULL,
  species_dispersal_kernel = HABITAT_AVAILABILITY_DEFAULT_KERNEL,
  min_patch_area_for_species = 0,
  patch_area_scaling = HABITAT_AVAILABILITY_DEFAULT_SCALING,
  mobility_detour_cap = 8,
  redundancy_method = "ime",
  metric_weights = NULL,
  weight_m = NULL,
  weight_lcc = NULL,
  weight_pc = NULL,
  weight_f = NULL,
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
    corridor_cell_assignment = corridor_cell_assignment,
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
    species_dispersal_distance = species_dispersal_distance,
    species_dispersal_kernel = species_dispersal_kernel,
    min_patch_area_for_species = min_patch_area_for_species,
    patch_area_scaling = patch_area_scaling,
    mobility_detour_cap = mobility_detour_cap,
    redundancy_method = redundancy_method,
    metric_weights = metric_weights,
    weight_m = weight_m,
    weight_lcc = weight_lcc,
    weight_pc = weight_pc,
    weight_f = weight_f,
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
#' Builds candidate corridors between polygon habitat patches and selects an
#' optimal corridor network under a budget constraint. This is the native
#' TerraLink workflow and is usually the better choice when planning inputs are
#' already polygon features.
#'
#' @param patches sf object with polygon geometry (one row per patch), or a
#'   file path to a GeoPackage / Shapefile. The CRS should be projected
#'   (e.g., UTM) so that area and distance calculations are meaningful.
#' @param budget Numeric. Total corridor area budget in hectares
#'   (\code{units = "metric"}) or acres (\code{units = "imperial"}).
#'   A reasonable starting point is 5--20 percent of your total patch area.
#' @param strategy Character string selecting the optimization objective.
#'   One of \code{"most_connected_habitat"} (default; maximizes total
#'   structurally connected habitat area), \code{"largest_single_network"}
#'   (maximizes the single largest connected component),
#'   \code{"landscape_fluidity"} (maximizes ease of movement and route
#'   redundancy), or \code{"reachable_habitat_advanced"} (maximizes
#'   dispersal-weighted reachable habitat; set
#'   \code{species_dispersal_distance}).
#' @param min_patch_size Numeric. Minimum patch area in hectares
#'   (\code{"metric"}) or acres (\code{"imperial"}). Patches smaller than
#'   this are dropped. Default: \code{NULL} (no filter).
#' @param min_corridor_width Numeric. Minimum corridor width in meters
#'   (\code{"metric"}) or feet (\code{"imperial"}). Controls the buffer
#'   applied to corridor center-lines. Typical values: 30--100 m for
#'   terrestrial mammals. Default: 100.
#' @param max_search_distance Numeric. Maximum edge-to-edge distance (meters
#'   or feet) between patches to consider a candidate corridor. Increase if
#'   few or no corridors are generated. Default: 5000.
#' @param obstacle_layers Optional sf object or file path to polygon barriers
#'   (roads, water bodies). Requires the \pkg{gdistance} package for
#'   shortest-path routing around obstacles.
#' @param obstacle_resolution Numeric. Raster cell size (in CRS units) used
#'   to rasterize obstacles for shortest-path routing. Smaller values give more
#'   accurate routing but increase computation time.
#' @param units Character string: \code{"metric"} (hectares / meters, default)
#'   or \code{"imperial"} (acres / feet).
#' @param max_pair_checks Integer. Upper limit on patch pairs evaluated.
#'   Default: 2,000,000.
#' @param max_candidates Integer. Upper limit on candidate corridors retained.
#'   Default: 200,000.
#' @param verbose Integer verbosity level: 0 = silent, 1 = progress, 2 =
#'   detailed. Default: 0.
#' @param progress Logical. Show progress bars. Default: \code{FALSE}.
#' @param obstacle_strategy Character string controlling behavior when
#'   obstacles are provided but \pkg{gdistance} is not installed. One of
#'   \code{"error"} (default; stop with an error),
#'   \code{"straight_line"} (fall back to straight-line corridors), or
#'   \code{"disable_obstacles"} (silently ignore obstacles).
#' @param return_crs Character string controlling the output CRS.
#'   \code{"input"} (default) returns outputs in the same CRS as the input
#'   patches. \code{"utm"} returns outputs in the UTM zone used internally.
#' @param species_dispersal_distance Numeric. Typical movement distance for
#'   the focal species in meters (\code{"metric"}) or feet
#'   (\code{"imperial"}). Used by \code{"reachable_habitat_advanced"} and
#'   habitat availability metrics. If \code{NULL} (default),
#'   \code{max_search_distance} is used as a proxy.
#' @param species_dispersal_kernel Character string. Dispersal probability
#'   kernel. Currently only \code{"exponential"} is supported (default).
#' @param min_patch_area_for_species Numeric. Minimum patch area (in analysis
#'   area units) for inclusion in species-level metrics. Default: 0.
#' @param patch_area_scaling Character string controlling how patch area is
#'   transformed before weighting. \code{"sqrt"} (default) applies square-root
#'   scaling, giving moderate weight to large patches. \code{"log"} applies
#'   logarithmic scaling, reducing the influence of very large patches.
#' @param patch_quality_field Optional character string naming a numeric
#'   column in \code{patches} that provides a 0--1 quality weight per patch
#'   (e.g., habitat suitability). Patches with higher quality contribute more
#'   to connectivity metrics.
#' @param mobility_detour_cap Numeric. Maximum detour ratio for fluidity
#'   metrics. Controls how much longer an indirect route can be before it is
#'   considered non-functional. Default: 8.
#' @param redundancy_method Character string selecting the flow redundancy
#'   method. \code{"ime"} (default) uses Inverse Mean Effective-resistance.
#'   \code{"fri"} uses the Flow Redundancy Index.
#' @param metric_weights Optional named numeric vector with elements
#'   \code{"mesh"}, \code{"lcc"}, \code{"pc"}, \code{"flow"} for the
#'   composite score. Normalized to sum to 1.
#' @param weight_m Optional numeric override for the mesh weight.
#' @param weight_lcc Optional numeric override for the LCC weight.
#' @param weight_pc Optional numeric override for the PC weight.
#' @param weight_f Optional numeric override for the flow weight.
#' @param output_dir Optional output directory for \code{write_outputs}.
#' @param output_prefix Optional name prefix for output files.
#' @param output_paths Optional named list of explicit output file paths.
#' @param write_outputs Logical. Write GeoPackage and CSV outputs to disk.
#'   Default: \code{FALSE}.
#' @param keep_candidates Logical. Include full candidate table in result.
#'   Default: \code{FALSE}.
#'
#' @return An object of class \code{"terralink_result"} (a list) with the
#'   following elements:
#'   \itemize{
#'     \item \code{corridors}: sf object of selected corridors with columns
#'       \code{patch1}, \code{patch2} (endpoint patch IDs),
#'       \code{corridor_area} (ha or ac), \code{corridor_length} (m or ft),
#'       \code{connected_area}, \code{network_area}, and geometry.
#'     \item \code{patches}: sf object of patches used in the analysis, with
#'       area and centroid attributes.
#'     \item \code{networks}: sf object of connected network polygons (one
#'       feature per component of patches + corridors).
#'     \item \code{summary}: Named list including \code{budget_total},
#'       \code{budget_used}, \code{corridors_used}, \code{candidate_edges},
#'       \code{patches}, \code{raw_patches}, \code{filtered_out},
#'       \code{primary_links}, \code{redundant_links}, \code{strategy},
#'       \code{units}.
#'     \item \code{metrics}: Named list of PRE/POST landscape connectivity
#'       metrics. Every metric has a \code{_pre} and \code{_post} value.
#'       Key metrics: \code{total_connected_habitat_area},
#'       \code{largest_network_area}, \code{habitat_availability},
#'       \code{mean_effective_resistance} (lower is better),
#'       \code{mesh_norm}, \code{lcc}, \code{pc}, \code{flow_redundancy},
#'       \code{strategic_mobility}, \code{landscape_fluidity},
#'       \code{composite_connectivity}.
#'     \item \code{metrics_report}: Character vector with a human-readable
#'       PRE/POST table. Print with
#'       \code{cat(result$metrics_report, sep = "\\n")}.
#'     \item \code{strategy_stats}: Named list of strategy-specific
#'       statistics.
#'     \item \code{mode}: Character string \code{"vector"}.
#'     \item \code{inputs}: Named list echoing key input parameters.
#'     \item \code{run_stats}: Named list with \code{elapsed_s},
#'       \code{candidate_edges}, \code{candidate_pairs}.
#'     \item \code{warnings}: Character vector of warnings.
#'     \item \code{diagnostics}: List of diagnostic messages.
#'   }
#'   The object has \code{print()}, \code{summary()}, and \code{plot()}
#'   methods.
#'
#' @section Parameter guidance:
#' \itemize{
#'   \item \strong{budget}: A practical starting point is often around 5--20
#'     percent of total patch area. Run multiple budgets and compare PRE/POST
#'     metrics to find the point of diminishing returns.
#'   \item \strong{min_corridor_width}: Depending on species and context,
#'     30--100 m can be a useful starting range for mammals and 10--30 m for
#'     some small birds.
#'   \item \strong{max_search_distance}: 500--5000 m is a common starting
#'     range. Increase if 0 corridors are generated.
#'   \item \strong{species_dispersal_distance}: Set to the focal species'
#'     typical natal or daily movement range. Directly affects
#'     \code{"reachable_habitat_advanced"} and habitat availability metrics.
#' }
#'
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
#'
#' # Access PRE/POST metrics
#' result$metrics$largest_network_area_pre
#' result$metrics$largest_network_area_post
#'
#' # Print the full metrics report
#' cat(result$metrics_report, sep = "\n")
#' @export
terralink_vector <- function(
  patches,
  budget,
  strategy = "most_connected_habitat",
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
  species_dispersal_distance = NULL,
  species_dispersal_kernel = HABITAT_AVAILABILITY_DEFAULT_KERNEL,
  min_patch_area_for_species = 0,
  patch_area_scaling = HABITAT_AVAILABILITY_DEFAULT_SCALING,
  patch_quality_field = NULL,
  mobility_detour_cap = 8,
  redundancy_method = "ime",
  metric_weights = NULL,
  weight_m = NULL,
  weight_lcc = NULL,
  weight_pc = NULL,
  weight_f = NULL,
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
    species_dispersal_distance = species_dispersal_distance,
    species_dispersal_kernel = species_dispersal_kernel,
    min_patch_area_for_species = min_patch_area_for_species,
    patch_area_scaling = patch_area_scaling,
    patch_quality_field = patch_quality_field,
    mobility_detour_cap = mobility_detour_cap,
    redundancy_method = redundancy_method,
    metric_weights = metric_weights,
    weight_m = weight_m,
    weight_lcc = weight_lcc,
    weight_pc = weight_pc,
    weight_f = weight_f,
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
