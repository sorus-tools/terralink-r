#' terralink: Connectivity Corridor Optimization
#'
#' Standalone R implementation of habitat connectivity corridor optimization
#' algorithms for raster and vector workflows, without desktop GIS
#' dependencies.
#'
#' Main user entry points:
#'
#' - `terralink_raster()` for raster inputs
#' - `terralink_vector()` for polygon patch inputs
#' - `terralink_run()` for a single generic wrapper
#'
#' Raster inputs are funneled through TerraLink's vector corridor pipeline so
#' the R package matches the current QGIS plugin workflow.
#'
#' @seealso \itemize{
#'   \item \code{\link{terralink_raster}}
#'   \item \code{\link{terralink_vector}}
#'   \item \code{\link{terralink_run}}
#'   \item \url{https://github.com/sorus-tools/terralink-r}
#'   \item Report bugs at \url{https://github.com/sorus-tools/terralink-r/issues}
#' }
#' @import R6
"_PACKAGE"
