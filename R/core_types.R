#' Create a patch descriptor
#'
#' @param id Patch identifier.
#' @param weight Numeric patch weight (area, quality).
#' @param geometry Optional geometry object.
#' @return A patch list with class 'terralink_patch'.
#' @export
new_patch <- function(id, weight, geometry = NULL) {
  structure(
    list(id = id, weight = weight, geometry = geometry),
    class = "terralink_patch"
  )
}

#' Create a candidate corridor descriptor
#'
#' @param patch_ids Integer vector of patch ids.
#' @param cost Numeric cost of the corridor.
#' @param weight Numeric benefit or ROI.
#' @param geometry Optional geometry object.
#' @return A candidate list with class 'terralink_candidate'.
#' @export
new_candidate <- function(patch_ids, cost, weight = NULL, geometry = NULL) {
  structure(
    list(patch_ids = patch_ids, cost = cost, weight = weight, geometry = geometry),
    class = "terralink_candidate"
  )
}
