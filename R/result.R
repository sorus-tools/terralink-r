# Result object helpers for TerraLink

terralink_as_result <- function(result, mode, inputs = list(), run_stats = list(), warnings = NULL, diagnostics = NULL) {
  result$mode <- mode
  result$inputs <- inputs
  result$run_stats <- run_stats
  result$warnings <- warnings %||% character(0)
  result$diagnostics <- diagnostics %||% list()
  class(result) <- unique(c("terralink_result", class(result)))
  result
}

#' @export
print.terralink_result <- function(x, ...) {
  cat("TerraLink result (", x$mode %||% "unknown", ")\n", sep = "")
  if (!is.null(x$summary)) {
    cat("  Corridors: ", x$summary$corridors_used %||% "?", "\n", sep = "")
    cat("  Patches:   ", x$summary$patches %||% "?", "\n", sep = "")
    cat("  Candidates:", x$summary$candidate_edges %||% "?", "\n", sep = "")
    if (!is.null(x$summary$budget_used)) {
      cat("  Budget used:", x$summary$budget_used, "\n")
    }
  }
  warn_count <- length(x$warnings %||% character(0))
  if (warn_count > 0) {
    cat("  Warnings:  ", warn_count, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.terralink_result <- function(object, ...) {
  list(
    mode = object$mode %||% "unknown",
    summary = object$summary,
    inputs = object$inputs,
    run_stats = object$run_stats,
    warnings = object$warnings,
    diagnostics = object$diagnostics,
    metrics_report = object$metrics_report
  )
}

#' @export
plot.terralink_result <- function(x, ...) {
  if (!is.null(x$corridor_raster) && requireNamespace("terra", quietly = TRUE)) {
    terra::plot(x$corridor_raster, main = "TerraLink corridors", ...)
    return(invisible(x))
  }
  if (!is.null(x$contiguous_raster) && requireNamespace("terra", quietly = TRUE)) {
    terra::plot(x$contiguous_raster, main = "TerraLink contiguous network", ...)
    return(invisible(x))
  }
  if (!is.null(x$corridors) && inherits(x$corridors, "sf")) {
    plot(sf::st_geometry(x$corridors), col = "#0072B2", border = NA, ...)
    if (!is.null(x$patches) && inherits(x$patches, "sf")) {
      plot(sf::st_geometry(x$patches), col = "#56B4E9", border = "#2C7FB8", add = TRUE)
    }
    return(invisible(x))
  }
  if (!is.null(x$patches) && inherits(x$patches, "SpatRaster") && requireNamespace("terra", quietly = TRUE)) {
    terra::plot(x$patches, main = "TerraLink patches", ...)
    return(invisible(x))
  }
  warning("No plottable outputs found in result.", call. = FALSE)
  invisible(x)
}
