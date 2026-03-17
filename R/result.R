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
  mode <- x$mode %||% "unknown"
  strategy <- x$summary$strategy %||% x$strategy %||% ""
  units_label <- x$summary$units %||% x$inputs$units %||% ""
  area_unit <- switch(units_label, metric = "ha", imperial = "ac", pixels = "px", "")
  dist_unit <- switch(units_label, metric = "m", imperial = "ft", pixels = "px", "")

  header <- paste0("TerraLink result (", mode, ")")
  if (nzchar(strategy)) header <- paste0(header, " -- strategy: ", strategy)
  cat(header, "\n", sep = "")

  if (!is.null(x$summary)) {
    cat("  Patches:    ", x$summary$patches %||% "?", "\n", sep = "")
    cat("  Corridors:  ", x$summary$corridors_used %||% "?", "\n", sep = "")
    cat("  Candidates: ", x$summary$candidate_edges %||% "?", "\n", sep = "")
    if (!is.null(x$summary$budget_total)) {
      cat("  Budget:     ", formatC(x$summary$budget_used %||% 0, format = "f", digits = 2),
          " / ", formatC(x$summary$budget_total, format = "f", digits = 2),
          if (nzchar(area_unit)) paste0(" ", area_unit) else "", "\n", sep = "")
    } else if (!is.null(x$summary$budget_used)) {
      cat("  Budget used:", formatC(x$summary$budget_used, format = "f", digits = 2),
          if (nzchar(area_unit)) paste0(" ", area_unit) else "", "\n", sep = "")
    }
    if (!is.null(x$summary$primary_links) || !is.null(x$summary$redundant_links)) {
      cat("  Links:      ", x$summary$primary_links %||% 0, " primary, ",
          x$summary$redundant_links %||% 0, " redundant\n", sep = "")
    }
  }

  # Show key metric deltas if available

  if (!is.null(x$metrics)) {
    m <- x$metrics
    cat("\n  Key metrics (PRE -> POST):\n")
    fmt <- function(val, digits = 2) {
      num <- suppressWarnings(as.numeric(val))
      if (is.null(val) || !is.finite(num)) return("--")
      abs_num <- abs(num)
      if (abs_num > 0 && abs_num < 10^(-digits)) {
        return(formatC(num, format = "e", digits = max(2, digits)))
      }
      formatC(num, format = "f", digits = digits)
    }
    show_metric <- function(label, pre_name, post_name, digits = 2, lower_better = FALSE) {
      pre <- m[[pre_name]]
      post <- m[[post_name]]
      if (is.null(pre) && is.null(post)) return(invisible(NULL))
      delta <- as.numeric(post %||% 0) - as.numeric(pre %||% 0)
      arrow <- if (!is.finite(delta) || abs(delta) < 1e-12) "="
               else if ((!lower_better && delta > 0) || (lower_better && delta < 0)) "+"
               else "-"
      cat(sprintf("    %-28s %s -> %s  [%s]\n", label,
                  fmt(pre, digits), fmt(post, digits), arrow))
    }
    show_metric("Connected habitat area", "total_connected_habitat_area_pre",
                "total_connected_habitat_area_post")
    show_metric("Largest network area", "largest_network_area_pre",
                "largest_network_area_post")
    show_metric("Habitat availability", "habitat_availability_pre",
                "habitat_availability_post", digits = 4)
    show_metric("Landscape fluidity", "landscape_fluidity_pre",
                "landscape_fluidity_post", digits = 4)
    show_metric("Composite connectivity", "composite_connectivity_pre",
                "composite_connectivity_post", digits = 4)
    show_metric("Mean eff. resistance", "mean_effective_resistance_pre",
                "mean_effective_resistance_post", digits = 4, lower_better = TRUE)
    cat("\n  Use cat(result$metrics_report, sep=\"\\n\") for the full metrics table.\n")
  }

  warn_count <- length(x$warnings %||% character(0))
  if (warn_count > 0) {
    cat("  Warnings:   ", warn_count, "\n", sep = "")
  }
  diag_msg <- x$diagnostics$message %||% NULL
  if (!is.null(diag_msg) && nzchar(diag_msg)) {
    cat("  Diagnostic: ", diag_msg, "\n", sep = "")
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
    metrics = object$metrics,
    strategy_stats = object$strategy_stats,
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
