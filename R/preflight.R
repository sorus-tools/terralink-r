# Preflight validation helpers for TerraLink

terralink_make_valid <- function(x) {
  if (is.function(sf::st_make_valid)) {
    return(sf::st_make_valid(x))
  }
  if (requireNamespace("lwgeom", quietly = TRUE)) {
    fn <- tryCatch(get("st_make_valid", envir = asNamespace("lwgeom"), inherits = FALSE), error = function(e) NULL)
    if (is.function(fn)) {
      return(fn(x))
    }
  }
  NULL
}

terralink_preflight_raster <- function(
  raster,
  allow_large = FALSE,
  ctx = NULL,
  size_threshold = 40000000,
  fine_threshold = 10.0,
  fine_critical_count = 2000000
) {
  if (!inherits(raster, "SpatRaster")) {
    terralink_abort("raster must be a terra::SpatRaster.", class = "terralink_error_input")
  }
  n_cells <- tryCatch(terra::ncell(raster), error = function(e) NA_real_)
  if (!is.finite(n_cells) || n_cells <= 0) {
    terralink_abort("Raster has no cells or cannot be read.", class = "terralink_error_input")
  }
  res_xy <- tryCatch(terra::res(raster), error = function(e) c(NA_real_, NA_real_))
  pixel_size <- suppressWarnings(max(res_xy, na.rm = TRUE))
  pixel_area <- suppressWarnings(abs(res_xy[1] * res_xy[2]))
  if (!is.finite(pixel_area) || pixel_area <= 0) {
    pixel_area <- NA_real_
  }

  details <- c(
    sprintf("Raster cells: %s", format(n_cells, big.mark = ",")),
    sprintf("Pixel size: %s", if (is.finite(pixel_size)) terralink_format_number(pixel_size, 3) else "unknown"),
    sprintf("Approx. memory (1 layer): %s", terralink_format_bytes(n_cells * 8))
  )
  fixes <- c(
    "Crop to a smaller extent",
    "Aggregate/resample to larger pixels",
    "Clip to study area tiles",
    "Increase pixel size in your input"
  )

  if (n_cells >= size_threshold && !isTRUE(allow_large)) {
    terralink_abort(
      "Raster is too large for reliable processing.",
      class = "terralink_error_compute_limit",
      details = details,
      fix = fixes
    )
  }
  if (is.finite(pixel_size) && pixel_size > 0 && pixel_size < fine_threshold) {
    msg <- sprintf("High-resolution data detected (~%s map units per pixel).", terralink_format_number(pixel_size, 3))
    if (n_cells >= fine_critical_count && !isTRUE(allow_large)) {
      terralink_abort(
        msg,
        class = "terralink_error_compute_limit",
        details = details,
        fix = c(
          "Resample to a coarser resolution",
          "Reduce extent or increase pixel size",
          "Set allow_large = TRUE to override"
        )
      )
    } else {
      terralink_warn(paste(msg, "Proceeding because raster size is modest."), ctx = ctx)
    }
  }

  if (n_cells >= size_threshold && isTRUE(allow_large)) {
    terralink_warn("Large raster detected; run may be slow or memory-intensive.", ctx = ctx)
  }

  list(
    n_cells = n_cells,
    res = res_xy,
    pixel_size = pixel_size,
    pixel_area = pixel_area
  )
}

terralink_raster_value_summary <- function(raster, max_values = 10, sample_size = 100000) {
  n_cells <- tryCatch(terra::ncell(raster), error = function(e) 0)
  if (!is.finite(n_cells) || n_cells <= 0) return(NULL)

  values <- NULL
  if (n_cells <= sample_size) {
    values <- tryCatch(terra::values(raster, mat = FALSE), error = function(e) NULL)
  } else {
    values <- tryCatch(
      terra::spatSample(raster, size = sample_size, method = "random", na.rm = FALSE, as.raster = FALSE)[[1]],
      error = function(e) NULL
    )
  }
  if (is.null(values) || length(values) == 0) return(NULL)

  na_count <- sum(is.na(values))
  values <- values[!is.na(values)]
  if (length(values) == 0) return(NULL)

  tbl <- sort(table(values), decreasing = TRUE)
  top <- utils::head(tbl, max_values)
  total <- sum(tbl)
  items <- vapply(names(top), function(name) {
    cnt <- as.integer(top[[name]])
    pct <- (cnt / total) * 100
    sprintf("%s: %s (%.1f%%)", name, format(cnt, big.mark = ","), pct)
  }, character(1))

  list(
    sampled = min(n_cells, sample_size),
    total_sample = length(values),
    na_count = na_count,
    items = items
  )
}

terralink_preflight_vector <- function(patches_sf, ctx = NULL) {
  if (!inherits(patches_sf, "sf") && !inherits(patches_sf, "sfc")) {
    terralink_abort("patches must be an sf object.", class = "terralink_error_input")
  }
  patches_sf <- sf::st_as_sf(patches_sf)
  if (is.na(sf::st_crs(patches_sf))) {
    terralink_abort(
      "Vector input has no CRS defined.",
      class = "terralink_error_input",
      fix = "Assign a CRS with sf::st_set_crs() before running TerraLink."
    )
  }
  geom_types <- as.character(sf::st_geometry_type(patches_sf))
  if (!all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    terralink_abort("patches must be polygon features.", class = "terralink_error_geometry")
  }

  empty_mask <- sf::st_is_empty(patches_sf)
  if (any(empty_mask)) {
    terralink_warn(sprintf("Removed %s empty geometry feature(s).", sum(empty_mask)), ctx = ctx)
    patches_sf <- patches_sf[!empty_mask, , drop = FALSE]
  }

  if (nrow(patches_sf) == 0) {
    terralink_abort("No polygon features found after removing empty geometries.", class = "terralink_error_geometry")
  }

  valid_mask <- sf::st_is_valid(patches_sf)
  if (any(!valid_mask)) {
    valid_fixed <- terralink_make_valid(patches_sf)
    if (!is.null(valid_fixed)) {
      terralink_warn(sprintf("Fixing %s invalid geometry feature(s).", sum(!valid_mask)), ctx = ctx)
      patches_sf <- valid_fixed
      valid_mask <- sf::st_is_valid(patches_sf)
      if (any(!valid_mask)) {
        terralink_abort(
          "Invalid polygon geometries remain after st_make_valid().",
          class = "terralink_error_geometry",
          fix = c("Inspect invalid geometries with sf::st_is_valid()", "Simplify or repair geometries before running")
        )
      }
    } else {
      terralink_abort(
        "Invalid polygon geometries detected.",
        class = "terralink_error_geometry",
        details = sprintf("Invalid features: %s", paste(which(!valid_mask), collapse = ", ")),
        fix = c("Install sf >= 1.0 or lwgeom for st_make_valid()", "Repair geometries before running")
      )
    }
  }

  patches_sf
}

terralink_check_candidate_limits <- function(
  n_patches,
  max_pair_checks = 2000000,
  ctx = NULL,
  scope = "candidate pairs"
) {
  if (!is.finite(n_patches) || n_patches < 2) return(invisible(TRUE))
  possible_pairs <- (n_patches * (n_patches - 1)) / 2
  if (!is.null(max_pair_checks) && is.finite(max_pair_checks) && possible_pairs > max_pair_checks) {
    terralink_abort(
      sprintf("%s estimated to exceed safe limits.", scope),
      class = "terralink_error_compute_limit",
      details = c(
        sprintf("Patches: %s", format(n_patches, big.mark = ",")),
        sprintf("Possible pairs: %s", format(round(possible_pairs), big.mark = ",")),
        sprintf("Limit: %s", format(max_pair_checks, big.mark = ","))
      ),
      fix = c(
        "Increase min_patch_size to reduce patch count",
        "Reduce max_search_distance",
        "Clip to a smaller extent",
        "Increase max_pair_checks to override"
      )
    )
  }
  invisible(TRUE)
}

terralink_check_candidate_count <- function(count, max_candidates = 200000, ctx = NULL, scope = "Candidate corridors", override_param = "max_candidates") {
  if (!is.null(max_candidates) && is.finite(max_candidates) && is.finite(count) && count > max_candidates) {
    terralink_abort(
      sprintf("%s exceed safe limits.", scope),
      class = "terralink_error_compute_limit",
      details = c(
        sprintf("%s: %s", scope, format(count, big.mark = ",")),
        sprintf("Limit: %s", format(max_candidates, big.mark = ","))
      ),
      fix = c(
        "Reduce max_search_distance",
        "Increase min_patch_size",
        sprintf("Increase %s to override", override_param)
      )
    )
  }
  invisible(TRUE)
}
