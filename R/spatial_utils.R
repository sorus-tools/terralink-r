# Spatial helper utilities for TerraLink

terralink_check_package <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- if (is.null(purpose)) {
      sprintf("Package '%s' is required.", pkg)
    } else {
      sprintf("Package '%s' is required for %s.", pkg, purpose)
    }
    terralink_abort(
      msg,
      class = "terralink_error_dependency",
      fix = sprintf("Install with install.packages('%s')", pkg)
    )
  }
  TRUE
}

terralink_resolve_raster <- function(raster) {
  if (inherits(raster, "SpatRaster")) {
    return(raster)
  }
  if (is.character(raster) && length(raster) == 1) {
    if (!file.exists(raster)) {
      terralink_abort(
        sprintf("Raster path does not exist: %s", raster),
        class = "terralink_error_input"
      )
    }
    return(terra::rast(raster))
  }
  terralink_abort("raster must be a terra::SpatRaster or file path.", class = "terralink_error_input")
}

terralink_resolve_vector <- function(patches, quiet = TRUE) {
  if (inherits(patches, "sf") || inherits(patches, "sfc")) {
    return(sf::st_as_sf(patches))
  }
  if (is.character(patches) && length(patches) == 1) {
    if (!file.exists(patches)) {
      terralink_abort(
        sprintf("Vector path does not exist: %s", patches),
        class = "terralink_error_input"
      )
    }
    return(sf::st_read(patches, quiet = quiet))
  }
  terralink_abort("patches must be an sf object or file path.", class = "terralink_error_input")
}

terralink_is_projected <- function(crs) {
  if (is.null(crs)) return(FALSE)
  crs_obj <- NULL

  if (inherits(crs, "SpatRaster") || inherits(crs, "SpatVector")) {
    wkt <- tryCatch(terra::crs(crs), error = function(e) "")
    if (is.character(wkt) && length(wkt) == 1 && nzchar(wkt)) {
      crs_obj <- tryCatch(sf::st_crs(wkt), error = function(e) NULL)
    }
  } else if (inherits(crs, "crs")) {
    crs_obj <- crs
  } else if (is.character(crs) && length(crs) == 1 && nzchar(crs)) {
    crs_obj <- tryCatch(sf::st_crs(crs), error = function(e) NULL)
  } else if (is.list(crs) && !is.null(crs$wkt)) {
    crs_obj <- tryCatch(sf::st_crs(crs$wkt), error = function(e) NULL)
  } else if (is.data.frame(crs) && "code" %in% names(crs) && nrow(crs) >= 1) {
    code <- suppressWarnings(as.integer(crs$code[[1]]))
    if (is.finite(code)) {
      crs_obj <- tryCatch(sf::st_crs(code), error = function(e) NULL)
    }
  }

  if (is.null(crs_obj) || is.na(crs_obj)) return(FALSE)
  out <- tryCatch(!sf::st_is_longlat(crs_obj), error = function(e) FALSE)
  isTRUE(out)
}

terralink_pick_utm_crs <- function(geom) {
  sf_obj <- sf::st_as_sf(geom)
  sf_obj <- sf::st_transform(sf_obj, 4326)
  coords <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf_obj)))
  if (nrow(coords) == 0) {
    terralink_abort("Unable to compute centroid for UTM selection.", class = "terralink_error_geometry")
  }
  lon <- coords[1, 1]
  lat <- coords[1, 2]
  zone <- floor((lon + 180) / 6) + 1
  epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
  sf::st_crs(epsg)
}

terralink_unit_labels <- function(units) {
  units <- match.arg(units, c("pixels", "metric", "imperial"))
  if (units == "pixels") {
    list(distance = "px", area = "px^2")
  } else if (units == "imperial") {
    list(distance = "ft", area = "ac")
  } else {
    list(distance = "m", area = "ha")
  }
}

terralink_convert_raster_units <- function(
  raster,
  units,
  budget,
  min_patch_size,
  min_corridor_width,
  max_search_distance
) {
  units <- match.arg(units, c("pixels", "metric", "imperial"))
  res_xy <- terra::res(raster)
  pixel_size <- mean(res_xy)
  pixel_area <- abs(res_xy[1] * res_xy[2])
  if (units == "pixels") {
    return(list(
      budget_pixels = as.numeric(budget),
      min_patch_size_px = as.numeric(min_patch_size),
      min_corridor_width_px = as.numeric(min_corridor_width),
      max_search_distance_px = as.numeric(max_search_distance),
      pixel_size = pixel_size,
      pixel_area = pixel_area
    ))
  }
  if (!terralink_is_projected(raster)) {
    terralink_abort(
      "Raster must be in a projected CRS to use metric/imperial units.",
      class = "terralink_error_scale",
      fix = c("Project your raster to a CRS in meters/feet", "Use units = 'pixels'")
    )
  }

  if (units == "imperial") {
    dist_factor <- 0.3048  # feet to meters
    area_factor <- 4046.8564224  # acres to square meters
  } else {
    dist_factor <- 1.0
    area_factor <- 10000.0  # hectares to square meters
  }

  budget_m2 <- as.numeric(budget) * area_factor
  min_patch_m2 <- as.numeric(min_patch_size) * area_factor
  min_corridor_m <- as.numeric(min_corridor_width) * dist_factor
  max_search_m <- as.numeric(max_search_distance) * dist_factor

  list(
    budget_pixels = budget_m2 / pixel_area,
    min_patch_size_px = min_patch_m2 / pixel_area,
    min_corridor_width_px = min_corridor_m / pixel_size,
    max_search_distance_px = max_search_m / pixel_size,
    pixel_size = pixel_size,
    pixel_area = pixel_area
  )
}

terralink_mask_from_values_ranges <- function(raster, values = NULL, ranges = NULL) {
  if ((is.null(values) || length(values) == 0) && (is.null(ranges) || length(ranges) == 0)) {
    return(NULL)
  }
  mask <- NULL
  if (!is.null(values) && length(values) > 0) {
    vals <- as.numeric(values)
    mask <- terra::app(raster, fun = function(x) x %in% vals)
  }
  if (!is.null(ranges) && length(ranges) > 0) {
    for (r in ranges) {
      if (length(r) < 2) next
      low <- min(r, na.rm = TRUE)
      high <- max(r, na.rm = TRUE)
      range_mask <- terra::app(raster, fun = function(x) x >= low & x <= high)
      if (is.null(mask)) {
        mask <- range_mask
      } else {
        mask <- mask | range_mask
      }
    }
  }
  if (is.null(mask)) return(NULL)
  mask[is.na(mask)] <- FALSE
  mask
}

terralink_erode_mask <- function(mask, width_px) {
  width_px <- as.integer(round(width_px))
  if (is.null(width_px) || width_px <= 1) {
    return(mask)
  }
  w <- matrix(1, nrow = width_px, ncol = width_px)
  terra::focal(mask, w = w, fun = min, na.policy = "omit", fillvalue = 0)
}

terralink_safe_name <- function(name, max_len = 64) {
  safe <- gsub("[^A-Za-z0-9_-]+", "_", name)
  safe <- gsub("^_+|_+$", "", safe)
  if (!nzchar(safe)) safe <- "layer"
  substr(safe, 1, max_len)
}

terralink_format_number <- function(value, digits = 2) {
  if (is.null(value) || is.na(value)) return("")
  sprintf(paste0("%.", digits, "f"), as.numeric(value))
}

terralink_default_output_dir <- function(input) {
  if (is.character(input) && length(input) == 1 && nzchar(input)) {
    return(dirname(normalizePath(input, mustWork = FALSE)))
  }
  if (inherits(input, "SpatRaster")) {
    srcs <- tryCatch(terra::sources(input), error = function(e) NULL)
    if (!is.null(srcs) && length(srcs) > 0 && nzchar(srcs[[1]])) {
      return(dirname(normalizePath(srcs[[1]], mustWork = FALSE)))
    }
  }
  getwd()
}
