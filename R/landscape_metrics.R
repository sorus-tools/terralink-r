# Landscape metrics for TerraLink (raster + vector)

terralink_metrics_from_mask <- function(mask, raster, units = "pixels") {
  units <- match.arg(units, c("pixels", "metric", "imperial"))
  if (!inherits(mask, "SpatRaster")) {
    terralink_abort("mask must be a terra::SpatRaster.", class = "terralink_error_input")
  }
  if (!inherits(raster, "SpatRaster")) {
    terralink_abort("raster must be a terra::SpatRaster.", class = "terralink_error_input")
  }

  res_xy <- terra::res(raster)
  pixel_area <- abs(res_xy[1] * res_xy[2])
  pixel_size <- mean(res_xy)

  if (units == "metric") {
    area_div <- 10000.0  # ha
    dist_mult <- 1.0
    area_unit <- "ha"
    dist_unit <- "m"
  } else if (units == "imperial") {
    area_div <- 4046.8564224  # ac from m^2
    dist_mult <- 3.28084      # m to ft
    area_unit <- "ac"
    dist_unit <- "ft"
  } else {
    area_div <- 1.0
    dist_mult <- 1.0
    area_unit <- "px^2"
    dist_unit <- "px"
  }

  mask_bin <- terra::ifel(!is.na(mask) & mask > 0, 1, 0)
  mask_vals <- terra::values(mask_bin, mat = FALSE)
  hab_pixels <- sum(mask_vals, na.rm = TRUE)
  if (is.null(mask_vals) || hab_pixels <= 0) {
    return(list(
      metrics = list(
        total_area = 0,
        num_patches = 0,
        mesh = 0,
        lps = 0,
        split = 0,
        pd = 0,
        lpi = 0,
        total_edge = 0,
        msi = 0,
        ed = 0,
        mean_para = 0,
        frac = 0,
        ai = 0,
        cai = 0,
        mean_gyrate = 0
      ),
      area_unit = area_unit,
      dist_unit = dist_unit,
      labels = terra::rast(raster)
    ))
  }

  labels <- terra::patches(terra::ifel(mask_bin > 0, 1, NA), directions = 8)
  freq <- terra::freq(labels)
  freq <- freq[!is.na(freq$value) & freq$value != 0, , drop = FALSE]
  if (nrow(freq) == 0) {
    return(list(
      metrics = list(
        total_area = 0,
        num_patches = 0,
        mesh = 0,
        lps = 0,
        split = 0,
        pd = 0,
        lpi = 0,
        total_edge = 0,
        msi = 0,
        ed = 0,
        mean_para = 0,
        frac = 0,
        ai = 0,
        cai = 0,
        mean_gyrate = 0
      ),
      area_unit = area_unit,
      dist_unit = dist_unit,
      labels = labels
    ))
  }

  patch_ids <- freq$value
  patch_counts <- freq$count
  patch_area_m2 <- patch_counts * pixel_area
  patch_area_unit <- patch_area_m2 / area_div
  total_area_unit <- sum(patch_area_unit)
  num_patches <- length(patch_ids)

  # Edge map using 4-neighbor kernel
  kernel <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, byrow = TRUE)
  neighbor_count <- terra::focal(mask_bin, w = kernel, fun = sum, na.policy = "omit", fillvalue = 0)
  edge_map <- (4 - neighbor_count) * mask_bin
  edge_vals <- terra::values(edge_map, mat = FALSE)
  edge_vals[edge_vals < 0] <- 0
  total_edge_m <- sum(edge_vals, na.rm = TRUE) * pixel_size
  total_edge_unit <- total_edge_m * dist_mult

  # Patch perimeter per patch (approx using edge map)
  label_vals <- terra::values(labels, mat = FALSE)
  per_patch_edge <- tapply(edge_vals, label_vals, sum, na.rm = TRUE)
  per_patch_edge <- per_patch_edge[as.character(patch_ids)]
  per_patch_edge[is.na(per_patch_edge)] <- 0
  patch_perim_base <- per_patch_edge * pixel_size
  patch_area_base <- patch_counts * pixel_area

  mesh <- if (total_area_unit > 0) sum(patch_area_unit^2) / total_area_unit else 0
  split <- if (sum(patch_area_unit^2) > 0) (total_area_unit^2) / sum(patch_area_unit^2) else 0
  lps <- if (num_patches > 0) max(patch_area_unit) else 0
  lpi <- if (total_area_unit > 0) (lps / total_area_unit) * 100 else 0
  pd <- if (total_area_unit > 0) (num_patches / total_area_unit) * 100 else 0
  ed <- if (total_area_unit > 0) total_edge_unit / total_area_unit else 0

  withCallingHandlers({
    msi <- if (num_patches > 0) mean(0.25 * patch_perim_base / sqrt(patch_area_base), na.rm = TRUE) else 0
  }, warning = function(w) {})

  mean_para <- if (num_patches > 0) mean(patch_perim_base / patch_area_base, na.rm = TRUE) else 0
  frac <- 0
  if (num_patches > 0) {
    frac_vals <- (2 * log(0.25 * patch_perim_base)) / log(patch_area_base)
    frac <- mean(frac_vals[is.finite(frac_vals)], na.rm = TRUE)
    if (is.nan(frac)) frac <- 0
  }

  # Aggregation index approximation
  denom_ai <- (2 * hab_pixels - 2 * sqrt(hab_pixels) - 2)
  ai <- if (denom_ai > 0) ((hab_pixels * 4 - sum(edge_vals, na.rm = TRUE)) / 2) / denom_ai * 100 else 0

  # Core area index (1-pixel erosion)
  core_kernel <- matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, byrow = TRUE)
  core_mask <- terra::focal(mask_bin, w = core_kernel, fun = min, na.policy = "omit", fillvalue = 0)
  core_pixels <- sum(terra::values(core_mask, mat = FALSE), na.rm = TRUE)
  cai <- if (hab_pixels > 0) (core_pixels / hab_pixels) * 100 else 0

  # Mean gyration (limit to 200 patches for speed)
  mean_gyrate <- 0
  if (num_patches > 0) {
    limit_ids <- patch_ids
    if (num_patches > 1000) {
      limit_ids <- patch_ids[seq_len(min(200, length(patch_ids)))]
    }
    pts <- terra::as.data.frame(labels, xy = TRUE, na.rm = TRUE)
    if (nrow(pts) > 0) {
      value_col <- setdiff(names(pts), c("x", "y"))[1]
      pts <- pts[pts[[value_col]] %in% limit_ids, , drop = FALSE]
      if (nrow(pts) > 0) {
        coords_by <- split(pts[, c("x", "y")], pts[[value_col]])
        gyr <- vapply(coords_by, function(df) {
          if (nrow(df) < 2) return(0)
          cx <- mean(df$x)
          cy <- mean(df$y)
          dist <- sqrt((df$x - cx)^2 + (df$y - cy)^2)
          mean(dist) * dist_mult
        }, numeric(1))
        mean_gyrate <- mean(gyr, na.rm = TRUE)
        if (is.nan(mean_gyrate)) mean_gyrate <- 0
      }
    }
  }

  list(
    metrics = list(
      total_area = total_area_unit,
      num_patches = num_patches,
      mesh = mesh,
      lps = lps,
      split = split,
      pd = pd,
      lpi = lpi,
      total_edge = total_edge_unit,
      msi = msi,
      ed = ed,
      mean_para = mean_para,
      frac = frac,
      ai = ai,
      cai = cai,
      mean_gyrate = mean_gyrate
    ),
    area_unit = area_unit,
    dist_unit = dist_unit,
    labels = labels
  )
}

terralink_landscape_report <- function(pre_mask, post_mask, raster, units = "pixels", label = "TerraLink") {
  pre <- terralink_metrics_from_mask(pre_mask, raster, units = units)
  post <- terralink_metrics_from_mask(post_mask, raster, units = units)
  connected_groups <- 0L
  if (inherits(pre$labels, "SpatRaster") && inherits(post$labels, "SpatRaster")) {
    pre_vals <- terra::values(pre$labels, mat = FALSE)
    post_vals <- terra::values(post$labels, mat = FALSE)
    post_ids <- sort(unique(post_vals[is.finite(post_vals) & post_vals > 0]))
    if (length(post_ids) > 0) {
      for (gid in post_ids) {
        overlap <- unique(pre_vals[post_vals == gid])
        overlap <- overlap[is.finite(overlap) & overlap > 0]
        if (length(overlap) >= 2) connected_groups <- connected_groups + 1L
      }
    }
  }

  header <- paste(rep("=", 100), collapse = "")
  results <- c(
    header,
    paste(" LANDSCAPE ANALYSIS:", label),
    header,
    sprintf("%-30s | %-15s | %-15s | %s", "METRIC NAME", "PRE", "POST", "INTERPRETATION"),
    paste(rep("-", 100), collapse = "")
  )

  fmt <- function(x) terralink_format_number(x, digits = 2)
  results <- c(
    results,
    sprintf("%-30s | %-15s | %-15s | %s", "Num. Patches (NP)", format(pre$metrics$num_patches, big.mark = ","), format(post$metrics$num_patches, big.mark = ","), "count"),
    sprintf("%-30s | %-15s | %-15s | %s", "Num. Connected Groups", "0", format(connected_groups, big.mark = ","), "Post = networks"),
    sprintf("%-30s | %-15s | %-15s | %s", "Total Area", fmt(pre$metrics$total_area), fmt(post$metrics$total_area), pre$area_unit),
    sprintf("%-30s | %-15s | %-15s | %s", "Eff. Mesh Size", fmt(pre$metrics$mesh), fmt(post$metrics$mesh), paste(pre$area_unit, "(effective connected habitat)", sep = " ")),
    sprintf("%-30s | %-15s | %-15s | %s", "Largest Patch Size", fmt(pre$metrics$lps), fmt(post$metrics$lps), pre$area_unit),
    sprintf("%-30s | %-15s | %-15s | %s", "Splitting Index", fmt(pre$metrics$split), fmt(post$metrics$split), "Higher = more fragmented"),
    sprintf("%-30s | %-15s | %-15s | %s", "Largest Patch Index (LPI)", fmt(pre$metrics$lpi), fmt(post$metrics$lpi), "% of landscape"),
    sprintf("%-30s | %-15s | %-15s | %s", "Mean Shape Index (MSI)", fmt(pre$metrics$msi), fmt(post$metrics$msi), "1 = compact, higher = irregular"),
    sprintf("%-30s | %-15s | %-15s | %s", "Edge Density (ED)", fmt(pre$metrics$ed), fmt(post$metrics$ed), paste(pre$dist_unit, "edge per", pre$area_unit)),
    sprintf("%-30s | %-15s | %-15s | %s", "Aggregation Index (AI)", fmt(pre$metrics$ai), fmt(post$metrics$ai), "0 = Dispersed, 100 = Clumped"),
    sprintf("%-30s | %-15s | %-15s | %s", "Core Area Index (CAI)", fmt(pre$metrics$cai), fmt(post$metrics$cai), "% core habitat")
  )
  results <- c(results, header)
  results
}

terralink_metrics_from_polygons <- function(sf_polys, units = "metric") {
  units <- match.arg(units, c("metric", "imperial"))
  if (!inherits(sf_polys, "sf")) {
    terralink_abort("sf_polys must be an sf object.", class = "terralink_error_input")
  }
  if (nrow(sf_polys) == 0) {
    return(list(metrics = list(total_area = 0, num_patches = 0, mesh = 0, lps = 0, split = 0, pd = 0, lpi = 0, total_edge = 0, msi = 0, ed = 0, mean_para = 0, frac = 0, ai = 0, cai = 0, mean_gyrate = 0), area_unit = ifelse(units == "imperial", "ac", "ha"), dist_unit = ifelse(units == "imperial", "ft", "m")))
  }

  area_div <- if (units == "imperial") 4046.8564224 else 10000.0
  dist_mult <- if (units == "imperial") 3.28084 else 1.0
  area_unit <- if (units == "imperial") "ac" else "ha"
  dist_unit <- if (units == "imperial") "ft" else "m"

  areas_m2 <- as.numeric(sf::st_area(sf_polys))
  perims_m <- as.numeric(sf::st_length(sf::st_cast(sf_polys, "MULTILINESTRING")))
  areas_unit <- areas_m2 / area_div
  perims_unit <- perims_m * dist_mult

  total_area <- sum(areas_unit)
  num_patches <- length(areas_unit)
  mesh <- if (total_area > 0) sum(areas_unit^2) / total_area else 0
  split <- if (sum(areas_unit^2) > 0) (total_area^2) / sum(areas_unit^2) else 0
  lps <- if (num_patches > 0) max(areas_unit) else 0
  lpi <- if (total_area > 0) (lps / total_area) * 100 else 0
  pd <- if (total_area > 0) (num_patches / total_area) * 100 else 0
  total_edge <- sum(perims_unit)
  ed <- if (total_area > 0) total_edge / total_area else 0
  msi <- if (num_patches > 0) mean(0.25 * perims_unit / sqrt(areas_unit), na.rm = TRUE) else 0
  mean_para <- if (num_patches > 0) mean(perims_unit / areas_unit, na.rm = TRUE) else 0
  frac_vals <- (2 * log(0.25 * perims_unit)) / log(areas_unit)
  frac <- mean(frac_vals[is.finite(frac_vals)], na.rm = TRUE)
  if (is.nan(frac)) frac <- 0

  list(
    metrics = list(
      total_area = total_area,
      num_patches = num_patches,
      mesh = mesh,
      lps = lps,
      split = split,
      pd = pd,
      lpi = lpi,
      total_edge = total_edge,
      msi = msi,
      ed = ed,
      mean_para = mean_para,
      frac = frac,
      ai = 0,
      cai = 0,
      mean_gyrate = 0
    ),
    area_unit = area_unit,
    dist_unit = dist_unit
  )
}

terralink_vector_report <- function(pre_polys, post_polys, units = "metric", label = "TerraLink", grid_resolution = 50) {
  units <- match.arg(units, c("metric", "imperial"))
  if (!inherits(pre_polys, "sf") || !inherits(post_polys, "sf")) {
    terralink_abort("pre_polys and post_polys must be sf objects.", class = "terralink_error_input")
  }

  if (nrow(pre_polys) == 0 && nrow(post_polys) == 0) {
    empty <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    pre_mask <- terra::ifel(empty > -Inf, 0, NA)
    post_mask <- pre_mask
    return(terralink_landscape_report(pre_mask, post_mask, empty, units = units, label = label))
  }

  work <- pre_polys
  if (nrow(work) == 0) work <- post_polys
  if (nrow(post_polys) > 0 && !isTRUE(sf::st_crs(post_polys) == sf::st_crs(work))) {
    post_polys <- sf::st_transform(post_polys, sf::st_crs(work))
  }

  geoms_all <- c(sf::st_geometry(pre_polys), sf::st_geometry(post_polys))
  all_bbox <- sf::st_bbox(sf::st_union(geoms_all))
  if (!all(is.finite(as.numeric(all_bbox)))) {
    all_bbox <- sf::st_bbox(work)
  }
  res_m <- max(1, as.numeric(grid_resolution))
  pad <- 2 * res_m
  template <- terra::rast(
    xmin = as.numeric(all_bbox["xmin"]) - pad,
    xmax = as.numeric(all_bbox["xmax"]) + pad,
    ymin = as.numeric(all_bbox["ymin"]) - pad,
    ymax = as.numeric(all_bbox["ymax"]) + pad,
    resolution = res_m,
    crs = sf::st_crs(work)$wkt
  )

  pre_mask <- template
  terra::values(pre_mask) <- NA
  if (nrow(pre_polys) > 0) {
    pre_mask <- terra::rasterize(terra::vect(pre_polys), template, field = 1, background = NA)
  }
  pre_mask <- terra::ifel(pre_mask > 0, 1, NA)

  post_mask <- template
  terra::values(post_mask) <- NA
  if (nrow(post_polys) > 0) {
    post_mask <- terra::rasterize(terra::vect(post_polys), template, field = 1, background = NA)
  }
  post_mask <- terra::ifel(post_mask > 0, 1, NA)

  terralink_landscape_report(pre_mask, post_mask, template, units = units, label = label)
}
