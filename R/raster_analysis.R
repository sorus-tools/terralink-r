terralink_candidate_col <- function(cand, name) {
  if (!name %in% names(cand)) return(NULL)
  val <- cand[[name]]
  if (is.list(val) && length(val) == 1) {
    inner <- val[[1]]
    if (is.atomic(inner) || inherits(inner, "sfc")) return(inner)
  }
  val
}

terralink_raster_endpoint_area_guard_ok <- function(patch1, patch2, corridor_cost, patch_sizes) {
  p1 <- suppressWarnings(as.integer(patch1))
  p2 <- suppressWarnings(as.integer(patch2))
  cost <- suppressWarnings(as.numeric(corridor_cost))
  if (!is.finite(p1) || !is.finite(p2) || p1 <= 0 || p2 <= 0 || p1 == p2) return(TRUE)
  if (!is.finite(cost) || cost <= 0) return(TRUE)
  p1_area <- suppressWarnings(as.numeric(patch_sizes[[as.character(p1)]] %||% 0))
  p2_area <- suppressWarnings(as.numeric(patch_sizes[[as.character(p2)]] %||% 0))
  if (is.finite(p1_area) && p1_area > 0 && p1_area + 1e-12 < cost) return(FALSE)
  if (is.finite(p2_area) && p2_area > 0 && p2_area + 1e-12 < cost) return(FALSE)
  TRUE
}

terralink_apply_hybrid_leftover_budget_raster <- function(candidates, selected_ids, patch_sizes, remaining_budget, roi_bias = 0.15) {
  if (!is.finite(remaining_budget) || remaining_budget <= 0 || nrow(candidates) == 0) {
    return(list(selected_ids = as.integer(selected_ids), budget_used = 0, low_value_added = 0L, redundancy_added = 0L))
  }

  patch_ids <- as.integer(names(patch_sizes))
  uf <- UnionFind$new()
  for (pid in patch_ids) {
    uf$find(pid)
    uf$size[[as.character(pid)]] <- as.numeric(patch_sizes[[as.character(pid)]])
    uf$count[[as.character(pid)]] <- 1L
  }

  g <- igraph::make_empty_graph(directed = FALSE)
  if (length(patch_ids) > 0) {
    g <- igraph::add_vertices(g, length(patch_ids), name = as.character(patch_ids))
  }

  selected_mask <- candidates$id %in% selected_ids
  if (any(selected_mask)) {
    for (i in which(selected_mask)) {
      p1 <- as.integer(candidates$patch1[[i]])
      p2 <- as.integer(candidates$patch2[[i]])
      uf$union(p1, p2)
      len <- as.numeric(candidates$length[[i]])
      if (!is.finite(len) || len <= 0) len <- as.numeric(candidates$cost[[i]])
      g <- terralink_graph_add_or_update(g, c(p1, p2), len)
    }
  }

  pick_main_root <- function() {
    if (length(patch_ids) == 0) return(NA_integer_)
    roots <- unique(vapply(patch_ids, function(pid) as.integer(uf$find(pid)), integer(1)))
    roots <- roots[is.finite(roots)]
    if (length(roots) == 0) return(NA_integer_)
    sizes <- vapply(roots, function(root) {
      val <- uf$size[[as.character(root)]]
      if (is.null(val) || !is.finite(val)) 0 else as.numeric(val)
    }, numeric(1))
    roots[[which.max(sizes)]]
  }

  main_root <- pick_main_root()
  budget_used <- 0
  low_value_added <- 0L
  redundancy_added <- 0L

  repeat {
    if (remaining_budget <= 0) break
    best_new <- NULL
    best_red <- NULL

    for (i in seq_len(nrow(candidates))) {
      if (selected_mask[[i]]) next
      cost <- as.numeric(candidates$cost[[i]])
      if (!is.finite(cost) || cost <= 0 || cost > remaining_budget) next
      p1 <- as.integer(candidates$patch1[[i]])
      p2 <- as.integer(candidates$patch2[[i]])
      r1 <- as.integer(uf$find(p1))
      r2 <- as.integer(uf$find(p2))

      if (r1 != r2) {
        if (!is.na(main_root) && xor(r1 == main_root, r2 == main_root)) {
          other <- if (r1 == main_root) r2 else r1
          gain <- uf$size[[as.character(other)]]
          if (is.null(gain) || !is.finite(gain) || gain <= 0) next
          roi <- as.numeric(gain) / cost
          if (is.null(best_new) || roi > best_new$roi) {
            best_new <- list(i = i, roi = roi)
          }
        }
        next
      }

      if (is.na(main_root) || r1 != main_root) next
      len <- as.numeric(candidates$length[[i]])
      if (!is.finite(len) || len <= 0) len <- cost
      score <- terralink_graph_shortcut_multiplier(
        g = g,
        p1 = p1,
        p2 = p2,
        length_val = len,
        diminishing_base = 0.5,
        shortcut_ratio_high = 3.0,
        shortcut_ratio_mid = 1.5,
        shortcut_ratio_low = 1.5,
        shortcut_mult_high = 0.9,
        shortcut_mult_mid = 0.5,
        shortcut_mult_low = 0.1
      )
      if (!is.finite(score) || score <= 0) next
      roi <- score / cost
      if (is.null(best_red) || roi > best_red$roi) {
        best_red <- list(i = i, roi = roi)
      }
    }

    if (is.null(best_new) && is.null(best_red)) break
    pick_red <- FALSE
    if (!is.null(best_red) && !is.null(best_new)) {
      pick_red <- best_red$roi > best_new$roi * (1 + as.numeric(roi_bias))
    } else if (!is.null(best_red)) {
      pick_red <- TRUE
    }
    pick <- if (pick_red) best_red else best_new
    if (is.null(pick)) break

    idx <- as.integer(pick$i)
    cost <- as.numeric(candidates$cost[[idx]])
    if (!is.finite(cost) || cost <= 0 || cost > remaining_budget) break
    p1 <- as.integer(candidates$patch1[[idx]])
    p2 <- as.integer(candidates$patch2[[idx]])
    uf$union(p1, p2)
    main_root <- pick_main_root()
    len <- as.numeric(candidates$length[[idx]])
    if (!is.finite(len) || len <= 0) len <- cost
    g <- terralink_graph_add_or_update(g, c(p1, p2), len)

    selected_mask[[idx]] <- TRUE
    remaining_budget <- remaining_budget - cost
    budget_used <- budget_used + cost
    if (pick_red) {
      redundancy_added <- redundancy_added + 1L
    } else {
      low_value_added <- low_value_added + 1L
    }
  }

  out_ids <- as.integer(candidates$id[selected_mask])
  list(
    selected_ids = out_ids,
    budget_used = as.numeric(budget_used),
    low_value_added = as.integer(low_value_added),
    redundancy_added = as.integer(redundancy_added)
  )
}

terralink_thicken_corridors_raster <- function(selected, remaining_budget, min_corridor_width_px, patch_sizes, labels, max_width_factor = 3L) {
  if (!is.finite(remaining_budget) || remaining_budget <= 0 || nrow(selected) == 0) {
    return(list(selected = selected, budget_used = 0))
  }
  if (!("path_cells" %in% names(selected) && "buffered_cells" %in% names(selected))) {
    return(list(selected = selected, budget_used = 0))
  }

  largest_patch <- NA_integer_
  if (length(patch_sizes) > 0) {
    nm <- names(patch_sizes)
    if (!is.null(nm) && length(nm) > 0) {
      largest_patch <- as.integer(nm[[which.max(as.numeric(patch_sizes))]])
    }
  }
  ordered <- seq_len(nrow(selected))
  if (is.finite(largest_patch)) {
    touches_largest <- (selected$patch1 == largest_patch) | (selected$patch2 == largest_patch)
    ordered <- c(which(touches_largest), which(!touches_largest))
  }

  nr <- terra::nrow(labels)
  nc <- terra::ncol(labels)
  label_vals <- terra::values(labels, mat = FALSE)
  inflate_ok <- is.na(label_vals) | (label_vals <= 0)
  budget_used <- 0

  if (!("thicken_factor" %in% names(selected))) selected$thicken_factor <- 1L

  for (idx in ordered) {
    if (remaining_budget <= 0) break
    path_cells <- selected$path_cells[[idx]]
    if (is.list(path_cells) && length(path_cells) == 1) path_cells <- path_cells[[1]]
    path_cells <- as.integer(path_cells)
    if (length(path_cells) == 0) next
    cur_buf <- selected$buffered_cells[[idx]]
    if (is.list(cur_buf) && length(cur_buf) == 1) cur_buf <- cur_buf[[1]]
    cur_buf <- as.integer(cur_buf)
    if (length(cur_buf) == 0) {
      cur_offsets <- terralink_corridor_offsets(min_corridor_width_px)
      cur_buf <- terralink_inflate_cells(path_cells, nrow = nr, ncol = nc, offsets = cur_offsets, passable_mask = inflate_ok)
    }
    cur_factor <- as.integer(selected$thicken_factor[[idx]] %||% 1L)
    if (!is.finite(cur_factor) || cur_factor < 1L) cur_factor <- 1L

    for (factor in seq.int(cur_factor + 1L, as.integer(max_width_factor))) {
      width <- max(1L, as.integer(round(as.numeric(min_corridor_width_px) * factor)))
      offsets <- terralink_corridor_offsets(width)
      new_buf <- terralink_inflate_cells(path_cells, nrow = nr, ncol = nc, offsets = offsets, passable_mask = inflate_ok)
      added <- length(setdiff(new_buf, cur_buf))
      if (added <= 0) {
        cur_buf <- new_buf
        cur_factor <- factor
        next
      }
      if (added > remaining_budget) {
        selected$buffered_cells[[idx]] <- as.integer(cur_buf)
        selected$cost[[idx]] <- as.numeric(length(cur_buf))
        selected$thicken_factor[[idx]] <- as.integer(cur_factor)
        return(list(selected = selected, budget_used = budget_used))
      }
      remaining_budget <- remaining_budget - added
      budget_used <- budget_used + added
      cur_buf <- new_buf
      cur_factor <- factor
    }
    selected$buffered_cells[[idx]] <- as.integer(cur_buf)
    selected$cost[[idx]] <- as.numeric(length(cur_buf))
    selected$thicken_factor[[idx]] <- as.integer(cur_factor)
  }

  list(selected = selected, budget_used = as.numeric(budget_used))
}

terralink_compute_raster_connected_sizes <- function(selected, patch_sizes) {
  if (is.null(selected) || nrow(selected) == 0) return(numeric(0))
  if (length(patch_sizes) == 0) return(rep(0, nrow(selected)))

  patch_ids <- suppressWarnings(as.integer(names(patch_sizes)))
  patch_ids <- patch_ids[is.finite(patch_ids)]
  uf <- UnionFind$new()
  for (pid in patch_ids) {
    uf$find(pid)
    uf$size[[as.character(pid)]] <- as.numeric(patch_sizes[[as.character(pid)]]) %||% 0
    uf$count[[as.character(pid)]] <- 1L
  }
  for (i in seq_len(nrow(selected))) {
    p1 <- suppressWarnings(as.integer(selected$patch1[[i]]))
    p2 <- suppressWarnings(as.integer(selected$patch2[[i]]))
    if (!is.finite(p1) || !is.finite(p2)) next
    uf$union(p1, p2)
  }

  comp_sizes <- list()
  for (pid in patch_ids) {
    root <- as.character(uf$find(pid))
    comp_sizes[[root]] <- as.numeric(comp_sizes[[root]] %||% 0) + as.numeric(patch_sizes[[as.character(pid)]] %||% 0)
  }
  for (i in seq_len(nrow(selected))) {
    p1 <- suppressWarnings(as.integer(selected$patch1[[i]]))
    cost <- as.numeric(selected$cost[[i]])
    if (!is.finite(p1) || !is.finite(cost) || cost <= 0) next
    root <- as.character(uf$find(p1))
    comp_sizes[[root]] <- as.numeric(comp_sizes[[root]] %||% 0) + cost
  }

  out <- numeric(nrow(selected))
  for (i in seq_len(nrow(selected))) {
    p1 <- suppressWarnings(as.integer(selected$patch1[[i]]))
    if (!is.finite(p1)) {
      out[[i]] <- 0
      next
    }
    root <- as.character(uf$find(p1))
    out[[i]] <- as.numeric(comp_sizes[[root]] %||% 0)
  }
  out
}

terralink_normalize_raster_backend <- function(backend, default = "delegate_vector") {
  raw <- default
  if (!is.null(backend) && length(backend) > 0) {
    candidate <- as.character(backend[[1]])
    if (!is.na(candidate) && nzchar(trimws(candidate))) raw <- candidate
  }
  key <- tolower(trimws(raw))
  key <- gsub("[[:space:]-]+", "_", key)
  aliases <- c(
    delegate = "delegate_vector",
    vector = "delegate_vector",
    vector_polygonize = "delegate_vector",
    qgis = "delegate_vector",
    current_qgis = "delegate_vector",
    native = "native",
    raster = "native"
  )
  if (key %in% names(aliases)) key <- unname(aliases[[key]])
  if (!key %in% c("delegate_vector", "native")) key <- default
  key
}

terralink_raster_bridge_pixel_metrics <- function(raster) {
  if (isTRUE(terra::is.lonlat(raster))) {
    terralink_abort(
      "Raster-to-vector bridge requires a projected raster CRS measured in meters.",
      class = "terralink_error_scale",
      fix = c("Reproject the raster to a projected CRS")
    )
  }
  wkt <- tryCatch(terra::crs(raster, proj = TRUE), error = function(e) "")
  if (!nzchar(wkt)) {
    terralink_abort(
      "Raster-to-vector bridge requires a raster CRS.",
      class = "terralink_error_input",
      fix = c("Assign a projected CRS to the raster")
    )
  }
  crs_info <- tryCatch(sf::st_crs(wkt), error = function(e) NULL)
  unit_name <- tolower(as.character(crs_info$units_gdal %||% ""))
  if (nzchar(unit_name) && !grepl("met", unit_name, fixed = FALSE)) {
    terralink_abort(
      "Raster-to-vector bridge requires projected coordinates measured in meters.",
      class = "terralink_error_scale",
      details = sprintf("CRS units: %s", unit_name),
      fix = c("Reproject the raster to a meter-based projected CRS")
    )
  }

  res_xy <- tryCatch(terra::res(raster), error = function(e) c(NA_real_, NA_real_))
  pixel_w_m <- abs(as.numeric(res_xy[[1]]))
  pixel_h_m <- abs(as.numeric(res_xy[[2]]))
  pixel_size_m <- suppressWarnings(max(pixel_w_m, pixel_h_m, na.rm = TRUE))
  pixel_area_m2 <- pixel_w_m * pixel_h_m
  if (!is.finite(pixel_size_m) || pixel_size_m <= 0 || !is.finite(pixel_area_m2) || pixel_area_m2 <= 0) {
    terralink_abort(
      "Unable to derive valid pixel size from the raster.",
      class = "terralink_error_scale"
    )
  }

  list(
    pixel_w_m = pixel_w_m,
    pixel_h_m = pixel_h_m,
    pixel_size_m = pixel_size_m,
    pixel_area_m2 = pixel_area_m2,
    pixel_area_ha = pixel_area_m2 / 10000.0
  )
}

terralink_raster_labels_to_sf <- function(labels, field_name = "patch_id") {
  vals <- tryCatch(terra::values(labels, mat = FALSE), error = function(e) numeric(0))
  if (length(vals) == 0 || !any(is.finite(vals) & vals > 0, na.rm = TRUE)) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = terra::crs(labels))))
  }
  polys <- terra::as.polygons(labels, values = TRUE, na.rm = TRUE)
  polys_sf <- sf::st_as_sf(polys)
  geom_col <- attr(polys_sf, "sf_column")
  value_cols <- setdiff(names(polys_sf), geom_col)
  if (length(value_cols) == 0) {
    polys_sf[[field_name]] <- seq_len(nrow(polys_sf))
  } else {
    names(polys_sf)[names(polys_sf) == value_cols[[1]]] <- field_name
  }
  polys_sf <- polys_sf[!is.na(polys_sf[[field_name]]) & polys_sf[[field_name]] > 0, , drop = FALSE]
  polys_sf[[field_name]] <- as.integer(polys_sf[[field_name]])
  polys_sf
}

terralink_rasterize_sf_to_template <- function(sf_obj, template, field = NULL) {
  out <- template
  terra::values(out) <- NA
  if (is.null(sf_obj) || !inherits(sf_obj, "sf") || nrow(sf_obj) == 0) return(out)
  vec <- terra::vect(sf_obj)
  if (!is.null(field) && field %in% names(sf_obj)) {
    return(terra::rasterize(vec, template, field = field, background = NA))
  }
  terra::rasterize(vec, template, field = 1, background = NA)
}

terralink_convert_named_numeric_units <- function(x, factor, patterns) {
  if (!is.list(x) || length(x) == 0 || !is.finite(factor) || factor <= 0) return(x)
  nm <- names(x)
  if (is.null(nm)) return(x)
  for (i in seq_along(x)) {
    key <- nm[[i]]
    if (!nzchar(key) || !is.numeric(x[[i]]) || length(x[[i]]) != 1 || !is.finite(x[[i]])) next
    if (any(vapply(patterns, function(pattern) grepl(pattern, key, perl = TRUE), logical(1)))) {
      x[[i]] <- as.numeric(x[[i]]) * factor
    }
  }
  x
}

terralink_raster_delegate_result <- function(
  raster,
  strategy_key,
  patch_values,
  patch_ranges,
  budget,
  min_patch_size,
  min_corridor_width,
  corridor_cell_assignment_key,
  max_search_distance,
  obstacle_values,
  obstacle_ranges,
  patch_connectivity,
  units,
  max_pair_checks,
  max_candidates,
  verbose,
  progress,
  obstacle_strategy,
  species_dispersal_distance,
  species_dispersal_kernel,
  min_patch_area_for_species,
  patch_area_scaling,
  mobility_detour_cap,
  redundancy_method,
  metric_weights,
  weight_m,
  weight_lcc,
  weight_pc,
  weight_f,
  keep_candidates,
  preflight,
  t_start
) {
  pixel_metrics <- terralink_raster_bridge_pixel_metrics(raster)
  unit_conv <- terralink_convert_raster_units(
    raster,
    units = units,
    budget = budget,
    min_patch_size = min_patch_size,
    min_corridor_width = min_corridor_width,
    max_search_distance = max_search_distance
  )

  budget_px <- unit_conv$budget_pixels
  min_patch_size_px <- unit_conv$min_patch_size_px
  min_corridor_width_px <- unit_conv$min_corridor_width_px
  max_search_distance_px <- unit_conv$max_search_distance_px

  habitat_mask <- terralink_mask_from_values_ranges(raster, patch_values, patch_ranges)
  if (is.null(habitat_mask)) {
    value_summary <- terralink_raster_value_summary(raster)
    details <- c(
      sprintf("patch_values: %s", paste(patch_values %||% character(0), collapse = ", ")),
      sprintf("patch_ranges: %s", if (is.null(patch_ranges)) "NULL" else paste(vapply(patch_ranges, function(r) paste(r, collapse = "-"), character(1)), collapse = "; "))
    )
    if (!is.null(value_summary)) {
      details <- c(details, sprintf("Sampled values (%s cells): %s", format(value_summary$sampled, big.mark = ","), paste(value_summary$items, collapse = "; ")))
    }
    terralink_abort(
      "Unable to create habitat mask from patch values/ranges.",
      class = "terralink_error_input",
      details = details,
      fix = c("Check categorical codes", "Ensure ranges match raster values", "Handle NA values explicitly")
    )
  }

  obstacle_mask <- terralink_mask_from_values_ranges(raster, obstacle_values, obstacle_ranges)
  patch_cells <- terra::global(habitat_mask, fun = "sum", na.rm = TRUE)
  patch_cells <- if (!is.null(patch_cells)) as.numeric(patch_cells[1, 1]) else NA_real_
  if (!is.null(obstacle_mask)) {
    habitat_mask[obstacle_mask == 1] <- FALSE
    obstacle_cells <- terra::global(obstacle_mask, fun = "sum", na.rm = TRUE)
    obstacle_cells <- if (!is.null(obstacle_cells)) as.numeric(obstacle_cells[1, 1]) else NA_real_
  } else {
    obstacle_cells <- 0
  }

  labels <- label_patches(habitat_mask, connectivity = patch_connectivity)
  patch_df <- patch_summary_from_labels(labels)
  if (nrow(patch_df) == 0) {
    terralink_abort(
      "No patches found for the given patch_values/ranges.",
      class = "terralink_error_input",
      fix = c("Verify raster codes", "Adjust patch_values/patch_ranges", "Check obstacle masks")
    )
  }

  habitat_sf <- terralink_raster_labels_to_sf(labels, field_name = "patch_id")
  obstacle_sf <- NULL
  if (!is.null(obstacle_mask)) {
    obstacle_labels <- label_patches(terra::ifel(obstacle_mask == 1, 1, NA), connectivity = patch_connectivity)
    obstacle_sf <- terralink_raster_labels_to_sf(obstacle_labels, field_name = "obstacle_id")
    if (nrow(obstacle_sf) == 0) obstacle_sf <- NULL
  }

  delegate_units <- if (units == "imperial") "imperial" else "metric"
  budget_delegate <- if (units == "pixels") budget_px * pixel_metrics$pixel_area_ha else budget
  min_patch_delegate <- if (units == "pixels") min_patch_size_px * pixel_metrics$pixel_area_ha else min_patch_size
  min_corridor_width_delegate <- if (units == "pixels") min_corridor_width_px * pixel_metrics$pixel_size_m else min_corridor_width
  max_search_distance_delegate <- if (units == "pixels") max_search_distance_px * pixel_metrics$pixel_size_m else max_search_distance
  species_dispersal_delegate <- species_dispersal_distance
  if (!is.null(species_dispersal_delegate) && units == "pixels") {
    species_dispersal_delegate <- as.numeric(species_dispersal_delegate) * pixel_metrics$pixel_size_m
  }
  min_patch_area_for_species_delegate <- if (units == "pixels") {
    as.numeric(min_patch_area_for_species %||% 0) * pixel_metrics$pixel_area_ha
  } else {
    min_patch_area_for_species
  }

  vector_result <- run_vector_analysis(
    patches = habitat_sf,
    budget = budget_delegate,
    strategy = strategy_key,
    min_patch_size = min_patch_delegate,
    min_corridor_width = min_corridor_width_delegate,
    max_search_distance = max_search_distance_delegate,
    obstacle_layers = obstacle_sf,
    obstacle_resolution = max(pixel_metrics$pixel_size_m, 1),
    units = delegate_units,
    max_pair_checks = max_pair_checks,
    max_candidates = max_candidates,
    verbose = verbose,
    progress = progress,
    obstacle_strategy = obstacle_strategy,
    return_crs = "input",
    species_dispersal_distance = species_dispersal_delegate,
    species_dispersal_kernel = species_dispersal_kernel,
    min_patch_area_for_species = min_patch_area_for_species_delegate,
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

  corridors_out <- vector_result$corridors
  networks_out <- vector_result$networks
  area_to_pixels <- 10000.0 / pixel_metrics$pixel_area_m2
  distance_to_pixels <- 1.0 / pixel_metrics$pixel_size_m

  if (!is.null(corridors_out) && inherits(corridors_out, "sf") && nrow(corridors_out) > 0) {
    if (units == "pixels") {
      corridors_out$corridor_area <- as.numeric(corridors_out$corridor_area_m2) * area_to_pixels
      corridors_out$connected_area <- as.numeric(corridors_out$connected_area_m2) * area_to_pixels
      corridors_out$corridor_length <- as.numeric(corridors_out$corridor_length_m) * distance_to_pixels
    }
    if (corridor_cell_assignment_key == "sum_direct_connected_patches") {
      patch_lookup <- stats::setNames(as.numeric(patch_df$cell_count), as.character(patch_df$patch_id))
      corridors_out$direct_patch_area <- vapply(seq_len(nrow(corridors_out)), function(i) {
        p1 <- as.character(corridors_out$patch1[[i]])
        p2 <- as.character(corridors_out$patch2[[i]])
        a1 <- as.numeric(patch_lookup[[p1]] %||% 0)
        a2 <- as.numeric(patch_lookup[[p2]] %||% 0)
        a1 + a2
      }, numeric(1))
    }
  }
  if (!is.null(networks_out) && inherits(networks_out, "sf") && nrow(networks_out) > 0 && units == "pixels") {
    networks_out$area <- as.numeric(networks_out$area_m2) * area_to_pixels
  }

  corridor_field <- NULL
  if (!is.null(corridors_out) && inherits(corridors_out, "sf") && nrow(corridors_out) > 0) {
    if (corridor_cell_assignment_key == "efficiency" && "efficiency" %in% names(corridors_out)) {
      corridor_field <- "efficiency"
    } else if (corridor_cell_assignment_key == "sum_direct_connected_patches" && "direct_patch_area" %in% names(corridors_out)) {
      corridor_field <- "direct_patch_area"
    } else if ("connected_area" %in% names(corridors_out)) {
      corridor_field <- "connected_area"
    }
  }
  corridor_raster <- terralink_rasterize_sf_to_template(corridors_out, labels, field = corridor_field)
  contiguous_field <- if (!is.null(networks_out) && inherits(networks_out, "sf") && "area" %in% names(networks_out)) "area" else NULL
  contiguous_raster <- if (!is.null(networks_out) && inherits(networks_out, "sf") && nrow(networks_out) > 0) {
    terralink_rasterize_sf_to_template(networks_out, labels, field = contiguous_field)
  } else {
    NULL
  }

  summary <- vector_result$summary %||% list()
  summary$budget_total <- budget
  if (is.numeric(summary$budget_used) && length(summary$budget_used) == 1 && is.finite(summary$budget_used) && units == "pixels") {
    summary$budget_used <- as.numeric(summary$budget_used) * area_to_pixels
  }
  summary$patch_values <- patch_values
  summary$patch_ranges <- patch_ranges
  summary$obstacle_values <- obstacle_values
  summary$obstacle_ranges <- obstacle_ranges
  summary$patch_value_cells <- patch_cells
  summary$obstacle_cells <- obstacle_cells
  summary$min_patch_size <- min_patch_size
  summary$min_corridor_width <- min_corridor_width
  summary$corridor_cell_assignment <- corridor_cell_assignment_key
  summary$max_search_distance <- max_search_distance
  summary$patch_connectivity <- patch_connectivity
  summary$units <- units
  summary$raster_delegate_backend <- "vector_polygonize"
  summary$source_raster_rows <- terra::nrow(raster)
  summary$source_raster_cols <- terra::ncol(raster)
  summary$source_raster_pixels_total <- terra::ncell(raster)
  summary$source_raster_pixel_width_m <- pixel_metrics$pixel_w_m
  summary$source_raster_pixel_height_m <- pixel_metrics$pixel_h_m
  summary$source_raster_pixel_area_ha <- pixel_metrics$pixel_area_ha

  metrics <- vector_result$metrics %||% list()
  strategy_stats <- vector_result$strategy_stats %||% list()
  if (units == "pixels") {
    convert_patterns <- c(
      "(^|_)budget(_|$)",
      "(^|_)area(_|$)",
      "_cluster(_|$)",
      "^habitat_availability",
      "^mean_reachable_area",
      "^median_reachable_area",
      "^largest_group_area",
      "^total_connected_area",
      "^total_patch_area",
      "^seed_area",
      "^final_patch_area",
      "^bigconnect_objective"
    )
    metrics <- terralink_convert_named_numeric_units(metrics, area_to_pixels, convert_patterns)
    strategy_stats <- terralink_convert_named_numeric_units(strategy_stats, area_to_pixels, convert_patterns)
  }

  terralink_as_result(
    list(
      patches = labels,
      patch_table = patch_df,
      corridors = corridors_out,
      networks = networks_out,
      corridor_raster = corridor_raster,
      contiguous_raster = contiguous_raster,
      strategy = strategy_key,
      summary = summary,
      metrics = metrics,
      metrics_report = vector_result$metrics_report,
      strategy_stats = strategy_stats,
      delegate_result = vector_result
    ),
    mode = "raster",
    inputs = list(
      units = units,
      pixel_size = unit_conv$pixel_size,
      pixel_area = unit_conv$pixel_area,
      raster_cells = preflight$n_cells,
      budget = budget,
      min_patch_size = min_patch_size,
      min_corridor_width = min_corridor_width,
      max_search_distance = max_search_distance,
      raster_backend = "delegate_vector"
    ),
    run_stats = list(
      elapsed_s = proc.time()[[3]] - t_start,
      candidate_edges = vector_result$run_stats$candidate_edges %||% summary$candidate_edges %||% 0,
      candidate_pairs = vector_result$run_stats$candidate_pairs %||% summary$candidate_pairs %||% 0,
      possible_pairs = vector_result$run_stats$possible_pairs %||% summary$possible_pairs %||% 0
    ),
    warnings = unique(c(vector_result$warnings %||% character(0))),
    diagnostics = vector_result$diagnostics %||% list()
  )
}

#' Run TerraLink raster workflow
#'
#' @param raster SpatRaster or path to raster.
#' @param patch_values Numeric values representing habitat.
#' @param patch_ranges Optional list of value ranges defining habitat.
#' @param budget Total corridor budget (units defined by `units`).
#' @param budget_pixels Back-compat alias for budget (pixels).
#' @param strategy Strategy name. Canonical TerraLink 1.7 values are
#'   "most_connected_habitat", "largest_single_network",
#'   "landscape_fluidity", and "reachable_habitat_advanced".
#' @param min_patch_size Minimum patch size (units defined by `units`).
#' @param min_corridor_width Minimum corridor width (units defined by `units`).
#' @param corridor_cell_assignment Corridor cell assignment mode.
#' @param max_search_distance Maximum search distance (units defined by `units`).
#' @param obstacle_values Optional impassable raster values.
#' @param obstacle_ranges Optional list of impassable ranges.
#' @param allow_bottlenecks Whether to allow corridors to squeeze through gaps.
#' @param patch_connectivity Connectivity for patch labeling (4 or 8).
#' @param units Unit system: "pixels", "metric", or "imperial".
#' @param allow_large Allow processing very large rasters.
#' @param max_pair_checks Limit for candidate pair checks (prevents O(n^2) blowups).
#' @param max_candidates Limit for candidate corridors.
#' @param verbose Verbosity level (0-2).
#' @param progress Show progress bars.
#' @param obstacle_strategy Behavior when gdistance is unavailable and obstacles are provided.
#' @param species_dispersal_distance Species movement distance used by
#'   "reachable_habitat_advanced" and connectivity reporting.
#' @param species_dispersal_kernel Dispersal kernel for habitat availability.
#' @param min_patch_area_for_species Minimum patch area eligible for species metrics.
#' @param patch_area_scaling Patch-area scaling for habitat availability ("sqrt" or "log").
#' @param mobility_detour_cap Cap used by graph-based mobility/fluidity metrics.
#' @param redundancy_method Flow redundancy method ("ime" or "fri").
#' @param metric_weights Named numeric vector for composite connectivity score.
#' @param weight_m Optional mesh weight override for composite score.
#' @param weight_lcc Optional LCC weight override for composite score.
#' @param weight_pc Optional PC weight override for composite score.
#' @param weight_f Optional flow weight override for composite score.
#' @param keep_candidates Whether to keep candidate list in the output.
#' @return List with patches, corridors, rasters, and summary.
#' @details
#' Raster inputs are funneled through TerraLink's vector corridor pipeline after
#' habitat and impassable classes are polygonized, matching the current QGIS
#' plugin workflow.
#' @export
run_raster_analysis <- function(
  raster,
  patch_values,
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
  patch_connectivity = 4,
  units = "pixels",
  patch_ranges = NULL,
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
  keep_candidates = FALSE
) {
  raster <- terralink_resolve_raster(raster)
  ctx <- terralink_new_run_context(verbose = verbose, progress = progress)
  terralink_progress_start(ctx, message = "Starting raster analysis")
  t_start <- proc.time()[[3]]
  preflight <- terralink_preflight_raster(raster, allow_large = allow_large, ctx = ctx)
  strategy_key <- terralink_normalize_strategy_key(strategy, default = "most_connected_habitat")
  corridor_cell_assignment_key <- terralink_normalize_corridor_cell_assignment(corridor_cell_assignment)

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
  if (!is.numeric(budget) || budget <= 0) {
    terralink_abort("budget must be a positive number.", class = "terralink_error_input")
  }
  if (!is.numeric(max_search_distance) || max_search_distance <= 0) {
    terralink_abort("max_search_distance must be a positive number.", class = "terralink_error_input")
  }
  if ((is.null(patch_values) || length(patch_values) == 0) && (is.null(patch_ranges) || length(patch_ranges) == 0)) {
    terralink_abort(
      "patch_values or patch_ranges must be provided.",
      class = "terralink_error_input",
      fix = c("Provide patch_values (e.g., c(1, 2))", "Provide patch_ranges (e.g., list(c(1, 3)))")
    )
  }

  result <- terralink_raster_delegate_result(
    raster = raster,
    strategy_key = strategy_key,
    patch_values = patch_values,
    patch_ranges = patch_ranges,
    budget = budget,
    min_patch_size = min_patch_size,
    min_corridor_width = min_corridor_width,
    corridor_cell_assignment_key = corridor_cell_assignment_key,
    max_search_distance = max_search_distance,
    obstacle_values = obstacle_values,
    obstacle_ranges = obstacle_ranges,
    patch_connectivity = patch_connectivity,
    units = units,
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
    keep_candidates = keep_candidates,
    preflight = preflight,
    t_start = t_start
  )
  terralink_progress_done(ctx)
  return(result)

  obstacle_strategy <- match.arg(obstacle_strategy)
  has_obstacles <- (!is.null(obstacle_values) && length(obstacle_values) > 0) || (!is.null(obstacle_ranges) && length(obstacle_ranges) > 0)
  has_gdistance <- requireNamespace("gdistance", quietly = TRUE) &&
    requireNamespace("raster", quietly = TRUE) &&
    requireNamespace("sp", quietly = TRUE)
  if (has_obstacles && !has_gdistance) {
    if (obstacle_strategy == "error") {
      terralink_abort(
        "Obstacle-aware routing requires the gdistance package.",
        class = "terralink_error_dependency",
        fix = c("Install gdistance, raster, and sp", "Set obstacle_strategy = 'straight_line' or 'disable_obstacles'")
      )
    } else if (obstacle_strategy == "disable_obstacles") {
      terralink_warn("Obstacle layers disabled because gdistance is not available.", ctx = ctx)
      obstacle_values <- NULL
      obstacle_ranges <- NULL
    } else {
      terralink_warn("gdistance not available; using straight-line routing without obstacle awareness.", ctx = ctx)
    }
  }

  terralink_progress_update(ctx, 10, "Converting units")
  unit_conv <- terralink_convert_raster_units(
    raster,
    units = units,
    budget = budget,
    min_patch_size = min_patch_size,
    min_corridor_width = min_corridor_width,
    max_search_distance = max_search_distance
  )

  budget_px <- unit_conv$budget_pixels
  min_patch_size_px <- unit_conv$min_patch_size_px
  min_corridor_width_px <- unit_conv$min_corridor_width_px
  max_search_distance_px <- unit_conv$max_search_distance_px
  area_div <- if (units == "imperial") 4046.8564224 else if (units == "metric") 10000.0 else 1.0
  dist_mult_report <- if (units == "imperial") 3.28084 else 1.0
  pixel_size_report <- if (units == "pixels") unit_conv$pixel_size else 1.0

  terralink_progress_update(ctx, 20, "Building habitat mask")
  habitat_mask <- terralink_mask_from_values_ranges(raster, patch_values, patch_ranges)
  if (is.null(habitat_mask)) {
    value_summary <- terralink_raster_value_summary(raster)
    details <- c(
      sprintf("patch_values: %s", paste(patch_values %||% character(0), collapse = ", ")),
      sprintf("patch_ranges: %s", if (is.null(patch_ranges)) "NULL" else paste(vapply(patch_ranges, function(r) paste(r, collapse = "-"), character(1)), collapse = "; "))
    )
    if (!is.null(value_summary)) {
      details <- c(details, sprintf("Sampled values (%s cells): %s", format(value_summary$sampled, big.mark = ","), paste(value_summary$items, collapse = "; ")))
    }
    terralink_abort(
      "Unable to create habitat mask from patch values/ranges.",
      class = "terralink_error_input",
      details = details,
      fix = c("Check categorical codes", "Ensure ranges match raster values", "Handle NA values explicitly")
    )
  }

  obstacle_mask <- terralink_mask_from_values_ranges(raster, obstacle_values, obstacle_ranges)
  patch_cells <- terra::global(habitat_mask, fun = "sum", na.rm = TRUE)
  patch_cells <- if (!is.null(patch_cells)) as.numeric(patch_cells[1, 1]) else NA_real_
  if (!is.null(obstacle_mask)) {
    habitat_mask[obstacle_mask == 1] <- FALSE
    obstacle_cells <- terra::global(obstacle_mask, fun = "sum", na.rm = TRUE)
    obstacle_cells <- if (!is.null(obstacle_cells)) as.numeric(obstacle_cells[1, 1]) else NA_real_
  } else {
    obstacle_cells <- 0
  }

  labels <- label_patches(habitat_mask, connectivity = patch_connectivity)
  freq <- terra::freq(labels)
  if (is.null(freq) || nrow(freq) == 0) {
    value_summary <- terralink_raster_value_summary(raster)
    details <- c(
      sprintf("patch_values: %s", paste(patch_values %||% character(0), collapse = ", ")),
      sprintf("patch_ranges: %s", if (is.null(patch_ranges)) "NULL" else paste(vapply(patch_ranges, function(r) paste(r, collapse = "-"), character(1)), collapse = "; ")),
      sprintf("obstacle_values: %s", paste(obstacle_values %||% character(0), collapse = ", "))
    )
    if (!is.null(value_summary)) {
      details <- c(details, sprintf("Sampled values (%s cells): %s", format(value_summary$sampled, big.mark = ","), paste(value_summary$items, collapse = "; ")))
      if (value_summary$na_count > 0) {
        details <- c(details, sprintf("NA values in sample: %s", format(value_summary$na_count, big.mark = ",")))
      }
    }
    terralink_abort(
      "No patches found for the given patch_values/ranges.",
      class = "terralink_error_input",
      details = details,
      fix = c("Verify raster codes", "Adjust patch_values/patch_ranges", "Check NA handling and obstacle masks")
    )
  }
  freq <- freq[!is.na(freq$value) & freq$value != 0, , drop = FALSE]
  raw_patch_count <- nrow(freq)
  valid_ids <- freq$value[freq$count >= min_patch_size_px]
  if (length(valid_ids) == 0) {
    detail_counts <- if (nrow(freq) > 0) {
      stats::quantile(freq$count, probs = c(0, 0.5, 0.9, 0.99), na.rm = TRUE)
    } else {
      NA
    }
    terralink_abort(
      "No patches meet min_patch_size.",
      class = "terralink_error_scale",
      details = c(
        sprintf("Raw patches: %s", format(raw_patch_count, big.mark = ",")),
        sprintf("min_patch_size (px): %s", terralink_format_number(min_patch_size_px, 2)),
        sprintf("Patch size quantiles (px): %s", paste(round(detail_counts, 2), collapse = ", "))
      ),
      fix = c("Lower min_patch_size", "Increase pixel size", "Adjust patch_values")
    )
  }
  filtered_out <- raw_patch_count - length(valid_ids)

  valid_ids <- as.numeric(valid_ids)
  valid_mask <- terra::app(labels, fun = function(x) x %in% valid_ids)
  labels[valid_mask == 0] <- NA
  patch_df <- patch_summary_from_labels(labels)
  possible_pairs <- if (nrow(patch_df) > 1) {
    as.integer((nrow(patch_df) * (nrow(patch_df) - 1)) / 2)
  } else {
    0L
  }
  terralink_check_candidate_count(
    possible_pairs,
    max_candidates = max_pair_checks,
    ctx = ctx,
    scope = "Candidate pairs",
    override_param = "max_pair_checks"
  )
  terralink_inform(
    sprintf("Patches labeled: %s (raw %s, filtered %s)", nrow(patch_df), raw_patch_count, filtered_out),
    ctx = ctx,
    level = 1
  )

  # Build passable and blocked masks using the same semantics as QGIS:
  # filtered patch pixels are blocked, while matrix cells remain passable.
  filtered_patch_mask <- terra::ifel(!is.na(labels) & labels > 0, 1, 0)
  obstacle_for_traversal <- filtered_patch_mask
  if (!is.null(obstacle_mask)) {
    obstacle_for_traversal <- terra::ifel(obstacle_mask == 1 | filtered_patch_mask == 1, 1, 0)
  }
  passable_mask <- terra::ifel(filtered_patch_mask == 1 | obstacle_for_traversal == 1, 0, 1)
  if (!allow_bottlenecks && min_corridor_width_px > 1) {
    true_obstacles <- terra::ifel(obstacle_for_traversal == 1 & filtered_patch_mask == 0, 1, 0)
    clearance_space <- terra::ifel(true_obstacles == 1, 0, 1)
    kernel <- terralink_offsets_kernel(terralink_corridor_offsets(min_corridor_width_px))
    clearance_ok <- terra::focal(clearance_space, w = kernel, fun = min, na.policy = "omit", fillvalue = 0)
    passable_mask <- terra::ifel(passable_mask == 1 & clearance_ok == 1, 1, 0)
  }

  terralink_progress_update(ctx, 55, "Generating corridor candidates")
  candidates <- build_raster_candidates(
    labels = labels,
    patch_df = patch_df,
    passable_mask = passable_mask,
    max_search_distance_px = max_search_distance_px,
    raster_ref = raster,
    min_corridor_width_px = min_corridor_width_px,
    pair_index = NULL,
    patch_connectivity = patch_connectivity,
    habitat_mask = habitat_mask,
    obstacle_mask = obstacle_for_traversal
  )

  if (nrow(candidates) > 0) {
    pair_keys <- paste(pmin(candidates$patch1, candidates$patch2), pmax(candidates$patch1, candidates$patch2), sep = "_")
    candidate_pairs_count <- length(unique(pair_keys))
  } else {
    candidate_pairs_count <- 0L
  }

  terralink_check_candidate_count(nrow(candidates), max_candidates = max_candidates, ctx = ctx)
  terralink_inform(sprintf("Candidate corridors generated: %s", nrow(candidates)), ctx = ctx, level = 1)

  if (nrow(candidates) == 0) {
    summary <- list(
      budget_total = budget,
      budget_used = 0,
      corridors_used = 0,
      candidate_edges = 0,
      candidate_pairs = candidate_pairs_count,
      possible_pairs = possible_pairs,
      patches = nrow(patch_df),
      raw_patches = raw_patch_count,
      filtered_out = filtered_out,
      patch_values = patch_values,
      patch_ranges = patch_ranges,
      obstacle_values = obstacle_values,
      obstacle_ranges = obstacle_ranges,
      patch_value_cells = patch_cells,
      obstacle_cells = obstacle_cells,
      min_patch_size = min_patch_size,
      min_corridor_width = min_corridor_width,
      corridor_cell_assignment = corridor_cell_assignment_key,
      max_search_distance = max_search_distance,
      patch_connectivity = patch_connectivity,
      units = units
    )
    result <- list(
      patches = labels,
      patch_table = patch_df,
      corridors = candidates,
      corridor_raster = build_corridor_raster(
        labels,
        patch_df,
        candidates,
        min_corridor_width_px = min_corridor_width_px,
        assignment_mode = corridor_cell_assignment_key
      ),
      contiguous_raster = NULL,
      strategy = strategy_key,
      summary = summary
    )
    result <- terralink_as_result(
      result,
      mode = "raster",
      inputs = list(
        units = units,
        pixel_size = unit_conv$pixel_size,
        pixel_area = unit_conv$pixel_area,
        raster_cells = preflight$n_cells,
        budget = budget,
        min_patch_size = min_patch_size,
        min_corridor_width = min_corridor_width,
        max_search_distance = max_search_distance
      ),
      run_stats = list(
        elapsed_s = proc.time()[[3]] - t_start,
        candidate_edges = 0,
        candidate_pairs = candidate_pairs_count,
        possible_pairs = possible_pairs
      ),
      warnings = ctx$warnings,
      diagnostics = list(message = "No candidates generated; try increasing max_search_distance or lowering min_patch_size.")
    )
    terralink_progress_done(ctx)
    return(result)
  }

  nodes <- stats::setNames(patch_df$cell_count, patch_df$patch_id)
  engine_edges <- data.frame(
    u = candidates$patch1,
    v = candidates$patch2,
    id = candidates$id,
    cost = candidates$cost
  )
  patch_metric_df <- data.frame(
    patch_id = patch_df$patch_id,
    area = if (units == "pixels") patch_df$cell_count else (patch_df$cell_count * unit_conv$pixel_area / area_div),
    x = if (units == "pixels") (patch_df$x / unit_conv$pixel_size) else (patch_df$x * dist_mult_report),
    y = if (units == "pixels") (patch_df$y / unit_conv$pixel_size) else (patch_df$y * dist_mult_report),
    quality_weight = rep(1, nrow(patch_df)),
    stringsAsFactors = FALSE
  )
  candidates$cost_metric <- if (units == "pixels") candidates$cost else (candidates$cost * unit_conv$pixel_area / area_div)
  candidates$length_metric <- if (units == "pixels") candidates$length else (candidates$length * unit_conv$pixel_size * dist_mult_report)

  strategy_stats <- list(strategy = strategy_key)
  if (strategy_key == "largest_single_network") {
    patch_optimizer_df <- data.frame(
      patch_id = patch_df$patch_id,
      area = patch_df$cell_count,
      x = patch_df$x / unit_conv$pixel_size,
      y = patch_df$y / unit_conv$pixel_size,
      stringsAsFactors = FALSE
    )
    opt <- terralink_optimize_largest_single_network_plugin_parity_df(
      candidates = candidates,
      patch_df = patch_optimizer_df,
      budget = budget_px,
      cost_col = "cost",
      length_col = "length"
    )
    selected_ids <- opt$selected
    budget_used <- as.numeric(opt$stats$budget_used %||% 0)
    strategy_stats <- opt$stats
  } else if (strategy_key == "reachable_habitat_advanced") {
    opt <- terralink_optimize_habitat_availability_df(
      candidates = candidates,
      patch_df = patch_metric_df,
      budget = if (units == "pixels") budget_px else budget,
      min_patch_area_for_species = min_patch_area_for_species,
      patch_area_scaling = patch_area_scaling,
      species_dispersal_distance = species_dispersal_distance %||% max_search_distance,
      species_dispersal_kernel = species_dispersal_kernel
    )
    selected_ids <- opt$selected
    budget_used <- if (units == "pixels") {
      as.numeric(opt$stats$budget_used %||% 0)
    } else {
      as.numeric(opt$stats$budget_used %||% 0) * area_div / unit_conv$pixel_area
    }
    strategy_stats <- opt$stats
  } else if (strategy_key == "landscape_fluidity") {
    opt <- terralink_optimize_landscape_fluidity_df(
      candidates = candidates,
      patch_df = patch_metric_df,
      budget = if (units == "pixels") budget_px else budget,
      shortcut_threshold = 3,
      detour_cap = mobility_detour_cap,
      redundancy_method = redundancy_method
    )
    selected_ids <- opt$selected
    budget_used <- if (units == "pixels") {
      as.numeric(opt$stats$budget_used %||% 0)
    } else {
      as.numeric(opt$stats$budget_used %||% 0) * area_div / unit_conv$pixel_area
    }
    strategy_stats <- opt$stats
  } else {
    patch_optimizer_df <- data.frame(
      patch_id = patch_df$patch_id,
      area = patch_df$cell_count,
      x = patch_df$x / unit_conv$pixel_size,
      y = patch_df$y / unit_conv$pixel_size,
      stringsAsFactors = FALSE
    )
    opt <- terralink_optimize_connected_area_plugin_parity_df(
      candidates = candidates,
      patch_df = patch_optimizer_df,
      budget = budget_px,
      cost_col = "cost",
      length_col = "length",
      scale = 1,
      merge_equiv_area = 32,
      merge_equiv_ratio = 0.01
    )
    selected_ids <- opt$selected
    budget_used <- opt$stats$budget_used %||% 0
    strategy_stats <- opt$stats
  }

  selected_ids <- as.integer(unique(selected_ids))
  remaining_budget <- as.numeric(max(0, budget_px - budget_used))
  # Preserve the strategy-selected corridor set for the parity-aligned selectors.
  pure_metric_modes <- c("largest_single_network", "most_connected_habitat", "reachable_habitat_advanced", "landscape_fluidity")
  if (!(strategy_key %in% pure_metric_modes) && remaining_budget > 0 && nrow(candidates) > length(selected_ids)) {
    unselected <- candidates[!(candidates$id %in% selected_ids), , drop = FALSE]
    min_extra_cost <- suppressWarnings(min(unselected$cost, na.rm = TRUE))
    if (is.finite(min_extra_cost) && remaining_budget >= min_extra_cost) {
      hybrid <- terralink_apply_hybrid_leftover_budget_raster(
        candidates = candidates,
        selected_ids = selected_ids,
        patch_sizes = nodes,
        remaining_budget = remaining_budget
      )
      if (hybrid$budget_used > 0) {
        selected_ids <- as.integer(unique(hybrid$selected_ids))
        budget_used <- as.numeric(budget_used) + as.numeric(hybrid$budget_used)
        remaining_budget <- as.numeric(max(0, budget_px - budget_used))
      }
    }
  }

  selected <- candidates[candidates$id %in% selected_ids, , drop = FALSE]
  if (nrow(selected) > 0) {
    selected <- selected[order(match(selected$id, selected_ids)), , drop = FALSE]
  }
  if (!(strategy_key %in% pure_metric_modes) && remaining_budget > 0 && nrow(selected) > 0) {
    thickened <- terralink_thicken_corridors_raster(
      selected = selected,
      remaining_budget = remaining_budget,
      min_corridor_width_px = min_corridor_width_px,
      patch_sizes = nodes,
      labels = labels,
      max_width_factor = 3L
    )
    if (thickened$budget_used > 0) {
      selected <- thickened$selected
      budget_used <- as.numeric(budget_used) + as.numeric(thickened$budget_used)
    }
  }
  if (nrow(selected) > 0) {
    endpoint_ok <- vapply(seq_len(nrow(selected)), function(i) {
      terralink_raster_endpoint_area_guard_ok(
        selected$patch1[[i]],
        selected$patch2[[i]],
        selected$cost[[i]],
        nodes
      )
    }, logical(1))
    if (!all(endpoint_ok)) {
      selected <- selected[endpoint_ok, , drop = FALSE]
      budget_used <- sum(as.numeric(selected$cost), na.rm = TRUE)
    }
  }
  if (nrow(selected) > 0) {
    selected$connected_size <- terralink_compute_raster_connected_sizes(selected, nodes)
  }

  terralink_inform(sprintf("Corridors selected: %s", nrow(selected)), ctx = ctx, level = 1)
  terralink_inform(sprintf("Budget used (px): %s", terralink_format_number(budget_used, 2)), ctx = ctx, level = 1)

  terralink_progress_update(ctx, 80, "Building outputs")
  pre_mask <- terra::ifel(!is.na(labels) & labels > 0, 1, NA)
  corridor_raster <- build_corridor_raster(
    labels,
    patch_df,
    selected,
    min_corridor_width_px = min_corridor_width_px,
    assignment_mode = corridor_cell_assignment_key
  )
  contiguous_raster <- build_contiguous_raster(pre_mask, corridor_raster, patch_connectivity)

  terralink_progress_update(ctx, 90, "Calculating metrics")
  metric_context <- terralink_metric_context(
    patch_df = patch_metric_df,
    selected_corridors = selected,
    label = terralink_safe_name("TerraLink Raster"),
    area_unit = if (units == "imperial") "ac" else if (units == "metric") "ha" else "px^2",
    distance_unit = if (units == "imperial") "ft" else if (units == "metric") "m" else "px",
    budget_used = if (units == "pixels") budget_used else (budget_used * unit_conv$pixel_area / area_div),
    species_dispersal_distance = species_dispersal_distance %||% max_search_distance,
    species_dispersal_kernel = species_dispersal_kernel,
    min_patch_area_for_species = min_patch_area_for_species,
    patch_area_scaling = patch_area_scaling,
    max_search_distance = max_search_distance,
    mobility_detour_cap = mobility_detour_cap,
    redundancy_method = redundancy_method,
    metric_weights = metric_weights,
    weight_m = weight_m,
    weight_lcc = weight_lcc,
    weight_pc = weight_pc,
    weight_f = weight_f
  )
  metrics_report <- metric_context$report
  metrics <- metric_context$metrics
  terralink_inform("Metrics calculated.", ctx = ctx, level = 2)

  summary <- list(
    budget_total = budget,
    budget_used = budget_used,
    corridors_used = nrow(selected),
    candidate_edges = nrow(candidates),
    candidate_pairs = candidate_pairs_count,
    possible_pairs = possible_pairs,
    patches = nrow(patch_df),
    raw_patches = raw_patch_count,
    filtered_out = filtered_out,
    patch_values = patch_values,
    patch_ranges = patch_ranges,
    obstacle_values = obstacle_values,
    obstacle_ranges = obstacle_ranges,
    patch_value_cells = patch_cells,
    obstacle_cells = obstacle_cells,
    min_patch_size = min_patch_size,
    min_corridor_width = min_corridor_width,
    corridor_cell_assignment = corridor_cell_assignment_key,
    max_search_distance = max_search_distance,
    patch_connectivity = patch_connectivity,
    strategy = strategy_key,
    units = units
  )

  result <- list(
    patches = labels,
    patch_table = patch_df,
    corridors = selected,
    corridor_raster = corridor_raster,
    contiguous_raster = contiguous_raster,
    strategy = strategy_key,
    summary = summary,
    metrics = metrics,
    metrics_report = metrics_report,
    strategy_stats = strategy_stats
  )
  if (keep_candidates) result$candidates <- candidates
  result <- terralink_as_result(
    result,
    mode = "raster",
    inputs = list(
      units = units,
      pixel_size = unit_conv$pixel_size,
      pixel_area = unit_conv$pixel_area,
      raster_cells = preflight$n_cells,
      budget = budget,
      min_patch_size = min_patch_size,
      min_corridor_width = min_corridor_width,
      max_search_distance = max_search_distance
    ),
    run_stats = list(
      elapsed_s = proc.time()[[3]] - t_start,
      candidate_edges = nrow(candidates),
      candidate_pairs = candidate_pairs_count,
      possible_pairs = possible_pairs
    ),
    warnings = ctx$warnings,
    diagnostics = ctx$diagnostics
  )
  terralink_progress_done(ctx)
  result
}
