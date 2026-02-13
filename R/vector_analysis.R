# Vector workflow for TerraLink

terralink_merge_intersecting_patches <- function(patches_sf) {
  if (!inherits(patches_sf, "sf")) {
    terralink_abort("patches_sf must be an sf object.", class = "terralink_error_input")
  }
  if (nrow(patches_sf) == 0) return(patches_sf)

  geom <- sf::st_geometry(patches_sf)
  n <- length(geom)
  uf <- UnionFind$new()
  for (i in seq_len(n)) uf$find(i)

  ints <- sf::st_intersects(geom, geom, sparse = TRUE)
  for (i in seq_len(n)) {
    nei <- ints[[i]]
    if (length(nei) == 0) next
    for (j in nei) {
      if (j <= i) next
      uf$union(i, j)
    }
  }

  roots <- vapply(seq_len(n), function(i) as.integer(uf$find(i)), integer(1))
  groups <- split(seq_len(n), roots)
  merged <- list()
  k <- 1L
  crs_obj <- sf::st_crs(patches_sf)
  for (idx in groups) {
    if (length(idx) == 0) next
    g <- tryCatch(sf::st_union(geom[idx]), error = function(e) sf::st_combine(geom[idx]))
    g <- tryCatch(sf::st_make_valid(g), error = function(e) g)
    g_types <- tryCatch(unique(as.character(sf::st_geometry_type(g))), error = function(e) character(0))
    if (any(grepl("GEOMETRYCOLLECTION", g_types, fixed = TRUE))) {
      g <- tryCatch(sf::st_collection_extract(g, "POLYGON"), error = function(e) g)
    }
    if (length(g) == 0) next
    g <- tryCatch(sf::st_union(g), error = function(e) g[[1]])
    if (inherits(g, "sfc")) {
      merged[[k]] <- g[[1]]
    } else {
      merged[[k]] <- g
    }
    k <- k + 1L
  }
  if (length(merged) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_obj))[0, ])
  }
  sf::st_sf(geometry = sf::st_sfc(merged, crs = crs_obj))
}

terralink_geom_parts <- function(geom) {
  if (is.null(geom) || length(geom) == 0) return(NULL)
  geom <- tryCatch(sf::st_make_valid(geom), error = function(e) geom)
  g_types <- tryCatch(unique(as.character(sf::st_geometry_type(geom))), error = function(e) character(0))
  parts <- geom
  if (any(grepl("GEOMETRYCOLLECTION", g_types, fixed = TRUE))) {
    parts <- tryCatch(sf::st_collection_extract(geom, "POLYGON"), error = function(e) geom)
  }
  if (inherits(parts, "sfc") && length(parts) > 0) {
    poly <- tryCatch(sf::st_cast(parts, "POLYGON"), error = function(e) parts)
    if (inherits(poly, "sfc") && length(poly) > 0) return(poly)
    return(parts)
  }
  if (inherits(parts, "sfg")) {
    return(sf::st_sfc(parts, crs = sf::st_crs(geom)))
  }
  NULL
}

terralink_geom_overlap_ratio <- function(g1, g2, proximity_dist = 0) {
  if (is.null(g1) || is.null(g2)) return(1)
  if (length(g1) == 0 || length(g2) == 0) return(1)
  if (isTRUE(sf::st_is_empty(g1)) || isTRUE(sf::st_is_empty(g2))) return(1)
  if (is.finite(proximity_dist) && proximity_dist > 0) {
    d <- suppressWarnings(tryCatch(as.numeric(sf::st_distance(g1, g2, by_element = TRUE)), error = function(e) Inf))
    if (is.finite(d) && d <= proximity_dist) return(1)
  }
  a1 <- suppressWarnings(tryCatch(as.numeric(sf::st_area(g1)), error = function(e) NA_real_))
  a2 <- suppressWarnings(tryCatch(as.numeric(sf::st_area(g2)), error = function(e) NA_real_))
  denom <- min(a1, a2, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(1)
  inter <- tryCatch(sf::st_intersection(g1, g2), error = function(e) NULL)
  if (is.null(inter) || length(inter) == 0 || isTRUE(sf::st_is_empty(inter))) return(0)
  ia <- suppressWarnings(tryCatch(as.numeric(sf::st_area(inter)), error = function(e) 0))
  if (!is.finite(ia)) ia <- 0
  ia / denom
}

terralink_pick_part_closest_to_patch <- function(geom, patch_geom) {
  parts <- terralink_geom_parts(geom)
  if (is.null(parts) || length(parts) == 0) return(NULL)
  d <- suppressWarnings(tryCatch(as.numeric(sf::st_distance(parts, patch_geom)), error = function(e) rep(Inf, length(parts))))
  if (length(d) == 0) return(parts[1])
  parts[which.min(d)]
}

terralink_finalize_vector_corridor <- function(pid1, pid2, corridor_geom, patch_geom, patch_union = NULL) {
  if (is.null(corridor_geom) || length(corridor_geom) == 0 || isTRUE(sf::st_is_empty(corridor_geom))) {
    return(list(geom = NULL, patch_ids = integer(0), extra_patch_ids = integer(0)))
  }
  touched <- tryCatch(sf::st_intersects(corridor_geom, patch_geom, sparse = TRUE)[[1]], error = function(e) integer(0))
  patch_ids <- sort(unique(c(as.integer(pid1), as.integer(pid2), as.integer(touched))))
  patch_ids <- patch_ids[is.finite(patch_ids)]
  final_geom <- corridor_geom
  if (!is.null(patch_union) && length(patch_union) > 0 && !isTRUE(sf::st_is_empty(patch_union))) {
    final_geom <- tryCatch(sf::st_difference(final_geom, patch_union), error = function(e) final_geom)
  } else if (length(patch_ids) > 0) {
    for (pid in patch_ids) {
      final_geom <- tryCatch(sf::st_difference(final_geom, patch_geom[pid]), error = function(e) final_geom)
      if (isTRUE(sf::st_is_empty(final_geom))) break
    }
  }
  final_geom <- tryCatch(sf::st_make_valid(final_geom), error = function(e) final_geom)
  if (is.null(final_geom) || length(final_geom) == 0 || isTRUE(sf::st_is_empty(final_geom))) {
    return(list(geom = NULL, patch_ids = integer(0), extra_patch_ids = integer(0)))
  }
  extra <- setdiff(patch_ids, c(as.integer(pid1), as.integer(pid2)))
  list(geom = final_geom, patch_ids = patch_ids, extra_patch_ids = extra)
}

terralink_push_vector_candidate <- function(candidates_by_pair, candidate, max_keep_per_pair = 8L, min_distinct_overlap_ratio = 0.75, proximity_dist = 0) {
  p1 <- as.integer(candidate$patch1)
  p2 <- as.integer(candidate$patch2)
  if (!is.finite(p1) || !is.finite(p2) || p1 == p2) return(candidates_by_pair)
  if (is.null(candidate$corridor) || isTRUE(sf::st_is_empty(candidate$corridor))) return(candidates_by_pair)
  key <- paste(sort(c(p1, p2)), collapse = "_")
  existing <- candidates_by_pair[[key]]
  if (is.null(existing)) existing <- list()
  for (prev in existing) {
    ov <- suppressWarnings(tryCatch(
      terralink_geom_overlap_ratio(candidate$corridor, prev$corridor, proximity_dist = proximity_dist),
      error = function(e) 1
    ))
    if (is.finite(ov) && ov >= min_distinct_overlap_ratio) {
      return(candidates_by_pair)
    }
  }
  existing[[length(existing) + 1L]] <- candidate
  ord <- order(vapply(existing, function(x) as.numeric(x$area), numeric(1)), vapply(existing, function(x) as.numeric(x$length), numeric(1)))
  existing <- existing[ord]
  if (length(existing) > as.integer(max_keep_per_pair)) {
    existing <- existing[seq_len(as.integer(max_keep_per_pair))]
  }
  candidates_by_pair[[key]] <- existing
  candidates_by_pair
}

build_vector_cost_surface <- function(patches, obstacles, resolution) {
  ext <- terra::ext(sf::st_bbox(patches))
  r <- terra::rast(ext, res = resolution, crs = sf::st_crs(patches)$wkt)
  terra::values(r) <- 1
  if (!is.null(obstacles) && nrow(obstacles) > 0) {
    obs_vect <- terra::vect(obstacles)
    r <- terra::rasterize(obs_vect, r, field = 0, background = 1)
    r[r == 0] <- NA
  }
  r
}

terralink_as_linestring_sfc <- function(path, crs) {
  if (is.null(path)) return(NULL)
  line <- NULL
  if (inherits(path, "sfc")) {
    line <- path
  } else if (inherits(path, "sf")) {
    line <- sf::st_geometry(path)
  } else if (inherits(path, "SpatialLines") || inherits(path, "SpatialLinesDataFrame")) {
    line <- tryCatch(sf::st_geometry(sf::st_as_sf(path)), error = function(e) NULL)
  } else {
    line <- tryCatch(sf::st_as_sfc(path), error = function(e) NULL)
  }
  if (is.null(line) || length(line) == 0) return(NULL)
  line <- tryCatch(sf::st_cast(line, "LINESTRING"), error = function(e) line)
  if (is.null(line) || length(line) == 0 || isTRUE(all(sf::st_is_empty(line)))) return(NULL)
  if (length(line) > 1) {
    line <- tryCatch(sf::st_line_merge(sf::st_union(line)), error = function(e) line)
    if (inherits(line, "sfc") && length(line) > 1) line <- line[1]
  }
  sf::st_set_crs(line, sf::st_crs(crs))
}

build_vector_candidates <- function(
  patches,
  patch_df,
  max_search_distance_m,
  width_m,
  obstacles = NULL,
  obstacle_resolution = NULL,
  pair_index = NULL,
  strategy_key = "circuit_utility"
) {
  if (nrow(patch_df) < 2) return(data.frame())
  patch_geom <- sf::st_geometry(patches)
  patch_union <- tryCatch(sf::st_union(patch_geom), error = function(e) NULL)
  strategy_key <- match.arg(tolower(strategy_key), c("circuit_utility", "largest_network"))
  largest_network_mode <- identical(strategy_key, "largest_network")
  max_keep_per_pair <- 8L
  min_distinct_overlap_ratio <- 0.75
  proximity_dist <- max(as.numeric(width_m) * 1.5, 0)

  use_gdistance <- requireNamespace("gdistance", quietly = TRUE) &&
    requireNamespace("raster", quietly = TRUE) &&
    requireNamespace("sp", quietly = TRUE)
  if (use_gdistance && !is.null(obstacles)) {
    if (is.null(obstacle_resolution) || !is.finite(obstacle_resolution)) {
      obstacle_resolution <- max(width_m / 2, max_search_distance_m / 200)
    }
  }

  route_obstacles <- obstacles
  if (!is.null(route_obstacles) && nrow(route_obstacles) > 0 && is.finite(width_m) && width_m > 0) {
    route_obstacles <- tryCatch(sf::st_buffer(route_obstacles, width_m / 2), error = function(e) route_obstacles)
  }
  obstacle_union <- NULL
  if (!is.null(obstacles) && nrow(obstacles) > 0) {
    obstacle_union <- tryCatch(sf::st_union(sf::st_geometry(obstacles)), error = function(e) NULL)
    obstacle_union <- tryCatch(sf::st_make_valid(obstacle_union), error = function(e) obstacle_union)
  }

  tr <- NULL
  if (use_gdistance && !is.null(route_obstacles)) {
    cost_surface <- build_vector_cost_surface(patches, route_obstacles, obstacle_resolution)
    cost_r <- tryCatch(raster::raster(cost_surface), error = function(e) NULL)
    if (!is.null(cost_r)) {
      tr <- tryCatch(gdistance::transition(cost_r, function(x) 1 / mean(x), directions = 8), error = function(e) NULL)
      if (!is.null(tr)) {
        tr <- tryCatch(gdistance::geoCorrection(tr, type = "c"), error = function(e) tr)
      }
    }
  }

  candidates_by_pair <- list()

  build_line_for_pair <- function(i, j) {
    line <- NULL
    length_m <- suppressWarnings(tryCatch(
      as.numeric(sf::st_distance(sf::st_geometry(patches)[i], sf::st_geometry(patches)[j], by_element = TRUE)),
      error = function(e) Inf
    ))
    if (!is.finite(length_m)) return(NULL)
    if (!is.null(tr)) {
      p1 <- matrix(c(patch_df$x[i], patch_df$y[i]), ncol = 2)
      p2 <- matrix(c(patch_df$x[j], patch_df$y[j]), ncol = 2)
      path <- tryCatch(gdistance::shortestPath(tr, p1, p2, output = "SpatialLines"), error = function(e) NULL)
      if (!is.null(path)) {
        line <- terralink_as_linestring_sfc(path, sf::st_crs(patches))
        length_m <- as.numeric(sf::st_length(line))
      }
    }
    if (is.null(line)) {
      geom_i <- sf::st_geometry(patches)[i]
      geom_j <- sf::st_geometry(patches)[j]
      line <- tryCatch({
        geom <- sf::st_nearest_points(geom_i, geom_j)
        sf::st_cast(geom, "LINESTRING")
      }, error = function(e) NULL)
      if (is.null(line)) {
        line <- sf::st_sfc(
          sf::st_linestring(matrix(c(patch_df$x[i], patch_df$y[i], patch_df$x[j], patch_df$y[j]), ncol = 2, byrow = TRUE)),
          crs = sf::st_crs(patches)
        )
      }
      length_m <- as.numeric(sf::st_length(line))
    }
    list(line = line, length_m = as.numeric(length_m))
  }

  emit_candidate <- function(p1, p2, line, corridor_geom, length_m, patch_ids = c(p1, p2), distance_m = NULL) {
    if (is.null(corridor_geom) || length(corridor_geom) == 0 || isTRUE(sf::st_is_empty(corridor_geom))) return(invisible(NULL))
    area <- suppressWarnings(tryCatch(as.numeric(sf::st_area(corridor_geom)), error = function(e) NA_real_))
    if (!is.finite(area) || area <= 0) return(invisible(NULL))
    p1_idx <- match(p1, patch_df$patch_id)
    p2_idx <- match(p2, patch_df$patch_id)
    if (!is.finite(p1_idx) || !is.finite(p2_idx)) return(invisible(NULL))
    cand <- list(
      patch1 = as.integer(p1),
      patch2 = as.integer(p2),
      patch_ids = as.integer(sort(unique(patch_ids))),
      cost = as.numeric(area),
      length = as.numeric(length_m),
      area = as.numeric(area),
      roi = sqrt(patch_df$area_m2[p1_idx] * patch_df$area_m2[p2_idx]) / max(as.numeric(area), 1e-6),
      line = line,
      corridor = corridor_geom,
      distance_m = if (is.null(distance_m)) as.numeric(length_m) else as.numeric(distance_m)
    )
    candidates_by_pair <<- terralink_push_vector_candidate(
      candidates_by_pair = candidates_by_pair,
      candidate = cand,
      max_keep_per_pair = max_keep_per_pair,
      min_distinct_overlap_ratio = min_distinct_overlap_ratio,
      proximity_dist = proximity_dist
    )
    invisible(NULL)
  }

  process_pair <- function(i, j) {
    dist_ij <- suppressWarnings(tryCatch(
      as.numeric(sf::st_distance(sf::st_geometry(patches)[i], sf::st_geometry(patches)[j], by_element = TRUE)),
      error = function(e) Inf
    ))
    if (!is.finite(dist_ij) || dist_ij <= 0 || dist_ij > max_search_distance_m) return(invisible(NULL))
    routed <- build_line_for_pair(i, j)
    if (is.null(routed)) return(invisible(NULL))
    corridor_geom_raw <- tryCatch(sf::st_buffer(routed$line, dist = width_m / 2, nQuadSegs = 16), error = function(e) NULL)
    if (is.null(corridor_geom_raw) || length(corridor_geom_raw) == 0 || isTRUE(sf::st_is_empty(corridor_geom_raw))) return(invisible(NULL))
    if (!is.null(obstacle_union) && length(obstacle_union) > 0 && !isTRUE(sf::st_is_empty(obstacle_union))) {
      corridor_geom_raw <- tryCatch(sf::st_difference(corridor_geom_raw, obstacle_union), error = function(e) corridor_geom_raw)
      corridor_geom_raw <- tryCatch(sf::st_make_valid(corridor_geom_raw), error = function(e) corridor_geom_raw)
      if (is.null(corridor_geom_raw) || length(corridor_geom_raw) == 0 || isTRUE(sf::st_is_empty(corridor_geom_raw))) return(invisible(NULL))
    }

    pid1 <- patch_df$patch_id[i]
    pid2 <- patch_df$patch_id[j]
    finalized <- terralink_finalize_vector_corridor(
      pid1 = pid1,
      pid2 = pid2,
      corridor_geom = corridor_geom_raw,
      patch_geom = patch_geom,
      patch_union = patch_union
    )
    if (is.null(finalized$geom)) return(invisible(NULL))

    extra_patches <- as.integer(finalized$extra_patch_ids)
    emitted_any <- FALSE
    corridor_parts <- terralink_geom_parts(finalized$geom)
    if (length(extra_patches) > 0 && !is.null(corridor_parts) && length(corridor_parts) > 1) {
      part_a <- terralink_pick_part_closest_to_patch(finalized$geom, patch_geom[i])
      part_b <- terralink_pick_part_closest_to_patch(finalized$geom, patch_geom[j])
      if (!is.null(part_a) && !isTRUE(sf::st_is_empty(part_a))) {
        tgt_a <- extra_patches[which.min(suppressWarnings(tryCatch(
          as.numeric(sf::st_distance(patch_geom[i], patch_geom[extra_patches])),
          error = function(e) rep(Inf, length(extra_patches))
        )))]
        if (is.finite(tgt_a)) {
          dist_a <- suppressWarnings(tryCatch(as.numeric(sf::st_distance(patch_geom[i], patch_geom[tgt_a], by_element = TRUE)), error = function(e) dist_ij))
          emit_candidate(pid1, patch_df$patch_id[tgt_a], routed$line, part_a, dist_a, patch_ids = c(pid1, patch_df$patch_id[tgt_a]), distance_m = dist_a)
          emitted_any <- TRUE
        }
      }
      if (!is.null(part_b) && !isTRUE(sf::st_is_empty(part_b))) {
        tgt_b <- extra_patches[which.min(suppressWarnings(tryCatch(
          as.numeric(sf::st_distance(patch_geom[j], patch_geom[extra_patches])),
          error = function(e) rep(Inf, length(extra_patches))
        )))]
        if (is.finite(tgt_b)) {
          dist_b <- suppressWarnings(tryCatch(as.numeric(sf::st_distance(patch_geom[j], patch_geom[tgt_b], by_element = TRUE)), error = function(e) dist_ij))
          emit_candidate(pid2, patch_df$patch_id[tgt_b], routed$line, part_b, dist_b, patch_ids = c(pid2, patch_df$patch_id[tgt_b]), distance_m = dist_b)
          emitted_any <- TRUE
        }
      }
    }

    if (isTRUE(emitted_any)) return(invisible(NULL))
    if (isTRUE(largest_network_mode) && length(extra_patches) > 0) return(invisible(NULL))
    emit_candidate(pid1, pid2, routed$line, finalized$geom, routed$length_m, patch_ids = finalized$patch_ids, distance_m = routed$length_m)
    invisible(NULL)
  }

  if (is.null(pair_index)) {
    for (i in seq_len(nrow(patch_df) - 1)) {
      for (j in (i + 1):nrow(patch_df)) {
        process_pair(i, j)
      }
    }
  } else {
    if (nrow(pair_index) == 0) return(data.frame())
    for (row in seq_len(nrow(pair_index))) {
      i <- pair_index[row, 1]
      j <- pair_index[row, 2]
      process_pair(i, j)
    }
  }

  edges <- unlist(candidates_by_pair, recursive = FALSE)
  if (length(edges) == 0) return(data.frame())
  out <- do.call(rbind, lapply(seq_along(edges), function(i) {
    e <- edges[[i]]
    data.frame(
      patch1 = as.integer(e$patch1),
      patch2 = as.integer(e$patch2),
      cost = as.numeric(e$cost),
      length = as.numeric(e$length),
      area = as.numeric(e$area),
      id = as.integer(i),
      roi = as.numeric(e$roi),
      distance_m = as.numeric(e$distance_m),
      stringsAsFactors = FALSE
    )
  }))
  out$patch_ids <- I(lapply(edges, function(e) as.integer(e$patch_ids)))
  out$line <- sf::st_sfc(lapply(edges, function(e) e$line[[1]]), crs = sf::st_crs(patches))
  out$corridor <- sf::st_sfc(lapply(edges, function(e) e$corridor[[1]]), crs = sf::st_crs(patches))
  rownames(out) <- NULL
  out
}

#' Run TerraLink vector workflow
#'
#' @param patches sf polygons (one feature per patch) or file path.
#' @param budget Corridor budget (ha/ac).
#' @param strategy Strategy name: "circuit_utility", "most_connectivity" (alias), or "largest_network".
#' @param min_patch_size Minimum patch size (ha/ac).
#' @param min_corridor_width Minimum corridor width (m/ft).
#' @param max_search_distance Maximum search distance (m/ft).
#' @param obstacle_layers Optional obstacle layers (sf or file paths).
#' @param obstacle_resolution Raster resolution for obstacle routing.
#' @param units "metric" or "imperial".
#' @param max_pair_checks Limit for candidate pair checks.
#' @param max_candidates Limit for candidate corridors.
#' @param verbose Verbosity level (0-2).
#' @param progress Show progress bars.
#' @param obstacle_strategy Behavior when gdistance is unavailable and obstacles are provided.
#' @param return_crs CRS for outputs ("input" or "utm").
#' @param keep_candidates Keep candidate list in output.
#' @return List with corridors, networks, and summary.
#' @export
run_vector_analysis <- function(
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
  keep_candidates = FALSE
) {
  strategy_key <- match.arg(tolower(strategy), c("most_connectivity", "largest_network", "circuit_utility"))
  if (strategy_key == "most_connectivity") strategy_key <- "circuit_utility"

  ctx <- terralink_new_run_context(verbose = verbose, progress = progress)
  terralink_progress_start(ctx, message = "Starting vector analysis")
  t_start <- proc.time()[[3]]

  patches_sf <- terralink_resolve_vector(patches, quiet = TRUE)
  patches_sf <- terralink_preflight_vector(patches_sf, ctx = ctx)
  input_crs <- sf::st_crs(patches_sf)

  if (!terralink_is_projected(input_crs)) {
    utm <- terralink_pick_utm_crs(patches_sf)
    patches_sf <- sf::st_transform(patches_sf, utm)
    work_crs <- utm
  } else {
    work_crs <- input_crs
  }

  patches_sf <- terralink_merge_intersecting_patches(patches_sf)
  raw_patch_count <- nrow(patches_sf)

  return_crs <- match.arg(return_crs)
  output_crs <- if (return_crs == "input") input_crs else work_crs

  if (!is.numeric(budget) || budget <= 0) {
    terralink_abort("budget must be a positive number.", class = "terralink_error_input")
  }
  if (!is.numeric(max_search_distance) || max_search_distance <= 0) {
    terralink_abort("max_search_distance must be a positive number.", class = "terralink_error_input")
  }
  if (!is.numeric(min_corridor_width) || min_corridor_width <= 0) {
    terralink_abort("min_corridor_width must be a positive number.", class = "terralink_error_input")
  }

  obstacle_strategy <- match.arg(obstacle_strategy)
  has_obstacles <- !is.null(obstacle_layers)
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
      obstacle_layers <- NULL
    } else {
      terralink_warn("gdistance not available; using straight-line routing without obstacle awareness.", ctx = ctx)
    }
  }

  terralink_progress_update(ctx, 10, "Converting units")
  units <- match.arg(units, c("metric", "imperial"))
  area_factor <- if (units == "imperial") 4046.8564224 else 10000.0
  dist_factor <- if (units == "imperial") 0.3048 else 1.0

  budget_m2 <- as.numeric(budget) * area_factor
  min_patch_m2 <- if (is.null(min_patch_size)) 0 else as.numeric(min_patch_size) * area_factor
  width_m <- as.numeric(min_corridor_width) * dist_factor
  max_search_m <- as.numeric(max_search_distance) * dist_factor

  patch_area_m2 <- as.numeric(sf::st_area(patches_sf))
  keep <- patch_area_m2 >= min_patch_m2
  patches_sf <- patches_sf[keep, , drop = FALSE]
  patch_area_m2 <- patch_area_m2[keep]
  filtered_out <- raw_patch_count - nrow(patches_sf)
  if (nrow(patches_sf) < 2) {
    terralink_abort("Vector workflow requires at least two patch features.", class = "terralink_error_input")
  }

  centroids <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(patches_sf)))
  patch_df <- data.frame(
    patch_id = seq_len(nrow(patches_sf)),
    area_m2 = patch_area_m2,
    x = centroids[, 1],
    y = centroids[, 2],
    stringsAsFactors = FALSE
  )
  possible_pairs <- if (nrow(patch_df) > 1) as.integer((nrow(patch_df) * (nrow(patch_df) - 1)) / 2) else 0L
  terralink_check_candidate_count(
    possible_pairs,
    max_candidates = max_pair_checks,
    ctx = ctx,
    scope = "Candidate pair checks",
    override_param = "max_pair_checks"
  )
  terralink_inform(
    sprintf("Patches labeled: %s (raw %s, filtered %s)", nrow(patch_df), raw_patch_count, filtered_out),
    ctx = ctx,
    level = 1
  )

  obstacles <- NULL
  if (!is.null(obstacle_layers)) {
    if (inherits(obstacle_layers, "list")) {
      obs_list <- lapply(obstacle_layers, function(x) terralink_resolve_vector(x, quiet = TRUE))
      obstacles <- do.call(rbind, obs_list)
    } else {
      obstacles <- terralink_resolve_vector(obstacle_layers, quiet = TRUE)
    }
    obstacles <- sf::st_transform(obstacles, sf::st_crs(patches_sf))
    invalid_obs <- sf::st_is_valid(obstacles)
    if (any(!invalid_obs, na.rm = TRUE)) {
      fixed_obs <- terralink_make_valid(obstacles)
      if (!is.null(fixed_obs)) {
        terralink_warn(sprintf("Fixing %s invalid obstacle geometries.", sum(!invalid_obs, na.rm = TRUE)), ctx = ctx)
        obstacles <- fixed_obs
      } else {
        terralink_warn("Invalid obstacle geometries detected; results may be unreliable. Install sf >= 1.0 or lwgeom to auto-fix.", ctx = ctx)
      }
    }
  }

  terralink_progress_update(ctx, 55, "Generating corridor candidates")
  pair_index <- terralink_vector_pair_index(
    patches_sf = patches_sf,
    max_distance = max_search_m,
    x = patch_df$x,
    y = patch_df$y
  )

  candidates <- build_vector_candidates(
    patches = patches_sf,
    patch_df = patch_df,
    max_search_distance_m = max_search_m,
    width_m = width_m,
    obstacles = obstacles,
    obstacle_resolution = obstacle_resolution,
    pair_index = pair_index,
    strategy_key = strategy_key
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
      patches = nrow(patches_sf),
      raw_patches = raw_patch_count,
      filtered_out = filtered_out,
      units = units
    )
    patches_out <- patches_sf
    if (!is.null(output_crs) && !is.na(output_crs) && !isTRUE(sf::st_crs(patches_out) == output_crs)) {
      patches_out <- sf::st_transform(patches_out, output_crs)
    }
    result <- list(
      patches = patches_out,
      corridors = candidates,
      networks = NULL,
      summary = summary
    )
    result <- terralink_as_result(
      result,
      mode = "vector",
      inputs = list(units = units, budget = budget, min_patch_size = min_patch_size, min_corridor_width = min_corridor_width, max_search_distance = max_search_distance),
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

  nodes <- stats::setNames(patch_df$area_m2, patch_df$patch_id)
  engine_edges <- data.frame(
    u = candidates$patch1,
    v = candidates$patch2,
    id = candidates$id,
    cost = candidates$cost
  )

  if (strategy_key == "largest_network") {
    opt <- optimize_largest_network(nodes, engine_edges, budget = budget_m2)
    selected_ids <- opt$selected
    budget_used <- opt$total_cost
    primary_links <- NA_integer_
    redundant_links <- NA_integer_
  } else {
    opt <- optimize_strategy(
      strategy = "circuit_utility",
      nodes = nodes,
      edges = engine_edges,
      candidates = candidates,
      budget = budget_m2,
      get_patch_ids = function(cand) {
        if ("patch_ids" %in% names(cand)) {
          ids <- cand$patch_ids
          if (is.list(ids) && length(ids) > 0) ids <- ids[[1]]
          ids <- suppressWarnings(as.integer(ids))
          ids <- ids[is.finite(ids)]
          if (length(ids) >= 2) return(ids)
        }
        suppressWarnings(as.integer(c(cand$patch1, cand$patch2)))
      },
      get_pair_key = function(cand) sort(c(cand$patch1, cand$patch2)),
      get_cost = function(cand) cand$cost,
      get_base_roi = function(cand) cand$roi,
      get_length = function(cand) cand$length,
      get_patch_size = function(pid) nodes[[as.character(pid)]],
      overlap_ratio = function(cand, prior) {
        if (length(prior) == 0) return(0)
        ratios <- vapply(prior, function(p) {
          inter <- tryCatch(sf::st_intersection(cand$corridor, p), error = function(e) NULL)
          if (is.null(inter) || length(inter) == 0) return(0)
          as.numeric(sf::st_area(inter)) / max(as.numeric(sf::st_area(cand$corridor)), 1e-6)
        }, numeric(1))
        max(ratios, na.rm = TRUE)
      },
      overlap_obj = function(cand) cand$corridor,
      redundancy_distance_ok = function(cand, prior) {
        if (length(prior) == 0 || max_search_m <= 0) return(TRUE)
        for (p in prior) {
          d <- suppressWarnings(tryCatch(as.numeric(sf::st_distance(cand$corridor, p, by_element = TRUE)), error = function(e) Inf))
          if (is.finite(d) && d < max_search_m) return(FALSE)
        }
        TRUE
      },
      max_links_per_pair = 1
    )
    selected_ids <- opt$selected
    budget_used <- opt$stats$budget_used %||% 0
    primary_links <- as.integer(opt$stats$primary_links %||% 0L)
    redundant_links <- as.integer(opt$stats$redundant_links %||% 0L)
  }

  selected <- candidates[candidates$id %in% selected_ids, , drop = FALSE]
  if (nrow(selected) > 0) {
    pair_keys_selected <- paste(pmin(selected$patch1, selected$patch2), pmax(selected$patch1, selected$patch2), sep = "_")
    keep_unique_pair <- !duplicated(pair_keys_selected)
    if (any(!keep_unique_pair)) {
      selected <- selected[keep_unique_pair, , drop = FALSE]
      budget_used <- sum(as.numeric(selected$cost), na.rm = TRUE)
      if (!is.na(primary_links)) {
        primary_links <- as.integer(nrow(selected))
      }
      if (!is.na(redundant_links)) {
        redundant_links <- 0L
      }
    }
  }
  terralink_inform(sprintf("Corridors selected: %s", nrow(selected)), ctx = ctx, level = 1)
  if (nrow(selected) == 0) {
    summary <- list(
      budget_total = budget,
      budget_used = budget_used,
      corridors_used = 0,
      candidate_edges = nrow(candidates),
      patches = nrow(patches_sf),
      units = units
    )
    patches_out <- patches_sf
    if (!is.null(output_crs) && !is.na(output_crs) && !isTRUE(sf::st_crs(patches_out) == output_crs)) {
      patches_out <- sf::st_transform(patches_out, output_crs)
    }
    result <- list(
      patches = patches_out,
      corridors = selected,
      networks = NULL,
      summary = summary
    )
    result <- terralink_as_result(
      result,
      mode = "vector",
      inputs = list(units = units, budget = budget, min_patch_size = min_patch_size, min_corridor_width = min_corridor_width, max_search_distance = max_search_distance),
      run_stats = list(elapsed_s = proc.time()[[3]] - t_start, candidate_edges = nrow(candidates)),
      warnings = ctx$warnings,
      diagnostics = list(message = "No corridors selected; try increasing budget or search distance.")
    )
    terralink_progress_done(ctx)
    return(result)
  }

  terralink_progress_update(ctx, 80, "Building outputs")
  # Build corridors sf
  corridors_sf <- sf::st_sf(
    patch1 = selected$patch1,
    patch2 = selected$patch2,
    patch_ids = vapply(selected$patch_ids, function(ids) paste(sort(as.integer(ids)), collapse = ","), character(1)),
    corridor_area_m2 = selected$area,
    corridor_length_m = selected$length,
    geometry = selected$corridor,
    crs = sf::st_crs(patches_sf)
  )

  # Component sizes and efficiency (ensure all patches are included)
  patch_ids_all <- as.character(patch_df$patch_id)
  g <- igraph::make_empty_graph(n = length(patch_ids_all), directed = FALSE)
  igraph::V(g)$name <- patch_ids_all
  for (row in seq_len(nrow(selected))) {
    ids <- selected$patch_ids[[row]]
    ids <- suppressWarnings(as.integer(ids))
    ids <- ids[is.finite(ids)]
    if (length(ids) < 2) {
      ids <- suppressWarnings(as.integer(c(selected$patch1[[row]], selected$patch2[[row]])))
      ids <- ids[is.finite(ids)]
    }
    if (length(ids) < 2) next
    anchor <- as.character(ids[[1]])
    for (other in ids[-1]) {
      g <- igraph::add_edges(g, c(anchor, as.character(other)))
    }
  }
  comps <- igraph::components(g)
  comp_area <- tapply(patch_area_m2, comps$membership, sum)
  comp_names <- names(comps$membership)
  corridors_sf$component_id <- comps$membership[match(as.character(corridors_sf$patch1), comp_names)]
  corridors_sf$connected_area_m2 <- comp_area[as.character(corridors_sf$component_id)]
  corridors_sf$efficiency <- corridors_sf$connected_area_m2 / pmax(corridors_sf$corridor_area_m2, 1e-6)
  geom_types <- as.character(sf::st_geometry_type(corridors_sf))
  corridors_sf$multipart <- geom_types %in% c("MULTIPOLYGON", "MULTILINESTRING")

  # Build network polygons per component
  patches_sf$patch_id <- patch_df$patch_id
  patches_sf$component_id <- comps$membership[match(as.character(patches_sf$patch_id), comp_names)]
  net_polys <- list()
  for (comp_id in sort(unique(patches_sf$component_id))) {
    patch_part <- patches_sf[patches_sf$component_id == comp_id, ]
    corr_part <- corridors_sf[corridors_sf$component_id == comp_id, ]
    patch_geom <- sf::st_union(patch_part)
    if (nrow(corr_part) > 0) {
      corr_geom <- sf::st_union(corr_part)
      geom <- sf::st_union(patch_geom, corr_geom)
    } else {
      geom <- patch_geom
    }
    dissolve_tolerance <- max(0, as.numeric(width_m) * 0.01)
    if (is.finite(dissolve_tolerance) && dissolve_tolerance > 0) {
      geom <- tryCatch(sf::st_buffer(geom, dissolve_tolerance), error = function(e) geom)
      geom <- tryCatch(sf::st_buffer(geom, -dissolve_tolerance), error = function(e) geom)
    }
    geom <- tryCatch(sf::st_make_valid(geom), error = function(e) geom)
    net_polys[[length(net_polys) + 1]] <- sf::st_sf(component_id = comp_id, area_m2 = as.numeric(sf::st_area(geom)), geometry = geom)
  }
  networks_sf <- do.call(rbind, net_polys)

  if (is.null(networks_sf) || !inherits(networks_sf, "sf")) {
    networks_sf <- sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(patches_sf)))[0, ]
  }
  terralink_progress_update(ctx, 90, "Calculating metrics")
  metrics_report <- terralink_vector_report(patches_sf, networks_sf, units = units, label = terralink_safe_name("TerraLink Vector"))
  terralink_inform("Metrics calculated.", ctx = ctx, level = 2)

  # Convert to desired units for reporting
  area_div <- if (units == "imperial") 4046.8564224 else 10000.0
  dist_mult <- if (units == "imperial") 3.28084 else 1.0
  corridors_sf$corridor_area <- corridors_sf$corridor_area_m2 / area_div
  corridors_sf$connected_area <- corridors_sf$connected_area_m2 / area_div
  corridors_sf$corridor_length <- corridors_sf$corridor_length_m * dist_mult
  terralink_inform(
    sprintf("Budget used (%s): %s", if (units == "imperial") "ac" else "ha", terralink_format_number(budget_used / area_div, 2)),
    ctx = ctx,
    level = 1
  )

  summary <- list(
    budget_total = budget,
    budget_used = budget_used / area_div,
    corridors_used = nrow(corridors_sf),
    candidate_edges = nrow(candidates),
    candidate_pairs = candidate_pairs_count,
    possible_pairs = possible_pairs,
    patches = nrow(patches_sf),
    raw_patches = raw_patch_count,
    filtered_out = filtered_out,
    primary_links = primary_links,
    redundant_links = redundant_links,
    units = units
  )

  patches_out <- patches_sf
  corridors_out <- corridors_sf
  networks_out <- networks_sf
  if (!is.null(output_crs) && !is.na(output_crs) && !isTRUE(sf::st_crs(patches_out) == output_crs)) {
    patches_out <- sf::st_transform(patches_out, output_crs)
    corridors_out <- sf::st_transform(corridors_out, output_crs)
    networks_out <- sf::st_transform(networks_out, output_crs)
  }

  result <- list(
    patches = patches_out,
    corridors = corridors_out,
    networks = networks_out,
    summary = summary,
    metrics_report = metrics_report
  )
  if (keep_candidates) result$candidates <- candidates
  result <- terralink_as_result(
    result,
    mode = "vector",
    inputs = list(
      units = units,
      budget = budget,
      min_patch_size = min_patch_size,
      min_corridor_width = min_corridor_width,
      max_search_distance = max_search_distance,
      crs_input = input_crs,
      crs_work = work_crs,
      crs_output = output_crs
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
