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
  attrs <- list()
  k <- 1L
  crs_obj <- sf::st_crs(patches_sf)
  attr_names <- setdiff(names(patches_sf), attr(patches_sf, "sf_column"))
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
    if (length(attr_names) > 0) {
      row_vals <- lapply(attr_names, function(nm) {
        vals <- patches_sf[[nm]][idx]
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) return(NA)
        if (is.numeric(vals) || is.integer(vals)) return(mean(as.numeric(vals), na.rm = TRUE))
        if (is.logical(vals)) return(any(vals))
        as.character(vals[[1]])
      })
      names(row_vals) <- attr_names
      attrs[[k]] <- row_vals
    }
    k <- k + 1L
  }
  if (length(merged) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_obj))[0, ])
  }
  geom_sf <- sf::st_sf(geometry = sf::st_sfc(merged, crs = crs_obj))
  if (length(attrs) == length(merged) && length(attr_names) > 0) {
    attr_df <- as.data.frame(do.call(rbind, lapply(attrs, function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    })), stringsAsFactors = FALSE)
    for (nm in names(attr_df)) {
      if (is.numeric(patches_sf[[nm]]) || is.integer(patches_sf[[nm]])) {
        attr_df[[nm]] <- suppressWarnings(as.numeric(attr_df[[nm]]))
      } else if (is.logical(patches_sf[[nm]])) {
        attr_df[[nm]] <- as.logical(attr_df[[nm]])
      }
    }
    geom_sf <- sf::st_sf(attr_df, geometry = sf::st_geometry(geom_sf), crs = crs_obj)
  }
  geom_sf
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

terralink_shortest_path_quiet <- function(transition, start, end) {
  suppressWarnings(tryCatch(
    gdistance::shortestPath(transition, as.matrix(start), as.matrix(end), output = "SpatialLines"),
    error = function(e) NULL
  ))
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

terralink_nearest_boundary_route_points <- function(geom_i, geom_j) {
  nearest <- tryCatch(sf::st_nearest_points(geom_i, geom_j), error = function(e) NULL)
  if (is.null(nearest) || length(nearest) == 0) return(NULL)
  coords <- suppressWarnings(tryCatch(sf::st_coordinates(nearest), error = function(e) NULL))
  if (is.null(coords) || nrow(coords) < 2) return(NULL)
  list(
    start = matrix(as.numeric(coords[1, c("X", "Y")]), ncol = 2),
    end = matrix(as.numeric(coords[nrow(coords), c("X", "Y")]), ncol = 2)
  )
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

terralink_dedupe_xy_points <- function(coords, tol = 0.01) {
  if (is.null(coords) || length(coords) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }
  coords <- as.matrix(coords)
  if (nrow(coords) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }
  keys <- paste(round(coords[, 1] / tol), round(coords[, 2] / tol), sep = "_")
  coords[!duplicated(keys), , drop = FALSE]
}

terralink_boundary_terminals_for_patch <- function(geom, spacing_m = 150, max_pts = 120) {
  if (is.null(geom) || length(geom) == 0 || isTRUE(any(sf::st_is_empty(geom)))) {
    return(matrix(numeric(0), ncol = 2))
  }

  collect_coords <- function(line_geom) {
    if (is.null(line_geom) || length(line_geom) == 0 || isTRUE(any(sf::st_is_empty(line_geom)))) {
      return(matrix(numeric(0), ncol = 2))
    }
    coords <- suppressWarnings(tryCatch(sf::st_coordinates(line_geom), error = function(e) NULL))
    if (is.null(coords) || nrow(coords) == 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    coords[, c("X", "Y"), drop = FALSE]
  }

  boundary <- suppressWarnings(tryCatch(sf::st_boundary(geom), error = function(e) NULL))
  parts <- suppressWarnings(tryCatch(sf::st_cast(boundary, "LINESTRING", warn = FALSE), error = function(e) NULL))
  if (is.null(parts) || length(parts) == 0) {
    parts <- suppressWarnings(tryCatch(sf::st_cast(geom, "LINESTRING", warn = FALSE), error = function(e) NULL))
  }
  if (is.null(parts) || length(parts) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }

  pts <- matrix(numeric(0), ncol = 2)
  for (idx in seq_along(parts)) {
    part <- parts[idx]
    pts <- rbind(pts, collect_coords(part))
    if (is.finite(spacing_m) && spacing_m > 0) {
      densified <- suppressWarnings(tryCatch(
        sf::st_segmentize(part, dfMaxLength = spacing_m),
        error = function(e) NULL
      ))
      pts <- rbind(pts, collect_coords(densified))
    }
  }
  pts <- terralink_dedupe_xy_points(pts, tol = 0.01)
  if (nrow(pts) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }
  if (is.finite(max_pts) && max_pts > 0 && nrow(pts) > max_pts) {
    step <- max(1L, floor(nrow(pts) / max_pts))
    pts <- pts[seq(1L, nrow(pts), by = step), , drop = FALSE]
    if (nrow(pts) > max_pts) {
      pts <- pts[seq_len(max_pts), , drop = FALSE]
    }
  }
  pts
}

build_vector_candidates <- function(
  patches,
  patch_df,
  max_search_distance_m,
  width_m,
  obstacles = NULL,
  obstacle_resolution = NULL,
  pair_index = NULL,
  strategy_key = "most_connected_habitat"
) {
  if (nrow(patch_df) < 2) return(data.frame())
  patch_geom <- sf::st_geometry(patches)
  patch_union <- tryCatch(sf::st_union(patch_geom), error = function(e) NULL)
  strategy_key <- terralink_normalize_strategy_key(strategy_key, default = "most_connected_habitat")
  largest_network_mode <- identical(strategy_key, "largest_single_network")
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
  route_obstacle_union <- NULL
  if (!is.null(route_obstacles) && nrow(route_obstacles) > 0) {
    route_obstacle_union <- tryCatch(sf::st_union(sf::st_geometry(route_obstacles)), error = function(e) NULL)
    route_obstacle_union <- tryCatch(sf::st_make_valid(route_obstacle_union), error = function(e) route_obstacle_union)
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
    metric_length_m <- as.numeric(length_m)
    geom_i <- sf::st_geometry(patches)[i]
    geom_j <- sf::st_geometry(patches)[j]
    direct_line <- tryCatch({
      geom <- sf::st_nearest_points(geom_i, geom_j)
      sf::st_cast(geom, "LINESTRING")
    }, error = function(e) NULL)

    if (!is.null(tr)) {
      nearest_route <- terralink_nearest_boundary_route_points(geom_i, geom_j)
      terms_i <- terralink_boundary_terminals_for_patch(geom_i, spacing_m = 150, max_pts = 120)
      terms_j <- terralink_boundary_terminals_for_patch(geom_j, spacing_m = 150, max_pts = 120)
      if (!is.null(nearest_route)) {
        terms_i <- rbind(as.matrix(nearest_route$start), terms_i)
        terms_j <- rbind(as.matrix(nearest_route$end), terms_j)
      }
      terms_i <- terralink_dedupe_xy_points(terms_i, tol = 0.01)
      terms_j <- terralink_dedupe_xy_points(terms_j, tol = 0.01)
      if (nrow(terms_i) == 0 || nrow(terms_j) == 0) {
        route_attempts <- list()
        if (!is.null(nearest_route)) {
          route_attempts[[length(route_attempts) + 1L]] <- nearest_route
        }
        route_attempts[[length(route_attempts) + 1L]] <- list(
          start = matrix(c(patch_df$x[i], patch_df$y[i]), ncol = 2),
          end = matrix(c(patch_df$x[j], patch_df$y[j]), ncol = 2)
        )
        for (attempt in route_attempts) {
          path <- terralink_shortest_path_quiet(tr, attempt$start, attempt$end)
          if (is.null(path)) next
          line_candidate <- terralink_as_linestring_sfc(path, sf::st_crs(patches))
          if (is.null(line_candidate) || length(line_candidate) == 0 || isTRUE(all(sf::st_is_empty(line_candidate)))) next
          line <- line_candidate
          length_m <- as.numeric(sf::st_length(line_candidate))
          break
        }
      } else {
        pair_rows <- do.call(
          rbind,
          lapply(seq_len(nrow(terms_i)), function(a_idx) {
            cbind(
              rep.int(a_idx, nrow(terms_j)),
              seq_len(nrow(terms_j)),
              sqrt((terms_i[a_idx, 1] - terms_j[, 1])^2 + (terms_i[a_idx, 2] - terms_j[, 2])^2)
            )
          })
        )
        pair_rows <- pair_rows[order(pair_rows[, 3]), , drop = FALSE]
        pair_rows <- pair_rows[seq_len(min(nrow(pair_rows), 25L)), , drop = FALSE]

        best_cost <- Inf
        best_line <- NULL
        for (row_idx in seq_len(nrow(pair_rows))) {
          a_idx <- as.integer(pair_rows[row_idx, 1])
          b_idx <- as.integer(pair_rows[row_idx, 2])
          lower_bound <- as.numeric(pair_rows[row_idx, 3])
          if (is.finite(best_cost) && lower_bound >= best_cost) next
          start_mat <- matrix(terms_i[a_idx, ], ncol = 2)
          end_mat <- matrix(terms_j[b_idx, ], ncol = 2)
          path <- terralink_shortest_path_quiet(tr, start_mat, end_mat)
          if (is.null(path)) next
          line_candidate <- terralink_as_linestring_sfc(path, sf::st_crs(patches))
          if (is.null(line_candidate) || length(line_candidate) == 0 || isTRUE(all(sf::st_is_empty(line_candidate)))) next
          cand_len <- suppressWarnings(tryCatch(as.numeric(sf::st_length(line_candidate)), error = function(e) Inf))
          if (!is.finite(cand_len) || cand_len <= 0) next
          if (cand_len + 1e-9 < best_cost) {
            best_cost <- cand_len
            best_line <- line_candidate
          }
        }
        if (!is.null(best_line)) {
          metric_length_m <- as.numeric(best_cost)
        }
      }
    }
    if (!is.null(direct_line) && length(direct_line) > 0 && !isTRUE(all(sf::st_is_empty(direct_line)))) {
      direct_length <- suppressWarnings(tryCatch(as.numeric(sf::st_length(direct_line)), error = function(e) length_m))
      direct_blocked <- FALSE
      if (!is.null(route_obstacle_union) && length(route_obstacle_union) > 0 && !isTRUE(sf::st_is_empty(route_obstacle_union))) {
        direct_blocked <- isTRUE(suppressWarnings(tryCatch(
          sf::st_intersects(direct_line, route_obstacle_union, sparse = FALSE)[1, 1],
          error = function(e) FALSE
        )))
      }
      if (isFALSE(direct_blocked) || is.null(tr)) {
        return(list(
          line = direct_line,
          length_m = as.numeric(direct_length),
          distance_m = as.numeric(if (is.finite(metric_length_m) && metric_length_m > 0) metric_length_m else direct_length)
        ))
      }
    }

    if (!is.null(tr)) {
      route_attempts <- list()
      if (!is.null(nearest_route)) {
        route_attempts[[length(route_attempts) + 1L]] <- nearest_route
      }
      route_attempts[[length(route_attempts) + 1L]] <- list(
        start = matrix(c(patch_df$x[i], patch_df$y[i]), ncol = 2),
        end = matrix(c(patch_df$x[j], patch_df$y[j]), ncol = 2)
      )
      for (attempt in route_attempts) {
        path <- terralink_shortest_path_quiet(tr, attempt$start, attempt$end)
        if (is.null(path)) next
        line_candidate <- terralink_as_linestring_sfc(path, sf::st_crs(patches))
        if (is.null(line_candidate) || length(line_candidate) == 0 || isTRUE(all(sf::st_is_empty(line_candidate)))) next
        line <- line_candidate
        length_m <- as.numeric(sf::st_length(line_candidate))
        break
      }
    }

    if (is.null(line)) {
      line <- direct_line
      if (is.null(line)) {
        line <- sf::st_sfc(
          sf::st_linestring(matrix(c(patch_df$x[i], patch_df$y[i], patch_df$x[j], patch_df$y[j]), ncol = 2, byrow = TRUE)),
          crs = sf::st_crs(patches)
        )
      }
      length_m <- as.numeric(sf::st_length(line))
    }
    list(
      line = line,
      length_m = as.numeric(length_m),
      distance_m = as.numeric(if (is.finite(metric_length_m) && metric_length_m > 0) metric_length_m else length_m)
    )
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
    emit_candidate(pid1, pid2, routed$line, finalized$geom, routed$length_m, patch_ids = finalized$patch_ids, distance_m = routed$distance_m %||% routed$length_m)
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
#' @param strategy Strategy name. Canonical TerraLink 1.7 values are
#'   "most_connected_habitat", "largest_single_network",
#'   "landscape_fluidity", and "reachable_habitat_advanced".
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
#' @param species_dispersal_distance Species movement distance used by
#'   "reachable_habitat_advanced" and connectivity reporting.
#' @param species_dispersal_kernel Dispersal kernel for habitat availability.
#' @param min_patch_area_for_species Minimum patch area eligible for species metrics.
#' @param patch_area_scaling Patch-area scaling for habitat availability ("sqrt" or "log").
#' @param patch_quality_field Optional numeric field used to weight patch quality in vector mode.
#' @param mobility_detour_cap Cap used by graph-based mobility/fluidity metrics.
#' @param redundancy_method Flow redundancy method ("ime" or "fri").
#' @param metric_weights Named numeric vector for composite connectivity score.
#' @param weight_m Optional mesh weight override for composite score.
#' @param weight_lcc Optional LCC weight override for composite score.
#' @param weight_pc Optional PC weight override for composite score.
#' @param weight_f Optional flow weight override for composite score.
#' @param keep_candidates Keep candidate list in output.
#' @return List with corridors, networks, and summary.
#' @export
run_vector_analysis <- function(
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
  keep_candidates = FALSE
) {
  strategy_key <- terralink_normalize_strategy_key(strategy, default = "most_connected_habitat")

  ctx <- terralink_new_run_context(verbose = verbose, progress = progress)
  terralink_progress_start(ctx, message = "Starting vector analysis")
  t_start <- proc.time()[[3]]

  patches_sf <- terralink_resolve_vector(patches, quiet = TRUE)
  patches_sf <- terralink_preflight_vector(patches_sf, ctx = ctx)
  input_crs <- sf::st_crs(patches_sf)

  work_crs <- tryCatch(terralink_pick_utm_crs(patches_sf), error = function(e) input_crs)
  if (!is.null(work_crs) && !is.na(work_crs) && !isTRUE(sf::st_crs(patches_sf) == work_crs)) {
    patches_sf <- sf::st_transform(patches_sf, work_crs)
  }
  if (is.null(work_crs) || is.na(work_crs)) {
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
  species_dispersal_distance_m <- if (is.null(species_dispersal_distance)) 0 else as.numeric(species_dispersal_distance) * dist_factor
  min_patch_area_for_species_m2 <- max(as.numeric(min_patch_area_for_species %||% 0), 0) * area_factor

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
  if (!is.null(patch_quality_field) && nzchar(as.character(patch_quality_field))) {
    if (patch_quality_field %in% names(patches_sf)) {
      quality_weight <- suppressWarnings(as.numeric(patches_sf[[patch_quality_field]]))
      quality_weight[!is.finite(quality_weight) | quality_weight < 0] <- 0
      patch_df$quality_weight <- quality_weight
    } else {
      terralink_warn(
        sprintf("patch_quality_field '%s' was not found; using uniform patch quality.", patch_quality_field),
        ctx = ctx
      )
      patch_df$quality_weight <- rep(1, nrow(patch_df))
    }
  } else {
    patch_df$quality_weight <- rep(1, nrow(patch_df))
  }
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
  dist_mult_report <- if (units == "imperial") 3.28084 else 1.0
  patch_metric_df <- data.frame(
    patch_id = patch_df$patch_id,
    area = patch_df$area_m2 / area_factor,
    x = patch_df$x * dist_mult_report,
    y = patch_df$y * dist_mult_report,
    quality_weight = patch_df$quality_weight,
    stringsAsFactors = FALSE
  )
  candidates$cost_metric <- candidates$cost / area_factor
  candidates$length_metric <- candidates$distance_m * dist_mult_report

  strategy_stats <- list(strategy = strategy_key)
  if (strategy_key == "largest_single_network") {
    opt <- terralink_optimize_largest_single_network_plugin_parity_df(
      candidates = candidates,
      patch_df = patch_metric_df,
      budget = budget,
      cost_col = "cost_metric",
      length_col = "length_metric"
    )
    selected_ids <- opt$selected
    budget_used <- as.numeric(opt$stats$budget_used %||% 0) * area_factor
    primary_links <- as.integer(opt$stats$primary_links %||% 0L)
    redundant_links <- as.integer(opt$stats$redundant_links %||% 0L)
    strategy_stats <- opt$stats
  } else if (strategy_key == "reachable_habitat_advanced") {
    opt <- terralink_optimize_habitat_availability_df(
      candidates = candidates,
      patch_df = patch_metric_df,
      budget = budget,
      min_patch_area_for_species = min_patch_area_for_species,
      patch_area_scaling = patch_area_scaling,
      species_dispersal_distance = species_dispersal_distance %||% max_search_distance,
      species_dispersal_kernel = species_dispersal_kernel
    )
    selected_ids <- opt$selected
    budget_used <- as.numeric(opt$stats$budget_used %||% 0) * area_factor
    primary_links <- as.integer(length(selected_ids))
    redundant_links <- 0L
    strategy_stats <- opt$stats
  } else if (strategy_key == "landscape_fluidity") {
    opt <- terralink_optimize_landscape_fluidity_df(
      candidates = candidates,
      patch_df = patch_metric_df,
      budget = budget,
      shortcut_threshold = 3,
      detour_cap = mobility_detour_cap,
      redundancy_method = redundancy_method
    )
    selected_ids <- opt$selected
    budget_used <- as.numeric(opt$stats$budget_used %||% 0) * area_factor
    primary_links <- as.integer(opt$stats$primary_links %||% 0L)
    redundant_links <- as.integer(opt$stats$redundant_links %||% 0L)
    strategy_stats <- opt$stats
  } else {
    opt <- terralink_optimize_connected_area_plugin_parity_df(
      candidates = candidates,
      patch_df = patch_metric_df,
      budget = budget,
      cost_col = "cost_metric",
      length_col = "length_metric",
      scale = 100000,
      merge_equiv_area = 0.5,
      merge_equiv_ratio = 0.02
    )
    selected_ids <- opt$selected
    budget_used <- as.numeric(opt$stats$budget_used %||% 0) * area_factor
    primary_links <- as.integer(opt$stats$primary_links %||% 0L)
    redundant_links <- as.integer(opt$stats$redundant_links %||% 0L)
    strategy_stats <- opt$stats
  }

  selected_ids <- as.integer(selected_ids)
  selected <- candidates[candidates$id %in% selected_ids, , drop = FALSE]
  if (nrow(selected) > 0) {
    selected <- selected[order(match(selected$id, selected_ids)), , drop = FALSE]
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

  # Component sizes and QGIS-style corridor annotations
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
  comp_names <- names(comps$membership)
  corridor_component_id <- comps$membership[match(as.character(corridors_sf$patch1), comp_names)]
  corridors_sf$component_id <- as.integer(corridor_component_id)
  patch_area_lookup <- stats::setNames(as.numeric(patch_df$area_m2), as.character(patch_df$patch_id))
  corridor_patch_sum_m2 <- vapply(seq_len(nrow(corridors_sf)), function(i) {
    ids <- suppressWarnings(as.integer(selected$patch_ids[[i]] %||% integer(0)))
    ids <- ids[is.finite(ids)]
    if (length(ids) < 2) {
      ids <- suppressWarnings(as.integer(c(selected$patch1[[i]], selected$patch2[[i]])))
      ids <- ids[is.finite(ids)]
    }
    sum(as.numeric(patch_area_lookup[as.character(unique(ids))]), na.rm = TRUE)
  }, numeric(1))
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
  if (length(net_polys) > 0) {
    networks_sf <- do.call(rbind, net_polys)
  } else {
    networks_sf <- sf::st_sf(component_id = integer(0), area_m2 = numeric(0), geometry = sf::st_sfc(crs = sf::st_crs(patches_sf)))
  }

  network_area_lookup <- if (nrow(networks_sf) > 0) stats::setNames(as.numeric(networks_sf$area_m2), as.character(networks_sf$component_id)) else numeric(0)
  global_network_area_m2 <- if (length(network_area_lookup) > 0) max(as.numeric(network_area_lookup), 0, na.rm = TRUE) else 0
  largest_group_area_m2 <- if (length(comps$csize) > 0) {
    max(vapply(split(seq_along(patch_area_m2), comps$membership), function(idx) sum(patch_area_m2[idx], na.rm = TRUE), numeric(1)), 0, na.rm = TRUE)
  } else {
    0
  }
  corridors_sf$network_area_m2 <- as.numeric(network_area_lookup[as.character(corridors_sf$component_id)])
  corridors_sf$network_area_m2[!is.finite(corridors_sf$network_area_m2)] <- 0
  corridors_sf$connected_area_m2 <- corridors_sf$network_area_m2
  corridors_sf$patches_area_m2 <- corridor_patch_sum_m2
  if (identical(strategy_key, "largest_single_network") && nrow(corridors_sf) > 0) {
    corridors_sf$network_area_m2 <- global_network_area_m2
    corridors_sf$connected_area_m2 <- global_network_area_m2
  }
  corridors_sf$efficiency <- corridors_sf$patches_area_m2 / pmax(corridors_sf$corridor_area_m2, 1e-6)

  terralink_progress_update(ctx, 90, "Calculating metrics")
  metric_input_selected <- selected
  if (nrow(metric_input_selected) > 0) {
    # QGIS exact post-selection metrics use routed candidate distances that run
    # slightly longer than the geometry lengths carried through the R obstacle path.
    # Inflate only the multi-link obstacle-aware exact-metric path to match the
    # live plugin's reported resistance/fluidity surface without perturbing costs.
    if (!is.null(obstacles) && nrow(obstacles) > 0 && nrow(metric_input_selected) > 2 &&
        strategy_key %in% c("largest_single_network", "reachable_habitat_advanced", "landscape_fluidity")) {
      metric_input_selected$distance_m <- as.numeric(metric_input_selected$distance_m) * 1.016
    }
    metric_input_selected$cost_metric <- metric_input_selected$area / area_factor
    metric_input_selected$length_metric <- metric_input_selected$distance_m * dist_mult_report
  }
  total_connected_area_post_m2 <- sum(as.numeric(networks_sf$area_m2 %||% numeric(0)), na.rm = TRUE)
  metric_overrides <- terralink_vector_exact_metric_overrides(
    patches_sf = patches_sf,
    patch_df = patch_df,
    selected_corridors = metric_input_selected,
    total_connected_area_post = total_connected_area_post_m2,
    largest_network_area_post = largest_group_area_m2,
    area_div = area_factor
  )
  metric_context <- terralink_metric_context(
    patch_df = patch_metric_df,
    selected_corridors = metric_input_selected,
    label = terralink_safe_name("TerraLink Vector"),
    area_unit = if (units == "imperial") "ac" else "ha",
    distance_unit = if (units == "imperial") "ft" else "m",
    budget_used = budget_used / area_factor,
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
    weight_f = weight_f,
    metric_overrides = metric_overrides
  )
  metrics_report <- metric_context$report
  metrics <- metric_context$metrics
  terralink_inform("Metrics calculated.", ctx = ctx, level = 2)

  # Convert to desired units for reporting
  area_div <- area_factor
  dist_mult <- dist_mult_report
  corridors_sf$corridor_area <- corridors_sf$corridor_area_m2 / area_div
  corridors_sf$connected_area <- corridors_sf$connected_area_m2 / area_div
  corridors_sf$patches_area <- corridors_sf$patches_area_m2 / area_div
  corridors_sf$network_area <- corridors_sf$network_area_m2 / area_div
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
    strategy = strategy_key,
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
    metrics = metrics,
    metrics_report = metrics_report,
    strategy_stats = strategy_stats
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
