#' Create a logical mask from raster values
#'
#' @param raster SpatRaster.
#' @param values Numeric values to keep.
#' @return Logical SpatRaster mask.
#' @export
raster_mask_from_values <- function(raster, values) {
  if (length(values) == 0) {
    terralink_abort("values must contain at least one value.", class = "terralink_error_input")
  }
  values <- as.numeric(values)
  mask <- terra::app(raster, fun = function(x) x %in% values)
  mask[is.na(mask)] <- FALSE
  mask
}

#' Label contiguous habitat patches
#'
#' @param mask Logical SpatRaster mask.
#' @param connectivity 4 or 8.
#' @return SpatRaster of patch labels.
#' @export
label_patches <- function(mask, connectivity = 8) {
  directions <- if (connectivity == 8) 8 else 4
  # Ensure background is NA so only habitat cells get labeled.
  m <- mask
  if (is.logical(terra::values(m, mat = FALSE))) {
    m <- terra::ifel(m == 1, 1, NA)
  } else {
    m <- terra::ifel(m > 0, 1, NA)
  }
  labels <- terra::patches(m, directions = directions)
  labels
}

#' Summarize patch labels into a data frame
#'
#' @param labels SpatRaster of patch labels.
#' @return Data frame with patch_id, cell_count, area, x, y.
#' @export
patch_summary_from_labels <- function(labels) {
  pts <- terra::as.data.frame(labels, xy = TRUE, na.rm = TRUE)
  if (nrow(pts) == 0) {
    return(data.frame(patch_id = integer(0), cell_count = integer(0), area = numeric(0), x = numeric(0), y = numeric(0)))
  }
  value_cols <- setdiff(names(pts), c("x", "y"))
  if (length(value_cols) != 1) {
    terralink_abort("Unable to determine patch label column.", class = "terralink_error_input")
  }
  value_col <- value_cols[[1]]
  pts <- pts[!is.na(pts[[value_col]]) & pts[[value_col]] != 0, , drop = FALSE]
  if (nrow(pts) == 0) {
    return(data.frame(patch_id = integer(0), cell_count = integer(0), area = numeric(0), x = numeric(0), y = numeric(0)))
  }

  ids <- pts[[value_col]]
  cell_count <- tapply(pts$x, ids, length)
  centroids_x <- tapply(pts$x, ids, mean)
  centroids_y <- tapply(pts$y, ids, mean)

  cell_area <- abs(prod(terra::res(labels)))
  patch_ids <- names(cell_count)

  data.frame(
    patch_id = as.integer(patch_ids),
    cell_count = as.integer(cell_count[patch_ids]),
    area = as.numeric(cell_count[patch_ids]) * cell_area,
    x = as.numeric(centroids_x[patch_ids]),
    y = as.numeric(centroids_y[patch_ids])
  )
}

#' Build candidate edges between patches (centroid distance)
#'
#' @param patch_df Patch summary data frame.
#' @param max_search_distance Maximum distance in pixels.
#' @param raster_ref Raster reference for pixel size.
#' @return Data frame with patch1, patch2, cost, distance_map, id.
#' @export
build_patch_candidates <- function(patch_df, max_search_distance, raster_ref) {
  if (nrow(patch_df) < 2) {
    return(data.frame(patch1 = integer(0), patch2 = integer(0), cost = numeric(0), distance_map = numeric(0), id = integer(0)))
  }
  res_xy <- terra::res(raster_ref)
  pixel_size <- max(res_xy)

  edges <- list()
  k <- 1
  for (i in seq_len(nrow(patch_df) - 1)) {
    for (j in (i + 1):nrow(patch_df)) {
      dx <- patch_df$x[i] - patch_df$x[j]
      dy <- patch_df$y[i] - patch_df$y[j]
      dist_map <- sqrt(dx * dx + dy * dy)
      dist_px <- dist_map / pixel_size
      if (dist_px <= max_search_distance) {
        edges[[k]] <- data.frame(
          patch1 = patch_df$patch_id[i],
          patch2 = patch_df$patch_id[j],
          cost = dist_px,
          distance_map = dist_map,
          id = k
        )
        k <- k + 1
      }
    }
  }

  if (length(edges) == 0) {
    return(data.frame(patch1 = integer(0), patch2 = integer(0), cost = numeric(0), distance_map = numeric(0), id = integer(0)))
  }
  do.call(rbind, edges)
}

#' Create corridor raster from selected edges
#'
#' @param labels SpatRaster labels.
#' @param patch_df Patch summary data frame.
#' @param corridors Data frame with patch1, patch2, and optional line geometry.
#' @param min_corridor_width_px Width (pixels) for buffering corridors.
#' @return SpatRaster with corridor cells set to 1.
#' @export
build_corridor_raster <- function(labels, patch_df, corridors, min_corridor_width_px = 1) {
  out <- labels
  vals <- terra::values(out)
  vals[] <- NA

  if (is.null(corridors) || nrow(corridors) == 0) {
    terra::values(out) <- vals
    return(out)
  }

  if ("buffered_cells" %in% names(corridors)) {
    for (i in seq_len(nrow(corridors))) {
      cells <- corridors$buffered_cells[[i]]
      if (is.list(cells) && length(cells) == 1) cells <- cells[[1]]
      if (is.null(cells) || length(cells) == 0) next
      cells <- as.integer(cells)
      cells <- unique(cells[is.finite(cells) & cells >= 1 & cells <= length(vals)])
      if (length(cells) == 0) next
      score <- 1
      if ("connected_size" %in% names(corridors)) {
        score <- suppressWarnings(as.numeric(corridors$connected_size[[i]]))
        if (!is.finite(score) || score <= 0) score <- 1
      }
      old <- vals[cells]
      old[is.na(old)] <- -Inf
      vals[cells] <- pmax(old, score)
    }
    terra::values(out) <- vals
    return(out)
  }

  pixel_size <- max(terra::res(labels))
  width_map <- max(1, min_corridor_width_px) * pixel_size

  if ("line" %in% names(corridors)) {
    lines <- corridors$line
    if (!inherits(lines, "sfc")) {
      lines <- sf::st_sfc(lines, crs = terra::crs(labels))
    }
    line_vec <- terra::vect(lines)
    if (min_corridor_width_px > 1) {
      line_vec <- terra::buffer(line_vec, width = width_map / 2)
    }
    out <- terra::rasterize(line_vec, out, field = 1, background = NA)
    return(out)
  }

  coords <- patch_df[, c("patch_id", "x", "y")]
  get_xy <- function(pid) {
    row <- coords[coords$patch_id == pid, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    c(row$x[1], row$y[1])
  }

  for (i in seq_len(nrow(corridors))) {
    p1 <- corridors$patch1[i]
    p2 <- corridors$patch2[i]
    xy1 <- get_xy(p1)
    xy2 <- get_xy(p2)
    if (is.null(xy1) || is.null(xy2)) {
      next
    }
    dx <- xy2[1] - xy1[1]
    dy <- xy2[2] - xy1[2]
    steps <- max(2, ceiling(max(abs(dx), abs(dy)) / pixel_size))
    t_vals <- seq(0, 1, length.out = steps)
    xs <- xy1[1] + dx * t_vals
    ys <- xy1[2] + dy * t_vals
    cells <- terra::cellFromXY(labels, cbind(xs, ys))
    cells <- unique(cells[!is.na(cells)])
    if (length(cells) > 0) {
      vals[cells] <- 1
    }
  }

  terra::values(out) <- vals
  if (min_corridor_width_px > 1) {
    w <- matrix(1, nrow = min_corridor_width_px, ncol = min_corridor_width_px)
    out <- terra::focal(out, w = w, fun = max, na.policy = "omit")
  }
  out
}

terralink_shift_bool_mask <- function(mask, dr, dc) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  out <- matrix(FALSE, nrow = nr, ncol = nc)
  if (abs(dr) >= nr || abs(dc) >= nc) return(out)

  if (dr >= 0) {
    src_r <- (1 + dr):nr
    dst_r <- 1:(nr - dr)
  } else {
    src_r <- 1:(nr + dr)
    dst_r <- (1 - dr):nr
  }
  if (dc >= 0) {
    src_c <- (1 + dc):nc
    dst_c <- 1:(nc - dc)
  } else {
    src_c <- 1:(nc + dc)
    dst_c <- (1 - dc):nc
  }

  out[dst_r, dst_c] <- mask[src_r, src_c]
  out
}

terralink_boundary_coords_by_label <- function(labels_mat, patch_ids, connectivity = 4) {
  patch_ids <- as.integer(patch_ids)
  out <- vector("list", length(patch_ids))
  names(out) <- as.character(patch_ids)
  offsets <- if (as.integer(connectivity) == 8L) {
    rbind(
      c(-1L, -1L), c(-1L, 0L), c(-1L, 1L),
      c(0L, -1L), c(0L, 1L),
      c(1L, -1L), c(1L, 0L), c(1L, 1L)
    )
  } else {
    rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  }

  for (pid in patch_ids) {
    inside <- !is.na(labels_mat) & (labels_mat == pid)
    if (!any(inside, na.rm = TRUE)) {
      out[[as.character(pid)]] <- matrix(integer(0), ncol = 2)
      next
    }
    interior <- inside
    for (k in seq_len(nrow(offsets))) {
      interior <- interior & terralink_shift_bool_mask(inside, offsets[k, 1], offsets[k, 2])
    }
    boundary <- inside & !interior
    out[[as.character(pid)]] <- which(boundary, arr.ind = TRUE)
  }
  out
}

terralink_subsample_coords_grid <- function(coords, block_size = 2L, max_points = 2000L) {
  if (is.null(coords) || nrow(coords) == 0) return(coords)
  bs <- as.integer(max(1L, block_size))
  rr <- (coords[, 1] - 1L) %/% bs
  cc <- (coords[, 2] - 1L) %/% bs
  key <- paste(rr, cc, sep = ":")
  keep <- !duplicated(key)
  sub <- coords[keep, , drop = FALSE]
  if (!is.null(max_points) && nrow(sub) > as.integer(max_points)) {
    stride <- as.integer(max(1L, floor(nrow(sub) / as.integer(max_points))))
    sub <- sub[seq(1L, nrow(sub), by = stride), , drop = FALSE]
    if (nrow(sub) > as.integer(max_points)) {
      sub <- sub[seq_len(as.integer(max_points)), , drop = FALSE]
    }
  }
  sub
}

terralink_pair_key <- function(a, b) {
  a <- as.integer(a)
  b <- as.integer(b)
  if (a <= b) paste(a, b, sep = "_") else paste(b, a, sep = "_")
}

terralink_generate_boundary_pair_seeds <- function(
  labels_mat,
  patch_ids,
  connectivity,
  max_search_distance,
  min_corridor_width,
  max_seeds_per_pair = 8L
) {
  max_d <- as.integer(max(1L, round(max_search_distance)))
  max_d2 <- as.numeric(max_d) * as.numeric(max_d)
  patch_ids <- as.integer(patch_ids)
  boundary_by_pid <- terralink_boundary_coords_by_label(labels_mat, patch_ids, connectivity = connectivity)

  pts <- vector("list", length(patch_ids))
  npts <- 0L
  for (i in seq_along(patch_ids)) {
    pid <- patch_ids[[i]]
    coords <- boundary_by_pid[[as.character(pid)]]
    if (is.null(coords) || nrow(coords) == 0) next
    sub <- terralink_subsample_coords_grid(coords, block_size = 2L, max_points = 2000L)
    if (nrow(sub) == 0) next
    pts[[i]] <- data.frame(
      pid = rep.int(as.integer(pid), nrow(sub)),
      row = as.integer(sub[, 1]),
      col = as.integer(sub[, 2])
    )
    npts <- npts + nrow(sub)
  }
  if (npts == 0L) return(list())
  pts <- do.call(rbind, pts[!vapply(pts, is.null, logical(1))])
  rownames(pts) <- NULL

  bin_size <- as.integer(max(8L, min(64L, max_d)))
  bx <- (pts$row - 1L) %/% bin_size
  by <- (pts$col - 1L) %/% bin_size
  bin_key <- paste(bx, by, sep = ":")
  bins <- split(seq_len(nrow(pts)), bin_key)
  bin_parts <- do.call(rbind, strsplit(names(bins), ":", fixed = TRUE))
  bin_x <- as.integer(bin_parts[, 1])
  bin_y <- as.integer(bin_parts[, 2])
  bin_index <- seq_along(bins)

  pair_cands <- list()
  for (bi in seq_along(bins)) {
    idx_a <- bins[[bi]]
    if (length(idx_a) == 0) next
    nb_mask <- abs(bin_x - bin_x[[bi]]) <= 1L & abs(bin_y - bin_y[[bi]]) <= 1L
    idx_b <- unique(unlist(bins[bin_index[nb_mask]], use.names = FALSE))
    if (length(idx_b) == 0) next

    for (ia in idx_a) {
      p1 <- pts$pid[[ia]]
      r1 <- pts$row[[ia]]
      c1 <- pts$col[[ia]]
      for (ib in idx_b) {
        p2 <- pts$pid[[ib]]
        if (p2 <= p1) next
        dr <- as.numeric(r1 - pts$row[[ib]])
        dc <- as.numeric(c1 - pts$col[[ib]])
        d2 <- dr * dr + dc * dc
        if (d2 > max_d2) next
        key <- terralink_pair_key(p1, p2)
        row_new <- c(d2, r1, c1, pts$row[[ib]], pts$col[[ib]])
        cur <- pair_cands[[key]]
        if (is.null(cur)) {
          cur <- matrix(row_new, ncol = 5)
        } else {
          cur <- rbind(cur, row_new)
          if (nrow(cur) > 200) {
            ord <- order(cur[, 1], decreasing = FALSE)
            cur <- cur[ord[seq_len(min(50L, length(ord)))], , drop = FALSE]
          }
        }
        pair_cands[[key]] <- cur
      }
    }
  }

  if (length(pair_cands) == 0) return(list())
  diversity_radius <- as.integer(max(4L, ceiling(as.numeric(min_corridor_width) * 1.5)))
  out <- list()
  out_k <- 1L
  for (key in names(pair_cands)) {
    vals <- pair_cands[[key]]
    if (is.null(vals) || nrow(vals) == 0) next
    vals <- vals[order(vals[, 1], decreasing = FALSE), , drop = FALSE]
    chosen <- matrix(integer(0), ncol = 4)
    chosen_d <- numeric(0)
    for (i in seq_len(nrow(vals))) {
      if (nrow(chosen) >= as.integer(max_seeds_per_pair)) break
      a <- as.integer(vals[i, 2:3])
      b <- as.integer(vals[i, 4:5])
      ok <- TRUE
      if (nrow(chosen) > 0) {
        for (j in seq_len(nrow(chosen))) {
          if (abs(a[[1]] - chosen[j, 1]) <= diversity_radius && abs(a[[2]] - chosen[j, 2]) <= diversity_radius) {
            ok <- FALSE
            break
          }
          if (abs(b[[1]] - chosen[j, 3]) <= diversity_radius && abs(b[[2]] - chosen[j, 4]) <= diversity_radius) {
            ok <- FALSE
            break
          }
        }
      }
      if (!ok) next
      chosen <- rbind(chosen, c(a[[1]], a[[2]], b[[1]], b[[2]]))
      chosen_d <- c(chosen_d, sqrt(as.numeric(vals[i, 1])))
    }
    if (nrow(chosen) == 0) next
    pair_ids <- as.integer(strsplit(key, "_", fixed = TRUE)[[1]])
    out[[out_k]] <- list(
      patch1 = pair_ids[[1]],
      patch2 = pair_ids[[2]],
      seeds = data.frame(
        r1 = as.integer(chosen[, 1]),
        c1 = as.integer(chosen[, 2]),
        r2 = as.integer(chosen[, 3]),
        c2 = as.integer(chosen[, 4]),
        d = as.numeric(chosen_d)
      )
    )
    out_k <- out_k + 1L
  }
  out
}

terralink_rc_to_cells <- function(ncol, rc) {
  ((as.integer(rc[, 1]) - 1L) * as.integer(ncol)) + as.integer(rc[, 2])
}

terralink_cells_to_rc <- function(ncol, cells) {
  cells <- as.integer(cells)
  r <- ((cells - 1L) %/% as.integer(ncol)) + 1L
  c <- ((cells - 1L) %% as.integer(ncol)) + 1L
  cbind(r, c)
}

terralink_bresenham_rc <- function(r1, c1, r2, c2) {
  x1 <- as.integer(c1)
  y1 <- as.integer(r1)
  x2 <- as.integer(c2)
  y2 <- as.integer(r2)
  dx <- abs(x2 - x1)
  sx <- if (x1 < x2) 1L else -1L
  dy <- -abs(y2 - y1)
  sy <- if (y1 < y2) 1L else -1L
  err <- dx + dy
  out <- matrix(integer(0), ncol = 2)
  repeat {
    out <- rbind(out, c(y1, x1))
    if (x1 == x2 && y1 == y2) break
    e2 <- 2L * err
    if (e2 >= dy) {
      err <- err + dy
      x1 <- x1 + sx
    }
    if (e2 <= dx) {
      err <- err + dx
      y1 <- y1 + sy
    }
  }
  out
}

terralink_corridor_offsets <- function(width_px) {
  width <- max(1L, as.integer(round(width_px)))
  if (width <= 1L) return(matrix(c(0L, 0L), ncol = 2))
  radius <- max(0, width / 2)
  max_off <- as.integer(ceiling(radius))
  radius_sq <- radius * radius
  out <- matrix(integer(0), ncol = 2)
  for (dr in seq.int(-max_off, max_off)) {
    for (dc in seq.int(-max_off, max_off)) {
      if ((dr * dr + dc * dc) <= radius_sq + 1e-9) {
        out <- rbind(out, c(as.integer(dr), as.integer(dc)))
      }
    }
  }
  if (nrow(out) == 0) out <- matrix(c(0L, 0L), ncol = 2)
  out
}

terralink_offsets_kernel <- function(offsets) {
  if (is.null(offsets) || nrow(offsets) == 0) return(matrix(1, nrow = 1, ncol = 1))
  max_off <- max(abs(as.integer(offsets)))
  n <- as.integer(max_off * 2L + 1L)
  k <- matrix(0, nrow = n, ncol = n)
  center <- as.integer(max_off + 1L)
  for (i in seq_len(nrow(offsets))) {
    rr <- center + as.integer(offsets[i, 1])
    cc <- center + as.integer(offsets[i, 2])
    if (rr >= 1L && rr <= n && cc >= 1L && cc <= n) {
      k[rr, cc] <- 1
    }
  }
  if (!any(k > 0, na.rm = TRUE)) k[center, center] <- 1
  k
}

terralink_nearest_passable_rc <- function(r, c, nrow, ncol, passable_mask, search_radius = 6L) {
  r <- as.integer(r)
  c <- as.integer(c)
  if (!is.finite(r) || !is.finite(c)) return(c(NA_integer_, NA_integer_))
  if (r >= 1L && r <= nrow && c >= 1L && c <= ncol) {
    idx <- terralink_rc_to_cells(ncol, matrix(c(r, c), ncol = 2))
    if (length(idx) == 1L && !is.na(idx) && isTRUE(passable_mask[[idx]])) {
      return(c(r, c))
    }
  }
  max_r <- as.integer(max(1L, search_radius))
  for (rad in seq_len(max_r)) {
    rr <- seq.int(max(1L, r - rad), min(nrow, r + rad))
    cc <- seq.int(max(1L, c - rad), min(ncol, c + rad))
    grid <- as.matrix(expand.grid(rr = rr, cc = cc))
    if (nrow(grid) == 0) next
    d2 <- (grid[, 1] - r)^2 + (grid[, 2] - c)^2
    ord <- order(d2)
    grid <- grid[ord, , drop = FALSE]
    cells <- terralink_rc_to_cells(ncol, grid)
    keep <- passable_mask[cells]
    if (any(keep, na.rm = TRUE)) {
      first <- which(keep)[1]
      return(c(as.integer(grid[first, 1]), as.integer(grid[first, 2])))
    }
  }
  c(NA_integer_, NA_integer_)
}

terralink_inflate_cells <- function(cells, nrow, ncol, offsets, passable_mask = NULL) {
  if (length(cells) == 0 || nrow(offsets) == 0) return(integer(0))
  rc <- terralink_cells_to_rc(ncol, cells)
  out <- integer(0)
  for (k in seq_len(nrow(offsets))) {
    rr <- rc[, 1] + offsets[k, 1]
    cc <- rc[, 2] + offsets[k, 2]
    keep <- rr >= 1L & rr <= nrow & cc >= 1L & cc <= ncol
    if (!any(keep)) next
    out <- c(out, ((rr[keep] - 1L) * ncol) + cc[keep])
  }
  out <- unique(as.integer(out))
  if (!is.null(passable_mask) && length(passable_mask) >= max(out)) {
    out <- out[passable_mask[out]]
  }
  out
}

terralink_expand_cells <- function(cells, nrow, ncol, radius) {
  rad <- as.integer(max(1L, ceiling(radius)))
  offsets <- as.matrix(expand.grid(dr = seq.int(-rad, rad), dc = seq.int(-rad, rad)))
  terralink_inflate_cells(cells, nrow = nrow, ncol = ncol, offsets = offsets, passable_mask = NULL)
}

terralink_line_from_cells <- function(raster_ref, cells) {
  xy <- terra::xyFromCell(raster_ref, cells)
  if (is.null(xy) || nrow(xy) == 0) {
    xy <- matrix(c(0, 0, 0, 0), ncol = 2, byrow = TRUE)
  }
  if (nrow(xy) >= 2) {
    dup <- c(FALSE, (diff(xy[, 1]) == 0 & diff(xy[, 2]) == 0))
    xy <- xy[!dup, , drop = FALSE]
  }
  if (nrow(xy) < 2) {
    xy <- rbind(xy[1, , drop = FALSE], xy[1, , drop = FALSE] + 1e-9)
  }
  sf::st_sfc(sf::st_linestring(as.matrix(xy)), crs = terra::crs(raster_ref))
}

terralink_path_cells_from_line <- function(line, raster_ref) {
  if (is.null(line) || length(line) == 0 || isTRUE(sf::st_is_empty(line))) return(integer(0))
  coords <- tryCatch(sf::st_coordinates(line), error = function(e) NULL)
  if (is.null(coords) || nrow(coords) == 0) return(integer(0))
  xy <- as.matrix(coords[, seq_len(min(2, ncol(coords))), drop = FALSE])
  if (ncol(xy) < 2) return(integer(0))
  cells <- terra::cellFromXY(raster_ref, xy)
  cells <- as.integer(cells[is.finite(cells) & !is.na(cells)])
  cells <- unique(cells)
  if (length(cells) <= 1) return(cells)

  nc <- terra::ncol(raster_ref)
  rc <- terralink_cells_to_rc(nc, cells)
  path_rc <- rc[1, , drop = FALSE]
  for (idx in seq.int(2, nrow(rc))) {
    seg <- terralink_bresenham_rc(rc[idx - 1, 1], rc[idx - 1, 2], rc[idx, 1], rc[idx, 2])
    if (nrow(seg) <= 1) next
    path_rc <- rbind(path_rc, seg[-1, , drop = FALSE])
  }
  terralink_rc_to_cells(nc, path_rc)
}

terralink_shortest_path_cells <- function(transition_obj, raster_ref, from_cell, to_cell, fallback_rc = NULL) {
  if (is.null(transition_obj)) {
    if (is.null(fallback_rc)) return(list(cells = integer(0), line = NULL))
    nc <- terra::ncol(raster_ref)
    return(list(cells = terralink_rc_to_cells(nc, fallback_rc), line = NULL))
  }
  xy <- terra::xyFromCell(raster_ref, c(from_cell, to_cell))
  if (is.null(xy) || nrow(xy) < 2) {
    if (is.null(fallback_rc)) return(list(cells = integer(0), line = NULL))
    nc <- terra::ncol(raster_ref)
    return(list(cells = terralink_rc_to_cells(nc, fallback_rc), line = NULL))
  }

  path <- tryCatch(
    gdistance::shortestPath(transition_obj, origin = xy[1, , drop = FALSE], goal = xy[2, , drop = FALSE], output = "SpatialLines"),
    error = function(e) NULL
  )

  line <- NULL
  if (!is.null(path)) {
    line <- tryCatch(sf::st_geometry(sf::st_as_sf(path)), error = function(e) NULL)
    if (!is.null(line) && length(line) > 1) {
      line <- tryCatch(sf::st_line_merge(sf::st_union(line)), error = function(e) line)
      if (inherits(line, "sfc") && length(line) > 1) {
        line <- line[1]
      }
    }
  }

  cells <- terralink_path_cells_from_line(line, raster_ref)
  if (length(cells) == 0 && !is.null(fallback_rc)) {
    cells <- terralink_rc_to_cells(terra::ncol(raster_ref), fallback_rc)
  }
  list(cells = as.integer(cells), line = line)
}

terralink_neighbor_moves <- function(connectivity = 4L) {
  if (as.integer(connectivity) == 4L) {
    return(matrix(c(
      -1L, 0L,
      1L, 0L,
      0L, -1L,
      0L, 1L
    ), ncol = 2, byrow = TRUE))
  }
  matrix(c(
    -1L, -1L,
    -1L, 0L,
    -1L, 1L,
    0L, -1L,
    0L, 1L,
    1L, -1L,
    1L, 0L,
    1L, 1L
  ), ncol = 2, byrow = TRUE)
}

terralink_compute_start_positions_by_patch <- function(
  labels_mat,
  habitat_mat,
  obstacle_block_mat,
  passable_mat,
  connectivity = 4L
) {
  nr <- nrow(labels_mat)
  nc <- ncol(labels_mat)
  moves <- terralink_neighbor_moves(connectivity)
  out <- list()
  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      if (isTRUE(habitat_mat[r, c])) next
      if (isTRUE(obstacle_block_mat[r, c])) next
      if (!isTRUE(passable_mat[r, c])) next
      found <- NA_integer_
      for (k in seq_len(nrow(moves))) {
        rr <- r + moves[k, 1]
        cc <- c + moves[k, 2]
        if (rr < 1L || rr > nr || cc < 1L || cc > nc) next
        pid <- labels_mat[rr, cc]
        if (is.finite(pid) && pid > 0) {
          found <- as.integer(pid)
          break
        }
      }
      if (!is.finite(found)) next
      key <- as.character(found)
      row <- matrix(c(as.integer(r), as.integer(c)), ncol = 2)
      if (is.null(out[[key]])) {
        out[[key]] <- row
      } else {
        out[[key]] <- rbind(out[[key]], row)
      }
    }
  }
  out
}

terralink_filter_start_positions_near <- function(starts_rc, center_rc, radius) {
  if (is.null(starts_rc) || nrow(starts_rc) == 0) return(matrix(integer(0), ncol = 2))
  rad <- as.integer(max(1L, radius))
  cr <- as.integer(center_rc[[1]])
  cc <- as.integer(center_rc[[2]])
  keep <- abs(starts_rc[, 1] - cr) <= rad & abs(starts_rc[, 2] - cc) <= rad
  starts_rc[keep, , drop = FALSE]
}

terralink_find_shortest_corridor_targeted <- function(
  start_patch_id,
  target_patch_id,
  start_positions_rc,
  labels_mat,
  habitat_mat,
  max_width,
  connectivity = 4L,
  obstacle_block_mat,
  passable_mat,
  target_boundary_center = NULL,
  target_boundary_radius = NULL
) {
  if (is.null(start_positions_rc) || nrow(start_positions_rc) == 0) return(NULL)
  nr <- nrow(labels_mat)
  nc <- ncol(labels_mat)
  ncell <- nr * nc
  sqrt2 <- sqrt(2)
  max_width <- as.numeric(max_width)
  moves <- terralink_neighbor_moves(connectivity)

  best_cost <- rep(Inf, ncell)
  parent_cell <- integer(ncell)
  state_kind <- integer(ncell) # 0 path, 1 habitat-bridge

  heap_cost <- numeric(0)
  heap_r <- integer(0)
  heap_c <- integer(0)
  heap_cell <- integer(0)
  heap_size <- 0L

  heap_less <- function(i, j) {
    if (heap_cost[[i]] < heap_cost[[j]]) return(TRUE)
    if (heap_cost[[i]] > heap_cost[[j]]) return(FALSE)
    if (heap_r[[i]] < heap_r[[j]]) return(TRUE)
    if (heap_r[[i]] > heap_r[[j]]) return(FALSE)
    heap_c[[i]] < heap_c[[j]]
  }

  heap_swap <- function(i, j) {
    tc <- heap_cost[[i]]
    tr <- heap_r[[i]]
    tcc <- heap_c[[i]]
    tcell <- heap_cell[[i]]
    heap_cost[[i]] <<- heap_cost[[j]]
    heap_r[[i]] <<- heap_r[[j]]
    heap_c[[i]] <<- heap_c[[j]]
    heap_cell[[i]] <<- heap_cell[[j]]
    heap_cost[[j]] <<- tc
    heap_r[[j]] <<- tr
    heap_c[[j]] <<- tcc
    heap_cell[[j]] <<- tcell
  }

  heap_push <- function(cost, r, c, cell) {
    heap_size <<- heap_size + 1L
    if (length(heap_cost) < heap_size) {
      heap_cost <<- c(heap_cost, as.numeric(cost))
      heap_r <<- c(heap_r, as.integer(r))
      heap_c <<- c(heap_c, as.integer(c))
      heap_cell <<- c(heap_cell, as.integer(cell))
    } else {
      heap_cost[[heap_size]] <<- as.numeric(cost)
      heap_r[[heap_size]] <<- as.integer(r)
      heap_c[[heap_size]] <<- as.integer(c)
      heap_cell[[heap_size]] <<- as.integer(cell)
    }
    i <- heap_size
    while (i > 1L) {
      p <- i %/% 2L
      if (!heap_less(i, p)) break
      heap_swap(i, p)
      i <- p
    }
  }

  heap_pop <- function() {
    if (heap_size <= 0L) return(NULL)
    out <- list(
      cost = heap_cost[[1]],
      r = heap_r[[1]],
      c = heap_c[[1]],
      cell = heap_cell[[1]]
    )
    if (heap_size == 1L) {
      heap_size <<- 0L
      return(out)
    }
    heap_cost[[1]] <<- heap_cost[[heap_size]]
    heap_r[[1]] <<- heap_r[[heap_size]]
    heap_c[[1]] <<- heap_c[[heap_size]]
    heap_cell[[1]] <<- heap_cell[[heap_size]]
    heap_size <<- heap_size - 1L

    i <- 1L
    repeat {
      l <- i * 2L
      r <- l + 1L
      if (l > heap_size) break
      m <- l
      if (r <= heap_size && heap_less(r, l)) m <- r
      if (!heap_less(m, i)) break
      heap_swap(i, m)
      i <- m
    }
    out
  }

  rc_to_cell <- function(r, c) as.integer((r - 1L) * nc + c)

  reconstruct <- function(state_cell) {
    path <- integer(0)
    habitat <- integer(0)
    cur <- as.integer(state_cell)
    while (is.finite(cur) && cur > 0L) {
      kind <- state_kind[[cur]]
      if (kind == 1L) {
        habitat <- c(habitat, cur)
      } else {
        path <- c(path, cur)
      }
      cur <- parent_cell[[cur]]
    }
    list(
      path_cells = unique(as.integer(path)),
      habitat_cells = unique(as.integer(habitat))
    )
  }

  for (i in seq_len(nrow(start_positions_rc))) {
    r <- as.integer(start_positions_rc[i, 1])
    c <- as.integer(start_positions_rc[i, 2])
    if (r < 1L || r > nr || c < 1L || c > nc) next
    if (isTRUE(obstacle_block_mat[r, c])) next
    if (!isTRUE(passable_mat[r, c])) next
    cell <- rc_to_cell(r, c)
    if (0 < best_cost[[cell]]) {
      best_cost[[cell]] <- 0
      parent_cell[[cell]] <- 0L
      state_kind[[cell]] <- 0L
      heap_push(0, r, c, cell)
    }
  }
  if (heap_size <= 0L) return(NULL)

  tr <- NULL
  tc <- NULL
  rad <- NULL
  if (!is.null(target_boundary_center) && length(target_boundary_center) >= 2 && !is.null(target_boundary_radius)) {
    tr <- as.integer(target_boundary_center[[1]])
    tc <- as.integer(target_boundary_center[[2]])
    rad <- as.integer(target_boundary_radius)
  }

  start_patch_id <- as.integer(start_patch_id)
  target_patch_id <- as.integer(target_patch_id)
  while (heap_size > 0L) {
    node <- heap_pop()
    if (is.null(node)) break
    cost <- as.numeric(node$cost)
    r <- as.integer(node$r)
    c <- as.integer(node$c)
    cell <- as.integer(node$cell)
    if (cost > best_cost[[cell]]) next
    if (cost > max_width) next

    for (k in seq_len(nrow(moves))) {
      rr <- r + moves[k, 1]
      cc <- c + moves[k, 2]
      if (rr < 1L || rr > nr || cc < 1L || cc > nc) next
      lbl <- labels_mat[rr, cc]
      if (is.finite(lbl) && lbl > 0 && as.integer(lbl) != start_patch_id) {
        if (as.integer(lbl) == target_patch_id) {
          if (!is.null(tr) && !is.null(tc) && !is.null(rad)) {
            if (abs(rr - tr) > rad || abs(cc - tc) > rad) next
          }
          recon <- reconstruct(cell)
          return(list(
            path_cells = recon$path_cells,
            habitat_cells = recon$habitat_cells,
            cost = as.numeric(cost)
          ))
        }
        next
      }

      if (isTRUE(obstacle_block_mat[rr, cc])) next
      if (isTRUE(habitat_mat[rr, cc])) {
        if (is.finite(lbl) && as.integer(lbl) == 0L) {
          move_cost <- 0
          include_kind <- 1L
        } else {
          next
        }
      } else {
        if (!isTRUE(passable_mat[rr, cc])) next
        move_cost <- if (moves[k, 1] != 0L && moves[k, 2] != 0L) sqrt2 else 1
        include_kind <- 0L
      }
      new_cost <- cost + as.numeric(move_cost)
      if (new_cost > max_width) next

      next_cell <- rc_to_cell(rr, cc)
      if (best_cost[[next_cell]] <= new_cost) next
      best_cost[[next_cell]] <- new_cost
      parent_cell[[next_cell]] <- cell
      state_kind[[next_cell]] <- include_kind
      heap_push(new_cost, rr, cc, next_cell)
    }
  }

  NULL
}

terralink_build_routing_graph <- function(
  labels_mat,
  habitat_mat,
  obstacle_block_mat,
  passable_mat,
  connectivity = 4L
) {
  nr <- nrow(labels_mat)
  nc <- ncol(labels_mat)
  bridge_mat <- habitat_mat & (is.na(labels_mat) | labels_mat <= 0)
  traversable_mat <- (!obstacle_block_mat) & (bridge_mat | passable_mat)
  rc <- which(traversable_mat, arr.ind = TRUE)
  ncell <- nr * nc
  if (is.null(rc) || nrow(rc) == 0) {
    return(list(
      graph = igraph::make_empty_graph(n = 0, directed = TRUE),
      cell_to_vid = integer(ncell),
      vid_to_cell = integer(0),
      bridge_cells = rep(FALSE, ncell)
    ))
  }

  cells <- terralink_rc_to_cells(nc, rc)
  ord <- order(cells)
  rc <- rc[ord, , drop = FALSE]
  cells <- as.integer(cells[ord])
  cell_to_vid <- integer(ncell)
  cell_to_vid[cells] <- seq_along(cells)
  bridge_cells <- rep(FALSE, ncell)
  bridge_cells[terralink_rc_to_cells(nc, which(bridge_mat, arr.ind = TRUE))] <- TRUE

  moves <- terralink_neighbor_moves(connectivity)
  sqrt2 <- sqrt(2)
  edges_from <- integer(0)
  edges_to <- integer(0)
  weights <- numeric(0)
  for (idx in seq_along(cells)) {
    r <- rc[idx, 1]
    c <- rc[idx, 2]
    for (k in seq_len(nrow(moves))) {
      rr <- r + moves[k, 1]
      cc <- c + moves[k, 2]
      if (rr < 1L || rr > nr || cc < 1L || cc > nc) next
      cell2 <- as.integer((rr - 1L) * nc + cc)
      vid2 <- cell_to_vid[[cell2]]
      if (vid2 <= 0) next
      is_diag <- (moves[k, 1] != 0L && moves[k, 2] != 0L)
      step <- if (isTRUE(bridge_mat[rr, cc])) 0 else if (is_diag) sqrt2 else 1
      edges_from <- c(edges_from, as.integer(idx))
      edges_to <- c(edges_to, as.integer(vid2))
      weights <- c(weights, as.numeric(step))
    }
  }

  if (length(edges_from) == 0) {
    g <- igraph::make_empty_graph(n = length(cells), directed = TRUE)
  } else {
    g <- igraph::graph_from_edgelist(cbind(edges_from, edges_to), directed = TRUE)
    igraph::E(g)$weight <- weights
  }

  list(
    graph = g,
    cell_to_vid = cell_to_vid,
    vid_to_cell = cells,
    bridge_cells = bridge_cells
  )
}

terralink_route_between_sets <- function(routing, start_cells, target_cells, max_cost = Inf) {
  if (length(start_cells) == 0 || length(target_cells) == 0) return(NULL)
  if (igraph::vcount(routing$graph) == 0) return(NULL)

  start_cells <- unique(as.integer(start_cells))
  target_cells <- unique(as.integer(target_cells))
  start_cells <- start_cells[start_cells >= 1L & start_cells <= length(routing$cell_to_vid)]
  target_cells <- target_cells[target_cells >= 1L & target_cells <= length(routing$cell_to_vid)]
  if (length(start_cells) == 0 || length(target_cells) == 0) return(NULL)

  start_vids <- unique(routing$cell_to_vid[start_cells])
  target_vids <- unique(routing$cell_to_vid[target_cells])
  start_vids <- start_vids[start_vids > 0]
  target_vids <- target_vids[target_vids > 0]
  if (length(start_vids) == 0 || length(target_vids) == 0) return(NULL)

  dmat <- suppressWarnings(
    igraph::distances(
      routing$graph,
      v = start_vids,
      to = target_vids,
      mode = "out",
      weights = igraph::E(routing$graph)$weight
    )
  )
  if (!is.matrix(dmat)) {
    dmat <- matrix(dmat, nrow = length(start_vids), ncol = length(target_vids))
  }
  finite_idx <- which(is.finite(dmat), arr.ind = TRUE)
  if (is.null(finite_idx) || nrow(finite_idx) == 0) return(NULL)
  costs <- dmat[finite_idx]
  best_pos <- which.min(costs)
  min_cost <- as.numeric(costs[[best_pos]])
  if (is.finite(max_cost) && min_cost > as.numeric(max_cost)) return(NULL)

  best_idx <- finite_idx[best_pos, , drop = FALSE]
  best_start <- start_vids[[best_idx[1, 1]]]
  best_target <- target_vids[[best_idx[1, 2]]]
  sp <- suppressWarnings(
    igraph::shortest_paths(
      routing$graph,
      from = best_start,
      to = best_target,
      mode = "out",
      output = "vpath",
      weights = igraph::E(routing$graph)$weight
    )
  )
  if (is.null(sp$vpath) || length(sp$vpath) == 0 || length(sp$vpath[[1]]) == 0) return(NULL)
  path_vids <- as.integer(sp$vpath[[1]])
  path_cells <- as.integer(routing$vid_to_cell[path_vids])
  is_bridge <- routing$bridge_cells[path_cells]

  list(
    path_cells = unique(path_cells[!is_bridge]),
    habitat_cells = unique(path_cells[is_bridge]),
    cost = min_cost
  )
}

#' Build raster candidates using shortest paths
#'
#' @param labels SpatRaster of filtered patch labels.
#' @param patch_df Patch summary data frame.
#' @param passable_mask SpatRaster with 1 for passable cells.
#' @param max_search_distance_px Max search distance in pixels.
#' @param raster_ref Raster reference for CRS/resolution.
#' @param min_corridor_width_px Corridor width in pixels (used for area-based candidate cost).
#' @param pair_index Optional two-column matrix of patch index pairs to evaluate.
#' @param patch_connectivity Patch connectivity (4 or 8).
#' @param habitat_mask SpatRaster of habitat (pre-filter), optional.
#' @param obstacle_mask SpatRaster of blocked pixels, optional.
#' @return Data frame of candidates with geometry.
build_raster_candidates <- function(
  labels,
  patch_df,
  passable_mask,
  max_search_distance_px,
  raster_ref,
  min_corridor_width_px = 1,
  pair_index = NULL,
  patch_connectivity = 4,
  habitat_mask = NULL,
  obstacle_mask = NULL
) {
  if (nrow(patch_df) < 2) {
    return(data.frame())
  }

  labels_mat <- tryCatch(terra::as.matrix(labels, wide = TRUE), error = function(e) NULL)
  if (is.null(labels_mat)) {
    return(data.frame())
  }

  patch_ids <- as.integer(sort(unique(patch_df$patch_id)))
  max_search_distance_px <- as.numeric(max_search_distance_px)
  seeds_by_pair <- terralink_generate_boundary_pair_seeds(
    labels_mat = labels_mat,
    patch_ids = patch_ids,
    connectivity = patch_connectivity,
    max_search_distance = max_search_distance_px,
    min_corridor_width = min_corridor_width_px,
    max_seeds_per_pair = 8L
  )

  allowed_pairs <- NULL
  if (!is.null(pair_index) && nrow(pair_index) > 0) {
    allowed_pairs <- new.env(parent = emptyenv())
    for (row in seq_len(nrow(pair_index))) {
      i <- pair_index[row, 1]
      j <- pair_index[row, 2]
      if (!is.finite(i) || !is.finite(j)) next
      if (i < 1 || i > nrow(patch_df) || j < 1 || j > nrow(patch_df)) next
      key <- terralink_pair_key(patch_df$patch_id[[i]], patch_df$patch_id[[j]])
      allowed_pairs[[key]] <- TRUE
    }
  }

  if (length(seeds_by_pair) == 0) return(data.frame())

  nr <- terra::nrow(raster_ref)
  nc <- terra::ncol(raster_ref)
  pass_vals <- terra::values(passable_mask, mat = FALSE)
  if (is.null(pass_vals) || length(pass_vals) != terra::ncell(raster_ref)) {
    pass_vals <- rep(0, terra::ncell(raster_ref))
  }
  passable_ok <- !is.na(pass_vals) & pass_vals > 0
  passable_mat <- matrix(passable_ok, nrow = nr, ncol = nc, byrow = TRUE)

  obstacle_vals <- NULL
  if (!is.null(obstacle_mask)) {
    obstacle_vals <- terra::values(obstacle_mask, mat = FALSE)
  }
  if (is.null(obstacle_vals) || length(obstacle_vals) != terra::ncell(raster_ref)) {
    obstacle_vals <- rep(0, terra::ncell(raster_ref))
  }
  obstacle_block <- !is.na(obstacle_vals) & obstacle_vals > 0
  obstacle_block_mat <- matrix(obstacle_block, nrow = nr, ncol = nc, byrow = TRUE)
  inflation_ok <- !obstacle_block

  habitat_mat <- NULL
  if (!is.null(habitat_mask)) {
    habitat_raw <- tryCatch(terra::as.matrix(habitat_mask, wide = TRUE), error = function(e) NULL)
    if (!is.null(habitat_raw)) {
      habitat_mat <- !is.na(habitat_raw) & habitat_raw > 0
    }
  }
  if (is.null(habitat_mat)) {
    habitat_mat <- !is.na(labels_mat) & labels_mat > 0
  }
  labels_num <- labels_mat
  labels_num[is.na(labels_num)] <- 0

  start_positions_by_patch <- terralink_compute_start_positions_by_patch(
    labels_mat = labels_num,
    habitat_mat = habitat_mat,
    obstacle_block_mat = obstacle_block_mat,
    passable_mat = passable_mat,
    connectivity = patch_connectivity
  )
  if (length(start_positions_by_patch) == 0) return(data.frame())

  corridor_offsets <- terralink_corridor_offsets(min_corridor_width_px)
  prox_radius <- as.integer(max(1L, ceiling(as.numeric(min_corridor_width_px) / 2)))
  pair_keep <- 4L
  overlap_keep_ratio <- 0.75
  pixel_size <- max(terra::res(raster_ref))
  size_by_patch <- stats::setNames(as.numeric(patch_df$cell_count), as.character(patch_df$patch_id))
  start_filter_radius <- as.integer(max(6L, ceiling(as.numeric(min_corridor_width_px) * 2.0)))
  target_filter_radius <- as.integer(max(6L, ceiling(as.numeric(min_corridor_width_px) * 2.0)))

  records <- list()
  rec_id <- 1L
  for (pair_rec in seeds_by_pair) {
    p1 <- as.integer(pair_rec$patch1)
    p2 <- as.integer(pair_rec$patch2)
    pkey <- terralink_pair_key(p1, p2)
    if (!is.null(allowed_pairs) && !exists(pkey, envir = allowed_pairs, inherits = FALSE)) next
    seeds <- pair_rec$seeds
    if (is.null(seeds) || nrow(seeds) == 0) next
    starts_all <- start_positions_by_patch[[as.character(p1)]]
    if (is.null(starts_all) || nrow(starts_all) == 0) next

    kept <- 0L
    prior_expanded <- list()
    for (si in seq_len(nrow(seeds))) {
      if (kept >= pair_keep) break
      seed_a <- c(as.integer(seeds$r1[[si]]), as.integer(seeds$c1[[si]]))
      seed_b <- c(as.integer(seeds$r2[[si]]), as.integer(seeds$c2[[si]]))
      starts_filtered <- terralink_filter_start_positions_near(starts_all, seed_a, start_filter_radius)
      if (nrow(starts_filtered) == 0) starts_filtered <- starts_all

      route <- terralink_find_shortest_corridor_targeted(
        start_patch_id = p1,
        target_patch_id = p2,
        start_positions_rc = starts_filtered,
        labels_mat = labels_num,
        habitat_mat = habitat_mat,
        max_width = max_search_distance_px,
        connectivity = patch_connectivity,
        obstacle_block_mat = obstacle_block_mat,
        passable_mat = passable_mat,
        target_boundary_center = seed_b,
        target_boundary_radius = target_filter_radius
      )
      if (is.null(route) && si == 1L && kept == 0L) {
        route <- terralink_find_shortest_corridor_targeted(
          start_patch_id = p1,
          target_patch_id = p2,
          start_positions_rc = starts_filtered,
          labels_mat = labels_num,
          habitat_mat = habitat_mat,
          max_width = max_search_distance_px,
          connectivity = patch_connectivity,
          obstacle_block_mat = obstacle_block_mat,
          passable_mat = passable_mat
        )
      }
      if (is.null(route)) next

      path_cells <- as.integer(route$path_cells)
      habitat_touched <- as.integer(route$habitat_cells)
      if (length(path_cells) == 0) next
      buffered <- terralink_inflate_cells(path_cells, nrow = nr, ncol = nc, offsets = corridor_offsets, passable_mask = inflation_ok)
      if (length(buffered) == 0) next
      expanded <- terralink_expand_cells(buffered, nrow = nr, ncol = nc, radius = prox_radius)
      overlap_ratio <- 0
      if (length(prior_expanded) > 0) {
        denom <- max(1, length(buffered))
        for (prev in prior_expanded) {
          overlap_ratio <- max(overlap_ratio, length(intersect(buffered, prev)) / denom)
        }
      }
      if (overlap_ratio >= overlap_keep_ratio) next
      prior_expanded[[length(prior_expanded) + 1L]] <- expanded
      kept <- kept + 1L

      line <- terralink_line_from_cells(raster_ref, path_cells)
      length_map <- as.numeric(route$cost) * pixel_size
      if (!is.finite(length_map) || length_map <= 0) length_map <- as.numeric(seeds$d[[si]]) * pixel_size
      area_px <- as.numeric(length(buffered))
      s1 <- size_by_patch[[as.character(p1)]]
      s2 <- size_by_patch[[as.character(p2)]]
      roi <- sqrt(max(0, s1) * max(0, s2)) / max(area_px, 1e-6)
      records[[length(records) + 1L]] <- list(
        patch1 = p1,
        patch2 = p2,
        cost = area_px,
        distance_map = as.numeric(length_map),
        length = as.numeric(length_map),
        id = rec_id,
        roi = as.numeric(roi),
        variant = paste0("seed_", si),
        line = line,
        path_cells = as.integer(path_cells),
        buffered_cells = as.integer(buffered),
        buffered_cells_expanded = as.integer(expanded),
        habitat_cells = as.integer(habitat_touched)
      )
      rec_id <- rec_id + 1L
    }
  }

  if (length(records) == 0) {
    return(data.frame())
  }

  out <- data.frame(
    patch1 = vapply(records, function(r) r$patch1, integer(1)),
    patch2 = vapply(records, function(r) r$patch2, integer(1)),
    cost = vapply(records, function(r) r$cost, numeric(1)),
    distance_map = vapply(records, function(r) r$distance_map, numeric(1)),
    length = vapply(records, function(r) r$length, numeric(1)),
    id = vapply(records, function(r) r$id, integer(1)),
    roi = vapply(records, function(r) r$roi, numeric(1)),
    variant = vapply(records, function(r) r$variant, character(1)),
    stringsAsFactors = FALSE
  )
  out$line <- sf::st_sfc(lapply(records, function(r) r$line[[1]]), crs = terra::crs(raster_ref))
  out$path_cells <- I(lapply(records, function(r) r$path_cells))
  out$buffered_cells <- I(lapply(records, function(r) r$buffered_cells))
  out$buffered_cells_expanded <- I(lapply(records, function(r) r$buffered_cells_expanded))
  out$habitat_cells <- I(lapply(records, function(r) r$habitat_cells))
  rownames(out) <- NULL
  out
}

#' Build contiguous network raster (patches + corridors)
#'
#' @param habitat_mask SpatRaster of habitat.
#' @param corridor_raster SpatRaster of corridors.
#' @param connectivity Connectivity for patches (4 or 8).
#' @return SpatRaster with component sizes.
build_contiguous_raster <- function(habitat_mask, corridor_raster, connectivity = 8) {
  if (is.null(corridor_raster)) {
    return(NULL)
  }
  combined <- terra::ifel(habitat_mask == 1 | corridor_raster > 0, 1, NA)
  labels <- terra::patches(combined, directions = ifelse(connectivity == 8, 8, 4))
  freq <- terra::freq(labels)
  freq <- freq[!is.na(freq$value) & freq$value != 0, , drop = FALSE]
  if (nrow(freq) == 0) {
    return(labels)
  }
  values_map <- stats::setNames(freq$count, freq$value)
  vals <- terra::values(labels)
  new_vals <- values_map[as.character(vals)]
  out <- labels
  terra::values(out) <- as.numeric(new_vals)
  out
}
