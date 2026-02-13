# Candidate pairing utilities

terralink_pairs_within_distance <- function(x, y, max_distance) {
  n <- length(x)
  if (n < 2) return(matrix(integer(0), ncol = 2))
  if (!is.finite(max_distance) || max_distance <= 0) return(matrix(integer(0), ncol = 2))

  bin_size <- max_distance
  bx <- floor(x / bin_size)
  by <- floor(y / bin_size)
  bins <- split(seq_len(n), paste(bx, by, sep = ":"))

  pairs <- vector("list", n)
  k <- 1
  for (i in seq_len(n)) {
    key <- paste(bx[i], by[i], sep = ":")
    neighbor_keys <- c(
      key,
      paste(bx[i] + 1, by[i], sep = ":"),
      paste(bx[i] - 1, by[i], sep = ":"),
      paste(bx[i], by[i] + 1, sep = ":"),
      paste(bx[i], by[i] - 1, sep = ":"),
      paste(bx[i] + 1, by[i] + 1, sep = ":"),
      paste(bx[i] + 1, by[i] - 1, sep = ":"),
      paste(bx[i] - 1, by[i] + 1, sep = ":"),
      paste(bx[i] - 1, by[i] - 1, sep = ":")
    )
    candidates <- unique(unlist(bins[neighbor_keys], use.names = FALSE))
    if (length(candidates) == 0) next
    candidates <- candidates[candidates > i]
    if (length(candidates) == 0) next
    dx <- x[candidates] - x[i]
    dy <- y[candidates] - y[i]
    dist <- sqrt(dx * dx + dy * dy)
    keep <- which(dist <= max_distance)
    if (length(keep) == 0) next
    idx <- candidates[keep]
    pairs[[k]] <- cbind(rep.int(i, length(idx)), idx)
    k <- k + 1
  }

  if (k == 1) return(matrix(integer(0), ncol = 2))
  do.call(rbind, pairs[seq_len(k - 1)])
}

terralink_vector_pair_index <- function(patches_sf, max_distance, x = NULL, y = NULL) {
  n <- nrow(patches_sf)
  if (n < 2) return(matrix(integer(0), ncol = 2))
  if (!is.finite(max_distance) || max_distance <= 0) return(matrix(integer(0), ncol = 2))

  geom <- tryCatch(sf::st_geometry(patches_sf), error = function(e) NULL)
  if (!is.null(geom) && length(geom) == n) {
    idx <- tryCatch(
      sf::st_is_within_distance(geom, geom, dist = as.numeric(max_distance), sparse = TRUE),
      error = function(e) NULL
    )
    if (!is.null(idx) && length(idx) == n) {
      pairs <- vector("list", n)
      k <- 1L
      for (i in seq_len(n - 1L)) {
        j <- as.integer(idx[[i]])
        j <- j[is.finite(j) & j > i & j <= n]
        if (length(j) == 0) next
        pairs[[k]] <- cbind(rep.int(i, length(j)), j)
        k <- k + 1L
      }
      if (k == 1L) return(matrix(integer(0), ncol = 2))
      return(do.call(rbind, pairs[seq_len(k - 1L)]))
    }
  }

  if (!is.null(x) && !is.null(y) && length(x) == n && length(y) == n) {
    return(terralink_pairs_within_distance(x = x, y = y, max_distance = max_distance))
  }

  out <- utils::combn(seq_len(n), 2)
  t(out)
}
