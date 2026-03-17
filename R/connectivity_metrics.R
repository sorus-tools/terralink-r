HABITAT_AVAILABILITY_DEFAULT_KERNEL <- "exponential"
HABITAT_AVAILABILITY_DEFAULT_SCALING <- "sqrt"
HABITAT_AVAILABILITY_CUTOFF_MULTIPLIER <- 3.0
HABITAT_AVAILABILITY_REDUNDANT_PROB_THRESHOLD <- 0.7
HABITAT_AVAILABILITY_REDUNDANT_GAIN_FACTOR <- 0.25

terralink_normalize_habitat_availability_kernel <- function(value) {
  key <- tolower(trimws(as.character(value %||% HABITAT_AVAILABILITY_DEFAULT_KERNEL)))
  key <- gsub("[[:space:]-]+", "_", key)
  if (!key %in% c("exponential")) {
    key <- HABITAT_AVAILABILITY_DEFAULT_KERNEL
  }
  key
}

terralink_normalize_patch_area_scaling <- function(value) {
  key <- tolower(trimws(as.character(value %||% HABITAT_AVAILABILITY_DEFAULT_SCALING)))
  key <- gsub("[[:space:]-]+", "_", key)
  aliases <- c(
    square_root = "sqrt",
    root = "sqrt",
    logarithmic = "log",
    ln = "log"
  )
  if (key %in% names(aliases)) {
    key <- aliases[[key]]
  }
  if (!key %in% c("sqrt", "log")) {
    key <- HABITAT_AVAILABILITY_DEFAULT_SCALING
  }
  key
}

terralink_habitat_availability_alpha <- function(dispersal_distance) {
  max(as.numeric(dispersal_distance %||% 0) / 3, 1e-9)
}

terralink_scale_patch_area <- function(area, quality_weight = 1, scaling = HABITAT_AVAILABILITY_DEFAULT_SCALING) {
  raw <- max(0, as.numeric(area %||% 0)) * max(0, as.numeric(quality_weight %||% 0))
  mode <- terralink_normalize_patch_area_scaling(scaling)
  if (raw <= 0) {
    return(0)
  }
  if (identical(mode, "log")) {
    return(log1p(raw))
  }
  sqrt(raw)
}

terralink_habitat_probability <- function(distance, alpha, kernel = HABITAT_AVAILABILITY_DEFAULT_KERNEL, cutoff = NULL) {
  d <- as.numeric(distance %||% Inf)
  if (!is.finite(d)) {
    return(0)
  }
  d <- max(d, 0)
  if (!is.null(cutoff) && is.finite(as.numeric(cutoff)) && d > as.numeric(cutoff)) {
    return(0)
  }
  kernel <- terralink_normalize_habitat_availability_kernel(kernel)
  if (identical(kernel, "exponential")) {
    return(exp(-d / max(as.numeric(alpha %||% 0), 1e-9)))
  }
  exp(-d / max(as.numeric(alpha %||% 0), 1e-9))
}

terralink_matrix_changed <- function(a, b, tol = 1e-12) {
  same_inf <- is.infinite(a) & is.infinite(b)
  same_val <- suppressWarnings(abs(a - b) <= tol)
  same_val[is.na(same_val)] <- FALSE
  !all(same_inf | same_val)
}

terralink_candidate_patch_ids <- function(cand, valid_node_ids = NULL) {
  valid <- NULL
  if (!is.null(valid_node_ids)) {
    valid <- unique(as.integer(valid_node_ids[is.finite(valid_node_ids)]))
  }
  ids <- integer(0)
  if ("patch_ids" %in% names(cand)) {
    raw <- cand$patch_ids
    if (is.list(raw) && length(raw) == 1) {
      raw <- raw[[1]]
    }
    add <- suppressWarnings(as.integer(raw))
    add <- add[is.finite(add)]
    ids <- c(ids, add)
  }
  if (length(ids) < 2) {
    for (nm in c("patch1", "p1", "patch2", "p2")) {
      if (!nm %in% names(cand)) next
      add <- suppressWarnings(as.integer(cand[[nm]][[1]] %||% cand[[nm]]))
      if (is.finite(add)) ids <- c(ids, add)
    }
  }
  ids <- unique(ids)
  if (!is.null(valid)) {
    ids <- ids[ids %in% valid]
  }
  ids
}

terralink_bbox_intersects <- function(bbox_a, bbox_b) {
  !(as.numeric(bbox_a[["xmax"]]) < as.numeric(bbox_b[["xmin"]]) ||
      as.numeric(bbox_b[["xmax"]]) < as.numeric(bbox_a[["xmin"]]) ||
      as.numeric(bbox_a[["ymax"]]) < as.numeric(bbox_b[["ymin"]]) ||
      as.numeric(bbox_b[["ymax"]]) < as.numeric(bbox_a[["ymin"]]))
}

terralink_vector_build_resistance_context <- function(
  patches_sf,
  patch_df,
  pair_samples = 300L
) {
  if (!inherits(patches_sf, "sf") || !is.data.frame(patch_df) || nrow(patch_df) < 2) {
    return(NULL)
  }

  node_ids <- suppressWarnings(as.integer(patch_df$patch_id))
  node_ids <- node_ids[is.finite(node_ids)]
  if (length(node_ids) < 2) {
    return(NULL)
  }

  geom <- sf::st_geometry(patches_sf)
  coords_xy <- stats::setNames(
    lapply(seq_len(nrow(patch_df)), function(i) c(as.numeric(patch_df$x[[i]]), as.numeric(patch_df$y[[i]]))),
    as.character(node_ids)
  )
  area_ha_by_id <- stats::setNames(as.numeric(patch_df$area_m2) / 10000, as.character(node_ids))
  bbox_list <- lapply(seq_along(geom), function(i) sf::st_bbox(geom[[i]]))

  graph <- igraph::make_empty_graph(n = length(node_ids), directed = FALSE)
  igraph::V(graph)$name <- as.character(node_ids)

  if (length(node_ids) >= 2) {
    for (i in seq_len(length(node_ids) - 1L)) {
      geom_i <- geom[[i]]
      if (is.null(geom_i)) next
      for (j in (i + 1L):length(node_ids)) {
        geom_j <- geom[[j]]
        if (is.null(geom_j) || !terralink_bbox_intersects(bbox_list[[i]], bbox_list[[j]])) next
        d_poly <- suppressWarnings(tryCatch(as.numeric(sf::st_distance(geom_i, geom_j)), error = function(e) Inf))
        if (!is.finite(d_poly) || d_poly > 1e-6) next
        xy_i <- coords_xy[[as.character(node_ids[[i]])]]
        xy_j <- coords_xy[[as.character(node_ids[[j]])]]
        w <- max(sqrt((xy_i[[1]] - xy_j[[1]])^2 + (xy_i[[2]] - xy_j[[2]])^2), 1e-6)
        graph <- terralink_graph_add_or_update(graph, c(node_ids[[i]], node_ids[[j]]), w)
      }
    }
  }

  pair_records <- data.frame(
    u = integer(0),
    v = integer(0),
    euclid = numeric(0),
    stringsAsFactors = FALSE
  )
  if (length(node_ids) >= 2) {
    for (i in seq_len(length(node_ids) - 1L)) {
      xy_i <- coords_xy[[as.character(node_ids[[i]])]]
      for (j in (i + 1L):length(node_ids)) {
        xy_j <- coords_xy[[as.character(node_ids[[j]])]]
        eu <- sqrt((xy_i[[1]] - xy_j[[1]])^2 + (xy_i[[2]] - xy_j[[2]])^2)
        if (!is.finite(eu) || eu <= 1e-9) next
        pair_records <- rbind(
          pair_records,
          data.frame(u = node_ids[[i]], v = node_ids[[j]], euclid = eu, stringsAsFactors = FALSE)
        )
      }
    }
  }
  if (nrow(pair_records) == 0) {
    return(NULL)
  }

  pair_target <- max(200L, min(500L, as.integer(pair_samples %||% 300L)))
  if (nrow(pair_records) > pair_target) {
    pair_records <- pair_records[order(-pair_records$euclid, pair_records$u, pair_records$v), , drop = FALSE]
    pair_records <- pair_records[seq_len(pair_target), , drop = FALSE]
  }

  pair_weights_raw <- sqrt(
    pmax(as.numeric(area_ha_by_id[as.character(pair_records$u)]), 0) *
      pmax(as.numeric(area_ha_by_id[as.character(pair_records$v)]), 0)
  )
  pair_weights_raw[!is.finite(pair_weights_raw) | pair_weights_raw <= 0] <- 1e-9
  pair_weight_sum <- sum(pair_weights_raw)
  if (!is.finite(pair_weight_sum) || pair_weight_sum <= 0) {
    pair_weight_sum <- nrow(pair_records)
  }
  pair_records$weight <- pair_weights_raw / pair_weight_sum

  list(
    graph = graph,
    pair_records = pair_records,
    disconnected_penalty = max(1, max(pair_records$euclid, na.rm = TRUE) * 20),
    total_habitat_area_ha = sum(as.numeric(area_ha_by_id), na.rm = TRUE)
  )
}

terralink_vector_weighted_mean_graph_distance <- function(
  graph,
  pair_records,
  disconnected_penalty
) {
  if (is.null(graph) || !inherits(graph, "igraph") || !is.data.frame(pair_records) || nrow(pair_records) == 0) {
    return(as.numeric(disconnected_penalty %||% 1))
  }
  terralink_matrix_pseudoinverse <- function(mat, tol = 1e-9) {
    sv <- tryCatch(svd(mat), error = function(e) NULL)
    if (is.null(sv)) return(NULL)
    if (length(sv$d) == 0) return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
    cutoff <- max(as.numeric(tol), 1e-12) * max(sv$d)
    d_inv <- ifelse(sv$d > cutoff, 1 / sv$d, 0)
    sv$v %*% (diag(d_inv, nrow = length(d_inv), ncol = length(d_inv)) %*% t(sv$u))
  }

  comp_membership <- tryCatch(igraph::components(graph)$membership, error = function(e) NULL)
  resistance_by_component <- list()
  if (!is.null(comp_membership)) {
    for (comp_id in sort(unique(as.integer(comp_membership)))) {
      vids <- names(comp_membership)[comp_membership == comp_id]
      if (length(vids) < 2) next
      subg <- tryCatch(igraph::induced_subgraph(graph, vids = vids), error = function(e) NULL)
      if (is.null(subg)) next
      adj <- suppressWarnings(tryCatch(
        as.matrix(igraph::as_adjacency_matrix(subg, attr = "weight", sparse = FALSE)),
        error = function(e) NULL
      ))
      if (is.null(adj)) next
      conduct <- matrix(0, nrow = nrow(adj), ncol = ncol(adj), dimnames = dimnames(adj))
      valid <- is.finite(adj) & adj > 0
      conduct[valid] <- 1 / adj[valid]
      lap <- diag(rowSums(conduct), nrow = nrow(conduct), ncol = ncol(conduct)) - conduct
      lap_plus <- terralink_matrix_pseudoinverse(lap)
      if (is.null(lap_plus)) next
      diag_vals <- diag(lap_plus)
      resistance <- outer(diag_vals, diag_vals, "+") - (2 * lap_plus)
      diag(resistance) <- 0
      rownames(resistance) <- igraph::V(subg)$name
      colnames(resistance) <- igraph::V(subg)$name
      resistance_by_component[[as.character(comp_id)]] <- resistance
    }
  }

  dist_mat <- suppressWarnings(tryCatch(
    igraph::distances(graph, v = igraph::V(graph), to = igraph::V(graph), weights = igraph::E(graph)$weight),
    error = function(e) NULL
  ))
  if (!is.null(dist_mat)) {
    rownames(dist_mat) <- igraph::V(graph)$name
    colnames(dist_mat) <- igraph::V(graph)$name
  }
  vals <- vapply(seq_len(nrow(pair_records)), function(i) {
    u <- as.character(pair_records$u[[i]])
    v <- as.character(pair_records$v[[i]])
    comp_id <- if (!is.null(comp_membership)) as.character(comp_membership[[u]] %||% NA_integer_) else NA_character_
    if (!is.na(comp_id) && identical(comp_id, as.character(comp_membership[[v]] %||% NA_integer_))) {
      resistance <- resistance_by_component[[comp_id]]
      if (!is.null(resistance)) {
        d_eff <- suppressWarnings(as.numeric(resistance[u, v]))
        if (is.finite(d_eff) && d_eff >= 0) {
          return(d_eff)
        }
      }
    }
    if (!is.null(dist_mat)) {
      d <- suppressWarnings(as.numeric(dist_mat[u, v]))
      if (is.finite(d)) return(max(d, 0))
    }
    as.numeric(disconnected_penalty %||% 1)
  }, numeric(1))
  weights <- as.numeric(pair_records$weight %||% rep(1 / max(nrow(pair_records), 1), nrow(pair_records)))
  weights[!is.finite(weights) | weights < 0] <- 0
  if (sum(weights) <= 0) {
    weights <- rep(1 / max(length(vals), 1), length(vals))
  }
  sum(vals * weights)
}

terralink_vector_exact_metric_overrides <- function(
  patches_sf,
  patch_df,
  selected_corridors,
  total_connected_area_post,
  largest_network_area_post,
  area_div
) {
  overrides <- list(
    total_connected_habitat_area_post = as.numeric(total_connected_area_post %||% 0) / max(as.numeric(area_div %||% 1), 1e-9),
    largest_network_area_post = as.numeric(largest_network_area_post %||% 0) / max(as.numeric(area_div %||% 1), 1e-9)
  )
  ctx <- terralink_vector_build_resistance_context(patches_sf = patches_sf, patch_df = patch_df)
  if (is.null(ctx)) {
    return(overrides)
  }

  graph_post <- ctx$graph
  if (is.data.frame(selected_corridors) && nrow(selected_corridors) > 0) {
    for (i in seq_len(nrow(selected_corridors))) {
      ids <- terralink_candidate_patch_ids(selected_corridors[i, , drop = FALSE], valid_node_ids = patch_df$patch_id)
      if (length(ids) < 2) next
      len_val <- as.numeric(
        selected_corridors$distance_m[[i]] %||%
          selected_corridors$length[[i]] %||%
          selected_corridors$length_metric[[i]] %||%
          0
      )
      if (!is.finite(len_val) || len_val <= 0) next
      graph_post <- terralink_graph_add_or_update(graph_post, ids, len_val)
    }
  }

  mean_pre <- terralink_vector_weighted_mean_graph_distance(
    graph = ctx$graph,
    pair_records = ctx$pair_records,
    disconnected_penalty = ctx$disconnected_penalty
  )
  mean_post <- terralink_vector_weighted_mean_graph_distance(
    graph = graph_post,
    pair_records = ctx$pair_records,
    disconnected_penalty = ctx$disconnected_penalty
  )
  total_habitat_area_ha <- max(as.numeric(ctx$total_habitat_area_ha %||% 0), 0)
  overrides$mean_effective_resistance_pre <- mean_pre
  overrides$mean_effective_resistance_post <- mean_post
  overrides$landscape_fluidity_pre <- if (mean_pre > 0) total_habitat_area_ha / mean_pre else 0
  overrides$landscape_fluidity_post <- if (mean_post > 0) total_habitat_area_ha / mean_post else 0
  overrides
}

terralink_ha_distance_to_probability_matrix <- function(dist, alpha, cutoff) {
  out <- matrix(0, nrow = nrow(dist), ncol = ncol(dist))
  valid <- is.finite(dist) & (dist <= cutoff + 1e-12)
  if (any(valid)) {
    out[valid] <- exp(-pmax(dist[valid], 0) / alpha)
  }
  if (nrow(out) == ncol(out) && nrow(out) > 0) {
    diag(out) <- 1
  }
  out
}

terralink_ha_cluster_area_within_threshold <- function(state, dist = state$dist) {
  n <- length(state$node_ids)
  if (n <= 0) {
    return(0)
  }
  parent <- seq_len(n)
  find_root <- function(idx) {
    root <- idx
    while (parent[[root]] != root) {
      root <- parent[[root]]
    }
    while (parent[[idx]] != idx) {
      nxt <- parent[[idx]]
      parent[[idx]] <<- root
      idx <- nxt
    }
    root
  }
  union_root <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) {
      parent[[rb]] <<- ra
    }
  }
  thr <- max(as.numeric(state$dispersal_distance %||% 0), 0)
  if (n >= 2) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        if (is.finite(dist[i, j]) && dist[i, j] <= thr + 1e-12) {
          union_root(i, j)
        }
      }
    }
  }
  cluster_area <- list()
  for (idx in seq_len(n)) {
    root <- as.character(find_root(idx))
    cluster_area[[root]] <- as.numeric(cluster_area[[root]] %||% 0) + max(as.numeric(state$raw_areas[[idx]] %||% 0), 0)
  }
  vals <- unlist(cluster_area, use.names = FALSE)
  if (length(vals) == 0) return(0)
  max(vals, 0, na.rm = TRUE)
}

terralink_ha_metrics_from_state <- function(
  state,
  dist = state$dist,
  prob = state$prob,
  reachable_area = state$reachable_area,
  habitat_availability = state$habitat_availability
) {
  reachable <- as.numeric(reachable_area %||% numeric(0))
  list(
    habitat_availability = as.numeric(habitat_availability %||% 0),
    mean_reachable_area = if (length(reachable) > 0) mean(reachable) else 0,
    median_reachable_area = if (length(reachable) > 0) stats::median(reachable) else 0,
    largest_reachable_habitat_cluster = terralink_ha_cluster_area_within_threshold(state, dist = dist)
  )
}

terralink_new_habitat_availability_state <- function(
  nodes,
  dispersal_distance,
  kernel = HABITAT_AVAILABILITY_DEFAULT_KERNEL,
  cutoff_multiplier = HABITAT_AVAILABILITY_CUTOFF_MULTIPLIER
) {
  if (!is.data.frame(nodes) || nrow(nodes) == 0) {
    return(list(
      node_ids = integer(0),
      index_by_id = integer(0),
      raw_areas = numeric(0),
      weights = numeric(0),
      dispersal_distance = max(as.numeric(dispersal_distance %||% 0), 1e-9),
      alpha = terralink_habitat_availability_alpha(dispersal_distance),
      kernel = terralink_normalize_habitat_availability_kernel(kernel),
      cutoff = max(as.numeric(cutoff_multiplier %||% 0), 1) * max(as.numeric(dispersal_distance %||% 0), 1e-9),
      graph_edges = list(),
      dist = matrix(numeric(0), 0, 0),
      prob = matrix(numeric(0), 0, 0),
      reachable_area = numeric(0),
      habitat_availability = 0
    ))
  }
  keep <- is.finite(nodes$effective_area) & (nodes$effective_area > 0)
  nodes <- nodes[keep, , drop = FALSE]
  if (nrow(nodes) == 0) {
    return(terralink_new_habitat_availability_state(data.frame(), dispersal_distance, kernel, cutoff_multiplier))
  }
  nodes <- nodes[order(as.integer(nodes$node_id)), , drop = FALSE]
  node_ids <- as.integer(nodes$node_id)
  size <- length(node_ids)
  dist <- matrix(Inf, nrow = size, ncol = size)
  diag(dist) <- 0
  alpha <- terralink_habitat_availability_alpha(dispersal_distance)
  cutoff <- max(as.numeric(cutoff_multiplier %||% 0), 1) * max(as.numeric(dispersal_distance %||% 0), 1e-9)
  prob <- terralink_ha_distance_to_probability_matrix(dist, alpha = alpha, cutoff = cutoff)
  weights <- as.numeric(nodes$effective_area)
  reachable_area <- as.numeric(prob %*% weights)
  list(
    node_ids = node_ids,
    index_by_id = stats::setNames(seq_along(node_ids), as.character(node_ids)),
    raw_areas = as.numeric(nodes$raw_area),
    weights = weights,
    dispersal_distance = max(as.numeric(dispersal_distance %||% 0), 1e-9),
    alpha = alpha,
    kernel = terralink_normalize_habitat_availability_kernel(kernel),
    cutoff = cutoff,
    graph_edges = list(),
    dist = dist,
    prob = prob,
    reachable_area = reachable_area,
    habitat_availability = sum(weights * reachable_area)
  )
}

terralink_apply_edge_to_distance_matrix <- function(dist, u_idx, v_idx, length, cutoff = NULL) {
  if (!is.finite(length) || length <= 0 || u_idx == v_idx) {
    return(dist)
  }
  out <- dist
  if (is.finite(out[u_idx, v_idx]) && out[u_idx, v_idx] <= length + 1e-12) {
    return(out)
  }
  out[u_idx, v_idx] <- min(out[u_idx, v_idx], length)
  out[v_idx, u_idx] <- min(out[v_idx, u_idx], length)
  via_uv <- outer(dist[, u_idx], dist[v_idx, ], "+") + length
  via_vu <- outer(dist[, v_idx], dist[u_idx, ], "+") + length
  out <- pmin(out, via_uv, via_vu)
  if (!is.null(cutoff) && is.finite(as.numeric(cutoff))) {
    out[out > as.numeric(cutoff)] <- Inf
  }
  diag(out) <- 0
  out
}

terralink_ha_normalize_candidate_edges <- function(patch_ids, length, valid_node_ids = NULL) {
  ids <- suppressWarnings(as.integer(patch_ids))
  ids <- ids[is.finite(ids)]
  ids <- unique(ids)
  if (!is.null(valid_node_ids)) {
    valid <- unique(as.integer(valid_node_ids[is.finite(valid_node_ids)]))
    ids <- ids[ids %in% valid]
  }
  if (length(ids) < 2) {
    return(matrix(numeric(0), nrow = 0, ncol = 3))
  }
  edge_len <- max(as.numeric(length %||% 0), 1e-9)
  out <- matrix(numeric(0), nrow = 0, ncol = 3)
  if (length(ids) >= 2) {
    for (i in seq_len(length(ids) - 1L)) {
      for (j in (i + 1L):length(ids)) {
        a <- ids[[i]]
        b <- ids[[j]]
        if (a == b) next
        if (a > b) {
          tmp <- a
          a <- b
          b <- tmp
        }
        out <- rbind(out, c(a, b, edge_len))
      }
    }
  }
  out
}

terralink_ha_evaluate_candidate <- function(state, edge_specs, corridor_cost) {
  cost <- max(as.numeric(corridor_cost %||% 0), 0)
  if (cost <= 0 || length(state$node_ids) < 2 || length(edge_specs) == 0) {
    return(NULL)
  }
  if (is.vector(edge_specs)) {
    edge_specs <- matrix(edge_specs, ncol = 3)
  }
  temp_dist <- state$dist
  applied_edges <- matrix(numeric(0), nrow = 0, ncol = 3)
  endpoint_probs_before <- numeric(0)
  connects_new_component <- FALSE
  for (k in seq_len(nrow(edge_specs))) {
    u_id <- as.integer(edge_specs[k, 1])
    v_id <- as.integer(edge_specs[k, 2])
    edge_len <- as.numeric(edge_specs[k, 3])
    u_idx <- state$index_by_id[[as.character(u_id)]]
    v_idx <- state$index_by_id[[as.character(v_id)]]
    if (is.null(u_idx) || is.null(v_idx) || u_idx == v_idx) next
    prev_dist <- as.numeric(temp_dist[u_idx, v_idx])
    if (!is.finite(prev_dist)) {
      connects_new_component <- TRUE
    }
    endpoint_probs_before <- c(
      endpoint_probs_before,
      terralink_habitat_probability(prev_dist, alpha = state$alpha, kernel = state$kernel, cutoff = state$cutoff)
    )
    new_temp <- terralink_apply_edge_to_distance_matrix(temp_dist, u_idx, v_idx, edge_len, cutoff = state$cutoff)
    if (!terralink_matrix_changed(new_temp, temp_dist)) next
    temp_dist <- new_temp
    applied_edges <- rbind(applied_edges, c(u_id, v_id, edge_len))
  }
  if (nrow(applied_edges) == 0) {
    return(NULL)
  }
  new_prob <- terralink_ha_distance_to_probability_matrix(temp_dist, alpha = state$alpha, cutoff = state$cutoff)
  delta_prob <- new_prob - state$prob
  if (all(abs(delta_prob) <= 1e-12)) {
    return(NULL)
  }
  new_reachable <- state$reachable_area + as.numeric(delta_prob %*% state$weights)
  new_ha <- sum(state$weights * new_reachable)
  gain <- as.numeric(new_ha - state$habitat_availability)
  if (!is.finite(gain) || gain <= 1e-12) {
    return(NULL)
  }
  penalized_gain <- gain
  if (length(endpoint_probs_before) > 0 && max(endpoint_probs_before, na.rm = TRUE) > HABITAT_AVAILABILITY_REDUNDANT_PROB_THRESHOLD) {
    penalized_gain <- penalized_gain * HABITAT_AVAILABILITY_REDUNDANT_GAIN_FACTOR
  }
  metrics <- terralink_ha_metrics_from_state(
    state,
    dist = temp_dist,
    prob = new_prob,
    reachable_area = new_reachable,
    habitat_availability = new_ha
  )
  list(
    score = penalized_gain / max(cost, 1e-12),
    gain = gain,
    penalized_gain = penalized_gain,
    cost = cost,
    applied_edges = applied_edges,
    dist = temp_dist,
    prob = new_prob,
    reachable_area = new_reachable,
    metrics = metrics,
    type = if (connects_new_component) "primary" else "redundant"
  )
}

terralink_ha_apply_candidate <- function(state, evaluation) {
  if (is.null(evaluation)) {
    return(state)
  }
  applied <- evaluation$applied_edges
  if (!is.null(applied) && nrow(applied) > 0) {
    for (k in seq_len(nrow(applied))) {
      key <- paste(sort(as.integer(applied[k, 1:2])), collapse = "_")
      prev <- state$graph_edges[[key]]
      edge_len <- as.numeric(applied[k, 3])
      if (is.null(prev) || edge_len + 1e-12 < as.numeric(prev)) {
        state$graph_edges[[key]] <- edge_len
      }
    }
  }
  state$dist <- evaluation$dist
  state$prob <- evaluation$prob
  state$reachable_area <- evaluation$reachable_area
  state$habitat_availability <- as.numeric(evaluation$metrics$habitat_availability %||% state$habitat_availability)
  state
}

terralink_pair_euclid_matrix <- function(node_df) {
  n <- nrow(node_df)
  if (n <= 0) {
    return(matrix(numeric(0), 0, 0))
  }
  x <- as.numeric(node_df$x)
  y <- as.numeric(node_df$y)
  dx <- outer(x, x, "-")
  dy <- outer(y, y, "-")
  sqrt(dx * dx + dy * dy)
}

terralink_ime_norm <- function(dist, pair_euclid) {
  if (length(dist) == 0 || length(pair_euclid) == 0) {
    return(0)
  }
  idx <- upper.tri(pair_euclid)
  eu <- pair_euclid[idx]
  d <- dist[idx]
  if (length(eu) == 0) {
    return(0)
  }
  eps <- 1e-9
  num <- sum(ifelse(is.finite(d), 1 / (d + eps), 0))
  den <- sum(1 / (pmax(eu, 0) + eps))
  if (!is.finite(den) || den <= 0) {
    return(0)
  }
  max(0, min(1, num / den))
}

terralink_strategic_mobility <- function(dist, pair_euclid, detour_cap = 8) {
  idx <- upper.tri(pair_euclid)
  eu <- pair_euclid[idx]
  d <- dist[idx]
  if (length(eu) == 0) {
    return(0)
  }
  cap <- max(as.numeric(detour_cap %||% 0), 1)
  ratios <- rep(cap, length(eu))
  finite_eu <- is.finite(eu) & eu > 0
  finite_d <- is.finite(d)
  ratios[finite_eu & finite_d] <- pmax(1, d[finite_eu & finite_d] / eu[finite_eu & finite_d])
  ratios <- pmin(ratios, cap)
  mean_ratio <- mean(ratios)
  if (!is.finite(mean_ratio) || mean_ratio <= 0) {
    return(0)
  }
  max(0, min(1, 1 / mean_ratio))
}

terralink_fri_norm <- function(dist, pair_euclid) {
  n <- nrow(dist)
  if (n < 2) {
    return(0)
  }
  finite_idx <- which(upper.tri(dist) & is.finite(dist) & dist > 0, arr.ind = TRUE)
  if (nrow(finite_idx) == 0) {
    return(0)
  }
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  igraph::V(g)$name <- as.character(seq_len(n))
  for (k in seq_len(nrow(finite_idx))) {
    i <- finite_idx[k, 1]
    j <- finite_idx[k, 2]
    g <- igraph::add_edges(g, c(as.character(i), as.character(j)), attr = list(weight = as.numeric(dist[i, j])))
  }
  comps <- igraph::components(g)$membership
  if (is.null(comps) || length(comps) != n) {
    return(0)
  }
  resist_vals <- numeric(0)
  eu_vals <- numeric(0)
  for (comp_id in unique(comps)) {
    nodes <- which(comps == comp_id)
    if (length(nodes) < 2) next
    sub <- igraph::induced_subgraph(g, vids = as.character(nodes))
    edf <- igraph::as_data_frame(sub, what = "edges")
    if (nrow(edf) == 0) next
    local_index <- stats::setNames(seq_along(nodes), as.character(nodes))
    n_local <- length(nodes)
    lap <- matrix(0, nrow = n_local, ncol = n_local)
    for (row in seq_len(nrow(edf))) {
      i <- local_index[[edf$from[[row]]]]
      j <- local_index[[edf$to[[row]]]]
      d <- as.numeric(edf$weight[[row]])
      if (!is.finite(d) || d <= 0) next
      cnd <- 1 / max(d, 1e-9)
      lap[i, i] <- lap[i, i] + cnd
      lap[j, j] <- lap[j, j] + cnd
      lap[i, j] <- lap[i, j] - cnd
      lap[j, i] <- lap[j, i] - cnd
    }
    ref <- n_local
    lap_r <- lap[seq_len(ref - 1L), seq_len(ref - 1L), drop = FALSE]
    if (nrow(lap_r) == 0) next
    inv_ok <- TRUE
    for (i in seq_len(n_local - 1L)) {
      for (j in (i + 1L):n_local) {
        b <- rep(0, ref - 1L)
        if (i != ref) b[[i]] <- b[[i]] + 1
        if (j != ref) b[[j]] <- b[[j]] - 1
        x <- tryCatch(solve(lap_r, b), error = function(e) {
          inv_ok <<- FALSE
          NULL
        })
        if (!inv_ok || is.null(x)) break
        vi <- if (i == ref) 0 else x[[i]]
        vj <- if (j == ref) 0 else x[[j]]
        r_eff <- abs(vi - vj)
        if (is.finite(r_eff)) {
          resist_vals <- c(resist_vals, r_eff)
          eu_vals <- c(eu_vals, pair_euclid[nodes[[i]], nodes[[j]]])
        }
      }
      if (!inv_ok) break
    }
    if (!inv_ok) next
  }
  if (length(resist_vals) == 0 || length(eu_vals) == 0) {
    return(0)
  }
  r_bar <- mean(resist_vals)
  eu_bar <- mean(eu_vals)
  if (!is.finite(eu_bar) || eu_bar <= 0) {
    return(0)
  }
  max(0, min(1, 1 / (1 + (r_bar / max(eu_bar, 1e-9)))))
}

terralink_mean_path_resistance <- function(dist, pair_euclid, detour_cap = 20) {
  idx <- upper.tri(pair_euclid)
  eu <- pair_euclid[idx]
  d <- dist[idx]
  if (length(eu) == 0) {
    return(0)
  }
  max_eu <- max(eu, na.rm = TRUE)
  penalty <- max(1, max_eu * max(as.numeric(detour_cap %||% 0), 1))
  vals <- ifelse(is.finite(d), pmax(d, 0), penalty)
  mean(vals)
}

terralink_new_fluidity_state <- function(node_df, detour_cap = 8, shortcut_threshold = 3, redundancy_method = "ime") {
  if (!is.data.frame(node_df) || nrow(node_df) == 0) {
    return(list(
      node_ids = integer(0),
      index_by_id = integer(0),
      pair_euclid = matrix(numeric(0), 0, 0),
      dist = matrix(numeric(0), 0, 0),
      detour_cap = max(as.numeric(detour_cap %||% 0), 1),
      shortcut_threshold = max(as.numeric(shortcut_threshold %||% 0), 1),
      redundancy_method = tolower(as.character(redundancy_method %||% "ime")),
      current_metrics = list(
        mean_effective_resistance = 0,
        landscape_fluidity = 0,
        ime = 0,
        strategic_mobility = 0,
        flow_redundancy = 0,
        flow_method = "IME"
      )
    ))
  }
  node_df <- node_df[order(as.integer(node_df$patch_id)), , drop = FALSE]
  node_ids <- as.integer(node_df$patch_id)
  n <- nrow(node_df)
  dist <- matrix(Inf, nrow = n, ncol = n)
  diag(dist) <- 0
  state <- list(
    node_ids = node_ids,
    index_by_id = stats::setNames(seq_along(node_ids), as.character(node_ids)),
    pair_euclid = terralink_pair_euclid_matrix(data.frame(x = node_df$x, y = node_df$y)),
    dist = dist,
    detour_cap = max(as.numeric(detour_cap %||% 0), 1),
    shortcut_threshold = max(as.numeric(shortcut_threshold %||% 0), 1),
    redundancy_method = tolower(as.character(redundancy_method %||% "ime"))
  )
  state$current_metrics <- terralink_fluidity_metrics_from_dist(state, dist)
  state
}

terralink_fluidity_metrics_from_dist <- function(state, dist) {
  mean_res <- terralink_mean_path_resistance(dist, state$pair_euclid, detour_cap = state$detour_cap)
  fluidity <- if (is.finite(mean_res) && mean_res > 0) 1 / mean_res else 0
  ime <- terralink_ime_norm(dist, state$pair_euclid)
  fri <- terralink_fri_norm(dist, state$pair_euclid)
  mobility <- terralink_strategic_mobility(dist, state$pair_euclid, detour_cap = state$detour_cap)
  use_fri <- identical(state$redundancy_method, "fri") && is.finite(fri) && fri > 0
  flow <- if (use_fri) fri else ime
  list(
    mean_effective_resistance = mean_res,
    landscape_fluidity = fluidity,
    ime = ime,
    fri = fri,
    strategic_mobility = mobility,
    flow_redundancy = flow,
    flow_method = if (use_fri) "FRI" else "IME"
  )
}

terralink_fluidity_evaluate_candidate <- function(state, edge_specs, corridor_cost) {
  cost <- max(as.numeric(corridor_cost %||% 0), 0)
  if (cost <= 0 || length(state$node_ids) < 2 || length(edge_specs) == 0) {
    return(NULL)
  }
  if (is.vector(edge_specs)) {
    edge_specs <- matrix(edge_specs, ncol = 3)
  }
  temp_dist <- state$dist
  applied_edges <- matrix(numeric(0), nrow = 0, ncol = 3)
  endpoint_ratios <- numeric(0)
  connects_new_component <- FALSE
  for (k in seq_len(nrow(edge_specs))) {
    u_id <- as.integer(edge_specs[k, 1])
    v_id <- as.integer(edge_specs[k, 2])
    edge_len <- max(as.numeric(edge_specs[k, 3]), 1e-9)
    u_idx <- state$index_by_id[[as.character(u_id)]]
    v_idx <- state$index_by_id[[as.character(v_id)]]
    if (is.null(u_idx) || is.null(v_idx) || u_idx == v_idx) next
    prev_dist <- as.numeric(temp_dist[u_idx, v_idx])
    if (!is.finite(prev_dist)) {
      connects_new_component <- TRUE
    } else {
      endpoint_ratios <- c(endpoint_ratios, prev_dist / edge_len)
    }
    new_temp <- terralink_apply_edge_to_distance_matrix(temp_dist, u_idx, v_idx, edge_len, cutoff = NULL)
    if (!terralink_matrix_changed(new_temp, temp_dist)) next
    temp_dist <- new_temp
    applied_edges <- rbind(applied_edges, c(u_id, v_id, edge_len))
  }
  if (nrow(applied_edges) == 0) {
    return(NULL)
  }
  metrics_new <- terralink_fluidity_metrics_from_dist(state, temp_dist)
  metrics_old <- state$current_metrics
  gain_fluidity <- max(as.numeric(metrics_new$landscape_fluidity %||% 0) - as.numeric(metrics_old$landscape_fluidity %||% 0), 0)
  gain_flow <- max(as.numeric(metrics_new$flow_redundancy %||% 0) - as.numeric(metrics_old$flow_redundancy %||% 0), 0)
  gain_mobility <- max(as.numeric(metrics_new$strategic_mobility %||% 0) - as.numeric(metrics_old$strategic_mobility %||% 0), 0)
  gain <- gain_fluidity + (0.5 * gain_flow) + (0.25 * gain_mobility)
  if (!is.finite(gain) || gain <= 1e-12) {
    return(NULL)
  }
  penalized_gain <- gain
  if (length(endpoint_ratios) > 0 && max(endpoint_ratios, na.rm = TRUE) < state$shortcut_threshold && !connects_new_component) {
    penalized_gain <- penalized_gain * 0.25
  }
  list(
    score = penalized_gain / max(cost, 1e-12),
    gain = gain,
    penalized_gain = penalized_gain,
    cost = cost,
    applied_edges = applied_edges,
    dist = temp_dist,
    metrics = metrics_new,
    type = if (connects_new_component) "primary" else "redundant"
  )
}

terralink_fluidity_apply_candidate <- function(state, evaluation) {
  if (is.null(evaluation)) {
    return(state)
  }
  state$dist <- evaluation$dist
  state$current_metrics <- evaluation$metrics
  state
}

terralink_component_stats <- function(node_ids, patch_areas, edge_specs) {
  node_ids <- unique(as.integer(node_ids[is.finite(node_ids)]))
  patch_areas <- stats::setNames(as.numeric(patch_areas[as.character(node_ids)]), as.character(node_ids))
  if (length(node_ids) == 0) {
    return(list(connected_area = 0, largest_network = 0, largest_patch = 0))
  }
  uf <- UnionFind$new()
  for (pid in node_ids) {
    uf$find(pid)
  }
  if (!is.null(edge_specs) && length(edge_specs) > 0) {
    if (is.data.frame(edge_specs)) {
      for (i in seq_len(nrow(edge_specs))) {
        ids <- terralink_candidate_patch_ids(edge_specs[i, , drop = FALSE], valid_node_ids = node_ids)
        if (length(ids) < 2) next
        anchor <- ids[[1]]
        for (other in ids[-1]) {
          uf$union(anchor, other)
        }
      }
    } else if (is.matrix(edge_specs) && nrow(edge_specs) > 0) {
      for (i in seq_len(nrow(edge_specs))) {
        uf$union(as.integer(edge_specs[i, 1]), as.integer(edge_specs[i, 2]))
      }
    }
  }
  comp_area <- list()
  comp_count <- list()
  for (pid in node_ids) {
    root <- as.character(uf$find(pid))
    comp_area[[root]] <- as.numeric(comp_area[[root]] %||% 0) + as.numeric(patch_areas[[as.character(pid)]] %||% 0)
    comp_count[[root]] <- as.integer(comp_count[[root]] %||% 0L) + 1L
  }
  area_vals <- unlist(comp_area, use.names = TRUE)
  count_vals <- unlist(comp_count, use.names = TRUE)
  list(
    connected_area = sum(area_vals[count_vals >= 2], na.rm = TRUE),
    largest_network = max(area_vals, na.rm = TRUE),
    largest_patch = max(patch_areas, na.rm = TRUE)
  )
}

terralink_probability_of_connectivity <- function(pair_euclid, patch_areas, alpha, dist_matrix = NULL, cutoff = NULL) {
  n <- length(patch_areas)
  if (n <= 0) {
    return(0)
  }
  total_area <- sum(patch_areas)
  if (!is.finite(total_area) || total_area <= 0) {
    return(0)
  }
  self_term <- sum(patch_areas^2)
  pair_term <- 0
  if (n >= 2) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        d <- as.numeric(pair_euclid[i, j])
        if (!is.null(dist_matrix)) {
          g_d <- as.numeric(dist_matrix[i, j])
          if (is.finite(g_d) && g_d < d) {
            d <- g_d
          }
        }
        if (!is.finite(d)) next
        if (!is.null(cutoff) && is.finite(cutoff) && d > cutoff) next
        pij <- exp(-max(d, 0) / max(alpha, 1e-9))
        pair_term <- pair_term + (patch_areas[[i]] * patch_areas[[j]] * pij)
      }
    }
  }
  max(0, min(1, (self_term + (2 * pair_term)) / max(total_area^2, 1e-12)))
}

terralink_normalize_metric_weights <- function(metric_weights = NULL, weight_m = NULL, weight_lcc = NULL, weight_pc = NULL, weight_f = NULL) {
  if (is.null(metric_weights)) {
    metric_weights <- c(mesh = weight_m %||% 0.25, lcc = weight_lcc %||% 0.25, pc = weight_pc %||% 0.25, flow = weight_f %||% 0.25)
  }
  if (is.null(names(metric_weights)) || !all(c("mesh", "lcc", "pc", "flow") %in% names(metric_weights))) {
    vals <- as.numeric(metric_weights)
    vals <- vals[seq_len(min(4, length(vals)))]
    if (length(vals) < 4) {
      vals <- c(vals, rep(0.25, 4 - length(vals)))
    }
    names(vals) <- c("mesh", "lcc", "pc", "flow")
    metric_weights <- vals
  }
  metric_weights <- pmax(as.numeric(metric_weights[c("mesh", "lcc", "pc", "flow")]), 0)
  total <- sum(metric_weights)
  if (!is.finite(total) || total <= 0) {
    metric_weights <- c(0.25, 0.25, 0.25, 0.25)
  } else {
    metric_weights <- metric_weights / total
  }
  names(metric_weights) <- c("mesh", "lcc", "pc", "flow")
  metric_weights
}

terralink_metric_context <- function(
  patch_df,
  selected_corridors,
  label,
  area_unit,
  distance_unit,
  budget_used = 0,
  species_dispersal_distance = NULL,
  species_dispersal_kernel = HABITAT_AVAILABILITY_DEFAULT_KERNEL,
  min_patch_area_for_species = 0,
  patch_area_scaling = HABITAT_AVAILABILITY_DEFAULT_SCALING,
  max_search_distance = NULL,
  mobility_detour_cap = 8,
  redundancy_method = "ime",
  metric_weights = NULL,
  weight_m = NULL,
  weight_lcc = NULL,
  weight_pc = NULL,
  weight_f = NULL,
  metric_overrides = NULL
) {
  patch_df <- patch_df[order(as.integer(patch_df$patch_id)), , drop = FALSE]
  patch_areas <- as.numeric(patch_df$area)
  pair_euclid <- terralink_pair_euclid_matrix(data.frame(x = patch_df$x, y = patch_df$y))
  alpha <- as.numeric(species_dispersal_distance %||% max_search_distance %||% 0)
  if (!is.finite(alpha) || alpha <= 0) {
    finite_pair <- pair_euclid[upper.tri(pair_euclid)]
    finite_pair <- finite_pair[is.finite(finite_pair) & finite_pair > 0]
    alpha <- if (length(finite_pair) > 0) stats::median(finite_pair) else 1
  }
  cutoff <- max(alpha * 3, alpha)

  quality_weight <- patch_df$quality_weight %||% rep(1, nrow(patch_df))
  ha_nodes <- data.frame(
    node_id = patch_df$patch_id,
    raw_area = patch_df$area,
    effective_area = vapply(seq_len(nrow(patch_df)), function(i) {
      if (!is.finite(patch_df$area[[i]]) || patch_df$area[[i]] + 1e-12 < as.numeric(min_patch_area_for_species %||% 0)) {
        return(0)
      }
      terralink_scale_patch_area(patch_df$area[[i]], quality_weight = quality_weight[[i]], scaling = patch_area_scaling)
    }, numeric(1))
  )
  ha_state <- terralink_new_habitat_availability_state(
    ha_nodes,
    dispersal_distance = max(as.numeric(species_dispersal_distance %||% 0), alpha),
    kernel = species_dispersal_kernel
  )
  ha_before <- terralink_ha_metrics_from_state(ha_state)

  fluidity_state <- terralink_new_fluidity_state(
    data.frame(patch_id = patch_df$patch_id, x = patch_df$x, y = patch_df$y),
    detour_cap = mobility_detour_cap,
    shortcut_threshold = 3,
    redundancy_method = redundancy_method
  )
  fluidity_before <- fluidity_state$current_metrics

  edge_specs <- matrix(numeric(0), nrow = 0, ncol = 3)
  if (!is.null(selected_corridors) && nrow(selected_corridors) > 0) {
    for (i in seq_len(nrow(selected_corridors))) {
      ids <- terralink_candidate_patch_ids(selected_corridors[i, , drop = FALSE], valid_node_ids = patch_df$patch_id)
      if (length(ids) < 2) next
      len_val <- as.numeric(selected_corridors$length_metric[[i]] %||% selected_corridors$length[[i]] %||% selected_corridors$distance_m[[i]] %||% selected_corridors$distance[[i]] %||% 1)
      e <- terralink_ha_normalize_candidate_edges(ids, length = max(len_val, 1e-9), valid_node_ids = patch_df$patch_id)
      if (nrow(e) == 0) next
      edge_specs <- rbind(edge_specs, e)
      ha_eval <- terralink_ha_evaluate_candidate(ha_state, e, corridor_cost = as.numeric(selected_corridors$cost_metric[[i]] %||% selected_corridors$cost[[i]] %||% 1))
      if (!is.null(ha_eval)) {
        ha_state <- terralink_ha_apply_candidate(ha_state, ha_eval)
      }
      fl_eval <- terralink_fluidity_evaluate_candidate(fluidity_state, e, corridor_cost = as.numeric(selected_corridors$cost_metric[[i]] %||% selected_corridors$cost[[i]] %||% 1))
      if (!is.null(fl_eval)) {
        fluidity_state <- terralink_fluidity_apply_candidate(fluidity_state, fl_eval)
      }
    }
  }

  comp_stats <- terralink_component_stats(patch_df$patch_id, stats::setNames(patch_df$area, patch_df$patch_id), edge_specs)
  mesh_norm_pre <- if (sum(patch_areas) > 0) sum(patch_areas^2) / (sum(patch_areas)^2) else 0
  mesh_norm_post <- {
    comps <- terralink_component_stats(patch_df$patch_id, stats::setNames(patch_df$area, patch_df$patch_id), edge_specs)
    uf <- UnionFind$new()
    for (pid in patch_df$patch_id) uf$find(pid)
    if (nrow(edge_specs) > 0) {
      for (i in seq_len(nrow(edge_specs))) uf$union(as.integer(edge_specs[i, 1]), as.integer(edge_specs[i, 2]))
    }
    comp_area <- list()
    for (i in seq_len(nrow(patch_df))) {
      root <- as.character(uf$find(as.integer(patch_df$patch_id[[i]])))
      comp_area[[root]] <- as.numeric(comp_area[[root]] %||% 0) + as.numeric(patch_df$area[[i]])
    }
    areas <- unlist(comp_area, use.names = FALSE)
    if (sum(patch_areas) > 0) sum(areas^2) / (sum(patch_areas)^2) else 0
  }
  lcc_pre <- if (sum(patch_areas) > 0) max(patch_areas) / sum(patch_areas) else 0
  lcc_post <- if (sum(patch_areas) > 0) comp_stats$largest_network / sum(patch_areas) else 0
  pc_pre <- terralink_probability_of_connectivity(pair_euclid, patch_areas, alpha = alpha, dist_matrix = NULL, cutoff = cutoff)
  pc_post <- terralink_probability_of_connectivity(pair_euclid, patch_areas, alpha = alpha, dist_matrix = fluidity_state$dist, cutoff = cutoff)
  ha_after <- terralink_ha_metrics_from_state(ha_state)
  fluidity_metrics <- fluidity_state$current_metrics
  weights <- terralink_normalize_metric_weights(metric_weights, weight_m, weight_lcc, weight_pc, weight_f)
  composite_pre <- (weights[["mesh"]] * mesh_norm_pre) + (weights[["lcc"]] * lcc_pre) + (weights[["pc"]] * pc_pre) + (weights[["flow"]] * as.numeric(fluidity_before$flow_redundancy %||% 0))
  composite_post <- (weights[["mesh"]] * mesh_norm_post) + (weights[["lcc"]] * lcc_post) + (weights[["pc"]] * pc_post) + (weights[["flow"]] * as.numeric(fluidity_metrics$flow_redundancy %||% 0))

  metrics <- list(
    total_connected_habitat_area_pre = 0,
    total_connected_habitat_area_post = comp_stats$connected_area,
    largest_network_area_pre = comp_stats$largest_patch,
    largest_network_area_post = comp_stats$largest_network,
    habitat_availability_pre = ha_before$habitat_availability,
    habitat_availability_post = ha_after$habitat_availability,
    mean_effective_resistance_pre = as.numeric(fluidity_before$mean_effective_resistance %||% 0),
    mean_effective_resistance_post = as.numeric(fluidity_metrics$mean_effective_resistance %||% 0),
    mean_reachable_area_pre = ha_before$mean_reachable_area,
    mean_reachable_area_post = ha_after$mean_reachable_area,
    median_reachable_area_pre = ha_before$median_reachable_area,
    median_reachable_area_post = ha_after$median_reachable_area,
    largest_reachable_habitat_cluster_pre = ha_before$largest_reachable_habitat_cluster,
    largest_reachable_habitat_cluster_post = ha_after$largest_reachable_habitat_cluster,
    mesh_norm_pre = mesh_norm_pre,
    mesh_norm_post = mesh_norm_post,
    lcc_pre = lcc_pre,
    lcc_post = lcc_post,
    pc_pre = pc_pre,
    pc_post = pc_post,
    ime_pre = as.numeric(fluidity_before$ime %||% 0),
    ime_post = as.numeric(fluidity_metrics$ime %||% 0),
    flow_redundancy_pre = as.numeric(fluidity_before$flow_redundancy %||% 0),
    flow_redundancy_post = as.numeric(fluidity_metrics$flow_redundancy %||% 0),
    flow_method = as.character(fluidity_metrics$flow_method %||% "IME"),
    strategic_mobility_pre = as.numeric(fluidity_before$strategic_mobility %||% 0),
    strategic_mobility_post = as.numeric(fluidity_metrics$strategic_mobility %||% 0),
    landscape_fluidity_pre = as.numeric(fluidity_before$landscape_fluidity %||% 0),
    landscape_fluidity_post = as.numeric(fluidity_metrics$landscape_fluidity %||% 0),
    composite_connectivity_pre = composite_pre,
    composite_connectivity_post = composite_post,
    budget_used = as.numeric(budget_used %||% 0),
    area_unit = area_unit,
    distance_unit = distance_unit,
    label = label
  )
  if (is.list(metric_overrides) && length(metric_overrides) > 0) {
    metrics[names(metric_overrides)] <- metric_overrides
  }

  metric_rows <- list(
    list(name = "Total Connected Habitat Area", pre = metrics$total_connected_habitat_area_pre, post = metrics$total_connected_habitat_area_post, desc = paste("Habitat area connected into multi-patch networks (", area_unit, ")", sep = ""), higher = TRUE, digits = 4),
    list(name = "Largest Network Area", pre = metrics$largest_network_area_pre, post = metrics$largest_network_area_post, desc = paste("Habitat area in the largest connected network (", area_unit, ")", sep = ""), higher = TRUE, digits = 4),
    list(name = "Habitat Availability", pre = metrics$habitat_availability_pre, post = metrics$habitat_availability_post, desc = "Reachable habitat weighted by dispersal probability", higher = TRUE, digits = 6),
    list(name = "Mean Effective Resistance", pre = metrics$mean_effective_resistance_pre, post = metrics$mean_effective_resistance_post, desc = "Average movement difficulty across the patch graph", higher = FALSE, digits = 6),
    list(name = "Habitat-Normalized Mesh", pre = metrics$mesh_norm_pre, post = metrics$mesh_norm_post, desc = "Effective mesh size normalized by total habitat area", higher = TRUE, digits = 6),
    list(name = "Largest Connected Component Proportion", pre = metrics$lcc_pre, post = metrics$lcc_post, desc = "Share of habitat area captured by the largest network", higher = TRUE, digits = 6),
    list(name = "Probability of Connectivity", pre = metrics$pc_pre, post = metrics$pc_post, desc = "Area-weighted dispersal connectivity among patches", higher = TRUE, digits = 6),
    list(name = paste0("Flow Redundancy (", metrics$flow_method, ")"), pre = metrics$flow_redundancy_pre, post = metrics$flow_redundancy_post, desc = "Redundancy / alternative-route efficiency", higher = TRUE, digits = 6),
    list(name = "Strategic Mobility", pre = metrics$strategic_mobility_pre, post = metrics$strategic_mobility_post, desc = "Inverse detour ratio across the patch graph", higher = TRUE, digits = 6),
    list(name = "Landscape Fluidity", pre = metrics$landscape_fluidity_pre, post = metrics$landscape_fluidity_post, desc = "Inverse mean path resistance on the patch graph", higher = TRUE, digits = 6),
    list(name = "Composite Connectivity Score", pre = metrics$composite_connectivity_pre, post = metrics$composite_connectivity_post, desc = "Weighted blend of mesh, LCC, PC, and flow redundancy", higher = TRUE, digits = 6)
  )

  fmt_num <- function(x, digits) {
    if (!is.finite(as.numeric(x))) return("")
    formatC(as.numeric(x), format = "f", digits = digits)
  }
  report <- c(
    paste(rep("=", 120), collapse = ""),
    paste("LANDSCAPE ANALYSIS:", label),
    paste(rep("=", 120), collapse = "")
  )
  col_header <- sprintf("%-40s | %-14s | %-14s | %-14s | %-14s | %s", "METRIC NAME", "PRE", "POST", "DELTA", "GAIN/BUDGET", "DESCRIPTION")
  report <- c(report, col_header, paste(rep("-", nchar(col_header)), collapse = ""))
  for (row in metric_rows) {
    delta <- as.numeric(row$post) - as.numeric(row$pre)
    beneficial_gain <- if (isTRUE(row$higher)) delta else -delta
    gain_per_budget <- if (is.finite(metrics$budget_used) && metrics$budget_used > 1e-12) beneficial_gain / metrics$budget_used else 0
    report <- c(
      report,
      sprintf(
        "%-40s | %-14s | %-14s | %-14s | %-14s | %s",
        row$name,
        fmt_num(row$pre, row$digits),
        fmt_num(row$post, row$digits),
        sprintf("%+.*f", row$digits, delta),
        sprintf("%+0.6f", gain_per_budget),
        row$desc
      )
    )
  }
  report <- c(report, paste(rep("=", 120), collapse = ""))
  list(report = report, metrics = metrics)
}

terralink_optimize_connected_area_df <- function(
  candidates,
  patch_df,
  budget,
  max_redundant_links = 1L,
  cost_col = "cost_metric",
  length_col = "length_metric"
) {
  budget_limit <- max(as.numeric(budget %||% 0), 0)
  if (budget_limit <= 0 || !is.data.frame(candidates) || nrow(candidates) == 0 || !is.data.frame(patch_df) || nrow(patch_df) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "most_connected_habitat", corridors_used = 0L, budget_used = 0)))
  }

  node_ids <- unique(as.integer(patch_df$patch_id[is.finite(patch_df$patch_id)]))
  if (length(node_ids) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "most_connected_habitat", corridors_used = 0L, budget_used = 0)))
  }

  area_lookup <- stats::setNames(as.numeric(patch_df$area[match(node_ids, patch_df$patch_id)]), as.character(node_ids))
  uf <- UnionFind$new()
  comp_area <- list()
  comp_count <- list()
  for (pid in node_ids) {
    root <- as.character(uf$find(pid))
    comp_area[[root]] <- as.numeric(area_lookup[[as.character(pid)]] %||% 0)
    comp_count[[root]] <- 1L
  }

  current_connected_area <- 0
  current_largest_network <- max(unlist(comp_area, use.names = FALSE), 0, na.rm = TRUE)
  remaining <- budget_limit
  budget_used <- 0
  selected <- integer(0)
  selected_types <- character(0)
  redundant_links_used <- 0L
  available <- seq_len(nrow(candidates))

  while (remaining > 1e-12 && length(available) > 0) {
    best_primary_idx <- NA_integer_
    best_primary_eval <- NULL
    best_primary_score <- -Inf
    best_primary_gain <- -Inf
    best_primary_cost <- Inf
    best_primary_length <- Inf
    best_redundant_idx <- NA_integer_
    best_redundant_eval <- NULL
    best_redundant_score <- -Inf
    best_redundant_gain <- -Inf
    best_redundant_cost <- Inf
    best_redundant_length <- Inf

    for (idx in available) {
      cand <- candidates[idx, , drop = FALSE]
      cost <- as.numeric(cand[[cost_col]][[1]] %||% cand$cost[[1]] %||% cand$area[[1]] %||% 0)
      if (!is.finite(cost) || cost <= 0 || cost > remaining + 1e-12) next
      ids <- terralink_candidate_patch_ids(cand, valid_node_ids = node_ids)
      if (length(ids) < 2) next
      roots <- unique(vapply(ids, function(pid) as.character(uf$find(pid)), character(1)))
      root_areas <- vapply(roots, function(root) as.numeric(comp_area[[root]] %||% 0), numeric(1))
      root_counts <- vapply(roots, function(root) as.integer(comp_count[[root]] %||% 0L), integer(1))
      length_val <- as.numeric(cand[[length_col]][[1]] %||% cand$length[[1]] %||% cand$distance[[1]] %||% 1)
      if (!is.finite(length_val) || length_val <= 0) {
        length_val <- cost
      }

      if (length(roots) > 1) {
        merged_area <- sum(root_areas)
        merged_count <- sum(root_counts)
        replaced_connected <- sum(root_areas[root_counts >= 2], na.rm = TRUE)
        connected_after <- current_connected_area - replaced_connected + if (merged_count >= 2) merged_area else 0
        largest_after <- max(current_largest_network, merged_area, na.rm = TRUE)
        gain_connected <- connected_after - current_connected_area
        gain_largest <- largest_after - current_largest_network
        gain <- max(gain_connected, 0) + max(gain_largest, 0)
        score <- gain / max(cost, 1e-12)
        better <-
          score > best_primary_score + 1e-12 ||
          (abs(score - best_primary_score) <= 1e-12 && gain > best_primary_gain + 1e-12) ||
          (abs(score - best_primary_score) <= 1e-12 && abs(gain - best_primary_gain) <= 1e-12 && cost < best_primary_cost - 1e-12) ||
          (abs(score - best_primary_score) <= 1e-12 && abs(gain - best_primary_gain) <= 1e-12 && abs(cost - best_primary_cost) <= 1e-12 && length_val < best_primary_length - 1e-12)
        if (isTRUE(better)) {
          best_primary_idx <- idx
          best_primary_eval <- list(
            ids = ids,
            roots = roots,
            merged_area = merged_area,
            merged_count = merged_count,
            connected_after = connected_after,
            largest_after = largest_after
          )
          best_primary_score <- score
          best_primary_gain <- gain
          best_primary_cost <- cost
          best_primary_length <- length_val
        }
      } else {
        if (redundant_links_used >= as.integer(max_redundant_links %||% 0L)) next
        comp_area_val <- root_areas[[1]]
        gain <- max(as.numeric(comp_area_val %||% 0), 0) / max(length(ids), 2)
        score <- (gain / max(cost, 1e-12)) * 0.05
        better <-
          score > best_redundant_score + 1e-12 ||
          (abs(score - best_redundant_score) <= 1e-12 && gain > best_redundant_gain + 1e-12) ||
          (abs(score - best_redundant_score) <= 1e-12 && abs(gain - best_redundant_gain) <= 1e-12 && cost < best_redundant_cost - 1e-12) ||
          (abs(score - best_redundant_score) <= 1e-12 && abs(gain - best_redundant_gain) <= 1e-12 && abs(cost - best_redundant_cost) <= 1e-12 && length_val < best_redundant_length - 1e-12)
        if (isTRUE(better)) {
          best_redundant_idx <- idx
          best_redundant_eval <- list(ids = ids)
          best_redundant_score <- score
          best_redundant_gain <- gain
          best_redundant_cost <- cost
          best_redundant_length <- length_val
        }
      }
    }

    if (is.finite(best_primary_idx) && !is.null(best_primary_eval) && best_primary_score > 1e-12) {
      chosen_idx <- best_primary_idx
      chosen_eval <- best_primary_eval
      chosen_cost <- best_primary_cost
      chosen_type <- "primary"
    } else if (is.finite(best_redundant_idx) && !is.null(best_redundant_eval) && best_redundant_score > 1e-12) {
      chosen_idx <- best_redundant_idx
      chosen_eval <- best_redundant_eval
      chosen_cost <- best_redundant_cost
      chosen_type <- "redundant"
    } else {
      break
    }

    if (identical(chosen_type, "primary")) {
      anchor <- chosen_eval$ids[[1]]
      for (other in chosen_eval$ids[-1]) {
        uf$union(anchor, other)
      }
      new_root <- as.character(uf$find(anchor))
      comp_area[[new_root]] <- as.numeric(chosen_eval$merged_area %||% 0)
      comp_count[[new_root]] <- as.integer(chosen_eval$merged_count %||% 0L)
      old_roots <- setdiff(unique(as.character(chosen_eval$roots)), new_root)
      for (old_root in old_roots) {
        comp_area[[old_root]] <- 0
        comp_count[[old_root]] <- 0L
      }
      current_connected_area <- as.numeric(chosen_eval$connected_after %||% current_connected_area)
      current_largest_network <- as.numeric(chosen_eval$largest_after %||% current_largest_network)
    } else {
      redundant_links_used <- redundant_links_used + 1L
    }

    selected <- c(selected, as.integer(candidates$id[[chosen_idx]] %||% chosen_idx))
    selected_types <- c(selected_types, chosen_type)
    budget_used <- budget_used + chosen_cost
    remaining <- remaining - chosen_cost
    available <- setdiff(available, chosen_idx)
  }

  list(
    selected = as.integer(selected),
    stats = list(
      strategy = "most_connected_habitat",
      corridors_used = length(selected),
      budget_used = budget_used,
      primary_links = sum(selected_types == "primary"),
      redundant_links = sum(selected_types == "redundant"),
      connected_area_after = current_connected_area,
      largest_network_after = current_largest_network
    )
  )
}

terralink_optimize_habitat_availability_df <- function(
  candidates,
  patch_df,
  budget,
  min_patch_area_for_species = 0,
  patch_area_scaling = HABITAT_AVAILABILITY_DEFAULT_SCALING,
  species_dispersal_distance,
  species_dispersal_kernel = HABITAT_AVAILABILITY_DEFAULT_KERNEL,
  max_redundant_links = 1L,
  cost_col = "cost_metric",
  length_col = "length_metric"
) {
  budget_limit <- max(as.numeric(budget %||% 0), 0)
  if (budget_limit <= 0 || !is.data.frame(candidates) || nrow(candidates) == 0 || !is.data.frame(patch_df) || nrow(patch_df) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "reachable_habitat_advanced", corridors_used = 0L, budget_used = 0)))
  }
  nodes <- data.frame(
    node_id = patch_df$patch_id,
    raw_area = patch_df$area,
    effective_area = vapply(seq_len(nrow(patch_df)), function(i) {
      if (patch_df$area[[i]] + 1e-12 < max(as.numeric(min_patch_area_for_species %||% 0), 0)) {
        return(0)
      }
      terralink_scale_patch_area(
        area = patch_df$area[[i]],
        quality_weight = (patch_df$quality_weight %||% rep(1, nrow(patch_df)))[[i]],
        scaling = patch_area_scaling
      )
    }, numeric(1))
  )
  state <- terralink_new_habitat_availability_state(
    nodes = nodes,
    dispersal_distance = species_dispersal_distance,
    kernel = species_dispersal_kernel
  )
  if (length(state$node_ids) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "reachable_habitat_advanced", corridors_used = 0L, budget_used = 0)))
  }
  before_metrics <- terralink_ha_metrics_from_state(state)
  remaining <- budget_limit
  budget_used <- 0
  selected <- integer(0)
  selected_types <- character(0)
  redundant_links_used <- 0L
  available <- seq_len(nrow(candidates))

  while (remaining > 1e-12 && length(available) > 0) {
    best_primary_idx <- NA_integer_
    best_primary_eval <- NULL
    best_primary_score <- -Inf
    best_primary_gain <- -Inf
    best_primary_cost <- Inf
    best_primary_length <- Inf
    best_redundant_idx <- NA_integer_
    best_redundant_eval <- NULL
    best_redundant_score <- -Inf
    best_redundant_gain <- -Inf
    best_redundant_cost <- Inf
    best_redundant_length <- Inf
    for (idx in available) {
      cand <- candidates[idx, , drop = FALSE]
      cost <- as.numeric(cand[[cost_col]][[1]] %||% cand$cost[[1]] %||% cand$area[[1]] %||% 0)
      if (!is.finite(cost) || cost <= 0 || cost > remaining + 1e-12) next
      ids <- terralink_candidate_patch_ids(cand, valid_node_ids = state$node_ids)
      if (length(ids) < 2) next
      length_val <- as.numeric(cand[[length_col]][[1]] %||% cand$length[[1]] %||% cand$distance[[1]] %||% 1)
      edge_specs <- terralink_ha_normalize_candidate_edges(ids, length = max(length_val, 1e-9), valid_node_ids = state$node_ids)
      eval_out <- terralink_ha_evaluate_candidate(state, edge_specs, corridor_cost = cost)
      if (is.null(eval_out)) next
      score <- as.numeric(eval_out$score %||% 0)
      gain <- as.numeric(eval_out$gain %||% 0)
      candidate_type <- as.character(eval_out$type %||% "primary")
      if (!identical(candidate_type, "primary") && redundant_links_used >= as.integer(max_redundant_links %||% 0L)) {
        next
      }
      if (identical(candidate_type, "primary")) {
        better <-
          score > best_primary_score + 1e-12 ||
          (abs(score - best_primary_score) <= 1e-12 && gain > best_primary_gain + 1e-12) ||
          (abs(score - best_primary_score) <= 1e-12 && abs(gain - best_primary_gain) <= 1e-12 && cost < best_primary_cost - 1e-12) ||
          (abs(score - best_primary_score) <= 1e-12 && abs(gain - best_primary_gain) <= 1e-12 && abs(cost - best_primary_cost) <= 1e-12 && length_val < best_primary_length - 1e-12)
        if (isTRUE(better)) {
          best_primary_idx <- idx
          best_primary_eval <- eval_out
          best_primary_score <- score
          best_primary_gain <- gain
          best_primary_cost <- cost
          best_primary_length <- length_val
        }
      } else {
        better <-
          score > best_redundant_score + 1e-12 ||
          (abs(score - best_redundant_score) <= 1e-12 && gain > best_redundant_gain + 1e-12) ||
          (abs(score - best_redundant_score) <= 1e-12 && abs(gain - best_redundant_gain) <= 1e-12 && cost < best_redundant_cost - 1e-12) ||
          (abs(score - best_redundant_score) <= 1e-12 && abs(gain - best_redundant_gain) <= 1e-12 && abs(cost - best_redundant_cost) <= 1e-12 && length_val < best_redundant_length - 1e-12)
        if (isTRUE(better)) {
          best_redundant_idx <- idx
          best_redundant_eval <- eval_out
          best_redundant_score <- score
          best_redundant_gain <- gain
          best_redundant_cost <- cost
          best_redundant_length <- length_val
        }
      }
    }
    if (is.finite(best_primary_idx) && !is.null(best_primary_eval) && best_primary_score > 1e-12) {
      chosen_idx <- best_primary_idx
      chosen_eval <- best_primary_eval
      chosen_cost <- best_primary_cost
      chosen_type <- "primary"
    } else if (is.finite(best_redundant_idx) && !is.null(best_redundant_eval) && best_redundant_score > 1e-12) {
      chosen_idx <- best_redundant_idx
      chosen_eval <- best_redundant_eval
      chosen_cost <- best_redundant_cost
      chosen_type <- "redundant"
    } else {
      break
    }
    state <- terralink_ha_apply_candidate(state, chosen_eval)
    selected <- c(selected, as.integer(candidates$id[[chosen_idx]] %||% chosen_idx))
    selected_types <- c(selected_types, chosen_type)
    if (identical(chosen_type, "redundant")) {
      redundant_links_used <- redundant_links_used + 1L
    }
    budget_used <- budget_used + chosen_cost
    remaining <- remaining - chosen_cost
    available <- setdiff(available, chosen_idx)
  }
  after_metrics <- terralink_ha_metrics_from_state(state)
  list(
    selected = as.integer(selected),
    stats = list(
      strategy = "reachable_habitat_advanced",
      corridors_used = length(selected),
      budget_used = budget_used,
      primary_links = sum(selected_types == "primary"),
      redundant_links = sum(selected_types == "redundant"),
      habitat_availability_before = before_metrics$habitat_availability,
      habitat_availability_after = after_metrics$habitat_availability,
      habitat_availability_gain = after_metrics$habitat_availability - before_metrics$habitat_availability,
      mean_reachable_area_before = before_metrics$mean_reachable_area,
      mean_reachable_area = after_metrics$mean_reachable_area,
      median_reachable_area_before = before_metrics$median_reachable_area,
      median_reachable_area = after_metrics$median_reachable_area,
      largest_reachable_habitat_cluster_before = before_metrics$largest_reachable_habitat_cluster,
      largest_reachable_habitat_cluster = after_metrics$largest_reachable_habitat_cluster,
      species_dispersal_distance = species_dispersal_distance,
      species_dispersal_kernel = terralink_normalize_habitat_availability_kernel(species_dispersal_kernel),
      patch_area_scaling = terralink_normalize_patch_area_scaling(patch_area_scaling)
    )
  )
}

terralink_optimize_landscape_fluidity_df <- function(
  candidates,
  patch_df,
  budget,
  shortcut_threshold = 3,
  detour_cap = 8,
  redundancy_method = "ime",
  max_redundant_links = 1L,
  cost_col = "cost_metric",
  length_col = "length_metric"
) {
  budget_limit <- max(as.numeric(budget %||% 0), 0)
  if (budget_limit <= 0 || !is.data.frame(candidates) || nrow(candidates) == 0 || !is.data.frame(patch_df) || nrow(patch_df) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "landscape_fluidity", corridors_used = 0L, budget_used = 0)))
  }
  state <- terralink_new_fluidity_state(
    data.frame(patch_id = patch_df$patch_id, x = patch_df$x, y = patch_df$y),
    detour_cap = detour_cap,
    shortcut_threshold = shortcut_threshold,
    redundancy_method = redundancy_method
  )
  remaining <- budget_limit
  budget_used <- 0
  selected <- integer(0)
  selected_types <- character(0)
  redundant_links_used <- 0L
  available <- seq_len(nrow(candidates))
  before_metrics <- state$current_metrics

  while (remaining > 1e-12 && length(available) > 0) {
    best_primary_idx <- NA_integer_
    best_primary_eval <- NULL
    best_primary_score <- -Inf
    best_primary_gain <- -Inf
    best_primary_cost <- Inf
    best_primary_length <- Inf
    best_redundant_idx <- NA_integer_
    best_redundant_eval <- NULL
    best_redundant_score <- -Inf
    best_redundant_gain <- -Inf
    best_redundant_cost <- Inf
    best_redundant_length <- Inf
    for (idx in available) {
      cand <- candidates[idx, , drop = FALSE]
      cost <- as.numeric(cand[[cost_col]][[1]] %||% cand$cost[[1]] %||% cand$area[[1]] %||% 0)
      if (!is.finite(cost) || cost <= 0 || cost > remaining + 1e-12) next
      ids <- terralink_candidate_patch_ids(cand, valid_node_ids = patch_df$patch_id)
      if (length(ids) < 2) next
      length_val <- as.numeric(cand[[length_col]][[1]] %||% cand$length[[1]] %||% cand$distance[[1]] %||% 1)
      edge_specs <- terralink_ha_normalize_candidate_edges(ids, length = max(length_val, 1e-9), valid_node_ids = patch_df$patch_id)
      eval_out <- terralink_fluidity_evaluate_candidate(state, edge_specs, corridor_cost = cost)
      if (is.null(eval_out)) next
      score <- as.numeric(eval_out$score %||% 0)
      gain <- as.numeric(eval_out$gain %||% 0)
      candidate_type <- as.character(eval_out$type %||% "primary")
      if (!identical(candidate_type, "primary") && redundant_links_used >= as.integer(max_redundant_links %||% 0L)) {
        next
      }
      if (identical(candidate_type, "primary")) {
        better <-
          score > best_primary_score + 1e-12 ||
          (abs(score - best_primary_score) <= 1e-12 && gain > best_primary_gain + 1e-12) ||
          (abs(score - best_primary_score) <= 1e-12 && abs(gain - best_primary_gain) <= 1e-12 && cost < best_primary_cost - 1e-12) ||
          (abs(score - best_primary_score) <= 1e-12 && abs(gain - best_primary_gain) <= 1e-12 && abs(cost - best_primary_cost) <= 1e-12 && length_val < best_primary_length - 1e-12)
        if (isTRUE(better)) {
          best_primary_idx <- idx
          best_primary_eval <- eval_out
          best_primary_score <- score
          best_primary_gain <- gain
          best_primary_cost <- cost
          best_primary_length <- length_val
        }
      } else {
        better <-
          score > best_redundant_score + 1e-12 ||
          (abs(score - best_redundant_score) <= 1e-12 && gain > best_redundant_gain + 1e-12) ||
          (abs(score - best_redundant_score) <= 1e-12 && abs(gain - best_redundant_gain) <= 1e-12 && cost < best_redundant_cost - 1e-12) ||
          (abs(score - best_redundant_score) <= 1e-12 && abs(gain - best_redundant_gain) <= 1e-12 && abs(cost - best_redundant_cost) <= 1e-12 && length_val < best_redundant_length - 1e-12)
        if (isTRUE(better)) {
          best_redundant_idx <- idx
          best_redundant_eval <- eval_out
          best_redundant_score <- score
          best_redundant_gain <- gain
          best_redundant_cost <- cost
          best_redundant_length <- length_val
        }
      }
    }
    if (is.finite(best_primary_idx) && !is.null(best_primary_eval) && best_primary_score > 1e-12) {
      chosen_idx <- best_primary_idx
      chosen_eval <- best_primary_eval
      chosen_cost <- best_primary_cost
      chosen_type <- "primary"
    } else if (is.finite(best_redundant_idx) && !is.null(best_redundant_eval) && best_redundant_score > 1e-12) {
      chosen_idx <- best_redundant_idx
      chosen_eval <- best_redundant_eval
      chosen_cost <- best_redundant_cost
      chosen_type <- "redundant"
    } else {
      break
    }
    state <- terralink_fluidity_apply_candidate(state, chosen_eval)
    selected <- c(selected, as.integer(candidates$id[[chosen_idx]] %||% chosen_idx))
    selected_types <- c(selected_types, chosen_type)
    if (identical(chosen_type, "redundant")) {
      redundant_links_used <- redundant_links_used + 1L
    }
    budget_used <- budget_used + chosen_cost
    remaining <- remaining - chosen_cost
    available <- setdiff(available, chosen_idx)
  }
  after_metrics <- state$current_metrics
  list(
    selected = as.integer(selected),
    stats = list(
      strategy = "landscape_fluidity",
      corridors_used = length(selected),
      budget_used = budget_used,
      primary_links = sum(selected_types %in% c("primary", "backbone")),
      redundant_links = sum(selected_types == "redundant"),
      mean_effective_resistance_before = before_metrics$mean_effective_resistance,
      mean_effective_resistance_after = after_metrics$mean_effective_resistance,
      landscape_fluidity_pre = before_metrics$landscape_fluidity,
      landscape_fluidity_post = after_metrics$landscape_fluidity,
      landscape_fluidity_gain = after_metrics$landscape_fluidity - before_metrics$landscape_fluidity,
      flow_redundancy_pre = before_metrics$flow_redundancy,
      flow_redundancy_post = after_metrics$flow_redundancy,
      strategic_mobility_pre = before_metrics$strategic_mobility,
      strategic_mobility_post = after_metrics$strategic_mobility,
      redundancy_method = after_metrics$flow_method
    )
  )
}

terralink_bigconnect_key <- function(value, scale = 1) {
  val <- suppressWarnings(as.numeric(value))
  if (!is.finite(val) || val <= 0) return(0L)
  as.integer(round(val * scale))
}

terralink_bigconnect_score_is_better <- function(candidate_score, incumbent_score) {
  if (is.null(incumbent_score)) return(TRUE)
  for (i in seq_along(candidate_score)) {
    cand <- as.numeric(candidate_score[[i]])
    inc <- as.numeric(incumbent_score[[i]])
    if (cand > inc + 1e-12) return(TRUE)
    if (cand < inc - 1e-12) return(FALSE)
  }
  FALSE
}

terralink_bigconnect_candidate_value <- function(cand, fields, default = NA_real_) {
  for (nm in fields) {
    if (!nm %in% names(cand)) next
    raw <- cand[[nm]]
    if (is.list(raw) && length(raw) == 1L) raw <- raw[[1]]
    if (length(raw) == 0) next
    val <- suppressWarnings(as.numeric(raw[[1]] %||% raw))
    if (is.finite(val)) return(val)
  }
  suppressWarnings(as.numeric(default))
}

terralink_bigconnect_candidate_id_local <- function(cand, fallback_id) {
  if (!"id" %in% names(cand)) return(as.integer(fallback_id))
  raw <- cand$id
  if (is.list(raw) && length(raw) == 1L) raw <- raw[[1]]
  out <- suppressWarnings(as.integer(raw[[1]] %||% raw))
  if (is.finite(out)) out else as.integer(fallback_id)
}

terralink_bigconnect_score_tuple <- function(
  connected_area_key,
  cohesion_key,
  budget_used_key,
  corridor_count,
  total_length
) {
  c(
    as.numeric(connected_area_key),
    -as.numeric(budget_used_key),
    -as.numeric(total_length),
    as.numeric(cohesion_key),
    -as.numeric(corridor_count)
  )
}

terralink_bigconnect_canonical_signature <- function(sig) {
  mapping <- list()
  out <- integer(length(sig))
  next_id <- 0L
  for (i in seq_along(sig)) {
    key <- as.character(sig[[i]])
    if (is.null(mapping[[key]])) {
      mapping[[key]] <- next_id
      next_id <- next_id + 1L
    }
    out[[i]] <- as.integer(mapping[[key]])
  }
  out
}

terralink_bigconnect_union_signature <- function(sig, members) {
  idx <- sort(unique(as.integer(members)))
  idx <- idx[is.finite(idx)]
  if (length(idx) < 2L) return(as.integer(sig))
  root_val <- sig[[idx[[1]] + 1L]]
  targets <- unique(sig[idx + 1L])
  updated <- as.integer(sig)
  for (i in seq_along(updated)) {
    if (updated[[i]] %in% targets) updated[[i]] <- as.integer(root_val)
  }
  terralink_bigconnect_canonical_signature(updated)
}

terralink_bigconnect_signature_groups <- function(sig) {
  split(seq_along(sig) - 1L, sig)
}

terralink_bigconnect_connected_area <- function(sig, patch_area_keys) {
  groups <- terralink_bigconnect_signature_groups(sig)
  total <- 0L
  for (members in groups) {
    if (length(members) < 2L) next
    total <- total + as.integer(sum(as.integer(patch_area_keys[members + 1L])))
  }
  total
}

terralink_bigconnect_cohesion_key <- function(sig, patch_area_keys) {
  groups <- terralink_bigconnect_signature_groups(sig)
  total <- 0
  for (members in groups) {
    if (length(members) < 2L) next
    comp_area <- sum(as.numeric(patch_area_keys[members + 1L]))
    total <- total + (comp_area * comp_area)
  }
  total
}

terralink_bigconnect_remaining_upper_bound <- function(sig, patch_area_keys, remaining_patch_sets) {
  current <- terralink_bigconnect_connected_area(sig, patch_area_keys)
  groups <- terralink_bigconnect_signature_groups(sig)
  counted <- integer(0)
  for (members in groups) {
    if (length(members) >= 2L) counted <- unique(c(counted, as.integer(members)))
  }
  incident_remaining <- integer(0)
  for (pset in remaining_patch_sets) {
    incident_remaining <- unique(c(incident_remaining, as.integer(pset)))
  }
  extra <- 0
  for (idx in incident_remaining) {
    if (!(idx %in% counted)) {
      extra <- extra + as.numeric(patch_area_keys[[idx + 1L]])
    }
  }
  current + extra
}

terralink_bigconnect_state_key <- function(sig, spend_key) {
  paste0(paste(sig, collapse = ","), "|", as.integer(spend_key))
}

terralink_bigconnect_keep_best_frontier_rows <- function(rows) {
  best_by_spend <- list()
  for (row in rows) {
    spend_key <- as.integer(row$budget_used_key %||% 0L)
    score <- terralink_bigconnect_score_tuple(
      connected_area_key = row$connected_area_key %||% 0,
      cohesion_key = row$cohesion_key %||% 0,
      budget_used_key = spend_key,
      corridor_count = row$corridor_count %||% 0L,
      total_length = row$total_length %||% 0
    )
    prev <- best_by_spend[[as.character(spend_key)]]
    prev_score <- NULL
    if (!is.null(prev)) {
      prev_score <- terralink_bigconnect_score_tuple(
        connected_area_key = prev$connected_area_key %||% 0,
        cohesion_key = prev$cohesion_key %||% 0,
        budget_used_key = prev$budget_used_key %||% 0L,
        corridor_count = prev$corridor_count %||% 0L,
        total_length = prev$total_length %||% 0
      )
    }
    if (terralink_bigconnect_score_is_better(score, prev_score)) {
      best_by_spend[[as.character(spend_key)]] <- row
    }
  }
  spend_vals <- sort(as.integer(names(best_by_spend)))
  unname(lapply(spend_vals, function(spend_key) best_by_spend[[as.character(spend_key)]]))
}

terralink_bigconnect_build_canonical_candidates_df <- function(
  candidates,
  patch_df,
  budget,
  cost_col = "cost_metric",
  length_col = "length_metric",
  scale = 1,
  merge_equiv_area = 0,
  merge_equiv_ratio = 0.01
) {
  if (!is.data.frame(candidates) || nrow(candidates) == 0 || !is.data.frame(patch_df) || nrow(patch_df) == 0) {
    return(list(canonical = list(), stats = list(raw_count = 0L, nondominated_count = 0L)))
  }
  valid_node_ids <- suppressWarnings(as.integer(patch_df$patch_id))
  valid_node_ids <- valid_node_ids[is.finite(valid_node_ids)]
  patch_area_by_id <- stats::setNames(as.numeric(patch_df$area), as.character(patch_df$patch_id))
  patch_area_keys <- lapply(patch_area_by_id, function(value) round(as.numeric(value) * scale))
  budget_key_limit <- terralink_bigconnect_key(budget, scale = scale)
  grouped <- list()
  raw_count <- 0L

  endpoint_pair_key <- function(cand, patch_ids) {
    p1 <- terralink_bigconnect_candidate_value(cand, c("patch1", "p1"), default = patch_ids[[1]])
    p2 <- terralink_bigconnect_candidate_value(cand, c("patch2", "p2"), default = patch_ids[[length(patch_ids)]])
    pair <- sort(c(as.integer(p1), as.integer(p2)))
    paste(pair, collapse = "_")
  }

  for (idx in seq_len(nrow(candidates))) {
    cand <- candidates[idx, , drop = FALSE]
    cost <- terralink_bigconnect_candidate_value(cand, c(cost_col, "cost", "area"), default = 0)
    cost_key <- terralink_bigconnect_key(cost, scale = scale)
    if (cost_key <= 0L || cost_key > budget_key_limit) next
    patch_ids <- terralink_candidate_patch_ids(cand, valid_node_ids = valid_node_ids)
    raw_count <- raw_count + 1L
    if (length(patch_ids) < 2L) next
    patch_ids <- sort(unique(as.integer(patch_ids)))
    length_val <- terralink_bigconnect_candidate_value(cand, c(length_col, "length", "distance_m", "distance"), default = cost)
    if (!is.finite(length_val) || length_val <= 0) length_val <- cost
    row <- list(
      candidate_index = idx,
      candidate_id = terralink_bigconnect_candidate_id_local(cand, idx),
      patch_ids = patch_ids,
      endpoint_pair = endpoint_pair_key(cand, patch_ids),
      cost = as.numeric(cost),
      cost_key = as.integer(cost_key),
      length = as.numeric(length_val),
      potential_area_key = sum(as.numeric(unlist(patch_area_keys[as.character(patch_ids)])))
    )
    key <- paste(patch_ids, collapse = "_")
    grouped[[key]] <- c(grouped[[key]], list(row))
  }

  canonical_raw <- list()
  canon_id <- 1L
  for (key in sort(names(grouped))) {
    rows <- grouped[[key]]
    order_idx <- order(
      vapply(rows, function(rec) rec$cost_key, integer(1)),
      vapply(rows, function(rec) rec$length, numeric(1)),
      vapply(rows, function(rec) rec$candidate_id, integer(1))
    )
    rows <- rows[order_idx]
    kept <- list()
    for (row in rows) {
      dominated <- FALSE
      for (prior in kept) {
        if (
          as.integer(prior$cost_key) <= as.integer(row$cost_key) &&
          as.numeric(prior$length) <= as.numeric(row$length) + 1e-12
        ) {
          dominated <- TRUE
          break
        }
      }
      if (isTRUE(dominated)) next
      row$canon_id <- canon_id
      kept[[length(kept) + 1L]] <- row
      canonical_raw[[length(canonical_raw) + 1L]] <- row
      canon_id <- canon_id + 1L
    }
  }

  collapsed <- list()
  endpoint_groups <- split(canonical_raw, vapply(canonical_raw, function(rec) rec$endpoint_pair, character(1)))
  merge_equiv_key <- terralink_bigconnect_key(merge_equiv_area, scale = scale)
  for (endpoint_key in sort(names(endpoint_groups))) {
    rows <- endpoint_groups[[endpoint_key]]
    rows <- rows[order(
      -vapply(rows, function(rec) rec$potential_area_key, numeric(1)),
      vapply(rows, function(rec) rec$cost_key, integer(1)),
      vapply(rows, function(rec) rec$length, numeric(1)),
      vapply(rows, function(rec) rec$candidate_id, integer(1))
    )]
    while (length(rows) > 0L) {
      ref_area_key <- as.numeric(rows[[1]]$potential_area_key)
      area_eps_key <- max(
        as.numeric(merge_equiv_key),
        round(as.numeric(ref_area_key) * as.numeric(merge_equiv_ratio))
      )
      tier_idx <- which(vapply(rows, function(rec) ref_area_key - as.numeric(rec$potential_area_key) <= area_eps_key, logical(1)))
      tier <- rows[tier_idx]
      best_idx <- order(
        vapply(tier, function(rec) rec$cost_key, integer(1)),
        vapply(tier, function(rec) rec$length, numeric(1)),
        vapply(tier, function(rec) rec$candidate_id, integer(1))
      )[[1]]
      collapsed[[length(collapsed) + 1L]] <- tier[[best_idx]]
      keep_ids <- vapply(tier, function(rec) rec$candidate_id, integer(1))
      rows <- rows[!(vapply(rows, function(rec) rec$candidate_id, integer(1)) %in% keep_ids)]
    }
  }

  canonical <- collapsed[order(
    vapply(collapsed, function(rec) paste(rec$patch_ids, collapse = "_"), character(1)),
    vapply(collapsed, function(rec) rec$cost_key, integer(1)),
    vapply(collapsed, function(rec) rec$length, numeric(1)),
    vapply(collapsed, function(rec) rec$candidate_id, integer(1))
  )]
  list(
    canonical = canonical,
    stats = list(raw_count = as.integer(raw_count), nondominated_count = as.integer(length(canonical)))
  )
}

terralink_bigconnect_build_clusters_df <- function(canonical, patch_df, scale = 1) {
  if (length(canonical) == 0L) return(list())
  area_lookup <- stats::setNames(as.numeric(patch_df$area), as.character(patch_df$patch_id))
  valid_ids <- suppressWarnings(as.integer(names(area_lookup)))
  valid_ids <- valid_ids[is.finite(valid_ids)]

  uf <- UnionFind$new()
  for (pid in valid_ids) uf$find(pid)
  for (row in canonical) {
    pids <- as.integer(row$patch_ids)
    if (length(pids) < 2L) next
    anchor <- pids[[1]]
    for (other in pids[-1]) uf$union(anchor, other)
  }

  clusters <- list()
  for (row in canonical) {
    pids <- as.integer(row$patch_ids)
    if (length(pids) < 2L) next
    root <- as.character(uf$find(pids[[1]]))
    if (is.null(clusters[[root]])) {
      clusters[[root]] <- list(patch_ids = integer(0), candidates = list())
    }
    clusters[[root]]$patch_ids <- unique(c(clusters[[root]]$patch_ids, pids))
    clusters[[root]]$candidates <- c(clusters[[root]]$candidates, list(row))
  }

  out <- list()
  for (root in sort(names(clusters))) {
    patch_ids <- sort(unique(as.integer(clusters[[root]]$patch_ids)))
    patch_index <- stats::setNames(seq_along(patch_ids) - 1L, as.character(patch_ids))
    patch_area_keys <- vapply(patch_ids, function(pid) round(as.numeric(area_lookup[[as.character(pid)]]) * scale), numeric(1))
    cluster_candidates <- list()
    rows <- clusters[[root]]$candidates
    rows <- rows[order(
      -vapply(rows, function(rec) rec$potential_area_key, numeric(1)),
      vapply(rows, function(rec) rec$cost_key, integer(1)),
      vapply(rows, function(rec) rec$length, numeric(1))
    )]
    for (row in rows) {
      local_indices <- sort(as.integer(unname(patch_index[as.character(row$patch_ids)])))
      if (length(local_indices) < 2L) next
      cluster_candidates[[length(cluster_candidates) + 1L]] <- list(
        canon_id = as.integer(row$canon_id),
        candidate_index = as.integer(row$candidate_index),
        patch_ids = as.integer(row$patch_ids),
        local_indices = as.integer(local_indices),
        cost = as.numeric(row$cost),
        cost_key = as.integer(row$cost_key),
        length = as.numeric(row$length),
        potential_area_key = as.integer(sum(as.integer(patch_area_keys[local_indices + 1L])))
      )
    }
    out[[length(out) + 1L]] <- list(
      patch_ids = patch_ids,
      patch_area_keys = patch_area_keys,
      candidates = cluster_candidates
    )
  }
  out
}

terralink_bigconnect_rank_states <- function(states, patch_area_keys, remaining_patch_sets, beam_width) {
  if (length(states) <= beam_width) return(states)
  meta <- lapply(states, function(state) {
    sig <- state$sig
    list(
      optimistic_area = terralink_bigconnect_remaining_upper_bound(sig, patch_area_keys, remaining_patch_sets),
      current_area = terralink_bigconnect_connected_area(sig, patch_area_keys),
      cohesion = terralink_bigconnect_cohesion_key(sig, patch_area_keys),
      spend_key = as.integer(state$budget_used_key),
      count = as.integer(state$corridor_count),
      total_length = as.numeric(state$total_length)
    )
  })
  order_idx <- order(
    -vapply(meta, function(x) x$optimistic_area, numeric(1)),
    -vapply(meta, function(x) x$current_area, numeric(1)),
    -vapply(meta, function(x) x$cohesion, numeric(1)),
    vapply(meta, function(x) x$spend_key, numeric(1)),
    vapply(meta, function(x) x$total_length, numeric(1)),
    vapply(meta, function(x) x$count, numeric(1))
  )
  states[order_idx[seq_len(min(length(order_idx), beam_width))]]
}

terralink_bigconnect_exact_frontier_df <- function(cluster, budget_key_limit, max_states = 250000L) {
  candidates <- cluster$candidates
  patch_area_keys <- cluster$patch_area_keys
  start_sig <- seq_along(patch_area_keys) - 1L
  states <- list()
  states[[terralink_bigconnect_state_key(start_sig, 0L)]] <- list(
    sig = start_sig,
    budget_used_key = 0L,
    total_length = 0,
    corridor_count = 0L,
    selected_canon_ids = integer(0)
  )
  remaining_patch_sets <- lapply(candidates, function(cand) as.integer(cand$local_indices))
  peak_states <- 1L

  for (idx in seq_along(candidates)) {
    cand <- candidates[[idx]]
    next_states <- states
    for (state in states) {
      new_spend <- as.integer(state$budget_used_key) + as.integer(cand$cost_key)
      if (new_spend > budget_key_limit) next
      upper_bound <- terralink_bigconnect_remaining_upper_bound(
        state$sig,
        patch_area_keys,
        remaining_patch_sets[idx:length(remaining_patch_sets)]
      )
      if (upper_bound <= 0L) next
      new_sig <- terralink_bigconnect_union_signature(state$sig, cand$local_indices)
      key <- terralink_bigconnect_state_key(new_sig, new_spend)
      new_state <- list(
        sig = new_sig,
        budget_used_key = new_spend,
        total_length = as.numeric(state$total_length) + as.numeric(cand$length),
        corridor_count = as.integer(state$corridor_count) + 1L,
        selected_canon_ids = c(as.integer(state$selected_canon_ids), as.integer(cand$canon_id))
      )
      prev <- next_states[[key]]
      if (is.null(prev)) {
        next_states[[key]] <- new_state
      } else {
        area_key <- terralink_bigconnect_connected_area(new_sig, patch_area_keys)
        cohesion_key <- terralink_bigconnect_cohesion_key(new_sig, patch_area_keys)
        prev_score <- terralink_bigconnect_score_tuple(area_key, cohesion_key, prev$budget_used_key, prev$corridor_count, prev$total_length)
        new_score <- terralink_bigconnect_score_tuple(area_key, cohesion_key, new_state$budget_used_key, new_state$corridor_count, new_state$total_length)
        if (terralink_bigconnect_score_is_better(new_score, prev_score)) {
          next_states[[key]] <- new_state
        }
      }
    }
    states <- next_states
    peak_states <- max(peak_states, length(states))
    if (length(states) > max_states) stop("state_limit")
  }

  rows <- lapply(states, function(state) {
    list(
      connected_area_key = terralink_bigconnect_connected_area(state$sig, patch_area_keys),
      cohesion_key = terralink_bigconnect_cohesion_key(state$sig, patch_area_keys),
      budget_used_key = as.integer(state$budget_used_key),
      total_length = as.numeric(state$total_length),
      corridor_count = as.integer(state$corridor_count),
      selected_canon_ids = as.integer(state$selected_canon_ids),
      exact_flag = TRUE
    )
  })
  list(frontier = terralink_bigconnect_keep_best_frontier_rows(rows), meta = list(exact = TRUE, peak_states = peak_states, abort_reason = ""))
}

terralink_bigconnect_beam_frontier_df <- function(cluster, budget_key_limit, beam_width = 256L) {
  candidates <- cluster$candidates
  patch_area_keys <- cluster$patch_area_keys
  states <- list()
  start_sig <- seq_along(patch_area_keys) - 1L
  states[[terralink_bigconnect_state_key(start_sig, 0L)]] <- list(
    sig = start_sig,
    budget_used_key = 0L,
    total_length = 0,
    corridor_count = 0L,
    selected_canon_ids = integer(0)
  )
  remaining_patch_sets <- lapply(candidates, function(cand) as.integer(cand$local_indices))
  peak_states <- 1L

  for (idx in seq_along(candidates)) {
    cand <- candidates[[idx]]
    merged <- states
    for (state in states) {
      new_spend <- as.integer(state$budget_used_key) + as.integer(cand$cost_key)
      if (new_spend > budget_key_limit) next
      new_sig <- terralink_bigconnect_union_signature(state$sig, cand$local_indices)
      key <- terralink_bigconnect_state_key(new_sig, new_spend)
      new_state <- list(
        sig = new_sig,
        budget_used_key = new_spend,
        total_length = as.numeric(state$total_length) + as.numeric(cand$length),
        corridor_count = as.integer(state$corridor_count) + 1L,
        selected_canon_ids = c(as.integer(state$selected_canon_ids), as.integer(cand$canon_id))
      )
      prev <- merged[[key]]
      if (is.null(prev)) {
        merged[[key]] <- new_state
      } else {
        area_key <- terralink_bigconnect_connected_area(new_sig, patch_area_keys)
        cohesion_key <- terralink_bigconnect_cohesion_key(new_sig, patch_area_keys)
        prev_score <- terralink_bigconnect_score_tuple(area_key, cohesion_key, prev$budget_used_key, prev$corridor_count, prev$total_length)
        new_score <- terralink_bigconnect_score_tuple(area_key, cohesion_key, new_state$budget_used_key, new_state$corridor_count, new_state$total_length)
        if (terralink_bigconnect_score_is_better(new_score, prev_score)) {
          merged[[key]] <- new_state
        }
      }
    }
    states <- terralink_bigconnect_rank_states(
      merged,
      patch_area_keys = patch_area_keys,
      remaining_patch_sets = if (idx < length(remaining_patch_sets)) remaining_patch_sets[(idx + 1L):length(remaining_patch_sets)] else list(),
      beam_width = beam_width
    )
    peak_states <- max(peak_states, length(states))
  }

  rows <- lapply(states, function(state) {
    list(
      connected_area_key = terralink_bigconnect_connected_area(state$sig, patch_area_keys),
      cohesion_key = terralink_bigconnect_cohesion_key(state$sig, patch_area_keys),
      budget_used_key = as.integer(state$budget_used_key),
      total_length = as.numeric(state$total_length),
      corridor_count = as.integer(state$corridor_count),
      selected_canon_ids = as.integer(state$selected_canon_ids),
      exact_flag = FALSE
    )
  })
  list(frontier = terralink_bigconnect_keep_best_frontier_rows(rows), meta = list(exact = FALSE, peak_states = peak_states, abort_reason = "candidate_limit"))
}

terralink_bigconnect_solve_cluster_df <- function(cluster, budget_key_limit, exact_max_candidates = 26L) {
  candidates <- cluster$candidates
  if (length(candidates) == 0L) {
    return(list(frontier = list(list(connected_area_key = 0L, cohesion_key = 0L, budget_used_key = 0L, total_length = 0, corridor_count = 0L, selected_canon_ids = integer(0), exact_flag = TRUE)), meta = list(exact = TRUE, peak_states = 1L, abort_reason = "")))
  }
  if (length(candidates) <= exact_max_candidates) {
    exact <- tryCatch(
      terralink_bigconnect_exact_frontier_df(cluster, budget_key_limit = budget_key_limit),
      error = function(e) NULL
    )
    if (!is.null(exact)) return(exact)
  }
  terralink_bigconnect_beam_frontier_df(cluster, budget_key_limit = budget_key_limit)
}

terralink_bigconnect_objective_from_selection_df <- function(selected_rows, patch_df, area_scale = 1) {
  valid_node_ids <- suppressWarnings(as.integer(patch_df$patch_id))
  valid_node_ids <- valid_node_ids[is.finite(valid_node_ids)]
  area_lookup <- stats::setNames(as.numeric(patch_df$area), as.character(patch_df$patch_id))
  area_keys <- stats::setNames(vapply(area_lookup, terralink_bigconnect_key, integer(1), scale = area_scale), names(area_lookup))
  uf <- UnionFind$new()
  for (pid in valid_node_ids) uf$find(pid)
  if (length(selected_rows) > 0L) {
    for (row in selected_rows) {
      pids <- terralink_candidate_patch_ids(row, valid_node_ids = valid_node_ids)
      if (length(pids) < 2L) next
      anchor <- pids[[1]]
      for (other in pids[-1]) uf$union(anchor, other)
    }
  }
  groups <- list()
  for (pid in valid_node_ids) {
    root <- as.character(uf$find(pid))
    groups[[root]] <- c(groups[[root]], pid)
  }
  connected_area_key <- 0
  connected_patches <- 0L
  for (members in groups) {
    members <- as.integer(members)
    if (length(members) < 2L) next
    connected_patches <- connected_patches + length(members)
    connected_area_key <- connected_area_key + sum(as.numeric(area_keys[as.character(members)]))
  }
  list(
    connected_area_key = connected_area_key,
    connected_area = as.numeric(connected_area_key) / as.numeric(area_scale),
    connected_patches = as.integer(connected_patches)
  )
}

terralink_bigconnect_refill_df <- function(
  selected_ids,
  candidates,
  patch_df,
  budget,
  cost_col = "cost_metric",
  length_col = "length_metric",
  area_scale = 1
) {
  valid_node_ids <- suppressWarnings(as.integer(patch_df$patch_id))
  valid_node_ids <- valid_node_ids[is.finite(valid_node_ids)]
  total_patch_count <- length(valid_node_ids)
  remaining_budget <- max(as.numeric(budget) - sum(vapply(selected_ids, function(idx) terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(cost_col, "cost", "area"), default = 0), numeric(1))), 0)
  if (remaining_budget <= 1e-12) {
    return(list(selected_ids = as.integer(selected_ids), redundant_added = 0L))
  }

  selected_rows <- lapply(selected_ids, function(idx) candidates[idx, , drop = FALSE])
  objective <- terralink_bigconnect_objective_from_selection_df(selected_rows, patch_df, area_scale = area_scale)
  current_area_key <- objective$connected_area_key
  current_patches <- objective$connected_patches
  redundant_added <- 0L

  repeat {
    feasible_idx <- integer(0)
    for (idx in seq_len(nrow(candidates))) {
      if (idx %in% selected_ids) next
      cost <- terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(cost_col, "cost", "area"), default = 0)
      if (!is.finite(cost) || cost <= 0 || cost > remaining_budget + 1e-12) next
      ids <- terralink_candidate_patch_ids(candidates[idx, , drop = FALSE], valid_node_ids = valid_node_ids)
      if (length(ids) < 2L) next
      feasible_idx <- c(feasible_idx, idx)
    }
    if (length(feasible_idx) == 0L) break

    best_direct <- NULL
    best_direct_score <- NULL
    for (idx in feasible_idx) {
      trial_rows <- c(selected_rows, list(candidates[idx, , drop = FALSE]))
      objective_after <- terralink_bigconnect_objective_from_selection_df(trial_rows, patch_df, area_scale = area_scale)
      area_gain_key <- as.numeric(objective_after$connected_area_key) - as.numeric(current_area_key)
      if (area_gain_key <= 0) next
      patch_gain <- as.integer(objective_after$connected_patches) - as.integer(current_patches)
      cost_key <- terralink_bigconnect_key(
        terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(cost_col, "cost", "area"), default = 0),
        scale = area_scale
      )
      length_val <- terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(length_col, "length", "distance_m", "distance"), default = cost_key)
      score <- c(as.numeric(area_gain_key), as.numeric(patch_gain), -as.numeric(cost_key), -as.numeric(length_val))
      if (terralink_bigconnect_score_is_better(score, best_direct_score)) {
        best_direct_score <- score
        best_direct <- list(idx = idx, objective = objective_after)
      }
    }

    chosen_idx <- NULL
    chosen_objective <- NULL
    if (!is.null(best_direct)) {
      chosen_idx <- best_direct$idx
      chosen_objective <- best_direct$objective
    } else if (current_patches >= total_patch_count) {
      uf <- UnionFind$new()
      for (pid in valid_node_ids) uf$find(pid)
      for (row in selected_rows) {
        ids <- terralink_candidate_patch_ids(row, valid_node_ids = valid_node_ids)
        if (length(ids) < 2L) next
        anchor <- ids[[1]]
        for (other in ids[-1]) uf$union(anchor, other)
      }
      best_redundant <- NULL
      best_score <- NULL
      for (idx in feasible_idx) {
        ids <- terralink_candidate_patch_ids(candidates[idx, , drop = FALSE], valid_node_ids = valid_node_ids)
        roots <- unique(vapply(ids, function(pid) as.integer(uf$find(pid)), integer(1)))
        if (length(roots) != 1L) next
        cost_key <- terralink_bigconnect_key(
          terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(cost_col, "cost", "area"), default = 0),
          scale = area_scale
        )
        length_val <- terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(length_col, "length", "distance_m", "distance"), default = cost_key)
        score <- c(-as.numeric(length_val), -as.numeric(cost_key), -as.numeric(length(ids)))
        if (terralink_bigconnect_score_is_better(score, best_score)) {
          best_score <- score
          best_redundant <- idx
        }
      }
      if (!is.null(best_redundant)) {
        chosen_idx <- best_redundant
        chosen_objective <- objective
        redundant_added <- redundant_added + 1L
      }
    }

    if (is.null(chosen_idx)) break
    selected_ids <- c(selected_ids, as.integer(chosen_idx))
    selected_rows <- c(selected_rows, list(candidates[chosen_idx, , drop = FALSE]))
    remaining_budget <- remaining_budget - terralink_bigconnect_candidate_value(candidates[chosen_idx, , drop = FALSE], c(cost_col, "cost", "area"), default = 0)
    current_area_key <- as.numeric(chosen_objective$connected_area_key)
    current_patches <- as.integer(chosen_objective$connected_patches)
    if (remaining_budget <= 1e-12) break
  }

  list(selected_ids = as.integer(selected_ids), redundant_added = as.integer(redundant_added))
}

terralink_optimize_connected_area_plugin_parity_df <- function(
  candidates,
  patch_df,
  budget,
  cost_col = "cost_metric",
  length_col = "length_metric",
  scale = 1,
  merge_equiv_area = 0,
  merge_equiv_ratio = 0.01
) {
  budget_limit <- max(as.numeric(budget %||% 0), 0)
  if (budget_limit <= 0 || !is.data.frame(candidates) || nrow(candidates) == 0 || !is.data.frame(patch_df) || nrow(patch_df) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "most_connected_habitat", corridors_used = 0L, budget_used = 0)))
  }

  build <- terralink_bigconnect_build_canonical_candidates_df(
    candidates = candidates,
    patch_df = patch_df,
    budget = budget_limit,
    cost_col = cost_col,
    length_col = length_col,
    scale = scale,
    merge_equiv_area = merge_equiv_area,
    merge_equiv_ratio = merge_equiv_ratio
  )
  canonical <- build$canonical
  if (length(canonical) == 0L) {
    return(list(selected = integer(0), stats = list(strategy = "most_connected_habitat", corridors_used = 0L, budget_used = 0)))
  }

  clusters <- terralink_bigconnect_build_clusters_df(canonical, patch_df, scale = scale)
  budget_key_limit <- terralink_bigconnect_key(budget_limit, scale = scale)
  cluster_frontiers <- list()
  exact_clusters <- 0L
  heuristic_clusters <- 0L
  frontier_states_total <- 0L
  for (cluster in clusters) {
    solved <- terralink_bigconnect_solve_cluster_df(cluster, budget_key_limit = budget_key_limit)
    zero_row <- list(
      connected_area_key = 0,
      cohesion_key = 0,
      budget_used_key = 0L,
      total_length = 0,
      corridor_count = 0L,
      selected_canon_ids = integer(0),
      exact_flag = solved$meta$exact
    )
    frontier <- terralink_bigconnect_keep_best_frontier_rows(c(list(zero_row), solved$frontier))
    cluster_frontiers[[length(cluster_frontiers) + 1L]] <- frontier
    frontier_states_total <- frontier_states_total + length(frontier)
    if (isTRUE(solved$meta$exact)) {
      exact_clusters <- exact_clusters + 1L
    } else {
      heuristic_clusters <- heuristic_clusters + 1L
    }
  }

  dp <- list()
  dp[["0"]] <- list(
    connected_area_key = 0,
    cohesion_key = 0,
    budget_used_key = 0L,
    total_length = 0,
    corridor_count = 0L,
    selected_canon_ids = integer(0),
    exact_flag = TRUE
  )
  for (frontier in cluster_frontiers) {
    next_dp <- list()
    for (row_a in dp) {
      for (row_b in frontier) {
        new_spend <- as.integer(row_a$budget_used_key) + as.integer(row_b$budget_used_key)
        if (new_spend > budget_key_limit) next
        combined <- list(
          connected_area_key = as.numeric(row_a$connected_area_key) + as.numeric(row_b$connected_area_key),
          cohesion_key = as.numeric(row_a$cohesion_key) + as.numeric(row_b$cohesion_key),
          budget_used_key = new_spend,
          total_length = as.numeric(row_a$total_length) + as.numeric(row_b$total_length),
          corridor_count = as.integer(row_a$corridor_count) + as.integer(row_b$corridor_count),
          selected_canon_ids = c(as.integer(row_a$selected_canon_ids), as.integer(row_b$selected_canon_ids)),
          exact_flag = isTRUE(row_a$exact_flag) && isTRUE(row_b$exact_flag)
        )
        prev <- next_dp[[as.character(new_spend)]]
        prev_score <- NULL
        if (!is.null(prev)) {
          prev_score <- terralink_bigconnect_score_tuple(prev$connected_area_key, prev$cohesion_key, prev$budget_used_key, prev$corridor_count, prev$total_length)
        }
        new_score <- terralink_bigconnect_score_tuple(combined$connected_area_key, combined$cohesion_key, combined$budget_used_key, combined$corridor_count, combined$total_length)
        if (terralink_bigconnect_score_is_better(new_score, prev_score)) {
          next_dp[[as.character(new_spend)]] <- combined
        }
      }
    }
    dp <- next_dp
  }

  best_row <- NULL
  best_score <- NULL
  for (row in dp) {
    score <- terralink_bigconnect_score_tuple(row$connected_area_key, row$cohesion_key, row$budget_used_key, row$corridor_count, row$total_length)
    if (terralink_bigconnect_score_is_better(score, best_score)) {
      best_score <- score
      best_row <- row
    }
  }
  if (is.null(best_row)) {
    return(list(selected = integer(0), stats = list(strategy = "most_connected_habitat", corridors_used = 0L, budget_used = 0)))
  }

  canon_by_id <- stats::setNames(canonical, vapply(canonical, function(rec) as.character(rec$canon_id), character(1)))
  selected_ids <- integer(0)
  for (canon_id in sort(unique(as.integer(best_row$selected_canon_ids)))) {
    selected_ids <- c(selected_ids, as.integer(canon_by_id[[as.character(canon_id)]]$candidate_index))
  }

  refill <- terralink_bigconnect_refill_df(
    selected_ids = selected_ids,
    candidates = candidates,
    patch_df = patch_df,
    budget = budget_limit,
    cost_col = cost_col,
    length_col = length_col,
    area_scale = scale
  )
  selected_ids <- as.integer(unique(refill$selected_ids))
  selected_rows <- lapply(selected_ids, function(idx) candidates[idx, , drop = FALSE])
  objective <- terralink_bigconnect_objective_from_selection_df(selected_rows, patch_df, area_scale = scale)
  budget_used <- sum(vapply(selected_ids, function(idx) terralink_bigconnect_candidate_value(candidates[idx, , drop = FALSE], c(cost_col, "cost", "area"), default = 0), numeric(1)))

  list(
    selected = as.integer(candidates$id[selected_ids] %||% selected_ids),
    stats = list(
      strategy = "most_connected_habitat",
      corridors_used = length(selected_ids),
      budget_used = as.numeric(budget_used),
      primary_links = as.integer(length(selected_ids) - as.integer(refill$redundant_added %||% 0L)),
      redundant_links = as.integer(refill$redundant_added %||% 0L),
      bigconnect_objective_post = as.numeric(objective$connected_area),
      bigconnect_connected_patches = as.integer(objective$connected_patches),
      bigconnect_candidate_count_raw = as.integer(build$stats$raw_count %||% 0L),
      bigconnect_candidate_count_nondominated = as.integer(build$stats$nondominated_count %||% 0L),
      bigconnect_clusters_total = as.integer(length(clusters)),
      bigconnect_clusters_exact = as.integer(exact_clusters),
      bigconnect_clusters_heuristic = as.integer(heuristic_clusters),
      bigconnect_proven_optimal = heuristic_clusters == 0L,
      bigconnect_frontier_states_total = as.integer(frontier_states_total)
    )
  )
}

terralink_lsn_shortcut_multiplier <- function(g, p1, p2, length_val) {
  if (!is.finite(length_val) || length_val <= 0) return(0.1)
  if (is.null(g) || igraph::ecount(g) == 0) return(0.1)
  dist_m <- tryCatch(
    igraph::distances(g, v = as.character(p1), to = as.character(p2), weights = igraph::E(g)$weight),
    error = function(e) matrix(NA_real_, nrow = 1, ncol = 1)
  )
  current_len <- suppressWarnings(as.numeric(dist_m[[1]]))
  if (!is.finite(current_len)) return(0.1)
  ratio <- current_len / max(length_val, 1e-9)
  if (ratio >= 3.0) return(0.9)
  if (ratio >= 1.5) return(0.5)
  0.1
}

terralink_optimize_largest_single_network_plugin_parity_df <- function(
  candidates,
  patch_df,
  budget,
  cost_col = "cost_metric",
  length_col = "length_metric",
  max_redundant_links = 1L
) {
  budget_limit <- max(as.numeric(budget %||% 0), 0)
  if (budget_limit <= 0 || !is.data.frame(candidates) || nrow(candidates) == 0 || !is.data.frame(patch_df) || nrow(patch_df) < 2) {
    return(list(selected = integer(0), stats = list(strategy = "largest_single_network", corridors_used = 0L, budget_used = 0)))
  }

  valid_node_ids <- suppressWarnings(as.integer(patch_df$patch_id))
  valid_node_ids <- valid_node_ids[is.finite(valid_node_ids)]
  area_lookup <- stats::setNames(as.numeric(patch_df$area), as.character(patch_df$patch_id))
  seed_patch <- as.integer(patch_df$patch_id[[which.max(as.numeric(patch_df$area))]])
  remaining <- budget_limit
  selected_idx <- integer(0)
  selected_types <- character(0)
  redundant_links <- 0L

  uf <- UnionFind$new()
  for (pid in valid_node_ids) {
    uf$find(pid)
    uf$size[[as.character(pid)]] <- as.numeric(area_lookup[[as.character(pid)]] %||% 0)
    uf$count[[as.character(pid)]] <- 1L
  }
  g <- igraph::make_empty_graph(directed = FALSE)
  if (length(valid_node_ids) > 0) {
    g <- igraph::add_vertices(g, length(valid_node_ids), name = as.character(valid_node_ids))
  }

  main_root <- function() as.integer(uf$find(seed_patch))
  touches_main_and_other <- function(ids) {
    mr <- main_root()
    roots <- unique(vapply(ids, function(pid) as.integer(uf$find(pid)), integer(1)))
    any(roots == mr) && any(roots != mr)
  }
  candidate_cost <- function(cand) terralink_bigconnect_candidate_value(cand, c(cost_col, "cost", "area"), default = 0)
  candidate_length <- function(cand, fallback_cost) {
    val <- terralink_bigconnect_candidate_value(cand, c(length_col, "length", "distance_m", "distance"), default = fallback_cost)
    if (is.finite(val) && val > 0) val else fallback_cost
  }
  base_roi <- function(cand, ids, cost) {
    if (!is.finite(cost) || cost <= 0 || length(ids) < 2L) return(0)
    p1 <- as.integer(ids[[1]])
    p2 <- as.integer(ids[[length(ids)]])
    w <- sqrt(max(as.numeric(area_lookup[[as.character(p1)]] %||% 0), 0) * max(as.numeric(area_lookup[[as.character(p2)]] %||% 0), 0))
    if (!is.finite(w) || w <= 0) return(0)
    w / cost
  }
  commit_candidate <- function(idx, corr_type) {
    cand <- candidates[idx, , drop = FALSE]
    cost <- candidate_cost(cand)
    ids <- terralink_candidate_patch_ids(cand, valid_node_ids = valid_node_ids)
    if (!is.finite(cost) || cost <= 0 || cost > remaining + 1e-12 || length(ids) < 2L) return(FALSE)
    anchor <- ids[[1]]
    for (other in ids[-1]) uf$union(anchor, other)
    len_val <- candidate_length(cand, cost)
    g <<- terralink_graph_add_or_update(g, ids, len_val)
    remaining <<- remaining - cost
    selected_idx <<- c(selected_idx, idx)
    selected_types <<- c(selected_types, corr_type)
    TRUE
  }

  repeat {
    best_idx <- NA_integer_
    best_score <- -Inf
    best_gain <- -Inf
    best_cost <- Inf
    best_length <- Inf
    for (idx in seq_len(nrow(candidates))) {
      if (idx %in% selected_idx) next
      cand <- candidates[idx, , drop = FALSE]
      cost <- candidate_cost(cand)
      if (!is.finite(cost) || cost <= 0 || cost > remaining + 1e-12) next
      ids <- terralink_candidate_patch_ids(cand, valid_node_ids = valid_node_ids)
      if (length(ids) < 2L || !touches_main_and_other(ids)) next
      roi <- base_roi(cand, ids, cost)
      if (!is.finite(roi) || roi <= 0) next
      gain <- min(vapply(ids, function(pid) as.numeric(area_lookup[[as.character(pid)]] %||% 0), numeric(1)))
      if (!is.finite(gain) || gain <= 0) gain <- 1
      score <- roi * gain
      len_val <- candidate_length(cand, cost)
      better <- score > best_score + 1e-12 ||
        (abs(score - best_score) <= 1e-12 && gain > best_gain + 1e-12) ||
        (abs(score - best_score) <= 1e-12 && abs(gain - best_gain) <= 1e-12 && cost < best_cost - 1e-12) ||
        (abs(score - best_score) <= 1e-12 && abs(gain - best_gain) <= 1e-12 && abs(cost - best_cost) <= 1e-12 && len_val < best_length - 1e-12)
      if (isTRUE(better)) {
        best_idx <- idx
        best_score <- score
        best_gain <- gain
        best_cost <- cost
        best_length <- len_val
      }
    }
    if (!is.finite(best_idx)) break
    if (!commit_candidate(best_idx, "primary")) break
  }

  repeat {
    best_idx <- NA_integer_
    best_score <- -Inf
    best_cost <- Inf
    best_length <- Inf
    best_type <- "primary"
    for (idx in seq_len(nrow(candidates))) {
      if (idx %in% selected_idx) next
      cand <- candidates[idx, , drop = FALSE]
      cost <- candidate_cost(cand)
      if (!is.finite(cost) || cost <= 0 || cost > remaining + 1e-12) next
      ids <- terralink_candidate_patch_ids(cand, valid_node_ids = valid_node_ids)
      if (length(ids) < 2L) next
      roi <- base_roi(cand, ids, cost)
      if (!is.finite(roi) || roi <= 0) next
      roots <- unique(vapply(ids, function(pid) as.integer(uf$find(pid)), integer(1)))
      len_val <- candidate_length(cand, cost)
      if (length(roots) > 1L) {
        if (!touches_main_and_other(ids)) next
        gain <- min(vapply(ids, function(pid) as.numeric(area_lookup[[as.character(pid)]] %||% 0), numeric(1)))
        if (!is.finite(gain) || gain <= 0) gain <- 1
        mult <- gain
        corr_type <- "primary"
      } else {
        if (as.integer(roots[[1]]) != main_root()) next
        if (redundant_links >= as.integer(max_redundant_links)) next
        mult <- terralink_lsn_shortcut_multiplier(g, ids[[1]], ids[[length(ids)]], len_val)
        if (!is.finite(mult) || mult <= 0.1 + 1e-12) next
        corr_type <- "redundant"
      }
      score <- roi * mult
      better <- score > best_score + 1e-12 ||
        (abs(score - best_score) <= 1e-12 && cost < best_cost - 1e-12) ||
        (abs(score - best_score) <= 1e-12 && abs(cost - best_cost) <= 1e-12 && len_val < best_length - 1e-12)
      if (isTRUE(better)) {
        best_idx <- idx
        best_score <- score
        best_cost <- cost
        best_length <- len_val
        best_type <- corr_type
      }
    }
    if (!is.finite(best_idx)) break
    if (!commit_candidate(best_idx, best_type)) break
    if (identical(best_type, "redundant")) redundant_links <- redundant_links + 1L
  }

  selected_ids <- as.integer(candidates$id[selected_idx] %||% selected_idx)
  budget_used <- sum(vapply(selected_idx, function(idx) candidate_cost(candidates[idx, , drop = FALSE]), numeric(1)))
  list(
    selected = selected_ids,
    stats = list(
      strategy = "largest_single_network",
      corridors_used = length(selected_ids),
      budget_used = as.numeric(budget_used),
      primary_links = as.integer(sum(selected_types == "primary")),
      redundant_links = as.integer(sum(selected_types == "redundant"))
    )
  )
}
