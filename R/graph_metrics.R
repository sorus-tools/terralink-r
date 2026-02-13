#' Build an igraph graph from patch ids and corridor edges
#'
#' @param patches Patch ids (vector) or data frame with column id.
#' @param corridors Data frame with patch1, patch2, and optional distance column.
#' @param distance_col Column name to use for edge distance/weight.
#' @return An igraph graph.
#' @export
build_graph_from_corridors <- function(patches, corridors, distance_col = "distance_m") {
  if (is.data.frame(patches)) {
    patch_ids <- patches$id
  } else {
    patch_ids <- patches
  }
  patch_ids <- as.character(patch_ids)

  if (!is.data.frame(corridors) || nrow(corridors) == 0) {
    g <- igraph::make_empty_graph(n = length(patch_ids), directed = FALSE)
    if (length(patch_ids) > 0) {
      igraph::V(g)$name <- patch_ids
    }
    return(g)
  }

  edge_mat <- as.matrix(corridors[, c("patch1", "patch2")])
  edge_mat <- apply(edge_mat, 2, as.character)
  g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)

  missing <- setdiff(patch_ids, igraph::V(g)$name)
  if (length(missing) > 0) {
    g <- igraph::add_vertices(g, length(missing), name = missing)
  }

  if (!is.null(distance_col) && distance_col %in% names(corridors)) {
    igraph::E(g)$distance <- as.numeric(corridors[[distance_col]])
    igraph::E(g)$weight <- igraph::E(g)$distance
  } else {
    igraph::E(g)$distance <- 1.0
    igraph::E(g)$weight <- 1.0
  }

  g
}

#' Calculate movement entropy for a graph
#'
#' @param graph igraph graph.
#' @param alpha Dispersal kernel parameter.
#' @return Numeric entropy value.
#' @export
calculate_movement_entropy <- function(graph, alpha = 0.002) {
  n <- igraph::vcount(graph)
  if (n == 0) {
    return(0)
  }
  total <- 0
  vertices <- igraph::V(graph)
  for (v in vertices) {
    edges <- igraph::incident(graph, v, mode = "all")
    if (length(edges) == 0) {
      next
    }
    distances <- igraph::E(graph)[edges]$distance
    distances <- ifelse(is.finite(distances), distances, 1.0)
    probs <- exp(-alpha * distances)
    prob_sum <- sum(probs)
    if (prob_sum > 0) {
      norm_probs <- probs / prob_sum
      total <- total - sum(norm_probs * log(norm_probs + 1e-10))
    }
  }
  total / n
}

#' Calculate topology penalty (cycle count)
#'
#' @param graph igraph graph.
#' @return Numeric penalty.
#' @export
calculate_topology_penalty <- function(graph) {
  if (igraph::vcount(graph) == 0) {
    return(0)
  }
  n_edges <- igraph::ecount(graph)
  n_nodes <- igraph::vcount(graph)
  n_components <- igraph::components(graph)$no
  min_edges_needed <- n_nodes - n_components
  max(0, n_edges - min_edges_needed)
}

#' Calculate disturbance penalty (diameter normalized)
#'
#' @param graph igraph graph.
#' @return Numeric penalty.
#' @export
calculate_disturbance_penalty <- function(graph) {
  n <- igraph::vcount(graph)
  if (n < 2) {
    return(0)
  }
  if (igraph::is_connected(graph)) {
    d <- igraph::diameter(graph, weights = NA)
    return(d / n)
  }
  1.0
}

#' Calculate total entropy summary
#'
#' @param graph igraph graph.
#' @param lambda_c Connectivity penalty multiplier.
#' @param lambda_f Topology penalty multiplier.
#' @param lambda_d Disturbance penalty multiplier.
#' @return Named list of entropy components.
#' @export
calculate_total_entropy <- function(graph, lambda_c = 1.0, lambda_f = 1.0, lambda_d = 1.0) {
  h_mov <- calculate_movement_entropy(graph)
  n_comp <- igraph::components(graph)$no
  c_val <- if (n_comp > 0) (n_comp - 1) else 0
  f_val <- calculate_topology_penalty(graph)
  d_val <- calculate_disturbance_penalty(graph)
  h_total <- h_mov + (lambda_c * c_val) + (lambda_f * f_val) + (lambda_d * d_val)

  list(
    H_total = h_total,
    H_mov = h_mov,
    C_connectivity = c_val,
    F_topology = f_val,
    D_disturbance = d_val
  )
}

#' Calculate the fraction of node pairs with two edge-disjoint paths
#'
#' @param graph igraph graph.
#' @return Numeric ratio.
#' @export
calculate_two_edge_connectivity <- function(graph) {
  n <- igraph::vcount(graph)
  if (n < 3) {
    return(0)
  }
  verts <- igraph::V(graph)
  total_pairs <- n * (n - 1) / 2
  two_connected <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      k <- igraph::edge_connectivity(graph, source = verts[i], target = verts[j])
      if (is.finite(k) && k >= 2) {
        two_connected <- two_connected + 1
      }
    }
  }
  two_connected / total_pairs
}

#' Estimate failure probability via edge removal
#'
#' @param graph igraph graph.
#' @param k_failures Number of edges removed each trial.
#' @param iterations Number of Monte Carlo iterations.
#' @return Failure probability.
#' @export
calculate_failure_probability <- function(graph, k_failures = 1, iterations = 100) {
  edges <- igraph::E(graph)
  m <- length(edges)
  if (m <= k_failures) {
    return(1.0)
  }
  base_components <- igraph::components(graph)$no
  failures <- 0
  for (i in seq_len(iterations)) {
    remove_ids <- sample(seq_len(m), k_failures, replace = FALSE)
    g_test <- igraph::delete_edges(graph, edges[remove_ids])
    if (!igraph::is_connected(g_test) && igraph::components(g_test)$no > base_components) {
      failures <- failures + 1
    }
  }
  failures / iterations
}

#' Score a loop edge for shortcut value
#'
#' @param graph igraph graph.
#' @param u First node id.
#' @param v Second node id.
#' @param weight Edge cost.
#' @return Numeric score.
#' @export
score_edge_for_loops <- function(graph, u, v, weight) {
  if (weight <= 0) {
    return(0)
  }
  if (!igraph::are.connected(graph, as.character(u), as.character(v))) {
    return(0)
  }
  path <- igraph::shortest_paths(graph, from = as.character(u), to = as.character(v), weights = igraph::E(graph)$weight)
  if (length(path$vpath[[1]]) == 0) {
    return(0)
  }
  length_val <- igraph::distances(graph, v = as.character(u), to = as.character(v), weights = igraph::E(graph)$weight)
  if (!is.finite(length_val)) {
    return(0)
  }
  as.numeric(length_val) / (weight + 1e-6)
}
