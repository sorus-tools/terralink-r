#' Union-Find data structure
#'
#' @param x Node id for lookup operations.
#' @param a Node id for union operations.
#' @param b Node id for union operations.
#' @field parent Environment mapping nodes to parents.
#' @field size Environment mapping roots to component sizes.
#' @field count Environment mapping roots to component counts.
#' @param ... Unused (included for R6 compatibility).
#' @section Methods:
#' \describe{
#'   \item{initialize}{Create a new UnionFind.}
#'   \item{find}{Find root of a node with path compression.}
#'   \item{union}{Union two nodes; returns TRUE if merged.}
#'   \item{get_size}{Get component size for a node.}
#'   \item{get_count}{Get component count for a node.}
#' }
#' @param x Node id for lookup operations.
#' @param a Node id for union operations.
#' @param b Node id for union operations.
#' @name UnionFind
#' @export
UnionFind <- R6::R6Class(
  "UnionFind",
  public = list(
    parent = NULL,
    size = NULL,
    count = NULL,
    initialize = function() {
      self$parent <- new.env(parent = emptyenv())
      self$size <- new.env(parent = emptyenv())
      self$count <- new.env(parent = emptyenv())
      invisible(self)
    },
    find = function(x) {
      private$ensure(x)
      path <- character(0)
      node <- x
      repeat {
        key <- as.character(node)
        private$ensure(node)
        path <- c(path, key)
        parent <- self$parent[[key]]
        if (identical(parent, node)) {
          break
        }
        node <- parent
      }
      root <- self$parent[[as.character(node)]]
      for (key in path) {
        self$parent[[key]] <- root
      }
      root
    },
    union = function(a, b) {
      ra <- self$find(a)
      rb <- self$find(b)
      if (identical(ra, rb)) {
        return(FALSE)
      }
      key_a <- as.character(ra)
      key_b <- as.character(rb)
      sa <- self$size[[key_a]]
      sb <- self$size[[key_b]]
      if (is.null(sa)) sa <- 0
      if (is.null(sb)) sb <- 0
      if (sa < sb) {
        tmp <- ra
        ra <- rb
        rb <- tmp
        tmp <- sa
        sa <- sb
        sb <- tmp
        key_a <- as.character(ra)
        key_b <- as.character(rb)
      }
      self$parent[[key_b]] <- ra
      self$size[[key_a]] <- sa + sb
      ca <- self$count[[key_a]]
      cb <- self$count[[key_b]]
      if (is.null(ca)) ca <- 0L
      if (is.null(cb)) cb <- 0L
      self$count[[key_a]] <- as.integer(ca + cb)
      TRUE
    },
    get_size = function(x) {
      root <- self$find(x)
      val <- self$size[[as.character(root)]]
      if (is.null(val)) 0 else val
    },
    get_count = function(x) {
      root <- self$find(x)
      val <- self$count[[as.character(root)]]
      if (is.null(val)) 0L else as.integer(val)
    }
  ),
  private = list(
    ensure = function(x) {
      key <- as.character(x)
      if (!exists(key, envir = self$parent, inherits = FALSE)) {
        self$parent[[key]] <- x
        self$size[[key]] <- 0
        self$count[[key]] <- 0L
      }
    }
  )
)

#' Network optimizer for corridor selection
#'
#' Implements a two-phase optimizer: MST backbone, then optional loop additions.
#'
#' @param nodes Named numeric vector of node weights.
#' @param u Candidate edge start node.
#' @param v Candidate edge end node.
#' @param cand_id Candidate edge id.
#' @param cost Candidate edge cost.
#' @param budget Numeric budget for corridor costs.
#' @param loop_fraction Fraction of budget reserved for loops.
#' @param max_redundancy Max redundant edges per component.
#' @field nodes Node weights.
#' @field edges Candidate edge list.
#' @field uf UnionFind instance.
#' @param ... Unused (included for R6 compatibility).
#' @section Methods:
#' \describe{
#'   \item{initialize}{Create a new optimizer with nodes.}
#'   \item{add_candidate}{Add a candidate edge.}
#'   \item{solve}{Run optimization.}
#' }
#' @param nodes Named numeric vector of node weights.
#' @param u Candidate edge start node.
#' @param v Candidate edge end node.
#' @param cand_id Candidate edge id.
#' @param cost Candidate edge cost.
#' @param budget Numeric budget for corridor costs.
#' @param loop_fraction Fraction of budget reserved for loops.
#' @param max_redundancy Max redundant edges per component.
#' @name NetworkOptimizer
#' @export
NetworkOptimizer <- R6::R6Class(
  "NetworkOptimizer",
  public = list(
    nodes = NULL,
    edges = NULL,
    uf = NULL,
    initialize = function(nodes) {
      if (is.null(nodes) || length(nodes) == 0) {
        terralink_abort("nodes must be a named vector or list of node weights.", class = "terralink_error_input")
      }
      self$nodes <- nodes
      self$edges <- list()
      self$uf <- UnionFind$new()
      node_ids <- names(nodes)
      if (is.null(node_ids)) {
        node_ids <- as.character(seq_along(nodes))
      }
      for (i in seq_along(nodes)) {
        node_id <- node_ids[[i]]
        self$uf$parent[[as.character(node_id)]] <- node_id
        self$uf$size[[as.character(node_id)]] <- as.numeric(nodes[[i]])
        self$uf$count[[as.character(node_id)]] <- 1L
      }
      invisible(self)
    },
    add_candidate = function(u, v, cand_id, cost) {
      if (isTRUE(u == v)) {
        return(invisible(FALSE))
      }
      edge <- list(
        u = as.integer(u),
        v = as.integer(v),
        id = as.integer(cand_id),
        cost = as.numeric(cost)
      )
      self$edges[[length(self$edges) + 1]] <- edge
      invisible(TRUE)
    },
    solve = function(budget, loop_fraction = 0.05, max_redundancy = 2) {
      if (!is.numeric(budget) || budget <= 0) {
        terralink_abort("budget must be a positive number.", class = "terralink_error_input")
      }
      if (length(self$edges) == 0) {
        return(list(
          selected = integer(0),
          component_sizes = numeric(0),
          component_counts = integer(0),
          total_cost = 0
        ))
      }
      costs <- vapply(self$edges, function(e) e$cost, numeric(1))
      order_idx <- order(costs)
      edges_sorted <- self$edges[order_idx]

      selected <- integer(0)
      current_cost <- 0

      for (edge in edges_sorted) {
        if (current_cost + edge$cost > budget) {
          break
        }
        if (self$uf$union(edge$u, edge$v)) {
          selected <- c(selected, edge$id)
          current_cost <- current_cost + edge$cost
        }
      }

      redundancy_count <- new.env(parent = emptyenv())
      get_redundancy <- function(root) {
        key <- as.character(root)
        if (!exists(key, envir = redundancy_count, inherits = FALSE)) {
          return(0L)
        }
        as.integer(redundancy_count[[key]])
      }
      set_redundancy <- function(root, value) {
        redundancy_count[[as.character(root)]] <- as.integer(value)
      }

      if (current_cost < budget) {
        loop_thresh <- if (is.null(loop_fraction) || loop_fraction == 0) NA_real_ else budget * loop_fraction
        for (edge in edges_sorted) {
          if (edge$id %in% selected) {
            next
          }
          remaining <- budget - current_cost
          if (edge$cost > remaining) {
            break
          }
          if (!is.na(loop_thresh) && edge$cost > loop_thresh) {
            next
          }
          root_u <- self$uf$find(edge$u)
          root_v <- self$uf$find(edge$v)
          if (identical(root_u, root_v)) {
            count <- get_redundancy(root_u)
            if (count >= max_redundancy) {
              next
            }
            set_redundancy(root_u, count + 1L)
          }
          selected <- c(selected, edge$id)
          current_cost <- current_cost + edge$cost
          self$uf$union(edge$u, edge$v)
        }
      }

      node_ids <- names(self$nodes)
      if (is.null(node_ids)) {
        node_ids <- as.character(seq_along(self$nodes))
      }
      final_sizes <- numeric(0)
      final_counts <- integer(0)
      for (node_id in node_ids) {
        root <- self$uf$find(node_id)
        key <- as.character(root)
        final_sizes[[key]] <- self$uf$size[[key]]
        final_counts[[key]] <- self$uf$count[[key]]
      }

      list(
        selected = as.integer(selected),
        component_sizes = final_sizes,
        component_counts = final_counts,
        total_cost = current_cost
      )
    }
  )
)

#' Optimize a network given nodes and candidate edges
#'
#' @param nodes Named numeric vector of node weights.
#' @param edges Data frame with columns u, v, id, cost.
#' @param budget Numeric budget for corridor costs.
#' @param loop_fraction Fraction of budget reserved for loops.
#' @param max_redundancy Max redundant edges per component.
#' @return List with selected edge ids, component sizes/counts, and total cost.
#' @export
optimize_network <- function(nodes, edges, budget, loop_fraction = 0.05, max_redundancy = 2) {
  if (!is.data.frame(edges)) {
    terralink_abort("edges must be a data frame with columns u, v, id, cost.", class = "terralink_error_input")
  }
  required <- c("u", "v", "id", "cost")
  if (!all(required %in% names(edges))) {
    terralink_abort("edges is missing required columns: u, v, id, cost.", class = "terralink_error_input")
  }
  optimizer <- NetworkOptimizer$new(nodes)
  for (i in seq_len(nrow(edges))) {
    optimizer$add_candidate(edges$u[[i]], edges$v[[i]], edges$id[[i]], edges$cost[[i]])
  }
  optimizer$solve(budget = budget, loop_fraction = loop_fraction, max_redundancy = max_redundancy)
}
