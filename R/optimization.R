# Optimization strategies for corridor selection

terralink_heap_new <- function(capacity = 1024L) {
  cap <- as.integer(max(16L, capacity))
  h <- new.env(parent = emptyenv())
  h$score <- rep(0, cap)
  h$seq <- rep(0L, cap)
  h$stamp <- rep(0L, cap)
  h$idx <- rep(0L, cap)
  h$n <- 0L
  h
}

terralink_heap_greater <- function(h, i, j) {
  si <- h$score[[i]]
  sj <- h$score[[j]]
  if (si > sj) return(TRUE)
  if (si < sj) return(FALSE)
  h$seq[[i]] < h$seq[[j]]
}

terralink_heap_swap <- function(h, i, j) {
  tmp <- h$score[[i]]
  h$score[[i]] <- h$score[[j]]
  h$score[[j]] <- tmp

  tmp <- h$seq[[i]]
  h$seq[[i]] <- h$seq[[j]]
  h$seq[[j]] <- tmp

  tmp <- h$stamp[[i]]
  h$stamp[[i]] <- h$stamp[[j]]
  h$stamp[[j]] <- tmp

  tmp <- h$idx[[i]]
  h$idx[[i]] <- h$idx[[j]]
  h$idx[[j]] <- tmp
}

terralink_heap_ensure_capacity <- function(h) {
  if (h$n < length(h$score)) return(invisible(h))
  new_cap <- as.integer(max(16L, length(h$score) * 2L))
  length(h$score) <- new_cap
  length(h$seq) <- new_cap
  length(h$stamp) <- new_cap
  length(h$idx) <- new_cap
  invisible(h)
}

terralink_heap_push <- function(h, score, seq, stamp, idx) {
  terralink_heap_ensure_capacity(h)
  h$n <- h$n + 1L
  i <- h$n
  h$score[[i]] <- as.numeric(score)
  h$seq[[i]] <- as.integer(seq)
  h$stamp[[i]] <- as.integer(stamp)
  h$idx[[i]] <- as.integer(idx)

  while (i > 1L) {
    p <- i %/% 2L
    if (!terralink_heap_greater(h, i, p)) break
    terralink_heap_swap(h, i, p)
    i <- p
  }
  invisible(h)
}

terralink_heap_pop <- function(h) {
  if (h$n <= 0L) return(NULL)
  out <- list(
    score = h$score[[1L]],
    seq = h$seq[[1L]],
    stamp = h$stamp[[1L]],
    idx = h$idx[[1L]]
  )
  if (h$n == 1L) {
    h$n <- 0L
    return(out)
  }
  h$score[[1L]] <- h$score[[h$n]]
  h$seq[[1L]] <- h$seq[[h$n]]
  h$stamp[[1L]] <- h$stamp[[h$n]]
  h$idx[[1L]] <- h$idx[[h$n]]
  h$n <- h$n - 1L

  i <- 1L
  repeat {
    left <- i * 2L
    right <- left + 1L
    if (left > h$n) break
    best <- left
    if (right <= h$n && terralink_heap_greater(h, right, left)) {
      best <- right
    }
    if (terralink_heap_greater(h, i, best)) break
    terralink_heap_swap(h, i, best)
    i <- best
  }
  out
}

terralink_pair_key_chr <- function(key) {
  if (length(key) == 0) return("0_0")
  vals <- suppressWarnings(as.integer(key))
  if (all(is.finite(vals))) {
    vals <- sort(vals)
    return(paste(vals, collapse = "_"))
  }
  txt <- sort(as.character(key))
  paste(txt, collapse = "_")
}

terralink_candidate_id <- function(cand, fallback_id) {
  if ("id" %in% names(cand)) {
    idv <- suppressWarnings(as.integer(cand$id[[1]]))
    if (is.finite(idv)) return(idv)
  }
  as.integer(fallback_id)
}

terralink_graph_shortcut_multiplier <- function(
  g,
  p1,
  p2,
  length_val,
  diminishing_base,
  shortcut_ratio_high,
  shortcut_ratio_mid,
  shortcut_ratio_low,
  shortcut_mult_high,
  shortcut_mult_mid,
  shortcut_mult_low
) {
  if (!is.finite(length_val) || length_val <= 0) return(as.numeric(shortcut_mult_low))
  if (is.null(g) || igraph::vcount(g) == 0) return(as.numeric(diminishing_base))
  dist_m <- tryCatch(
    igraph::distances(g, v = as.character(p1), to = as.character(p2), weights = igraph::E(g)$weight),
    error = function(e) matrix(NA_real_, nrow = 1, ncol = 1)
  )
  sp_len <- suppressWarnings(as.numeric(dist_m[[1]]))
  if (!is.finite(sp_len)) return(as.numeric(diminishing_base))
  ratio <- sp_len / max(length_val, 1e-9)
  if (ratio >= shortcut_ratio_high) return(as.numeric(shortcut_mult_high))
  if (ratio >= shortcut_ratio_mid) return(as.numeric(shortcut_mult_mid))
  if (ratio <= shortcut_ratio_low) return(as.numeric(shortcut_mult_low))
  as.numeric(shortcut_mult_low)
}

terralink_graph_add_or_update <- function(g, pids, weight) {
  if (is.null(g) || length(pids) < 2) return(g)
  anchor <- as.character(pids[[1]])
  for (other in pids[-1]) {
    other_chr <- as.character(other)
    has_edge <- tryCatch(igraph::are_adjacent(g, anchor, other_chr), error = function(e) FALSE)
    if (isTRUE(has_edge)) {
      eid <- tryCatch(igraph::get.edge.ids(g, c(anchor, other_chr)), error = function(e) 0L)
      if (is.finite(eid) && eid > 0) {
        igraph::E(g)[eid]$weight <- as.numeric(weight)
      }
    } else {
      g <- igraph::add_edges(g, c(anchor, other_chr), attr = list(weight = as.numeric(weight)))
    }
  }
  g
}

#' Select corridors for the "Most Connectivity" strategy
#'
#' Greedy ROI-based selector with dynamic rescoring, bridge seeding, and optional overlap checks.
#'
#' @param candidates Iterable of candidate objects (data.frame rows or lists).
#' @param budget Total corridor budget.
#' @param get_patch_ids Function that returns patch ids for a candidate.
#' @param get_pair_key Function that returns a sorted pair key for a candidate.
#' @param get_cost Function that returns candidate cost.
#' @param get_base_roi Function that returns candidate base ROI.
#' @param get_length Function that returns candidate length for shortcut scoring.
#' @param get_patch_size Function that returns patch size by id.
#' @param overlap_ratio Function that returns overlap ratio vs prior objects.
#' @param overlap_obj Function that returns overlap object representation.
#' @param redundancy_distance_ok Optional callback that can reject near-duplicate redundant corridors.
#' @param overlap_reject_ratio Overlap ratio threshold for heavy redundancy penalty.
#' @param max_prior_per_pair Maximum overlap objects retained per patch pair.
#' @param diminishing_base Base for redundancy penalty when no shortcut context is available.
#' @param max_links_per_pair Optional hard limit of selected corridors per patch pair.
#' @param enable_bridge_pairs Whether to pre-seed bridge corridor pairs.
#' @param bridge_max_per_patch Max candidates retained per bridge midpoint patch.
#' @param shortcut_ratio_high High shortcut ratio threshold.
#' @param shortcut_ratio_mid Mid shortcut ratio threshold.
#' @param shortcut_ratio_low Low shortcut ratio threshold.
#' @param shortcut_mult_high Multiplier when shortcut ratio is high.
#' @param shortcut_mult_mid Multiplier when shortcut ratio is mid.
#' @param shortcut_mult_low Multiplier when shortcut ratio is low.
#' @return List with picks, selected_ids, and summary stats.
#' @export
select_circuit_utility <- function(
  candidates,
  budget,
  get_patch_ids,
  get_pair_key,
  get_cost,
  get_base_roi,
  get_length,
  get_patch_size,
  overlap_ratio,
  overlap_obj,
  redundancy_distance_ok = NULL,
  overlap_reject_ratio = 0.30,
  max_prior_per_pair = 3,
  diminishing_base = 0.5,
  max_links_per_pair = Inf,
  enable_bridge_pairs = TRUE,
  bridge_max_per_patch = 25,
  shortcut_ratio_high = 3.0,
  shortcut_ratio_mid = 1.5,
  shortcut_ratio_low = 1.5,
  shortcut_mult_high = 0.9,
  shortcut_mult_mid = 0.5,
  shortcut_mult_low = 0.1
) {
  budget_total <- as.numeric(budget)
  budget_guard <- max(1e-3, abs(budget_total) * 1e-12)
  remaining <- budget_total - budget_guard
  if (is.na(remaining) || remaining <= 0) {
    return(list(
      picks = list(),
      selected_ids = integer(0),
      stats = list(
        budget_used = 0,
        budget_remaining = max(0, budget_total),
        primary_links = 0L,
        redundant_links = 0L,
        wasteful_links = 0L,
        corridors_used = 0L
      )
    ))
  }

  cand_list <- if (is.data.frame(candidates)) {
    split(candidates, seq_len(nrow(candidates)))
  } else {
    as.list(candidates)
  }

  base_roi <- vapply(cand_list, get_base_roi, numeric(1))
  cost_vec <- vapply(cand_list, get_cost, numeric(1))
  valid <- is.finite(base_roi) & is.finite(cost_vec) & base_roi > 0 & cost_vec > 0
  cand_list <- cand_list[valid]
  base_roi <- base_roi[valid]
  cost_vec <- cost_vec[valid]

  if (length(cand_list) == 0) {
    return(list(
      picks = list(),
      selected_ids = integer(0),
      stats = list(
        budget_used = 0,
        budget_remaining = remaining,
        primary_links = 0L,
        redundant_links = 0L,
        wasteful_links = 0L,
        corridors_used = 0L
      )
    ))
  }

  patch_ids_all <- unique(unlist(lapply(cand_list, get_patch_ids)))
  patch_ids_all <- suppressWarnings(as.integer(patch_ids_all))
  patch_ids_all <- patch_ids_all[is.finite(patch_ids_all)]
  uf <- UnionFind$new()
  for (pid in patch_ids_all) {
    uf$find(pid)
    try({
      uf$size[[as.character(pid)]] <- as.numeric(get_patch_size(pid))
      uf$count[[as.character(pid)]] <- 1L
    }, silent = TRUE)
  }

  g <- igraph::make_empty_graph(directed = FALSE)
  if (length(patch_ids_all) > 0) {
    g <- igraph::add_vertices(g, length(patch_ids_all), name = as.character(patch_ids_all))
  }

  selected_overlap_by_pair <- list()
  selected_count_by_pair <- list()
  selected_overlap_by_component <- list()
  selected_overlap_global <- list()
  picks <- list()
  selected_ids <- integer(0)
  primary_links <- 0L
  redundant_links <- 0L
  wasteful_links <- 0L

  n <- length(cand_list)
  selected_mask <- rep(FALSE, n)
  stamps <- integer(n)
  heap <- terralink_heap_new(capacity = n + 16L)
  counter <- 0L

  for (i in seq_len(n)) {
    counter <- counter + 1L
    stamps[[i]] <- stamps[[i]] + 1L
    terralink_heap_push(heap, base_roi[[i]], counter, stamps[[i]], i)
  }

  commit_candidate <- function(i, corr_type, score, overlap_r = 0) {
    if (isTRUE(selected_mask[[i]])) return(FALSE)
    cand <- cand_list[[i]]
    cost <- as.numeric(cost_vec[[i]])
    if (!is.finite(cost) || cost <= 0 || cost > remaining) return(FALSE)
    pids <- suppressWarnings(as.integer(get_patch_ids(cand)))
    pids <- pids[is.finite(pids)]
    if (length(pids) < 2) return(FALSE)

    roots_before <- unique(vapply(pids, function(pid) suppressWarnings(as.integer(uf$find(pid))), integer(1)))
    roots_before <- roots_before[is.finite(roots_before)]

    for (pid in pids) uf$find(pid)
    anchor <- pids[[1]]
    for (other in pids[-1]) {
      uf$union(anchor, other)
    }

    remaining <<- remaining - cost
    selected_mask[[i]] <<- TRUE
    selected_ids <<- c(selected_ids, terralink_candidate_id(cand, i))
    picks[[length(picks) + 1L]] <<- list(
      candidate = cand,
      corr_type = corr_type,
      score = as.numeric(score),
      overlap_ratio = as.numeric(overlap_r)
    )

    pair_key_chr <- terralink_pair_key_chr(get_pair_key(cand))
    obj <- tryCatch(overlap_obj(cand), error = function(e) NULL)
    if (!is.null(obj)) {
      prior_pair <- selected_overlap_by_pair[[pair_key_chr]] %||% list()
      selected_overlap_by_pair[[pair_key_chr]] <<- c(prior_pair, list(obj))
      if (length(selected_overlap_by_pair[[pair_key_chr]]) > as.integer(max_prior_per_pair)) {
        selected_overlap_by_pair[[pair_key_chr]] <<- utils::tail(selected_overlap_by_pair[[pair_key_chr]], as.integer(max_prior_per_pair))
      }
      selected_overlap_global <<- c(selected_overlap_global, list(obj))

      merged <- list(obj)
      for (root in roots_before) {
        root_chr <- as.character(root)
        if (!is.null(selected_overlap_by_component[[root_chr]])) {
          merged <- c(merged, selected_overlap_by_component[[root_chr]])
        }
      }
      new_root <- suppressWarnings(as.integer(uf$find(anchor)))
      new_root_chr <- as.character(new_root)
      selected_overlap_by_component[[new_root_chr]] <<- merged
      for (root in roots_before) {
        root_chr <- as.character(root)
        if (root_chr != new_root_chr) {
          selected_overlap_by_component[[root_chr]] <<- NULL
        }
      }
    }
    selected_count_by_pair[[pair_key_chr]] <<- (selected_count_by_pair[[pair_key_chr]] %||% 0L) + 1L

    len_val <- suppressWarnings(as.numeric(get_length(cand)))
    if (!is.finite(len_val) || len_val <= 0) len_val <- cost
    g <<- terralink_graph_add_or_update(g, pids, len_val)
    TRUE
  }

  if (isTRUE(enable_bridge_pairs) && remaining > 0) {
    bridge_by_mid <- list()
    for (i in seq_len(n)) {
      pids <- suppressWarnings(as.integer(get_patch_ids(cand_list[[i]])))
      pids <- pids[is.finite(pids)]
      if (length(pids) < 2) next
      base <- as.numeric(base_roi[[i]])
      cost <- as.numeric(cost_vec[[i]])
      if (!is.finite(base) || !is.finite(cost) || base <= 0 || cost <= 0) next
      p1 <- pids[[1]]
      p2 <- pids[[2]]
      bridge_by_mid[[as.character(p1)]] <- c(bridge_by_mid[[as.character(p1)]], list(list(base = base, idx = i, other = p2, cost = cost)))
      bridge_by_mid[[as.character(p2)]] <- c(bridge_by_mid[[as.character(p2)]], list(list(base = base, idx = i, other = p1, cost = cost)))
    }

    if (length(bridge_by_mid) > 0) {
      for (mid_key in names(bridge_by_mid)) {
        items <- bridge_by_mid[[mid_key]]
        if (length(items) <= bridge_max_per_patch) next
        base_vals <- vapply(items, function(x) x$base, numeric(1))
        keep_idx <- order(base_vals, decreasing = TRUE)[seq_len(bridge_max_per_patch)]
        bridge_by_mid[[mid_key]] <- items[keep_idx]
      }

      bridge_pairs <- list()
      for (mid_key in names(bridge_by_mid)) {
        items <- bridge_by_mid[[mid_key]]
        if (length(items) < 2) next
        for (ii in seq_len(length(items) - 1L)) {
          for (jj in (ii + 1L):length(items)) {
            a <- items[[ii]]
            b <- items[[jj]]
            if (a$other == b$other) next
            total_cost <- as.numeric(a$cost + b$cost)
            if (!is.finite(total_cost) || total_cost <= 0 || total_cost > remaining) next
            s1 <- as.numeric(get_patch_size(a$other))
            s2 <- as.numeric(get_patch_size(b$other))
            if (!is.finite(s1) || !is.finite(s2) || s1 <= 0 || s2 <= 0) next
            score <- sqrt(s1 * s2) / total_cost
            bridge_pairs[[length(bridge_pairs) + 1L]] <- list(score = score, total_cost = total_cost, idx1 = a$idx, idx2 = b$idx, a = a$other, c = b$other)
          }
        }
      }

      if (length(bridge_pairs) > 0) {
        bridge_scores <- vapply(bridge_pairs, function(x) x$score, numeric(1))
        bridge_pairs <- bridge_pairs[order(bridge_scores, decreasing = TRUE)]
        for (bp in bridge_pairs) {
          if (bp$total_cost > remaining) next
          if (bp$a == bp$c) next
          root_a <- suppressWarnings(as.integer(uf$find(bp$a)))
          root_c <- suppressWarnings(as.integer(uf$find(bp$c)))
          if (is.finite(root_a) && is.finite(root_c) && root_a == root_c) next
          if (isTRUE(selected_mask[[bp$idx1]]) || isTRUE(selected_mask[[bp$idx2]])) next
          ok1 <- commit_candidate(bp$idx1, corr_type = "primary", score = bp$score, overlap_r = 0)
          ok2 <- commit_candidate(bp$idx2, corr_type = "primary", score = bp$score, overlap_r = 0)
          if (ok1) primary_links <- primary_links + 1L
          if (ok2) primary_links <- primary_links + 1L
          if (remaining <= 0) break
        }
      }
    }
  }

  repeat {
    if (heap$n <= 0L || remaining <= 0) break
    entry <- terralink_heap_pop(heap)
    if (is.null(entry)) break
    i <- as.integer(entry$idx)
    if (!is.finite(i) || i < 1 || i > n) next
    if (isTRUE(selected_mask[[i]])) next
    if (stamps[[i]] != as.integer(entry$stamp)) next

    cand <- cand_list[[i]]
    cost <- as.numeric(cost_vec[[i]])
    if (!is.finite(cost) || cost <= 0 || cost > remaining) next
    base <- as.numeric(base_roi[[i]])
    if (!is.finite(base) || base <= 0) next
    pids <- suppressWarnings(as.integer(get_patch_ids(cand)))
    pids <- pids[is.finite(pids)]
    if (length(pids) < 2) next

    roots <- unique(vapply(pids, function(pid) suppressWarnings(as.integer(uf$find(pid))), integer(1)))
    roots <- roots[is.finite(roots)]
    pair_key_chr <- terralink_pair_key_chr(get_pair_key(cand))
    selected_for_pair <- selected_count_by_pair[[pair_key_chr]] %||% 0L
    if (is.finite(max_links_per_pair) && selected_for_pair >= max_links_per_pair) next

    overlap_r <- 0
    if (length(roots) > 1) {
      root_sizes <- vapply(roots, function(root) {
        val <- uf$size[[as.character(root)]]
        if (is.null(val) || !is.finite(val)) 0 else as.numeric(val)
      }, numeric(1))
      target_importance <- suppressWarnings(min(root_sizes, na.rm = TRUE))
      if (!is.finite(target_importance) || target_importance <= 0) target_importance <- 1.0
      multiplier <- target_importance
    } else {
      if (!is.null(redundancy_distance_ok) && length(selected_overlap_global) > 0) {
        ok <- tryCatch(isTRUE(redundancy_distance_ok(cand, selected_overlap_global)), error = function(e) TRUE)
        if (!ok) next
      }
      prior <- selected_overlap_by_pair[[pair_key_chr]] %||% list()
      overlap_r <- suppressWarnings(as.numeric(tryCatch(overlap_ratio(cand, prior), error = function(e) 0)))
      if (!is.finite(overlap_r)) overlap_r <- 0
      if (overlap_r > overlap_reject_ratio) {
        multiplier <- 0.01
      } else {
        len_val <- suppressWarnings(as.numeric(get_length(cand)))
        if (!is.finite(len_val) || len_val <= 0) len_val <- cost
        multiplier <- terralink_graph_shortcut_multiplier(
          g = g,
          p1 = pids[[1]],
          p2 = pids[[length(pids)]],
          length_val = len_val,
          diminishing_base = diminishing_base,
          shortcut_ratio_high = shortcut_ratio_high,
          shortcut_ratio_mid = shortcut_ratio_mid,
          shortcut_ratio_low = shortcut_ratio_low,
          shortcut_mult_high = shortcut_mult_high,
          shortcut_mult_mid = shortcut_mult_mid,
          shortcut_mult_low = shortcut_mult_low
        )
      }
    }

    new_score <- base * multiplier
    if (!is.finite(new_score) || new_score <= 0) next
    old_score <- as.numeric(entry$score)
    if (abs(new_score - old_score) > 1e-12) {
      counter <- counter + 1L
      stamps[[i]] <- stamps[[i]] + 1L
      terralink_heap_push(heap, new_score, counter, stamps[[i]], i)
      next
    }

    corr_type <- if (length(roots) > 1) "primary" else "redundant"
    if (corr_type == "primary") {
      primary_links <- primary_links + 1L
    } else {
      redundant_links <- redundant_links + 1L
      if (multiplier <= shortcut_mult_low || multiplier <= 0.01) {
        corr_type <- "wasteful"
        wasteful_links <- wasteful_links + 1L
      }
    }

    if (!commit_candidate(i, corr_type = corr_type, score = new_score, overlap_r = overlap_r)) next
    if (remaining <= 0) break
  }

  budget_used <- budget_total - remaining - budget_guard
  if (!is.finite(budget_used)) budget_used <- 0
  if (budget_used < 0) budget_used <- 0
  budget_remaining_out <- remaining + budget_guard
  if (!is.finite(budget_remaining_out) || budget_remaining_out < 0) budget_remaining_out <- 0

  list(
    picks = picks,
    selected_ids = as.integer(selected_ids),
    stats = list(
      budget_used = budget_used,
      budget_remaining = budget_remaining_out,
      primary_links = as.integer(primary_links),
      redundant_links = as.integer(redundant_links),
      wasteful_links = as.integer(wasteful_links),
      corridors_used = length(picks)
    )
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Optimize for largest connected network (MST backbone + loops)
#'
#' @param nodes Named numeric vector of patch sizes.
#' @param edges Data frame with u, v, id, cost columns.
#' @param budget Numeric budget for corridor cost.
#' @param loop_fraction Fraction of budget reserved for loops.
#' @param max_redundancy Max redundant edges per component.
#' @return List with selected ids and summary.
#' @export
optimize_largest_network <- function(nodes, edges, budget, loop_fraction = 0.05, max_redundancy = 2) {
  optimize_network(nodes, edges, budget, loop_fraction = loop_fraction, max_redundancy = max_redundancy)
}

#' Choose optimization strategy
#'
#' @param strategy Strategy name.
#' @param nodes Named numeric vector of patch sizes.
#' @param edges Data frame with u, v, id, cost.
#' @param candidates Candidate list (for circuit utility).
#' @param budget Numeric budget.
#' @param ... Additional args forwarded to circuit utility selector.
#' @return List with selected ids and stats.
#' @export
optimize_strategy <- function(strategy, nodes, edges, candidates, budget, ...) {
  strategy <- match.arg(strategy, c("largest_network", "most_connectivity", "circuit_utility"))
  if (strategy == "most_connectivity") strategy <- "circuit_utility"
  if (strategy == "largest_network") {
    res <- optimize_largest_network(nodes, edges, budget)
    return(list(selected = res$selected, stats = list(budget_used = res$total_cost)))
  }
  res <- select_circuit_utility(candidates, budget = budget, ...)
  list(selected = res$selected_ids, stats = res$stats, picks = res$picks)
}
