terralink_candidate_col <- function(cand, name) {
  if (!name %in% names(cand)) return(NULL)
  val <- cand[[name]]
  if (is.list(val) && length(val) == 1) {
    inner <- val[[1]]
    if (is.atomic(inner) || inherits(inner, "sfc")) return(inner)
  }
  val
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

#' Run TerraLink raster workflow
#'
#' @param raster SpatRaster or path to raster.
#' @param patch_values Numeric values representing habitat.
#' @param patch_ranges Optional list of value ranges defining habitat.
#' @param budget Total corridor budget (units defined by `units`).
#' @param budget_pixels Back-compat alias for budget (pixels).
#' @param strategy Strategy name: "circuit_utility", "most_connectivity" (alias), or "largest_network".
#' @param min_patch_size Minimum patch size (units defined by `units`).
#' @param min_corridor_width Minimum corridor width (units defined by `units`).
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
#' @param keep_candidates Whether to keep candidate list in the output.
#' @return List with patches, corridors, rasters, and summary.
#' @export
run_raster_analysis <- function(
  raster,
  patch_values,
  budget = NULL,
  budget_pixels = NULL,
  strategy = "circuit_utility",
  min_patch_size = 10,
  min_corridor_width = 3,
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
  keep_candidates = FALSE
) {
  raster <- terralink_resolve_raster(raster)
  ctx <- terralink_new_run_context(verbose = verbose, progress = progress)
  terralink_progress_start(ctx, message = "Starting raster analysis")
  t_start <- proc.time()[[3]]
  preflight <- terralink_preflight_raster(raster, allow_large = allow_large, ctx = ctx)
  strategy_key <- match.arg(tolower(strategy), c("most_connectivity", "largest_network", "circuit_utility"))
  if (strategy_key == "most_connectivity") strategy_key <- "circuit_utility"

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
      max_search_distance = max_search_distance,
      patch_connectivity = patch_connectivity,
      units = units
    )
    result <- list(
      patches = labels,
      patch_table = patch_df,
      corridors = candidates,
      corridor_raster = build_corridor_raster(labels, patch_df, candidates),
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

  if (strategy_key == "largest_network") {
    opt <- optimize_largest_network(nodes, engine_edges, budget = budget_px)
    selected_ids <- opt$selected
    budget_used <- opt$total_cost
  } else {
    opt <- optimize_strategy(
      strategy = "circuit_utility",
      nodes = nodes,
      edges = engine_edges,
      candidates = candidates,
      budget = budget_px,
      get_patch_ids = function(cand) c(cand$patch1, cand$patch2),
      get_pair_key = function(cand) sort(c(cand$patch1, cand$patch2)),
      get_cost = function(cand) cand$cost,
      get_base_roi = function(cand) cand$roi,
      get_length = function(cand) cand$length,
      get_patch_size = function(pid) nodes[[as.character(pid)]],
      enable_bridge_pairs = TRUE,
      overlap_ratio = function(cand, prior) {
        new_buf <- terralink_candidate_col(cand, "buffered_cells")
        new_buf <- as.integer(new_buf)
        if (length(new_buf) == 0 || length(prior) == 0) return(0)
        denom <- max(1, length(new_buf))
        ratios <- vapply(prior, function(p) {
          prev <- if (is.list(p) && !is.null(p$cells)) as.integer(p$cells) else as.integer(p)
          if (length(prev) == 0) return(0)
          length(intersect(new_buf, prev)) / denom
        }, numeric(1))
        if (length(ratios) == 0) return(0)
        suppressWarnings(max(ratios, na.rm = TRUE))
      },
      overlap_obj = function(cand) {
        list(
          cells = as.integer(terralink_candidate_col(cand, "buffered_cells_expanded")),
          line = terralink_candidate_col(cand, "line")
        )
      },
      redundancy_distance_ok = function(cand, prior) {
        if (max_search_distance_px <= 0 || length(prior) == 0) return(TRUE)
        FALSE
      }
    )
    selected_ids <- opt$selected
    budget_used <- opt$stats$budget_used %||% 0
  }

  selected_ids <- as.integer(unique(selected_ids))
  remaining_budget <- as.numeric(max(0, budget_px - budget_used))
  if (remaining_budget > 0 && nrow(candidates) > length(selected_ids)) {
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
  if (remaining_budget > 0 && nrow(selected) > 0) {
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
    selected$connected_size <- terralink_compute_raster_connected_sizes(selected, nodes)
  }

  terralink_inform(sprintf("Corridors selected: %s", nrow(selected)), ctx = ctx, level = 1)
  terralink_inform(sprintf("Budget used (px): %s", terralink_format_number(budget_used, 2)), ctx = ctx, level = 1)

  terralink_progress_update(ctx, 80, "Building outputs")
  pre_mask <- terra::ifel(!is.na(labels) & labels > 0, 1, NA)
  corridor_raster <- build_corridor_raster(labels, patch_df, selected, min_corridor_width_px)
  contiguous_raster <- build_contiguous_raster(pre_mask, corridor_raster, patch_connectivity)

  terralink_progress_update(ctx, 90, "Calculating metrics")
  metrics_report <- terralink_landscape_report(
    pre_mask = pre_mask,
    post_mask = terra::ifel(contiguous_raster > 0, 1, NA),
    raster = raster,
    units = units,
    label = terralink_safe_name("TerraLink Raster")
  )
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
    max_search_distance = max_search_distance,
    patch_connectivity = patch_connectivity,
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
    metrics_report = metrics_report
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
