test_that("terralink_vector returns corridors and networks", {
  skip_if_not_installed("sf")

  p1 <- sf::st_polygon(list(rbind(c(0, 0), c(0, 10), c(10, 10), c(10, 0), c(0, 0))))
  p2 <- sf::st_polygon(list(rbind(c(30, 0), c(30, 10), c(40, 10), c(40, 0), c(30, 0))))
  patches <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(p1, p2), crs = 32618)

  result <- terralink_vector(
    patches = patches,
    budget = 10,
    min_patch_size = 0.001,
    min_corridor_width = 5,
    max_search_distance = 100,
    units = "metric"
  )

  expect_true(inherits(result, "terralink_result"))
  expect_true(inherits(result$corridors, "sf"))
  expect_true(inherits(result$networks, "sf") || is.null(result$networks))
  expect_true(is.list(result$summary))
})

test_that("synthetic vector fixture stays close to current QGIS vector outputs", {
  skip_if_not_installed("sf")

  expect_rel_close <- function(actual, expected, rel_tol, abs_tol = 1e-6) {
    scale <- max(abs(expected), abs_tol)
    expect_lte(abs(actual - expected), max(scale * rel_tol, abs_tol))
  }

  vec <- terralink_sample_data("synthetic_vector")
  obs <- terralink_sample_data("synthetic_obstacle")
  expect_true(length(vec) == 1)
  expect_true(file.exists(vec))

  runs <- list(
    most_connected_habitat = terralink_vector(
      patches = vec,
      obstacle_layers = obs,
      budget = 18,
      min_patch_size = 0.01,
      min_corridor_width = 60,
      max_search_distance = 900,
      units = "metric",
      strategy = "most_connected_habitat"
    ),
    largest_single_network = terralink_vector(
      patches = vec,
      obstacle_layers = obs,
      budget = 18,
      min_patch_size = 0.01,
      min_corridor_width = 60,
      max_search_distance = 900,
      units = "metric",
      strategy = "largest_single_network"
    ),
    reachable_habitat_advanced = terralink_vector(
      patches = vec,
      obstacle_layers = obs,
      budget = 18,
      min_patch_size = 0.01,
      min_corridor_width = 60,
      max_search_distance = 900,
      units = "metric",
      strategy = "reachable_habitat_advanced",
      species_dispersal_distance = 800,
      patch_quality_field = "quality"
    ),
    landscape_fluidity = terralink_vector(
      patches = vec,
      obstacle_layers = obs,
      budget = 18,
      min_patch_size = 0.01,
      min_corridor_width = 60,
      max_search_distance = 900,
      units = "metric",
      strategy = "landscape_fluidity"
    )
  )

  expected <- list(
    most_connected_habitat = list(
      pairs = c("1-2", "3-4"),
      corridors_used = 2L,
      budget_used = 7.049097299531684,
      total_connected = 22.973001481098706,
      largest_network = 7.9619715847678485,
      mean_effective_resistance_pre = 32696.508479327207,
      mean_effective_resistance_post = 21996.66719796463,
      landscape_fluidity_pre = 0.0004870215482376388,
      landscape_fluidity_post = 0.000723923494330108,
      tol = list(budget = 0.011, total_connected = 0.005, largest_network = 0.005, resistance = 0.01, fluidity = 0.01)
    ),
    largest_single_network = list(
      pairs = c("1-2", "1-3", "2-3", "3-4"),
      corridors_used = 4L,
      budget_used = 13.589674830344116,
      total_connected = 29.51357901191114,
      largest_network = 15.923904181567021,
      mean_effective_resistance_pre = 32696.508479327207,
      mean_effective_resistance_post = 652.6075621809816,
      landscape_fluidity_pre = 0.0004870215482376388,
      landscape_fluidity_post = 0.024562498121177504,
      tol = list(budget = 0.006, total_connected = 0.005, largest_network = 0.005, resistance = 0.01, fluidity = 0.01)
    ),
    reachable_habitat_advanced = list(
      pairs = c("1-2", "1-3", "2-3", "3-4"),
      corridors_used = 4L,
      budget_used = 13.589674830344117,
      total_connected = 29.51357901191114,
      largest_network = 15.923904181567021,
      habitat_availability_post = 19.681070141617113,
      mean_reachable_area_post = 2.4612737742998294,
      largest_reachable_habitat_cluster_post = 15.923904181567021,
      mean_effective_resistance_pre = 32696.508479327207,
      mean_effective_resistance_post = 652.6075621809816,
      landscape_fluidity_pre = 0.0004870215482376388,
      landscape_fluidity_post = 0.024513021296251578,
      tol = list(
        budget = 0.006,
        total_connected = 0.005,
        largest_network = 0.005,
        habitat = 0.01,
        mean_reachable = 0.01,
        largest_reachable = 0.005,
        resistance = 0.01,
        fluidity = 0.01
      )
    ),
    landscape_fluidity = list(
      pairs = c("1-2", "1-3", "2-3", "3-4"),
      corridors_used = 4L,
      budget_used = 13.58967483034412,
      total_connected = 29.51357901191114,
      largest_network = 15.923904181567021,
      mean_effective_resistance_pre = 32696.508479327207,
      mean_effective_resistance_post = 652.6075621809816,
      landscape_fluidity_pre = 0.0004870215482376388,
      landscape_fluidity_post = 0.024550551806054553,
      tol = list(budget = 0.006, total_connected = 0.005, largest_network = 0.005, resistance = 0.01, fluidity = 0.01)
    )
  )

  for (name in names(runs)) {
    result <- runs[[name]]
    pair_sig <- unname(sort(apply(as.data.frame(result$corridors)[, c("patch1", "patch2"), drop = FALSE], 1, function(x) paste(x, collapse = "-"))))
    expect_equal(result$summary$strategy, name)
    expect_equal(result$summary$corridors_used, expected[[name]]$corridors_used)
    expect_equal(pair_sig, expected[[name]]$pairs)
    expect_rel_close(result$summary$budget_used, expected[[name]]$budget_used, expected[[name]]$tol$budget)
    expect_rel_close(result$metrics$total_connected_habitat_area_post, expected[[name]]$total_connected, expected[[name]]$tol$total_connected)
    expect_rel_close(result$metrics$largest_network_area_post, expected[[name]]$largest_network, expected[[name]]$tol$largest_network)
    expect_rel_close(result$metrics$mean_effective_resistance_pre, expected[[name]]$mean_effective_resistance_pre, expected[[name]]$tol$resistance)
    expect_rel_close(result$metrics$mean_effective_resistance_post, expected[[name]]$mean_effective_resistance_post, expected[[name]]$tol$resistance)
    expect_rel_close(result$metrics$landscape_fluidity_pre, expected[[name]]$landscape_fluidity_pre, expected[[name]]$tol$fluidity)
    expect_rel_close(result$metrics$landscape_fluidity_post, expected[[name]]$landscape_fluidity_post, expected[[name]]$tol$fluidity)
    if (!is.null(expected[[name]]$habitat_availability_post)) {
      expect_rel_close(result$metrics$habitat_availability_post, expected[[name]]$habitat_availability_post, expected[[name]]$tol$habitat)
    }
    if (!is.null(expected[[name]]$mean_reachable_area_post)) {
      expect_rel_close(result$metrics$mean_reachable_area_post, expected[[name]]$mean_reachable_area_post, expected[[name]]$tol$mean_reachable)
    }
    if (!is.null(expected[[name]]$largest_reachable_habitat_cluster_post)) {
      expect_rel_close(
        result$metrics$largest_reachable_habitat_cluster_post,
        expected[[name]]$largest_reachable_habitat_cluster_post,
        expected[[name]]$tol$largest_reachable
      )
    }
  }
})
