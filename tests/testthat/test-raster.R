test_that("terralink_raster returns corridor output", {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000, crs = "EPSG:3857")
  terra::values(r) <- 0
  terra::values(r)[1:10] <- 1
  terra::values(r)[300:320] <- 1

  result <- terralink_raster(
    raster = r,
    patch_values = 1,
    budget = 200,
    min_patch_size = 3,
    max_search_distance = 50,
    units = "pixels"
  )

  expect_true(inherits(result, "terralink_result"))
  expect_true(inherits(result$patches, "SpatRaster"))
  expect_true(is.data.frame(result$corridors))
  expect_true(is.list(result$summary))
  expect_equal(result$summary$raster_delegate_backend, "vector_polygonize")
})

test_that("build_contiguous_raster includes weighted corridor cells", {
  habitat <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5)
  terra::values(habitat) <- NA
  habitat_vals <- terra::values(habitat, mat = FALSE)
  habitat_vals[c(1, 2, 24, 25)] <- 1
  terra::values(habitat) <- habitat_vals

  corridor <- habitat
  terra::values(corridor) <- NA
  corridor_vals <- terra::values(corridor, mat = FALSE)
  corridor_vals[13] <- 10
  terra::values(corridor) <- corridor_vals

  contiguous <- build_contiguous_raster(habitat, corridor, connectivity = 8)
  contiguous_vals <- terra::values(contiguous, mat = FALSE)

  expect_false(is.na(contiguous_vals[13]))
  expect_gt(contiguous_vals[13], 0)
})

test_that("synthetic raster fixture is stable across TerraLink 1.7 strategies", {
  synth <- terralink_sample_data("synthetic_raster")
  expect_true(length(synth) == 1)
  expect_true(file.exists(synth))

  runs <- suppressWarnings(list(
    most_connected_habitat = terralink_raster(
      raster = synth,
      patch_values = 1,
      budget = 120,
      min_patch_size = 4,
      min_corridor_width = 1,
      max_search_distance = 30,
      units = "pixels",
      strategy = "most_connected_habitat"
    ),
    largest_single_network = terralink_raster(
      raster = synth,
      patch_values = 1,
      budget = 120,
      min_patch_size = 4,
      min_corridor_width = 1,
      max_search_distance = 30,
      units = "pixels",
      strategy = "largest_single_network"
    ),
    reachable_habitat_advanced = terralink_raster(
      raster = synth,
      patch_values = 1,
      budget = 120,
      min_patch_size = 4,
      min_corridor_width = 1,
      max_search_distance = 30,
      units = "pixels",
      strategy = "reachable_habitat_advanced",
      species_dispersal_distance = 12
    ),
    landscape_fluidity = terralink_raster(
      raster = synth,
      patch_values = 1,
      budget = 120,
      min_patch_size = 4,
      min_corridor_width = 1,
      max_search_distance = 30,
      units = "pixels",
      strategy = "landscape_fluidity"
    )
  ))

  expected <- list(
    most_connected_habitat = list(
      corridors_used = 4,
      budget_used = 40.3603345958348,
      total_connected = 164.767630863503,
      largest_network = 74.6441467882577,
      fluidity_post = 3.7348044904674e-05
    ),
    largest_single_network = list(
      corridors_used = 4,
      budget_used = 27.7010838783668,
      total_connected = 152.108326495541,
      largest_network = 124.406924253536,
      fluidity_post = 0.00121354282468765
    ),
    reachable_habitat_advanced = list(
      corridors_used = 4,
      budget_used = 27.7010838783668,
      total_connected = 152.108326495541,
      largest_network = 124.406924253536,
      fluidity_post = 0.00121354282468765
    ),
    landscape_fluidity = list(
      corridors_used = 4,
      budget_used = 27.7010838783668,
      total_connected = 152.108326495541,
      largest_network = 124.406924253536,
      fluidity_post = 0.00121354282468765
    )
  )

  for (name in names(runs)) {
    result <- runs[[name]]
    expect_equal(result$summary$strategy, name)
    expect_equal(result$summary$corridors_used, expected[[name]]$corridors_used)
    expect_equal(result$summary$budget_used, expected[[name]]$budget_used, tolerance = 1e-6)
    expect_equal(result$summary$raster_delegate_backend, "vector_polygonize")
    expect_equal(result$metrics$total_connected_habitat_area_post, expected[[name]]$total_connected, tolerance = 1e-6)
    expect_equal(result$metrics$largest_network_area_post, expected[[name]]$largest_network, tolerance = 1e-6)
    expect_equal(as.numeric(result$metrics$landscape_fluidity_post), as.numeric(expected[[name]]$fluidity_post), tolerance = 1e-6)
  }
})
