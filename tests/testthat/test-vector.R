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
