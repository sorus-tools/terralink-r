test_that("budget_pixels requires pixel units", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(r) <- 1

  expect_error(
    terralink_raster(
      raster = r,
      patch_values = 1,
      budget_pixels = 50,
      units = "metric",
      max_search_distance = 5,
      min_patch_size = 1
    ),
    class = "terralink_error_scale"
  )
})

test_that("candidate pair guard triggers for raster", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(r) <- 0
  terra::values(r)[1] <- 1
  terra::values(r)[50] <- 1
  terra::values(r)[100] <- 1

  expect_error(
    terralink_raster(
      raster = r,
      patch_values = 1,
      budget = 10,
      min_patch_size = 1,
      max_search_distance = 20,
      max_pair_checks = 1
    ),
    class = "terralink_error_compute_limit"
  )
})

test_that("invalid geometry is fixed or reported", {
  skip_if_not_installed("sf")

  bowtie <- sf::st_polygon(list(rbind(c(0, 0), c(1, 1), c(1, 0), c(0, 1), c(0, 0))))
  p2 <- sf::st_polygon(list(rbind(c(3, 0), c(3, 1), c(4, 1), c(4, 0), c(3, 0))))
  patches <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(bowtie, p2), crs = 32618)

  has_make_valid <- is.function(sf::st_make_valid) || (requireNamespace("lwgeom", quietly = TRUE) && "st_make_valid" %in% getNamespaceExports("lwgeom"))
  if (has_make_valid) {
    expect_error(
      terralink_vector(
        patches = patches,
        budget = 1,
        min_patch_size = 0,
        min_corridor_width = 5,
        max_search_distance = 100,
        units = "metric"
      ),
      NA
    )
  } else {
    expect_error(
      terralink_vector(
        patches = patches,
        budget = 1,
        min_patch_size = 0,
        min_corridor_width = 5,
        max_search_distance = 100,
        units = "metric"
      ),
      class = "terralink_error_geometry"
    )
  }
})

test_that("vector outputs return to input CRS by default", {
  skip_if_not_installed("sf")

  p1 <- sf::st_polygon(list(rbind(c(-73.0, 40.0), c(-73.0, 40.01), c(-72.99, 40.01), c(-72.99, 40.0), c(-73.0, 40.0))))
  p2 <- sf::st_polygon(list(rbind(c(-72.98, 40.0), c(-72.98, 40.01), c(-72.97, 40.01), c(-72.97, 40.0), c(-72.98, 40.0))))
  patches <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(p1, p2), crs = 4326)

  result <- terralink_vector(
    patches = patches,
    budget = 1,
    min_patch_size = 0.00001,
    min_corridor_width = 50,
    max_search_distance = 2000,
    units = "metric",
    return_crs = "input"
  )

  expect_true(inherits(result, "terralink_result"))
  expect_true(isTRUE(sf::st_crs(result$patches) == sf::st_crs(4326)))
})

test_that("metric raster units work with projected CRS", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:3857")
  v <- rep(0, terra::ncell(r))
  v[c(1:5, 96:100)] <- 1
  terra::values(r) <- v

  expect_error(
    terralink_raster(
      raster = r,
      patch_values = 1,
      budget = 0.01,
      min_patch_size = 0.0001,
      min_corridor_width = 1,
      max_search_distance = 30,
      units = "metric"
    ),
    NA
  )
})
