test_that("raster impassable routing avoids obstacle cells", {
  skip_if_not_installed("gdistance")
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")

  r <- terra::rast(nrows = 80, ncols = 80, xmin = 0, xmax = 80, ymin = 0, ymax = 80)
  r[,] <- 0
  r[20:30, 10:18] <- 1
  r[20:30, 62:70] <- 1
  r[1:80, 40:41] <- 2
  r[37:44, 40:41] <- 0

  result <- terralink_raster(
    raster = r,
    patch_values = 1,
    budget = 260,
    min_patch_size = 5,
    min_corridor_width = 3,
    max_search_distance = 120,
    obstacle_values = 2,
    allow_bottlenecks = FALSE,
    units = "pixels",
    verbose = 0,
    progress = FALSE
  )

  expect_gte(result$summary$corridors_used %||% 0, 1)
  corr_vals <- terra::values(result$corridor_raster, mat = FALSE)
  obs_vals <- terra::values(r, mat = FALSE)
  overlap <- sum((corr_vals > 0) & (obs_vals == 2), na.rm = TRUE)
  expect_equal(overlap, 0)
})

test_that("vector impassable routing avoids obstacle polygons", {
  skip_if_not_installed("gdistance")
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")
  skip_if_not_installed("sf")

  poly <- function(xmin, ymin, xmax, ymax) {
    sf::st_polygon(list(rbind(
      c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin)
    )))
  }

  patches <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      poly(60, 320, 220, 480),
      poly(780, 320, 940, 480)
    ),
    crs = 3857
  )

  obstacles <- sf::st_sf(
    oid = 1:2,
    geometry = sf::st_sfc(
      poly(430, 0, 570, 380),
      poly(430, 420, 570, 800)
    ),
    crs = 3857
  )

  result <- terralink_vector(
    patches = patches,
    budget = 8,
    min_patch_size = 0.01,
    min_corridor_width = 20,
    max_search_distance = 2000,
    obstacle_layers = obstacles,
    obstacle_resolution = 10,
    units = "metric",
    return_crs = "input",
    verbose = 0,
    progress = FALSE
  )

  expect_true(inherits(result$corridors, "sf"))
  expect_gte(nrow(result$corridors), 1)
  inter <- suppressWarnings(tryCatch(
    sf::st_intersection(sf::st_union(sf::st_geometry(result$corridors)), sf::st_union(sf::st_geometry(obstacles))),
    error = function(e) NULL
  ))
  overlap_area <- 0
  if (!is.null(inter) && length(inter) > 0 && !all(sf::st_is_empty(inter))) {
    overlap_area <- as.numeric(sf::st_area(inter))
  }
  expect_equal(overlap_area, 0, tolerance = 1e-6)
})
