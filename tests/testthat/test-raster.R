test_that("terralink_raster returns corridor output", {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
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
