test_that("UnionFind tracks components and sizes", {
  uf <- UnionFind$new()
  uf$find(1)
  uf$find(2)
  uf$size[["1"]] <- 5
  uf$size[["2"]] <- 3
  uf$count[["1"]] <- 1L
  uf$count[["2"]] <- 1L

  expect_false(identical(uf$find(1), uf$find(2)))
  uf$union(1, 2)
  expect_true(identical(uf$find(1), uf$find(2)))
  expect_equal(uf$get_size(1), 8)
  expect_equal(uf$get_count(1), 2L)
})

test_that("NetworkOptimizer selects MST under budget", {
  nodes <- c("1" = 1, "2" = 1, "3" = 1, "4" = 1)
  edges <- data.frame(
    u = c(1, 2, 3, 1),
    v = c(2, 3, 4, 4),
    id = c(1, 2, 3, 4),
    cost = c(1, 1, 1, 10)
  )

  result <- optimize_network(nodes, edges, budget = 3)
  expect_equal(sort(result$selected), c(1L, 2L, 3L))
  expect_equal(result$total_cost, 3)
  expect_true(length(result$component_sizes) >= 1)
})

test_that("NetworkOptimizer respects loop_fraction and redundancy", {
  nodes <- c("1" = 1, "2" = 1, "3" = 1)
  edges <- data.frame(
    u = c(1, 2, 1),
    v = c(2, 3, 3),
    id = c(1, 2, 3),
    cost = c(1, 1, 1)
  )

  result <- optimize_network(nodes, edges, budget = 2, loop_fraction = 0.0)
  expect_equal(length(result$selected), 2)
})
