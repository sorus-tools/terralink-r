test_that("Graph metrics compute without error", {
  patches <- data.frame(id = c(1, 2, 3))
  corridors <- data.frame(
    patch1 = c(1, 2),
    patch2 = c(2, 3),
    distance_m = c(10, 20)
  )

  g <- build_graph_from_corridors(patches, corridors)
  expect_true(igraph::vcount(g) == 3)
  expect_true(igraph::ecount(g) == 2)

  h_mov <- calculate_movement_entropy(g)
  expect_true(is.finite(h_mov))

  top_pen <- calculate_topology_penalty(g)
  expect_equal(top_pen, 0)

  dist_pen <- calculate_disturbance_penalty(g)
  expect_true(dist_pen >= 0)

  totals <- calculate_total_entropy(g)
  expect_true(is.list(totals))
  expect_true(!is.null(totals$H_total))

  rho2 <- calculate_two_edge_connectivity(g)
  expect_true(rho2 >= 0)

  fail_prob <- calculate_failure_probability(g, k_failures = 1, iterations = 10)
  expect_true(fail_prob >= 0)
})
