test_that("circuit utility selector matches plugin-style greedy picks", {
  candidates <- data.frame(
    id = 1:6,
    patch1 = c(1, 2, 3, 1, 1, 2),
    patch2 = c(2, 3, 4, 4, 3, 4),
    cost = c(3, 2, 2, 5, 4, 4.5),
    length = c(3, 2, 2, 5, 4, 4.5),
    roi = c(
      sqrt(10 * 8) / 3,
      sqrt(8 * 6) / 2,
      sqrt(6 * 4) / 2,
      sqrt(10 * 4) / 5,
      sqrt(10 * 6) / 4,
      sqrt(8 * 4) / 4.5
    )
  )
  nodes <- c("1" = 10, "2" = 8, "3" = 6, "4" = 4)

  res <- select_circuit_utility(
    candidates = candidates,
    budget = 9,
    get_patch_ids = function(cand) c(cand$patch1, cand$patch2),
    get_pair_key = function(cand) sort(c(cand$patch1, cand$patch2)),
    get_cost = function(cand) cand$cost,
    get_base_roi = function(cand) cand$roi,
    get_length = function(cand) cand$length,
    get_patch_size = function(pid) nodes[[as.character(pid)]],
    overlap_ratio = function(cand, prior) 0,
    overlap_obj = function(cand) NULL
  )

  expect_equal(res$selected_ids, c(1L, 2L, 3L))
  expect_equal(res$stats$budget_used, 7)
  expect_equal(res$stats$budget_remaining, 2)
  expect_equal(res$stats$corridors_used, 3L)
})
