test_that("diff history", {
  expect_basic_history <- function(ranks, history){
    expect_equal(length(history), length(ranks))
    expect_true(all(0 <= history & history <= 1))
  }

  expect_converge_to_exact <- function(ranks, history) {
    expect_basic_history(ranks, history)
    #expect_true(all(diff(history) <= 0))
    expect_equal(history[length(ranks)], 0)
  }

  test_converge_to_exact <- function(ranks, max_rank, p) {
    history <- compute_diff_history_norm(ranks, max_rank, p)
    expect_converge_to_exact(ranks, history)
  }

  test_converge_to_exact(0:100, max_rank = 100, p = 1)
  test_converge_to_exact(0:66, max_rank = 66, p = 2)
  test_converge_to_exact(0:213, max_rank = 213, p = 3)

  for(i in 1:5) {
    test_converge_to_exact(sample(0:84), max_rank = 84, p = 1)
    test_converge_to_exact(sample(0:331), max_rank = 331, p = 2)
  }

  history_all_0 <- compute_diff_history_norm(rep(0,5), max_rank = 3, p = 1)
  expect_basic_history(rep(0,5), history_all_0)
  expect_equal(history_all_0, rep(mean(1 - (1:4) / 4), 5))

  history_all_max <- compute_diff_history_norm(rep(6,4), max_rank = 6, p = 1)
  expect_basic_history(rep(6,4), history_all_max)
  expect_equal(history_all_max, rep(mean(1 - (1:7) / 7), 4))

})
