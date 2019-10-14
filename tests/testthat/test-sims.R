context("test-sims")

test_that("Test sim_durations is sane", {
  if (FALSE){
    devtools::load_all()
  }
  durations <- sim_durations(n = 30,
                             rates = 'fiebig',
                             eclipse_distribution <- list(
                               distr = function(n, shape, scale, location) {
                                 location + rweibull(n, shape = shape, scale = scale)
                               },
                               params = list(shape = 1.35, scale = 9, location = 4.8)))
  expect_equal(class(durations), 'data.frame')
  expect_equal(names(durations), c('ptid', 'eclipse', 'l1', 'l2', 'l3', 'l4', 'l5'))
  expect_equal(nrow(durations), 30)
  for (c_col in c('eclipse', 'l1', 'l2', 'l3', 'l4', 'l5')){
    expect_equal(class(durations[, c_col]), 'numeric')
    expect_true(all(durations[, c_col] > 0))
  }
})

