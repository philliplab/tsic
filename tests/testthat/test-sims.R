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

test_that("the numbers from sim_durations are similar to the distributions they were drawn from", {
  durations <- sim_durations(n = 1000,
                             rates = 'fiebig',
                             eclipse_distribution <- list(
                               distr = function(n, shape, scale, location) {
                                 location + rweibull(n, shape = shape, scale = scale)
                               },
                               params = list(shape = 1.35, scale = 9, location = 4.8)))
  expect_gt(mean(durations$eclipse), 0.9*13)
  expect_lt(mean(durations$eclipse), 1.1*13)

  expect_gt(mean(durations$l1), 0.9*5.0)
  expect_lt(mean(durations$l1), 1.1*5.0)
  
  expect_gt(mean(durations$l2), 0.9*5.3)
  expect_lt(mean(durations$l2), 1.1*5.3)
  
  expect_gt(mean(durations$l3), 0.9*3.2)
  expect_lt(mean(durations$l3), 1.1*3.2)
  
  expect_gt(mean(durations$l4), 0.9*5.6)
  expect_lt(mean(durations$l4), 1.1*5.6)
  
  expect_gt(mean(durations$l5), 0.9*69.5)
  expect_lt(mean(durations$l5), 1.1*69.5)
})
