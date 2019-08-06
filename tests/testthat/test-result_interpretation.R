context("test-result_interpretation")

test_that("assay_result_interpreter gets build correctly", {
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  expect_type(foo, 'closure')
  expect_equal(foo(20), 1)
  expect_equal(foo(2), 0)
  
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  expect_type(foo, 'closure')
  expect_equal(foo(20), 0)
  expect_equal(foo(2), 1)
})


test_that("mucking about in the outer scope does not mess with the closure", {
  # These tests are only necessary because I am not that familiar with closures
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  expect_type(foo, 'closure')
  expect_equal(foo(20), 1)
  expect_equal(foo(2), 0)
  
  assay_dynamics$parameters$diagnostic_delay <- 100
  expect_equal(foo(20), 1)
  assay_dynamics$parameters$diagnostic_delay <- 1.5
  expect_equal(foo(2), 0)
  
  result <- '-'
  expect_equal(foo(20), 1)
  expect_equal(foo(2), 0)
  
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  expect_type(foo, 'closure')
  expect_equal(foo(20), 0)
  expect_equal(foo(2), 1)
  
  assay_dynamics$parameters$diagnostic_delay <- 100
  expect_equal(foo(20), 0)
  assay_dynamics$parameters$diagnostic_delay <- 1.5
  expect_equal(foo(2), 1)
  
  result <- '+'
  expect_equal(foo(20), 0)
  expect_equal(foo(2), 1)
})

test_that("get_scatterpoints work in normal circumstances", {
  x <- 0:10/1000
  fun <- function(x){return(x)}
  scatterpoints <- get_scatterpoints(fun = fun, seedpoints = x)
  expect_true(all(scatterpoints$x == x))
  expect_true(all(scatterpoints$y == x))
})
