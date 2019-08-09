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

test_that("get_scatterpoints adds in enough intervals for a sloped line between two points", {
  x <- 0:1/40
  fun <- function(x){return(x)}
  scatterpoints <- get_scatterpoints(fun = fun, seedpoints = x, 
                                     max_delta = 0.02, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)
  delta <- abs(scatterpoints$y[1:(length(scatterpoints$y)-1)] - scatterpoints$y[2:length(scatterpoints$y)])
  int_lengths <- abs(scatterpoints$x[1:(length(scatterpoints$x)-1)] - scatterpoints$x[2:length(scatterpoints$x)])
  expect_true(all(sort(names(scatterpoints)) == c('x', 'y')))
  expect_equal(length(scatterpoints$x), length(unique(scatterpoints$x)))
  expect_equal(length(scatterpoints$x), length(scatterpoints$y))
  expect_true(all(delta < 0.02)) 
  expect_true(all(int_lengths > 0.01/20)) 
  expect_equal(length(scatterpoints$x), 20)
})

test_that("get_scatterpoints adds in enough intervals for a sloped line between many points", {
  x <- 0:10/40
  fun <- function(x){return(x)}
  scatterpoints <- get_scatterpoints(fun = fun, seedpoints = x, 
                                     max_delta = 0.02, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)
  delta <- abs(scatterpoints$y[1:(length(scatterpoints$y)-1)] - scatterpoints$y[2:length(scatterpoints$y)])
  int_lengths <- abs(scatterpoints$x[1:(length(scatterpoints$x)-1)] - scatterpoints$x[2:length(scatterpoints$x)])
  expect_true(all(sort(names(scatterpoints)) == c('x', 'y')))
  expect_equal(length(scatterpoints$x), length(unique(scatterpoints$x)))
  expect_equal(length(scatterpoints$x), length(scatterpoints$y))
  expect_true(all(delta < 0.02)) 
  expect_true(all(int_lengths > 0.01/20)) 
  expect_equal(length(scatterpoints$x), 191)
})

test_that("min_length prevents get_scatterpoints from adding in too many points", {
  x <- 0:1/40
  fun <- function(x){return(x)}
  scatterpoints <- get_scatterpoints(fun = fun, seedpoints = x, 
                                     max_delta = 0.0001, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)
  delta <- abs(scatterpoints$y[1:(length(scatterpoints$y)-1)] - scatterpoints$y[2:length(scatterpoints$y)])
  int_lengths <- abs(scatterpoints$x[1:(length(scatterpoints$x)-1)] - scatterpoints$x[2:length(scatterpoints$x)])
  expect_true(all(sort(names(scatterpoints)) == c('x', 'y')))
  expect_equal(length(scatterpoints$x), length(unique(scatterpoints$x)))
  expect_equal(length(scatterpoints$x), length(scatterpoints$y))
  expect_true(all(int_lengths > 0.01/20)) 
  expect_equal(length(scatterpoints$x), 20)
})

test_that("get_scatterpoints can recurse more than once", {
  x <- 0:1/40
  fun <- function(x){return(x)}
  scatterpoints <- get_scatterpoints(fun = fun, seedpoints = x, 
                                     max_delta = 0.0001, min_length = 0.001,
                                     n_new_segments = 20, verbose = FALSE)
  delta <- abs(scatterpoints$y[1:(length(scatterpoints$y)-1)] - scatterpoints$y[2:length(scatterpoints$y)])
  int_lengths <- abs(scatterpoints$x[1:(length(scatterpoints$x)-1)] - scatterpoints$x[2:length(scatterpoints$x)])
  expect_true(all(sort(names(scatterpoints)) == c('x', 'y')))
  expect_equal(length(scatterpoints$x), length(unique(scatterpoints$x)))
  expect_equal(length(scatterpoints$x), length(scatterpoints$y))
  expect_true(all(delta < 0.0001)) 
  expect_true(all(int_lengths > 0.001/20)) 
  expect_equal(length(scatterpoints$x), 362)
})

test_that("get_scatterpoints zoom in on the step", {
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  x <- 1:20
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  expect_equal(sum(scatterpoints$x <= 9), length(1:9))
  expect_equal(sum(scatterpoints$x >= 11), length(11:20))
  expect_true(sum(scatterpoints$x>=9 & scatterpoints$x <= 11) > 19)
  
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  x <- 1:20
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  expect_equal(sum(scatterpoints$x <= 9), length(1:9))
  expect_equal(sum(scatterpoints$x >= 11), length(11:20))
  expect_true(sum(scatterpoints$x>=9 & scatterpoints$x <= 11) > 19)
})

test_that('reduce_x_points works in a trivial case', {
  x <- 1:10
  y <- rep(0, 10)
  min_delta <- 0.0001
  result <- reduce_x_points(x = x, y = y, min_delta = min_delta)
  expect_equal(result$x, c(1, 10))
  expect_equal(result$y, c(0, 0))

  x_matches <- match(result$x, x)
  expect_true(all(result$y == y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = x)
  expect_true(all(abs(interpolated$y - y) < min_delta))
})

test_that('reduce_x_points handles changes in the y values correctly', {
  x <- 1:10
  y <- c(rep(0, 8), 1, 0)
  min_delta <- 0.0001
  result <- reduce_x_points(x = x, y = y, min_delta = min_delta)

  x_matches <- match(result$x, x)
  expect_true(all(result$y == y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = x)
  expect_true(all(abs(interpolated$y - y) < min_delta))
})

test_that('reduce_x_points handles a step function correctly', {
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  x <- 1:20
  y <- NULL
  for (i in x){y <- c(y, foo(i))}
  min_delta <- 0.0001
  result <- reduce_x_points(x = x, y = y, min_delta = min_delta)

  x_matches <- match(result$x, x)
  expect_true(all(result$y == y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = x)
  expect_true(all(abs(interpolated$y - y) < min_delta))
})

test_that('INTEGRATION: get_scatterpoints and reduce_x_points on Linear', {
  assay_dynamics <- get_assay_dynamics(assay = 'linear_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  x <- 1:20
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  min_delta <- 0.0001
  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = min_delta)

  x_matches <- match(result$x, scatterpoints$x)
  expect_true(all(result$y == scatterpoints$y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = scatterpoints$x)
  expect_true(all(abs(interpolated$y - scatterpoints$y) < min_delta))
})

test_that('INTEGRATION: get_scatterpoints and reduce_x_points on Weib3', {
  assay_dynamics <- get_assay_dynamics(assay = 'weib3_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  x <- 1:50
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  min_delta <- 0.0001
  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = min_delta)

  x_matches <- match(result$x, scatterpoints$x)
  expect_true(all(result$y == scatterpoints$y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = scatterpoints$x)
  expect_true(all(abs(interpolated$y - scatterpoints$y) < min_delta))
})



