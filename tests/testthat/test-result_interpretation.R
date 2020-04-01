context("test-result_interpretation")

test_that("assay_result_interpreter gets build correctly", {
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = as.numeric(as.Date('2018-06-01'))+0.5)
  expect_type(foo, 'closure')
  expect_equal(foo(as.numeric(as.Date('2010-05-01'))), 1)
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 1)
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 0)
  expect_equal(foo(as.numeric(as.Date('2018-06-30'))), 0)
  
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = as.numeric(as.Date('2018-06-01'))+0.5)
  expect_type(foo, 'closure')
  expect_equal(foo(as.numeric(as.Date('2010-05-01'))), 0)
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 0)
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 1)
  expect_equal(foo(as.numeric(as.Date('2018-06-30'))), 1)
})

test_that("mucking about in the outer scope does not mess with the closure", {
  # These tests are only necessary because I am not that familiar with closures
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = as.numeric(as.Date('2018-06-01'))+0.5)
  expect_type(foo, 'closure')
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 1)
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 0)
  
  assay_dynamics$parameters$diagnostic_delay <- 100
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 1)
  assay_dynamics$parameters$diagnostic_delay <- 1.5
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 0)
  
  result <- '-'
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 1)
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 0)
  
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = as.numeric(as.Date('2018-06-01'))+0.5)
  expect_type(foo, 'closure')
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 0)
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 1)
  
  assay_dynamics$parameters$diagnostic_delay <- 100
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 0)
  assay_dynamics$parameters$diagnostic_delay <- 1.5
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 1)
  
  result <- '+'
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 0)
  expect_equal(foo(as.numeric(as.Date('2018-05-30'))), 1)
})

test_that("aggregate_interpreter gets build correctly in a basic case", {
  diagnostic_history <- data.frame(
    ptid = c('p0', 'p0'),
    sample_date = c(as.numeric(as.Date('2016-03-01')), as.numeric(as.Date('2016-09-01'))),
    test = c('step_unit_testing', 'step_unit_testing'),
    result = c('-', '+'),
    stringsAsFactors = FALSE
  )
  diagnostic_history$sample_date <- diagnostic_history$sample_date + 0.5
  foo <- construct_aggregate_interpreter(diagnostic_history)
  expect_type(foo, 'closure')
  expect_equal(foo(as.numeric(as.Date('2010-05-01'))), 0)
  expect_equal(foo(as.numeric(as.Date('2018-05-01'))), 0)
  expect_equal(foo(as.numeric(as.Date('2016-05-29'))), 1)
  expect_equal(foo(as.numeric(as.Date('2016-03-11'))), 1)
  expect_equal(foo(as.numeric(as.Date('2017-03-11'))), 0)

#  Somehow this is now working. I think it is a lazy evaluation thing - I now reassign the important variables in the environment and this fixes the issue. Unsure if it is the actual assignment or just the forcing of the thing to be evaluated that does the trick. LEARN MORE.

  # Find the jumps and make sure that there are no values that are not zero or one
  expect_equal(foo(as.numeric(as.Date('2016-02-19'))), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20'))+0.5), 1)
  expect_equal(foo(as.numeric(as.Date('2016-02-20'))), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.4)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.49)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.499)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.4999)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.49999)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.499999)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.4999999)), 0)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.49999999)), 0)

  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.6)), 1)
  expect_equal(foo(as.numeric(as.Date('2016-02-20')+0.500001)), 1)
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
  sample_date <- as.numeric(as.Date('2015-01-01')) + 0.5
  range_start <- as.numeric(as.Date('2014-12-01'))
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:sample_date
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  if (FALSE){
    plot(scatterpoints$y ~ scatterpoints$x)
    scatterpoints <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y)
    plot(scatterpoints$y ~ scatterpoints$x)
  }

  date_less_diagnostic_delay <- sample_date - assay_dynamics$params$diagnostic_delay
  day_before <- date_less_diagnostic_delay - 1
  day_after <- date_less_diagnostic_delay + 1

  expect_equal( sum(scatterpoints$x <= day_before), 
                length(range_start:day_before) )
  expect_equal( sum(scatterpoints$x >= day_after), 
                length(day_after:(sample_date-0.5)) )
  expect_true( sum( scatterpoints$x >= day_before & 
                    scatterpoints$x <= day_after) > 19 )
  
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:sample_date
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  expect_equal( sum(scatterpoints$x <= day_before), 
                length(range_start:day_before) )
  expect_equal( sum(scatterpoints$x >= day_after), 
                length(day_after:(sample_date-0.5)) )
  expect_true( sum( scatterpoints$x >= day_before & 
                    scatterpoints$x <= day_after) > 19 )
})

test_that('get_scatterpoints + reduce_x_points do not produce a little "trail" for aptima', {
  # funny thing happening where the aptima biomarker is computed with a long trail of points going from 0.994 to 1 over 12 days
  # try to figure out the interaction between get_scatterpoints and reduce_x_points that is causing this to happen

  sample_date <- as.numeric(as.Date('2015-01-01')) + 0.5
  range_start <- as.numeric(as.Date('2014-12-01'))
  assay_dynamics <- get_assay_dynamics(assay = 'aptima_weib3_delaney')
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:(sample_date+20)
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.001, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = 0.0001)
  result <- data.frame(x = result$x, y = result$y)

  for (i in 1:(nrow(result)-1)){
    if (any(result$y[i+0:1] %in% c(0,1)) & any(!(result$y[i+0:1] %in% c(0,1)))){
      expect_lte(abs(result$x[i] - result$x[i+1]), 1)
    }
  }

  expect_false(any(duplicated(result)))
})

test_that('reduce_x_points works in a trivial case', {
  x <- 1:10
  y <- rep(0, 10)
  min_delta <- 0.0001
  result <- reduce_x_points(x = x, y = y, min_delta = min_delta)
  expect_false(any(duplicated(result)))
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
  expect_false(any(duplicated(result)))

  x_matches <- match(result$x, x)
  expect_true(all(result$y == y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = x)
  expect_true(all(abs(interpolated$y - y) < min_delta))
})

test_that('reduce_x_points handles a step function correctly', {
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  sample_date <- as.numeric(as.Date('2015-01-01')) + 0.5
  range_start <- as.numeric(as.Date('2014-12-01'))
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:sample_date
  y <- NULL
  for (i in x){y <- c(y, foo(i))}
  min_delta <- 0.0001
  result <- reduce_x_points(x = x, y = y, min_delta = min_delta)

  x_matches <- match(result$x, x)
  expect_true(all(result$y == y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = x)
  expect_true(all(abs(interpolated$y - y) < min_delta))
})

test_that('interpret_ihist works', {
  if (FALSE) {
    devtools::load_all()
  }
  ihist <- data.frame(
    ptid = c('p0', 'p0'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-09-01')) + 0.5),
    test = c('step_unit_testing', 'step_unit_testing'),
    result = c('-', '+'),
    stringsAsFactors = FALSE
  )

  iihist <- interpret_ihist(ihist = ihist,
                            range_start = as.numeric(as.Date('2016-01-01')),
                            range_end = as.numeric(as.Date('2016-11-30')))

  expect_true('data.frame' %in% class(iihist))
  expect_equal(sort(names(iihist)),
               sort(c("ptid", "sample_date", "test_details", "prob_val")))

  for (c_assay_name in unique(iihist$test_details)){
    c_dat <- subset(iihist, test_details == c_assay_name)
    expect_equal(length(c_dat$sample_date),
                 length(unique(c_dat$sample_date)))
  }

#interpret_ihist <- function(ihist, range_start, range_end, verbose = FALSE){
})

test_that('trim_range works', {
  if (FALSE){
    devtools::load_all()
  }
  fun <- function(x){ifelse(x<400, 0, ifelse(x>450, 0, 1))}
  res <- trim_range(fun = fun, range_start = 0, range_end = 2000)
  expect_lte(res$range_start, 400)
  expect_gte(res$range_start, 400-10)
  
  expect_gte(res$range_end, 450)
  expect_lte(res$range_end, 450+10)
})

test_that('INTEGRATION: get_scatterpoints and reduce_x_points on Linear', {
  assay_dynamics <- get_assay_dynamics(assay = 'linear_unit_testing')
  result <- '+'
  sample_date <- as.numeric(as.Date('2015-01-01')) + 0.5
  range_start <- as.numeric(as.Date('2014-12-01'))
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:sample_date
  max_delta <- 0.01
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = max_delta, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  min_delta <- 0.0001
  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = min_delta)

  x_matches <- match(result$x, scatterpoints$x)
  expect_true(all(result$y == scatterpoints$y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = scatterpoints$x)
  expect_true(all(abs(interpolated$y - scatterpoints$y) < max_delta))
})

test_that('INTEGRATION: get_scatterpoints and reduce_x_points on Weib3', {
  if (FALSE){
    devtools::load_all()
  }
  assay_dynamics <- get_assay_dynamics(assay = 'weib3_unit_testing')
  result <- '+'
  sample_date <- as.numeric(as.Date('2015-01-01')) + 0.5
  range_start <- as.numeric(as.Date('2014-11-01'))
  range_end <- as.numeric(as.Date('2015-02-01'))
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:range_end
  max_delta <- 0.01
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = max_delta, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  min_delta <- 0.0001
  max_length <- 14
  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = min_delta, max_length = max_length)
  int_lengths <- abs(result$x[1:(length(result$x)-1)] - result$x[2:(length(result$x))])
  expect_equal(sum(result$x <= range_start), 1)
  expect_equal(sum(result$x >= range_end), 1)
  expect_false(any(int_lengths > max_length))

  x_matches <- match(result$x, scatterpoints$x)
  expect_true(all(result$y == scatterpoints$y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = scatterpoints$x)
  expect_true(all(abs(interpolated$y - scatterpoints$y) < max_delta))
})

test_that('PERFORMANCE: get_scatterpoints and reduce_x_points on Weib3 large interval', {
  assay_dynamics <- get_assay_dynamics(assay = 'weib3_unit_testing')
  result <- '+'
  sample_date <- as.numeric(as.Date('2015-01-01')) + 0.5
  range_start <- as.numeric(as.Date('2014-03-11'))
  range_end <- as.numeric(as.Date('2016-02-01'))
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result,
                                            sample_date = sample_date)
  x <- range_start:range_end
  scatterpoints <- get_scatterpoints(fun = foo, seedpoints = x, 
                                     max_delta = 0.01, min_length = 0.01,
                                     n_new_segments = 20, verbose = FALSE)

  min_delta <- 0.0001
  max_length <- 14
  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = min_delta, max_length = max_length)
  int_lengths <- abs(result$x[1:(length(result$x)-1)] - result$x[2:(length(result$x))])
  expect_equal(sum(result$x <= range_start), 1)
  expect_equal(sum(result$x >= range_end), 1)
  expect_false(any(int_lengths > (max_length+0.5)))

  x_matches <- match(result$x, scatterpoints$x)
  expect_true(all(result$y == scatterpoints$y[x_matches]))

  interpolated <- approx(x = result$x, y = result$y, xout = scatterpoints$x)
  expect_true(all(abs(interpolated$y - scatterpoints$y) < min_delta))
})

