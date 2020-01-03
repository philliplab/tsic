context("test-result_interpretation")

if (FALSE) {
  devtools::load_all()
}

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
            if (FALSE) {
              devtools::load_all()
            }
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
#  data.frame(x = scatterpoints$x, y=scatterpoints$y)

  result <- reduce_x_points(x = scatterpoints$x, y = scatterpoints$y, min_delta = 0.0001)
  result <- data.frame(x = result$x, y = result$y)

  if (FALSE){
    devtools::load_all()
    reduce_x_points <- reduce_x_points_new
  }

  for (i in 1:(nrow(result)-1)){
    if (any(result$y[i+0:1] %in% c(0,1)) & any(!(result$y[i+0:1] %in% c(0,1)))){
      expect_lte(abs(result$x[i] - result$x[i+1]), 1)
    }
  }

#TODO: fix these
  expect_false(any(duplicated(result)))
#  # duplication issues only when the assay do not start at zero at range start
#  
#  # Just such horrible hacky algoritms, just make them pretty.
#  expect_equal(1, "elegant", label = "get_scatterpoints and reduce_x_points are")
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

test_that('which_is_faster works', {
  check_which_is_faster <- function(faster, slower){
    if (FALSE){
      faster <- "iscav2_weib3_delaney_and_tosiano"
      slower <- "geenius_fr_weib3_delaney"
    }
    devtools::load_all()
    faster <- get_assay_dynamics(faster)
    slower <- get_assay_dynamics(slower)
    resultfs <- which_is_faster(assay1 = faster, assay2 = slower)
    resultsf <- which_is_faster(assay1 = slower, assay2 = faster)
    expect_equal(resultfs$diff, -resultsf$diff)
    expect_equal(resultfs$faster, faster)
    expect_equal(resultfs$slower, slower)
    expect_equal(resultsf$faster, faster)
    expect_equal(resultsf$slower, slower)
  }
  check_which_is_faster(faster = "iscav2_weib3_delaney_and_tosiano", 
                        slower = "geenius_fr_weib3_delaney")
  check_which_is_faster(faster = "architect_weib3_delaney", 
                        slower = "geenius_fr_weib3_delaney")
  check_which_is_faster(faster = "gs_combo_weib3_delaney", 
                        slower = "geenius_fr_weib3_delaney")
  check_which_is_faster(faster = "aptima_weib3_delaney", 
                        slower = "gs_combo_weib3_delaney")
  check_which_is_faster(faster = "aptima_weib3_delaney", 
                        slower = "architect_weib3_delaney")

  check_which_is_faster(faster = "iscav2_weib3_delaney_and_tosiano", 
                        slower = "abbott_real_time_weib3_delaney_and_manufacturer")
  check_which_is_faster(faster = "iscav2_weib3_delaney_and_tosiano", 
                        slower = "taqman_weib3_delaney_and_manufacturer")
  check_which_is_faster(faster = "taqman_weib3_delaney_and_manufacturer", 
                        slower = "gs_combo_weib3_delaney")
  check_which_is_faster(faster = "abbott_real_time_weib3_delaney_and_manufacturer", 
                        slower = "architect_weib3_delaney")
  check_which_is_faster(faster = "iscav2_weib3_delaney_and_tosiano", 
                        slower = "abbott_real_time_weib3_delaney_and_manufacturer")
  check_which_is_faster(faster = "iscav2_weib3_delaney_and_tosiano", 
                        slower = "taqman_weib3_delaney_and_manufacturer")
})

if (FALSE) {
  devtools::load_all()
}

test_that('select_most_informative_results works', {
  ihist <- data.frame(
    ptid = c('p0', 'p1'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-09-01')) + 0.5),
    test = c('step_unit_testing', 'step_unit_testing'),
    result = c('-', '+'),
    stringsAsFactors = FALSE
  )
  expect_error(select_most_informative_results(ihist, NULL))

  ihist <- data.frame(
    ptid = c('p0', 'p0'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-09-01')) + 0.5),
    test = c('step_unit_testing', 'step_unit_testing'),
    result = c('-', '+'),
    stringsAsFactors = FALSE
  )
  expect_error(select_most_informative_results(ihist, NULL))

  ihist <- data.frame(
    ptid = c("p314", "p314", "p314", "p314", "p314", "p314", "p314", "p314", "p314"), 
    sample_date = c(16860.5, 16910.5, 16921.5, 16921.5, 16910.5, 16921.5, 16921.5, 16910.5, 16860.5), 
    test = c("architect_weib3_delaney", "architect_weib3_delaney", "architect_weib3_delaney", 
             "geenius_fr_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_indet_weib3_delaney", 
             "taqman_weib3_delaney_and_manufacturer", "taqman_weib3_delaney_and_manufacturer", 
             "taqman_weib3_delaney_and_manufacturer"), 
    result = c("-", "+", "+", "+", "-", "+", "+", "+", "-"),
    stringsAsFactors = FALSE)
  expected_result <- list(kept_ihist = structure(list(ptid = c("p314", "p314", "p314",
"p314"), sample_date = c(16860.5, 16910.5, 16910.5, 16921.5),
    test = c("taqman_weib3_delaney_and_manufacturer", "geenius_indet_weib3_delaney",
    "architect_weib3_delaney", "geenius_fr_weib3_delaney"), result = c("-",
    "-", "+", "+")), row.names = c(9L, 5L, 2L, 4L), class = "data.frame"),
    rm_ihist = structure(list(ptid = c("p314", "p314", "p314",
    "p314", "p314"), sample_date = c(16860.5, 16910.5, 16921.5,
    16921.5, 16921.5), test = c("architect_weib3_delaney", "taqman_weib3_delaney_and_manufacturer",
    "architect_weib3_delaney", "geenius_indet_weib3_delaney",
    "taqman_weib3_delaney_and_manufacturer"), result = c("-",
    "+", "+", "+", "+")), row.names = c(1L, 8L, 3L, 6L, 7L), class = "data.frame"))

  inf_ihist <- select_most_informative_results(ihist)
  expect_true(class(inf_ihist) == 'list')
  expect_true(all(sort(names(inf_ihist)) == c('kept_ihist', 'rm_ihist')))
  expect_true(all(sort(unique(ihist$sample_date)) == sort(unique(inf_ihist$kept_ihist$sample_date))))
  expect_equal(nrow(ihist), nrow(inf_ihist$kept_ihist) + nrow(inf_ihist$rm_ihist))
  counts <- with(inf_ihist$kept_ihist, tapply(test, list(sample_date, result), length))
  expect_lte(max(counts, na.rm = TRUE), 1)
  expect_gte(min(counts, na.rm = TRUE), 0)
  expect_equal(inf_ihist, expected_result)
})







