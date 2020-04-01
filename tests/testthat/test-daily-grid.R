context('daily grid')

test_that('Daily Grid works on known distributions', {
  agg_fun <- function(x){dnorm(x, 17500, 10)}
  x <- compute_daily_grid(agg_fun, 1, 17400, 17600)
  expect_equal(sum(x$mass), 1)
  expect_equal(17400, min(x$interval_start))
  expect_equal(17600, max(x$interval_end))

  agg_fun <- function(x){dweibull(x, 10, 10)}
  x <- compute_daily_grid(agg_fun, 1, 0, 20)
  expect_equal(sum(x$mass), 1)
  expect_equal(0, min(x$interval_start))
  expect_equal(20, max(x$interval_end))
})


test_that('Daily Grid works on known distributions', {
  ihist <- structure(list(ptid = c("p01", "p01", "p01", "p01", "p01", "p01", "p01", "p01", "p01", "p01", "p01", "p01"), 
                          sample_date = c("2017-06-01", "2017-07-01", "2017-08-01", "2017-08-01", "2017-09-01", "2017-09-01"), 
                          test = c("aptima_weib3_delaney", "aptima_weib3_delaney", "aptima_weib3_delaney", "architect_weib3_delaney", 
                                   "aptima_weib3_delaney", "architect_weib3_delaney"), 
                          result = c("-", "-", "+", "+", "+", "+")), 
                     row.names = c(NA, 12L), 
                     class = "data.frame")
  ihist$sample_date <- as.numeric(as.Date(ihist$sample_date)) + 0.5

  inf_ihist <- select_most_informative_results(ihist)$kept_ihist
  agg_inter <-  construct_aggregate_interpreter(inf_ihist)
  agg_fun <- agg_inter
  range_start <- floor(min(inf_ihist$sample_date) - 60)
  range_end <- ceiling(max(inf_ihist$sample_date) + 30)
  lb_med_ub <- estimate_lb_med_ub(agg_fun, range_start, range_end)
  tauc <- lb_med_ub$aoc

  x <- compute_daily_grid(agg_fun, tauc, range_start, range_end)
  expect_equal(sum(x$mass), 1)
  expect_equal(range_start, min(x$interval_start))
  expect_equal(range_end, max(x$interval_end))
  expect_true(all(x$mass >= 0))
  # agg function always lte 1
  expect_true(all(x$mass <= 1))
})

test_that('Daily grid expands automatically', {
  agg_fun <- function(x){dnorm(x, 17500, 10)}
  x <- compute_daily_grid(agg_fun, 1, 17500, 17600)
  expect_equal(sum(x$mass), 1)
  expect_gte(17480, min(x$interval_start))
  expect_equal(17600, max(x$interval_end))

  agg_fun <- function(x){dweibull(x, 10, 10)}
  x <- compute_daily_grid(agg_fun, 1, 0, 10)
  expect_equal(sum(x$mass), 1)
  expect_equal(0, min(x$interval_start))
  expect_lte(11, max(x$interval_end))
})


