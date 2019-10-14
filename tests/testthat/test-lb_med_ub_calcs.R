context('lb / med / ub computation')

test_that('check_lb_med_ub works using known distributions', {
            if (FALSE){
              devtools::load_all()
            }

  res <- check_lb_med_ub(lb = -1.96, 
                         med = 0, 
                         ub = 1.96, 
                         fun = dnorm, 
                         range_start = -1000,
                         range_end = 1000)
  expect_lte( (res$area_left_of_lb - 0.025)^2, 0.00001)
  expect_lte( (res$area_left_of_med - 0.5)^2 , 0.00001)
  expect_lte( (res$area_left_of_ub - 0.975)^2, 0.00001)

  res <- check_lb_med_ub(lb  = qbeta(0.025, 2, 5), 
                         med = qbeta(0.5,   2, 5), 
                         ub  = qbeta(0.975, 2, 5), 
                         fun = function(x){dbeta(x, 2, 5)}, 
                         range_start = -100,
                         range_end = 100)
  expect_lte( (res$area_left_of_lb - 0.025)^2, 0.00001)
  expect_lte( (res$area_left_of_med - 0.5)^2 , 0.00001)
  expect_lte( (res$area_left_of_ub - 0.975)^2, 0.00001)
})

test_that('check_lb_med_ub works using a trivial distribution', {
  res <- check_lb_med_ub(lb = 2.5, 
                         med = 50, 
                         ub = 97.5, 
                         fun = function(x){ifelse(x<0, 0, ifelse(x>100, 0, 1/100))}, 
                         range_start = -10000,
                         range_end = 10000)
  expect_lte( (res$area_left_of_lb - 0.025)^2, 0.00001)
  expect_lte( (res$area_left_of_med - 0.5)^2 , 0.00001)
  expect_lte( (res$area_left_of_ub - 0.975)^2, 0.00001)
})

#estimate_lb_med_ub <- function(fun, range_start, range_end){
test_that('estimate_lb_med_ub works on basic functions', {
  if (FALSE){
  devtools::load_all()
  }
  fun <- function(x){ifelse(x<0, 0, ifelse(x>100, 0, 1/100))}
  range_to_int <- trim_range(fun = fun, range_start = -1000, range_end = 1000)
  res <- estimate_lb_med_ub(fun = fun,
                            range_start = range_to_int$range_start,
                            range_end = range_to_int$range_end)
  expect_lte((res$lb  - 2.5 )^2, 0.0001)
  expect_lte((res$med - 50  )^2, 0.0001)
  expect_lte((res$ub  - 97.5)^2, 0.0001)

  range_to_int <- trim_range(fun = dnorm, range_start = -10000, range_end = 10000)
  res <- estimate_lb_med_ub(fun = dnorm,
                            range_start = range_to_int$range_start,
                            range_end = range_to_int$range_end)
  expect_lte((res$lb  - qnorm(0.025))^2, 0.0001)
  expect_lte((res$med - qnorm(0.5)  )^2, 0.0001)
  expect_lte((res$ub  - qnorm(0.975))^2, 0.0001)

  range_to_int <- trim_range(fun = dexp, range_start = -10000, range_end = 10000)
  res <- estimate_lb_med_ub(fun = dexp,
                            range_start = range_to_int$range_start,
                            range_end = range_to_int$range_end)
  expect_lte((res$lb  - qexp(0.025))^2, 0.0001)
  expect_lte((res$med - qexp(0.5)  )^2, 0.0001)
  expect_lte((res$ub  - qexp(0.975))^2, 0.0001)
})

test_that('estimate_lb_med_ub works with diagnostic histories', {
  ihists <- read.csv('/fridge/data/tsic/test_data.csv', stringsAsFactors = F)
  ihist <- subset(ihists, ptid == 'p01')
  # weibull dynamics aptima + architect
  ihist <- structure(list(ptid = c("p01", "p01", "p01", "p01", "p01", "p01",
"p01", "p01", "p01", "p01", "p01", "p01"), sample_date = c("2017-02-01",
"2017-03-01", "2017-04-01", "2017-05-01", "2017-06-01", "2017-07-01",
"2017-08-01", "2017-08-01", "2017-09-01", "2017-09-01", "2017-10-01",
"2017-10-01"), test = c("aptima_weib3_delaney", "aptima_weib3_delaney",
"aptima_weib3_delaney", "aptima_weib3_delaney", "aptima_weib3_delaney",
"aptima_weib3_delaney", "aptima_weib3_delaney", "architect_weib3_delaney",
"aptima_weib3_delaney", "architect_weib3_delaney", "aptima_weib3_delaney",
"architect_weib3_delaney"), result = c("-", "-", "-", "-", "-",
"-", "+", "+", "+", "+", "+", "+")), row.names = c(NA, 12L), class = "data.frame")
  ihist$sample_date <- as.numeric(as.Date(ihist$sample_date))
  
  agg_interpreter <- construct_aggregate_interpreter(ihist)
  range_start <- min(ihist$sample_date)
  range_end <- max(ihist$sample_date)

  range_to_int <- trim_range(fun = agg_interpreter, range_start = range_start, range_end = range_end)
  res <- estimate_lb_med_ub(fun = agg_interpreter,
                            range_start = range_to_int$range_start,
                            range_end = range_to_int$range_end,
                            verbose = FALSE)
  expect_lte(res$lb, res$med)
  expect_lte(res$med, res$ub)
  expect_lte(min(ihist$sample_date) - 90, res$lb)
  expect_gte(max(ihist$sample_date) + 90, res$ub)
})




