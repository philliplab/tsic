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
test_that('estimate_lb_med_ub works', {
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




