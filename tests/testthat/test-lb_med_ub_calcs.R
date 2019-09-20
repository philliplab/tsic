context('lb / med / ub computation')

test_that('check_lb_med_ub works using known distributions', {
            if (FALSE){
              devtools::load_all()
            }

  res <- check_lb_med_ub(lb = -1.96, 
                         med = 0, 
                         ub = 1.96, 
                         fun = dnorm, 
                         range_start = -1000)
  expect_lte( (res$area_left_of_lb$value - 0.025)^2, 0.00001)
  expect_lte( (res$area_left_of_med$value - 0.5)^2 , 0.00001)
  expect_lte( (res$area_left_of_ub$value - 0.975)^2, 0.00001)

  res <- check_lb_med_ub(lb  = qbeta(0.025, 2, 5), 
                         med = qbeta(0.5,   2, 5), 
                         ub  = qbeta(0.975, 2, 5), 
                         fun = function(x){dbeta(x, 2, 5)}, 
                         range_start = 0,
                         range_end = 1)
  expect_lte( (res$area_left_of_lb$value - 0.025)^2, 0.00001)
  expect_lte( (res$area_left_of_med$value - 0.5)^2 , 0.00001)
  expect_lte( (res$area_left_of_ub$value - 0.975)^2, 0.00001)


#check_lb_med_ub <- function(lb, med, ub, fun, range_start){
})
