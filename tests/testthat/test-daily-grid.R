context('daily grid')

if(FALSE){
  devtools::load_all()
}

test_that('Daily Grid works on known distributions', {
  agg_fun <- function(x){dnorm(x, 17500, 10)}
  x <- compute_daily_grid(agg_fun, 1, 17400, 17600)
  expect_equal(sum(x$daily_probs), 1)
  expect_equal(17400, min(x$starts_of_daily_intervals))
  expect_equal(17599, max(x$starts_of_daily_intervals))

  agg_fun <- function(x){dweibull(x, 10, 10)}
  x <- compute_daily_grid(agg_fun, 1, 0, 20)
  expect_equal(sum(x$daily_probs), 1)
  expect_equal(0, min(x$starts_of_daily_intervals))
  expect_equal(19, max(x$starts_of_daily_intervals))
})





