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

test_that('sim_visit_dates is sane', {
  durations <- 
    sim_durations(n = 30,
                  rates = 'fiebig',
                  eclipse_distribution <- list(
                    distr = function(n, shape, scale, location) {
                      location + rweibull(n, shape = shape, scale = scale)
                    },
                    params = list(shape = 1.35, scale = 9, location = 4.8)))
  visit_dates <-
    sim_visit_dates(durations = durations, 
                    inf_date = as.numeric(as.Date('2015-05-01')),
                    n_before_infection = 2,
                    gap_distribution = function(n){rnorm(n, 30, 3)})
  expect_equal(class(visit_dates), 'data.frame')
  expect_equal(names(visit_dates), c('ptid', 'visit_id', 'visit_date'))
  expect_equal(class(visit_dates$visit_date), 'numeric')
  expect_true( class(visit_dates$visit_id) %in% c('numeric', 'integer'))

  for (c_ptid in sort(unique(durations$ptid))){
    total_duration <- sum(durations[durations$ptid == c_ptid, -1])
    c_visit_dates <- subset(visit_dates, ptid == c_ptid)
    follow_up <- max(c_visit_dates$visit_date) - min(c_visit_dates$visit_date)
    expect_gte(follow_up, total_duration)
    expect_equal(c_visit_dates$visit_id, 1:nrow(c_visit_dates))
    for (i in 1:(nrow(c_visit_dates) - 1)){
      expect_lte(c_visit_dates$visit_date[i], c_visit_dates$visit_date[i+1])
    }
  }
})






