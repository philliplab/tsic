context("test-sims")

if (FALSE){
  library(testthat)
}

test_that("sim_dx_results is sane", {
  list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                      "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                      "geenius_fr_weib3_delaney")

  expect_error(sim_dx_results(10, c(list_of_assays[5], list_of_assays[1:4]), skip_order_check = FALSE))

  x <- sim_dx_results(10, list_of_assays)

  expect_equal(class(x), 'list')
  expect_true(all(names(x) %in% list_of_assays))
  expect_true(all(unlist(x) %in% c('-','+')))

  x <- sim_dx_results(0, list_of_assays)
  expect_true(all(unlist(x) %in% c('-')))

  x <- sim_dx_results(10000, list_of_assays)
  expect_true(all(unlist(x) %in% c('+')))
})






