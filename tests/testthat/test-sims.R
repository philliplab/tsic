context("test-sims")

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

test_that("sim_sc_times is sane", {
  list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                      "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                      "geenius_fr_weib3_delaney")

  expect_error(sim_sc_times(c(list_of_assays[5], list_of_assays[1:4]), skip_order_check = FALSE))

  # stochasticity can cause issues in larger runs - do a whole bunch of checks
  for (i in 1:100){
    x <- sim_sc_times(list_of_assays)
    expect_equal(class(x), 'list')
    expect_true(all(names(x) %in% list_of_assays))
    expect_true(all(unlist(x) > 0 & unlist(x) < 150))
    for (j in 1:(length(x)-1)){
      expect_lte(x[[j]], x[[j+1]])
    }
  }

})

test_that('fix_draw parameter of sim_sc_times works', {
  list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                      "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                      "geenius_fr_weib3_delaney")
  x <- sim_sc_times(list_of_assays, fix_draw = 0.5)
  expect_lte(max(x$iscav2_weib3_delaney_and_tosiano/8.760281, 8.760281/x$iscav2_weib3_delaney_and_tosiano),
             1.01)
  expect_lte(max(x$taqman_weib3_delaney_and_manufacturer/12.1334, 12.1334/x$taqman_weib3_delaney_and_manufacturer),
             1.01)
  expect_lte(max(x$architect_weib3_delaney/17.87017, 17.87017/x$architect_weib3_delaney),
             1.01)
  expect_lte(max(x$geenius_indet_weib3_delaney/27.0469, 27.0469/x$geenius_indet_weib3_delaney),
             1.01)
  expect_lte(max(x$geenius_fr_weib3_delaney/32.87321, 32.87321/x$geenius_fr_weib3_delaney),
             1.01)
})

test_that('sim_sc_times, combine_sc_and_visit_times and select_most_informative_results works together', {
  list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                      "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                      "geenius_fr_weib3_delaney")
  sc_times <- sim_sc_times(list_of_assays, fix_draw = 0.5)
  visit_times <- ((0:2)*28) - 14
  ihist <- combine_sc_and_visit_times(sc_times, visit_times)
  mihist <- select_most_informative_results(ihist)

  expected_mihist <- structure(list(
    ptid = c("?", "?", "?", "?"), 
    sample_date = c(16986.5, 17014.5, 17014.5, 17042.5), 
    test = c("iscav2_weib3_delaney_and_tosiano", "architect_weib3_delaney", 
             "taqman_weib3_delaney_and_manufacturer", "geenius_fr_weib3_delaney"), 
    result = c("-", "-", "+", "+")), 
    row.names = c( 1L, 8L, 7L, 15L), 
    class = "data.frame")
  expect_equal(mihist$kept_ihist, expected_mihist)
})




