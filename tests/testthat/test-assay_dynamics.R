context('assay dynamics')

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

test_that('check_assay_order works', {
  list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
                      "architect_weib3_delaney", "geenius_indet_weib3_delaney",
                      "geenius_fr_weib3_delaney")
  expect_true(check_assay_order(list_of_assays))
  expect_false(check_assay_order(c(list_of_assays[5], list_of_assays[1:4]), verbose = FALSE))
  expect_warning(check_assay_order(c(list_of_assays[5], list_of_assays[1:4]), verbose = TRUE))
  expect_false(check_assay_order(rev(list_of_assays), verbose = FALSE))

  expect_true(check_assay_order(rev(list_of_assays), short_window_period_first = FALSE))
  expect_false(check_assay_order(c(list_of_assays[5], list_of_assays[1:4]), short_window_period_first = FALSE, verbose = FALSE))
  expect_warning(check_assay_order(c(list_of_assays[5], list_of_assays[1:4]), short_window_period_first = FALSE, verbose = TRUE))
  expect_false(check_assay_order(list_of_assays, short_window_period_first = FALSE, verbose = FALSE))
})



