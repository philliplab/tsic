#' Removes strongly dependent results at the same timepoint
#'
#' The construction of the aggregate curve assumes that the contributions from the different tests are independent. For some tests, for example fully reactive and indeterminate versions of the same test this assumption is strongly violated. To offset this, remove the least informative version of the in case that they give the same result.
#'
#' @param ihist The diagnostic history from which to remove such cases
#' @param more_sensitive_test The test which detects infection faster. In case both are positive, this one will be removed.
#' @param less_sensitive_test The test which detects infection slower. In case both are negative, this one will be removed.
#' @export

remove_strongly_dependent_results <- function(ihist, more_sensitive_test, less_sensitive_test){
  if (FALSE){
    devtools::load_all('/home/phillipl/projects/tsic/code/tsic')
    dat <- load_dsmb_nov_2019_data(file_name = '/fridge/data/AMP/DSMB_timing_nov_2019/AMP_diagnostic_testing_history_DSMB_2019_Nov.csv')
    ihist <- subset(dat, ptid == 'p_703-0109')
    less_sensitive_test <- 'geenius_fr_weib3_delaney'
    more_sensitive_test <- 'geenius_indet_weib3_delaney'
  }
  ihist_names <- names(ihist)
  if (more_sensitive_test %in% ihist$test & less_sensitive_test %in% ihist$test){
    ihist$remove_flag <- 0
    c_date <- unique(ihist$sample_date)[3] ####
    for (c_date in unique(ihist$sample_date)){
      c_tests <- subset(ihist, sample_date == c_date & test %in% c(more_sensitive_test, less_sensitive_test))
      if (nrow(c_tests) %in% c(0, 1) ) {next}
      stopifnot(nrow(c_tests) == 2)
      if (all(c_tests$result == '-')){
        ihist <- c_tests$sample_date == c_date & test == less_sensitive_test
        stopifnot(sum(indx) == 1)
        ihist[indx, 'remove_flag'] <- 1
      } else if (all(c_tests$result == '+')) {
        indx <- ihist$sample_date == c_date & ihist$test == more_sensitive_test
        stopifnot(sum(indx) == 1)
        ihist[indx, 'remove_flag'] <- 1
      }
    }
  } else {
    return(ihist)
  }
  ihist <- subset(ihist, remove_flag == 0)
  return(ihist[, ihist_names])
}

#' Selects most informative test results only
#'
#' Given a diagnostic history and an ordered list of tests by window period, this function will return a restricted diagnostic history that only contains the most informative test results for each sample date. The most informative results are the test with the shortest window period that produced a negative result and the test with the longest window period that produced a positive result. The motivation for doing this is to ensure that the aggregate function is constructed using independent results.
#'
#' @param ihist The diagnostic test history of the individual.
#' @param fastest_to_slowest_tests A vector listing the tests in the order from the test with the shortest window period to the test with the longest window period. This ordering will be consulted in conjunction with the test result to decide which tests are the most informative.
#' @export

select_most_informative_results <- function(ihist, fastest_to_slowest_tests = NULL){
  if (FALSE){
    devtools::load_all()
    fastest_to_slowest_tests <- NULL
    for (i in 1:(length(fastest_to_slowest_tests)-1)){
      res <- which_is_faster(all_assay_dynamics[[fastest_to_slowest_tests[i]]], all_assay_dynamics[[fastest_to_slowest_tests[i+1]]])
      print(c("IN", fastest_to_slowest_tests[i], fastest_to_slowest_tests[i+1]))
      print(c("OUT", res$faster$full_assayname, res$slower$full_assayname))
    }
    dat <- load_dsmb_nov_2019_data(file_name = '/fridge/data/AMP/DSMB_timing_nov_2019/AMP_diagnostic_testing_history_DSMB_2019_Nov.csv')
    ihist <- subset(dat, ptid == "p_703-0013")
  } # end of debugging tools
  if (is.null(fastest_to_slowest_tests)){
    fastest_to_slowest_tests <- c("iscav2_weib3_delaney_and_tosiano",
      "aptima_weib3_delaney",
      "taqman_weib3_delaney_and_manufacturer",
      "abbott_real_time_weib3_delaney_and_manufacturer",
      "architect_weib3_delaney",
      "gs_combo_weib3_delaney",
      "determine_weib3_delaney",
      "geenius_indet_weib3_delaney",
      "geenius_fr_weib3_delaney",
      "oraquick_weib3_delaney")
  }
  stopifnot(length(unique(ihist$ptid))==1)
  stopifnot(all(ihist$test %in% fastest_to_slowest_tests))

  kept_ihist <- NULL
  rm_ihist   <- NULL

  for (c_date in sort(unique(ihist$sample_date))){
    c_ihist <- subset(ihist, sample_date == c_date)
    for (c_res in sort(unique(c_ihist$result))){
      if (c_res == '-') {
        test_order <- fastest_to_slowest_tests
      } else {
        test_order <- rev(fastest_to_slowest_tests)
      }
      cc_ihist <- subset(c_ihist, result == c_res)
      for (c_test in test_order){
        if (c_test %in% cc_ihist$test){
          kept_ihist <- rbind(kept_ihist, subset(cc_ihist, test == c_test))
          rm_ihist   <- rbind(rm_ihist,   subset(cc_ihist, test != c_test))
          break
        }
      }
    }
  }
  kept_ihist <- kept_ihist[with(kept_ihist, order(sample_date, test, result)), ]
  rm_ihist <- rm_ihist[with(rm_ihist, order(sample_date, test, result)), ]
  return(list(kept_ihist = kept_ihist,
              rm_ihist   = rm_ihist))
}

