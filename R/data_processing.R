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

#' Annotates an ihist with visit_numbers and rel_sample_dates
#'
#' Given an ihist with the four main columns, add a visit_number column with the first visit at which any test was positive as 0. The visit immediately before the first positive visit is visit_number -1 while the visit immediately after the first positive visit is visit_number 1.
#'
#' @param ihist An ihist.
#' @export

visit_labeller <- function(ihist){
  ihist$rel_sample_date <- NA_real_
  ihist$visit_number <- NA_real_
  fp_date <- min(subset(ihist, result == '+')$sample_date)
  uniq_sample_dates <- sort(unique(ihist$sample_date))
  rel_sample_dates <- uniq_sample_dates - fp_date
  visit_numbers <- order(rel_sample_dates)
  visit_numbers <- visit_numbers - which(rel_sample_dates == 0)
  for (indx in 1:length(uniq_sample_dates)){
    ihist_indx <- ihist$sample_date == uniq_sample_dates[indx]
    ihist$rel_sample_date[ihist_indx] <- rel_sample_dates[indx]
    ihist$visit_number[ihist_indx] <- visit_numbers[indx]
  }
  return(ihist)
}

#' Computes a summary of the results of a single assay for an ihist
#'
#' Given an ihist and an assay of interest, compute the observed LN and FP dates for that assay. Optionally, slower (faster) assays are inspected to see if earlier (later) LN (FP) dates can be deduced (only useful in the case where some tests are not performed at some visits).
#'
#' @return A single row of a data.frame with the columns:
#' \itemize{
#'   \item ptid - The participant identifier read from the ihist
#'   \item assay_name - The name of the assay this row describes. From the assay_name input to the function.
#'   \item obs_rel_LN - Observed LN date for the assay relative to the overall first positive for the participant. 
#'   \item ded_rel_LN - Deduced LN date (from faster assays) for this assay relative to the overall first positive for the participant. 
#'   \item LN_ded_from - The assay name from which the LN was deduced. Will always be the fastest assay that produced the earliest LN date.
#'   \item obs_rel_FP - Observed FP date for the assay relative to the overall first positive for the participant. 
#'   \item ded_rel_FP - Deduced FP date (from slower assays) for this assay relative to the overall first positive for the participant. 
#'   \item FP_ded_from - The assay name from which the FP was deduced. Will always be the fastest assay that produced the earliest FP date.
#'   \item rel_LN - The maximum of the observed and deduced LN dates. This will ensure that the interval during which the assay transitioned from negative to positive is as short as possible by factoring in the structural knowledge we assume about the assays.
#'   \item rel_FP - The minimun of the observed and deduced FP dates. This will ensure that the interval during which the assay transitioned from negative to positive is as short as possible by factoring in the structural knowledge we assume about the assays.
#'   \item LN_FP_gap - The width of the interval between which the assay transitioned from negative to positive.
#' }
#'
#' @param ihist The individual diagnostic history from which an assay's results must be summarized.
#' @param assay_name The name of the assay to summarize.
#' @param faster_assays The assays with shorter window periods than this assay that should be checked if a later LN result can be deduced.
#' @return

compute_result_summary <- function(ihist, assay_name, faster_assays = NULL){
  if (FALSE){
    devtools::load_all('..')
    list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                        "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                        "geenius_fr_weib3_delaney")
    sc_times <- sim_sc_times(list_of_assays)
    ihist_src <- combine_sc_and_visit_times(sc_times = sc_times,
                                        visit_times = c(-14, 14, 22, 30),
                                        true_infection_date = as.numeric(as.Date('2018-05-05')),
                                        ptid = 'p001')
    ihist_src <- structure(list(ptid = c("p001", "p001", "p001", "p001", "p001",
"p001", "p001", "p001", "p001", "p001", "p001", "p001", "p001",
"p001", "p001", "p001", "p001", "p001", "p001", "p001"), sample_date = c(17642,
17642, 17642, 17642, 17642, 17670, 17670, 17670, 17670, 17670,
17678, 17678, 17678, 17678, 17678, 17686, 17686, 17686, 17686,
17686), test = c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
"architect_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney",
"iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
"architect_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney",
"iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
"architect_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney",
"iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
"architect_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney"
), result = c("-", "-", "-", "-", "-", "+", "+", "-", "-", "-",
"+", "+", "-", "-", "-", "+", "+", "+", "+", "+")), row.names = c(NA,
-20L), class = "data.frame")
    ihist_src <- visit_labeller(ihist_src)

    ihist <- ihist_src

    assay_name <- 'gs_combo_weib3_delaney'
    slower_assays <- c("geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney")
    faster_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer")

    ihist <- ihist_src
    assay_name <- 'like_geenius_indet'
    slower_assays <- c("geenius_fr_weib3_delaney")
    faster_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 'architect_weib3_delaney')

  }

  warning('in progress - some things still need to be figured out')
  fp_date <- min(subset(ihist, test == assay_name & result == '+')$sample_date)
  ln_date <- max(subset(ihist, test == assay_name & sample_date < fp_date)$sample_date)
  rel_fp_date <- min(subset(ihist, test == assay_name & result == '+')$rel_sample_date)
  rel_ln_date <- max(subset(ihist, test == assay_name & rel_sample_date < rel_fp_date)$rel_sample_date)
  
  # compute the lastest faster LN date
  rel_latest_faster_ln_date <- -Inf
  for (c_assay in faster_assays){
    rel_c_ln_date <- max(subset(ihist, test == c_assay & result == '-')$rel_sample_date)
    if (rel_latest_faster_ln_date < rel_c_ln_date){
      rel_latest_faster_ln_date <- rel_c_ln_date
    }
  }

  # check that fp date is reasonable given other fp dates
  faster_assays <- c(faster_assays, assay_name)
  assay_sane <- TRUE
  for (c_assay in faster_assays){
    c_fp_date <- min(subset(ihist, test == c_assay & result == '+')$sample_date)
    assay_sane <- assay_sane & (max(subset(ihist, test == c_assay & result == '-')$sample_date) < c_fp_date)
  }

  assay_not_inf <- fp_date != Inf & ln_date != -Inf
  return(
         data.frame(ptid = unique(ihist$ptid),
                    test = assay_name,
                    assay_sane = assay_sane,
                    ln_date = ln_date,
                    fp_date = fp_date,
                    ln_fp_gap = rel_fp_date - max(rel_ln_date, rel_latest_faster_ln_date),
                    rel_ln_date = rel_ln_date,
                    rel_fp_date = rel_fp_date,
                    rel_latest_faster_ln_date = rel_latest_faster_ln_date,
                    useful_rel_ln_date = max(rel_ln_date, rel_latest_faster_ln_date),
                    assay_not_inf = assay_not_inf,
                    stringsAsFactors = FALSE)
         )
}
