#' Recognize a scenario from an ihist
#'
#' Given a scenario mapping and an ihist, this function will detect what scenario the given ihist represents.
#'
#' Note that this function will apply the 'select_most_informative_results' function to the ihist, so that places some limits on how exotic the scenario description can be.
#'
#' TODO: This function will be much more useful if it is extended to handle different assays that may produce the same result. For example, allowing the scenario 1 case to be defined by a negative on with Taqman or Abbott Real Time.
#'
#' @param ihist The ihist from which the scenario must be detected.
#' @param scenario_mapping A list of lists where the names in the outer list names the scenarios and the inner lists provides the test patterns that are associated with that scenario. If no scenario_mapping is specified, the AMP Taqman/Architect mapping will be used.
#' @export

recognize_scenario <- function(ihist, scenario_mapping = NULL){
  if (FALSE){
    list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
                        "architect_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney")
    sc_times <- sim_sc_times(list_of_assays, fix_draw = 0.5)
    
    f_ihist_s1 <- combine_sc_and_visit_times(sc_times,
                                             ((0:2)*28-46),
                                             true_infection_date = as.numeric(as.Date('2019-03-01'))+46.5)
    m_ihist_s1 <- select_most_informative_results(f_ihist_s1)$kept_ihist
    ihist <- f_ihist_s1
  }
  if (is.null(scenario_mapping)){
    scenario_mapping <- list()
    scenario_mapping[['Scenario 1']] <- list("iscav2_weib3_delaney_and_tosiano" = '+',
                                             "taqman_weib3_delaney_and_manufacturer" = '-')
    scenario_mapping[['Scenario 2']] <- list("taqman_weib3_delaney_and_manufacturer" = '+',
                                             "architect_weib3_delaney" = '-')
    scenario_mapping[['Scenario 3']] <- list("architect_weib3_delaney" = '+',
                                             "geenius_indet_weib3_delaney" = '-')
    scenario_mapping[['Scenario 4']] <- list("geenius_indet_weib3_delaney" = '+',
                                             "geenius_fr_weib3_delaney" = '-')
    scenario_mapping[['Scenario 5']] <- list("geenius_fr_weib3_delaney" = '+')
  }
  # get FP
  fp_date <- min(subset(ihist, result == '+')$sample_date)
  ln_date <- max(ihist$sample_date[ihist$sample_date < fp_date])

  # get most inform only
  m_ihist <- select_most_informative_results(ihist)$kept_ihist

  # scan mapping
  matching_scenarios <- NULL
  fpm_ihist <- subset(m_ihist, sample_date == fp_date)
  for (c_scenario_name in names(scenario_mapping)){
    c_scenario_spec <- scenario_mapping[[c_scenario_name]]
    it_is_this_scenario <- TRUE
    for (c_test in names(c_scenario_spec)){
      if (c_test %in% fpm_ihist$test){
        if (c_scenario_spec[[c_test]] == fpm_ihist$result[fpm_ihist$test == c_test]){
          it_is_this_scenario <- it_is_this_scenario & TRUE
        } else {
          it_is_this_scenario <- FALSE
        }
      } else {
        it_is_this_scenario <- FALSE
      }
    }
    if (it_is_this_scenario) {
      matching_scenarios <- c(matching_scenarios, c_scenario_name)
    }
  }
  stopifnot(length(matching_scenarios) < 2)
  return(list('matching_scenarios' = matching_scenarios,
              'lnfp_gap' = fp_date - ln_date))
}
