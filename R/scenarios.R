#' Recognize a scenario from an ihist
#'
#' Given a scenario mapping and an ihist, this function will detect what scenario the given ihist represents.
#'
#' Note that this function will apply the 'select_most_informative_results' function to the ihist, so that places some limits on how exotic the scenario description can be.
#'
#' @param ihist The ihist from which the scenario must be detected.
#' @param scenario_mapping A list of lists where the names in the outer list names the scenarios and the inner lists provides the test patterns that are associated with that scenario. If no scenario_mapping is specified, the AMP Taqman/Architect mapping will be used.
#' @export

recognize_scenario <- function(ihist, scenario_mapping = NULL){
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
  return(scenario_mapping)
}
