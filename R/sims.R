#' Simulates a visit's diagnostic results
#'
#' Given the time since infection, this function will draw from the window period distributions of the specified assays to obtain a set of diagnostic test result.
#'
#' The first negative test result will stop the procedure since then it is known that all other assays will also produce negative results.
#'
#' A key shortcoming of this function is that successive calls are independent. Thus if you want to construct a diagnostic history this function is inappropriate. See the sim-dx-history vignette.
#'
#' @param tsi The time since infection in days. tsi = 0 is the day of infection.
#' @param list_of_assays ORDERED list of assays for which test results should be produced. They must be ordered with the fastest assay first (fastest = assay with the shortest window period)
#' @param skip_order_check When set to TRUE (the default), it is assumed that the order of the list_of_assays is correct. It is EXTREMELY important that this order is correct. This option defaults to TRUE since this step is very slow, so run it the first time only to check that your list is ordered correctly (by setting skip_order_check = FALSE).
#' @export

sim_dx_results <- function(tsi, list_of_assays, skip_order_check = TRUE){

  dx_results <- list()
  
  # check that assay list is in correct order
  if (!skip_order_check){
    for (indx in 1:(length(list_of_assays)-1)){

#      assay1 <- all_assay_dynamics[[list_of_assays[indx]]]
#      assay2 <- all_assay_dynamics[[list_of_assays[indx+1]]]

      assay1 <- get_assay_dynamics(list_of_assays[indx])
      assay2 <- get_assay_dynamics(list_of_assays[indx+1])
      res <- which_is_faster(assay1, assay2)
      if (assay1$short_assayname != res$faster$short_assayname){
        stop('assay list must be from fastest to slowest')
      }
    }
  }
  
  # draw result for each assay until first negative result
  for (c_assay in list_of_assays){
    assay_dynamics <- get_assay_dynamics(c_assay)
  
    evaluate_dynamics <- function(x){ #x is time since infection event
      assay_dynamics$params$x <- x
      return(do.call(assay_dynamics$fun, assay_dynamics$params))
    }
  
    pos_prob <- evaluate_dynamics(tsi)
    if (runif(1) < pos_prob){
      dx_results[[c_assay]] <- '+'
    } else {
      dx_results[[c_assay]] <- '-'
      return(dx_results)
    }
  }
  return(dx_results)
}

#' Simulate time to seroconversion
#'
#' Given a list of assays (ordered according to window period - from short to long) assume that the infecting exposure occurred at time zero and simulate at what time each assay will start to produce positive results. This output can then be combined with a visit schedule to simulate test results for a diagnostic history.
#'
#' @param list_of_assays ORDERED list of assays for which test results should be produced. They must be ordered with the fastest assay first (fastest = assay with the shortest window period)
#' @param skip_order_check When set to TRUE (the default), it is assumed that the order of the list_of_assays is correct. It is EXTREMELY important that this order is correct. This option defaults to TRUE since this step is very slow, so run it the first time only to check that your list is ordered correctly (by setting skip_order_check = FALSE).
#' @export

sim_sc_times <- function(list_of_assays, skip_order_check = TRUE){

}

















