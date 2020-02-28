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
    order_ok <- check_assay_order(list_of_assays, verbose = FALSE)
    stopifnot(order_ok)
  }
  
  # draw result for each assay until first negative result
  for (c_assay in list_of_assays){
    assay_dynamics <- all_assay_dynamics[[c_assay]]
  
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
  if (FALSE){ # DEBUGGING notes
    assay_dynamics <- all_assay_dynamics[['geenius_fr_weib3_delaney']]
    list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                        "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                        "geenius_fr_weib3_delaney")
    skip_order_check = TRUE
  }
  sc_times <- list()
  # check that assay list is in correct order
  if (!skip_order_check){
    order_ok <- check_assay_order(list_of_assays, verbose = FALSE)
    stopifnot(order_ok)
  }

  # draw result for each assay
  for (c_assay in list_of_assays){
    assay_dynamics <- all_assay_dynamics[[c_assay]]
  
    evaluate_dynamics <- function(x){ #x is time since infection event
      assay_dynamics$params$x <- x
      return(do.call(assay_dynamics$fun, assay_dynamics$params))
    }

    min_sc <- 0
    sc <- -1
    counter <- 0

    while (sc < min_sc){ # enforce expected assay seroconversion time order
      counter <- counter + 1
      stopifnot(counter < 100)
      u <- runif(1)
      if (evaluate_dynamics(150) < u){
        warning('Simulated seroconversion date more than 150 days after infection date - this is not supported - returning 150. Proceed with caution')
        sc <- 150
      } else {
        sc <- optimize(f = function(y){(evaluate_dynamics(y) - u)^2},
                       interval = c(0, 150),
                       tol = 2.220446e-16)
        sc <- sc$minimum
        min_sc <- sc
        stopifnot(abs(evaluate_dynamics(sc) - u) < 0.02)
      }
    }
  
    sc_times[[c_assay]] <- sc
  }
  return(sc_times)
}

















