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
#' @param fix_draw Provides the option to fix the draws so that it is no longer random. Useful for constructing examples and systematic evaluations.
#' @export

sim_sc_times <- function(list_of_assays, skip_order_check = TRUE, fix_draw = NULL){
  if (FALSE){ # DEBUGGING notes
    devtools::load_all()
    assay_dynamics <- all_assay_dynamics[['geenius_fr_weib3_delaney']]
    list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                        "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                        "geenius_fr_weib3_delaney")
    skip_order_check = TRUE
  }
  # check that assay list is in correct order
  if (!skip_order_check){
    order_ok <- check_assay_order(list_of_assays, verbose = FALSE)
    stopifnot(order_ok)
  }
  if (!is.null(fix_draw)){
    stopifnot(fix_draw<1 & fix_draw>0)
  }

  # draw result for each assay
  sc_times <- list()
  min_sc <- 0
  for (c_assay in list_of_assays){
    assay_dynamics <- all_assay_dynamics[[c_assay]]
  
    evaluate_dynamics <- function(x){ #x is time since infection event
      assay_dynamics$params$x <- x
      return(do.call(assay_dynamics$fun, assay_dynamics$params))
    }

    sc <- -1
    counter <- 0
    r_unif_lb <- 0

    while (sc < min_sc){ # enforce expected assay seroconversion time order
      counter <- counter + 1
      stopifnot(counter < 100)
      if (is.null(fix_draw)){
        u <- runif(1, min = r_unif_lb) #use r_unif_lb to reduce number of iterations before valid sc pattern.
      } else {
        u <- fix_draw
      }
      if (evaluate_dynamics(150) < u){
        warning('Simulated seroconversion date more than 150 days after infection date - this is not supported - returning 150. Proceed with caution')
        sc <- 150
      } else {
        sc <- optimize(f = function(y){(evaluate_dynamics(y) - u)^2},
                       interval = c(0, 150),
                       tol = 2.220446e-16)
        sc <- sc$minimum
        stopifnot(abs(evaluate_dynamics(sc) - u) < 0.02)
      }
      r_unif_lb <- u
    }
    if (is.null(fix_draw)){
      min_sc <- sc
    } else {
      min_sc <- 0
    }
  
    sc_times[[c_assay]] <- sc
  }
  return(sc_times)
}

#' Produce ihist from sc and visit times
#'
#' Given a list of sc times as produced by sim_sc_times, and a vector of visit times relative to the time of true infection (t = 0 is the time of true infection) produce an ihist. This is done by comparing the visit time with the seroconversion times and deducing the result for the relevant test at the visit. Lastly, the dates are converted so that the true infection time is at an appropriate calendar date and not t = 0.
#'
#' TODO: Should there be an option that places the first positive visit at t = 0?
#'
#' @param sc_times A list with names giving the assay names and elements the seroconversion time of that assay relative to the true infection time as produced by sim_sc_times.
#' @param visit_times A vector of visit times relative to the time of true infection, with true infection at time = 0. Thus -5 indicates a visit 5 days before true infection.
#' @param true_infection_date The true of the true infection expressed as the number of days since 1970-01-01. This will be used to shift all the dates in the resulting ihist.
#' @param ptid The value to assign to the ptid column of the resulting ihist.
#' @export

combine_sc_and_visit_times <- function(sc_times, visit_times, true_infection_date = 17000.5, ptid = '?'){
  if (FALSE){ #debugging notes
    devtools::load_all()
    visit_times <- (-3:10)*5 + 2
    list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer", 
                        "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                        "geenius_fr_weib3_delaney")
    sc_times <- sim_sc_times(list_of_assays)
    true_infection_date = 17000.5
    ptid = '?'
  }
  ihist <- NULL
  for (c_visit_time in visit_times){
    c_ihist <- NULL
    for (c_test in names(sc_times)){
      if (c_visit_time > sc_times[[c_test]]){
        res <- '+'
      } else {
        res <- '-'
      }
      c_ihist <- rbind(c_ihist,
        data.frame(ptid = ptid,
                   sample_date = c_visit_time + true_infection_date,
                   test = c_test,
                   result = res,
                   stringsAsFactors = FALSE),
        stringsAsFactors = FALSE)
    }
    ihist <- rbind(ihist, c_ihist)
  }
  return(ihist)
}
















