#' Simulates a visit's diagnostic results
#'
#' Given the time since infection, this function will draw from the window period distributions of the specified assays to obtain a set of diagnostic test result.
#'
#' The first negative test result will stop the procedure since then it is known that all other assays will also produce negative results.
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






#' Simulates the durations times between fiebig stages
#'
#' For each patients, draw an eclispe duration from the supplied eclipse_distribution,
#' and stage durations from exponential distributions with parameters specified by the rates argument.
#'
#' Only the durations in stages I to V is drawn since the exit form stage VI is not defined / measured.
#' 
#' @return data.frame with columns for ptid and the drations of the eclipse phase and stages I to V.
#'
#' @param n The number of people to simulate
#' @param rates The parameters for the exponential distributions of each of the six stages. ( default = 'fiebig' which gets replaced with list(l1 = 5.0, l2 = 5.3, l3 = 3.2, l4 = 5.6, l5 = 69.5) )
#' @param eclipse_distribution A list with elements distr pointing to the function allowing draws from the eclipse distribution ( default = function(n, shape, scale, location) {location + rweibull(n, shape = shape, scale = scale)} ) and the parameters of said distribution ( default = list(shape = 1.35, scale = 9, location = 4.8) )
#' @export

sim_durations <- function(n = 30, rates = 'fiebig',
                            eclipse_distribution = list(
                              distr = function(n, shape, scale, location) {location + rweibull(n, shape = shape, scale = scale)},
                              params = list(shape = 1.35, scale = 9, location = 4.8))
                            ){
  if (FALSE){
    devtools::load_all()
    n <- 30
    rates <- 'fiebig'
    eclipse_distribution <- list(
      distr = function(n, shape, scale, location) {location + rweibull(n, shape = shape, scale = scale)},
      params = list(shape = 1.35, scale = 9, location = 4.8))
  }
  durations <- data.frame(ptid = 1:n)
  
  eclipse_distribution$params$n <- n
  durations$eclipse <- do.call(eclipse_distribution$distr, eclipse_distribution$params)

  if (class(rates) == 'character'){
    stopifnot(rates %in% c('fiebig'))
    if (rates == 'fiebig'){
      rates <- list(l1 = 5.0, l2 = 5.3, l3 = 3.2, l4 = 5.6, l5 = 69.5)
    }
  }
  stopifnot(all(sort(names(rates)) == c('l1', 'l2', 'l3', 'l4', 'l5')))

  for (c_stage in sort(names(rates))){
    durations[, gsub('l', 'f', c_stage)] <- rexp(n, 1/rates[[c_stage]])
  }
  return(durations)
}

#' Simulates the visit dates around a set of transitions
#'
#' @param durations The data.frame of transitions as produced by sim_durations.
#' @param inf_date The infection date. In this version, all infections occur on the same date. This makes comparisions easy also.
#' @param n_before_infection The number of visits before infection.
#' @param gap_distribution A function that generates the n gaps between the next n visits.
#' @export

sim_visit_dates <- function(durations, 
                            inf_date = as.numeric(as.Date('2015-05-01')),
                            n_before_infection = 2,
                            gap_distribution = function(n){rnorm(n, 30, 3)}){
  if (FALSE){
    durations <- sim_durations()
    inf_date <- as.numeric(as.Date('2015-05-01'))
    n_before_infection <- 2
    gap_distribution <- function(n, mu = 30, sig = 3){rnorm(n, mu, sig)}
    c_ptid <- 1
  }
  all_visit_dates <- NULL

  for (c_ptid in sort(unique(durations$ptid))){
    pot_gaps <- c(0, cumsum(gap_distribution(max(n_before_infection + 5, 10))))
    inf_offset <- runif(1, 1, pot_gaps[n_before_infection] - 1)
    v1_date <- inf_date - (inf_offset + pot_gaps[n_before_infection])
    visit_dates <- v1_date + pot_gaps
    last_transition <- inf_date + sum(durations[durations$ptid == c_ptid, -1])

    transitions_into <- c(0, as.numeric(durations[durations$ptid == c_ptid, -1]))
    names(transitions_into) <- c('eclipse', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6')
    
    while (last_transition > max(visit_dates)){
      pot_gaps <- cumsum(gap_distribution(max(n_before_infection + 5, 10)))
      visit_dates <- c(visit_dates, max(visit_dates) + pot_gaps)
    }
    all_visit_dates <- rbind(all_visit_dates,
      data.frame(ptid = c_ptid,
                 visit_id = 1:length(visit_dates),
                 visit_date = sort(visit_dates)))
  }
  return(all_visit_dates)
}

#' Reads test results from a duration and visit schedule
#'
#' @param durations The data.frame of transitions as produced by sim_durations.
#' @param visit_schedule The data.frame with the date of each visit as produced by sim_visit_dates.
#' @param tests_and_properties A list with each element corresponding to a test giving the Fiebig stage in which the test first starts returning positive results. Additionally, each test has an offset that allows the user to specify that the test starts testing positive x days before or after person enters the stage.
#' @export

sim_test_results <- function(durations, visit_schedule, 
  tests_and_properties = list(aptima_weib3_delaney = list(name = aptima_weib3_delaney,
                                                          stage = 1,
                                                          offset = -2),
                              architech_weib3_delaney = list(name = architech_weib3_delaney,
                                                             stage = 2,
                                                             offset = -1),
                              geenius_indet_weib3_delaney = list(name = geenius_indet_weib3_delaney,
                                                              stage = 3,
                                                              offset = 0),
                              geenius_fr_weib3_delaney = list(name = geenius_fr_weib3_delaney,
                                                              stage = 5,
                                                              offset = 0)
                              )){
    return(1)
}





