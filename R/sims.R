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
    durations[, c_stage] <- rexp(n, 1/rates[[c_stage]])
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






