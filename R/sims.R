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
