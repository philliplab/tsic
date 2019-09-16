#' Construct an assay result interpretation function
#'
#' Given the dynamics of an assay and a result, produce a function that when evaluated for a potential exposure date will produce the probability of observing the given result at the relevant sample date. This function will only take x as input, which is the potential exposure date, and will already have the particulars of the result and assay loaded into it.
#'
#' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
#' @param result The result of the test. Either '+' or '-'.
#' @param sample_date The date on which the sample was drawn.
#' @export

construct_assay_result_interpreter <- function(assay_dynamics, result, sample_date){
  evaluate_dynamics <- function(x){
    tmp_x <- lubridate::as.duration(lubridate::ymd(sample_date) - lubridate::as_date(x))
    tmp_x <- as.numeric(tmp_x)/(60*60*24)
    assay_dynamics$params$x <- tmp_x
    return(do.call(assay_dynamics$fun, assay_dynamics$params))
  }
  if (result == '+'){
    result_interpreter <- function(x){ return(evaluate_dynamics(x)) }
  } else if (result == "-"){
    result_interpreter <- function(x){ return(1 - evaluate_dynamics(x)) }
  } else {
    stop('result must be "+" or "-"')
  }

  return(result_interpreter)
}

#' Evaluates a function over a range
#'
#' Given a function and either the endpoints of a range of a number of seedpoints, evaluate the function at points in the range so that the change in y does not exceed a certain threshold between any two points, subject to the constraint that no interval between to x points go below some interval.
#'
#' Note that for performance reasons, you want to set up your seed points so that this function will not go through many rounds of recursion.
#'
#' @param fun The function to evaluate.
#' @param seedpoints The points at which the function should be evaluated initially. It is very important to provide the function with a sensible set of seedpoints. Due to the design of the tsic package, such a set of points should always be available.
#' @param max_delta The largest increase in y that is tolerated in a single interval. Default = 0.02. Warning, setting this small relative to your y range can lead to extreme performance requirements.
#' @param min_length The minimum interval length that will still be divided into smaller segments. If the length of an interval is below this number, then it will not be divided into shorter segments no matter how much y changes in it. Warning, depending on how wiggly your function is, setting this number small relative to y can cause extreme performance requirements.
#' @param n_new_segments Into how many segments should each segment that contains too much y movement be divided? You probably want to set this high due to R's terrible recursion performance.
#' @export

get_scatterpoints <- function(fun, seedpoints, max_delta = 0.02, min_length = 0.01, n_new_segments = 20, verbose = FALSE){
  stopifnot(length(seedpoints)>=2)
  x <- seedpoints
  y <- rep(-1, length(x))
  for (i in 1:length(x)){
    y[i] <- fun(x[i])
  }
  if (verbose){
    print('x, y')
    print(x)
    print(y)
  }

  delta <- abs(y[1:(length(y)-1)] - y[2:length(y)])
  int_lengths <- abs(x[1:(length(x)-1)] - x[2:length(x)])
  too_steep <- (delta > max_delta) & (int_lengths > min_length)

  if (verbose){
    print('delta and int_lengths and too_steep')
    print(delta)
    print(int_lengths)
    print(too_steep)
  }

  if (any(too_steep)){
    new_x <- NULL
    new_y <- NULL
    which_too_steep <- which(too_steep)
    prev_indx <- 0
    for (indx in which_too_steep){
      if (verbose){
        print('prev_indx, indx')
        print(prev_indx)
        print(indx)
      }
      new_x <- c(new_x, x[(prev_indx+1):(indx)])
      new_y <- c(new_y, y[(prev_indx+1):(indx)])
      results <- get_scatterpoints(fun = fun, 
                                   seedpoints = seq(from = x[indx], to = x[indx+1], length.out = n_new_segments), 
                                   max_delta = max_delta, 
                                   min_length = min_length, 
                                   n_new_segments = n_new_segments,
                                   verbose = verbose)
      new_x <- c(new_x, results$x[2:(length(results$x)-1)])
      new_y <- c(new_y, results$y[2:(length(results$y)-1)])
      prev_indx <- indx
    }
    last_few <- min(which(x > max(new_x)))
    if (verbose){
      print('which(x > max(new_x)) and last_few and length(x)')
      print(which(x > max(new_x)))
      print(last_few)
      print(length(x))
    }
    if (last_few < Inf){
      new_x <- c(new_x, x[last_few:length(x)])
      new_y <- c(new_y, y[last_few:length(y)])
    }
  } else {
    new_x <- x
    new_y <- y
  }
  return(list(x = new_x, y = new_y))
}

#' Reduce the number of scatter points
#'
#' Given a set of scatterpoints, remove those points that does not lead to an improvement in the accuracy of a plot. For now this is just if you take three points, and all their y values are identical, then remove the middle point. Identical is tested by ensuring that the absolute difference is below a tolerance parameter called min_delta.
#'
#' TODO: If three points lie on the same line, then the middle point adds no information. Expand this function to test this instead of just checking that the y values are the same. How will this affect the aggregated function? First implement that before considering making this improvement.
#'
#' NOTE: This algorithm can be further improved. If any of three points do not share the same y value, all three points will be kept. Usually, only two needs to be kept. This will result in a minor improvement that is not worth the effort at present.
#'
#' @param min_delta The smallest increase in y that is required before it is flagged as a real change that is worth keeping.
#' @param max_length The max length that an interval between two points on the x-axis will be allowed to grow to. If an interval is longer than this, then no other intervals will be to appended to it anymore.
#' @export

reduce_x_points <- function(x, y, min_delta = 0.0001, max_length = 14){
  stopifnot(length(x) == length(y))
  indx <- 1
  new_x <- rep(NA_real_, length(x))
  new_y <- rep(NA_real_, length(y))
  new_indx <- 2
  new_x[1] <- x[1]
  new_y[1] <- y[1]
  while (indx <= length(x)-2){
    if ((abs(y[indx] - y[indx+1]) < min_delta) & 
        (abs(y[indx] - y[indx+2]) < min_delta) & 
        (abs(y[indx+1] - y[indx+2]) < min_delta) &
        (abs(x[indx] - new_x[new_indx-1] + 1) < max_length)){
      indx <- indx + 1
    } else {
      new_x[new_indx:(new_indx+2)] <- x[indx:(indx+2)]
      new_y[new_indx:(new_indx+2)] <- y[indx:(indx+2)]
      new_indx <- new_indx + 3
      indx <- indx + 2
    }
    indx <- indx + 1
  }
  if (new_x[new_indx-1] != x[length(x)]){
    new_x[new_indx] <- x[length(x)]
    new_y[new_indx] <- y[length(y)]
  } else {
    new_indx <- new_indx-1
  }
  new_x <- new_x[1:new_indx]
  new_y <- new_y[1:new_indx]
  return(list(x = new_x, y = new_y))
}

#' Constructs an interpreter for a diagnostic history
#'
#' Construct a closure that gives the likelihood of a diagnostic history if the infection event was on the date supplied.
#'
#' @param ihist A data.frame with columns sample_date, test, and result that contains a diagnostic history (or a subset thereof) for a single person.
#' @export

construct_aggregate_interpreter <- function(ihist){
  stopifnot(class(ihist) == 'data.frame')
  stopifnot(all(names(ihist) == c('ptid', 'sample_date', 'test', 'result')))
  stopifnot(nrow(ihist) >= 1)
  assay_interpreters <- list()
  for (i in 1:nrow(ihist)){
    c_dynamics <- get_assay_dynamics(assay = ihist$test[i])
    c_interpreter <- construct_assay_result_interpreter(assay_dynamics = c_dynamics, 
                                                        result = ihist$result[i], 
                                                        sample_date = ihist$sample_date[i] )
    assay_interpreters[[i]] <- c_interpreter
    environment(assay_interpreters[[i]])$result <- ihist$result[i]
    environment(assay_interpreters[[i]])$assay_dynamics <- c_dynamics
    environment(assay_interpreters[[i]])$sample_date <- ihist$sample_date[i]
    #    print(ls(environment(assay_interpreters[[i]])))
  }
#  print(assay_interpreters)
#  for (i in 1:length(assay_interpreters)){
#    for (c_var in ls(environment(assay_interpreters[[i]]))){
#      print (c('Assay interpter', i, ' variable ', c_var))
#      print (get(c_var, envir = environment(assay_interpreters[[i]])))
#    }
#  }
  evaluate_proposed_infection_date <- function(x){
    overall_result <- 1
    for (i in 1:length(assay_interpreters)){
      c_result <- assay_interpreters[[i]](x)
#      print(c('Assay result # ', i, ' at time ', x, ' is equal to ', c_result))
      overall_result <- overall_result * c_result
    }
    return(overall_result)
  }
  return(evaluate_proposed_infection_date)
}

#' Interprets an ihist into daily likelihoods
#'
#' Given a range of interest and an individual's diagnostic history, compute the likelihood that each day in the range was the day of infection. This produces the output in long format.
#'
#' @param ihist The diagnostic history for an individual.
#' @param range_start The earliest day of interest.
#' @param range_end The latest day of interest.
#' @param verbose Should verbose output be printed during the computation?
#' @example
#' ihist <- data.frame(
#'   ptid = c('p0', 'p0'),
#'   sample_date = c('2016-03-01', '2016-09-01'),
#'   test = c('step_unit_testing', 'step_unit_testing'),
#'   result = c('-', '+'),
#'   stringsAsFactors = FALSE
#' )
#' range_start <- as.Date('2016-01-01')
#' range_end <- as.Date('2016-12-01')
#' interpret_ihist(ihist = ihist, range_start = range_start, range_end = range_end, verbose = TRUE)
#' @export

interpret_ihist <- function(ihist, range_start, range_end, verbose = FALSE){
  stopifnot(class(range_start) == 'Date')
  stopifnot(class(range_end) == 'Date')
  stopifnot(length(unique(ihist$ptid)) == 1)

  # Debugging stuff
  if (FALSE){
    ihist <- data.frame(
      ptid = c('p0', 'p0'),
      sample_date = c('2016-03-01', '2016-09-01'),
      test = c('step_unit_testing', 'step_unit_testing'),
      result = c('-', '+'),
      stringsAsFactors = FALSE
    )
    range_start <- as.Date('2016-01-01')
    range_end <- as.Date('2016-12-01')
    verbose <- FALSE
    i <- 1
  }


  all_xy_points <- data.frame(ptid = character(0),
                              sample_date = numeric(0),
                              test_details = character(0),
                              prob_val = numeric(0),
                              stringsAsFactors = TRUE)

  for (i in 1:nrow(ihist)){
    if (verbose){
      cat(paste0(ihist[i, 'test'], '_', ihist[i, 'sample_date'], '_', ihist[i, 'result'], ': Setting up'))
    }
    assay_dynamics <- get_assay_dynamics(ihist[i, 'test'])
    assay_interpreter <-
      construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                         result = ihist[i, 'result'],
                                         sample_date = ihist[i, 'sample_date'])

    if (verbose){ cat(paste0('. Getting')) }
    c_xy_points <- get_scatterpoints(fun = assay_interpreter,
                                          seedpoints = range_start:range_end)
    if (verbose){ cat(paste0('. Reducing')) }
    c_xy_points <- reduce_x_points(x = c_xy_points[['x']],
                                        y = c_xy_points[['y']])
    if (verbose){ cat(paste0('. Binding')) }
    all_xy_points <- rbind(
      all_xy_points,
      data.frame(ptid = ihist[i, 'ptid'],
                 sample_date = c_xy_points[['x']],
                 test_details = paste0(ihist[i, 'test'], '\n', ihist[i, 'sample_date'], '\n', ihist[i, 'result']),
                 prob_val = c_xy_points[['y']],
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
    if (verbose){cat('.\n')}
  }
#construct_aggregate_interpreter <- function(ihist){
  if (verbose){
    cat(paste0('Aggregate (Slow step)', ': Setting up'))
  }
  aggregate_interpreter <- construct_aggregate_interpreter(ihist)
  aggregate_seedpoints <- sort(unique(all_xy_points[['sample_date']]))
  if (verbose){ cat(paste0('. Getting')) }
  c_xy_points <- get_scatterpoints(fun = aggregate_interpreter,
                                   seedpoints = aggregate_seedpoints)
  if (verbose){ cat(paste0('. Reducing')) }
  c_xy_points <- reduce_x_points(x = c_xy_points[['x']],
                                      y = c_xy_points[['y']])
  if (verbose){ cat(paste0('. Binding')) }
  all_xy_points <- rbind(
    all_xy_points,
    data.frame(ptid = ihist[1, 'ptid'],
               sample_date = c_xy_points[['x']],
               test_details = paste0('Aggregate'),
               prob_val = c_xy_points[['y']],
               stringsAsFactors = FALSE),
    stringsAsFactors = FALSE)
  if (verbose){cat('.\n')}

  return(all_xy_points)
}


#################### OBSOLETE

##' Probabilities of certain test results
##'
##' Turns assay dynamics into the probability of observing a specified test result a given number of days since XXX
##'
##' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
##' @param result '+' or '-' indicating the test result.
##' @param days_to_visit The number of days between XXX and the test result.
##' @export
#
#prob_test_result_given_days_to_visit <- function(assay_dynamics, result, days_to_visit){
#  assay_dynamics$params$x <- days_to_visit
#  x <- do.call(assay_dynamics$fun, assay_dynamics$params)
#
#  if (result == "+"){
#    return(x)
#  } else if (result == "-"){
#    return(1-x)
#  } else {
#    stop('result must be "+" or "-"')
#  }
#}
#
##' Probability of DDI_1 curve implied by test result
##'
##' Compute a curve of probabilities of observing a specified test result on a specified date assuming the day of DDI_1. Thus assume individually that each day in a range of dates is the DDI_1 and then compute the probability of observing the specified result.
##'
##' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
##' @param result '+' or '-' indicating the test result.
##' @param range_start The date furthest into the past that should be considered as a plausible day for DDI_1.
##' @param range_end The most recent date that should be considered as a plausible day for DDI_1.
##' @param visit_date The date of the visit on whose sample the test result was observed
##' @param min_prob The most 'certain' you are that the given test results excludes dates far from the actual visit. Think of this as a chance that the test may give an incorrect result.
##' @param max_prob The most 'certain' you are that the given test results excludes dates far from the actual visit. Think of this as a chance that the test may give an incorrect result. TODO: Think more about what this actually means.
##' @export
#
#interpret_result <- function(assay_dynamics, result, range_start, range_end, visit_date, 
#                             min_prob = 0.00, max_prob = 1.00){
#  days_vec <- NULL
#  prob_vec <- NULL
#  d2v_vec <- NULL
#  for (i in ymd(range_start):ymd(range_end)){
#    time_to_visit = as.duration(ymd(visit_date) - as_date(i))
#    days_to_visit = as.numeric(time_to_visit)/(60*60*24)
#    stopifnot(days_to_visit == floor(days_to_visit))
#
#    days_vec <- c(days_vec, as_date(i))
#    prob_vec <- c(prob_vec, 
#                  prob_test_result_given_days_to_visit(assay_dynamics = assay_dynamics, 
#                                                       result = result, 
#                                                       days_to_visit = days_to_visit)
#                  )
#    d2v_vec <- c(d2v_vec, days_to_visit)
#  }
#  prob_vec[prob_vec < min_prob] <- min_prob
#  prob_vec[prob_vec > max_prob] <- max_prob
#  return(data.frame(days = as_date(days_vec),
#                    prob = prob_vec,
#                    d2v = d2v_vec)
#        )
#}
#
##' Given a set of results, construct DDI_1 probability curves
##'
##' Calls interpret_result on a set of test results and aggregate the results into a data.frame
##'
##' @param result_set A data.frame of all the test results and dates for a patient
##' @param range_start The date furthest into the past that should be considered as a plausible day for DDI_1. Can be NULL, then defaults to 2 months before the earliest result.
##' @param range_end The most recent date that should be considered as a plausible day for DDI_1. Can be NULL, then defaults to 1 week after the most recent result.
##' @export
#
#
#
#interpret_result_set <- function(result_set, range_start = NULL, range_end = NULL){
#  if (is.null(range_start)){
#    range_start <- min(ymd(result_set$visit_date)) %m-% months(2)
#  }
#  if (is.null(range_end)){
#    range_end <- max(ymd(result_set$visit_date)) %m+% weeks(1)
#  }
#
#  dat <- data.frame(date = as_date(ymd(range_start):ymd(range_end)))
#  for (r_indx in 1:nrow(result_set)){
#    c_result <- result_set[r_indx,]
#    assay_dynamics <- get_assay_dynamics(assay = c_result$assay)
#    i_result <- interpret_result(assay_dynamics = assay_dynamics,
#                                 result = c_result$result,
#                                 range_start = range_start,
#                                 range_end = range_end,
#                                 visit_date = c_result$visit_date)
#    p_result <- i_result[,c('days', 'prob')]
#    tmpx <- paste(c_result$assay, '_', 
#                  gsub("-", "", c_result$visit_date), 
#                  '_', c_result$result, sep = '')
#    names(p_result) <- c( 'date', tmpx )
#    dat <- merge(dat, p_result)
#  }
#  return(dat)
#}
#
