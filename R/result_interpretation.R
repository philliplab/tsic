#' Construct an assay result interpretation function
#'
#' Given the dynamics of an assay and a result with a date (in days since 1970-01-01), produce a function that when evaluated for a potential exposure date (in days since 1970-01-01) will produce the probability of observing the given result at the relevant sample date. This function will only take x as input, which is the potential exposure date (in days since 1970-01-01), and will already have the particulars of the result and assay loaded into it.
#'
#' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
#' @param result The result of the test. Either '+' or '-'.
#' @param sample_date The date (in days since 1970-01-01) on which the sample was drawn.
#' @export

construct_assay_result_interpreter <- function(assay_dynamics, result, sample_date){
  if (sample_date %% 1 != 0.5) {warning('Sample date not at midday')}
  evaluate_dynamics <- function(x){
    assay_dynamics$params$x <- sample_date - x
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

  #!DO NOT MOVE! - this infix function must be declared here to have the correct value for min_delta
  `%~=%` <- function(a, b){
    result <- NULL
    if ( (a %in% c(0,1)) | (b %in% c(0,1)) ){
      if (a != b){
        result <- FALSE
      } else {
        result <- TRUE
      }
    } else {
      if (abs(a - b) < min_delta){
        result <- TRUE
      } else {
        result <- FALSE
      }
    }
    return(result)
  } # end %~=%

  prev_is_end <- TRUE # to ensure indx == 1 get's flagged as a start
  new_x <- rep(NA_real_, length(x))
  new_y <- rep(NA_real_, length(x))
  new_indx <- 1

  for (indx in 1:length(x)){
    this_is_start <- FALSE
    this_is_end <- FALSE
    too_long <- FALSE
    # checking: is start?
    if (prev_is_end) { this_is_start <- TRUE }
    else if ( !(y[indx] %~=% y[indx - 1]) ) { this_is_start <- TRUE }
    # is end?
    if (indx == length(x)) { this_is_end <- TRUE }
    else if ( !(y[indx] %~=% y[indx + 1]) ){ this_is_end <- TRUE }
    # too long?
    if (!(indx %in% c(1, length(x)) ) ){
      if (abs(new_x[new_indx-1] - x[indx + 1]) > max_length) {too_long <- TRUE}
    }
    # doing
    if (this_is_start | this_is_end | too_long) {
      new_x[new_indx] <- x[indx]
      new_y[new_indx] <- y[indx]
      new_indx <- new_indx + 1
    }
    if (this_is_end){
      prev_is_end <- TRUE
    } else {
      prev_is_end <- FALSE
    }
  }
  return(list(x = new_x[1:(new_indx - 1)],
              y = new_y[1:(new_indx - 1)]))
}


#' Trims range on which to evaluate a function
#'
#' Because the duration of the trial is long and often the region in which the infection could plausibly have happened is small, the integratoin functions often produces zero. Since we know that any aggregate function that produces a probable infection interval shorted than 1 day is problematic, we can devise an easy scheme to remove large chunks from the range where the interpreter function is zero.
#'
#' Simply walk from the start to end in one day increments, if the likelihood is still zero, adjust the range_start to the current day. Likewise from range_end.
#'
#' @param fun The function that will dictate which portions of the range will get trimmed.
#' @param range_start The lower bound of the initial range.
#' @param range_end The upper bound of the initial range.
#' @param tol How close to zero before we consider it equal to zero? Default = 0.
#' @export

trim_range <- function(fun, range_start, range_end, tol = 0){
  new_range_start <- range_start
  new_range_end <- range_end
  for (i in range_start:range_end){
    if (fun(i) <= tol){
      new_range_start <- i
    } else {
      break
    }
  }
  if (new_range_start < range_end){
    for (i in range_end:new_range_start){
      if (fun(i) <= tol){
        new_range_end <- i
      } else {
        break
      }
    }
  }
  return(list(range_start = new_range_start,
              range_end = new_range_end))
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
    c_dynamics <- all_assay_dynamics[[assay = ihist$test[i]]]
    c_interpreter <- construct_assay_result_interpreter(assay_dynamics = c_dynamics, 
                                                        result = ihist$result[i], 
                                                        sample_date = ihist$sample_date[i] )
    assay_interpreters[[i]] <- c_interpreter
    environment(assay_interpreters[[i]])$result <- ihist$result[i]
    environment(assay_interpreters[[i]])$assay_dynamics <- c_dynamics
    environment(assay_interpreters[[i]])$sample_date <- ihist$sample_date[i]
  }
  evaluate_proposed_infection_date <- function(x){
    overall_result <- 1
    for (i in 1:length(assay_interpreters)){
      c_result <- assay_interpreters[[i]](x)
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
#' @examples
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
  stopifnot(class(range_start) == 'numeric')
  stopifnot(class(range_end) == 'numeric')
  stopifnot(length(unique(ihist$ptid)) == 1)

  if (FALSE){  # Debugging stuff
    ihist <- data.frame(
      ptid = c('p0', 'p0'),
      sample_date = c(as.numeric(as.Date('2016-03-01')), as.numeric(as.Date('2016-09-01'))),
      test = c('step_unit_testing', 'step_unit_testing'),
      result = c('-', '+'),
      stringsAsFactors = FALSE
    )
    range_start <- as.numeric(as.Date('2016-01-01'))
    range_end <- as.numeric(as.Date('2016-12-01'))
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
    assay_dynamics <- all_assay_dynamics[[ihist[i, 'test']]]
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
    test_details <- paste0(assay_dynamics$short_assayname, '\n', 
                           as.character(as.Date(floor(ihist[i, 'sample_date']), origin = '1970-01-01')), '\n', 
                           ihist[i, 'result'])
    all_xy_points <- rbind(
      all_xy_points,
      data.frame(ptid = ihist[i, 'ptid'],
                 sample_date = c_xy_points[['x']],
                 test_details = test_details,
                 prob_val = c_xy_points[['y']],
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
    if (verbose){cat('.\n')}
  }
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


