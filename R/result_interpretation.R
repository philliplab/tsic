#' Probabilities of certain test results
#'
#' Turns assay dynamics into the probability of observing a specified test result a given number of days since XXX
#'
#' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
#' @param result '+' or '-' indicating the test result.
#' @param days_to_visit The number of days between XXX and the test result.
#' @export

prob_test_result_given_days_to_visit <- function(assay_dynamics, result, days_to_visit){
  assay_dynamics$params$x <- days_to_visit
  x <- do.call(assay_dynamics$fun, assay_dynamics$params)

  if (result == "+"){
    return(x)
  } else if (result == "-"){
    return(1-x)
  } else {
    stop('result must be "+" or "-"')
  }
}

#' Probability of DDI_1 curve implied by test result
#'
#' Compute a curve of probabilities of observing a specified test result on a specified date assuming the day of DDI_1. Thus assume individually that each day in a range of dates is the DDI_1 and then compute the probability of observing the specified result.
#'
#' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
#' @param result '+' or '-' indicating the test result.
#' @param range_start The date furthest into the past that should be considered as a plausible day for DDI_1.
#' @param range_end The most recent date that should be considered as a plausible day for DDI_1.
#' @param visit_date The date of the visit on whose sample the test result was observed
#' @param min_prob The most 'certain' you are that the given test results excludes dates far from the actual visit. Think of this as a chance that the test may give an incorrect result.
#' @param max_prob The most 'certain' you are that the given test results excludes dates far from the actual visit. Think of this as a chance that the test may give an incorrect result. TODO: Think more about what this actually means.
#' @export

interpret_result <- function(assay_dynamics, result, range_start, range_end, visit_date, 
                             min_prob = 0.00, max_prob = 1.00){
  days_vec <- NULL
  prob_vec <- NULL
  d2v_vec <- NULL
  for (i in ymd(range_start):ymd(range_end)){
    time_to_visit = as.duration(ymd(visit_date) - as_date(i))
    days_to_visit = as.numeric(time_to_visit)/(60*60*24)
    stopifnot(days_to_visit == floor(days_to_visit))

    days_vec <- c(days_vec, as_date(i))
    prob_vec <- c(prob_vec, 
                  prob_test_result_given_days_to_visit(assay_dynamics = assay_dynamics, 
                                                       result = result, 
                                                       days_to_visit = days_to_visit)
                  )
    d2v_vec <- c(d2v_vec, days_to_visit)
  }
  prob_vec[prob_vec < min_prob] <- min_prob
  prob_vec[prob_vec > max_prob] <- max_prob
  return(data.frame(days = as_date(days_vec),
                    prob = prob_vec,
                    d2v = d2v_vec)
        )
}

#' Given a set of results, construct DDI_1 probability curves
#'
#' Calls interpret_result on a set of test results and aggregate the results into a data.frame
#'
#' @param result_set A data.frame of all the test results and dates for a patient
#' @param range_start The date furthest into the past that should be considered as a plausible day for DDI_1. Can be NULL, then defaults to 2 months before the earliest result.
#' @param range_end The most recent date that should be considered as a plausible day for DDI_1. Can be NULL, then defaults to 1 week after the most recent result.
#' @export

interpret_result_set <- function(result_set, range_start = NULL, range_end = NULL){
  if (is.null(range_start)){
    range_start <- min(ymd(result_set$visit_date)) %m-% months(2)
  }
  if (is.null(range_end)){
    range_end <- max(ymd(result_set$visit_date)) %m+% weeks(1)
  }

  dat <- data.frame(date = as_date(ymd(range_start):ymd(range_end)))
  for (r_indx in 1:nrow(result_set)){
    c_result <- result_set[r_indx,]
    assay_dynamics <- get_assay_dynamics(assay = c_result$assay)
    i_result <- interpret_result(assay_dynamics = assay_dynamics,
                                 result = c_result$result,
                                 range_start = range_start,
                                 range_end = range_end,
                                 visit_date = c_result$visit_date)
    p_result <- i_result[,c('days', 'prob')]
    tmpx <- paste(c_result$assay, '_', 
                  gsub("-", "", c_result$visit_date), 
                  '_', c_result$result, sep = '')
    names(p_result) <- c( 'date', tmpx )
    dat <- merge(dat, p_result)
  }
  return(dat)
}

