# Note that this predicts the time since a rna assay with a 50% chance of testing positive at a viral load of 1 and not the true infection date.


# Each diagnostic test result provides information about an interval in which the patient could have been infected.
# A negative result provides a short (recent) interval in which the patient may have been infected (eclipse phase) and a much longer interval during which the patient could not have been infected.
# A positive result provides a short (recent) interval in which the patint may not have been infected (eclipse phase) and a much longer interval during which the patient could have been infected.


# TODO:
# Package the plotting into a single function
# Code to process the input data
# Code to produce a document with the results for each patient






##############
# Start here
##############

# TODO: delete this
if (FALSE){
  # debugging / dev scribbles
  assay_dynamics <- list(fun = linear_assay_dynamics,
                         params = list(diagnostic_delay = 10))
  days_since_ddi_1 <- 9
  days_since_ddi_1 <- 10
  days_since_ddi_1 <- 3
  result <- '+'
  result <- '-'

  range_start <- '2017-02-01'
  range_end <- '2017-10-31'
  visit_date <- '2017-04-01'
  min_prob <- 0
  max_prob <- 1

  assay <- 'Abbott Architect HIV Ag/Ab Combo'

  result_set <- structure(list(ptid = c(102967486L, 102967486L, 102967486L, 102967486L,
    102967486L, 102967486L, 102967486L, 102967486L, 102967486L, 102967486L,
    102967486L, 102967486L, 102967486L, 102967486L), visit_date = c("2017-05-02",
    "2017-05-02", "2017-05-30", "2017-05-30", "2017-06-30", "2017-06-30",
    "2017-06-30", "2017-07-11", "2017-07-11", "2017-07-11", "2017-07-28",
    "2017-08-10", "2017-08-27", "2017-09-24"), assay = c("elisa",
    "rnapcr", "elisa", "rnapcr", "elisa", "geenius", "rnapcr", "elisa",
    "geenius", "rnapcr", "rnapcr", "rnapcr", "rnapcr", "rnapcr"),
        o_result = c(2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
        1L, 1L, 1L), result = c("-", "-", "-", "+", "+", "+", "+",
        "+", "+", "+", "+", "+", "+", "+")), row.names = c(NA, 14L
    ), class = "data.frame")


}

#' Probabilities of certain test results
#'
#' Turns assay dynamics into the probability of observing a specified test result a given number of days since DDI_1
#'
#' @param assay_dynamics A list providing the assay_dynamics, (TODO: elaborate)
#' @param result '+' or '-' indicating the test result.
#' @param days_to_visit The number of days between DDI_1 and the test result.
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
#' @param range_start The date furthest into the past that should be considered as a plausible day for DDI_1.
#' @param range_end The most recent date that should be considered as a plausible day for DDI_1.
# @export

interpret_result_set <- function(result_set, range_start = NULL, range_end = NULL){
  if (is.null(range_start)){
    range_start <- min(ymd(result$visit_date)) %m-% months(1)
  }

  if (is.null(range_end)){
    range_end <- max(ymd(result$visit_date)) %m+% months(1)
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

#' Parse and prepare the data for the first mock dataset
#'
#' This script will load the input file and prepare it for processing. It is currently specifically designed to be used with the mock datasets. A more general version of this function must be written.
#'
#' @param file_name Full path the the csv file
#' @export

parse_data_first_mock <- function(file_name){
  dat <- read.csv(file_name, stringsAsFactors = FALSE)
  ldat <- gather(dat, key = 'test', value = 'result', elisa, geenius, rnapcr, totalnucleicacid, westernblot)
  ldat <- ldat %>% arrange(ptid, drawdt, test) %>% filter(!is.na(result))
  ldat$p_result <- ''
  
  # 1: Positive = 1
  # 8,9: positive (HIV2)
  # 10,11: positive but not quantifiable
  # 15: HIV untypable
  # 18: cross reactive hiv1 and 2
  # 19,20,21 reactive antibody / antigen / both
  ldat$p_result[ldat$result %in% c(1, 8, 9, 10, 11, 19, 20, 21)] <- '+'
  
  # Negative = 2
  ldat$p_result[ldat$result == 2] <- '-'
  
  # Indeterminate
  ldat <- ldat %>% filter(!(result %in% c(3,4,5,6,7, 12, 13, 14, 16, 17, 22)))

  names(ldat)[names(ldat) == 'drawdt'] <- 'visit_date'
  names(ldat)[names(ldat) == 'test'] <- 'assay'
  names(ldat)[names(ldat) == 'result'] <- 'o_result'
  names(ldat)[names(ldat) == 'p_result'] <- 'result'
  return(ldat)
}

#remove_non_informative_results <- function(in_dat){ #{{{
## for each patient, for each test
## remove all entries where:
## the result is negative and the next result is also negative
## the result is positive and the previous result is also positive
#  trimmed_dat <- in_dat[0,]
#  c_ptid <- unique(in_dat$ptid)[5]
#  for (c_ptid in unique(in_dat$ptid)){
#    c_dat_pt <- subset(in_dat, ptid == c_ptid)
#    c_assay <- 'rnapcr'
#    c_assay <- unique(c_dat_pt$assay)[1]
#    for (c_assay in unique(c_dat_pt$assay)){
#      c_dat_as <- subset(c_dat_pt, assay == c_assay)
#      if (nrow(c_dat_as) == 1){
#        trimmed_dat <- rbind(trimmed_dat, c_dat_as)
#      } else {
#        c_dat_as <- c_dat_as[order(c_dat_as$visit_date),]
#        if (length(unique(c_dat_as$result)) == 1){
#          if (unique(c_dat_as$result) == '-'){
#            # If only negative results present, only the last one is informative
#            trimmed_dat <- rbind(trimmed_dat,
#                                 c_dat_as[c_dat_as$visit_dat == max(c_dat_as$visit_dat),])
#          } else {
#            # If only positive results present, only the first one is informative
#            trimmed_dat <- rbind(trimmed_dat,
#                                 c_dat_as[c_dat_as$visit_dat == min(c_dat_as$visit_dat),])
#          }
#        } else {
#          for (indx in 1:(nrow(c_dat_as)-1)){
#            if (c_dat_as$result[indx] == '-' & c_dat_as$result[indx+1] != '-'){
#              # a negative test immediately followed by another negative test is not informative
#              trimmed_dat <- rbind(trimmed_dat, c_dat_as[indx,])
#            }
#          }
#          for (indx in nrow(c_dat_as):2){
#            if (c_dat_as$result[indx] == '+' & c_dat_as$result[indx-1] != '+'){
#              # a positive test immediately preceeded by another positive test is not informative
#              trimmed_dat <- rbind(trimmed_dat, c_dat_as[indx,])
#            }
#          }
#        }
#      }
#    }
#  }
#  return(trimmed_dat)
#} #}}}
#
#make_vlines_dat <- function(in_dat){ #{{{
#  vlines <- in_dat[,c('assay', 'visit_date', 'result')]
#  vlines %>% map_if(is.factor, as.character) %>% as_tibble -> vlines
#  vlines$facet_lab <- paste(vlines$assay, '\n', gsub('-', '', vlines$visit_date), '\n', vlines$result, sep = '')
#  vlines$visit_date <- as_date(vlines$visit_date)
#  return(vlines)
#} #}}}
#
#make_lrs_dat <- function(in_dat){ #{{{
#  rs <- interpret_result_set(in_dat)
#  lrs <- gather(rs, key = 'test', value = 'prob', -date)
#  
#  all_tests <- lrs %>% 
#    group_by(date) %>%
#    summarize(test = 'Aggregate',
#              prob = prod(prob))
#  lrs <- rbind(lrs, all_tests)
#  lrs$facet_lab <- gsub('_', '\n', lrs$test)
#  return(lrs)
#} #}}}
#
## Plotting
#patient_plot <- function(lrs, vlines){ #{{{
#  x <- ggplot(lrs, aes(x = date, y = prob, col = facet_lab)) + 
#    facet_grid(rows = vars(facet_lab)) + 
#    geom_line() +
#    theme(legend.position = 'none') +
#    labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis',
#         x = 'Date of intial infection') +
#    geom_vline(data = vlines, aes(xintercept = visit_date), col = 'black') +
#    scale_x_date(date_breaks = "months", date_labels = "%b-%y")
#  return(x)
#} #}}}
#
#wrapped_patient_plot <- function(dat){ #{{{
#  vlines <- make_vlines_dat(dat)
#  lrs <- make_lrs_dat(dat)
#  x <- patient_plot(lrs, vlines)
#  return(x)
#} #}}}
#
#if (FALSE){
#  # scribbles to make picture for DnE 2019
#  str(ct_dat)
#  x_dat <- ct_dat
#  x_dat <- x_dat[!(x_dat$assay == 'elisa' & x_dat$result == '+'),]
#  x_dat <- x_dat[!(x_dat$assay == 'rnapcr' & x_dat$result == '-'),]
#  wrapped_patient_plot(x_dat)
#
#
#  print('hello')
#
#
#  # end of DnE 2019 presentation scribbles
#
#
#  file_name <- '/home/phillipl/projects/hiv-founder-id/data/mock_amp_infection_timing_db/mock_list27Aug2018.csv'
#  in_dat <- parse_data(file_name)
#  c_pat <- in_dat$ptid[10]
#  vlines <- make_vlines_dat(subset(in_dat, ptid == c_pat))
#  lrs <- make_lrs_dat(subset(in_dat, ptid == c_pat))
#  x <- patient_plot(lrs, vlines)
#  print('hello')
#
#
#  ###################################################################
#    # old stuff
#  ggplot(all_assay_dynamics, aes(col = assay, x = days, y = prob)) +
#    geom_line() +
#    labs(y = "Probability of testing positive",
#         x = "Days since infection",
#         col = "Assay")
#
#  ggplot(all_assay_dynamics, aes(col = assay, x = days, y = 1-prob)) +
#    geom_line() +
#    labs(y = "Probability of testing negative",
#         x = "Days since infection",
#         col = "Assay")
#
#  assay_dynamics <- rnapcr_dynamics
#  result <- "-"
#  days_to_visit <- 9
#
#
#  file_name <- '/home/phillipl/projects/hiv-founder-id/data/mock_amp_infection_timing_db/mock_list27Aug2018.csv'
#
#  range_start = "2017-05-25"
#  range_end = "2017-10-03"
#
#  result <- list('r1' = list(assay = 'rnapcr',
#                             visit_date = '2017-08-01',
#                             result = '+'),
#                 'r2' = list(assay = 'elisa',
#                             visit_date = '2017-08-01',
#                             result = '-'),
#                 'r3' = list(assay = 'rnapcr',
#                             visit_date = '2017-09-05',
#                             result = '+'),
#                 'r4' = list(assay = 'elisa',
#                             visit_date = '2017-09-05',
#                             result = '+'))
#
#  result_df <-
#  do.call(rbind, lapply(result, as.data.frame))
#
#  rs <-
#  interpret_result_set(result_df, range_start, range_end)
#
#  vlines <-
#  do.call("rbind", lapply(result, as.data.frame))
#  vlines %>% map_if(is.factor, as.character) %>% as_tibble -> vlines
#  vlines$facet_lab <- paste(vlines$assay, '\n', gsub('-', '', vlines$visit_date), '\n', vlines$result, sep = '')
#  vlines$visit_date <- as_date(vlines$visit_date)
#
#  lrs <-
#  gather(rs, key = 'test', value = 'prob', -date)
#
#  all_tests <-
#  lrs %>% 
#    group_by(date) %>%
#    summarize(test = 'Aggregate',
#              prob = prod(prob))
#  lrs <- rbind(lrs, all_tests)
#  lrs$facet_lab <- gsub('_', '\n', lrs$test)
#
#  patient_plot(lrs, vlines)
#
#
#  x <-
#  interpret_result(assay_dynamics = rnapcr_dynamics, 
#                   result = '+', 
#                   range_start = "2017-05-25", 
#                   range_end = "2017-10-03", 
#                   visit_date = "2017-08-01")
#  i_result <- x
#
#
#  prob_test_result_if_initial_infection_on_given_day(assay_dynamics = rnapcr_dynamics,
#                                                     result = '+',
#                                                     days_to_visit = 10)
#
#  #prob_test_result_if_initial_infection_on_given_day <- function(assay_dynamics, result, days_to_visit){
#
#  interpret_result_set_deprecated <- function(result, range_start, range_end){ #{{{
#    dat <- data.frame(date = as_date(ymd(range_start):ymd(range_end)))
#    for (r_indx in 1:length(result)){
#      c_result <- result[[r_indx]]
#      assay_dynamics <- get(paste(c_result[['assay']], '_dynamics', sep = ''))
#      i_result <- interpret_result(assay_dynamics = assay_dynamics,
#                                   result = c_result$result,
#                                   range_start = range_start,
#                                   range_end = range_end,
#                                   visit_date = c_result$visit_date)
#      p_result <- i_result[,c('days', 'prob')]
#      names(p_result) <- c('date', paste(c_result[['assay']], '_', gsub("-", "", c_result$visit_date), '_', c_result$result, sep = '') )
#      dat <- merge(dat, p_result)
#    }
#    return(dat)
#  } #}}}
#}
