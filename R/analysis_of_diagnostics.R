# Note that this predicts the time since a rna assay with a 50% chance of testing positive at a viral load of 1 and not the true infection date.

# Each diagnostic test result provides information about an interval in which the patient could have been infected.
# A negative result provides a short (recent) interval in which the patient may have been infected (eclipse phase) and a much longer interval during which the patient could not have been infected.
# A positive result provides a short (recent) interval in which the patint may not have been infected (eclipse phase) and a much longer interval during which the patient could have been infected.

# TODO:
# Package the plotting into a single function
# Code to process the input data
# Code to produce a document with the results for each patient

library(ggplot2)
library(lubridate)
library(tidyr)
library(dplyr)
library(purrr)

# Assay details
totalnucleicacid_dynamics <- matrix(data = c(0, 0, #{{{
                                             1, 0,
                                             2, 0.01,
                                             3, 0.05,
                                             4, 0.275,
                                             5, 0.5,
                                             6, 0.725,
                                             7, 0.95,
                                             8, 0.99,
                                             9, 1,
                                             10, 1,
                                             11, 1,
                                             12, 1,
                                             13, 1,
                                             14, 1,
                                             15, 1,
                                             16, 1,
                                             17, 1,
                                             18, 1,
                                             19, 1,
                                             20, 1,
                                             21, 1,
                                             22, 1,
                                             23, 1,
                                             24, 1,
                                             25, 1,
                                             26, 1,
                                             27, 1,
                                             28, 1,
                                             29, 1,
                                             30, 1),
                                             ncol = 2,
                                             byrow = TRUE) #}}}

rnapcr_dynamics <- matrix(data = c(0, 0, #{{{
                                   1, 0,
                                   2, 0,
                                   3, 0,
                                   4, 0,
                                   5, 0,
                                   6, 0.01,
                                   7, 0.05,
                                   8, 0.1785714,
                                   9, 0.3071429,
                                   10, 0.4357143,
                                   11, 0.5642857,
                                   12, 0.6928571,
                                   13, 0.8214286,
                                   14, 0.95,
                                   15, 0.99,
                                   16, 1,
                                   17, 1,
                                   18, 1,
                                   19, 1,
                                   20, 1,
                                   21, 1,
                                   22, 1,
                                   23, 1,
                                   24, 1,
                                   25, 1,
                                   26, 1,
                                   27, 1,
                                   28, 1,
                                   29, 1,
                                   30, 1),
                                   ncol = 2,
                                   byrow = TRUE) #}}}

elisa_dynamics <- matrix(data = c(0, 0, #{{{
                                   1, 0,
                                   2, 0,
                                   3, 0,
                                   4, 0,
                                   5, 0,
                                   6, 0,
                                   7, 0,
                                   8, 0,
                                   9, 0,
                                   10, 0,
                                   11, 0,
                                   12, 0.01,
                                   13, 0.05,
                                   14, 0.1318182,
                                   15, 0.2136364,
                                   16, 0.2954545,
                                   17, 0.3772727,
                                   18, 0.4590909,
                                   19, 0.5409091,
                                   20, 0.6227273,
                                   21, 0.7045455,
                                   22, 0.7863636,
                                   23, 0.8681818,
                                   24, 0.95,
                                   25, 0.99,
                                   26, 1,
                                   27, 1,
                                   28, 1,
                                   29, 1,
                                   30, 1),
                                   ncol = 2,
                                   byrow = TRUE) #}}}

westernblot_dynamics <- matrix(data = c(0, 0, #{{{
                                        1, 0,
                                        2, 0,
                                        3, 0,
                                        4, 0,
                                        5, 0,
                                        6, 0,
                                        7, 0,
                                        8, 0,
                                        9, 0,
                                        10, 0,
                                        11, 0.01,
                                        12, 0.05,
                                        13, 0.1318182,
                                        14, 0.2136364,
                                        15, 0.2954545,
                                        16, 0.3772727,
                                        17, 0.4590909,
                                        18, 0.5409091,
                                        19, 0.6227273,
                                        20, 0.7045455,
                                        21, 0.7863636,
                                        22, 0.8681818,
                                        23, 0.95,
                                        24, 0.99,
                                        25, 1,
                                        26, 1,
                                        27, 1,
                                        28, 1,
                                        29, 1,
                                        30, 1),
                                        ncol = 2,
                                        byrow = TRUE) #}}}

geenius_dynamics <- matrix(data = c(0, 0, #{{{
                                    1, 0,
                                    2, 0,
                                    3, 0,
                                    4, 0,
                                    5, 0,
                                    6, 0,
                                    7, 0,
                                    8, 0,
                                    9, 0,
                                    10, 0.025,
                                    11, 0.050,
                                    12, 0.075,
                                    13, 0.100,
                                    14, 0.125,
                                    15, 0.150,
                                    16, 0.175,
                                    17, 0.200,
                                    18, 0.225,
                                    19, 0.250,
                                    20, 0.275,
                                    21, 0.300,
                                    22, 0.325,
                                    23, 0.350,
                                    24, 0.375,
                                    25, 0.400,
                                    26, 0.425,
                                    27, 0.450,
                                    28, 0.475,
                                    29, 0.500,
                                    30, 0.525,
                                    31, 0.550,
                                    32, 0.575,
                                    33, 0.600,
                                    34, 0.625,
                                    35, 0.650,
                                    36, 0.675,
                                    37, 0.700,
                                    38, 0.725,
                                    39, 0.750,
                                    40, 0.775,
                                    41, 0.800,
                                    42, 0.825,
                                    43, 0.850,
                                    44, 0.875,
                                    45, 0.900,
                                    46, 0.925,
                                    47, 0.950,
                                    48, 0.975,
                                    49, 1,
                                    50, 1,
                                    51, 1,
                                    52, 1,
                                    53, 1,
                                    54, 1,
                                    55, 1,
                                    56, 1,
                                    57, 1,
                                    58, 1,
                                    59, 1,
                                    60, 1),
                                    ncol = 2,
                                    byrow = TRUE) #}}}

geenius_dynamics_old <- matrix(data = c(0, 0, #{{{
                                        1, 0,
                                        2, 0,
                                        3, 0,
                                        4, 0,
                                        5, 0,
                                        6, 0,
                                        7, 0,
                                        8, 0,
                                        9, 0,
                                        10, 0,
                                        11, 0,
                                        12, 0,
                                        13, 0.01,
                                        14, 0.05,
                                        15, 0.1318182,
                                        16, 0.2136364,
                                        17, 0.2954545,
                                        18, 0.3772727,
                                        19, 0.4590909,
                                        20, 0.5409091,
                                        21, 0.6227273,
                                        22, 0.7045455,
                                        23, 0.7863636,
                                        24, 0.8681818,
                                        25, 0.95,
                                        26, 0.99,
                                        27, 1,
                                        28, 1,
                                        29, 1,
                                        30, 1),
                                        ncol = 2,
                                        byrow = TRUE) #}}}

all_assay_dynamics <- rbind( #{{{
  data.frame(assay = 'totalnucleicacid',
             days = totalnucleicacid_dynamics[,1],
             prob = totalnucleicacid_dynamics[,2]),
  data.frame(assay = 'rnapcr',
             days = rnapcr_dynamics[,1],
             prob = rnapcr_dynamics[,2]),
  data.frame(assay = 'elisa',
             days = elisa_dynamics[,1],
             prob = elisa_dynamics[,2]),
  data.frame(assay = 'westernblot',
             days = westernblot_dynamics[,1],
             prob = westernblot_dynamics[,2]),
  data.frame(assay = 'geenius',
             days = geenius_dynamics[,1],
             prob = geenius_dynamics[,2])) #}}}

# Processing a patient's results
prob_test_result_if_initial_infection_on_given_day <- function(assay_dynamics, result, days_to_visit){ #{{{
  if (days_to_visit > max(assay_dynamics[,1])){
    if (result == "+"){
      return(1)
    } else if (result == "-"){
      return(0)
    }
  } else if (days_to_visit < 0) {
    if (result == "+"){
      return(0)
    } else if (result == "-"){
      return(1)
    }
  } else {
    if (result == "+"){
      return(assay_dynamics[assay_dynamics[,1] == days_to_visit, 2])
    } else if (result == "-"){
      return(1-assay_dynamics[assay_dynamics[,1] == days_to_visit, 2])
    } else {
      stop("not implemented")
    }
  }
} #}}}

interpret_result <- function(assay_dynamics, result, range_start, range_end, visit_date, 
                             min_prob = 0.00, max_prob = 1.00){ #{{{
  days_vec <- NULL
  prob_vec <- NULL
  d2v_vec <- NULL
  for (i in ymd(range_start):ymd(range_end)){
    time_to_visit = as.duration(ymd(visit_date) - as_date(i))
    days_to_visit = as.numeric(time_to_visit)/(60*60*24)
    stopifnot(days_to_visit == floor(days_to_visit))

    days_vec <- c(days_vec, as_date(i))
    prob_vec <- c(prob_vec, 
                  prob_test_result_if_initial_infection_on_given_day(assay_dynamics = assay_dynamics, 
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
} #}}}

interpret_result_set <- function(result, range_start = NULL, range_end = NULL){ #{{{
  if (is.null(range_start)){
    range_start <- min(ymd(result$visit_date)) %m-% months(1)
  }

  if (is.null(range_end)){
    range_end <- max(ymd(result$visit_date)) %m+% months(1)
  }

  dat <- data.frame(date = as_date(ymd(range_start):ymd(range_end)))
  for (r_indx in 1:nrow(result)){
    c_result <- result[r_indx,]
    assay_dynamics <- get(paste(c_result$assay, '_dynamics', sep = ''))
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
} #}}}

# Parse and prepare data
parse_data <- function(file_name){ #{{{
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
} #}}}

remove_non_informative_results <- function(in_dat){ #{{{
# for each patient, for each test
# remove all entries where:
# the result is negative and the next result is also negative
# the result is positive and the previous result is also positive
  trimmed_dat <- in_dat[0,]
  c_ptid <- unique(in_dat$ptid)[5]
  for (c_ptid in unique(in_dat$ptid)){
    c_dat_pt <- subset(in_dat, ptid == c_ptid)
    c_assay <- 'rnapcr'
    c_assay <- unique(c_dat_pt$assay)[1]
    for (c_assay in unique(c_dat_pt$assay)){
      c_dat_as <- subset(c_dat_pt, assay == c_assay)
      if (nrow(c_dat_as) == 1){
        trimmed_dat <- rbind(trimmed_dat, c_dat_as)
      } else {
        c_dat_as <- c_dat_as[order(c_dat_as$visit_date),]
        if (length(unique(c_dat_as$result)) == 1){
          if (unique(c_dat_as$result) == '-'){
            # If only negative results present, only the last one is informative
            trimmed_dat <- rbind(trimmed_dat,
                                 c_dat_as[c_dat_as$visit_dat == max(c_dat_as$visit_dat),])
          } else {
            # If only positive results present, only the first one is informative
            trimmed_dat <- rbind(trimmed_dat,
                                 c_dat_as[c_dat_as$visit_dat == min(c_dat_as$visit_dat),])
          }
        } else {
          for (indx in 1:(nrow(c_dat_as)-1)){
            if (c_dat_as$result[indx] == '-' & c_dat_as$result[indx+1] != '-'){
              # a negative test immediately followed by another negative test is not informative
              trimmed_dat <- rbind(trimmed_dat, c_dat_as[indx,])
            }
          }
          for (indx in nrow(c_dat_as):2){
            if (c_dat_as$result[indx] == '+' & c_dat_as$result[indx-1] != '+'){
              # a positive test immediately preceeded by another positive test is not informative
              trimmed_dat <- rbind(trimmed_dat, c_dat_as[indx,])
            }
          }
        }
      }
    }
  }
  return(trimmed_dat)
} #}}}

make_vlines_dat <- function(in_dat){ #{{{
  vlines <- in_dat[,c('assay', 'visit_date', 'result')]
  vlines %>% map_if(is.factor, as.character) %>% as_tibble -> vlines
  vlines$facet_lab <- paste(vlines$assay, '\n', gsub('-', '', vlines$visit_date), '\n', vlines$result, sep = '')
  vlines$visit_date <- as_date(vlines$visit_date)
  return(vlines)
} #}}}

make_lrs_dat <- function(in_dat){ #{{{
  rs <- interpret_result_set(in_dat)
  lrs <- gather(rs, key = 'test', value = 'prob', -date)
  
  all_tests <- lrs %>% 
    group_by(date) %>%
    summarize(test = 'Aggregate',
              prob = prod(prob))
  lrs <- rbind(lrs, all_tests)
  lrs$facet_lab <- gsub('_', '\n', lrs$test)
  return(lrs)
} #}}}

# Plotting
patient_plot <- function(lrs, vlines){ #{{{
  x <- ggplot(lrs, aes(x = date, y = prob, col = facet_lab)) + 
    facet_grid(rows = vars(facet_lab)) + 
    geom_line() +
    theme(legend.position = 'none') +
    labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis',
         x = 'Date of intial infection') +
    geom_vline(data = vlines, aes(xintercept = visit_date), col = 'black') +
    scale_x_date(date_breaks = "months", date_labels = "%b-%y")
  return(x)
} #}}}

wrapped_patient_plot <- function(dat){ #{{{
  vlines <- make_vlines_dat(dat)
  lrs <- make_lrs_dat(dat)
  x <- patient_plot(lrs, vlines)
  return(x)
} #}}}

if (FALSE){
  # scribbles to make picture for DnE 2019
  str(ct_dat)
  x_dat <- ct_dat
  x_dat <- x_dat[!(x_dat$assay == 'elisa' & x_dat$result == '+'),]
  x_dat <- x_dat[!(x_dat$assay == 'rnapcr' & x_dat$result == '-'),]
  wrapped_patient_plot(x_dat)


  print('hello')


  # end of DnE 2019 presentation scribbles


  file_name <- '/home/phillipl/projects/hiv-founder-id/data/mock_amp_infection_timing_db/mock_list27Aug2018.csv'
  in_dat <- parse_data(file_name)
  c_pat <- in_dat$ptid[10]
  vlines <- make_vlines_dat(subset(in_dat, ptid == c_pat))
  lrs <- make_lrs_dat(subset(in_dat, ptid == c_pat))
  x <- patient_plot(lrs, vlines)
  print('hello')


  ###################################################################
    # old stuff
  ggplot(all_assay_dynamics, aes(col = assay, x = days, y = prob)) +
    geom_line() +
    labs(y = "Probability of testing positive",
         x = "Days since infection",
         col = "Assay")

  ggplot(all_assay_dynamics, aes(col = assay, x = days, y = 1-prob)) +
    geom_line() +
    labs(y = "Probability of testing negative",
         x = "Days since infection",
         col = "Assay")

  assay_dynamics <- rnapcr_dynamics
  result <- "-"
  days_to_visit <- 9


  file_name <- '/home/phillipl/projects/hiv-founder-id/data/mock_amp_infection_timing_db/mock_list27Aug2018.csv'

  range_start = "2017-05-25"
  range_end = "2017-10-03"

  result <- list('r1' = list(assay = 'rnapcr',
                             visit_date = '2017-08-01',
                             result = '+'),
                 'r2' = list(assay = 'elisa',
                             visit_date = '2017-08-01',
                             result = '-'),
                 'r3' = list(assay = 'rnapcr',
                             visit_date = '2017-09-05',
                             result = '+'),
                 'r4' = list(assay = 'elisa',
                             visit_date = '2017-09-05',
                             result = '+'))

  result_df <-
  do.call(rbind, lapply(result, as.data.frame))

  rs <-
  interpret_result_set(result_df, range_start, range_end)

  vlines <-
  do.call("rbind", lapply(result, as.data.frame))
  vlines %>% map_if(is.factor, as.character) %>% as_tibble -> vlines
  vlines$facet_lab <- paste(vlines$assay, '\n', gsub('-', '', vlines$visit_date), '\n', vlines$result, sep = '')
  vlines$visit_date <- as_date(vlines$visit_date)

  lrs <-
  gather(rs, key = 'test', value = 'prob', -date)

  all_tests <-
  lrs %>% 
    group_by(date) %>%
    summarize(test = 'Aggregate',
              prob = prod(prob))
  lrs <- rbind(lrs, all_tests)
  lrs$facet_lab <- gsub('_', '\n', lrs$test)

  patient_plot(lrs, vlines)


  x <-
  interpret_result(assay_dynamics = rnapcr_dynamics, 
                   result = '+', 
                   range_start = "2017-05-25", 
                   range_end = "2017-10-03", 
                   visit_date = "2017-08-01")
  i_result <- x


  prob_test_result_if_initial_infection_on_given_day(assay_dynamics = rnapcr_dynamics,
                                                     result = '+',
                                                     days_to_visit = 10)

  #prob_test_result_if_initial_infection_on_given_day <- function(assay_dynamics, result, days_to_visit){

  interpret_result_set_deprecated <- function(result, range_start, range_end){ #{{{
    dat <- data.frame(date = as_date(ymd(range_start):ymd(range_end)))
    for (r_indx in 1:length(result)){
      c_result <- result[[r_indx]]
      assay_dynamics <- get(paste(c_result[['assay']], '_dynamics', sep = ''))
      i_result <- interpret_result(assay_dynamics = assay_dynamics,
                                   result = c_result$result,
                                   range_start = range_start,
                                   range_end = range_end,
                                   visit_date = c_result$visit_date)
      p_result <- i_result[,c('days', 'prob')]
      names(p_result) <- c('date', paste(c_result[['assay']], '_', gsub("-", "", c_result$visit_date), '_', c_result$result, sep = '') )
      dat <- merge(dat, p_result)
    }
    return(dat)
  } #}}}
}
