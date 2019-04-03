#' Converts an interpretated results set to long format
#'
#' Converts an interpretated result set as produced by interpret_result_set to long format
#'
#' @param interpretation_set An interpretated results set as produced by interpret_result_set
#' @export

convert_interpretation_to_long <- function(interpretation_set){
#  rs <- interpret_result_set(in_dat)
  rs <- interpretation_set
  lrs <- gather(rs, key = 'test', value = 'prob', -date)
  
  all_tests <- lrs %>% 
    group_by(date) %>%
    summarize(test = 'Aggregate',
              prob = prod(prob))
  lrs <- rbind(lrs, all_tests)
  lrs$facet_lab <- gsub('_', '\n', lrs$test)
  return(lrs)
}

#' Convert Aggregate prob curve into an ecdf.
#'
#' @param dat An interpreted resultset in long format
#' @export

convert_aggregate_into_ecdf <- function(dat){
  dat <- subset(dat, test == 'Aggregate')
  dat <- dat[order(dat$date),]
  if (length(unique(dat$prob)) == 1){
    return(NULL)
  } else {
    y <- cumsum(dat$prob)
    y <- y/max(y)
    x <- dat$date
    return(list(x = x, y = y))
  }
}

#' Estimate the LB, UB and median for DDI_1
#'
#' @param dat An interpreted resultset in long format
#' @export

estimate_lb_med_ub <- function(dat){
  aggre_ecdf <- convert_aggregate_into_ecdf(dat)
  if (is.null(aggre_ecdf)){
    return(NULL)
  } else {
    result <- approx(x = aggre_ecdf$y, y = as_datetime(aggre_ecdf$x), xout = c(0.025, 0.5, 0.975))
    dates <- as_datetime(result$y)
    return(dates)
  }
}

#' Remove non-informative test results
#'
#' For each assay, if this result is negative and the next one is also negative, remove this result. For each assys, this this result is positive and the previous result is also positive, then remove this result. Obviously this needs refinement.
#'
#' @param in_dat An interpreted resultset in long format
#' @export

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
}
