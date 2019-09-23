#' Estimates the 2.5, 50 and 97.5 percentiles
#'
#' Various hacks required to work around the inadequacy of R's numerical analysis tools.
#'
#' @param fun The function whose percentiles are required
#' @param range_start Start of interval containing the percentiles
#' @param range_end End of interval containing the percentiles
#' @export

estimate_lb_med_ub <- function(fun, range_start, range_end){
  if (FALSE){
    range_start <- -100
    range_end <- 100
    fun <- dexp
    fun <- dnorm

    tiles <- c(0.025, 0.5, 0.975)
    qnorm(tiles)
    qexp(tiles)

  }
  ranges <- trim_range(fun, range_start, range_end, tol = 0.1^50)
  range_start <- ranges$range_start
  range_end <- ranges$range_end
  xy_points <- get_scatterpoints(fun, range_start:range_end, max_delta = 0.001, min_length = 0.001)
  xy_points <- reduce_x_points(xy_points$x, xy_points$y)

  range_start <- min(xy_points$x)
  range_end <- max(xy_points$x)
  total_aoc <- pracma::integral(fun = fun,
                                xmin = range_start,
                                xmax = range_end,
                                no_intervals = 1000)

  integrated_fun <- function(x){
    pracma::integral(fun = function(x){fun(x)/total_aoc},
                     xmin = range_start, 
                     xmax = min(x, range_end),
                     no_intervals = 1000)
  }

  midpoint_heights <- (xy_points$y[1:(length(xy_points$y)-1)] + xy_points$y[2:(length(xy_points$y))]) / 2
  int_lengths <- xy_points$x[2:(length(xy_points$x))] - xy_points$x[1:(length(xy_points$x)-1)]
  riemans <- midpoint_heights * int_lengths
  midpoints <- (xy_points$x[2:(length(xy_points$x))] + xy_points$x[1:(length(xy_points$x)-1)]) / 2
  probs <- cumsum(riemans)
  probs <- probs / max(probs)
  x <- midpoints
#  plot(probs ~ x)

  find_perc <- function(value, width_toggle){
#  value <- 0.025
#  width_toggle <- 0.01
    start_indx <- max(which(probs < value - width_toggle))
    end_indx <- min(which(probs > value + width_toggle))

    tight_range_start <- x[start_indx]
    tight_range_end <- x[end_indx]

    c_perc_w <- optimize(f = function(x){abs(integrated_fun(x) - value)},
                         interval = c(range_start, range_end),
                         tol = 2.220446e-16)
    c_perc_t <- optimize(f = function(x){abs(integrated_fun(x) - value)},
                         interval = c(tight_range_start, tight_range_end),
                         tol = 2.220446e-16)
    if (c_perc_w$objective < c_perc_t$objective){
      c_perc <- c_perc_w
    } else {
      c_perc <- c_perc_t
    }
    return(c_perc)
  }

  lb  <- find_perc(value = 0.025, width_toggle = 0.01)
  med <- find_perc(value = 0.500, width_toggle = 0.01)
  ub  <- find_perc(value = 0.975, width_toggle = 0.01)

  return(list(lb = lb$minimum,
              med = med$minimum,
              ub = ub$minimum))
}

estimate_lb_med_ub_failed <- function(fun, range_start, range_end){
  if (FALSE){
    fun <- function(x){
      return(ifelse(x < 0 | x > 10, 0, x))
    }
    range_start <- -5
    range_end <- 15
    x <- 4


    fun <- dnorm
    range_start <- -39
    range_end <- 39
  
    fun <- dexp
    range_start <- -10000
    range_end <- 10000
    range_start <- -1
    range_end <- 746
  }

  total_aoc <- pracma::integral(fun = fun,
                                xmin = range_start,
                                xmax = range_end,
                                no_intervals = 1000)

  integrated_fun <- function(x){
    pracma::integral(fun = function(x){fun(x)/total_aoc},
                     xmin = range_start, 
                     xmax = min(x, range_end),
                     no_intervals = 1000)
  }

  # such a horrible hack to make this work for the normal distribution
  n_precuts <- 20
  precuts <- seq(from = range_start, to = range_end, length.out = n_precuts)
  lb_start <- precuts[1]
  med_start <- precuts[1]
  ub_start <- precuts[1]

  lb_end <- precuts[n_precuts]
  med_end <- precuts[n_precuts]
  ub_end <- precuts[n_precuts]

  start_from_here <- precuts[1]
  for (i in 1:n_precuts){
    area_to_left <- integrated_fun(precuts[i])
    if (area_to_left < 0.02){
      lb_start <- precuts[i]
    }
    if (area_to_left < 0.49){
      med_start <- precuts[i]
    }
    if (area_to_left < 0.097){
      ub_start <- precuts[i]
    }
    
    if (area_to_left > 0.03 & lb_end == precuts[n_precuts]){
      lb_end <- precuts[i]
    }
    if (area_to_left > 0.51 & med_end == precuts[n_precuts]){
      med_end <- precuts[i]
    }
    if (area_to_left > 0.098 & ub_end == precuts[n_precuts]){
      ub_end <- precuts[i]
    }
  }

  lb1 <- optimize(f = function(x){abs(integrated_fun(x) - 0.025)},
           interval = c(range_start, range_end),
           tol = 2.220446e-16)
  lb2 <- optimize(f = function(x){abs(integrated_fun(x) - 0.025)},
           interval = c(lb_start, lb_end),
           tol = 2.220446e-16)
  if (lb1$objective < lb2$objective){
    lb <- lb1
  } else {
    lb <- lb2
  }

  ub1 <- optimize(f = function(x){abs(integrated_fun(x) - 0.975)},
           interval = c(range_start, range_end),
           tol = 2.220446e-16)
  ub2 <- optimize(f = function(x){abs(integrated_fun(x) - 0.975)},
           interval = c(ub_start, ub_end),
           tol = 2.220446e-16)
  if (ub1$objective < ub2$objective){
    ub <- ub1
  } else {
    ub <- ub2
  }

  med1 <- optimize(f = function(x){abs(integrated_fun(x) - 0.5)},
           interval = c(range_start, range_end),
           tol = 2.220446e-16)
  med2 <- optimize(f = function(x){abs(integrated_fun(x) - 0.5)},
           interval = c(med_start, med_end),
           tol = 2.220446e-16)
  if (med1$objective < med2$objective){
    med <- med1
  } else {
    med <- med2
  }

  return(list(lb = lb$minimum,
              med = med$minimum,
              ub = ub$minimum))

}

check_lb_med_ub <- function(lb, med, ub, fun, range_start, range_end){
  if (FALSE){

    lb  = qbeta(0.025, 2, 5)
    med = qbeta(0.5,   2, 5)
    ub  = qbeta(0.975, 2, 5)
    fun = function(x){dbeta(x, 2, 5)}
    range_start = 0
    range_end = 1



    lb <- lb$minimum
    med <- med$minimum
    ub <- ub$minimum
  }
  total_aoc <- pracma::integral(f = fun,
                         xmin = range_start,
                         xmax = range_end,
                         no_intervals = 1000)
  area_left_of_lb <- pracma::integral(f = function(x){fun(x)/total_aoc},
                               xmin = range_start,
                               xmax = lb,
                               no_intervals = 1000)
  area_left_of_med <- pracma::integral(f = function(x){fun(x)/total_aoc},
                               xmin = range_start,
                               xmax = med,
                               no_intervals = 1000)
  area_left_of_ub <- pracma::integral(f = function(x){fun(x)/total_aoc},
                               xmin = range_start,
                               xmax = ub,
                               no_intervals = 1000)
  return(list(area_left_of_lb = area_left_of_lb,
              area_left_of_med = area_left_of_med,
              area_left_of_ub = area_left_of_ub))
}

## obsolete
#
##' Converts an interpretated results set to long format
##'
##' Converts an interpretated result set as produced by interpret_result_set to long format
##'
##' @param interpretation_set An interpretated results set as produced by interpret_result_set
##' @export
#
#convert_interpretation_to_long <- function(interpretation_set){
##  rs <- interpret_result_set(in_dat)
#  rs <- interpretation_set
#  lrs <- gather(rs, key = 'test', value = 'prob', -date)
#  
#  all_tests <- lrs %>% 
#    group_by(date) %>%
#    summarize(test = 'Aggregate',
#              prob = prod(prob))
#  lrs <- rbind(lrs, all_tests)
#  lrs$facet_lab <- gsub('_', '\n', lrs$test)
#  return(lrs)
#}
#
##' Convert Aggregate prob curve into an ecdf.
##'
##' @param dat An interpreted resultset in long format
##' @export
#
#convert_aggregate_into_ecdf <- function(dat){
#  dat <- subset(dat, test == 'Aggregate')
#  dat <- dat[order(dat$date),]
#  if (length(unique(dat$prob)) == 1){
#    return(NULL)
#  } else {
#    y <- cumsum(dat$prob)
#    y <- y/max(y)
#    x <- dat$date
#    return(list(x = x, y = y))
#  }
#}
#
##' Interpolates - single point
##'
##' Crappy function needed because internal approx fails
##' @param x, y vectors giving the coordinates of the points to be interpolated.
##' @param xout Single value specifying where interpolation is to take place.
#
#manual_approx_one <- function(x, y, xout){
#  ord_vec <- order(x)
#  x <- x[ord_vec]
#  y <- y[ord_vec]
#  for (i in 1:length(x)){
#    if (x[i] >= xout){
#      x2 <- x[i]
#      x1 <- x[i-1]
#      break
#    }
#  }
#  if (x2 == xout){
#    return(y[i])
#  } else {
#    top_gap <- x2 - xout
#    bot_gap <- xout - x1
#    tot_gap <- x2 - x1
#    return((bot_gap/tot_gap)*y[i] + (top_gap/tot_gap)*y[i-1])
#  }
#}
#
##' Interpolates - multiple points
##'
##' Crappy function needed because internal approx fails
##' @param x, y vectors giving the coordinates of the points to be interpolated.
##' @param xout Vector of values specifying where interpolation is to take place.
#
#manual_approx <- function(x, y, xout){
#  result <- NULL
#  for (c_xout in xout){
#    result <- c(result, manual_approx_one(x, y, c_xout))
#  }
#  return(result)
#}
#
#
##' Estimate the LB, UB and median for DDI_1
##'
##' @param dat An interpreted resultset in long format
##' @export
#
#estimate_lb_med_ub <- function(dat){
#  aggre_ecdf <- convert_aggregate_into_ecdf(dat)
#  if (is.null(aggre_ecdf)){
#    return(NULL)
#  } else {
#    result <- manual_approx(x = aggre_ecdf$y, y = as.numeric(as_datetime(aggre_ecdf$x)), xout = c(0.025, 0.5, 0.975))
#    dates <- as_datetime(result)
##    result <- approx(x = aggre_ecdf$y, y = as_datetime(aggre_ecdf$x), xout = c(0.025, 0.5, 0.975))
##    dates <- as_datetime(result$y)
#    return(dates)
#  }
#}
#
##' Remove non-informative test results
##'
##' For each assay, if this result is negative and the next one is also negative, remove this result. For each assys, this this result is positive and the previous result is also positive, then remove this result. Obviously this needs refinement.
##'
##' @param in_dat An interpreted resultset in long format
##' @export
#
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
#}
#
