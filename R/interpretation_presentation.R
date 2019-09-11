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
  stopifnot(length(unique(ihist$ptid)) == 1)
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
                 test_details = paste0(ihist[i, 'test'], '_', ihist[i, 'sample_date'], '_', ihist[i, 'result']),
                 prob_val = c_xy_points[['y']],
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
    if (verbose){cat('.\n')}
  }
#construct_aggregate_interpreter <- function(ihist){
  aggregate_interpreter <- construct_aggregate_interpreter(ihist)
  aggregate_seedpoints <- sort(unique(all_xy_points[['sample_date']]))
  c_xy_points <- get_scatterpoints(fun = aggregate_interpreter,
                                   seedpoints = aggregate_seedpoints)
  all_xy_points <- rbind(
    all_xy_points,
    data.frame(ptid = ihist[1, 'ptid'],
               sample_date = c_xy_points[['x']],
               test_details = paste0('Aggregate'),
               prob_val = c_xy_points[['y']],
               stringsAsFactors = FALSE),
    stringsAsFactors = FALSE)
  return(all_xy_points)
}













## obsolete
#
##' Create Dataset of test dates
##'
##' Prepares a dataset with the dates of the tests for easy inclusion in plots.
##'
##' @param in_dat An interpretated datasets in long format
##' @export
#
#make_vlines_dat <- function(in_dat){
#  vlines <- unique(in_dat$test)
#  vlines <- vlines[!grepl("Aggregate", vlines)]
#  vlines <- strsplit(vlines, '_')
#  vlines <- lapply(vlines, function(x){data.frame(assay = x[1], visit_date = x[2], result = x[3], stringsAsFactors = FALSE)})
#  all_dat <- NULL
#  for (i in vlines){
#    all_dat <- rbind(all_dat, i)
#  }
#  vlines <- all_dat
#  vlines$facet_lab <- paste(vlines$assay, '\n', gsub('-', '', vlines$visit_date), '\n', vlines$result, sep = '')
#  vlines$visit_date <- as_date(vlines$visit_date)
#
#  lb_med_ub <- round_date(estimate_lb_med_ub(in_dat), unit = 'day')
#  if (!is.null(lb_med_ub)){
#    vlines <- rbind(vlines, 
#      data.frame(assay = "Aggregate",
#                 visit_date = lb_med_ub,
#                 result = NA,
#                 facet_lab = "Aggregate"))
#  }
#
#  return(vlines)
#}
#
##' Plots a patients timelines
##'
##' Plots dates and the likelihoods that these dates are the DDI_1 for a patient.
##'
##' @param lrs An interpretated result set in long format.
##' @param vlines Vertical lines indicating the test dates as produced by make_vlines_dat.
##' @export
#
#patient_plot <- function(lrs, vlines){
#  lrs$facet_lab <- factor(lrs$facet_lab, 
#         levels = c('Aggregate', setdiff(sort(unique(lrs$facet_lab)), 'Aggregate')),
#         ordered = TRUE)
#  vlines$facet_lab <- factor(vlines$facet_lab, 
#         levels = c('Aggregate', setdiff(sort(unique(lrs$facet_lab)), 'Aggregate')),
#         ordered = TRUE)
#  x <- ggplot2::ggplot(lrs, aes(x = date, y = prob, col = facet_lab)) + 
#    facet_grid(rows = vars(facet_lab)) + 
#    geom_line() +
#    theme(legend.position = 'none') +
#    labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis',
#         x = 'Date of intial infection') +
#    geom_vline(data = vlines, aes(xintercept = visit_date), col = 'black') +
##    geom_text(data = vlines, 
##              aes(x = visit_date, y = 0.05, label = result), 
##              col = 'black', size = 10, 
##              hjust = 0, nudge_x = 0.05,
##              vjust = 0) +
#    scale_x_date(date_breaks = "months", date_labels = "%b-%y")
#  return(x)
#}
#
##' Plots a patients timelines compactly
##'
##' Plots dates and the likelihoods that these dates are the DDI_1 for a patient.
##'
##' @param lrs An interpretated result set in long format.
##' @param vlines Vertical lines indicating the test dates as produced by make_vlines_dat.
##' @export
#
#patient_plot_compact <- function(lrs, vlines){
#  lrs$assay <- sapply(strsplit(lrs$test, '_'), function(x){x[1]})
#  lrs$result <- sapply(strsplit(lrs$test, '_'), function(x){x[3]})
#  lrs$result[is.na(lrs$result)] <- '*'
#  x <- ggplot(lrs, aes(x = date, y = prob, col = result, group = facet_lab)) + 
#    facet_grid(rows = vars(assay)) + 
#    geom_line() +
#    theme(legend.position = 'none') +
#    labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis',
#         x = 'Date of intial infection') +
#    geom_vline(data = vlines, aes(xintercept = visit_date), col = 'black') +
#    scale_x_date(date_breaks = "months", date_labels = "%b-%y")
#}
#
##' Plots and annotates a patients timelines
##'
##' Plots the dates and the likelihoods that these dates are the DDI_1 for a
##' patient. Plots are annotated with the test dates as vertical lines.
##'
##' @param dat An interpretated result set in long format.
##' @export
#
#wrapped_patient_plot <- function(dat){
#  vlines <- make_vlines_dat(dat)
#  x <- patient_plot(dat, vlines)
#  return(x)
#}
#
