

plot_iihist <- function(iihist, lb_med_ub = NULL, produce_plot = TRUE){
  # debugging stuff
  if (FALSE) {
    ihist <- data.frame(
      ptid = c('p0', 'p0'),
      sample_date = c('2016-03-01', '2016-09-01'),
      test = c('step_unit_testing', 'step_unit_testing'),
      result = c('-', '+'),
      stringsAsFactors = FALSE
    )
    iihist <- interpret_ihist(ihist = ihist,
                              range_start = as.Date('2016-01-01'),
                              range_end = as.Date('2016-11-30'))
  }
  stopifnot('Aggregate' %in% iihist$test_details)

  iihist$result <- gsub('(.*)(.)$', '\\2', iihist$test_details)
  iihist$test_details <- 
    factor(iihist$test_details, 
           levels = c('Aggregate', setdiff(sort(unique(iihist$test_details)), 'Aggregate')),
           ordered = TRUE)

  x_breaks <- quantile(seq(from = min(iihist$sample_date), 
                           to = max(iihist$sample_date), 
                           length.out = 100), 
                       c(.1, .3, .5, .7, .9))

  x_tick_labels <- as.character(as.Date(round(x_breaks, 0), origin = '1970-01-01'))
  #x_tick_labels <- strftime(as.Date(round(x_breaks, 0), origin = '1970-01-01'), format = "%b-%y")

  x <- 
  ggplot2::ggplot(iihist, aes(x = sample_date, y = prob_val, col = result)) +
    geom_line() +
    facet_grid(test_details ~ .) +
    theme(legend.position = 'none') +
    labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis') +
    scale_x_continuous('Date of intial infection', breaks = x_breaks, labels = x_tick_labels)
  if (produce_plot) {print(x)}
  return(x)
}

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
