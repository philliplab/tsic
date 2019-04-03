#' Create Dataset of test dates
#'
#' Prepares a dataset with the dates of the tests for easy inclusion in plots.
#'
#' @param in_dat An interpretated datasets in long format
#' @export

make_vlines_dat <- function(in_dat){
  vlines <- in_dat[,c('assay', 'visit_date', 'result')]
  vlines %>% map_if(is.factor, as.character) %>% as_tibble -> vlines
  vlines$facet_lab <- paste(vlines$assay, '\n', gsub('-', '', vlines$visit_date), '\n', vlines$result, sep = '')
  vlines$visit_date <- as_date(vlines$visit_date)
  return(vlines)
}


#' Plots a patients timelines
#'
#' Plots dates and the likelihoods that these dates are the DDI_1 for a patient.
#'
#' @param lrs An interpretated result set in long format.
#' @param vlines Vertical lines indicating the test dates as produced by make_vlines_dat.
#' @export

patient_plot <- function(lrs, vlines){
  x <- ggplot(lrs, aes(x = date, y = prob, col = facet_lab)) + 
    facet_grid(rows = vars(facet_lab)) + 
    geom_line() +
    theme(legend.position = 'none') +
    labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis',
         x = 'Date of intial infection') +
    geom_vline(data = vlines, aes(xintercept = visit_date), col = 'black') +
    scale_x_date(date_breaks = "months", date_labels = "%b-%y")
  return(x)
}

#' Plots and annotates a patients timelines
#'
#' Plots the dates and the likelihoods that these dates are the DDI_1 for a
#' patient. Plots are annotated with the test dates as vertical lines.
#'
#' @param dat An interpretated result set in long format.
#' @export

wrapped_patient_plot <- function(dat){
  vlines <- make_vlines_dat(dat)
  lrs <- make_lrs_dat(dat)
  x <- patient_plot(lrs, vlines)
  return(x)
}

