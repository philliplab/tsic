#' Plots an individual history
#'
#' Will produce the interpretation itself.
#'
#' @export

plot_iihist <- function(ihist, lb_med_ub, range_start, range_end, 
                        x_breaks = NULL,
                        produce_plot = TRUE, save_plot = FALSE, 
                        verbose = FALSE, plot_aggregate = TRUE,
                        show_test_dates = FALSE,
                        custom_aggregate_label = NULL,
                        scales = 'fixed'){
  # debugging stuff
  if (FALSE) {
    lb_med_ub <- TRUE
    plot_aggregate <- TRUE
    produce_plot <- TRUE
    verbose <- TRUE
    show_test_dates <- TRUE
    x_breaks <- NULL
   
    custom_aggregate_label <- 'Aggregate\n5 repeats'
    plot_aggregate <- FALSE
    lb_med_ub <- NULL
    devtools::load_all()
    library(profvis)
    #1
    ihist <- data.frame(
      ptid = c('p0', 'p0'),
      sample_date = c(as.numeric(as.Date('2016-03-01')), as.numeric(as.Date('2016-09-01'))),
      test = c('step_unit_testing', 'step_unit_testing'),
      result = c('-', '+'),
      stringsAsFactors = FALSE
    )
    iihist <- interpret_ihist(ihist = ihist,
                              range_start = as.numeric(as.Date('2016-01-01')),
                              range_end = as.numeric(as.Date('2016-11-30')))

    # 2
    ihist <- read.csv('/fridge/data/tsic/test_data.csv', stringsAsFactors = FALSE)
    ihist$sample_date <- as.numeric(as.Date(ihist$sample_date))
    ihist <- subset(ihist, ptid == 'p01')
    range_start <- as.numeric(as.Date('2017-01-01'))
    range_end <- as.numeric(as.Date('2017-11-30'))
#    profvis({
    iihist <- interpret_ihist(ihist = ihist,
                              range_start = as.numeric(as.Date('2017-01-01')),
                              range_end = as.numeric(as.Date('2017-11-30')),
                              verbose = TRUE)
#    })
    agg_interpreter <- construct_aggregate_interpreter(ihist)
    range_start <- min(ihist$sample_date) - 100
    range_end <- max(ihist$sample_date) + 100
    range_to_int <- trim_range(fun = agg_interpreter, range_start = range_start, range_end = range_end)
    lb_med_ub <- estimate_lb_med_ub(fun = agg_interpreter,                            
                                    range_start = range_to_int$range_start,
                                    range_end = range_to_int$range_end,
                                    verbose = TRUE)
  } # end of debugging stuff

  stopifnot(scales %in% c('fixed', 'free', 'free_y', 'free_x'))

  iihist <- interpret_ihist(ihist = ihist,
                            range_start = range_start,
                            range_end = range_end,
                            verbose = FALSE)
  stopifnot('Aggregate' %in% iihist$test_details)
  if (verbose){print(str(iihist))}

  if (!is.null(lb_med_ub) & plot_aggregate){
    vlines_dat <- data.frame(ptid = unique(iihist$ptid))
  }

  iihist$result <- gsub('(.*)(.)$', '\\2', iihist$test_details)
  iihist$test_details <- 
    factor(iihist$test_details, 
           levels = c('Aggregate', setdiff(sort(unique(iihist$test_details)), 'Aggregate')),
           ordered = TRUE)
  if (verbose){print(str(iihist))}

  if (is.null(x_breaks)){
    x_breaks <- quantile(seq(from = min(iihist$sample_date), 
                             to = max(iihist$sample_date), 
                             length.out = 100), 
                         c(.1, .3, .5, .7, .9))
  } else {
    stopifnot(min(x_breaks) >= range_start)
    stopifnot(max(x_breaks) <= range_end)
  }

  x_tick_labels <- as.character(as.Date(round(x_breaks, 0), origin = '1970-01-01'))
  
  if (!is.null(lb_med_ub) & plot_aggregate){
    vlines_dat <- data.frame(ptid = unique(iihist$ptid),
                             sample_date = c(lb_med_ub$lb, lb_med_ub$med, lb_med_ub$ub),
                             test_details = 'Aggregate',
                             prob_val = NA,
                             result = c('lb', 'med', 'ub'),
                             stringsAsFactors = FALSE)
    vlines_dat$test_details <- 
      factor(vlines_dat$test_details, 
             levels = c('Aggregate', setdiff(sort(unique(vlines_dat$test_details)), 'Aggregate')),
             ordered = TRUE)
  }
  if (verbose){print(str(iihist))}

  if (show_test_dates){
    stopifnot(show_test_dates > 0)
    test_date_arrows <- data.frame(test_details = setdiff(levels(iihist$test_details), 'Aggregate'))
    test_date_arrows$sample_date <- as.numeric(as.Date(gsub('.*\n([0-9]{4}-[0-9]{2}-[0-9]{2})\n.*', '\\1', test_date_arrows$test_details)))
    test_date_arrows$result <- gsub('.*\n[0-9]{4}-[0-9]{2}-[0-9]{2}\n(.*)$', '\\1', test_date_arrows$test_details)
  }
  cols <- c("+" = rgb(1,0,0), "-" = rgb(0, 176/255, 240/255), "e" = "black")

  if (plot_aggregate){
    iihist_tmp <- iihist
    if (!is.null(custom_aggregate_label)){
      levels(iihist_tmp$test_details)[levels(iihist_tmp$test_details)=='Aggregate'] <- custom_aggregate_label
    }

    if(verbose){print(str(iihist_tmp))}

    x <- 
    ggplot2::ggplot(iihist_tmp, ggplot2::aes(x = sample_date, y = prob_val, col = result)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(test_details ~ ., scales = scales) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis') +
      ggplot2::scale_x_continuous('Date of intial infection', breaks = x_breaks, labels = x_tick_labels) +
      ggplot2::scale_colour_manual(values = cols)
    if (!is.null(lb_med_ub)){
      vlines_dat_tmp <- vlines_dat
      if (!is.null(custom_aggregate_label)){
        levels(vlines_dat_tmp$test_details)[levels(vlines_dat_tmp$test_details)=='Aggregate'] <- custom_aggregate_label
      }
      x <- x + ggplot2::geom_vline(data = vlines_dat_tmp, ggplot2::aes(xintercept = sample_date), col = 'black')
    }
  } else {
    x <- 
    ggplot2::ggplot(subset(iihist, test_details != 'Aggregate'), ggplot2::aes(x = sample_date, y = prob_val, col = result)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(test_details ~ ., scales = scales) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::labs(y = 'Probability of observed result given initial\ninfection on day indicated by x-axis') +
      ggplot2::scale_x_continuous('Date of intial infection', breaks = x_breaks, labels = x_tick_labels) +
      ggplot2::scale_colour_manual(values = cols)
  }
  if (show_test_dates){
    x <- x + ggplot2::geom_segment(data = test_date_arrows, 
                                   ggplot2::aes(x = sample_date, xend = sample_date, y = 0.333, yend = 0),
                                   arrow = ggplot2::arrow(length = ggplot2::unit(0.10, 'npc')),
                                   size = show_test_dates)
  }
  if (produce_plot) {print(x)}
  return(x)
}

