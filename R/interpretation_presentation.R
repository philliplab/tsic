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
