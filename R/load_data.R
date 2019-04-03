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
