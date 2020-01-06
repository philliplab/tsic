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

#' Loads dataset used for DSMB Nov 2019
#'
#' @param file_name Full path to the csv file
#' @export

load_dsmb_nov_2019_data <- function(file_name){
  if (FALSE){
    file_name <- '/fridge/data/AMP/DSMB_timing_nov_2019/AMP_diagnostic_testing_history_DSMB_2019_Nov.csv'
    file_name <- '/fridge/data/AMP/misc_timing_runs/erika_2019-12-30/v704_nonmitt_diagnostic_testing_history.csv'
  }
  dat <- read.csv(file_name, stringsAsFactors = FALSE)
  names(dat) <- c('ptid', 'sample_date', 'test', 'result')
  dat$ptid <- paste('p_', dat$ptid, sep = '')
  dat$sample_date <- as.numeric(as.Date(dat$sample_date)) + 0.5
  stopifnot(all(sort(unique(dat$result)) == c("Negative", "Positive")))
  dat$result <- ifelse(dat$result == 'Negative', '-', '+')

  valid_tests <- c(
    "Abbott ARCHITECT HIV Ag/Ab Combo", "Abbott Real Time HIV01 v1.0 m2000sp/m2000rt",
    "BioRad Geenius Fully Reactive", "BioRad Geenius Indeterminate",
    "BioRad GS HIV Combo Ag/Ab EIA", "Roche Taqman v2.0",
    "Alere Determine Rapid-HIV 1/2 Test", "OraSure OraQuick ADVANCE whole blood"
  )
  input_rows <- nrow(dat)
  dat <- subset(dat, test %in% valid_tests)
  valid_test_rows <- nrow(dat)

  test_mapping <- matrix(
c(
"Abbott ARCHITECT HIV Ag/Ab Combo", 'architect_weib3_delaney',
"Abbott Real Time HIV01 v1.0 m2000sp/m2000rt", 'abbott_real_time_weib3_delaney_and_manufacturer',
"BioRad Geenius Fully Reactive", 'geenius_fr_weib3_delaney',
"BioRad Geenius Indeterminate", 'geenius_indet_weib3_delaney',
"BioRad GS HIV Combo Ag/Ab EIA", 'gs_combo_weib3_delaney',
"Roche Taqman v2.0", 'taqman_weib3_delaney_and_manufacturer',
"Alere Determine Rapid-HIV 1/2 Test", "determine_weib3_delaney",
"OraSure OraQuick ADVANCE whole blood", "oraquick_weib3_delaney"
),
    ncol = 2,
    byrow = TRUE)
  test_mapping <- data.frame(
    AMP_name = test_mapping[,1],
    tsic_name = test_mapping[,2],
    stringsAsFactors = FALSE)

  dat <- merge(dat, test_mapping,
    by.x = 'test', by.y = 'AMP_name')
  test_mapped_rows <- nrow(dat)
  stopifnot(test_mapped_rows == valid_test_rows)

  dat <- dat[,c('ptid', 'sample_date', 'tsic_name', 'result')]
  names(dat) <- c('ptid', 'sample_date', 'test', 'result')
  return(dat)
}




