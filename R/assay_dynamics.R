#' Linear Assay Dynamics
#'
#' Returns the probability of testing positive x days after infection
#'
#' A very simple scheme is used to compute these probabilities. Each assay is characterized by a diagnostic delay after which time 50 percent of patients will test positive. The assay also has a spread parameter. This parameter controls the window during which this probability is not 0 or 1. The default is 1.5 meaning that the period of time during which the assay result is not perfect is 50% longer than the diagnostic delay.
#'
#' @param x The time since infection
#' @param diagnostic_delay The number of days since infection after which time 50 percent of patients will test positive.
#' @param spread The multiplier applied to the diagnostic_delay to obtain the length of window centered at diagnostic_delay days after infection during some patients will have different results.
#' @param abs_spread The width of the range over which the probabilities are both not zero and not one. This range will be centered around the diagnostic_delay. Setting this will cause the function to ignore whatever was specified for spread. Defaults to NULL.
#' @export

linear_assay_dynamics <- function(x, diagnostic_delay, spread = 1.5, abs_spread = NULL){
  if (is.null(abs_spread)){ # using spread parameter
    y <- x/(spread*diagnostic_delay) + (0.5-(1/spread))
  } else { # using abs_spread parameter and ignoring spread parameter.
    y <- x/abs_spread + 0.5 - diagnostic_delay/abs_spread
  }
  y <- ifelse(y < 0, 0, y)
  y <- ifelse(y > 1, 1, y)
  return(y)
}

#' Step Assay Dynamics
#'
#' Returns the probability of testing positive x days after infection
#'
#' The simplest scheme is used to compute these probabilities. Each assay is characterized by a diagnostic delay after which time all of the results will be positive. Before this time, all results will be negative.
#'
#' @param x The time since infection
#' @param diagnostic_delay The number of days since infection after which time all the results will be positive.
#' @export

step_assay_dynamics <- function(x, diagnostic_delay){
  return(ifelse(x > diagnostic_delay, 1, 0))
}

#' Weibull3 Assay Dynamics
#'
#' Returns the probability of testing positive x days after infection
#'
#' Assumes that the probability of testing positive as a function of time given that the person is infected has the form of a Weibull distribution.
#'
#' @param x The time since infection
#' @param location The location parameter.
#' @param scale The scale parameter.
#' @param shape The shape parameter.
#' @export

weib3_assay_dynamics <- function(x, location, scale, shape){
  return(ifelse(x <= location, 0, 1-exp(-((x - location)/scale)^shape)))
}

#' Get dynamics functions for an assay
#'
#' Constructs the assay dynamics list with suitable parameter set given the name of an assay
#'
#' @param assay The name of the assay whose dynamics are required. If NULL (default), list of available assays will be returned.
#' @export

get_assay_dynamics <- function(assay = NULL){

  # see /fridge/data/tsic/characterization for details: TODO make this a vignette

  # Build library of assay parameters
  all_assays <- list()

  # BASICs for UNIT TESTING
  all_assays[[tolower('step_unit_testing')]] <- list(
    class = 'unit testing',
    full_assayname = 'Unit Testing',
    short_assayname = 'Unit_Testing',
    form = 'step',
    source = 'made-up',
    fun = 'step_assay_dynamics',
    params = list(diagnostic_delay = 10)
  )
  all_assays[[tolower('linear_unit_testing')]] <- list(
    class = 'unit testing',
    full_assayname = 'Unit Testing',
    short_assayname = 'Unit_Testing',
    form = 'linear_abs_spread',
    source = 'made-up',
    fun = 'linear_assay_dynamics',
    params = list(diagnostic_delay = 10, abs_spread = 10)
  )
  all_assays[[tolower('weib3_unit_testing')]] <- list(
    class = 'unit testing',
    full_assayname = 'Unit Testing',
    short_assayname = 'Unit_Testing',
    form = 'weib3',
    source = 'aptima_in_delaney_2017',
    fun = 'weib3_assay_dynamics',
    params = list(location = 4.8, shape = 1.35, scale = 9)
  )

  # ARCHITECT
  all_assays[[tolower('architect_step_delaney')]] <- list(
    class = 'Ag/Ab',
    full_assayname = 'Abbott Architect HIV Ag/Ab Combo',
    short_assayname = 'Architect',
    form = 'step',
    source = 'delaney_2017',
    fun = 'step_assay_dynamics',
    params = list(diagnostic_delay = 17.9)
  )
  all_assays[[tolower('architect_linear_abs_spread_delaney')]] <- list(
    class = 'Ag/Ab',
    full_assayname = 'Abbott Architect HIV Ag/Ab Combo',
    short_assayname = 'Architect',
    form = 'linear_abs_spread',
    source = 'delaney_2017',
    fun = 'linear_assay_dynamics',
    params = list(diagnostic_delay = 17.9, abs_spread = 2*9.5)
  )
  all_assays[[tolower('architect_weib3_delaney')]] <- list(
    class = 'Ag/Ab',
    full_assayname = 'Abbott Architect HIV Ag/Ab Combo',
    short_assayname = 'Architect',
    form = 'weib3',
    source = 'delaney_2017',
    fun = 'weib3_assay_dynamics',
    params = list(location = 7.209, shape = 1.725, scale = 13.185)
  )

  # APTIMA
  all_assays[[tolower('aptima_step_delaney')]] <- list(
    class = 'RNA',
    full_assayname = 'Aptima HIV-1 RNA Qualitative Assay',
    short_assayname = 'Aptima RNA',
    form = 'step',
    source = 'delaney_2017',
    fun = 'step_assay_dynamics',
    params = list(diagnostic_delay = 11.679)
  )
  all_assays[[tolower('aptima_linear_abs_spread_delaney')]] <- list(
    class = 'RNA',
    full_assayname = 'Aptima HIV-1 RNA Qualitative Assay',
    short_assayname = 'Aptima RNA',
    form = 'linear_abs_spread',
    source = 'delaney_2017',
    fun = 'linear_assay_dynamics',
    params = list(diagnostic_delay = 11.679, abs_spread = 2*7.875)
  )
  all_assays[[tolower('aptima_weib3_delaney')]] <- list(
    class = 'RNA',
    full_assayname = 'Aptima HIV-1 RNA Qualitative Assay',
    short_assayname = 'Aptima RNA',
    form = 'weib3',
    source = 'delaney_2017',
    fun = 'weib3_assay_dynamics',
    params = list(location = 4.8, shape = 1.35, scale = 9)
  )

  # Abbott real time
  all_assays[[tolower('abbott_real_time_weib3_delaney_and_manufacturer')]] <- list(
    class = 'RNA',
    full_assayname = 'Abbott Real Time HIV01 v1.0 m2000sp/m2000rt',
    short_assayname = 'Abbott Real Time',
    form = 'weib3',
    source = "manufacturer's details and delaney_2017",
    fun = 'weib3_assay_dynamics',
    params = list(location = 5.1252, shape = 1.350, scale = 9.660)
  )

  # Taqman
  all_assays[[tolower('taqman_weib3_delaney_and_manufacturer')]] <- list(
    class = 'RNA',
    full_assayname = 'Roche Taqman v2.0',
    short_assayname = 'Taqman',
    form = 'weib3',
    source = "manufacturer's details and delaney_2017",
    fun = 'weib3_assay_dynamics',
    params = list(location = 4.995, shape = 1.350, scale = 9.365)
  )

  # iSCAv2 (Mellor's low copy)
  all_assays[[tolower('iSCAv2_weib3_delaney_and_tosiano')]] <- list(
    class = 'RNA',
    full_assayname = 'iSCA v2.0',
    short_assayname = 'iSCAv2',
    form = 'weib3',
    source = "tosiano_2019 and delaney_2017",
    fun = 'weib3_assay_dynamics',
    params = list(location = 3.606, shape = 1.350, scale = 6.762)
  )

  # Geenius Fully Reactive
  all_assays[[tolower('geenius_fr_step_delaney')]] <- list(
    class = 'IgG_Supp',
    full_assayname = 'BioRad Geenius Fully Reactive',
    short_assayname = 'Geenius_FR',
    form = 'step',
    source = 'delaney_2017',
    fun = 'step_assay_dynamics',
    params = list(diagnostic_delay = 32.9)
  )
  all_assays[[tolower('geenius_fr_linear_abs_spread_delaney')]] <- list(
    class = 'IgG_Supp',
    full_assayname = 'BioRad Geenius Fully Reactive',
    short_assayname = 'Geenius_FR',
    form = 'linear_abs_spread',
    source = 'delaney_2017',
    fun = 'linear_assay_dynamics',
    params = list(diagnostic_delay = 32.9, abs_spread = 2*10.4)
  )
  all_assays[[tolower('geenius_fr_weib3_delaney')]] <- list(
    class = 'IgG_Supp',
    full_assayname = 'BioRad Geenius Fully Reactive',
    short_assayname = 'Geenius_FR',
    form = 'weib3',
    source = 'delaney_2017',
    fun = 'weib3_assay_dynamics',
    params = list(location = 21.151, shape = 1.733, scale = 14.483)
  )

  # Geenius indeterminate
  all_assays[[tolower('geenius_indet_weib3_delaney')]] <- list(
    class = 'IgG_Rapid',
    full_assayname = 'BioRad Geenius Indeterminate',
    short_assayname = 'Geenius_Indet',
    form = 'weib3',
    source = 'delaney_2017_relationship_to_multispot',
    fun = 'weib3_assay_dynamics',
    params = list(location = 16.556, shape = 1.773, scale = 12.9)
  )

  # GS COMBO
  all_assays[[tolower('gs_combo_weib3_delaney')]] <- list(
    class = 'Ag/Ab',
    full_assayname = 'BioRad GS HIV Combo Ag/Ab EIA',
    short_assayname = 'GS_Combo',
    form = 'weib3',
    source = 'delaney_2017',
    fun = 'weib3_assay_dynamics',
    params = list(location = 9.9675, shape = 0.7206, scale = 7.6388)
  )

  # lookup
  if (is.null(assay)){
    return(names(all_assays))
  }
  if (assay == 'all'){
    return(all_assays)
  }
  if (assay %in% names(all_assays)){
    return(all_assays[[assay]])
  } else {
    print('Valid Assay Names:')
    print(names(all_assays))
    stop('Requested Assay Do Not Exists')
  }
}

