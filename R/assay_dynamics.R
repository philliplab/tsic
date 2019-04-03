#' Assay dynamics (assuming linear behaviour)
#'
#' Returns the probability of testing positive x days after DDI_1
#'
#' A very simple scheme is used to compute these probabilities. Each assay is characterized by a diagnostic delay after which time 50% of patients will test positive. The assay also has a spread parameter. This parameter controls the window during which this probability is not 0 or 1. The default is 1.5 meaning that the period of time during which the assay result is not perfect is 50% longer than the diagnostic delay.
#'
#' @param x The time since DDI_1
#' @param diagnostic_delay The number of days since DDI_1 after which time 50% of patients will test positive.
#' @param spread The multiplier applied to the diagnostic_delay to obtain the length of window centered at diagnostic_delay days after DDI_1 during some patients will have different results.
#' @export

linear_assay_dynamics <- function(x, diagnostic_delay, spread = 1.5){
  # first model:
  # probability of positive test diagnostic_delay days after ddi1 = 0.5
  # probability of positive test (1-spread/2)*diagnostic_delay days after ddi1 = 0.
  # A straight line through these points gives the probabilities:
  #  if prob > 1, then return 1
  #  if prob < 0, then return 0

  y <- x/(spread*diagnostic_delay) + (0.5-(1/spread))
  if (y < 0) {y <- 0}
  if (y > 1) {y <- 1}
  return(y)
}

#' Get dynamics functions for an assay
#'
#' Constructs the assay dynamics list with suitable parameter set given the name of an assay
#'
#' @param assay The name of the assay whose dynamics are required.
#' @export

get_assay_dynamics <- function(assay){
  if (assay %in% c('Abbott Architect HIV Ag/Ab Combo', 'elisa')){
    return(list(fun = get('linear_assay_dynamics'),
                params = list(diagnostic_delay = 10.8)))
  } else if (assay %in% c('rnapcr')){
    return(list(fun = get('linear_assay_dynamics'),
                params = list(diagnostic_delay = 10.8)))
  } else if (assay %in% c('geenius')){
    return(list(fun = get('linear_assay_dynamics'),
                params = list(diagnostic_delay = 10.8)))
  } else if (assay %in% c('totalnucleicacid')){
    return(list(fun = get('linear_assay_dynamics'),
                params = list(diagnostic_delay = 10.8)))
  } else {
    stop(paste(assay, ' is not available'))
  }
}


## old assay dynamics
#
#if (FALSE){
#
#totalnucleicacid_dynamics <- matrix(data = c(0, 0, #{{{
#                                             1, 0,
#                                             2, 0.01,
#                                             3, 0.05,
#                                             4, 0.275,
#                                             5, 0.5,
#                                             6, 0.725,
#                                             7, 0.95,
#                                             8, 0.99,
#                                             9, 1,
#                                             10, 1,
#                                             11, 1,
#                                             12, 1,
#                                             13, 1,
#                                             14, 1,
#                                             15, 1,
#                                             16, 1,
#                                             17, 1,
#                                             18, 1,
#                                             19, 1,
#                                             20, 1,
#                                             21, 1,
#                                             22, 1,
#                                             23, 1,
#                                             24, 1,
#                                             25, 1,
#                                             26, 1,
#                                             27, 1,
#                                             28, 1,
#                                             29, 1,
#                                             30, 1),
#                                             ncol = 2,
#                                             byrow = TRUE) #}}}
#
#rnapcr_dynamics <- matrix(data = c(0, 0, #{{{
#                                   1, 0,
#                                   2, 0,
#                                   3, 0,
#                                   4, 0,
#                                   5, 0,
#                                   6, 0.01,
#                                   7, 0.05,
#                                   8, 0.1785714,
#                                   9, 0.3071429,
#                                   10, 0.4357143,
#                                   11, 0.5642857,
#                                   12, 0.6928571,
#                                   13, 0.8214286,
#                                   14, 0.95,
#                                   15, 0.99,
#                                   16, 1,
#                                   17, 1,
#                                   18, 1,
#                                   19, 1,
#                                   20, 1,
#                                   21, 1,
#                                   22, 1,
#                                   23, 1,
#                                   24, 1,
#                                   25, 1,
#                                   26, 1,
#                                   27, 1,
#                                   28, 1,
#                                   29, 1,
#                                   30, 1),
#                                   ncol = 2,
#                                   byrow = TRUE) #}}}
#
#elisa_dynamics <- matrix(data = c(0, 0, #{{{
#                                   1, 0,
#                                   2, 0,
#                                   3, 0,
#                                   4, 0,
#                                   5, 0,
#                                   6, 0,
#                                   7, 0,
#                                   8, 0,
#                                   9, 0,
#                                   10, 0,
#                                   11, 0,
#                                   12, 0.01,
#                                   13, 0.05,
#                                   14, 0.1318182,
#                                   15, 0.2136364,
#                                   16, 0.2954545,
#                                   17, 0.3772727,
#                                   18, 0.4590909,
#                                   19, 0.5409091,
#                                   20, 0.6227273,
#                                   21, 0.7045455,
#                                   22, 0.7863636,
#                                   23, 0.8681818,
#                                   24, 0.95,
#                                   25, 0.99,
#                                   26, 1,
#                                   27, 1,
#                                   28, 1,
#                                   29, 1,
#                                   30, 1),
#                                   ncol = 2,
#                                   byrow = TRUE) #}}}
#
#westernblot_dynamics <- matrix(data = c(0, 0, #{{{
#                                        1, 0,
#                                        2, 0,
#                                        3, 0,
#                                        4, 0,
#                                        5, 0,
#                                        6, 0,
#                                        7, 0,
#                                        8, 0,
#                                        9, 0,
#                                        10, 0,
#                                        11, 0.01,
#                                        12, 0.05,
#                                        13, 0.1318182,
#                                        14, 0.2136364,
#                                        15, 0.2954545,
#                                        16, 0.3772727,
#                                        17, 0.4590909,
#                                        18, 0.5409091,
#                                        19, 0.6227273,
#                                        20, 0.7045455,
#                                        21, 0.7863636,
#                                        22, 0.8681818,
#                                        23, 0.95,
#                                        24, 0.99,
#                                        25, 1,
#                                        26, 1,
#                                        27, 1,
#                                        28, 1,
#                                        29, 1,
#                                        30, 1),
#                                        ncol = 2,
#                                        byrow = TRUE) #}}}
#
#geenius_dynamics <- matrix(data = c(0, 0, #{{{
#                                    1, 0,
#                                    2, 0,
#                                    3, 0,
#                                    4, 0,
#                                    5, 0,
#                                    6, 0,
#                                    7, 0,
#                                    8, 0,
#                                    9, 0,
#                                    10, 0.025,
#                                    11, 0.050,
#                                    12, 0.075,
#                                    13, 0.100,
#                                    14, 0.125,
#                                    15, 0.150,
#                                    16, 0.175,
#                                    17, 0.200,
#                                    18, 0.225,
#                                    19, 0.250,
#                                    20, 0.275,
#                                    21, 0.300,
#                                    22, 0.325,
#                                    23, 0.350,
#                                    24, 0.375,
#                                    25, 0.400,
#                                    26, 0.425,
#                                    27, 0.450,
#                                    28, 0.475,
#                                    29, 0.500,
#                                    30, 0.525,
#                                    31, 0.550,
#                                    32, 0.575,
#                                    33, 0.600,
#                                    34, 0.625,
#                                    35, 0.650,
#                                    36, 0.675,
#                                    37, 0.700,
#                                    38, 0.725,
#                                    39, 0.750,
#                                    40, 0.775,
#                                    41, 0.800,
#                                    42, 0.825,
#                                    43, 0.850,
#                                    44, 0.875,
#                                    45, 0.900,
#                                    46, 0.925,
#                                    47, 0.950,
#                                    48, 0.975,
#                                    49, 1,
#                                    50, 1,
#                                    51, 1,
#                                    52, 1,
#                                    53, 1,
#                                    54, 1,
#                                    55, 1,
#                                    56, 1,
#                                    57, 1,
#                                    58, 1,
#                                    59, 1,
#                                    60, 1),
#                                    ncol = 2,
#                                    byrow = TRUE) #}}}
#
#geenius_dynamics_old <- matrix(data = c(0, 0, #{{{
#                                        1, 0,
#                                        2, 0,
#                                        3, 0,
#                                        4, 0,
#                                        5, 0,
#                                        6, 0,
#                                        7, 0,
#                                        8, 0,
#                                        9, 0,
#                                        10, 0,
#                                        11, 0,
#                                        12, 0,
#                                        13, 0.01,
#                                        14, 0.05,
#                                        15, 0.1318182,
#                                        16, 0.2136364,
#                                        17, 0.2954545,
#                                        18, 0.3772727,
#                                        19, 0.4590909,
#                                        20, 0.5409091,
#                                        21, 0.6227273,
#                                        22, 0.7045455,
#                                        23, 0.7863636,
#                                        24, 0.8681818,
#                                        25, 0.95,
#                                        26, 0.99,
#                                        27, 1,
#                                        28, 1,
#                                        29, 1,
#                                        30, 1),
#                                        ncol = 2,
#                                        byrow = TRUE) #}}}
#
#all_assay_dynamics <- rbind( #{{{
#  data.frame(assay = 'totalnucleicacid',
#             days = totalnucleicacid_dynamics[,1],
#             prob = totalnucleicacid_dynamics[,2]),
#  data.frame(assay = 'rnapcr',
#             days = rnapcr_dynamics[,1],
#             prob = rnapcr_dynamics[,2]),
#  data.frame(assay = 'elisa',
#             days = elisa_dynamics[,1],
#             prob = elisa_dynamics[,2]),
#  data.frame(assay = 'westernblot',
#             days = westernblot_dynamics[,1],
#             prob = westernblot_dynamics[,2]),
#  data.frame(assay = 'geenius',
#             days = geenius_dynamics[,1],
#             prob = geenius_dynamics[,2])) #}}}
#
#}
