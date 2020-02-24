#' Window Period Distributions for Assays
#'
#' A list of the window period distributions, class, data source and various names of the assays included in tsic. The data is mostly sourced from Delaney 2017, but augmented with some extra computations and data sources.
#'
#' @format A list of lists each of which contains the details of a version of an assay. Each assay can appear multiple times because different assumptions can be made for the window period distributions. 
#' \describe{
#'   \item{class}{The class of the assay, mostly following the conventions of Delaney 2017}
#'   \item{full_assayname}{The full name of the assay}
#'   \item{short_assayname}{A Shortened name used to annotate figures}
#'   \item{form}{The name and sometimes description of the distribution assumed}
#'   \item{source}{Source(s) of the data}
#'   \item{fun}{The name of the function that implements the distribution assumed for this assay}
#'   \item{params}{The parameters of the window period distribution. A list whose element names depend on the distribution assumed}
#' }
"all_assay_dynamics"
