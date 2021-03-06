% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/result_interpretation.R
\name{trim_range}
\alias{trim_range}
\title{Trims range on which to evaluate a function}
\usage{
trim_range(fun, range_start, range_end, tol = 0)
}
\arguments{
\item{fun}{The function that will dictate which portions of the range will get trimmed.}

\item{range_start}{The lower bound of the initial range.}

\item{range_end}{The upper bound of the initial range.}

\item{tol}{How close to zero before we consider it equal to zero? Default = 0.}
}
\description{
Because the duration of the trial is long and often the region in which the infection could plausibly have happened is small, the integratoin functions often produces zero. Since we know that any aggregate function that produces a probable infection interval shorted than 1 day is problematic, we can devise an easy scheme to remove large chunks from the range where the interpreter function is zero.
}
\details{
Simply walk from the start to end in one day increments, if the likelihood is still zero, adjust the range_start to the current day. Likewise from range_end.
}
