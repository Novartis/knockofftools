% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{find_ratio}
\alias{find_ratio}
\title{Ratio from Ren et al. (2021) used to derandomized knockoffs.}
\usage{
find_ratio(M, eta)
}
\arguments{
\item{M}{number of knockoff runs}

\item{eta}{the proportion of runs in which a feature
must be selected to be selected in the overall derandomized
procedure.}
}
\value{
returns pre-calculated ratios to be used in the multiple knockoff variable selection with error.type = "pfer" or "kfwer"
}
\description{
Do not call this function on its own.
}
\keyword{internal}
