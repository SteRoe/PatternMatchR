% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\docType{data}
\name{originalWd}
\alias{originalWd}
\title{structure of sessionVariables:}
\format{
An object of class \code{character} of length 1.
}
\usage{
originalWd
}
\arguments{
\item{libname}{library name}

\item{pkgname}{package name}
}
\value{
nothing
}
\description{
very first function during package load
}
\details{
general settings
session$userData$sessionVariables$MaxProbes - maximum number of probes to use in clustering
session$userData$sessionVariables$P_ValMinBorder - minimum p-value to consider in clustering and to visualize
session$userData$sessionVariables$P_ValMaxBorder - maximum p-value to consider in clustering and to visualize
general data
session$userData$globalVariables$Beta_tDF - general beta values from Illumina 450k array
very first function
}
\keyword{##@noRd}
\keyword{internal}
\keyword{noRd}
