% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResultData.R
\name{loadResultDF}
\alias{loadResultDF}
\title{structure of sessionVariables:}
\usage{
loadResultDF(session, folder, loadRDS = FALSE)
}
\arguments{
\item{session}{session object}

\item{folder}{folder containing all files to read to list}

\item{loadRDS}{loads RDS file or not (if only RDS files should be written)}
}
\value{
result list()
data.frame P_Val,
data.frame DM (DeltaMethylation),
data.frame N,
data.frame PHENOdata,
chr PHENOFileName
}
\description{
result data to be analyzed loaded from regression results
session$userData$sessionVariables$resultDFListTrait1()$resultDFP_Val - list of resulting p-values
session$userData$sessionVariables$resultDFListTrait1()$resultDFDM - list of resulting delta methylations
session$userData$sessionVariables$resultDFListTrait1()$resultDFN - list of resulting n
session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[\link{x}]$PHENODF - original data, the regression is based on
session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[\link{x}]$PHENOFileName - file name of the original data
PHENOdata <- list(PHENODF = getPHENODF(PHENOFileName),PHENOFileName = PHENOFileName)
loadResultDF
loads data.frame from folder
}
