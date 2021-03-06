% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_RegionSummaryMethod.R
\name{RegionSummaryMethod}
\alias{RegionSummaryMethod}
\alias{MaxSites}
\alias{MeanSites}
\alias{MedianSites}
\alias{PC1Sites}
\title{Methods to summarize RNA editing levels from multiple sites within a 
  single region.}
\usage{
MaxSites(rnaEditMatrix, ...)

MeanSites(rnaEditMatrix, ...)

MedianSites(rnaEditMatrix, ...)

PC1Sites(rnaEditMatrix, ...)
}
\arguments{
\item{rnaEditMatrix}{A matrix (or data frame) of RNA editing level values on
individual sites, with row names as site IDs in the form of
"chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
follow the format of example dataset (\code{data(rnaedit_df)}).}

\item{...}{Dots for additional internal arguments (currently unused).}
}
\value{
A named numeric vector of summarized RNA editing levels with sample
  IDs as names.
}
\description{
Summarize RNA editing sites in a single region by taking 
  maximum, mean, median or first principal component.
}
\examples{
  data(rnaedit_df)
  MedianSites(rnaEditMatrix = rnaedit_df)[1:3]

}
\keyword{internal}
