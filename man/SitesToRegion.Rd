% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_SitesToRegion.R
\name{SitesToRegion}
\alias{SitesToRegion}
\title{Create output data in the format of GRanges.}
\usage{
SitesToRegion(
  sitesSubregion_df,
  sitesAreOrdered = TRUE,
  keepminPairwiseCor_df,
  returnAllSites = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{sitesSubregion_df}{An output data frame from function
\code{FindCorrelatedRegions} with variables \code{site, subregion}. Please
see \code{\link{FindCorrelatedRegions}} for details.}

\item{sitesAreOrdered}{Are the sites in \code{sitesSubregion_df} ordered by 
location? Defaults to FALSE.}

\item{keepminPairwiseCor_df}{An output data frame from function
\code{GetMinPairwiseCor} with variables \code{subregion},
\code{keepminPairwiseCor} and \code{minPairwiseCor}. Please see
\code{\link{GetMinPairwiseCor}} for details.}

\item{returnAllSites}{When no contiguous co-edited regions are found in
a input genomic region, \code{returnAllSites = TRUE} indicates
outputting all the sites in this input region, while
\code{returnAllSites = FALSE} indicates not returning any site in this
input region. Defaults to FALSE.}

\item{verbose}{Should messages and warnings be displayed? Defaults to TRUE.}
}
\value{
A GRanges object with \code{seqnames}, \code{ranges} and
  \code{strand} of the contiguous co-edited regions.
}
\description{
Output contiguous co-edited subregions found by
  \code{\link{FindCorrelatedRegions}} function and filtered by
  \code{\link{GetMinPairwiseCor}} function.
}
\examples{
  data(t_rnaedit_df)
  
  ordered_cols <- OrderSitesByLocation(
    sites_char = colnames(t_rnaedit_df),
    output = "vector"
  )
  exm_data <- t_rnaedit_df[, ordered_cols]
  
  exm_sites <- MarkCoeditedSites(
    rnaEditCluster_mat = exm_data,
    method = "spearman"
  )
  
  exm_regions <- FindCorrelatedRegions(
    sites_df = exm_sites,
    featureType = "site"
  )
  
  exm_sites <- split(
    x = exm_regions$site,
    f = exm_regions$subregion
  )
  
  exm_cor <- GetMinPairwiseCor(
    rnaEditCluster_mat = exm_data,
    minPairCorr = 0.1,
    probes_ls = exm_sites,
    method = "spearman"
  )
  
  SitesToRegion(
    sitesSubregion_df = exm_regions,
    keepminPairwiseCor_df = exm_cor$keepminPairwiseCor_df
  )
   
}
\keyword{internal}
