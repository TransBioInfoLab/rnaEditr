#' Create output data in the format of GRanges.
#' 
#' @description Output contiguous co-edited subregions found by
#'   \code{\link{FindCorrelatedRegions}} function and filtered by
#'   \code{\link{GetMinPairwiseCor}} function. 
#'
#' @param sitesSubregion_df An output data frame from function
#'   \code{FindCorrelatedRegions} with variables \code{site, subregion}. Please
#'   see \code{\link{FindCorrelatedRegions}} for details.
#' @param sitesAreOrdered Are the sites in \code{sitesSubregion_df} ordered by 
#'   location? Defaults to FALSE.
#' @param keepminPairwiseCor_df An output data frame from function
#'   \code{GetMinPairwiseCor} with variables \code{subregion},
#'   \code{keepminPairwiseCor} and \code{minPairwiseCor}. Please see
#'   \code{\link{GetMinPairwiseCor}} for details.
#' @param returnAllSites When no contiguous co-edited regions are found in
#'   a input genomic region, \code{returnAllSites = TRUE} indicates
#'   outputting all the sites in this input region, while
#'   \code{returnAllSites = FALSE} indicates not returning any site in this
#'   input region. Defaults to FALSE.
#' @param verbose Should messages and warnings be displayed? Defaults to TRUE.
#'
#' @return A GRanges object with \code{seqnames}, \code{ranges} and
#'   \code{strand} of the contiguous co-edited regions.
#' 
#' @importFrom stats aggregate
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   data(t_rnaedit_df)
#'   
#'   ordered_cols <- OrderSitesByLocation(
#'     sites_char = colnames(t_rnaedit_df),
#'     output = "vector"
#'   )
#'   exm_data <- t_rnaedit_df[, ordered_cols]
#'   
#'   exm_sites <- MarkCoeditedSites(
#'     rnaEditCluster_mat = exm_data,
#'     method = "spearman"
#'   )
#'   
#'   exm_regions <- FindCorrelatedRegions(
#'     sites_df = exm_sites,
#'     featureType = "site"
#'   )
#'   
#'   exm_sites <- split(
#'     x = exm_regions$site,
#'     f = exm_regions$subregion
#'   )
#'   
#'   exm_cor <- GetMinPairwiseCor(
#'     rnaEditCluster_mat = exm_data,
#'     minPairCorr = 0.1,
#'     probes_ls = exm_sites,
#'     method = "spearman"
#'   )
#'   
#'   SitesToRegion(
#'     sitesSubregion_df = exm_regions,
#'     keepminPairwiseCor_df = exm_cor$keepminPairwiseCor_df
#'   )
#'    
SitesToRegion <- function(sitesSubregion_df,
                          sitesAreOrdered = TRUE,
                          keepminPairwiseCor_df,
                          returnAllSites = FALSE, 
                          verbose = TRUE){
  
    if (!returnAllSites & (sum(sitesSubregion_df$subregion) == 0)) {
      
      if (verbose) {
        message(
          "No co-edited subregions found. Returning empty GRanges object."
        )
      }
      return(GRanges())
      
    }
    
    # Separate site into chromosome and position. Eg.
    #   "chr22:41327462" --> c("chr22", "41327462")
    orderedSites_df <- OrderSitesByLocation(
      sites_char = sitesSubregion_df$site,
      output = "dataframe"
    )
    
    # If input sites are ordered then just combine variable subregion to 
    #   original data, if not ordered then we have to merge two datasets 
    #   instead of easily combine.
    if (!sitesAreOrdered) {
      
      orderedSites_df <- merge(
        orderedSites_df, sitesSubregion_df,
        by = "site"
      )
      
    } else {
      orderedSites_df$subregion <- sitesSubregion_df$subregion
    }
    
    orderedSites_df <- merge(
      orderedSites_df, keepminPairwiseCor_df,
      by = "subregion",
      all.x = TRUE
    )
    
    orderedSites_df <- orderedSites_df[
      !is.na(orderedSites_df$keepminPairwiseCor) &
        orderedSites_df$keepminPairwiseCor == 1,
    ]
    
    if (nrow(orderedSites_df) == 0) {
      
      if (verbose) {
        message(
          "No co-edited subregions found. Returning empty GRanges object."
        )
      }
      return(GRanges())
      
    } 
    
    # Find start and end position for each subregion ###
    start_df <- aggregate(pos ~ chr + subregion, data = orderedSites_df, min)
    end_df <- aggregate(pos ~ chr + subregion, data = orderedSites_df, max)
    
    # Create GRanges
    GRanges(
      seqnames = start_df$chr,
      ranges = IRanges(
        start = start_df$pos,
        end = end_df$pos
      )
    )
  
}
