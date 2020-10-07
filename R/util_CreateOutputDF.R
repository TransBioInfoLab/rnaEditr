#' Create output data in the format of data frame.
#'
#' @description Output all the contiguous co-edited subregions found by
#'   \code{\link{FindCorrelatedRegions}} function and filtered by
#'   \code{\link{GetMinPairwiseCor}} function.
#'    
#' @param keepSites_df An output data frame from function 
#'   \code{MarkCoeditedSites}, with variables \code{site, keep, ind, r_drop}. 
#'   Please see \code{\link{MarkCoeditedSites}} for details.
#' @param keepContiguousSites_df An output data frame from function
#'   \code{FindCorrelatedRegions} with variables \code{site, subregion}. Please
#'   see \code{\link{FindCorrelatedRegions}} for details.
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
#' @return A data frame with following columns:
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{chr} : }{chromosome number.}
#'     \item{\code{pos} : }{genomic location.}
#'     \item{\code{r_drop} : }{the correlation between RNA editing levels of 
#'     one site and the mean RNA editing levels of the rest of the sites.}
#'     \item{\code{keep} : }{indicator for co-edited sites, The sites with
#'     \code{keep = 1} belong to the contiguous and co-edited region.}
#'     \item{\code{keep_contiguous} : }{contiguous co-edited region number.}
#'     \item{\code{regionMinPairwiseCor} : }{the minimum pairwise correlation 
#'     between sites within a sub-region.}
#'     \item{\code{keep_regionMinPairwiseCor} : }{indicator for contiguous
#'     co-edited subregions, The regions with \code{keepminPairwiseCor = 1}
#'     are the ones that passed \code{regionMinPairwiseCor} filter and will be 
#'     returned as a contiguous co-edited sub-region.}
#'   }
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
#'   exm_probes <- split(
#'     x = exm_regions$site,
#'     f = exm_regions$subregion
#'   )
#'   
#'   exm_cor <- GetMinPairwiseCor(
#'     rnaEditCluster_mat = exm_data,
#'     minPairCorr = 0.1,
#'     probes_ls = exm_probes,
#'     method = "spearman"
#'   )
#'   
#'   CreateOutputDF(
#'     keepSites_df = exm_sites,
#'     keepContiguousSites_df = exm_regions,
#'     keepminPairwiseCor_df = exm_cor$keepminPairwiseCor_df
#'   )
#'    
CreateOutputDF <- function(keepSites_df,
                           keepContiguousSites_df,
                           keepminPairwiseCor_df,
                           returnAllSites = FALSE,
                           verbose = TRUE){

    if (!returnAllSites & (sum(keepContiguousSites_df$subregion) == 0)) {
  
      if (verbose) {
        message(
          "No co-edited subregions found. Returning empty GRanges object."
        )
      }
      return(data.frame())
  
    }
    
    output_df <- merge(
      x = keepSites_df,
      y = keepContiguousSites_df,
      by = "site",
      all.x = TRUE
    )
    
    output_df[is.na(output_df)] <- 0
    
    # to make sure final dataset does not have dis-continuous subregion 
    #   numbers, we will only include subregions with "keepminPairwiseCor = 1"
    #   here, and then re-number subregions.
    keepminPairwiseCor_df <- keepminPairwiseCor_df[
      keepminPairwiseCor_df$keepminPairwiseCor == 1,
    ]
    
    keepminPairwiseCor_df$subregion <- seq_len(
      nrow(keepminPairwiseCor_df)
    )
    
    output_df <- merge(
      output_df, keepminPairwiseCor_df,
      by = "subregion",
      all.x = TRUE
    )
    
    output_df <- output_df[
      !is.na(output_df$keepminPairwiseCor) &
        output_df$keepminPairwiseCor == 1,
    ]
    
    # dataset is not ordered after line 108-112 where we merge by "subregion" 
    #   only, but since each "subregion" contains more than one site, so their
    #   index number well not be ordered anymore. So we are ordering here.
    output_df <- output_df[
      order(output_df$subregion, output_df$ind),
    ]
    
    if (nrow(output_df) == 0){
      
      if (verbose) {
        message(
          "No co-edited subregions found. Returning empty GRanges object."
        )
      }
      return(data.frame())
      
    }
    
    output2_df <- OrderSitesByLocation(output_df$site)
    
    data.frame(
      site = output_df$site,
      chr = output2_df$chr,
      pos = output2_df$pos,
      r_drop = output_df$r_drop,
      keep = output_df$keep,
      keep_contiguous = output_df$subregion,
      regionMinPairwiseCor = output_df$minPairwiseCor,
      keep_regionMinPairwiseCor = output_df$keepminPairwiseCor
    )
  
}
