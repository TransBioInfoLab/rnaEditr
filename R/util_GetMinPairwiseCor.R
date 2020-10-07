#' Calculate minimum pairwise correlation for sub-regions.
#'
#' @description Filter the contiguous co-edited subregions found from
#'   \code{\link{FindCorrelatedRegions}}, by calculating pairwise correlations 
#'   and then selecting subregions passing the minimum correlation filter.
#'    
#' @param rnaEditCluster_mat A matrix of RNA editing level values on individual 
#'   sites, with row names as sample IDs and column names as site IDs in the 
#'   form of "chrAA:XXXXXXXX".
#' @param minPairCorr Minimum pairwise correlation coefficient of sites within 
#'   a cluster, used as a filter. To use this filter, set a number between -1 
#'   and 1 (defaults to 0.1). To turn it off, please set the number to -1.
#' @param probes_ls A list of regions with sites. Please note that probes in
#'   each list need to be ordered by their locations.
#' @param method Method for computing correlation. Defaults to 
#'   \code{"spearman"}.
#'
#' @return A list with a list of probes passing the minPairCorr and a data 
#'   frame with the following columns:
#'   \itemize{
#'     \item{\code{subregion} : }{index for each output contiguous co-edited
#'     region.}
#'     \item{\code{keepminPairwiseCor} : }{indicator for contiguous co-edited
#'     subregions, The regions with \code{keepminPairwiseCor = 1} passed the
#'     minimum correlation and will be returned as a contiguous co-edited
#'     subregion.}
#'     \item{\code{minPairwiseCor} : }{the minimum pairwise correlation of 
#'     sites within a subregion.}
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
#'   exm_sites <- list(
#'     "1" = c("chr1:28661656", "chr1:28661718", "chr1:28662148")
#'   )
#'   
#'   GetMinPairwiseCor(
#'     rnaEditCluster_mat = exm_data,
#'     minPairCorr = 0.1,
#'     probes_ls = exm_sites,
#'     method = "spearman"
#'   )
#'                
GetMinPairwiseCor <- function (rnaEditCluster_mat,
                               minPairCorr = 0.1,
                               probes_ls,
                               method = c("spearman", "pearson")) {
  
    method <- match.arg(method)
    
    # Calculate min pairwise correlation for each subregion
    # If users want all regions regardless of correlation, return all regions
    if (minPairCorr == -1) {
      
      nProbes <- length(probes_ls)
      minPairwiseCor_num <- rep(NA_real_, nProbes)
      keepminPairwiseCor <- rep(TRUE, nProbes)
      
    } else {
      
      minPairwiseCor_num <- vapply(
        X = probes_ls,
        FUN = function(probes){
          dat <- rnaEditCluster_mat[, probes]
          cor_matrix <- cor(
            dat, method = method, use = "pairwise.complete.obs"
          )
          min(cor_matrix)
        },
        FUN.VALUE = numeric(1)
      )
      
      keepminPairwiseCor <- minPairwiseCor_num > minPairCorr
      
    }
    
    # Create data frame of coedited regions
    df <- data.frame(
      subregion = names(probes_ls),
      keepminPairwiseCor = as.numeric(keepminPairwiseCor),
      minPairwiseCor = minPairwiseCor_num,
      stringsAsFactors = FALSE
    )
    
    if (all(keepminPairwiseCor == FALSE)) {
      probes_ls <- NULL
    } else {
      probes_ls <- probes_ls[keepminPairwiseCor]
    }
    
    # Return
    return(
      list(
        keepminPairwiseCor_df = df,
        probesFiltered_ls     = probes_ls
      )
    )
  
}
