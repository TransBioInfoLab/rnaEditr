#' Find contiguous co-edited subregions.
#'
#' @description Find contiguous co-edited subregions based on the output file
#'   from function \code{\link{MarkCoeditedSites}}.
#'
#' @param sites_df An output data frame from function \code{MarkCoeditedSites},
#'   with variables \code{site, keep, ind, r_drop}. Please see
#'   \code{\link{MarkCoeditedSites}} for details.
#' @param featureType Feature type, Defaults to \code{"site"}.
#' @param minSites_int An integer indicates the minimum number of sites to be
#'   considered a contiguous co-edited region.
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{subregion} : }{index for each output contiguous co-edited
#'     region.}
#'   }
#'
#' @importFrom bumphunter getSegments
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
#'   FindCorrelatedRegions(
#'     sites_df = exm_sites,
#'     featureType = "site"
#'   )
#'    
FindCorrelatedRegions <- function(sites_df,
                                  featureType = c("site", "cpg"),
                                  minSites_int = 3){

    featureType <- match.arg(featureType)
    
    # Get contiguous regions of sites
    contiguousRegion_ls <- getSegments(sites_df$keep, cutoff = 1)
    nSegs_int <- length(contiguousRegion_ls$upIndex)
  
    if (nSegs_int > 0){
  
      # Select segments with number of sites >= minSites
      contiguous_int <- lengths(contiguousRegion_ls$upIndex)
      contiguousMinSites_idx <- which(contiguous_int >= minSites_int)
      nSegsMinSites_int <- length(contiguousMinSites_idx)
  
      # Create output dataframe with sites and contiguous coedited subregion
      #   number
      ind<-NULL
      
      if (nSegsMinSites_int > 0){
        
        inner_ls <- lapply(
          seq_len(nSegsMinSites_int),
          function(u){
            
            data.frame(
              site = subset(
                sites_df,
                ind %in% contiguousRegion_ls$upIndex[[
                  contiguousMinSites_idx[u]
                ]],
                select = featureType
              ),
              subregion = rep(
                u, length(
                  contiguousRegion_ls$upIndex[[contiguousMinSites_idx[u]]]
                )
              )
            )
          }
        )
        
        contiguousRegionsSites <- do.call(rbind, inner_ls)
        
        
      } else {
        
        contiguousRegionsSites <- cbind(
          as.data.frame(sites_df[, featureType], stringsAsFactors = FALSE),
          rep(0,length(sites_df[, featureType]))
        )
        
      }
  
    } else {
      
        contiguousRegionsSites <- cbind(
          as.data.frame(sites_df[, featureType], stringsAsFactors = FALSE),
          rep(0,length(sites_df[, featureType]))
        )
        
    }
  
    colnames(contiguousRegionsSites) <- c(featureType,"subregion")
    
    contiguousRegionsSites
  
}
