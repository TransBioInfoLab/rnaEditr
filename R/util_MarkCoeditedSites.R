#' Mark RNA editing sites in contiguous and co-edited region.
#' 
#' @description Mark RNA editing sites in contiguous and co-edited region by
#'   selecting sites for which \code{r_drop} values calculated from inner 
#'   function \code{\link{CreateRdrop}} is greater than \code{rDropThresh_num}.
#'
#' @param rnaEditCluster_mat A matrix of RNA editing level values on individual 
#'   sites, with row names as sample IDs and column names as site IDs in the 
#'   form of "chrAA:XXXXXXXX".
#' @param rDropThresh_num Threshold for minimum correlation between RNA editing
#'   levels of one site and the mean RNA editing levels of the rest of the 
#'   sites. Please set a number between 0 and 1. Defaults to 0.4.
#' @param method Method for computing correlation. Defaults to 
#'   \code{"spearman"}.
#' @param minEditFreq Threshold for minimum percentage of edited samples for a 
#'   given site. The \code{r_drop} value of the sites with frequency lower than 
#'   \code{minEditFreq} will be set as NA. Please set a number between 0 and 1. 
#'   Defaults to 0.05.
#' @param verbose Should messages and warnings be displayed? Defaults to TRUE.
#' 
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{r_drop} : }{The correlation between RNA editing levels of 
#'     one site and the mean RNA editing levels of the rest of the sites.}
#'     \item{\code{keep} : }{indicator for co-edited sites, The sites with
#'     \code{keep = 1} belong to the contiguous and co-edited region.}
#'     \item{\code{keep_contiguous} : }{contiguous co-edited region number}
#'   }
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{keep} : }{indicator for co-edited sites, The sites with
#'     \code{keep = 1} belong to the contiguous and co-edited region.}
#'     \item{\code{ind} : }{index for the sites.}
#'     \item{\code{r_drop} : }{the correlation between RNA editing levels of 
#'     one site and the mean RNA editing levels of the rest of the sites.}
#'   }
#'
#' @details \code{r_drop} statistic is used to identify co-edited sites. An 
#'   outlier site (\code{keep = 0}) in a genomic region typically has low 
#'   correlation with the rest of the sites in a genomic region. The sites with 
#'   \code{r_drop} value greater than \code{rDropThresh_num} are marked to have 
#'   \code{keep = 1}. Please see \code{\link{CreateRdrop}} for more details.
#'
#' @seealso \code{\link{CreateRdrop}} 
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
#'   MarkCoeditedSites(
#'     rnaEditCluster_mat = exm_data,
#'     method = "spearman"
#'   )
#'    
MarkCoeditedSites <- function (rnaEditCluster_mat,
                               rDropThresh_num = 0.4,
                               method = c("spearman", "pearson"),
                               minEditFreq = 0.05,
                               verbose = TRUE) {
  
    method <- match.arg(method)
    
    # Calculate r_drop
    clusterRdrop_df <- CreateRdrop(
      data = rnaEditCluster_mat,
      method = method,
      minEditFreq = minEditFreq
    )
    
    sites_char <- clusterRdrop_df$site
    
    # Drop sites with missing r_drop or r_drop < rDropThresh_num
    dropSites_char <- sites_char[
      is.na(clusterRdrop_df$r_drop) |
        clusterRdrop_df$r_drop < rDropThresh_num
    ]
  
    # Create Output Data Frame
    data.frame(
      site = sites_char,
      keep = ifelse(sites_char %in% dropSites_char, 0, 1), #(drop=0, keep=1)
      ind = seq_len(ncol(rnaEditCluster_mat)),
      r_drop = clusterRdrop_df$r_drop,
      stringsAsFactors = FALSE
    )
  
}
