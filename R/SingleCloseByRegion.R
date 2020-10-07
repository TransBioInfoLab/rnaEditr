#' Extracts clusters of RNA editing sites located closely in a single 
#'   genomic region.
#'
#' @description Extracts clusters of RNA editing sites located closely in an 
#' input genomic region.
#'
#' @param region_df A data frame with the input genomic region. Please make 
#'   sure columns \code{seqnames}, \code{start}, and \code{end} are included in 
#'   the data frame.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values 
#'   for individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param maxGap An integer, genomic locations within \code{maxGap} from each
#'   other are placed into the same cluster. Defaults to 50.
#' @param minSites An integer, minimum number of edited sites for a
#'   cluster to be selected for output. Defaults to 3.
#'
#' @return A GRanges object containing genomic locations of RNA editing sites
#'   located closely within the single input pre-defined genomic region.
#'
#' @details The algorithm of this function is based on the
#'   \code{\link[bumphunter]{clusterMaker}} function in the \code{bumphunter}
#'   R package. Each cluster is essentially a group of sites such that
#'   two consecutive sites in the cluster are separated by less than
#'   \code{maxGap}.
#'
#' @importFrom bumphunter clusterMaker
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @export
#' @keywords internal
#'
#' @examples
#'   data(rnaedit_df)
#'   
#'   exm_region <- data.frame(
#'     seqnames = "chr1",
#'     start =  28691093,
#'     end = 28826881, 
#'     stringsAsFactors = FALSE
#'   )
#'   
#'   SingleCloseByRegion(
#'     region_df = exm_region,
#'     rnaEditMatrix = rnaedit_df,
#'     maxGap = 50,
#'     minSites = 3
#'   )  
#'   
SingleCloseByRegion <- function(region_df,
                                rnaEditMatrix,
                                maxGap = 50,
                                minSites = 3){
  
    sitesOrdered_df <- GetSitesLocations(
      region_df = region_df,
      rnaEditMatrix = rnaEditMatrix,
      output = "locationsOnly"
    )
    
    if (is.null(sitesOrdered_df)) {
     
      return(NULL)
      
    } else {
      
      # Find close by clusters
      sitesOrdered_df$cluster <- clusterMaker(
        chr = sitesOrdered_df$chr,
        pos = sitesOrdered_df$pos,
        maxGap = maxGap
      )
      
      # Count number of sites in each cluster
      clusterLength <- data.frame(table(sitesOrdered_df$cluster))
      
      # Include clusters with number of sites >= minSites
      clusterLength_flt <- clusterLength[clusterLength$Freq >= minSites, ]
      
      sitesOrdered_flt_df <- sitesOrdered_df[
        sitesOrdered_df$cluster %in% clusterLength_flt$Var1,
      ]
      
      if (nrow(sitesOrdered_flt_df) == 0) {
        
        return(NULL)
        
      } else {
        
        # Find start and end position for each cluster
        start_df <- aggregate(
          pos ~ chr + cluster,
          data = sitesOrdered_flt_df,
          min)
        
        end_df <- aggregate(
          pos ~ chr + cluster,
          data = sitesOrdered_flt_df,
          max)
        
        # Convert data frame to GRanges
        GRanges(
          seqnames = start_df$chr,
          ranges = IRanges(
            start = start_df$pos,
            end = end_df$pos
          )
        )
        
      }
    }
}
