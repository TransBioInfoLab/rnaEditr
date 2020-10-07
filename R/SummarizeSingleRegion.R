#' Summarizes RNA editing levels from multiple sites in a single region.
#'
#' @description Summarizes RNA editing levels from multiple sites in an input
#'   region.
#'
#' @param region_df A data frame with the input genomic region. Please make 
#'   sure columns \code{seqnames}, \code{start}, and \code{end} are included in
#'   the data frame.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values 
#'   for individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param selectMethod Method for summarizing regions. Available options are 
#'   \code{"MaxSites", "MeanSites", "MedianSites", "PC1Sites"}. Please see
#'   \code{\link{RegionSummaryMethod}} for more details.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return A named numeric vector of summarized RNA editing levels with sample
#'   IDs as column names.
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
#'   SummarizeSingleRegion(
#'     region_df = exm_region,
#'     rnaEditMatrix = rnaedit_df
#'   )[1:3]
#'
SummarizeSingleRegion <- function(region_df,
                                  rnaEditMatrix,
                                  selectMethod = MedianSites, ...){
  
    rnaEdit_df <- GetSitesLocations(
      region_df = region_df,
      rnaEditMatrix = rnaEditMatrix,
      output = "locationsAndValues"
    )
  
    selectMethod(rnaEditMatrix = rnaEdit_df, ...)

}


