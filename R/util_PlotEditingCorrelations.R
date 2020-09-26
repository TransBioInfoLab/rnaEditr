#' Plotting correlations of RNA editing levels within a region.
#'
#' @description Plotting correlations of RNA editing levels within a region.
#'
#' @param region_gr A GRanges object of a region.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values on
#'   individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param ... Dots for additional internal arguments, see 
#'   \code{\link[corrplot]{corrplot}} for details.
#'
#' @return (Invisibly) returns a reordered correlation matrix.
#'
#' @importFrom corrplot corrplot
#'
#' @export
#' @keywords internal
#'
#' @examples
#'   data(rnaedit_df)
#'   
#'   genes_gr <- TransformToGR(
#'     genes_char = c("PHACTR4", "CCR5", "METTL7A"),
#'     type = "symbol",
#'     genome = "hg19"
#'   )
#'   
#'   exm_regions <- AllCoeditedRegions(
#'     regions_gr = genes_gr,
#'     rnaEditMatrix = rnaedit_df,
#'     output = "GRanges",
#'     method = "spearman"
#'   )
#'   
#'   PlotEditingCorrelations(
#'     region_gr = exm_regions[1],
#'     rnaEditMatrix = rnaedit_df
#'   )
#'    
PlotEditingCorrelations <- function(region_gr, rnaEditMatrix, ...){
  
  region_df <- as.data.frame(
    region_gr,
    stringsAsFactors = FALSE
  )
  
  region_df$seqnames <- as.character(region_df$seqnames)
  
  rnaedit_df <- GetSitesLocations(
    region_df = region_df,
    rnaEditMatrix = rnaEditMatrix,
    output = "locationsAndValues"
  )
  
  rnaedit_t_mat <- t(rnaedit_df)
  
  corr <- cor(
    rnaedit_t_mat,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  corrplot(corr, method="number", number.cex = 1, tl.cex = 0.7, ...)
  
}
