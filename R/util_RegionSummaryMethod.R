#' Methods to summarize RNA editing levels from multiple sites within a 
#'   single region.
#' 
#' @description Summarize RNA editing sites in a single region by taking 
#'   maximum, mean, median or first principal component.
#'
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values on
#'   individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return A named numeric vector of summarized RNA editing levels with sample
#'   IDs as names.
#' 
#' @importFrom stats median
#' @importFrom stats prcomp
#' 
#'
#' @examples
#'   data(rnaedit_df)
#'   MedianSites(rnaEditMatrix = rnaedit_df)[1:3]
#'
#' @name RegionSummaryMethod
NULL

#' @rdname RegionSummaryMethod
#' @export
#' @keywords internal
MaxSites <- function(rnaEditMatrix, ...){
  
    apply(
      X = rnaEditMatrix,
      MARGIN = 2,
      FUN = function(x){
        max(x, na.rm = TRUE, ...)
      }
    )
  
}


#' @rdname RegionSummaryMethod
#' @export
#' @keywords internal
MeanSites <- function(rnaEditMatrix, ...){
  
    colMeans(
      x = rnaEditMatrix,
      na.rm = TRUE,
      ...
    )
  
}


#' @rdname RegionSummaryMethod
#' @export
#' @keywords internal
MedianSites <- function(rnaEditMatrix, ...){
  
    apply(
      X = rnaEditMatrix,
      MARGIN = 2,
      FUN = function(x){
        median(x, na.rm = TRUE, ...)
      }
    )
  
}


#' @rdname RegionSummaryMethod
#' @export
#' @keywords internal
PC1Sites <- function(rnaEditMatrix, ...){
  
    site_prcomp <- prcomp(t(rnaEditMatrix), scale. = TRUE, ...)
    site_prcomp$x[, 1, drop = TRUE]
  
}

