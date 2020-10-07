#' Convert RNA editing matrix into a special data frame with class
#'   \code{rnaEdit_df}.
#'
#' @description Convert RNA editing matrix to a special data frame with class 
#'   \code{rnaEdit_df}, which is then used to identify differentially co-edited
#'   regions with function \code{\link{TestAssociations}}.
#' 
#' @param rnaEditMatrix A matrix of RNA editing level values on individual 
#'   sites, with row names as site IDs in the form of "chrAA:XXXXXXXX", and 
#'   column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' 
#' @return A dataset of class \code{rnaEdit_df}, includes variables 
#'   \code{seqnames, start, end, width} and summarized RNA editing levels in 
#'   each sample.
#' 
#' @export
#'
#' @seealso \code{\link{TransformToGR}}, \code{\link{AllCloseByRegions}}, 
#'   \code{\link{AllCoeditedRegions}}, \code{\link{SummarizeAllRegions}}, 
#'   \code{\link{TestAssociations}}, \code{\link{AnnotateResults}}
#'   
#' @examples
#'   data(rnaedit_df)
#'   CreateEditingTable(rnaEditMatrix = rnaedit_df)[1:3, 1:5]
#'   
CreateEditingTable <- function(rnaEditMatrix){
  
    sites_mat <- do.call(
      rbind,
      strsplit(
        row.names(rnaEditMatrix),
        split = ":"
      )
    )
    
    # We call the chromosome column "seqnames" to match the default GRanges
    #   output.
    metaCols_df <- data.frame(
      seqnames = sites_mat[, 1],
      start    = as.integer(sites_mat[, 2]),
      end      = as.integer(sites_mat[, 2]),
      width    = 1L,
      stringsAsFactors = FALSE
    )
    
    dat_df <- cbind(metaCols_df, rnaEditMatrix)
    row.names(dat_df) <- NULL
    
    class(dat_df) <- c("rnaEdit_df", "data.frame")
    
    dat_df
  
}
